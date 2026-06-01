#include "da_filter.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <chrono>

#ifdef __linux__
#include <sched.h>
#include <pthread.h>
#endif

// SIMD intrinsics
#if defined(__SSE2__) || defined(__AVX__)
#include <immintrin.h>
#include <xmmintrin.h>
#endif

// ARM NEON optimized functions (from ThreadedDSP compatibility layer)
#if defined(__arm__) || defined(__aarch64__)

inline IppStatus ippsAdd_32f_I(const Ipp32f* pSrc, Ipp32f* pSrcDst, int len) {
    while(len >= 4) { 
        vst1q_f32(pSrcDst, vaddq_f32(vld1q_f32(pSrcDst), vld1q_f32(pSrc))); 
        pSrcDst += 4; pSrc += 4; len -= 4; 
    }
    while(len > 0) { *pSrcDst++ += *pSrc++; len--; }
    return ippStsNoErr; 
}

inline IppStatus ippsMaxAbs_32f(const Ipp32f* pSrc, int len, Ipp32f* pMaxAbs) {
    *pMaxAbs = 0; 
    while(len--) { 
        float abs_val = std::abs(*pSrc); 
        if (abs_val > *pMaxAbs) *pMaxAbs = abs_val; 
        pSrc++; 
    }
    return ippStsNoErr; 
}

inline IppStatus ippsAddProduct_32fc(const Ipp32fc* pSrc1, const Ipp32fc* pSrc2, Ipp32fc* pSrcDst, int len) {   
    assert(len % 4 == 0);
    float32x4x2_t A, B, C;
    while (len >= 4) {
        A = vld2q_f32((const float*)pSrcDst);    
        B = vld2q_f32((const float*)pSrc1);    
        C = vld2q_f32((const float*)pSrc2);
        // Complex multiply-accumulate: (a+bi)(c+di) = (ac-bd) + (ad+bc)i
        A.val[0] = vmlaq_f32(A.val[0], B.val[0], C.val[0]);  // real += real1 * real2
        A.val[0] = vmlsq_f32(A.val[0], B.val[1], C.val[1]);  // real -= imag1 * imag2
        A.val[1] = vmlaq_f32(A.val[1], B.val[0], C.val[1]);  // imag += real1 * imag2
        A.val[1] = vmlaq_f32(A.val[1], B.val[1], C.val[0]);  // imag += imag1 * real2
        vst2q_f32((float*)pSrcDst, A);
        pSrcDst += 8; pSrc1 += 8; pSrc2 += 8; len -= 4;
    }
    return ippStsNoErr;
}

inline IppStatus ippsFFTGetSize_R_32f(int order, int flag, IppHintAlgorithm hint, int*pSpecSize, int* pSpecBufferSize, int* pBufferSize) {
    *pSpecSize = 0; *pSpecBufferSize = 0; *pBufferSize = (int)sizeof(float32_t) * ((1 << order) + 2); 
    return ippStsNoErr; 
}

inline IppStatus ippsFFTInit_R_32f(IppsFFTSpec_R_32f** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer) { 
    *ppFFTSpec = (IppsFFTSpec_R_32f *)ne10_fft_alloc_r2c_float32(1 << order); 
    return ippStsNoErr; 
}

inline IppStatus ippsFFTFwd_RToPerm_32f(Ipp32f* pSrc, Ipp32f* pDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer) {
    ne10_fft_r2c_cfg_float32_t FFT = (ne10_fft_r2c_cfg_float32_t)pFFTSpec;
    float32_t *Temp = (float32_t *)pBuffer;
    int32_t N = FFT->nfft;
    ne10_fft_r2c_1d_float32_c((ne10_fft_cpx_float32_t *)Temp, pSrc, FFT);
    pDst[0] = Temp[0]; pDst[1] = Temp[N]; 
    memmove(pDst + 2, Temp + 2, (N - 2) * sizeof(Ipp32f));
    return ippStsNoErr; 
}

inline IppStatus ippsFFTInv_PermToR_32f(const Ipp32f* pSrc, Ipp32f* pDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer) {
    ne10_fft_r2c_cfg_float32_t FFT = (ne10_fft_r2c_cfg_float32_t)pFFTSpec;
    float32_t *Temp = (float32_t *)pBuffer;
    int32_t N = FFT->nfft;
    Temp[0] = pSrc[0]; Temp[1] = 0; Temp[N] = pSrc[1]; Temp[N + 1] = 0; 
    memmove(Temp + 2, (void *)(pSrc + 2), (N - 2) * sizeof(Ipp32f));
    ne10_fft_c2r_1d_float32_c(pDst, (ne10_fft_cpx_float32_t *)Temp, FFT);
    return ippStsNoErr; 
}

inline IppStatus ippsMulC_32f_I(Ipp32f val, Ipp32f* pSrcDst, int len) {
    while(len >= 4) {
        float32x4_t v = vld1q_f32(pSrcDst);
        vst1q_f32(pSrcDst, vmulq_n_f32(v, val));
        pSrcDst += 4; len -= 4;
    }
    while(len > 0) { *pSrcDst++ *= val; len--; }
    return ippStsNoErr;
}

inline IppStatus ippsMul_32f_I(const Ipp32f* pSrc, Ipp32f* pSrcDst, int len) {
    while(len >= 4) {
        float32x4_t a = vld1q_f32(pSrcDst);
        float32x4_t b = vld1q_f32(pSrc);
        vst1q_f32(pSrcDst, vmulq_f32(a, b));
        pSrcDst += 4; pSrc += 4; len -= 4;
    }
    while(len > 0) { *pSrcDst++ *= *pSrc++; len--; }
    return ippStsNoErr;
}

#endif

DAFilter::DAFilter() {
    stats_.last_process_time = std::chrono::high_resolution_clock::now();
}

DAFilter::~DAFilter() {
    shutdown();
}

bool DAFilter::initialize(int block_size, int num_input_channels, int num_output_channels, 
                         int num_threads, int overlap_factor, WindowType window_type) {
    if (running_) {
        return false; // Already initialized
    }
    
    // Validate parameters
    if (block_size < 16 || block_size > MAX_BLOCK_SIZE || 
        num_input_channels < 0 || num_input_channels > MAX_CHANNELS ||
        num_output_channels < 0 || num_output_channels > MAX_CHANNELS ||
        overlap_factor < 2 || overlap_factor > 8) {
        return false;
    }
    
    // Set core parameters
    block_size_ = block_size;
    num_input_channels_ = num_input_channels;
    num_output_channels_ = num_output_channels;
    window_type_ = window_type;
    
    // Calculate FFT size (next power of 2 that accommodates block + overlap)
    fft_size_ = calculate_optimal_fft_size(block_size, overlap_factor);
    overlap_size_ = fft_size_ - block_size;
    
    // Set number of threads (default to hardware concurrency)
    if (num_threads <= 0) {
        num_threads_ = std::max(1, (int)std::thread::hardware_concurrency());
    } else {
        num_threads_ = std::min(num_threads, MAX_THREADS);
    }
    
    try {
        // Allocate audio buffers with improved memory layout
        size_t input_time_size = num_input_channels * (overlap_size_ + block_size_);
        size_t output_time_size = num_output_channels * (overlap_size_ + block_size_);
        size_t input_freq_size = num_input_channels * fft_size_ * 2; // Complex
        size_t output_freq_size = num_output_channels * fft_size_ * 2; // Complex
        size_t overlap_buffer_size = num_output_channels * overlap_size_;
        size_t temp_buffer_size = num_threads_ * fft_size_ * 2; // Complex per thread
        
        input_time_buffer_ = std::make_unique<float[]>(input_time_size);
        output_time_buffer_ = std::make_unique<float[]>(output_time_size);
        input_freq_buffer_ = std::make_unique<float[]>(input_freq_size);
        output_freq_buffer_ = std::make_unique<float[]>(output_freq_size);
        overlap_buffer_ = std::make_unique<float[]>(overlap_buffer_size);
        window_function_ = std::make_unique<float[]>(fft_size_);
        temp_buffers_ = std::make_unique<float[]>(temp_buffer_size);
        
        // Clear all buffers
        std::memset(input_time_buffer_.get(), 0, input_time_size * sizeof(float));
        std::memset(output_time_buffer_.get(), 0, output_time_size * sizeof(float));
        std::memset(input_freq_buffer_.get(), 0, input_freq_size * sizeof(float));
        std::memset(output_freq_buffer_.get(), 0, output_freq_size * sizeof(float));
        std::memset(overlap_buffer_.get(), 0, overlap_buffer_size * sizeof(float));
        std::memset(temp_buffers_.get(), 0, temp_buffer_size * sizeof(float));
        
        // Initialize window function
        generate_window_function();
        
        // Initialize FFT contexts for each thread
        initialize_fft_contexts();
        
        // Initialize statistics
        std::memset(stats_.peak_input, 0, sizeof(stats_.peak_input));
        std::memset(stats_.peak_output, 0, sizeof(stats_.peak_output));
        stats_.blocks_processed = 0;
        stats_.cpu_load = 0.0f;
        
        // Create worker threads
        running_ = true;
        threads_ready_ = 0;
        threads_completed_ = 0;
        
        worker_threads_.reserve(num_threads_);
        for (int i = 0; i < num_threads_; i++) {
            worker_threads_.emplace_back(
                std::make_unique<std::thread>(&DAFilter::worker_thread_proc, this, i)
            );
        }
        
        // Wait for all threads to be ready
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        
        return true;
        
    } catch (const std::exception& e) {
        shutdown();
        return false;
    }
}

void DAFilter::shutdown() {
    if (!running_) {
        return;
    }
    
    // Signal threads to stop
    running_ = false;
    start_processing_.notify_all();
    
    // Wait for all worker threads to finish
    for (auto& thread : worker_threads_) {
        if (thread && thread->joinable()) {
            thread->join();
        }
    }
    worker_threads_.clear();
    
    // Cleanup FFT contexts
    cleanup_fft_contexts();
    
    // Reset state
    filter_banks_.clear();
    active_bank_index_ = -1;
    target_bank_index_ = -1;
}

void DAFilter::process_block(const float* const* inputs, float* const* outputs) {
    if (!running_ || processing_) {
        return; // Not ready or already processing
    }
    
    auto process_start = std::chrono::high_resolution_clock::now();
    
    // Copy input data to time buffers (with overlap preservation)
    for (int ch = 0; ch < num_input_channels_; ch++) {
        float* time_buf = get_input_time_buffer(ch);
        // Shift existing overlap data
        std::memmove(time_buf, time_buf + block_size_, overlap_size_ * sizeof(float));
        // Copy new input block
        std::memcpy(time_buf + overlap_size_, inputs[ch], block_size_ * sizeof(float));
    }
    
    // Signal worker threads to start processing
    {
        std::lock_guard<std::mutex> lock(process_mutex_);
        processing_ = true;
        threads_ready_ = 0;
        threads_completed_ = 0;
    }
    start_processing_.notify_all();
    
    // Wait for all threads to complete
    std::unique_lock<std::mutex> lock(process_mutex_);
    processing_complete_.wait(lock, [this] { return threads_completed_ == num_threads_; });
    processing_ = false;
    
    // Copy output data from time buffers
    for (int ch = 0; ch < num_output_channels_; ch++) {
        float* time_buf = get_output_time_buffer(ch);
        std::memcpy(outputs[ch], time_buf + overlap_size_, block_size_ * sizeof(float));
    }
    
    // Update statistics
    auto process_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(process_end - process_start);
    float process_time_ms = duration.count() / 1000.0f;
    float block_time_ms = (block_size_ * 1000.0f) / 48000.0f; // Assume 48kHz for now
    stats_.cpu_load = 0.9f * stats_.cpu_load + 0.1f * (process_time_ms / block_time_ms) * 100.0f;
    stats_.blocks_processed++;
    stats_.last_process_time = process_end;
}

void DAFilter::worker_thread_proc(int thread_id) {
    // Set thread priority on Linux
#ifdef __linux__
    struct sched_param param{};
    param.sched_priority = 50; // High priority but not maximum
    pthread_setschedparam(pthread_self(), SCHED_FIFO, &param);
#endif

    ippSetDenormAreZeros(true);
    
    std::unique_lock<std::mutex> lock(process_mutex_);
    
    while (running_) {
        // Wait for processing to start
        start_processing_.wait(lock, [this] { return processing_ || !running_; });
        
        if (!running_) break;
        
        lock.unlock();
        
        // Divide work among threads
        int channels_per_thread = num_input_channels_ / num_threads_;
        int start_input_ch = thread_id * channels_per_thread;
        int end_input_ch = (thread_id == num_threads_ - 1) ? num_input_channels_ : (thread_id + 1) * channels_per_thread;
        
        // Phase 1: Forward FFT on input channels
        process_input_channels(thread_id, start_input_ch, end_input_ch);
        
        // Synchronization barrier
        lock.lock();
        if (++threads_ready_ == num_threads_) {
            threads_ready_ = 0;
            start_processing_.notify_all();
        } else {
            start_processing_.wait(lock, [this] { return threads_ready_ == 0; });
        }
        lock.unlock();
        
        // Phase 2: Filter processing (frequency domain multiply-accumulate)
        if (active_bank_index_ >= 0) {
            int filters_per_thread = filter_banks_[active_bank_index_].filters.size() / num_threads_;
            int start_filter = thread_id * filters_per_thread;
            int end_filter = (thread_id == num_threads_ - 1) ? 
                filter_banks_[active_bank_index_].filters.size() : (thread_id + 1) * filters_per_thread;
            
            process_filters(thread_id, start_filter, end_filter);
        }
        
        // Synchronization barrier
        lock.lock();
        if (++threads_ready_ == num_threads_) {
            threads_ready_ = 0;
            start_processing_.notify_all();
        } else {
            start_processing_.wait(lock, [this] { return threads_ready_ == 0; });
        }
        lock.unlock();
        
        // Phase 3: Inverse FFT on output channels
        channels_per_thread = num_output_channels_ / num_threads_;
        int start_output_ch = thread_id * channels_per_thread;
        int end_output_ch = (thread_id == num_threads_ - 1) ? num_output_channels_ : (thread_id + 1) * channels_per_thread;
        
        process_output_channels(thread_id, start_output_ch, end_output_ch);
        
        // Final synchronization
        lock.lock();
        if (++threads_completed_ == num_threads_) {
            processing_complete_.notify_one();
        }
    }
}

void DAFilter::process_input_channels(int thread_id, int start_channel, int end_channel) {
    FFTContext& fft_ctx = fft_contexts_[thread_id];
    float* temp_buf = get_temp_buffer(thread_id);
    
    for (int ch = start_channel; ch < end_channel; ch++) {
        float* time_buf = get_input_time_buffer(ch);
        float* freq_buf = get_input_freq_buffer(ch);
        
        // Apply window function and perform forward FFT
        for (int i = 0; i < fft_size_; i++) {
            temp_buf[i] = time_buf[i] * window_function_[i];
        }
        
        ippsFFTFwd_RToPerm_32f(temp_buf, freq_buf, fft_ctx.fft_spec, fft_ctx.fft_buffer);
        
        // Update peak level
        float peak;
        ippsMaxAbs_32f(time_buf + overlap_size_, block_size_, &peak);
        stats_.peak_input[ch] = 0.99f * stats_.peak_input[ch] + 0.01f * peak;
    }
}

void DAFilter::process_filters(int thread_id, int start_filter, int end_filter) {
    // Handle crossfading between filter banks
    bool is_crossfading = (bank_mode_ == FilterBankMode::CROSSFADE && 
                          active_bank_index_ != target_bank_index_ &&
                          target_bank_index_ >= 0 && target_bank_index_ < filter_banks_.size());
    
    if (active_bank_index_ < 0 || active_bank_index_ >= filter_banks_.size()) {
        return;
    }
    
    // Calculate crossfade progress (0.0 = old bank only, 1.0 = new bank only)
    float crossfade_progress = 0.0f;
    if (is_crossfading) {
        crossfade_progress = get_crossfade_progress();
    }
    
    // Process old/active bank with fade-out
    if (!is_crossfading || crossfade_progress < 1.0f) {
        const FilterBank& active_bank = filter_banks_[active_bank_index_];
        float active_gain = is_crossfading ? (1.0f - crossfade_progress) : 1.0f;
        
        process_bank_filters(thread_id, start_filter, end_filter, active_bank, active_gain);
    }
    
    // Process new/target bank with fade-in (during crossfade)
    if (is_crossfading && crossfade_progress > 0.0f) {
        const FilterBank& target_bank = filter_banks_[target_bank_index_];
        float target_gain = crossfade_progress;
        
        process_bank_filters(thread_id, start_filter, end_filter, target_bank, target_gain);
        
        // Update crossfade progress (only in one thread to avoid race conditions)
        if (thread_id == 0) {
            update_crossfade_progress();
        }
    }
}

void DAFilter::process_bank_filters(int thread_id, int start_filter, int end_filter, 
                                   const FilterBank& bank, float bank_gain) {
    for (int f = start_filter; f < end_filter; f++) {
        if (f >= bank.filters.size()) break;
        
        const FilterConfig& filter = bank.filters[f];
        if (!filter.enabled || filter.input_channel < 0 || filter.output_channel < 0 ||
            filter.input_channel >= num_input_channels_ || filter.output_channel >= num_output_channels_) {
            continue;
        }
        
        float* input_freq = get_input_freq_buffer(filter.input_channel);
        float* output_freq = get_output_freq_buffer(filter.output_channel);
        
        // Get filter frequency response (in a real implementation, this would be pre-computed)
        float* filter_freq = get_filter_freq_response(f, bank);
        
        // Apply complex convolution with crossfade gain
        float effective_gain = filter.gain * bank_gain * bank.mix_level;
        if (effective_gain != 0.0f) {
            // Use IPP optimized complex multiply-accumulate if available, 
            // otherwise fall back to our custom SIMD implementation
#if defined(IPP_VERSION_MAJOR) || defined(__arm__) || defined(__aarch64__)
            if (effective_gain == 1.0f) {
                // Direct complex multiply-accumulate using IPP (from ThreadedDSP)
                ippsAddProduct_32fc((const Ipp32fc*)input_freq, (const Ipp32fc*)filter_freq, 
                                   (Ipp32fc*)output_freq, fft_size_);
            } else {
                // Apply gain and then accumulate
                complex_multiply_accumulate_optimized(input_freq, filter_freq, output_freq, 
                                                    fft_size_, effective_gain);
            }
#else
            // Use our custom optimized implementation
            complex_multiply_accumulate_optimized(input_freq, filter_freq, output_freq, 
                                                fft_size_, effective_gain);
#endif
        }
    }
}

// Optimized complex multiply-accumulate using SIMD and proper complex math
void DAFilter::complex_multiply_accumulate_optimized(const float* input, const float* filter, 
                                                    float* output, int complex_samples, float gain) {
#if defined(__arm__) || defined(__aarch64__)
    // ARM NEON optimized version
    const float* in_ptr = input;
    const float* filt_ptr = filter;
    float* out_ptr = output;
    
    // Process 4 complex samples at a time (8 floats)
    int simd_samples = complex_samples & ~3; // Round down to multiple of 4
    
    for (int i = 0; i < simd_samples; i += 4) {
        // Load 4 complex input samples (8 floats)
        float32x4x2_t in_complex = vld2q_f32(in_ptr);
        // Load 4 complex filter samples (8 floats)  
        float32x4x2_t filt_complex = vld2q_f32(filt_ptr);
        // Load 4 complex output samples (8 floats)
        float32x4x2_t out_complex = vld2q_f32(out_ptr);
        
        // Complex multiply: (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
        float32x4_t real_result = vmulq_f32(in_complex.val[0], filt_complex.val[0]); // a*c
        real_result = vmlsq_f32(real_result, in_complex.val[1], filt_complex.val[1]); // a*c - b*d
        
        float32x4_t imag_result = vmulq_f32(in_complex.val[0], filt_complex.val[1]); // a*d  
        imag_result = vmlaq_f32(imag_result, in_complex.val[1], filt_complex.val[0]); // a*d + b*c
        
        // Apply gain
        if (gain != 1.0f) {
            real_result = vmulq_n_f32(real_result, gain);
            imag_result = vmulq_n_f32(imag_result, gain);
        }
        
        // Accumulate to output
        out_complex.val[0] = vaddq_f32(out_complex.val[0], real_result);
        out_complex.val[1] = vaddq_f32(out_complex.val[1], imag_result);
        
        // Store result
        vst2q_f32(out_ptr, out_complex);
        
        in_ptr += 8; filt_ptr += 8; out_ptr += 8;
    }
    
    // Handle remaining samples (scalar fallback)
    for (int i = simd_samples; i < complex_samples; i++) {
        float in_real = in_ptr[0], in_imag = in_ptr[1];
        float filt_real = filt_ptr[0], filt_imag = filt_ptr[1];
        
        float result_real = (in_real * filt_real - in_imag * filt_imag) * gain;
        float result_imag = (in_real * filt_imag + in_imag * filt_real) * gain;
        
        out_ptr[0] += result_real;
        out_ptr[1] += result_imag;
        
        in_ptr += 2; filt_ptr += 2; out_ptr += 2;
    }
    
#elif defined(__SSE2__) || defined(__AVX__)
    // x86 SSE/AVX optimized version
    const float* in_ptr = input;
    const float* filt_ptr = filter;
    float* out_ptr = output;
    
    // Process 2 complex samples at a time with SSE (4 floats)
    int simd_samples = complex_samples & ~1; // Round down to multiple of 2
    
    for (int i = 0; i < simd_samples; i += 2) {
        __m128 in_vec = _mm_load_ps(in_ptr);     // [in_r0, in_i0, in_r1, in_i1]
        __m128 filt_vec = _mm_load_ps(filt_ptr); // [filt_r0, filt_i0, filt_r1, filt_i1]
        __m128 out_vec = _mm_load_ps(out_ptr);
        
        // Rearrange for complex multiply
        __m128 in_real = _mm_shuffle_ps(in_vec, in_vec, _MM_SHUFFLE(2,0,2,0)); // [r0,r0,r1,r1]
        __m128 in_imag = _mm_shuffle_ps(in_vec, in_vec, _MM_SHUFFLE(3,1,3,1)); // [i0,i0,i1,i1]
        __m128 filt_real = _mm_shuffle_ps(filt_vec, filt_vec, _MM_SHUFFLE(2,0,2,0));
        __m128 filt_imag = _mm_shuffle_ps(filt_vec, filt_vec, _MM_SHUFFLE(3,1,3,1));
        
        // Complex multiply
        __m128 ac = _mm_mul_ps(in_real, filt_real);
        __m128 bd = _mm_mul_ps(in_imag, filt_imag);
        __m128 ad = _mm_mul_ps(in_real, filt_imag);
        __m128 bc = _mm_mul_ps(in_imag, filt_real);
        
        __m128 result_real = _mm_sub_ps(ac, bd); // ac - bd
        __m128 result_imag = _mm_add_ps(ad, bc); // ad + bc
        
        // Interleave real and imaginary parts
        __m128 result_low = _mm_unpacklo_ps(result_real, result_imag);
        __m128 result_high = _mm_unpackhi_ps(result_real, result_imag);
        
        // Apply gain
        if (gain != 1.0f) {
            __m128 gain_vec = _mm_set1_ps(gain);
            result_low = _mm_mul_ps(result_low, gain_vec);
            result_high = _mm_mul_ps(result_high, gain_vec);
        }
        
        // Accumulate
        out_vec = _mm_add_ps(out_vec, result_low);
        _mm_store_ps(out_ptr, out_vec);
        
        in_ptr += 4; filt_ptr += 4; out_ptr += 4;
    }
    
    // Handle remaining samples
    for (int i = simd_samples; i < complex_samples; i++) {
        float in_real = in_ptr[0], in_imag = in_ptr[1];
        float filt_real = filt_ptr[0], filt_imag = filt_ptr[1];
        
        float result_real = (in_real * filt_real - in_imag * filt_imag) * gain;
        float result_imag = (in_real * filt_imag + in_imag * filt_real) * gain;
        
        out_ptr[0] += result_real;
        out_ptr[1] += result_imag;
        
        in_ptr += 2; filt_ptr += 2; out_ptr += 2;
    }
    
#else
    // Scalar fallback version with proper complex math
    for (int i = 0; i < complex_samples; i++) {
        float in_real = input[i*2], in_imag = input[i*2 + 1];
        float filt_real = filter[i*2], filt_imag = filter[i*2 + 1];
        
        // Complex multiply: (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
        float result_real = (in_real * filt_real - in_imag * filt_imag) * gain;
        float result_imag = (in_real * filt_imag + in_imag * filt_real) * gain;
        
        // Accumulate to output
        output[i*2] += result_real;
        output[i*2 + 1] += result_imag;
    }
#endif
}

void DAFilter::process_output_channels(int thread_id, int start_channel, int end_channel) {
    FFTContext& fft_ctx = fft_contexts_[thread_id];
    float* temp_buf = get_temp_buffer(thread_id);
    
    for (int ch = start_channel; ch < end_channel; ch++) {
        float* freq_buf = get_output_freq_buffer(ch);
        float* time_buf = get_output_time_buffer(ch);
        float* overlap_buf = get_overlap_buffer(ch);
        
        // Perform inverse FFT
        ippsFFTInv_PermToR_32f(freq_buf, temp_buf, fft_ctx.fft_spec, fft_ctx.fft_buffer);
        
        // Apply overlap-add
        for (int i = 0; i < overlap_size_; i++) {
            time_buf[i] = overlap_buf[i] + temp_buf[i];
        }
        
        // Copy new block data and save overlap for next iteration
        for (int i = 0; i < block_size_; i++) {
            time_buf[overlap_size_ + i] = temp_buf[overlap_size_ + i];
        }
        for (int i = 0; i < overlap_size_; i++) {
            overlap_buf[i] = temp_buf[block_size_ + overlap_size_ + i];
        }
        
        // Clear frequency buffer for next iteration
        std::memset(freq_buf, 0, fft_size_ * 2 * sizeof(float));
        
        // Update peak level
        float peak;
        ippsMaxAbs_32f(time_buf + overlap_size_, block_size_, &peak);
        stats_.peak_output[ch] = 0.99f * stats_.peak_output[ch] + 0.01f * peak;
    }
}

void DAFilter::initialize_fft_contexts() {
    fft_contexts_.resize(num_threads_);
    
    int fft_order = (int)(std::log2(fft_size_) + 0.5);
    
    for (int i = 0; i < num_threads_; i++) {
        FFTContext& ctx = fft_contexts_[i];
        
        int spec_size, init_size, buf_size;
        ippsFFTGetSize_R_32f(fft_order, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, 
                            &spec_size, &init_size, &buf_size);
        
        ctx.fft_spec_buffer = ippsMalloc_8u(spec_size);
        ctx.fft_init_buffer = ippsMalloc_8u(init_size);
        ctx.fft_buffer = ippsMalloc_8u(buf_size);
        
        ippsFFTInit_R_32f(&ctx.fft_spec, fft_order, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, 
                          ctx.fft_spec_buffer, ctx.fft_init_buffer);
    }
}

void DAFilter::cleanup_fft_contexts() {
    for (FFTContext& ctx : fft_contexts_) {
        if (ctx.fft_buffer) ippsFree(ctx.fft_buffer);
        if (ctx.fft_init_buffer) ippsFree(ctx.fft_init_buffer);
        if (ctx.fft_spec_buffer) ippsFree(ctx.fft_spec_buffer);
    }
    fft_contexts_.clear();
}

void DAFilter::generate_window_function() {
    for (int i = 0; i < fft_size_; i++) {
        switch (window_type_) {
            case WindowType::HANN:
                window_function_[i] = 0.5f * (1.0f - std::cos(2.0f * M_PI * i / (fft_size_ - 1)));
                break;
            case WindowType::HAMMING:
                window_function_[i] = 0.54f - 0.46f * std::cos(2.0f * M_PI * i / (fft_size_ - 1));
                break;
            case WindowType::BLACKMAN:
                window_function_[i] = 0.42f - 0.5f * std::cos(2.0f * M_PI * i / (fft_size_ - 1)) + 
                                     0.08f * std::cos(4.0f * M_PI * i / (fft_size_ - 1));
                break;
            case WindowType::KAISER:
            default:
                window_function_[i] = 1.0f; // Rectangular window as fallback
                break;
        }
    }
}

int DAFilter::calculate_optimal_fft_size(int block_size, int overlap_factor) {
    int min_fft_size = block_size * overlap_factor;
    
    // Find next power of 2
    int fft_size = 1;
    while (fft_size < min_fft_size) {
        fft_size <<= 1;
    }
    
    return fft_size;
}

// Filter bank management functions
int DAFilter::create_filter_bank(const std::string& name) {
    FilterBank bank;
    bank.name = name;
    bank.active = false;
    bank.mix_level = 1.0f;
    bank.crossfade_progress = 0.0f;
    
    filter_banks_.push_back(bank);
    return filter_banks_.size() - 1;
}

bool DAFilter::load_filter_bank(int bank_index, const std::vector<FilterConfig>& filters) {
    if (bank_index < 0 || bank_index >= filter_banks_.size()) {
        return false;
    }
    
    filter_banks_[bank_index].filters = filters;
    return true;
}

bool DAFilter::activate_filter_bank(int bank_index, FilterBankMode mode) {
    if (bank_index < 0 || bank_index >= filter_banks_.size()) {
        return false;
    }
    
    bank_mode_ = mode;
    
    switch (mode) {
        case FilterBankMode::REPLACE:
            active_bank_index_ = bank_index;
            target_bank_index_ = bank_index;
            break;
            
        case FilterBankMode::CROSSFADE:
            target_bank_index_ = bank_index;
            // Initialize crossfade timing
            crossfade_start_block_ = stats_.blocks_processed;
            // Convert duration from ms to blocks (assuming 48kHz sample rate)
            float sample_rate = 48000.0f; // TODO: Make this configurable
            float blocks_per_second = sample_rate / block_size_;
            crossfade_duration_blocks_ = (crossfade_duration_ms_ / 1000.0f) * blocks_per_second;
            break;
            
        case FilterBankMode::PARALLEL:
            // Enable multiple banks - implementation would need modification
            active_bank_index_ = bank_index;
            break;
    }
    
    return true;
}

// Utility functions
std::vector<float> DAFilter::design_window(WindowType type, int size, float beta) {
    std::vector<float> window(size);
    
    for (int i = 0; i < size; i++) {
        switch (type) {
            case WindowType::HANN:
                window[i] = 0.5f * (1.0f - std::cos(2.0f * M_PI * i / (size - 1)));
                break;
            case WindowType::HAMMING:
                window[i] = 0.54f - 0.46f * std::cos(2.0f * M_PI * i / (size - 1));
                break;
            case WindowType::BLACKMAN:
                window[i] = 0.42f - 0.5f * std::cos(2.0f * M_PI * i / (size - 1)) + 
                           0.08f * std::cos(4.0f * M_PI * i / (size - 1));
                break;
            case WindowType::KAISER:
            default:
                window[i] = 1.0f; // Rectangular window as fallback
                break;
        }
    }
    
    return window;
}

float DAFilter::calculate_filter_delay(int filter_length, int block_size) {
    return static_cast<float>(filter_length) / 2.0f;
}

// Additional implementation stubs for remaining methods...
void DAFilter::process_block(const int32_t* const* inputs, int32_t* const* outputs, int scale_factor) {
    // Implementation for integer processing - would convert to float, process, convert back
}

bool DAFilter::add_filter_to_bank(int bank_index, const FilterConfig& filter) {
    if (bank_index < 0 || bank_index >= filter_banks_.size()) {
        return false;
    }
    filter_banks_[bank_index].filters.push_back(filter);
    return true;
}

// Additional stub implementations for completeness...
bool DAFilter::update_filter_in_bank(int bank_index, int filter_index, const FilterConfig& filter) { return false; }
bool DAFilter::remove_filter_from_bank(int bank_index, int filter_index) { return false; }
bool DAFilter::remove_filter_bank(int bank_index) { return false; }
void DAFilter::clear_all_filter_banks() { filter_banks_.clear(); }
void DAFilter::set_bank_mix_level(int bank_index, float level) {}
void DAFilter::set_input_gain(int channel, float gain) {}
void DAFilter::set_output_gain(int channel, float gain) {}
void DAFilter::update_crossfade_progress() {
    if (bank_mode_ != FilterBankMode::CROSSFADE || 
        active_bank_index_ == target_bank_index_) {
        return;
    }
    
    uint64_t current_block = stats_.blocks_processed;
    uint64_t elapsed_blocks = current_block - crossfade_start_block_;
    float duration_blocks = crossfade_duration_blocks_;
    
    if (duration_blocks <= 0.0f) {
        // Instant switch if no duration specified
        active_bank_index_ = target_bank_index_;
        return;
    }
    
    float progress = static_cast<float>(elapsed_blocks) / duration_blocks;
    
    if (progress >= 1.0f) {
        // Crossfade complete - switch to target bank
        active_bank_index_ = target_bank_index_;
        if (active_bank_index_ < filter_banks_.size()) {
            filter_banks_[active_bank_index_].crossfade_progress = 1.0f;
        }
    } else {
        // Update progress for both banks
        if (active_bank_index_ < filter_banks_.size()) {
            filter_banks_[active_bank_index_].crossfade_progress = 1.0f - progress;
        }
        if (target_bank_index_ < filter_banks_.size()) {
            filter_banks_[target_bank_index_].crossfade_progress = progress;
        }
    }
}

float DAFilter::get_crossfade_progress() const {
    if (bank_mode_ != FilterBankMode::CROSSFADE || 
        active_bank_index_ == target_bank_index_) {
        return 1.0f; // No crossfade active
    }
    
    uint64_t current_block = stats_.blocks_processed;
    uint64_t elapsed_blocks = current_block - crossfade_start_block_;
    float duration_blocks = crossfade_duration_blocks_;
    
    if (duration_blocks <= 0.0f) {
        return 1.0f;
    }
    
    float progress = static_cast<float>(elapsed_blocks) / duration_blocks;
    return std::min(1.0f, std::max(0.0f, progress));
}

float* DAFilter::get_filter_freq_response(int filter_index, const FilterBank& bank) {
    // In a real implementation, this would return pre-computed frequency domain 
    // filter coefficients. For now, return a simple all-pass filter (gain only)
    static thread_local std::vector<float> temp_filter_freq;
    
    if (temp_filter_freq.size() != fft_size_ * 2) {
        temp_filter_freq.resize(fft_size_ * 2);
    }
    
    // Create a simple gain filter in frequency domain
    if (filter_index < bank.filters.size()) {
        float gain = bank.filters[filter_index].gain;
        for (int i = 0; i < fft_size_ * 2; i += 2) {
            temp_filter_freq[i] = gain;     // Real part
            temp_filter_freq[i + 1] = 0.0f; // Imaginary part (pure gain has no phase)
        }
    } else {
        // Unity filter
        for (int i = 0; i < fft_size_ * 2; i += 2) {
            temp_filter_freq[i] = 1.0f;     // Real part
            temp_filter_freq[i + 1] = 0.0f; // Imaginary part
        }
    }
    
    return temp_filter_freq.data();
}
void DAFilter::enable_parallel_bank_processing(bool enable) {}
void DAFilter::set_thread_affinity(const std::vector<int>& cpu_cores) {}
void DAFilter::export_filter_bank(int bank_index, const std::string& filename) const {}
bool DAFilter::import_filter_bank(const std::string& filename) { return false; }