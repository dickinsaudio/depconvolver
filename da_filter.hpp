#pragma once

#include <cstdint>
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <memory>
#include <vector>
#include <cstring>
#include <complex>
#include <cassert>

#ifdef _WIN32 
#include "../extern/ipp/include/ipp.h"
#else
#if defined(__arm__) || defined(__aarch64__)
#include <arm_neon.h>
#include "../Ne10/inc/NE10.h"
// Include the same IPP compatibility layer from ThreadedDSP
#define ippStsNoErr          0
#define IPP_FFT_DIV_INV_BY_N 2
#define ippAlgHintNone       0
#define ippRndZero           0

typedef float    Ipp32f, Ipp32fc;
typedef int32_t  Ipp32s;
typedef uint8_t  Ipp8u;
typedef int      IppStatus;
typedef void     IppsFFTSpec_R_32f;
typedef int      IppRoundMode;
typedef int      IppHintAlgorithm;

// IPP compatibility functions (same as ThreadedDSP)
inline Ipp8u*    ippsMalloc_8u         (int len) { if (len==0) return 0; else return (Ipp8u*)calloc(len,1); };
inline void      ippsFree              (void *p) { if (p) free(p); };
inline IppStatus ippsSet_32f           (Ipp32f val, Ipp32f* pDst, int len) { memset (pDst, 0, len*sizeof(Ipp32f)); return ippStsNoErr; };
inline IppStatus ippsMove_32f          (const Ipp32f* pSrc, Ipp32f* pDst, int len) { memmove(pDst, pSrc, len*sizeof(Ipp32f)); return ippStsNoErr; };
inline IppStatus ippsAdd_32f_I         (const Ipp32f* pSrc, Ipp32f* pSrcDst, int len);
inline IppStatus ippsMaxAbs_32f        (const Ipp32f* pSrc, int len, Ipp32f* pMaxAbs);
inline IppStatus ippsAddProduct_32fc   (const Ipp32fc* pSrc1, const Ipp32fc* pSrc2, Ipp32fc* pSrcDst, int len);
inline IppStatus ippsFFTGetSize_R_32f  (int order, int flag, IppHintAlgorithm hint, int*pSpecSize, int* pSpecBufferSize, int* pBufferSize);
inline IppStatus ippsFFTInit_R_32f     (IppsFFTSpec_R_32f** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer);
inline IppStatus ippSetDenormAreZeros  (int value) { return ippStsNoErr; };
inline IppStatus ippsFFTFwd_RToPerm_32f(Ipp32f* pSrc, Ipp32f* pDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer);
inline IppStatus ippsFFTInv_PermToR_32f(const Ipp32f* pSrc, Ipp32f* pDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer);
inline IppStatus ippsConvert_32f32s_Sfs(const Ipp32f* pSrc, Ipp32s* pDst, int len, IppRoundMode rndMode, int scaleFactor);
inline IppStatus ippsMulC_32f_I        (Ipp32f val, Ipp32f* pSrcDst, int len);
inline IppStatus ippsMul_32f_I         (const Ipp32f* pSrc, Ipp32f* pSrcDst, int len);
#else 
#include "ipp.h"
#endif
#endif

/**
 * DAFilter - Advanced multi-threaded convolution filter processor
 * 
 * Key improvements over ThreadedDSP:
 * - Enhanced filter bank management with smooth crossfading
 * - Optimized overlap-add processing for reduced latency
 * - Better thread load balancing
 * - Support for dynamic filter bank switching
 * - Improved memory layout for cache efficiency
 */
class DAFilter
{
public:
    static constexpr int MAX_CHANNELS = 1024;    // Increased from 512
    static constexpr int MAX_BLOCK_SIZE = 2048;  // Increased from 1024  
    static constexpr int MAX_FILTER_LENGTH = 524288; // Increased from 262144
    static constexpr int MAX_THREADS = 64;       // Increased from 32
    static constexpr int MAX_FILTER_BANKS = 16;  // New: support multiple filter banks
    static constexpr int MAX_FILTERS_PER_BANK = 512;
    static constexpr int DEFAULT_OVERLAP_FACTOR = 4; // For crossfading

    enum class FilterBankMode {
        REPLACE,        // Immediate replacement
        CROSSFADE,      // Smooth crossfade between banks
        PARALLEL        // Run banks in parallel and mix
    };

    enum class WindowType {
        HANN,
        HAMMING, 
        BLACKMAN,
        KAISER
    };

    struct FilterConfig {
        int input_channel = -1;
        int output_channel = -1;
        int length = 0;
        float gain = 1.0f;
        bool enabled = true;
        std::vector<float> coefficients;
    };

    struct FilterBank {
        std::string name;
        std::vector<FilterConfig> filters;
        float mix_level = 1.0f;
        bool active = false;
        float crossfade_progress = 0.0f; // 0.0 = fade out, 1.0 = fade in
    };

    struct ProcessingStats {
        float cpu_load = 0.0f;
        float peak_input[MAX_CHANNELS] = {0};
        float peak_output[MAX_CHANNELS] = {0};
        uint64_t blocks_processed = 0;
        std::chrono::high_resolution_clock::time_point last_process_time;
    };

private:
    // Core configuration
    int block_size_ = 0;
    int num_input_channels_ = 0;
    int num_output_channels_ = 0;
    int fft_size_ = 0;
    int overlap_size_ = 0;
    int num_threads_ = 0;
    WindowType window_type_ = WindowType::HANN;
    
    // Threading infrastructure
    std::vector<std::unique_ptr<std::thread>> worker_threads_;
    std::mutex process_mutex_;
    std::condition_variable start_processing_;
    std::condition_variable processing_complete_;
    std::atomic<bool> running_{false};
    std::atomic<bool> processing_{false};
    std::atomic<int> threads_ready_{0};
    std::atomic<int> threads_completed_{0};
    
    // Filter banks and crossfading
    std::vector<FilterBank> filter_banks_;
    std::atomic<int> active_bank_index_{-1};
    std::atomic<int> target_bank_index_{-1};
    std::atomic<float> crossfade_duration_ms_{100.0f};
    FilterBankMode bank_mode_ = FilterBankMode::CROSSFADE;
    
    // Audio buffers (improved memory layout)
    std::unique_ptr<float[]> input_time_buffer_;     // [channels][overlap_history + block_size]
    std::unique_ptr<float[]> output_time_buffer_;    // [channels][overlap_history + block_size] 
    std::unique_ptr<float[]> input_freq_buffer_;     // [channels][fft_size * 2] (complex)
    std::unique_ptr<float[]> output_freq_buffer_;    // [channels][fft_size * 2] (complex)
    std::unique_ptr<float[]> filter_freq_buffer_;    // [bank][filter][block][fft_size * 2]
    std::unique_ptr<float[]> overlap_buffer_;        // [channels][overlap_size]
    std::unique_ptr<float[]> window_function_;       // [fft_size]
    std::unique_ptr<float[]> temp_buffers_;          // [threads][fft_size * 2]
    
    // FFT infrastructure (per-thread)
    struct FFTContext {
        IppsFFTSpec_R_32f* fft_spec = nullptr;
        Ipp8u* fft_buffer = nullptr;
        Ipp8u* fft_spec_buffer = nullptr;
        Ipp8u* fft_init_buffer = nullptr;
    };
    std::vector<FFTContext> fft_contexts_;
    
    // Statistics and monitoring
    ProcessingStats stats_;
    std::chrono::high_resolution_clock::time_point process_start_time_;
    
    // Crossfade state tracking
    std::atomic<uint64_t> crossfade_start_block_{0};
    std::atomic<float> crossfade_duration_blocks_{0.0f};
    
    // Internal helper functions
    void initialize_fft_contexts();
    void cleanup_fft_contexts(); 
    void generate_window_function();
    void worker_thread_proc(int thread_id);
    void process_input_channels(int thread_id, int start_channel, int end_channel);
    void process_filters(int thread_id, int start_filter, int end_filter);
    void process_bank_filters(int thread_id, int start_filter, int end_filter, 
                             const FilterBank& bank, float bank_gain);
    void complex_multiply_accumulate_optimized(const float* input, const float* filter, 
                                             float* output, int complex_samples, float gain);
    float* get_filter_freq_response(int filter_index, const FilterBank& bank);
    void process_output_channels(int thread_id, int start_channel, int end_channel);
    void update_crossfade_progress();
    void apply_overlap_add(int channel);
    
    // Buffer access helpers (cache-friendly indexing)
    float* get_input_time_buffer(int channel) {
        return input_time_buffer_.get() + channel * (overlap_size_ + block_size_);
    }
    
    float* get_output_time_buffer(int channel) {
        return output_time_buffer_.get() + channel * (overlap_size_ + block_size_);
    }
    
    float* get_input_freq_buffer(int channel) {
        return input_freq_buffer_.get() + channel * fft_size_ * 2;
    }
    
    float* get_output_freq_buffer(int channel) {
        return output_freq_buffer_.get() + channel * fft_size_ * 2;
    }
    
    float* get_overlap_buffer(int channel) {
        return overlap_buffer_.get() + channel * overlap_size_;
    }
    
    float* get_temp_buffer(int thread_id) {
        return temp_buffers_.get() + thread_id * fft_size_ * 2;
    }

public:
    DAFilter();
    ~DAFilter();
    
    // Core lifecycle
    bool initialize(int block_size, int num_input_channels, int num_output_channels, 
                   int num_threads = 0, int overlap_factor = DEFAULT_OVERLAP_FACTOR,
                   WindowType window_type = WindowType::HANN);
    void shutdown();
    
    // Audio processing
    void process_block(const float* const* inputs, float* const* outputs);
    void process_block(const int32_t* const* inputs, int32_t* const* outputs, int scale_factor = 8);
    
    // Filter bank management  
    int create_filter_bank(const std::string& name);
    bool load_filter_bank(int bank_index, const std::vector<FilterConfig>& filters);
    bool activate_filter_bank(int bank_index, FilterBankMode mode = FilterBankMode::CROSSFADE);
    bool remove_filter_bank(int bank_index);
    void clear_all_filter_banks();
    
    // Individual filter management within banks
    bool add_filter_to_bank(int bank_index, const FilterConfig& filter);
    bool update_filter_in_bank(int bank_index, int filter_index, const FilterConfig& filter);
    bool remove_filter_from_bank(int bank_index, int filter_index);
    
    // Configuration and control
    void set_crossfade_duration(float duration_ms) { crossfade_duration_ms_ = duration_ms; }
    void set_bank_mix_level(int bank_index, float level);
    void set_input_gain(int channel, float gain);
    void set_output_gain(int channel, float gain);
    
    // Status and monitoring
    const ProcessingStats& get_stats() const { return stats_; }
    bool is_running() const { return running_; }
    int get_active_bank() const { return active_bank_index_; }
    int get_target_bank() const { return target_bank_index_; }
    float get_crossfade_progress() const;
    int get_latency_samples() const { return overlap_size_; }
    
    // Advanced features
    void enable_parallel_bank_processing(bool enable);
    void set_thread_affinity(const std::vector<int>& cpu_cores);
    void export_filter_bank(int bank_index, const std::string& filename) const;
    bool import_filter_bank(const std::string& filename);
    
    // Utility functions
    static int calculate_optimal_fft_size(int block_size, int overlap_factor = DEFAULT_OVERLAP_FACTOR);
    static float calculate_filter_delay(int filter_length, int block_size);
    static std::vector<float> design_window(WindowType type, int size, float beta = 0.0f);
};
