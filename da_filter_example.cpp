#include "da_filter.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>

/**
 * Example usage of DAFilter - demonstrates multi-channel processing
 * with filter banks and crossfading capabilities
 */

void generate_test_signal(float* buffer, int samples, float frequency, float sample_rate) {
    for (int i = 0; i < samples; i++) {
        buffer[i] = 0.5f * std::sin(2.0f * M_PI * frequency * i / sample_rate);
    }
}

void add_noise(float* buffer, int samples, float level) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    
    for (int i = 0; i < samples; i++) {
        buffer[i] += level * dist(gen);
    }
}

int main() {
    std::cout << "DAFilter Example - Multi-channel convolution with filter banks\n";
    std::cout << "============================================================\n\n";
    
    // Configuration
    const int BLOCK_SIZE = 512;
    const int NUM_INPUT_CHANNELS = 8;
    const int NUM_OUTPUT_CHANNELS = 8;
    const int NUM_THREADS = 4;
    const int OVERLAP_FACTOR = 4;
    const float SAMPLE_RATE = 48000.0f;
    const int NUM_BLOCKS_TO_PROCESS = 100;
    
    // Create and initialize the filter
    DAFilter filter;
    
    std::cout << "Initializing DAFilter...\n";
    std::cout << "  Block size: " << BLOCK_SIZE << " samples\n";
    std::cout << "  Input channels: " << NUM_INPUT_CHANNELS << "\n";
    std::cout << "  Output channels: " << NUM_OUTPUT_CHANNELS << "\n";
    std::cout << "  Threads: " << NUM_THREADS << "\n";
    std::cout << "  Overlap factor: " << OVERLAP_FACTOR << "\n";
    
    if (!filter.initialize(BLOCK_SIZE, NUM_INPUT_CHANNELS, NUM_OUTPUT_CHANNELS, 
                          NUM_THREADS, OVERLAP_FACTOR, DAFilter::WindowType::HANN)) {
        std::cerr << "Failed to initialize DAFilter!\n";
        return -1;
    }
    
    std::cout << "  FFT size: " << DAFilter::calculate_optimal_fft_size(BLOCK_SIZE, OVERLAP_FACTOR) << "\n";
    std::cout << "  Latency: " << filter.get_latency_samples() << " samples (" 
              << (filter.get_latency_samples() * 1000.0f / SAMPLE_RATE) << " ms)\n\n";
    
    // Create filter banks
    std::cout << "Creating filter banks...\n";
    
    // Filter Bank 1: Simple gain filters
    int bank1 = filter.create_filter_bank("Gain Bank");
    std::vector<DAFilter::FilterConfig> gain_filters;
    
    for (int i = 0; i < NUM_INPUT_CHANNELS; i++) {
        DAFilter::FilterConfig config;
        config.input_channel = i;
        config.output_channel = i;
        config.gain = 0.5f; // -6dB gain
        config.enabled = true;
        gain_filters.push_back(config);
    }
    filter.load_filter_bank(bank1, gain_filters);
    std::cout << "  Bank 1: " << gain_filters.size() << " gain filters (-6dB)\n";
    
    // Filter Bank 2: Cross-channel routing
    int bank2 = filter.create_filter_bank("Cross Route Bank");
    std::vector<DAFilter::FilterConfig> route_filters;
    
    for (int i = 0; i < NUM_INPUT_CHANNELS; i++) {
        DAFilter::FilterConfig config;
        config.input_channel = i;
        config.output_channel = (i + 1) % NUM_OUTPUT_CHANNELS; // Route to next channel
        config.gain = 0.7f; // -3dB gain
        config.enabled = true;
        route_filters.push_back(config);
    }
    filter.load_filter_bank(bank2, route_filters);
    std::cout << "  Bank 2: " << route_filters.size() << " cross-routing filters (-3dB)\n\n";
    
    // Activate first filter bank
    filter.activate_filter_bank(bank1, DAFilter::FilterBankMode::REPLACE);
    
    // Allocate audio buffers
    std::vector<std::vector<float>> input_buffers(NUM_INPUT_CHANNELS, std::vector<float>(BLOCK_SIZE));
    std::vector<std::vector<float>> output_buffers(NUM_OUTPUT_CHANNELS, std::vector<float>(BLOCK_SIZE));
    
    // Create pointer arrays for processing
    std::vector<const float*> input_ptrs(NUM_INPUT_CHANNELS);
    std::vector<float*> output_ptrs(NUM_OUTPUT_CHANNELS);
    
    for (int ch = 0; ch < NUM_INPUT_CHANNELS; ch++) {
        input_ptrs[ch] = input_buffers[ch].data();
    }
    for (int ch = 0; ch < NUM_OUTPUT_CHANNELS; ch++) {
        output_ptrs[ch] = output_buffers[ch].data();
    }
    
    std::cout << "Processing " << NUM_BLOCKS_TO_PROCESS << " blocks...\n";
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    for (int block = 0; block < NUM_BLOCKS_TO_PROCESS; block++) {
        // Generate test signals for each input channel
        for (int ch = 0; ch < NUM_INPUT_CHANNELS; ch++) {
            float frequency = 440.0f + ch * 55.0f; // A4 + harmonics
            generate_test_signal(input_buffers[ch].data(), BLOCK_SIZE, frequency, SAMPLE_RATE);
            add_noise(input_buffers[ch].data(), BLOCK_SIZE, 0.01f); // Add small amount of noise
        }
        
        // Switch filter banks halfway through for demonstration
        if (block == NUM_BLOCKS_TO_PROCESS / 2) {
            std::cout << "  Switching to filter bank 2 (crossfade)...\n";
            filter.activate_filter_bank(bank2, DAFilter::FilterBankMode::CROSSFADE);
        }
        
        // Process the block
        filter.process_block(input_ptrs.data(), output_ptrs.data());
        
        // Print statistics every 20 blocks
        if (block % 20 == 0) {
            const auto& stats = filter.get_stats();
            std::cout << "  Block " << block << ": CPU Load = " 
                      << std::fixed << std::setprecision(1) << stats.cpu_load << "%";
            
            // Show peak levels for first few channels
            std::cout << ", Peaks: ";
            for (int ch = 0; ch < std::min(4, NUM_INPUT_CHANNELS); ch++) {
                std::cout << "I" << ch << "=" << std::fixed << std::setprecision(3) 
                          << stats.peak_input[ch] << " ";
            }
            for (int ch = 0; ch < std::min(4, NUM_OUTPUT_CHANNELS); ch++) {
                std::cout << "O" << ch << "=" << std::fixed << std::setprecision(3) 
                          << stats.peak_output[ch] << " ";
            }
            std::cout << "\n";
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    std::cout << "\nProcessing completed!\n";
    std::cout << "  Total time: " << duration.count() << " ms\n";
    std::cout << "  Real-time factor: " << std::fixed << std::setprecision(2) 
              << (NUM_BLOCKS_TO_PROCESS * BLOCK_SIZE * 1000.0f / SAMPLE_RATE) / duration.count() << "x\n";
    
    const auto& final_stats = filter.get_stats();
    std::cout << "  Final CPU load: " << std::fixed << std::setprecision(1) 
              << final_stats.cpu_load << "%\n";
    std::cout << "  Blocks processed: " << final_stats.blocks_processed << "\n";
    
    // Demonstrate advanced features
    std::cout << "\nAdvanced Features Demo:\n";
    std::cout << "  Active filter bank: " << filter.get_active_bank() << "\n";
    std::cout << "  Target filter bank: " << filter.get_target_bank() << "\n";
    std::cout << "  Crossfade progress: " << std::fixed << std::setprecision(1) 
              << filter.get_crossfade_progress() * 100.0f << "%\n";
    
    // Test window function generation
    auto hann_window = DAFilter::design_window(DAFilter::WindowType::HANN, 512);
    auto blackman_window = DAFilter::design_window(DAFilter::WindowType::BLACKMAN, 512);
    
    std::cout << "\nWindow Functions:\n";
    std::cout << "  Hann window (512 points): ";
    for (int i = 0; i < 5; i++) {
        std::cout << std::fixed << std::setprecision(3) << hann_window[i * 100] << " ";
    }
    std::cout << "...\n";
    
    std::cout << "  Blackman window (512 points): ";
    for (int i = 0; i < 5; i++) {
        std::cout << std::fixed << std::setprecision(3) << blackman_window[i * 100] << " ";
    }
    std::cout << "...\n";
    
    // Performance comparison
    std::cout << "\nPerformance Analysis:\n";
    std::cout << "  Theoretical latency for " << BLOCK_SIZE << " samples: " 
              << (BLOCK_SIZE * 1000.0f / SAMPLE_RATE) << " ms\n";
    std::cout << "  Actual latency with overlap: " 
              << (filter.get_latency_samples() * 1000.0f / SAMPLE_RATE) << " ms\n";
    std::cout << "  Memory usage estimate: " 
              << ((NUM_INPUT_CHANNELS + NUM_OUTPUT_CHANNELS) * 
                  DAFilter::calculate_optimal_fft_size(BLOCK_SIZE, OVERLAP_FACTOR) * 
                  sizeof(float) / 1024) << " KB for audio buffers\n";
    
    // Shutdown
    std::cout << "\nShutting down...\n";
    filter.shutdown();
    
    std::cout << "Example completed successfully!\n";
    return 0;
}