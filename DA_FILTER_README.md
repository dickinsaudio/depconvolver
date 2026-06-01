# DAFilter - Advanced Multi-threaded Convolution Filter

DAFilter is an enhanced implementation inspired by ThreadedDSP, designed for high-performance multi-channel audio convolution with advanced features including filter banks, crossfading, and optimized parallel processing.

## Key Features

### Enhanced over ThreadedDSP
- **Filter Banks**: Support for multiple filter sets with smooth switching
- **Crossfading**: Seamless transitions between filter banks
- **Improved Threading**: Better load balancing and reduced latency
- **Cache Optimization**: Memory layout optimized for modern CPUs
- **Larger Capacity**: Support for up to 1024 channels, 2048 sample blocks
- **Advanced Windowing**: Multiple window functions (Hann, Hamming, Blackman, Kaiser)

### Core Capabilities
- **Multi-threaded Processing**: Parallel FFT, multiply-accumulate, and iFFT
- **Overlap-Add Processing**: Efficient convolution with configurable overlap
- **Real-time Performance**: Optimized for low-latency audio processing
- **Cross-platform**: ARM NEON and Intel IPP acceleration
- **Memory Efficient**: Shared buffers and optimized data structures

## Architecture

### Threading Model
```
Input Channels → [Thread Pool] → FFT Forward
                      ↓
Filter Processing → [Thread Pool] → Frequency Domain MAC
                      ↓  
Output Channels → [Thread Pool] → FFT Inverse + Overlap-Add
```

### Filter Bank Management
- **Multiple Banks**: Up to 16 filter banks with 512 filters each
- **Switching Modes**: 
  - `REPLACE`: Immediate switching
  - `CROSSFADE`: Smooth transition with configurable duration
  - `PARALLEL`: Multiple banks active simultaneously
- **Dynamic Loading**: Runtime filter coefficient updates

### Memory Layout
```
Input Time Buffer:  [CH0][CH1]...[CHn] → [Overlap + Block Size]
Input Freq Buffer:  [CH0][CH1]...[CHn] → [FFT Size * 2 (complex)]
Filter Freq Data:   [Bank][Filter][Block] → [FFT Size * 2]
Output Freq Buffer: [CH0][CH1]...[CHn] → [FFT Size * 2 (complex)]
Output Time Buffer: [CH0][CH1]...[CHn] → [Overlap + Block Size]
```

## Usage Example

```cpp
#include "da_filter.hpp"

// Initialize filter
DAFilter filter;
filter.initialize(
    512,    // block_size
    8,      // input_channels  
    8,      // output_channels
    4,      // num_threads
    4,      // overlap_factor
    DAFilter::WindowType::HANN
);

// Create filter bank
int bank_id = filter.create_filter_bank("My Filters");

// Add filters to bank
std::vector<DAFilter::FilterConfig> filters;
for (int i = 0; i < 8; i++) {
    DAFilter::FilterConfig config;
    config.input_channel = i;
    config.output_channel = i;
    config.gain = 0.7f;
    config.enabled = true;
    filters.push_back(config);
}

filter.load_filter_bank(bank_id, filters);
filter.activate_filter_bank(bank_id);

// Process audio blocks
const float* inputs[8];   // Input channel pointers
float* outputs[8];        // Output channel pointers

filter.process_block(inputs, outputs);

// Monitor performance
const auto& stats = filter.get_stats();
std::cout << "CPU Load: " << stats.cpu_load << "%\n";
```

## Building

### Prerequisites
- C++17 compatible compiler
- CMake 3.8+
- pthread (Linux/macOS)
- Optional: Intel IPP (x86_64)
- Optional: ARM Ne10 (ARM/AARCH64)

### Build Steps
```bash
# Add to existing CMakeLists.txt or build standalone
mkdir build && cd build
cmake ..
make da_filter_example

# Run example
./da_filter_example
```

### Integration with DepConvolver
```cmake
# In main CMakeLists.txt, add:
include(da_filter_CMakeLists.txt)
```

## Performance Characteristics

### Latency
- **Algorithmic Latency**: `(FFT_Size - Block_Size) samples`
- **Typical Values**: 
  - 512 sample blocks → 1536 samples latency (32ms @ 48kHz)
  - 256 sample blocks → 768 samples latency (16ms @ 48kHz)

### Throughput
- **Multi-core Scaling**: Near-linear with thread count
- **Memory Bandwidth**: Optimized for cache hierarchy
- **SIMD Optimization**: ARM NEON and x86 SSE/AVX

### Resource Usage
- **Memory**: ~4KB per channel per block for buffers
- **CPU**: Scales with number of active filters and FFT size
- **Real-time Performance**: Typically <50% CPU for moderate loads

## Advanced Features

### Filter Bank Crossfading
```cpp
// Set crossfade duration
filter.set_crossfade_duration(100.0f); // 100ms

// Switch with crossfade
filter.activate_filter_bank(new_bank_id, DAFilter::FilterBankMode::CROSSFADE);

// Monitor crossfade progress
float progress = filter.get_crossfade_progress(); // 0.0 to 1.0
```

### Window Function Design
```cpp
// Generate custom windows
auto hann = DAFilter::design_window(DAFilter::WindowType::HANN, 512);
auto blackman = DAFilter::design_window(DAFilter::WindowType::BLACKMAN, 1024);
```

### Performance Monitoring
```cpp
const auto& stats = filter.get_stats();
std::cout << "Blocks processed: " << stats.blocks_processed << "\n";
std::cout << "Input peaks: ";
for (int i = 0; i < 8; i++) {
    std::cout << stats.peak_input[i] << " ";
}
```

## Comparison with ThreadedDSP

| Feature | ThreadedDSP | DAFilter |
|---------|-------------|----------|
| Max Channels | 512 | 1024 |
| Max Block Size | 1024 | 2048 |
| Filter Banks | 2 groups | 16 banks |
| Crossfading | Basic | Advanced with timing |
| Memory Layout | Basic | Cache-optimized |
| Window Functions | Fixed | Multiple types |
| Thread Scaling | Good | Improved |

## Future Enhancements

- **GPU Acceleration**: CUDA/OpenCL backends
- **Variable Block Sizes**: Dynamic block size adaptation  
- **Advanced Filters**: IIR, parametric EQ, dynamics
- **Network Distribution**: Multi-machine processing
- **Real-time Visualization**: Spectrum analysis and monitoring

## License

Compatible with DepConvolver project licensing.

## Contributing

Improvements welcome, especially:
- Additional window functions
- Platform-specific optimizations  
- Filter design utilities
- Performance benchmarks