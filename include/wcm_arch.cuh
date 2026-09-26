// Author(s): Ron Acda (using an iterative LLM-guided workflow, https://github.com/quarkron/iterative-hillclimber/tree/main)
// GPU-architecture portability for the fused CG / BD kernels (sm_70 and newer).
// Programmatic dependent launch (PDL) needs compute capability 9.0+ (H100, B200): there the kernels of a chain are launched
// with programmatic stream serialization and each starts with cudaGridDependencySynchronize(). On older GPUs (V100, A100, ...)
// the kernels are launched normally and the device-side synchronisation compiles to nothing (stream order already serialises).
#pragma once
#include <cstdlib>
#include <cuda_runtime.h>

__device__ __forceinline__ void wcm_grid_dependency_sync() {
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 900)
    cudaGridDependencySynchronize();
#endif
}

// warp-wide integer max over the full warp: the sm_80+ reduction instruction, a shuffle tree before sm_80 (same result)
__device__ __forceinline__ int wcm_warp_max(int v) {
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 800)
    return __reduce_max_sync(0xffffffffu, v);
#else
    for (int o = 16; o > 0; o >>= 1) v = max(v, __shfl_xor_sync(0xffffffffu, v, o));
    return v;
#endif
}

// PDL for launches on the current device: compute capability >= 9.0 and not disabled with WCM_PDL_OFF=1
static inline bool wcm_pdl_enabled() {
    static const bool on = [] {
        if (std::getenv("WCM_PDL_OFF") != nullptr) return false;
        int dev = 0, major = 0;
        if (cudaGetDevice(&dev) != cudaSuccess || cudaDeviceGetAttribute(&major, cudaDevAttrComputeCapabilityMajor, dev) != cudaSuccess) return false;
        return major >= 9;
    }();
    return on;
}
