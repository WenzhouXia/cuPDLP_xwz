#ifndef CUPDLP_CUDA_KERNALS_H
#define CUPDLP_CUDA_KERNALS_H

#include "cuda_runtime.h"
#define CUPDLP_BLOCK_SIZE 512

#ifndef SFLOAT
#ifdef DLONG
typedef long long cupdlp_int;
#else
typedef int cupdlp_int;
#endif
typedef double cupdlp_float;
#define CudaComputeType CUDA_R_64F
#else
#define CudaComputeType CUDA_R_32F
#endif

#define CHECK_CUDA(func)                                               \
  {                                                                    \
    cudaError_t status = (func);                                       \
    if (status != cudaSuccess)                                         \
    {                                                                  \
      printf("CUDA API failed at line %d of %s with error: %s (%d)\n", \
             __LINE__, __FILE__, cudaGetErrorString(status), status);  \
      return EXIT_FAILURE;                                             \
    }                                                                  \
  }

#define CHECK_CUSPARSE(func)                                               \
  {                                                                        \
    cusparseStatus_t status = (func);                                      \
    if (status != CUSPARSE_STATUS_SUCCESS)                                 \
    {                                                                      \
      printf("CUSPARSE API failed at line %d of %s with error: %s (%d)\n", \
             __LINE__, __FILE__, cusparseGetErrorString(status), status);  \
      return EXIT_FAILURE;                                                 \
    }                                                                      \
  }

#define CHECK_CUBLAS(func)                                               \
  {                                                                      \
    cublasStatus_t status = (func);                                      \
    if (status != CUBLAS_STATUS_SUCCESS)                                 \
    {                                                                    \
      printf("CUBLAS API failed at line %d of %s with error: %s (%d)\n", \
             __LINE__, __FILE__, cublasGetStatusString(status), status); \
      return EXIT_FAILURE;                                               \
    }                                                                    \
  }

#define CUPDLP_FREE_VEC(x) \
  {                        \
    cudaFree(x);           \
    x = cupdlp_NULL;       \
  }

#define CUPDLP_COPY_VEC(dst, src, type, size) \
  cudaMemcpy(dst, src, sizeof(type) * (size), cudaMemcpyDefault)

#define CUPDLP_INIT_VEC(var, size)                                             \
  {                                                                            \
    cusparseStatus_t status =                                                  \
        cudaMalloc((void **)&var, (size) * sizeof(typeof(*var)));              \
    if (status != CUSPARSE_STATUS_SUCCESS)                                     \
    {                                                                          \
      printf("CUSPARSE API failed at line %d with error: %s (%d)\n", __LINE__, \
             cusparseGetErrorString(status), status);                          \
      goto exit_cleanup;                                                       \
    }                                                                          \
  }
#define CUPDLP_INIT_ZERO_VEC(var, size)                                        \
  {                                                                            \
    cusparseStatus_t status =                                                  \
        cudaMalloc((void **)&var, (size) * sizeof(typeof(*var)));              \
    if (status != CUSPARSE_STATUS_SUCCESS)                                     \
    {                                                                          \
      printf("CUSPARSE API failed at line %d with error: %s (%d)\n", __LINE__, \
             cusparseGetErrorString(status), status);                          \
      goto exit_cleanup;                                                       \
    }                                                                          \
    status = cudaMemset(var, 0, (size) * sizeof(typeof(*var)));                \
    if (status != CUSPARSE_STATUS_SUCCESS)                                     \
    {                                                                          \
      printf("CUSPARSE API failed at line %d with error: %s (%d)\n", __LINE__, \
             cusparseGetErrorString(status), status);                          \
      goto exit_cleanup;                                                       \
    }                                                                          \
  }
#define CUPDLP_ZERO_VEC(var, type, size) \
  cudaMemset(var, 0, sizeof(type) * (size))

dim3 cuda_gridsize(cupdlp_int n);

__global__ void element_wise_dot_kernel(cupdlp_float *x, const cupdlp_float *y,
                                        const cupdlp_int len);

__global__ void element_wise_div_kernel(cupdlp_float *x, const cupdlp_float *y,
                                        const cupdlp_int len);

__global__ void element_wise_projlb_kernel(cupdlp_float *x,
                                           const cupdlp_float *lb,
                                           const cupdlp_int len);

__global__ void element_wise_projub_kernel(cupdlp_float *x,
                                           const cupdlp_float *ub,
                                           const cupdlp_int len);

__global__ void element_wise_projSamelb_kernel(cupdlp_float *x,
                                               const cupdlp_float lb,
                                               const cupdlp_int len);

__global__ void element_wise_projSameub_kernel(cupdlp_float *x,
                                               const cupdlp_float ub,
                                               const cupdlp_int len);

__global__ void element_wise_initHaslb_kernal(cupdlp_float *haslb,
                                              const cupdlp_float *lb,
                                              const cupdlp_float bound,
                                              const cupdlp_int len);

__global__ void element_wise_initHasub_kernal(cupdlp_float *hasub,
                                              const cupdlp_float *ub,
                                              const cupdlp_float bound,
                                              const cupdlp_int len);

__global__ void element_wise_filterlb_kernal(cupdlp_float *x,
                                             const cupdlp_float *lb,
                                             const cupdlp_float bound,
                                             const cupdlp_int len);

__global__ void element_wise_filterub_kernal(cupdlp_float *x,
                                             const cupdlp_float *ub,
                                             const cupdlp_float bound,
                                             const cupdlp_int len);

__global__ void init_cuda_vec_kernal(cupdlp_float *x, const cupdlp_float val,
                                     const cupdlp_int len);

__global__ void primal_grad_step_kernal(cupdlp_float *xUpdate,
                                        const cupdlp_float *x,
                                        const cupdlp_float *cost,
                                        const cupdlp_float *ATy,
                                        const cupdlp_float dPrimalStep,
                                        const cupdlp_int len);

__global__ void dual_grad_step_kernal(
    cupdlp_float *yUpdate, const cupdlp_float *y, const cupdlp_float *b,
    const cupdlp_float *Ax, const cupdlp_float *AxUpdate,
    const cupdlp_float dDualStep, const cupdlp_int len);

__global__ void dual_grad_step_kernal_AdapTheta(
    cupdlp_float *yUpdate, const cupdlp_float *y, const cupdlp_float *b,
    const cupdlp_float *Ax, const cupdlp_float *AxUpdate,
    const cupdlp_float dDualStep, const cupdlp_int len, const cupdlp_float theta);

__global__ void pdtest_x_update_with_beta_kernal(cupdlp_float *x_new,
                                                 const cupdlp_float *x_1,
                                                 const cupdlp_float *x_2,
                                                 const cupdlp_float beta,
                                                 const cupdlp_int len);

__global__ void pdtest_y_update_with_beta_kernal(cupdlp_float *y_new,
                                                 const cupdlp_float *y_1,
                                                 const cupdlp_float *y_2,
                                                 const cupdlp_float beta_inv,
                                                 const cupdlp_int len);

__global__ void pdtest_x_update_with_theta_kernal(cupdlp_float *x_barUpdate,
                                                  const cupdlp_float *xUpdate,
                                                  const cupdlp_float *x,
                                                  const cupdlp_float theta,
                                                  const cupdlp_int len);

__global__ void pdtest_primal_grad_step_kernal(cupdlp_float *xUpdate,
                                               const cupdlp_float *x,
                                               const cupdlp_float *cost,
                                               const cupdlp_float *ATy,
                                               const cupdlp_float dPrimalStep,
                                               const cupdlp_int len);

__global__ void pdtest_dual_grad_step_kernal(
    cupdlp_float *yUpdate, const cupdlp_float *y, const cupdlp_float *b,
    const cupdlp_float *Ax_bar,
    const cupdlp_float dDualStep, const cupdlp_int len);

__global__ void naive_sub_kernal(cupdlp_float *z, const cupdlp_float *x,
                                 const cupdlp_float *y, const cupdlp_int len);

__device__ double atomicAddDouble(double *address, double val);
__global__ void compute_dualOT_inf_kernal(const double *x, int vec_len, int n_coarse, double scale_const, double *d_c_norm, double *d_diff_norm);
__global__ void countZero_and_ckeckConstraint_before_cudaMalloc_kernal(long long *keep_local_len_array, const double *x, const double *y, int resolution_now, int resolution_last, double thr, double violate_degree);
__global__ void countZero_and_checkConstraint_kernal(long long *keep_fine_redundancy, long long *keep_len_UpToNow, const double *x, const double *y, int resolution_now, int resolution_last, double thr, double violate_degree);
__global__ void ckeckConstraint_before_cudaMalloc_kernal(long long *keep_local_len_array, const double *x, int resolution_now, int resolution_last, double thr, double violate_degree);
__global__ void checkConstraint_kernal(long long *keep_fine_redundancy, long long *keep_len_UpToNow, const double *x, int resolution_now, int resolution_last, double thr, double violate_degree);
#endif