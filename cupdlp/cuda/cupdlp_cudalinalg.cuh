#ifndef CUPDLP_CUDA_LINALG_H
#define CUPDLP_CUDA_LINALG_H

#include <cublas_v2.h>        // cublas
#include <cuda_runtime_api.h> // cudaMalloc, cudaMemcpy, etc.
#include <cusparse.h>         // cusparseSpMV

#include "cupdlp_cuda_kernels.cuh"

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>  // printf
#include <stdlib.h> // EXIT_FAILURE

// #include "../cupdlp_defs.h"
// #include "../glbopts.h"
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

extern "C" cupdlp_int cuda_alloc_MVbuffer(
    cusparseHandle_t handle, cusparseSpMatDescr_t cuda_csc,
    cusparseDnVecDescr_t vecX, cusparseDnVecDescr_t vecAx,
    cusparseSpMatDescr_t cuda_csr, cusparseDnVecDescr_t vecY,
    cusparseDnVecDescr_t vecATy, void **dBuffer);

extern "C" cupdlp_int cuda_csc_Ax(cusparseHandle_t handle,
                                  cusparseSpMatDescr_t cuda_csc,
                                  cusparseDnVecDescr_t vecX,
                                  cusparseDnVecDescr_t vecAx, void *dBuffer,
                                  const cupdlp_float alpha,
                                  const cupdlp_float beta);
extern "C" cupdlp_int cuda_csr_Ax(cusparseHandle_t handle,
                                  cusparseSpMatDescr_t cuda_csr,
                                  cusparseDnVecDescr_t vecX,
                                  cusparseDnVecDescr_t vecAx, void *dBuffer,
                                  const cupdlp_float alpha,
                                  const cupdlp_float beta);
extern "C" cupdlp_int cuda_csc_ATy(cusparseHandle_t handle,
                                   cusparseSpMatDescr_t cuda_csc,
                                   cusparseDnVecDescr_t vecY,
                                   cusparseDnVecDescr_t vecATy, void *dBuffer,
                                   const cupdlp_float alpha,
                                   const cupdlp_float beta);
extern "C" cupdlp_int cuda_csr_ATy(cusparseHandle_t handle,
                                   cusparseSpMatDescr_t cuda_csr,
                                   cusparseDnVecDescr_t vecY,
                                   cusparseDnVecDescr_t vecATy, void *dBuffer,
                                   const cupdlp_float alpha,
                                   const cupdlp_float beta);

extern "C" void cupdlp_projSameub_cuda(cupdlp_float *x, const cupdlp_float ub,
                                       const cupdlp_int len);
extern "C" void cupdlp_projSamelb_cuda(cupdlp_float *x, const cupdlp_float lb,
                                       const cupdlp_int len);
extern "C" void cupdlp_projub_cuda(cupdlp_float *x, const cupdlp_float *ub,
                                   const cupdlp_int len);
extern "C" void cupdlp_projlb_cuda(cupdlp_float *x, const cupdlp_float *lb,
                                   const cupdlp_int len);
extern "C" void cupdlp_ediv_cuda(cupdlp_float *x, const cupdlp_float *y,
                                 const cupdlp_int len);
extern "C" void cupdlp_edot_cuda(cupdlp_float *x, const cupdlp_float *y,
                                 const cupdlp_int len);
extern "C" void cupdlp_haslb_cuda(cupdlp_float *haslb, const cupdlp_float *lb,
                                  const cupdlp_float bound,
                                  const cupdlp_int len);
extern "C" void cupdlp_hasub_cuda(cupdlp_float *hasub, const cupdlp_float *ub,
                                  const cupdlp_float bound,
                                  const cupdlp_int len);
extern "C" void cupdlp_filterlb_cuda(cupdlp_float *x, const cupdlp_float *lb,
                                     const cupdlp_float bound,
                                     const cupdlp_int len);
extern "C" void cupdlp_filterub_cuda(cupdlp_float *x, const cupdlp_float *ub,
                                     const cupdlp_float bound,
                                     const cupdlp_int len);
extern "C" void cupdlp_initvec_cuda(cupdlp_float *x, const cupdlp_float val,
                                    const cupdlp_int len);

extern "C" void cupdlp_pgrad_cuda(cupdlp_float *xUpdate, const cupdlp_float *x,
                                  const cupdlp_float *cost,
                                  const cupdlp_float *ATy,
                                  const cupdlp_float dPrimalStep,
                                  const cupdlp_int len);

extern "C" void cupdlp_dgrad_cuda(cupdlp_float *yUpdate, const cupdlp_float *y,
                                  const cupdlp_float *b, const cupdlp_float *Ax,
                                  const cupdlp_float *AxUpdate,
                                  const cupdlp_float dDualStep,
                                  const cupdlp_int len);

extern "C" void cupdlp_dgrad_cuda_AdapTheta(cupdlp_float *yUpdate, const cupdlp_float *y, const cupdlp_float *b,
                                            const cupdlp_float *Ax, const cupdlp_float *AxUpdate,
                                            const cupdlp_float dDualStep, const cupdlp_int len, const cupdlp_float theta);

extern "C" void pdtest_x_md_update_cuda(cupdlp_float *x_agUpdate, const cupdlp_float *x_ag, const cupdlp_float *xUpdate, const cupdlp_float beta, const cupdlp_int len);

extern "C" void pdtest_x_ag_update_cuda(cupdlp_float *x_agUpdate, const cupdlp_float *x_ag, const cupdlp_float *xUpdate, const cupdlp_float beta, const cupdlp_int len);

extern "C" void pdtest_y_ag_update_cuda(cupdlp_float *y_agUpdate, const cupdlp_float *y_ag, const cupdlp_float *yUpdate, const cupdlp_float beta, const cupdlp_int len);

extern "C" void pdtest_x_bar_update_cuda(cupdlp_float *x_barUpdate, const cupdlp_float *xUpdate, const cupdlp_float *x, const cupdlp_float theta, const cupdlp_int len);

extern "C" void pdtest_pgrad_cuda(cupdlp_float *xUpdate, const cupdlp_float *x,
                                  const cupdlp_float *cost,
                                  const cupdlp_float *ATy,
                                  const cupdlp_float dPrimalStep,
                                  const cupdlp_int len);

extern "C" void pdtest_dgrad_cuda(cupdlp_float *yUpdate, const cupdlp_float *y, const cupdlp_float *b,
                                  const cupdlp_float *Ax_bar,
                                  const cupdlp_float dDualStep, const cupdlp_int len);

extern "C" void cupdlp_sub_cuda(cupdlp_float *z, const cupdlp_float *x,
                                const cupdlp_float *y, const cupdlp_int len);
extern "C" void compute_dualOT_inf_cuda(cupdlp_float *h_c_norm, cupdlp_float *h_diff_norm, cupdlp_float *x, int vec_len, int n_coarse, double scale_const);
extern "C" void countZero_and_checkConstraint_cuda(long long **h_keep_fine_redundancy, long long *keep_fine_redundancy_len, const double *h_x, const double *h_y, long long x_len, long long y_len, int resolution_now, int resolution_last, double thr, double violate_degree);
#endif
#endif