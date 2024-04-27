#include "cupdlp_cudalinalg.cuh"

extern "C" cupdlp_int cuda_alloc_MVbuffer(
    cusparseHandle_t handle, cusparseSpMatDescr_t cuda_csc,
    cusparseDnVecDescr_t vecX, cusparseDnVecDescr_t vecAx,
    cusparseSpMatDescr_t cuda_csr, cusparseDnVecDescr_t vecY,
    cusparseDnVecDescr_t vecATy, void **dBuffer) {
  cudaDataType computeType = CUDA_R_32F;
#ifndef SFLOAT
  computeType = CUDA_R_64F;
#endif

  size_t AxBufferSize = 0;
  size_t ATyBufferSize = 0;
  cupdlp_float alpha = 1.0;
  cupdlp_float beta = 0.0;
  // cusparseSpSVAlg_t alg = CUSPARSE_SPSV_ALG_DEFAULT;
  cusparseSpMVAlg_t alg = CUSPARSE_SPMV_CSR_ALG2; //deterministic

  // get the buffer size needed by csr Ax
  CHECK_CUSPARSE(cusparseSpMV_bufferSize(
      handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, cuda_csr, vecX, &beta,
      vecAx, computeType, alg, &AxBufferSize))

  // get the buffer size needed by csc ATy
  CHECK_CUSPARSE(cusparseSpMV_bufferSize(
      handle, CUSPARSE_OPERATION_TRANSPOSE, &alpha, cuda_csc, vecY, &beta,
      vecATy, computeType, alg, &ATyBufferSize))

  size_t bufferSize =
      (AxBufferSize > ATyBufferSize) ? AxBufferSize : ATyBufferSize;

  // allocate an external buffer if needed
  CHECK_CUDA(cudaMalloc(dBuffer, bufferSize))

  return EXIT_SUCCESS;
}

extern "C" cupdlp_int cuda_csc_Ax(cusparseHandle_t handle,
                                  cusparseSpMatDescr_t cuda_csc,
                                  cusparseDnVecDescr_t vecX,
                                  cusparseDnVecDescr_t vecAx, void *dBuffer,
                                  const cupdlp_float alpha,
                                  const cupdlp_float beta) {
  // hAx = alpha * Acsc * hX + beta * hAx

  cusparseOperation_t op = CUSPARSE_OPERATION_NON_TRANSPOSE;
  cudaDataType computeType = CUDA_R_32F;
#ifndef SFLOAT
  computeType = CUDA_R_64F;
#endif

  CHECK_CUSPARSE(cusparseSpMV(handle, op, &alpha, cuda_csc, vecX, &beta, vecAx,
                              // computeType, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer))
                              computeType, CUSPARSE_SPMV_CSR_ALG2, dBuffer))

  return EXIT_SUCCESS;
}

extern "C" cupdlp_int cuda_csr_Ax(cusparseHandle_t handle,
                                  cusparseSpMatDescr_t cuda_csr,
                                  cusparseDnVecDescr_t vecX,
                                  cusparseDnVecDescr_t vecAx, void *dBuffer,
                                  const cupdlp_float alpha,
                                  const cupdlp_float beta) {
  // hAx = alpha * Acsc * hX + beta * hAx

  cusparseOperation_t op = CUSPARSE_OPERATION_NON_TRANSPOSE;
  cudaDataType computeType = CUDA_R_32F;
#ifndef SFLOAT
  computeType = CUDA_R_64F;
#endif

  CHECK_CUSPARSE(cusparseSpMV(handle, op, &alpha, cuda_csr, vecX, &beta, vecAx,
                              // computeType, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer))
                              computeType, CUSPARSE_SPMV_CSR_ALG2, dBuffer))

  return EXIT_SUCCESS;
}

extern "C" cupdlp_int cuda_csc_ATy(cusparseHandle_t handle,
                                   cusparseSpMatDescr_t cuda_csc,
                                   cusparseDnVecDescr_t vecY,
                                   cusparseDnVecDescr_t vecATy, void *dBuffer,
                                   const cupdlp_float alpha,
                                   const cupdlp_float beta) {
  // hATy = alpha * Acsr^T * hY + beta * hATy
  cusparseOperation_t op = CUSPARSE_OPERATION_TRANSPOSE;
  cudaDataType computeType = CUDA_R_32F;
#ifndef SFLOAT
  computeType = CUDA_R_64F;
#endif

  CHECK_CUSPARSE(cusparseSpMV(handle, op, &alpha, cuda_csc, vecY, &beta, vecATy,
                              // computeType, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer))
                              computeType, CUSPARSE_SPMV_CSR_ALG2, dBuffer))

  return EXIT_SUCCESS;
}

extern "C" cupdlp_int cuda_csr_ATy(cusparseHandle_t handle,
                                   cusparseSpMatDescr_t cuda_csr,
                                   cusparseDnVecDescr_t vecY,
                                   cusparseDnVecDescr_t vecATy, void *dBuffer,
                                   const cupdlp_float alpha,
                                   const cupdlp_float beta) {
  // hATy = alpha * Acsr^T * hY + beta * hATy
  cusparseOperation_t op = CUSPARSE_OPERATION_TRANSPOSE;
  cudaDataType computeType = CUDA_R_32F;
#ifndef SFLOAT
  computeType = CUDA_R_64F;
#endif

  CHECK_CUSPARSE(cusparseSpMV(handle, op, &alpha, cuda_csr, vecY, &beta, vecATy,
                              // computeType, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer))
                              computeType, CUSPARSE_SPMV_CSR_ALG2, dBuffer))

  return EXIT_SUCCESS;
}

extern "C" void cupdlp_projSameub_cuda(cupdlp_float *x, const cupdlp_float ub,
                                       const cupdlp_int len) {
  element_wise_projSameub_kernel<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
      x, ub, len);
}

extern "C" void cupdlp_projSamelb_cuda(cupdlp_float *x, const cupdlp_float lb,
                                       const cupdlp_int len) {
  element_wise_projSamelb_kernel<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
      x, lb, len);
}

extern "C" void cupdlp_projub_cuda(cupdlp_float *x, const cupdlp_float *ub,
                                   const cupdlp_int len) {
  element_wise_projub_kernel<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(x, ub,
                                                                        len);
}

extern "C" void cupdlp_projlb_cuda(cupdlp_float *x, const cupdlp_float *lb,
                                   const cupdlp_int len) {
  element_wise_projlb_kernel<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(x, lb,
                                                                        len);
}

extern "C" void cupdlp_ediv_cuda(cupdlp_float *x, const cupdlp_float *y,
                                 const cupdlp_int len) {
  element_wise_div_kernel<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(x, y, len);
}

extern "C" void cupdlp_edot_cuda(cupdlp_float *x, const cupdlp_float *y,
                                 const cupdlp_int len) {
  element_wise_dot_kernel<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(x, y, len);
}

extern "C" void cupdlp_haslb_cuda(cupdlp_float *haslb, const cupdlp_float *lb,
                                  const cupdlp_float bound,
                                  const cupdlp_int len) {
  element_wise_initHaslb_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
      haslb, lb, bound, len);
}

extern "C" void cupdlp_hasub_cuda(cupdlp_float *hasub, const cupdlp_float *ub,
                                  const cupdlp_float bound,
                                  const cupdlp_int len) {
  element_wise_initHasub_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
      hasub, ub, bound, len);
}

extern "C" void cupdlp_filterlb_cuda(cupdlp_float *x, const cupdlp_float *lb,
                                     const cupdlp_float bound,
                                     const cupdlp_int len) {
  element_wise_filterlb_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
      x, lb, bound, len);
}

extern "C" void cupdlp_filterub_cuda(cupdlp_float *x, const cupdlp_float *ub,
                                     const cupdlp_float bound,
                                     const cupdlp_int len) {
  element_wise_filterub_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
      x, ub, bound, len);
}

extern "C" void cupdlp_initvec_cuda(cupdlp_float *x, const cupdlp_float val,
                                    const cupdlp_int len) {
  init_cuda_vec_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(x, val, len);
}

extern "C" void cupdlp_pgrad_cuda(cupdlp_float *xUpdate,
                                        const cupdlp_float *x,
                                        const cupdlp_float *cost,
                                        const cupdlp_float *ATy,
                                        const cupdlp_float dPrimalStep,
                                        const cupdlp_int len) {
    primal_grad_step_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
        xUpdate, x, cost, ATy, dPrimalStep, len);
}

extern "C" void cupdlp_dgrad_cuda(cupdlp_float *yUpdate, const cupdlp_float *y, const cupdlp_float *b,
    const cupdlp_float *Ax, const cupdlp_float *AxUpdate,
    const cupdlp_float dDualStep, const cupdlp_int len) {
      dual_grad_step_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
          yUpdate, y, b, Ax, AxUpdate, dDualStep, len);
}

extern "C" void cupdlp_dgrad_cuda_AdapTheta(cupdlp_float *yUpdate, const cupdlp_float *y, const cupdlp_float *b,
    const cupdlp_float *Ax, const cupdlp_float *AxUpdate,
    const cupdlp_float dDualStep, const cupdlp_int len, const cupdlp_float theta) {
      dual_grad_step_kernal_AdapTheta<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
          yUpdate, y, b, Ax, AxUpdate, dDualStep, len, theta);
}

extern "C" void pdtest_x_md_update_cuda(cupdlp_float *x_md, const cupdlp_float *x_ag, const cupdlp_float *x, const cupdlp_float beta, const cupdlp_int len) {
      pdtest_x_update_with_beta_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
          x_md, x_ag, x, 1 / beta, len);
}


extern "C" void pdtest_x_ag_update_cuda(cupdlp_float *x_agUpdate, const cupdlp_float *x_ag, const cupdlp_float *xUpdate, const cupdlp_float beta, const cupdlp_int len) {
      pdtest_x_update_with_beta_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
          x_agUpdate, x_ag, xUpdate, 1 / beta, len);
}

extern "C" void pdtest_y_ag_update_cuda(cupdlp_float *y_agUpdate, const cupdlp_float *y_ag, const cupdlp_float *yUpdate, const cupdlp_float beta, const cupdlp_int len) {
      pdtest_y_update_with_beta_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
          y_agUpdate, y_ag, yUpdate, 1 / beta, len);
}

extern "C" void pdtest_x_bar_update_cuda(cupdlp_float *x_barUpdate, const cupdlp_float *xUpdate, const cupdlp_float *x, const cupdlp_float theta, const cupdlp_int len) {
      pdtest_x_update_with_theta_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
          x_barUpdate, xUpdate, x, theta, len);
}

extern "C" void pdtest_pgrad_cuda(cupdlp_float *xUpdate,
                                        const cupdlp_float *x,
                                        const cupdlp_float *cost,
                                        const cupdlp_float *ATyUpdate,
                                        const cupdlp_float dPrimalStep,
                                        const cupdlp_int len) {
    pdtest_primal_grad_step_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
        xUpdate, x, cost, ATyUpdate, dPrimalStep, len);
}

extern "C" void pdtest_dgrad_cuda(cupdlp_float *yUpdate, const cupdlp_float *y, const cupdlp_float *b,
    const cupdlp_float *Ax_bar,
    const cupdlp_float dDualStep, const cupdlp_int len) {
      pdtest_dual_grad_step_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
          yUpdate, y, b, Ax_bar, dDualStep, len);
}

extern "C" void cupdlp_sub_cuda(cupdlp_float *z, const cupdlp_float *x,
                                  const cupdlp_float *y, const cupdlp_int len)
{
   naive_sub_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(z, x, y, len);
}

// extern "C" void compute_dualOT_inf_cuda(cupdlp_float *h_c_norm, cupdlp_float *h_diff_norm, cupdlp_float *x, int vec_len, int n_coarse, double scale_const){
//   int blocks_x = 16;
//   int blocks_y = 16;
//   dim3 blocks(blocks_x, blocks_y);
//   dim3 threads(16, 16);
//   int num_blocks = blocks_x * blocks_y;
//   double *d_block_c_norm = NULL;
//   double *d_block_diff_norm = NULL;
//   double *d_c_norm = NULL;
//   double *d_diff_norm = NULL;
//   cudaMalloc(&d_block_c_norm, num_blocks * sizeof(double));
//   cudaMalloc(&d_block_diff_norm, num_blocks * sizeof(double));
//   cudaMalloc(&d_c_norm, sizeof(double));
//   cudaMalloc(&d_diff_norm, sizeof(double));
//   cudaMemset(d_block_c_norm, 0, num_blocks * sizeof(double));
//   cudaMemset(d_block_diff_norm, 0, num_blocks * sizeof(double));
//   cudaMemset(d_c_norm, 0, sizeof(double));
//   cudaMemset(d_diff_norm, 0, sizeof(double));
//   // <<<cuda_gridsize(vec_len), CUPDLP_BLOCK_SIZE>>>(x, vec_len, n_coarse, scale_const, block_c_norm, block_diff_norm);
//   compute_dualOT_inf_kernal<<<blocks, threads>>>(x, vec_len, n_coarse, scale_const, d_block_c_norm, d_block_diff_norm);

// // Call the reduction kernel to summarize the blocks' results
//   int reduce_threads = 256; // Make sure this is a power of 2
//   reduceFinal_kernal<<<1, reduce_threads, reduce_threads * sizeof(double)>>>(d_block_c_norm, d_c_norm, num_blocks);
//   reduceFinal_kernal<<<1, reduce_threads, reduce_threads * sizeof(double)>>>(d_block_diff_norm, d_diff_norm, num_blocks);

//   cudaMemcpy(h_c_norm, d_c_norm, sizeof(double), cudaMemcpyDeviceToHost);
//   cudaMemcpy(h_diff_norm, d_diff_norm, sizeof(double), cudaMemcpyDeviceToHost);
//   printf("compute_dualOT_inf_cuda, c_norm: %f, diff_norm: %f\n", *h_c_norm, *h_diff_norm);
//   cudaFree(d_block_c_norm);
//   cudaFree(d_block_diff_norm);
//   cudaFree(d_c_norm);
//   cudaFree(d_diff_norm);

// }

extern "C" void compute_dualOT_inf_cuda(cupdlp_float *h_c_norm, cupdlp_float *h_diff_norm, cupdlp_float *x, int vec_len, int n_coarse, double scale_const){
    // int grid_size = 16;
    // int block_size = 16;
    // if (n_coarse < block_size){
    //   block_size = n_coarse;
    // }
    // if (n_coarse < grid_size){
    //   grid_size = n_coarse;
    // }
    int grid_size = 0;
    if (n_coarse < 16)
    {
      grid_size = 1;
    }
    else
    {
      grid_size = n_coarse / 16;
    }
    int block_size = 16;
    int blocks_x = grid_size;
    int blocks_y = grid_size;
    dim3 blocks(blocks_x, blocks_y);
    dim3 threads(block_size, block_size);

    double *d_c_norm = NULL;
    double *d_diff_norm = NULL;
    cudaMalloc(&d_c_norm, sizeof(double));
    cudaMalloc(&d_diff_norm, sizeof(double));
    cudaMemset(d_c_norm, 0, sizeof(double));
    cudaMemset(d_diff_norm, 0, sizeof(double));

    compute_dualOT_inf_kernal<<<blocks, threads>>>(x, vec_len, n_coarse, scale_const, d_c_norm, d_diff_norm);

    cudaMemcpy(h_c_norm, d_c_norm, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_diff_norm, d_diff_norm, sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_c_norm);
    cudaFree(d_diff_norm);
}

extern "C" void countZero_and_checkConstraint_cuda(long long **h_keep_fine_redundancy, long long *keep_fine_redundancy_len, const double *h_x, const double *h_y, long long x_len, long long y_len, int resolution_now, int resolution_last, double thr, double violate_degree){
  int grid_size = 0;
  if (resolution_last < 16)
  {
    grid_size = 1;
  }
  else
  {
    grid_size = resolution_last / 16;
  }
  int block_size = 16;
  int blocks_x = grid_size;
  int blocks_y = grid_size;
  dim3 blocks(blocks_x, blocks_y);
  dim3 threads(block_size, block_size);
  int num_threads = block_size * block_size * grid_size * grid_size;
  // printf("num_threads: %d\n", num_threads);
  
  long long *keep_local_len_array = NULL;
  cudaMalloc(&keep_local_len_array, num_threads * sizeof(long long));
  double *d_x = NULL;
  double *d_y = NULL;
  cudaMalloc(&d_x, x_len * sizeof(double));
  cudaMalloc(&d_y, y_len * sizeof(double));
  cudaMemcpy(d_x, h_x, x_len * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, h_y, y_len * sizeof(double), cudaMemcpyHostToDevice);
  // printf("countZero_and_ckeckConstraint_before_cudaMalloc_kernal开始\n");
  countZero_and_ckeckConstraint_before_cudaMalloc_kernal<<<blocks, threads>>>(keep_local_len_array, d_x, d_y, resolution_now, resolution_last, thr, violate_degree);
  // printf("countZero_and_ckeckConstraint_before_cudaMalloc_kernal结束\n");
  long long *h_keep_local_len_array = (long long *)malloc(num_threads * sizeof(long long));
  cudaMemcpy(h_keep_local_len_array, keep_local_len_array, num_threads * sizeof(long long), cudaMemcpyDeviceToHost);
  long long *h_keep_len_UpToNow = (long long *)malloc(num_threads * sizeof(long long));
  h_keep_len_UpToNow[0] = 0;
  for (int i = 1; i < num_threads; i++){
    h_keep_len_UpToNow[i] = h_keep_len_UpToNow[i - 1] + h_keep_local_len_array[i - 1];
  }
  long long *keep_len_UpToNow = NULL;
  cudaMalloc(&keep_len_UpToNow, num_threads * sizeof(long long));
  cudaMemcpy(keep_len_UpToNow, h_keep_len_UpToNow, num_threads * sizeof(long long), cudaMemcpyHostToDevice);
  // printf("keep_len_UpToNow计算完毕\n");
  *keep_fine_redundancy_len = h_keep_len_UpToNow[num_threads - 1] + h_keep_local_len_array[num_threads - 1];
  long long *d_keep_fine_reduandancy = NULL;
  cudaMalloc(&d_keep_fine_reduandancy, *keep_fine_redundancy_len * sizeof(long long));
  // printf("countZero_and_ckeckConstraint_kernal开始\n");
  countZero_and_checkConstraint_kernal<<<blocks, threads>>>(d_keep_fine_reduandancy, keep_len_UpToNow, d_x, d_y, resolution_now, resolution_last, thr, violate_degree);

  *h_keep_fine_redundancy = (long long *)malloc(*keep_fine_redundancy_len * sizeof(long long));
  cudaMemcpy(*h_keep_fine_redundancy, d_keep_fine_reduandancy, *keep_fine_redundancy_len * sizeof(long long), cudaMemcpyDeviceToHost);

  cudaFree(keep_local_len_array);
  cudaFree(keep_len_UpToNow);
  cudaFree(d_keep_fine_reduandancy);
  cudaFree(d_x);
  cudaFree(d_y);
  free(h_keep_local_len_array);
  free(h_keep_len_UpToNow);

  
}
