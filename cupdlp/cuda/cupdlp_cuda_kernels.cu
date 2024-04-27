#include "cupdlp_cuda_kernels.cuh"

dim3 cuda_gridsize(cupdlp_int n) {
  cupdlp_int k = (n - 1) / CUPDLP_BLOCK_SIZE + 1;
  cupdlp_int x = k;
  cupdlp_int y = 1;
  if (x > 65535) {
    x = ceil(sqrt(k));
    y = (n - 1) / (x * CUPDLP_BLOCK_SIZE) + 1;
  }
  dim3 d = {x, y, 1};
  return d;
}

__global__ void element_wise_dot_kernel(cupdlp_float *x, const cupdlp_float *y,
                                        const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] *= y[i];
}

__global__ void element_wise_div_kernel(cupdlp_float *x, const cupdlp_float *y,
                                        const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] /= y[i];
}

__global__ void element_wise_projlb_kernel(cupdlp_float *x,
                                           const cupdlp_float *lb,
                                           const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] = x[i] < lb[i] ? lb[i] : x[i];
}

__global__ void element_wise_projub_kernel(cupdlp_float *x,
                                           const cupdlp_float *ub,
                                           const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] = x[i] > ub[i] ? ub[i] : x[i];
}

__global__ void element_wise_projSamelb_kernel(cupdlp_float *x,
                                               const cupdlp_float lb,
                                               const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] = x[i] <= lb ? lb : x[i];
}

__global__ void element_wise_projSameub_kernel(cupdlp_float *x,
                                               const cupdlp_float ub,
                                               const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] = x[i] >= ub ? ub : x[i];
}

__global__ void element_wise_initHaslb_kernal(cupdlp_float *haslb,
                                              const cupdlp_float *lb,
                                              const cupdlp_float bound,
                                              const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) haslb[i] = lb[i] > bound ? 1.0 : 0.0;
}

__global__ void element_wise_initHasub_kernal(cupdlp_float *hasub,
                                              const cupdlp_float *ub,
                                              const cupdlp_float bound,
                                              const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) hasub[i] = ub[i] < bound ? 1.0 : 0.0;
}

__global__ void element_wise_filterlb_kernal(cupdlp_float *x,
                                             const cupdlp_float *lb,
                                             const cupdlp_float bound,
                                             const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] = lb[i] > bound ? lb[i] : 0.0;
}

__global__ void element_wise_filterub_kernal(cupdlp_float *x,
                                             const cupdlp_float *ub,
                                             const cupdlp_float bound,
                                             const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] = ub[i] < bound ? ub[i] : 0.0;
}

__global__ void init_cuda_vec_kernal(cupdlp_float *x, const cupdlp_float val,
                                     const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] = val;
}

//xUpdate = x - dPrimalStep * (cost - ATy)
__global__ void primal_grad_step_kernal(cupdlp_float *xUpdate,
                                        const cupdlp_float *x,
                                        const cupdlp_float *cost,
                                        const cupdlp_float *ATy,
                                        const cupdlp_float dPrimalStep,
                                        const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) xUpdate[i] = x[i] - dPrimalStep * (cost[i] - ATy[i]);
}

//yUpdate = y + dDualStep * (b -2AxUpdate + Ax)
__global__ void dual_grad_step_kernal(cupdlp_float *yUpdate,
                                      const cupdlp_float *y,
                                      const cupdlp_float *b,
                                      const cupdlp_float *Ax,
                                      const cupdlp_float *AxUpdate,
                                      const cupdlp_float dDualStep,
                                      const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x; 
  if (i < len) yUpdate[i] = y[i] + dDualStep * (b[i] - 2 * AxUpdate[i] + Ax[i]);
}

//yUpdate = y + dDualStep * (b -(theat+1)*AxUpdate + theta*Ax)
__global__ void dual_grad_step_kernal_AdapTheta(cupdlp_float *yUpdate,
                                      const cupdlp_float *y,
                                      const cupdlp_float *b,
                                      const cupdlp_float *Ax,
                                      const cupdlp_float *AxUpdate,
                                      const cupdlp_float dDualStep,
                                      const cupdlp_int len, const cupdlp_float theta) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x; 
  if (i < len) yUpdate[i] = y[i] + dDualStep * (b[i] - (theta+1)* AxUpdate[i] + theta *Ax[i]);
}



__global__ void pdtest_x_update_with_beta_kernal(cupdlp_float *x_new,
                                      const cupdlp_float *x_1,
                                      const cupdlp_float *x_2,
                                      const cupdlp_float beta_inv,
                                      const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x; 
  if (i < len) x_new[i] = (1-beta_inv) * x_1[i] + beta_inv * x_2[i];
}

__global__ void pdtest_y_update_with_beta_kernal(cupdlp_float *y_new,
                                      const cupdlp_float *y_1,
                                      const cupdlp_float *y_2,
                                      const cupdlp_float beta_inv,
                                      const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x; 
  if (i < len) y_new[i] = (1-beta_inv) * y_1[i] + beta_inv * y_2[i];
}

__global__ void pdtest_x_update_with_theta_kernal(cupdlp_float *x_barUpdate,
                                      const cupdlp_float *xUpdate,
                                      const cupdlp_float *x,
                                      const cupdlp_float theta,
                                      const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x; 
  if (i < len) x_barUpdate[i] = theta * (xUpdate[i] - x[i]) + xUpdate[i];
}
__global__ void pdtest_primal_grad_step_kernal(cupdlp_float *xUpdate,
                                        const cupdlp_float *x,
                                        const cupdlp_float *cost,
                                        const cupdlp_float *ATyUpdate,
                                        const cupdlp_float dPrimalStep,
                                        const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) xUpdate[i] = x[i] - dPrimalStep * (cost[i] - ATyUpdate[i]);
}

//yUpdate = y + dDualStep * (b -2AxUpdate + Ax)
__global__ void pdtest_dual_grad_step_kernal(cupdlp_float *yUpdate,
                                      const cupdlp_float *y,
                                      const cupdlp_float *b,
                                      const cupdlp_float *Ax_bar,
                                      const cupdlp_float dDualStep,
                                      const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x; 
  if (i < len) yUpdate[i] = y[i] + dDualStep * (b[i] - Ax_bar[i]);
}

// z = x - y
__global__ void naive_sub_kernal(cupdlp_float *z, const cupdlp_float *x,
                                  const cupdlp_float *y, const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) z[i] = x[i] - y[i];
}

__device__ double atomicAddDouble(double* address, double val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
            __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}

__global__ void compute_dualOT_inf_kernal(const double *x, int vec_len, int n_coarse, double scale_const, double *d_c_norm, double *d_diff_norm) {
  // 总共有gridDim.x * gridDim.y个block，每个block有blockDim.x * blockDim.y个线程
  // n_coarse = resolution_now，是当前处理的图片的分辨率
  // 总共要计算resolution^4个数据的求和，相当于每个块的每个线程负责resolution^4/(gridDim.x * gridDim.y * blockDim.x * blockDim.y)个数据
  // 记int x_num = (resolution * resolution) / (gridDim.x * blockDim.x), y_num同理定义
  //  x方向上的第thread_idx_x = blockIdx.x * blockDim.x + threadIdx.x个线程，y方向上的第blockIdx.y * blockDim.y + threadIdx.y个线程，负责计算的数据索引为：
  // i = thread_idx_x * x_num; i < (thread_idx_x + 1) * x_num
  // j同理
  int thread_idx_x = blockIdx.x * blockDim.x + threadIdx.x;
  int thread_idx_y = blockIdx.y * blockDim.y + threadIdx.y;
  int x_num = vec_len / (gridDim.x * blockDim.x);
  int y_num = vec_len / (gridDim.y * blockDim.y);
  double c_norm = 0;
  double diff_norm = 0;
  for (int i = thread_idx_x * x_num; i < (thread_idx_x + 1) * x_num; i++){
    for (int j = thread_idx_y * y_num; j < (thread_idx_y + 1) * y_num; j++){
      int idx1 = i / n_coarse;
      int idx2 = i % n_coarse;
      int idx3 = j / n_coarse;
      int idx4 = j % n_coarse;

      double c = ((idx1 - idx3) * (idx1 - idx3) + (idx2 - idx4) * (idx2 - idx4)) / scale_const;
      double temp = (x[i] + x[vec_len + j]) - c;
      temp = temp > 0 ? temp : 0; // 使用条件运算符替代 fmax

      c_norm += c * c;
      diff_norm += temp * temp;
    }
  }
  atomicAddDouble(d_c_norm, c_norm);
  atomicAddDouble(d_diff_norm, diff_norm);
}


__global__ void countZero_and_ckeckConstraint_before_cudaMalloc_kernal(long long *keep_local_len_array,  const double *x, const double *y, int resolution_now, int resolution_last, double thr, double violate_degree){
  // 总共有gridDim.x * gridDim.y个block，每个block有blockDim.x * blockDim.y个线程
  // 希望这些线程处理resolution_last^4个循环，那么每个线程就要处理resolution_last^4/(gridDim.x * gridDim.y * blockDim.x * blockDim.y)个循环
  // 记int x_num = (resolution_last * resolution_last) / (gridDim.x * blockDim.x), y_num同理定义
  //  x方向上的第thread_idx_x = blockIdx.x * blockDim.x + threadIdx.x个线程，y方向上的第blockIdx.y * blockDim.y + threadIdx.y个线程，负责计算的数据索引为：
  // i = thread_idx_x * x_num; i < (thread_idx_x + 1) * x_num
  // j同理
  int vec_len_last = resolution_last * resolution_last;
  int thread_idx_x = blockIdx.x * blockDim.x + threadIdx.x;
  int thread_idx_y = blockIdx.y * blockDim.y + threadIdx.y;
  int x_num = vec_len_last / (gridDim.x * blockDim.x);
  int y_num = vec_len_last / (gridDim.y * blockDim.y);
  int scale = resolution_now / resolution_last;
  long long pow_resolution_last_2 = resolution_last * resolution_last;
  long long pow_resolution_now_2 = resolution_now * resolution_now;
  double scale_constant = 2.0 * pow_resolution_now_2;
  long long keep_local_len = 0;
  // i和j的取值范围是[0, resolution_last * resolution_last)
  for (int i = thread_idx_x * x_num; i < (thread_idx_x + 1) * x_num; i++){
    for (int j = thread_idx_y * y_num; j < (thread_idx_y + 1) * y_num; j++){
      for (int k1 = 0; k1 < scale; k1++){
        for (int k2 = 0; k2 < scale; k2++){
          for (int l1 = 0; l1 < scale; l1++){
            for (int l2 = 0; l2 < scale; l2++){
              long long i1 = i / resolution_last;
              long long i2 = i % resolution_last;
              long long j1 = j / resolution_last;
              long long j2 = j % resolution_last;
              long long idx_coarse = (i1 * resolution_last + j1) * pow_resolution_last_2 + (i2 * resolution_last + j2);
              long long idx_i1 = i1 * scale + k1;
              long long idx_i2 = i2 * scale + k2;
              long long idx_j1 = j1 * scale + l1;
              long long idx_j2 = j2 * scale + l2;
              long long idx_1 = idx_i1 * resolution_now + idx_j1;
              long long idx_2 = idx_i2 * resolution_now + idx_j2;
              if (fabs(y[idx_coarse]) >=thr){
                keep_local_len++;
              }
              else{
                if (x[idx_1] + x[pow_resolution_now_2 + idx_2] > (1 + violate_degree) * ((idx_i1 - idx_i2) * (idx_i1 - idx_i2) + (idx_j1 - idx_j2) * (idx_j1 - idx_j2)) / scale_constant){
                  keep_local_len++;
                }
              }
            }
          }
        }
      }
    }
  }
  int thread_idx = thread_idx_x * gridDim.y * blockDim.y + thread_idx_y;
  keep_local_len_array[thread_idx] = keep_local_len;
}


__global__ void countZero_and_checkConstraint_kernal(long long *keep_fine_redundancy, long long *keep_len_UpToNow, const double *x, const double *y, int resolution_now, int resolution_last, double thr, double violate_degree){
  // 总共有gridDim.x * gridDim.y个block，每个block有blockDim.x * blockDim.y个线程
  // resolution_now，是当前处理的图片的分辨率
  // 总共要完成resolution_now^4个数据的判定，相当于每个块的每个线程负责resolution_now^4/(gridDim.x * gridDim.y * blockDim.x * blockDim.y)个数据
  // 记int x_num = (resolution_last * resolution_last) / (gridDim.x * blockDim.x), y_num同理定义
  //  x方向上的第thread_idx_x = blockIdx.x * blockDim.x + threadIdx.x个线程，y方向上的第blockIdx.y * blockDim.y + threadIdx.y个线程，负责计算的数据索引为：
  // i = thread_idx_x * x_num; i < (thread_idx_x + 1) * x_num
  // j同理
  int vec_len = resolution_last * resolution_last;
  int thread_idx_x = blockIdx.x * blockDim.x + threadIdx.x;
  int thread_idx_y = blockIdx.y * blockDim.y + threadIdx.y;
  int x_num = vec_len / (gridDim.x * blockDim.x);
  int y_num = vec_len / (gridDim.y * blockDim.y);
  int scale = resolution_now / resolution_last;
  long long pow_resolution_last_2 = resolution_last * resolution_last;
  long long pow_resolution_now_2 = resolution_now * resolution_now;
  double scale_constant = 2.0 * pow_resolution_now_2;
  long long keep_part_len = 0;
  int thread_idx = thread_idx_x * gridDim.y * blockDim.y + thread_idx_y;
  long long begin = keep_len_UpToNow[thread_idx];
  long long count = 0;
  for (int i = thread_idx_x * x_num; i < (thread_idx_x + 1) * x_num; i++){
    for (int j = thread_idx_y * y_num; j < (thread_idx_y + 1) * y_num; j++){
      for (int k1 = 0; k1 < scale; k1++){
        for (int k2 = 0; k2 < scale; k2++){
          for (int l1 = 0; l1 < scale; l1++){
            for (int l2 = 0; l2 < scale; l2++){
              long long i1 = i / resolution_last;
              long long i2 = i % resolution_last;
              long long j1 = j / resolution_last;
              long long j2 = j % resolution_last;
              long long idx_coarse = (i1 * resolution_last + j1) * pow_resolution_last_2 + (i2 * resolution_last + j2);
              long long idx_i1 = i1 * scale + k1;
              long long idx_i2 = i2 * scale + k2;
              long long idx_j1 = j1 * scale + l1;
              long long idx_j2 = j2 * scale + l2;
              long long idx_1 = idx_i1 * resolution_now + idx_j1;
              long long idx_2 = idx_i2 * resolution_now + idx_j2;
              if (fabs(y[idx_coarse]) >=thr){
                long long idx_fine = idx_1 * pow_resolution_now_2 + idx_2;
                keep_fine_redundancy[count + begin] = idx_fine;
                count++;
              }
              else{
                if (x[idx_1] + x[pow_resolution_now_2 + idx_2] > (1 + violate_degree) * ((idx_i1 - idx_i2) * (idx_i1 - idx_i2) + (idx_j1 - idx_j2) * (idx_j1 - idx_j2)) / scale_constant){
                  long long idx_fine = idx_1 * pow_resolution_now_2 + idx_2;
                  keep_fine_redundancy[count + begin] = idx_fine;
                  count++;
                }
              }
            }
          }
        }
      }
    }
  }
}
