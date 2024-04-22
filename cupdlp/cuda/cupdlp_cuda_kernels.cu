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
  // int i = blockIdx.x * blockDim.x + threadIdx.x;
  // int j = blockIdx.y * blockDim.y + threadIdx.y;

  // if (i < vec_len && j < vec_len)
  // {
  //   int idx1 = i / n_coarse;
  //   int idx2 = i % n_coarse;
  //   int idx3 = j / n_coarse;
  //   int idx4 = j % n_coarse;

  //   double c = ((idx1 - idx3) * (idx1 - idx3) + (idx2 - idx4) * (idx2 - idx4)) / scale_const;
  //   double temp = (x[i] + x[vec_len + j]) - c;
  //   temp = temp > 0 ? temp : 0; // 使用条件运算符替代 fmax

  //   atomicAddDouble(d_c_norm, c * c);
  //   atomicAddDouble(d_diff_norm, temp * temp);
  //   }
}

// __global__ void compute_dualOT_inf_kernal(const double *x_solution, int vec_len, int n_coarse, double scale_const, double *block_c_norm, double *block_diff_norm) {
//     extern __shared__ double shared_data[];
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     int j = blockIdx.y * blockDim.y + threadIdx.y;
//     int local_tid = threadIdx.x + blockDim.x * threadIdx.y;
//     int block_tid = threadIdx.x;

//     // Initialize shared memory
//     shared_data[local_tid] = 0;
//     shared_data[blockDim.x * blockDim.y + local_tid] = 0;
//     __syncthreads();

//     if (i < vec_len && j < vec_len) {
//         int idx1 = i / n_coarse;
//         int idx2 = i % n_coarse;
//         int idx3 = j / n_coarse;
//         int idx4 = j % n_coarse;
        
//         double c = ((idx1 - idx3) * (idx1 - idx3) + (idx2 - idx4) * (idx2 - idx4)) / scale_const;
//         // double temp = fmax((x_solution[i] + x_solution[vec_len + j]) - c, 0);
//         double temp = (x_solution[i] + x_solution[vec_len + j]) - c;
//         temp = temp > 0 ? temp : 0;  // 使用条件运算符替代 fmax
//         // Local reduction in shared memory
//         shared_data[local_tid] += c * c;
//         shared_data[blockDim.x * blockDim.y + local_tid] += temp * temp;
//     }
//     __syncthreads();

//     // Reduction within the block
//     int num_threads = blockDim.x * blockDim.y;
//     while (num_threads > 1) {
//         int half_point = (num_threads + 1) / 2;
//         if (local_tid < half_point && (local_tid + half_point) < num_threads) {
//             shared_data[local_tid] += shared_data[local_tid + half_point];
//             shared_data[blockDim.x * blockDim.y + local_tid] += shared_data[blockDim.x * blockDim.y + local_tid + half_point];
//         }
//         __syncthreads();
//         num_threads = half_point;
//     }

//     // Write result for this block to global memory
//     if (block_tid == 0) {
//         block_c_norm[blockIdx.x + gridDim.x * blockIdx.y] = shared_data[0];
//         block_diff_norm[blockIdx.x + gridDim.x * blockIdx.y] = shared_data[blockDim.x * blockDim.y];
//     }
// }

// __global__ void reduceFinal_kernal(double *block_results, double *final_result, int num_blocks) {
//     extern __shared__ double shared_data[];
//     int tid = threadIdx.x;

//     shared_data[tid] = (tid < num_blocks) ? block_results[tid] : 0;
//     __syncthreads();

//     for (int s = blockDim.x / 2; s > 0; s >>= 1) {
//         if (tid < s) {
//             shared_data[tid] += shared_data[tid + s];
//         }
//         __syncthreads();
//     }

//     if (tid == 0) {
//         *final_result = shared_data[0];
//     }
// }
