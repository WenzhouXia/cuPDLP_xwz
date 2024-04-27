#include "wrapper_clp.h"

// #ifdef __cplusplus
#include <ClpPresolve.hpp>

#include "ClpModel.hpp"
#include "ClpSimplex.hpp"
#include "omp.h"
using namespace std;
// extern "C" void *createModel() { return new ClpModel(); }

// extern "C" void deleteModel(void *model) {
//   if (model != NULL) delete (ClpModel *)model;
//   // free(model);
// }

// extern "C" void loadMps(void *model, const char *filename) {
//   // model is lhs <= Ax <= rhs, l <= x <= u
//   cout << std::string(filename) << endl;
//   ((ClpModel *)model)->readMps(filename, true);
// }

extern "C" void *createModel() { return new ClpSimplex(); }

extern "C" void deleteModel(void *model) {
  if (model != NULL) delete (ClpSimplex *)model;
  // free(model);
}

extern "C" void loadMps(void *model, const char *filename) {
  // model is lhs <= Ax <= rhs, l <= x <= u
  cout << "--------------------------------------------------" << endl;
  cout << "reading file..." << endl;
  cout << "\t" << std::string(filename) << endl;
  cout << "--------------------------------------------------" << endl;
  ((ClpSimplex *)model)->readMps(filename, true);
}

extern "C" void loadProblemWrapper(void* model, int numCols, int numRows, 
                            const int* start, const int* index, 
                            const double* value, const double* colLower, 
                            const double* colUpper, const double* obj, 
                            const double* rowLower, const double* rowUpper) {
        ClpSimplex* simplex = static_cast<ClpSimplex*>(model);
        simplex->loadProblem(numCols, numRows, start, index, value, 
                             colLower, colUpper, obj, rowLower, rowUpper);
}

extern "C" void loadProblem_byMatrix_Wrapper(void* model, int resolution, 
                            const double* colLower, 
                            const double* colUpper, const double* obj, 
                            const double* rowLower, const double* rowUpper) {
        ClpSimplex* simplex = static_cast<ClpSimplex*>(model);
#pragma region 按列构造稀疏矩阵
        CoinPackedMatrix mat;
        int nCols = 2 * pow(resolution,2);
        int nRows = pow(resolution, 4);
        printf("nCols: %d, nRows: %d\n", nCols, nRows);
        int vec_len = pow(resolution, 2);
        // i代表列数，j代表行数
        // 前一半的列
        for (int i = 0; i < nCols /2 ; i++){
          std::vector<int> indices;
          for (int j = i * vec_len; j < (i + 1) * vec_len; j++){
            indices.push_back(j);
          }
          std::vector<double> elements;
          for (int i = 0; i < indices.size(); i++){
            elements.push_back(1);
          }
            mat.appendCol(indices.size(), &indices[0], &elements[0]);
        }
        // 后一半的列
        for (int i = 0; i < nCols /2 ; i++){
          std::vector<int> indices;
          for (int j = 0; j < vec_len; j++){
            indices.push_back(i + j * vec_len);
          }
          std::vector<double> elements;
          for (int i = 0; i < indices.size(); i++){
            elements.push_back(1);
          }
          mat.appendCol(indices.size(), &indices[0], &elements[0]);
        }
        printf("约束矩阵构造完成\n");
#pragma endregion

        simplex->loadProblem(mat, 
                             colLower, colUpper, obj, rowLower, rowUpper);
}

extern "C" void loadProblem_byMatrix_Wrapper_parallel(void* model, int resolution, 
                            const double* colLower, 
                            const double* colUpper, const double* obj, 
                            const double* rowLower, const double* rowUpper) {
        ClpSimplex* simplex = static_cast<ClpSimplex*>(model);
#pragma region 按列构造稀疏矩阵
        int nCols = 2 * pow(resolution, 2);
        int nRows = pow(resolution, 4);
        printf("nCols: %d, nRows: %d\n", nCols, nRows);
        CoinPackedMatrix mat;
        mat.setDimensions(nRows, 0);
        int vec_len = pow(resolution, 2);
        // i代表列数，j代表行数
        // 前一半的列
        double generate_mat_time = omp_get_wtime();
        int num_threads = 32;
        if (nCols / 2 < num_threads){
          num_threads = nCols / 2;
        }
        printf("多线程开始，线程数为%d\n", num_threads);
        CoinPackedMatrix *mat_array_1 = new CoinPackedMatrix[num_threads];
        for (int i = 0; i < num_threads; i++){
          mat_array_1[i].setDimensions(nRows, 0);
        }
#pragma omp parallel num_threads(num_threads)
{
        int thread_id = omp_get_thread_num();
        int start = thread_id * nCols / (2 * num_threads);
        int end = (thread_id+1) * nCols / (2 * num_threads);
        for (long long i = start; i < end; i++)
        {
          long long idx = 0;
          std::vector<int> indices;
          std::vector<double> elements;
          for (long long j = 0; j < vec_len; j++)
          {
            idx = i * vec_len + j;
            indices.push_back(idx);
            elements.push_back(1.0);
          }
          mat_array_1[thread_id].appendCol(indices.size(), &indices[0], &elements[0]);
        }
}
        generate_mat_time = omp_get_wtime() - generate_mat_time;
        printf("生成前一半矩阵时间为：%f\n", generate_mat_time);

        for (int i = 0; i < num_threads; i++)
        {
          mat.rightAppendPackedMatrix(mat_array_1[i]);
        }


        CoinPackedMatrix *mat_array_2 = new CoinPackedMatrix[num_threads];
        for (int i = 0; i < num_threads; i++){
          mat_array_2[i].setDimensions(nRows, 0);
        } 
#pragma omp parallel num_threads(num_threads)
{
        int thread_id = omp_get_thread_num();
        int start = thread_id * nCols / (2 * num_threads);
        int end = (thread_id+1) * nCols / (2 * num_threads);
        for (long long i = start; i < end ; i++)
        {
          long long idx = 0;
          std::vector<int> indices;
          std::vector<double> elements;
          for (long long j = 0; j < vec_len; j++)
          {
            idx = i + j * vec_len;
            indices.push_back(idx);
            elements.push_back(1.0);
          }
          mat_array_2[thread_id].appendCol(indices.size(), &indices[0], &elements[0]);
        }
}        
        for (int i = 0; i < num_threads; i++){
          mat.rightAppendPackedMatrix(mat_array_2[i]);
        }
        printf("约束矩阵构造完成\n");
#pragma endregion

        simplex->loadProblem(mat, 
                             colLower, colUpper, obj, rowLower, rowUpper);
}


extern "C" void loadProblem_delete_byMatrix_Wrapper(void* model, int resolution, int* zero_idx, int* zero_idx_len,
                            const double* colLower, 
                            const double* colUpper, const double* obj, 
                            const double* rowLower, const double* rowUpper) {
        ClpSimplex* simplex = static_cast<ClpSimplex*>(model);
#pragma region 构造稀疏矩阵
        // CoinPackedMatrix mat;
        // int nCols = 2 * pow(resolution,2);
        // int nRows = pow(resolution, 4);
        // printf("未修剪的矩阵维度：nCols: %d, nRows: %d\n", nCols, nRows);
        // int vec_len = pow(resolution, 2);
        // // i代表列数，j代表行数
        // // 前一半的列
        // std::vector<double> elements;
        // for (int i = 0; i < vec_len; i++){
        //     elements.push_back(1.0);
        //   }
        // for (int i = 0; i < nCols /2 ; i++){
        //   std::vector<int> indices;
          
        //   for (int j = i * vec_len; j < (i + 1) * vec_len; j++){
        //     indices.push_back(j);
        //   }
        //     mat.appendCol(indices.size(), &indices[0], &elements[0]);
        // }
        // // 后一半的列
        // for (int i = 0; i < nCols /2 ; i++){
        //   std::vector<int> indices;
        //   for (int j = 0; j < vec_len; j++){
        //     indices.push_back(i + j * vec_len);
        //   }
        //   // std::vector<double> elements;
        //   // for (int i = 0; i < indices.size(); i++){
        //   //   elements.push_back(1.0);
        //   // }
        //   mat.appendCol(indices.size(), &indices[0], &elements[0]);
        // }
        // printf("开始修剪\n");
        // // 删除对应的列
        // mat.deleteRows(*zero_idx_len, zero_idx);

#pragma region 按行构造稀疏矩阵，不需要的直接不添加  
        // colordered = false代表行主元
        CoinPackedMatrix mat = CoinPackedMatrix(false, 0, 0);   
        int nCols = 2 * pow(resolution,2);
        int nRows = pow(resolution, 4);
        printf("如果未修剪，矩阵尺寸为, nCols: %d, nRows: %d\n", nCols, nRows);
        int vec_len = pow(resolution, 2);
        int row_idx = 0;
        int count = 0;
        // i代表该行中第一个1的列数，j代表第二个1的列数（-vec_len)
        // for (int i = 0; i < *zero_idx_len; i++){
        //   printf("zero_idx: %d\n", zero_idx[i]);
        // }
          for (int i = 0; i < vec_len; i++)
          {
            for (int j = 0; j < vec_len; j++)
            {
              row_idx = i * vec_len + j;
              if (row_idx == zero_idx[count])
              {
                count += 1;
              }
              else
              {
                std::vector<int> indices;
                indices.push_back(i);
                indices.push_back(j + vec_len);
                std::vector<double> elements;
                elements.push_back(1.0);
                elements.push_back(1.0);
                mat.appendRow(indices.size(), &indices[0], &elements[0]);
                // printf("添加成功：zero_idx: %d, row_idx: %d\n", zero_idx[count], row_idx);

              }
            }
          }
#pragma endregion

        printf("修剪后的矩阵维度：nCols: %d, nRows: %d\n", mat.getNumCols(), mat.getNumRows());
#pragma endregion
        simplex->loadProblem(mat, 
                             colLower, colUpper, obj, rowLower, rowUpper);
}

extern "C" void loadProblem_delete_byMatrix_Wrapper_longlong(void* model, int resolution, long long* zero_idx, long long* zero_idx_len,
                            const double* colLower, 
                            const double* colUpper, const double* obj, 
                            const double* rowLower, const double* rowUpper) {
        ClpSimplex* simplex = static_cast<ClpSimplex*>(model);
#pragma region 按行构造稀疏矩阵，不需要的直接不添加  
        // // colordered = false代表行主元
        // CoinPackedMatrix mat = CoinPackedMatrix(false, 0, 0);   
        // int nCols = 2 * pow(resolution,2);
        // long long nRows = pow(resolution, 4);
        // printf("如果未修剪，矩阵尺寸为, nCols: %d, nRows: %lld\n", nCols, nRows);
        // int vec_len = pow(resolution, 2);
        // long long row_idx = 0;
        // long long count = 0;
        // // i代表该行中第一个1的列数，j代表第二个1的列数（-vec_len)
        // for (int i = 0; i < vec_len; i++)
        // {
        //   if (i % 100 == 0){
        //     printf("i: %d\n", i);
        //   }
        //   for (int j = 0; j < vec_len; j++)
        //   {
        //     row_idx = i * vec_len + j;
        //     if (row_idx == zero_idx[count])
        //     {
        //       count += 1;
        //     }
        //     else
        //     {
        //       std::vector<int> indices;
        //       indices.push_back(i);
        //       indices.push_back(j + vec_len);
        //       std::vector<double> elements;
        //       elements.push_back(1.0);
        //       elements.push_back(1.0);
        //       mat.appendRow(indices.size(), &indices[0], &elements[0]);
        //       // printf("添加成功：zero_idx: %d, row_idx: %d\n", zero_idx[count], row_idx);
        //     }
        //   }
        // }
#pragma endregion
#pragma region 按列构造稀疏矩阵，每列中(不需要的行对应的元素)直接不添加  
        // colordered = false代表行主元
        // CoinPackedMatrix mat;
        // int nCols = 2 * pow(resolution,2);
        // long long nRows = pow(resolution, 4);
        // printf("如果未修剪，矩阵尺寸为, nCols: %d, nRows: %lld\n", nCols, nRows);
        // int vec_len = pow(resolution, 2);
        // // 定义一个keep来标记每个元素是否保存
        // bool *keep = (bool *)malloc(sizeof(bool) * nRows);
        // int *keep_true_idx = (int *)malloc(sizeof(int) * nRows);
        // for (long long i = 0; i < nRows; i++){
        //   keep[i] = true;
        //   keep_true_idx[i] = 0;
        // }
        // for (long long i = 0; i < *zero_idx_len; i++){
        //   keep[zero_idx[i]] = false;
        // }
        // int true_count = 0;
        // for (long long i = 0; i < nRows; i++)
        // {
        //   if (keep[i]){
        //     keep_true_idx[i] = true_count;
        //     true_count += 1;
        //   }
        // }
        // // 再定义一个keep_true_idx，记录keep中的true是第几个true
        // long long idx = 0;
        // int print_const = 1;
        // for (int i = 0; i < nCols / 2; i++)
        // {
        //   if (i % print_const == 0){
        //     printf("i: %d\n", i);
        //   }
        //   std::vector<int> indices;
        //   std::vector<double> elements;
        //   for (int j = 0; j < vec_len; j++)
        //   {
        //     idx = i * vec_len + j;
        //     if (keep[idx]) 
        //     {
        //       // printf("idx: %lld\n", idx);
        //       indices.push_back(keep_true_idx[idx]);
        //       elements.push_back(1.0);
        //     }
        //   }
        //   mat.appendCol(indices.size(), &indices[0], &elements[0]);
          // printf("indices_len: %lld\n", indices.size());
          // printf("修剪后的矩阵维度：nCols: %d, nRows: %d\n", mat.getNumCols(), mat.getNumRows());
        // }
        // free(keep);
        // keep = NULL;
        // free(keep_true_idx);
        // keep_true_idx = NULL;
#pragma endregion
#pragma region 按列构造稀疏矩阵，每列中(不需要的行对应的元素)直接不添加，只用keep_true_idx不用keep
        CoinPackedMatrix mat;
        int nCols = 2 * pow(resolution,2);
        long long nRows = pow(resolution, 4);
        printf("如果未修剪，矩阵尺寸为, nCols: %d, nRows: %lld\n", nCols, nRows);
        int vec_len = pow(resolution, 2);
        // 定义一个keep来标记每个元素是否保存
        int *keep_true_idx = (int *)malloc(sizeof(int) * nRows);
        for (long long i = 0; i < nRows; i++){
          keep_true_idx[i] = 0;
        }
        for (long long i = 0; i < *zero_idx_len; i++){
          keep_true_idx[zero_idx[i]] = -1;
        }
        int true_count = 0;
        for (long long i = 0; i < nRows; i++)
        {
          if (keep_true_idx[i] != -1){
            keep_true_idx[i] = true_count;
            true_count += 1;
          }
        }
        // 再定义一个keep_true_idx，记录keep中的true是第几个true
        long long idx = 0;
        int print_const = floor(nCols / 10);
        for (long long i = 0; i < nCols / 2; i++)
        {
          std::vector<int> indices;
          std::vector<double> elements;
          for (long long j = 0; j < vec_len; j++)
          {
            idx = i * vec_len + j;
            if (keep_true_idx[idx] != -1) 
            {
              // printf("idx: %lld\n", idx);
              indices.push_back(keep_true_idx[idx]);
              elements.push_back(1.0);
            }
          }
          mat.appendCol(indices.size(), &indices[0], &elements[0]);
          if (i % print_const == 0){
            printf("i: %lld\n", i);
          }
          // printf("indices_len: %lld\n", indices.size());
          // printf("修剪后的矩阵维度：nCols: %d, nRows: %d\n", mat.getNumCols(), mat.getNumRows());
        }
        // 后一半的列
        for (long long i = 0; i < nCols / 2 ; i++){
          std::vector<int> indices;
          std::vector<double> elements;
          for (long long j = 0; j < vec_len; j++){
            idx = i + j * vec_len;
            if (keep_true_idx[idx] != -1){
              indices.push_back(keep_true_idx[idx]);
              elements.push_back(1.0);
            }
          }
          mat.appendCol(indices.size(), &indices[0], &elements[0]);
          if (i % print_const == 0){
            printf("i: %lld\n", i + nCols / 2);
          }

        }
          // printf("indices_len: %lld\n", indices.size());
          // printf("修剪后的矩阵维度：nCols: %d, nRows: %d\n", mat.getNumCols(), mat.getNumRows());
        // }
        free(keep_true_idx);
        keep_true_idx = NULL;
#pragma endregion
        printf("修剪后的矩阵维度：nCols: %d, nRows: %d\n", mat.getNumCols(), mat.getNumRows());
        simplex->loadProblem(mat, colLower, colUpper, obj, rowLower, rowUpper);
        printf("loadProblem_delete_byMatrix_Wrapper_longlong\n");
}

extern "C" void loadProblem_delete_byMatrix_Wrapper_longlong_parallel(void* model, int resolution, long long* zero_idx, long long* zero_idx_len,
                            const double* colLower, 
                            const double* colUpper, const double* obj, 
                            const double* rowLower, const double* rowUpper) {
        ClpSimplex* simplex = static_cast<ClpSimplex*>(model);
#pragma region 按列构造稀疏矩阵，不用添加列的方式而是用修改列的方式，方便并行，每列中(不需要的行对应的元素)直接不添加，只用keep_true_idx不用keep
        // CoinPackedMatrix matTest1;
        // CoinPackedMatrix matTest2;
        // std::vector<int> indices_test;
        // std::vector<double> elements_test;
        // indices_test.push_back(0);
        // indices_test.push_back(1);
        // elements_test.push_back(1.0);
        // elements_test.push_back(1.0);
        // matTest1.appendCol(indices_test.size(), &indices_test[0], &elements_test[0]);
        // matTest1.appendCol(indices_test.size(), &indices_test[0], &elements_test[0]);
        // matTest2.appendCol(indices_test.size(), &indices_test[0], &elements_test[0]);
        // matTest2.appendCol(indices_test.size(), &indices_test[0], &elements_test[0]);
        // matTest1.rightAppendPackedMatrix(matTest2);
        // printf("matTest1: %d, %d\n", matTest1.getNumCols(), matTest1.getNumRows());
        
        // int sum = 0;
        // int n_threads = 4;
        // int range = 2500; // 每个线程计算的范围大小
        // #pragma omp parallel num_threads(n_threads)
        // {
        //     int thread_id = omp_get_thread_num();
        //     int start = thread_id * range;
        //     int end = start + range - 1;
        //     int thread_sum = 0;
        //     for (int i = start; i <= end; i++) {
        //         thread_sum += i;
        //     }
        //     printf("Thread %d: %d\n", thread_id, thread_sum);
        //     #pragma omp atomic
        //     sum += thread_sum;
        // }
        // printf("Total sum is %d\n", sum);
        int nCols = 2 * pow(resolution,2);
        long long nRows_before_delete = pow(resolution, 4);
        int nRows = nRows_before_delete - *zero_idx_len;
        CoinPackedMatrix mat;
        mat.setDimensions(nRows, 0);
        printf("loadProblem_delete_byMatrix_Wrapper_longlong_parallel开始\n");
        printf("如果未修剪，矩阵尺寸为, nCols: %d, nRows: %lld\n", nCols, nRows_before_delete);
        // mat.setDimensions(nRows, nCols);
        int vec_len = pow(resolution, 2); 
        double keep_true_idx_time = omp_get_wtime();
        // 定义一个keep_true_idx来标记每个元素是否保存，如果保存，是第几个true
        int *keep_true_idx = (int *)malloc(sizeof(int) * nRows_before_delete);
        for (long long i = 0; i < nRows_before_delete; i++){
          keep_true_idx[i] = 0;
        }
        for (long long i = 0; i < *zero_idx_len; i++){
          keep_true_idx[zero_idx[i]] = -1;
        }
        int true_count = 0;
        for (long long i = 0; i < nRows_before_delete; i++)
        {
          if (keep_true_idx[i] != -1){
            keep_true_idx[i] = true_count;
            true_count += 1;
          }
        }
        keep_true_idx_time = omp_get_wtime() - keep_true_idx_time;
        printf("在loadProblem内部生成keep_true_idx耗时为：%f\n", keep_true_idx_time);
        double generate_mat_time = omp_get_wtime();
        // 多线程构造前一半的列，并存储在子矩阵中，最后通过引用，进行合并
        int num_threads = 32;
        printf("多线程开始，线程数为%d\n", num_threads);
        // 设置线程数
        CoinPackedMatrix *mat_array_1 = new CoinPackedMatrix[num_threads];
        for (int i = 0; i < num_threads; i++){
          mat_array_1[i].setDimensions(nRows, 0);
        } 
#pragma omp parallel num_threads(num_threads)
{
        int thread_id = omp_get_thread_num();
        int start = thread_id * nCols / (2 * num_threads);
        int end = (thread_id+1) * nCols / (2 * num_threads);
        // 前一半的列      
        for (long long i = start; i < end; i++)
        {
          long long idx = 0;
          std::vector<int> indices;
          std::vector<double> elements;
          for (long long j = 0; j < vec_len; j++)
          {
            idx = i * vec_len + j;
            if (keep_true_idx[idx] != -1) 
            {
              // printf("idx: %lld\n", idx);
              indices.push_back(keep_true_idx[idx]);
              elements.push_back(1.0);
            }
          }
          mat_array_1[thread_id].appendCol(indices.size(), &indices[0], &elements[0]);
        }
        generate_mat_time = omp_get_wtime() - generate_mat_time;
        
}
        printf("生成前一半矩阵时间为：%f\n", generate_mat_time);  
        printf("多线程结束，开始合并\n");
        for (int i = 0; i < num_threads; i++)
        {
          // printf("第%d个线程的矩阵维度：nCols: %d, nRows: %d\n", i, mat_array_1[i].getNumCols(), mat_array_1[i].getNumRows());
          mat.rightAppendPackedMatrix(mat_array_1[i]);
        }

        CoinPackedMatrix *mat_array_2 = new CoinPackedMatrix[num_threads];
        for (int i = 0; i < num_threads; i++){
          mat_array_2[i].setDimensions(nRows, 0);
        } 
#pragma omp parallel num_threads(num_threads)
{
        int thread_id = omp_get_thread_num();
        int start = thread_id * nCols / (2 * num_threads);
        int end = (thread_id+1) * nCols / (2 * num_threads);
        for (long long i = start; i < end ; i++)
        {
          long long idx = 0;
          std::vector<int> indices;
          std::vector<double> elements;
          for (long long j = 0; j < vec_len; j++)
          {
            idx = i + j * vec_len;
            if (keep_true_idx[idx] != -1){
              indices.push_back(keep_true_idx[idx]);
              elements.push_back(1.0);
            }
          }
          mat_array_2[thread_id].appendCol(indices.size(), &indices[0], &elements[0]);
        }
}
        printf("多线程结束，开始合并\n");
        for (int i = 0; i < num_threads; i++){
          // printf("第%d个线程的矩阵维度：nCols: %d, nRows: %d\n", i, mat_array_2[i].getNumCols(), mat_array_2[i].getNumRows());
          mat.rightAppendPackedMatrix(mat_array_2[i]);
        }
        

        free(keep_true_idx);
        keep_true_idx = NULL;

#pragma endregion
        printf("修剪后的矩阵维度：nCols: %d, nRows: %d\n", mat.getNumCols(), mat.getNumRows());
        double simplex_loadProblem_time = omp_get_wtime(); 
        simplex->loadProblem(mat, colLower, colUpper, obj, rowLower, rowUpper);
        delete[] mat_array_1;
        delete[] mat_array_2;
        printf("loadProblem_delete_byMatrix_Wrapper_longlong_parallel\n");
        simplex_loadProblem_time = omp_get_wtime() - simplex_loadProblem_time;
        printf("simplex_loadProblem_time耗时为：%f\n", simplex_loadProblem_time);

        
}

extern "C" void loadProblem_delete_byMatrix_byKeepIdx_Wrapper_longlong_parallel(void* model, int resolution, long long* keep_idx, long long *keep_nnz,
                            const double* colLower, 
                            const double* colUpper, const double* obj, 
                            const double* rowLower, const double* rowUpper) {
        ClpSimplex* simplex = static_cast<ClpSimplex*>(model);
#pragma region 按列构造稀疏矩阵，不用添加列的方式而是用修改列的方式，方便并行，每列中(不需要的行对应的元素)直接不添加，只用keep_idx不用keep
        int nCols = 2 * pow(resolution,2);
        long long nRows_before_delete = pow(resolution, 4);
        int nRows = *keep_nnz;
        CoinPackedMatrix mat;
        mat.setDimensions(nRows, 0);
        printf("loadProblem_delete_byMatrix_Wrapper_longlong_parallel开始\n");
        printf("如果未修剪，矩阵尺寸为, nCols: %d, nRows: %lld\n", nCols, nRows_before_delete);
        // mat.setDimensions(nRows, nCols);
        int vec_len = pow(resolution, 2); 

        double generate_mat_time = omp_get_wtime();
        double generate_mat_all_time = omp_get_wtime();
        // 多线程构造前一半的列，并存储在子矩阵中，最后通过引用，进行合并
        int num_threads = 128;
        if (nCols / 2 < num_threads){
          num_threads = nCols / 2;
        }
        printf("多线程开始，线程数为%d\n", num_threads);
        CoinPackedMatrix *mat_array_1 = new CoinPackedMatrix[num_threads];
        for (int i = 0; i < num_threads; i++){
          mat_array_1[i].setDimensions(nRows, 0);
        } 
#pragma omp parallel num_threads(num_threads)
{
        int thread_id = omp_get_thread_num();
        int start = thread_id * nCols / (2 * num_threads);
        int end = (thread_id+1) * nCols / (2 * num_threads);
        // 前一半的列      
        for (long long i = start; i < end; i++)
        {
          if(i % 1000 == 0){
            printf("i: %lld\n", i);
          }
          long long idx = 0;
          std::vector<int> indices;
          std::vector<double> elements;
          for (long long j = 0; j < vec_len; j++)
          {
            idx = i * vec_len + j;
            if (keep_idx[idx] != -1) 
            {
              // printf("idx: %lld\n", idx);
              indices.push_back(keep_idx[idx]);
              elements.push_back(1.0);
            }
          }
          mat_array_1[thread_id].appendCol(indices.size(), &indices[0], &elements[0]);
          
        }        
}       
        generate_mat_time = omp_get_wtime() - generate_mat_time;
        printf("生成前一半矩阵时间为：%f\n", generate_mat_time);  
        printf("多线程结束，开始合并\n");
        double merge_mat_time = omp_get_wtime();
        for (int i = 0; i < num_threads; i++)
        {
          // printf("第%d个线程的矩阵维度：nCols: %d, nRows: %d\n", i, mat_array_1[i].getNumCols(), mat_array_1[i].getNumRows());
          mat.rightAppendPackedMatrix(mat_array_1[i]);
        }
        merge_mat_time = omp_get_wtime() - merge_mat_time;
        printf("合并前一半矩阵完成，耗时%f\n", merge_mat_time);

        CoinPackedMatrix *mat_array_2 = new CoinPackedMatrix[num_threads];
        for (int i = 0; i < num_threads; i++){
          mat_array_2[i].setDimensions(nRows, 0);
        } 
#pragma omp parallel num_threads(num_threads)
{
        int thread_id = omp_get_thread_num();
        int start = thread_id * nCols / (2 * num_threads);
        int end = (thread_id+1) * nCols / (2 * num_threads);
        for (long long i = start; i < end ; i++)
        {
          if(i % 1000 == 0){
            printf("i: %lld\n", i);
          }
          long long idx = 0;
          std::vector<int> indices;
          std::vector<double> elements;
          for (long long j = 0; j < vec_len; j++)
          {
            idx = i + j * vec_len;
            if (keep_idx[idx] != -1){
              indices.push_back(keep_idx[idx]);
              elements.push_back(1.0);
            }
          }
          mat_array_2[thread_id].appendCol(indices.size(), &indices[0], &elements[0]);
        }
}
        printf("多线程结束，开始合并\n");
        for (int i = 0; i < num_threads; i++){
          // printf("第%d个线程的矩阵维度：nCols: %d, nRows: %d\n", i, mat_array_2[i].getNumCols(), mat_array_2[i].getNumRows());
          mat.rightAppendPackedMatrix(mat_array_2[i]);
        }
        generate_mat_all_time = omp_get_wtime() - generate_mat_all_time;
        printf("生成矩阵完成，总共耗时%f\n", generate_mat_all_time);
        


#pragma endregion
        printf("修剪后的矩阵维度：nCols: %d, nRows: %d\n", mat.getNumCols(), mat.getNumRows());
        // 打印矩阵进行检查
        // const int numRows = mat.getNumRows();
        // const int numCols = mat.getNumCols();
        // printf("numRows: %d, numCols: %d\n", numRows, numCols);

        // // 获取矩阵中的元素和它们的索引
        // const double* elements_print = mat.getElements();
        // const int* rowIndices = mat.getIndices();
        // const CoinBigIndex* columnStarts = mat.getVectorStarts();
        // const int* lengths = mat.getVectorLengths();

        // // 按列遍历矩阵
        // for (int col = 0; col < numCols; ++col) {
        //     std::cout << "Column " << col << ":" << std::endl;
        //     // 获取这一列的起始位置和长度
        //     CoinBigIndex start = columnStarts[col];
        //     int length = lengths[col];

        //     // 遍历这一列的每个元素
        //     for (int i = 0; i < length; ++i) {
        //         // 计算当前元素的索引
        //         CoinBigIndex index = start + i;
        //         // 输出行索引和元素的值
        //         std::cout << "  Row " << rowIndices[index] << ": " << elements_print[index] << std::endl;
        //     }
        // }

        double simplex_loadProblem_time = omp_get_wtime(); 
        simplex->loadProblem(mat, colLower, colUpper, obj, rowLower, rowUpper);
        delete[] mat_array_1;
        delete[] mat_array_2;
        simplex_loadProblem_time = omp_get_wtime() - simplex_loadProblem_time;
        printf("simplex_loadProblem_time耗时为：%f\n", simplex_loadProblem_time);
        printf("loadProblem_delete_byMatrix_Wrapper_longlong_parallel完成\n");


        
}


extern "C" void generate_ConstraintMatrix_byKeepIdx_Wrapper_longlong_parallel(void* model, int resolution, long long* keep_idx, long long *keep_nnz){
        ClpSimplex* simplex = static_cast<ClpSimplex*>(model);
#pragma region 按列构造稀疏矩阵，不用添加列的方式而是用修改列的方式，方便并行，每列中(不需要的行对应的元素)直接不添加，只用keep_idx不用keep
        int nCols = 2 * pow(resolution,2);
        long long nRows_before_delete = pow(resolution, 4);
        int nRows = *keep_nnz;
        CoinPackedMatrix mat;
        mat.setDimensions(nRows, 0);
        printf("loadProblem_delete_byMatrix_Wrapper_longlong_parallel开始\n");
        printf("如果未修剪，矩阵尺寸为, nCols: %d, nRows: %lld\n", nCols, nRows_before_delete);
        // mat.setDimensions(nRows, nCols);
        int vec_len = pow(resolution, 2); 

        double generate_mat_time = omp_get_wtime();
        double generate_mat_all_time = omp_get_wtime();
        // 多线程构造前一半的列，并存储在子矩阵中，最后通过引用，进行合并
        int num_threads = 64;
        if (nCols / 2 < num_threads){
          num_threads = nCols / 2;
        }
        printf("多线程开始，线程数为%d\n", num_threads);
        CoinPackedMatrix *mat_array_1 = new CoinPackedMatrix[num_threads];
        for (int i = 0; i < num_threads; i++){
          mat_array_1[i].setDimensions(nRows, 0);
        } 
#pragma omp parallel num_threads(num_threads)
{
        int thread_id = omp_get_thread_num();
        int start = thread_id * nCols / (2 * num_threads);
        int end = (thread_id+1) * nCols / (2 * num_threads);
        // 前一半的列      
        for (long long i = start; i < end; i++)
        {
          if(i % 1000 == 0){
            printf("i: %lld\n", i);
          }
          long long idx = 0;
          std::vector<int> indices;
          std::vector<double> elements;
          for (long long j = 0; j < vec_len; j++)
          {
            idx = i * vec_len + j;
            if (keep_idx[idx] != -1) 
            {
              // printf("idx: %lld\n", idx);
              indices.push_back(keep_idx[idx]);
              elements.push_back(1.0);
            }
          }
          mat_array_1[thread_id].appendCol(indices.size(), &indices[0], &elements[0]);
          
        }        
}       
        generate_mat_time = omp_get_wtime() - generate_mat_time;
        printf("生成前一半矩阵时间为：%f\n", generate_mat_time);  
        printf("多线程结束，开始合并\n");
        double merge_mat_time = omp_get_wtime();
        for (int i = 0; i < num_threads; i++)
        {
          // printf("第%d个线程的矩阵维度：nCols: %d, nRows: %d\n", i, mat_array_1[i].getNumCols(), mat_array_1[i].getNumRows());
          mat.rightAppendPackedMatrix(mat_array_1[i]);
        }
        merge_mat_time = omp_get_wtime() - merge_mat_time;
        printf("合并前一半矩阵完成，耗时%f\n", merge_mat_time);

        CoinPackedMatrix *mat_array_2 = new CoinPackedMatrix[num_threads];
        for (int i = 0; i < num_threads; i++){
          mat_array_2[i].setDimensions(nRows, 0);
        } 
#pragma omp parallel num_threads(num_threads)
{
        int thread_id = omp_get_thread_num();
        int start = thread_id * nCols / (2 * num_threads);
        int end = (thread_id+1) * nCols / (2 * num_threads);
        for (long long i = start; i < end ; i++)
        {
          if(i % 1000 == 0){
            printf("i: %lld\n", i);
          }
          long long idx = 0;
          std::vector<int> indices;
          std::vector<double> elements;
          for (long long j = 0; j < vec_len; j++)
          {
            idx = i + j * vec_len;
            if (keep_idx[idx] != -1){
              indices.push_back(keep_idx[idx]);
              elements.push_back(1.0);
            }
          }
          mat_array_2[thread_id].appendCol(indices.size(), &indices[0], &elements[0]);
        }
}
        printf("多线程结束，开始合并\n");
        for (int i = 0; i < num_threads; i++){
          // printf("第%d个线程的矩阵维度：nCols: %d, nRows: %d\n", i, mat_array_2[i].getNumCols(), mat_array_2[i].getNumRows());
          mat.rightAppendPackedMatrix(mat_array_2[i]);
        }
        generate_mat_all_time = omp_get_wtime() - generate_mat_all_time;
        printf("生成矩阵完成，总共耗时%f\n", generate_mat_all_time);
        


#pragma endregion
        printf("修剪后的矩阵维度：nCols: %d, nRows: %d\n", mat.getNumCols(), mat.getNumRows());
        void *mat_ptr = &mat;
        printSparseMatrix(mat_ptr);
        free(mat_ptr);
        mat_ptr = NULL;
        // // 打印矩阵进行检查
        // const int numRows = mat.getNumRows();
        // const int numCols = mat.getNumCols();
        // printf("numRows: %d, numCols: %d\n", numRows, numCols);

        // // 获取矩阵中的元素和它们的索引
        // const double* elements_print = mat.getElements();
        // const int* rowIndices = mat.getIndices();
        // const CoinBigIndex* columnStarts = mat.getVectorStarts();
        // const int* lengths = mat.getVectorLengths();

        // // 按列遍历矩阵
        // for (int col = 0; col < numCols; ++col) {
        //     std::cout << "Column " << col << ":" << std::endl;
        //     // 获取这一列的起始位置和长度
        //     CoinBigIndex start = columnStarts[col];
        //     int length = lengths[col];

        //     // 遍历这一列的每个元素
        //     for (int i = 0; i < length; ++i) {
        //         // 计算当前元素的索引
        //         CoinBigIndex index = start + i;
        //         // 输出行索引和元素的值
        //         std::cout << "  Row " << rowIndices[index] << ": " << elements_print[index] << std::endl;
        //     }
        // }

          delete[] mat_array_1;
        delete[] mat_array_2;
        printf("loadProblem_delete_byMatrix_Wrapper_longlong_parallel完成\n");


        
}

extern "C" void printSparseMatrix(void *mat_ptr) {

    CoinPackedMatrix* mat_ptr_temp = static_cast<CoinPackedMatrix*>(mat_ptr);
    CoinPackedMatrix mat = *mat_ptr_temp;
    const int numRows = mat.getNumRows();
    const int numCols = mat.getNumCols();
    const double *elements_print = mat.getElements();
    const int *rowIndices = mat.getIndices();
    const CoinBigIndex* columnStarts = mat.getVectorStarts();
    const int* lengths = mat.getVectorLengths();

    // 创建并初始化二维数组
    std::vector<std::vector<double>> matrix(numRows, std::vector<double>(numCols, 0.0));

    // 填充非零元素
    for (int col = 0; col < numCols; ++col) {
        int start = columnStarts[col];
        int length = lengths[col];
        for (int i = 0; i < length; ++i) {
            int row = rowIndices[start + i];
            double value = elements_print[start + i];
            matrix[row][col] = value;
        }
    }

    // 打印矩阵
    for (const auto &row : matrix) {
        for (double val : row) {
            std::cout << val << " ";
        }
        std::cout << "\n";
    }
    free(mat_ptr_temp);
    mat_ptr_temp = NULL;
}

extern "C" void loadProblem_delete_by_keep_byMatrix_Wrapper_longlong(void* model, int resolution, bool *keep, long long *keep_true_idx, long long *len_after_delete,
                            const double* colLower, 
                            const double* colUpper, const double* obj, 
                            const double* rowLower, const double* rowUpper) {
        ClpSimplex* simplex = static_cast<ClpSimplex*>(model);
#pragma region 按行构造稀疏矩阵，不需要的直接不添加  
        // // colordered = false代表行主元
        // CoinPackedMatrix mat = CoinPackedMatrix(false, 0, 0);   
        // int nCols = 2 * pow(resolution,2);
        // long long nRows = pow(resolution, 4);
        // printf("如果未修剪，矩阵尺寸为, nCols: %d, nRows: %lld\n", nCols, nRows);
        // int vec_len = pow(resolution, 2);
        // long long row_idx = 0;
        // long long count = 0;
        // // i代表该行中第一个1的列数，j代表第二个1的列数（-vec_len)
        // for (int i = 0; i < vec_len; i++)
        // {
        //   if (i % 100 == 0){
        //     printf("i: %d\n", i);
        //   }
        //   for (int j = 0; j < vec_len; j++)
        //   {
        //     row_idx = i * vec_len + j;
        //     if (row_idx == zero_idx[count])
        //     {
        //       count += 1;
        //     }
        //     else
        //     {
        //       std::vector<int> indices;
        //       indices.push_back(i);
        //       indices.push_back(j + vec_len);
        //       std::vector<double> elements;
        //       elements.push_back(1.0);
        //       elements.push_back(1.0);
        //       mat.appendRow(indices.size(), &indices[0], &elements[0]);
        //       // printf("添加成功：zero_idx: %d, row_idx: %d\n", zero_idx[count], row_idx);
        //     }
        //   }
        // }
#pragma endregion
#pragma region 按列构造稀疏矩阵，每列中(不需要的行对应的元素)直接不添加  
        // colordered = false代表行主元
        // for (long long i = 0; i < *len_after_delete; i++){
        //   if (keep[i]){
        //     printf("keep_true_idx: %lld\n", keep_true_idx[i]);
        //   }
        // }
          CoinPackedMatrix mat;
        int nCols = 2 * pow(resolution,2);
        long long nRows = pow(resolution, 4);
        printf("如果未修剪，矩阵尺寸为, nCols: %d, nRows: %lld\n", nCols, nRows);
        int vec_len = pow(resolution, 2);
        // 定义一个keep来标记每个元素是否保存
        // bool *keep = (bool *)malloc(sizeof(bool) * nRows);
        // int *keep_true_idx = (int *)malloc(sizeof(int) * nRows);
        // for (long long i = 0; i < nRows; i++){
        //   keep[i] = true;
        //   keep_true_idx[i] = 0;
        // }
        // for (long long i = 0; i < *zero_idx_len; i++){
        //   keep[zero_idx[i]] = false;
        // }
        // int true_count = 0;
        // for (long long i = 0; i < nRows; i++)
        // {
        //   if (keep[i]){
        //     keep_true_idx[i] = true_count;
        //     true_count += 1;
        //   }
        // }
        // 再定义一个keep_true_idx，记录keep中的true是第几个true
        long long idx = 0;
        for (int i = 0; i < nCols / 2; i++)
        {
          if (i % 1000 == 0){
            printf("i: %d\n", i);
          }
          std::vector<int> indices;
          std::vector<double> elements;
          for (int j = 0; j < vec_len; j++)
          {
            idx = i * vec_len + j;
            if (keep[idx]) 
            {
              // printf("idx: %lld\n", idx);
              indices.push_back(keep_true_idx[idx]);
              elements.push_back(1.0);
            }
          }
          mat.appendCol(indices.size(), &indices[0], &elements[0]);
          // printf("修剪后的矩阵维度：nCols: %d, nRows: %d\n", mat.getNumCols(), mat.getNumRows());
        }
        // 后一半的列
        for (int i = 0; i < nCols /2 ; i++){
          if (i % 1000 == 0){
            printf("i: %d\n", i);
          }
          std::vector<int> indices;
          std::vector<double> elements;
          for (int j = 0; j < vec_len; j++){
            long long idx = i + j * vec_len;
            if (keep[idx]){
              indices.push_back(keep_true_idx[idx]);
              elements.push_back(1.0);
            }
          }
          mat.appendCol(indices.size(), &indices[0], &elements[0]);
          // printf("indices_len: %lld\n", indices.size());
          // printf("修剪后的矩阵维度：nCols: %d, nRows: %d\n", mat.getNumCols(), mat.getNumRows());
        }
#pragma endregion

        printf("修剪后的矩阵维度：nCols: %d, nRows: %d\n", mat.getNumCols(), mat.getNumRows());
        simplex->loadProblem(mat, colLower, colUpper, obj, rowLower, rowUpper);
}

extern "C" void writeLpWrapper(void* model, const char* filename) {
    ClpSimplex* simplex = static_cast<ClpSimplex*>(model);
    simplex->writeLp(filename);
}

extern "C" void *createPresolve() { return new ClpPresolve(); }

extern "C" void deletePresolve(void *presolve) {
  if (presolve != NULL) delete (ClpPresolve *)presolve;
}

extern "C" void *presolvedModel(void *presolve, void *model) {
  cout << "--------------------------------------------------" << endl;
  cout << "running presolve" << endl;
  cout << "--------------------------------------------------" << endl;
  ClpSimplex *newModel =
      ((ClpPresolve *)presolve)->presolvedModel(*((ClpSimplex *)model));
  return newModel;
}
/*
 * formulate lhs <= Ax <= rhs, l <= x <= u
 * as Ax - z = 0, l <= x <= u, lhs <= z <= rhs
 * do not pre-allocate pointers except model, nCols, nRows, nnz and nEqs
 * set them to NULL is a better practice
 * but do remember to free them
 */

extern "C" int formulateLP(void *model, double **cost, int *nCols, int *nRows,
                           int *nnz, int *nEqs, int **csc_beg, int **csc_idx,
                           double **csc_val, double **rhs, double **lower,
                           double **upper, double *offset, int *nCols_origin) {
  int retcode = 0;

  // problem size for malloc
  int nCols_clp = ((ClpModel *)model)->getNumCols();
  int nRows_clp = ((ClpModel *)model)->getNumRows();
  int nnz_clp = ((ClpModel *)model)->getNumElements();
  *nCols_origin = nCols_clp;
  *nRows = nRows_clp;                                // need not recalculate
  *nCols = nCols_clp;                                // need recalculate
  *nEqs = nRows_clp;                                 // need not recalculate
  *nnz = nnz_clp;                                    // need recalculate
  *offset = ((ClpModel *)model)->objectiveOffset();  // need not recalculate
  
  if (*offset != 0.0) {
    printf("Has obj offset\n");
  } else {
    printf("No obj offset\n");
  }
  // allocate buffer memory
  bool *is_eq = NULL;

  assert(((ClpModel *)model)->matrix()->isColOrdered());

  const double *lhs_clp = ((const ClpModel *)model)->getRowLower();
  const double *rhs_clp = ((const ClpModel *)model)->getRowUpper();
  const CoinBigIndex *A_csc_beg =
      ((const ClpModel *)model)->matrix()->getVectorStarts();
  const int *A_csc_idx = ((const ClpModel *)model)->matrix()->getIndices();
  const double *A_csc_val = ((const ClpModel *)model)->matrix()->getElements();
  int has_lower, has_upper;

  // cupdlp_printf("------------------------------------------------\n");
  // vecIntPrint("A_csc_beg", A_csc_beg, nCols_clp + 1);
  // vecIntPrint("A_csc_idx", A_csc_idx, nnz_clp);
  // vecPrint("A_csc_val", A_csc_val, nnz_clp);

  // cupdlp_printf("%s (Trans):\n", "A_clp");
  // cupdlp_int deltaRow = 0;
  // for (cupdlp_int iCol = 0; iCol < nCols_clp; ++iCol)
  // {
  //     for (cupdlp_int iElem = A_csc_beg[iCol]; iElem < A_csc_beg[iCol + 1];
  //     ++iElem)
  //     {
  //         if (iElem == A_csc_beg[iCol])
  //             deltaRow = A_csc_idx[iElem];
  //         else
  //             deltaRow = A_csc_idx[iElem] - A_csc_idx[iElem - 1] - 1;
  //         for (cupdlp_int i = 0; i < deltaRow; ++i)
  //         {
  //             cupdlp_printf("       ");
  //         }
  //         cupdlp_printf("%6.3f ", A_csc_val[iElem]);
  //     }
  //     cupdlp_printf("\n");
  // }
  // cupdlp_printf("------------------------------------------------\n");

  CUPDLP_INIT(is_eq, nRows_clp);

  // recalculate nRows and nnz for Ax - z = 0
  for (int i = 0; i < nRows_clp; i++) {
    has_lower = lhs_clp[i] > -1e20;
    has_upper = rhs_clp[i] < 1e20;

    // count number of equations and rows
    if (has_lower && has_upper && lhs_clp[i] == rhs_clp[i]) {
      is_eq[i] = true;
    } else {
      is_eq[i] = false;
      (*nCols)++;
      (*nnz)++;
    }
  }

  // allocate memory
  CUPDLP_INIT(*cost, *nCols);
  CUPDLP_INIT(*lower, *nCols);
  CUPDLP_INIT(*upper, *nCols);
  CUPDLP_INIT(*csc_beg, *nCols + 1);
  CUPDLP_INIT(*csc_idx, *nnz);
  CUPDLP_INIT(*csc_val, *nnz);
  CUPDLP_INIT(*rhs, *nRows);

  // formulate LP matrix
  for (int i = 0; i < nCols_clp + 1; i++) (*csc_beg)[i] = A_csc_beg[i];
  for (int i = 0; i < nnz_clp; i++) {
    // can add sort here.
    (*csc_idx)[i] = A_csc_idx[i];
    (*csc_val)[i] = A_csc_val[i];
  }

  for (int i = 0, j = 0; i < nRows_clp; i++) {
    if (!is_eq[i]) {
      (*csc_beg)[nCols_clp + j + 1] = (*csc_beg)[nCols_clp + j] + 1;
      (*csc_idx)[(*csc_beg)[nCols_clp + j]] = i;
      (*csc_val)[(*csc_beg)[nCols_clp + j]] = -1.0;
      j++;
    }
  }

  // rhs
  for (int i = 0; i < *nRows; i++) {
    if (is_eq[i])
      (*rhs)[i] = lhs_clp[i];
    else
      (*rhs)[i] = 0.0;
  }

  // cost, lower, upper
  for (int i = 0; i < nCols_clp; i++) {
    (*cost)[i] = ((ClpModel *)model)->getObjCoefficients()[i] *
                 ((ClpModel *)model)->getObjSense();
    (*lower)[i] = ((ClpModel *)model)->getColLower()[i];

    (*upper)[i] = ((ClpModel *)model)->getColUpper()[i];
  }
  for (int i = nCols_clp; i < *nCols; i++) {
    (*cost)[i] = 0.0;
  }

  for (int i = 0, j = nCols_clp; i < *nRows; i++) {
    if (!is_eq[i]) {
      (*lower)[j] = lhs_clp[i];
      (*upper)[j] = rhs_clp[i];
      j++;
    }
  }

  for (int i = 0; i < *nCols; i++) {
    if ((*lower)[i] < -1e20) (*lower)[i] = -INFINITY;
    if ((*upper)[i] > 1e20) (*upper)[i] = INFINITY;
  }
  //    vecPrint("rhs", *rhs, *nRows);
  //    vecPrint("lower", *lower, *nCols);
  //    vecPrint("upper", *upper, *nCols);

exit_cleanup:
  // free buffer memory
  if (is_eq != NULL) {
    free(is_eq);
    is_eq = NULL;
  }

  return retcode;
}

/*
 * formulate
 *                A x =  b
 *         l1 <= G1 x
 *               G2 x <= u2
 *         l3 <= G3 x <= u3
 * with bounds
 *             l <= x <= u
 * as
 *                A x =  b
 *               G3 x - z = 0
 *               G1 x >= l1
 *              -G2 x >= -u2
 * with bounds
 *             l <= x <= u
 *            l3 <= z <= u3
 * do not pre-allocate pointers except model, nCols, nRows, nnz and nEqs
 * set them to NULL is a better practice
 * but do remember to free them
 */

extern "C" int formulateLP_new(void *model, double **cost, int *nCols,
                               int *nRows, int *nnz, int *nEqs, int **csc_beg,
                               int **csc_idx, double **csc_val, double **rhs,
                               double **lower, double **upper, double *offset,
                               double *sign_origin, int *nCols_origin,
                               int **constraint_new_idx) {
  int retcode = 0;

  // problem size for malloc
  int nCols_clp = ((ClpModel *)model)->getNumCols();
  int nRows_clp = ((ClpModel *)model)->getNumRows();
  int nnz_clp = ((ClpModel *)model)->getNumElements();
  *nCols_origin = nCols_clp;
  *nRows = nRows_clp;                                // need not recalculate
  *nCols = nCols_clp;                                // need recalculate
  *nEqs = 0;                                         // need recalculate
  *nnz = nnz_clp;                                    // need recalculate
  *offset = ((ClpModel *)model)->objectiveOffset();  // need not recalculate
  *sign_origin = ((ClpModel *)model)->getObjSense();

  // printf("offset: %f\n", *offset);
  // printf("nCols_clp: %d, nRows_clp: %d, nnz_clp: %d\n", nCols_clp, nRows_clp,
  //        nnz_clp);
  // printf("nCols: %d, nRows: %d, nnz: %d, nEqs: %d\n", *nCols, *nRows, *nnz,
  //        *nEqs);

  if (*offset != 0.0) {
    printf("Has obj offset %f\n", *offset);
  } else {
    printf("No obj offset\n");
  }
  // allocate buffer memory
  constraint_type *constraint_type_clp = NULL;  // the ONLY one need to free
  // int *constraint_original_idx = NULL;  // pass by user is better, for
  // postsolve recovering dual

  assert(((ClpModel *)model)->matrix()->isColOrdered());

  const double *lhs_clp = ((const ClpModel *)model)->getRowLower();
  const double *rhs_clp = ((const ClpModel *)model)->getRowUpper();
  const CoinBigIndex *A_csc_beg =
      ((const ClpModel *)model)->matrix()->getVectorStarts();
  const int *A_csc_idx = ((const ClpModel *)model)->matrix()->getIndices();
  const double *A_csc_val = ((const ClpModel *)model)->matrix()->getElements();
  int has_lower, has_upper;
  // 打印三个数组
  // for (int i = 0; i < nCols_clp + 1; i++) {
  //   printf("A_csc_beg[%d]: %d\n", i, A_csc_beg[i]);
  // }
  // for (int i = 0; i < nnz_clp; i++) {
  //   printf("A_csc_idx[%d]: %d\n", i, A_csc_idx[i]);
  // }
  // for (int i = 0; i < nnz_clp; i++) {
  //   printf("A_csc_val[%d]: %f\n", i, A_csc_val[i]);
  // }

  CUPDLP_INIT(constraint_type_clp, nRows_clp);
  CUPDLP_INIT(*constraint_new_idx, *nRows);

  // recalculate nRows and nnz for Ax - z = 0
  for (int i = 0; i < nRows_clp; i++) {
    has_lower = lhs_clp[i] > -1e20;
    has_upper = rhs_clp[i] < 1e20;

    // count number of equations and rows
    if (has_lower && has_upper && lhs_clp[i] == rhs_clp[i]) {
      // 既有上界又有下界，且上下界相等
      constraint_type_clp[i] = EQ;
      (*nEqs)++;
    } else if (has_lower && !has_upper) {
      // 只有下界没有上界
      constraint_type_clp[i] = GEQ;
    } else if (!has_lower && has_upper) {
      // 只有上界没有下界
      constraint_type_clp[i] = LEQ;
    } else if (has_lower && has_upper) {
      // 既有上界又有下界，但是二者不相等
      constraint_type_clp[i] = BOUND;
      (*nCols)++;
      (*nnz)++;
      (*nEqs)++;
    } else {
      // printf("Error: constraint %d has no lower and upper bound\n", i);
      // retcode = 1;
      // goto exit_cleanup;

      // what if regard free as bounded
      constraint_type_clp[i] = BOUND;
      (*nCols)++;
      (*nnz)++;
      (*nEqs)++;
    }
  }

  // allocate memory
  CUPDLP_INIT(*cost, *nCols);
  CUPDLP_INIT(*lower, *nCols);
  CUPDLP_INIT(*upper, *nCols);
  CUPDLP_INIT(*csc_beg, *nCols + 1);
  CUPDLP_INIT(*csc_idx, *nnz);
  CUPDLP_INIT(*csc_val, *nnz);
  CUPDLP_INIT(*rhs, *nRows);

  // cost, lower, upper
  for (int i = 0; i < nCols_clp; i++) {
    (*cost)[i] = ((ClpModel *)model)->getObjCoefficients()[i] *
                 ((ClpModel *)model)->getObjSense();
    (*lower)[i] = ((ClpModel *)model)->getColLower()[i];

    (*upper)[i] = ((ClpModel *)model)->getColUpper()[i];
  }
  // slack costs
  for (int i = nCols_clp; i < *nCols; i++) {
    (*cost)[i] = 0.0;
  }
  // slack bounds
  for (int i = 0, j = nCols_clp; i < *nRows; i++) {
    if (constraint_type_clp[i] == BOUND) {
      (*lower)[j] = lhs_clp[i];
      (*upper)[j] = rhs_clp[i];
      j++;
    }
  }

  for (int i = 0; i < *nCols; i++) {
    if ((*lower)[i] < -1e20) (*lower)[i] = -INFINITY;
    if ((*upper)[i] > 1e20) (*upper)[i] = INFINITY;
  }

  // permute LP rhs
  // EQ or BOUND first
  for (int i = 0, j = 0; i < *nRows; i++) {
    if (constraint_type_clp[i] == EQ) {
      (*rhs)[j] = lhs_clp[i];
      (*constraint_new_idx)[i] = j;
      j++;
    } else if (constraint_type_clp[i] == BOUND) {
      (*rhs)[j] = 0.0;
      (*constraint_new_idx)[i] = j;
      j++;
    }
  }
  // then LEQ or GEQ
  for (int i = 0, j = *nEqs; i < *nRows; i++) {
    if (constraint_type_clp[i] == LEQ) {
      (*rhs)[j] = -rhs_clp[i];  // multiply -1
      (*constraint_new_idx)[i] = j;
      j++;
    } else if (constraint_type_clp[i] == GEQ) {
      (*rhs)[j] = lhs_clp[i];
      (*constraint_new_idx)[i] = j;
      j++;
    }
  }

  // formulate and permute LP matrix
  // beg remains the same
  for (int i = 0; i < nCols_clp + 1; i++) (*csc_beg)[i] = A_csc_beg[i];
  for (int i = nCols_clp + 1; i < *nCols + 1; i++)
    (*csc_beg)[i] = (*csc_beg)[i - 1] + 1;

  // row idx changes
  for (int i = 0, k = 0; i < nCols_clp; i++) {
    // same order as in rhs
    // EQ or BOUND first
    for (int j = (*csc_beg)[i]; j < (*csc_beg)[i + 1]; j++) {
      if (constraint_type_clp[A_csc_idx[j]] == EQ ||
          constraint_type_clp[A_csc_idx[j]] == BOUND) {
        (*csc_idx)[k] = (*constraint_new_idx)[A_csc_idx[j]];
        (*csc_val)[k] = A_csc_val[j];
        k++;
      }
    }
    // then LEQ or GEQ
    for (int j = (*csc_beg)[i]; j < (*csc_beg)[i + 1]; j++) {
      if (constraint_type_clp[A_csc_idx[j]] == LEQ) {
        (*csc_idx)[k] = (*constraint_new_idx)[A_csc_idx[j]];
        (*csc_val)[k] = -A_csc_val[j];  // multiply -1
        k++;
      } else if (constraint_type_clp[A_csc_idx[j]] == GEQ) {
        (*csc_idx)[k] = (*constraint_new_idx)[A_csc_idx[j]];
        (*csc_val)[k] = A_csc_val[j];
        k++;
      }
    }
  }

  // slacks for BOUND
  for (int i = 0, j = nCols_clp; i < *nRows; i++) {
    if (constraint_type_clp[i] == BOUND) {
      (*csc_idx)[(*csc_beg)[j]] = (*constraint_new_idx)[i];
      (*csc_val)[(*csc_beg)[j]] = -1.0;
      j++;
    }
  }
  // // 打印三个数组
  // for (int i = 0; i < *nCols + 1; i++) {
  //   printf("csc_beg[%d]: %d\n", i, (*csc_beg)[i]);
  // }
  // for (int i = 0; i < *nnz; i++) {
  //   printf("csc_idx[%d]: %d\n", i, (*csc_idx)[i]);
  // }
  // for (int i = 0; i < *nnz; i++) {
  //   printf("csc_val[%d]: %f\n", i, (*csc_val)[i]);
  // }
  // for (int i = 0; i < *nCols; i++){
  //   printf("cost[%d]: %f\n", i, (*cost)[i]);
  // }
  // for (int i = 0; i < *nRows; i++){
  //   printf("rhs[%d]: %f\n", i, (*rhs)[i]);
  // }
  // for (int i = 0; i < *nCols; i++){
  //   printf("lower[%d]: %f\n", i, (*lower)[i]);
  // }
  // for (int i = 0; i < *nCols; i++){
  //   printf("upper[%d]: %f\n", i, (*upper)[i]);
  // }
  exit_cleanup:
    // free buffer memory
    if (constraint_type_clp != NULL)
    {
      free(constraint_type_clp);
      constraint_type_clp = NULL;
    }

  return retcode;
}

// #endif

extern "C" void Construct_dualOT_Matrix(void* matrix, int resolution){
  printf("Construct_dualOT_Matrix\n");
  CoinPackedMatrix *matrix_C = static_cast<CoinPackedMatrix*>(matrix);

  // true是按列添加数据，false是按行添加数据
  bool colordered = false;
  CoinPackedMatrix mat = CoinPackedMatrix(colordered, 0, 0);
  #pragma region example
  // 由于我们是按列添加数据，indices表示行索引
  // std::vector<int> indices;
  // indices.push_back(0); // 第一个约束的索引
  // indices.push_back(1); // 第二个约束的索引
  int *indices = NULL;
  int indices_len = 2;
  indices = (int *)malloc(indices_len * sizeof(int));
  indices[0] = 0;
  indices[1] = 1;

  // 定义两个变量在每个约束中的系数
  // std::vector<double> elements_x; // x的系数
  // elements_x.push_back(2); // 第一个约束中x的系数
  // elements_x.push_back(1); // 第二个约束中x的系数
  double *elements_x = NULL;
  elements_x = (double *)malloc(indices_len * sizeof(double));
  elements_x[0] = 2;
  elements_x[1] = 1;
  double *elements_y = NULL;
  elements_y = (double *)malloc(indices_len * sizeof(double));
  elements_y[0] = 3;
  elements_y[1] = 6;


  // std::vector<double> elements_y; // y的系数
  // elements_y.push_back(1); // 第一个约束中y的系数
  // elements_y.push_back(3); // 第二个约束中y的系数

  // // 添加第一列（x的系数）
  // mat.appendCol(indices_len, &indices[0], &elements_x[0]);

  // // 添加第二列（y的系数）
  // mat.appendCol(indices_len, &indices[0], &elements_y[0]);
  mat.appendRow(indices_len, indices, elements_x);
  mat.appendRow(indices_len, indices, elements_y);
  #pragma endregion
  
  // int nCols = 2 * pow(resolution,2);
  // int nRows = pow(resolution, 4);
  // printf("nCols: %d, nRows: %d\n", nCols, nRows);
  // int vec_len = pow(resolution, 2);
  
  
  // // i代表列数，j代表行数
  // // 前一半的列
  // for (int i = 0; i < nCols /2 ; i++){
  //   std::vector<int> indices;
  //   for (int j = i * vec_len; j < (i + 1) * vec_len; j++){
  //     indices.push_back(j);
  //   }
  //   std::vector<double> elements;
  //   for (int i = 0; i < indices.size(); i++){
  //     elements.push_back(1);
  //   }
  //     mat.appendCol(indices.size(), &indices[0], &elements[0]);
  // }
  // // 后一半的列
  // for (int i = 0; i < nCols /2 ; i++){
  //   std::vector<int> indices;
  //   for (int j = 0; j < vec_len; j++){
  //     indices.push_back(i + j * vec_len);
  //   }
  //   std::vector<double> elements;
  //   for (int i = 0; i < indices.size(); i++){
  //     elements.push_back(1);
  //   }
  //   mat.appendCol(indices.size(), &indices[0], &elements[0]);
  // }
  //   printf("Construct_dualOT_Matrix end\n");
#pragma region print
  // 列主元时的打印函数
  const int numRows = mat.getNumRows();
  const int numCols = mat.getNumCols();
  printf("numRows: %d, numCols: %d\n", numRows, numCols);

  // 获取矩阵中的元素和它们的索引
  const double* elements_print = mat.getElements();
  const int* rowIndices = mat.getIndices();
  const CoinBigIndex* columnStarts = mat.getVectorStarts();
  const int* lengths = mat.getVectorLengths();

  // 按列遍历矩阵
  for (int col = 0; col < numCols; ++col) {
      std::cout << "Column " << col << ":" << std::endl;
      // 获取这一列的起始位置和长度
      CoinBigIndex start = columnStarts[col];
      int length = lengths[col];

      // 遍历这一列的每个元素
      for (int i = 0; i < length; ++i) {
          // 计算当前元素的索引
          CoinBigIndex index = start + i;
          // 输出行索引和元素的值
          std::cout << "  Row " << rowIndices[index] << ": " << elements_print[index] << std::endl;
      }
  }
#pragma endregion
  
  int *zero_idx = NULL;
  int zero_idx_len = 2;
  zero_idx = (int *)malloc(zero_idx_len * sizeof(int));
  zero_idx[0] = 0;
  zero_idx[1] = 1;
  mat.deleteCols(zero_idx_len, zero_idx);
  printf("删除后的矩阵维度：nRows: %d, nCols: %d\n", mat.getNumRows(), mat.getNumCols());
  

  matrix_C = &mat;
}

extern "C" void countZero_and_CheckConstraint_Keep_Wrapper(long long **keep, long long *keep_nnz, double *y, double *x, int resolution_y, double thr, double violate_degree){
  printf("countZero_and_CheckConstraint_Keep_Wrapper开始\n");
  int num_threads = 128;
  if (num_threads > resolution_y)
  {
      num_threads = resolution_y;
  }
  std::vector<long long> nnz_array(num_threads);

  std::vector<std::vector<long long>> keep_array;
  for (int i = 0; i < num_threads; i++){
    std::vector<long long> keep_sub;
    keep_array.push_back(keep_sub);
  }
  
  int pow_resolution_2 = pow(resolution_y, 2);
  double scale = 2.0 * pow_resolution_2;
printf("countZero_and_CheckConstraint_Keep_Wrapper omp parallel开始\n");
#pragma omp parallel num_threads(num_threads)
{
  int thread_id = omp_get_thread_num();
  int start = thread_id * resolution_y / num_threads;
  int end = (thread_id + 1) * resolution_y / num_threads;
  nnz_array[thread_id] = 0;
  for (int i1 = start; i1 < end; i1++)
  {
    long long idx_temp = 0;
    long long idx_temp_min = 0;
    long long idx_1 = 0;
    long long idx_2 = 0;
    for (int j1 = 0; j1 < resolution_y; j1++)
    {
      idx_1 = i1 * resolution_y + j1;
      for (int i2 = 0; i2 < resolution_y; i2++)
      {
        for (int j2 = 0; j2 < resolution_y; j2++)
        {
          idx_2 = i2 * resolution_y + j2;
          idx_temp = idx_1 * pow_resolution_2 + idx_2;
          idx_temp_min = start * pow_resolution_2 + idx_2;
          idx_temp_min = idx_temp_min * resolution_y;
          if (fabs(y[idx_temp]) < thr)
          {
            if (x[idx_1] + x[pow_resolution_2 + idx_2] > (1+violate_degree) * ((i1 - i2) * (i1 - i2) + (j1 - j2) * (j1 - j2) + 1e-8) / scale)
            {
              keep_array[thread_id].push_back(idx_temp);
              nnz_array[thread_id] += 1;
            }
          }
          else
          {
            keep_array[thread_id].push_back(idx_temp);
            nnz_array[thread_id] += 1;
          }
        }
      }
    }
  }
}
  // 先前的nnz_array是每个线程的非零元素个数，现在需要计算到每个线程为止有多少个非零元素
  for (int i = 1; i < num_threads; i++)
  {
    nnz_array[i] += nnz_array[i - 1];
  }
  *keep_nnz = nnz_array[num_threads - 1];
  // for (int i = 0; i < num_threads; i++){
  //   printf("nnz_array[%d]: %lld\n", i, nnz_array[i]);
  // }
  // for (int i = 0; i < num_threads; i++){
  //   for (int j = 0; j < keep_array[i].size(); j++){
  //     printf("keep_array[%d][%d]: %lld\n", i, j, keep_array[i][j]);
  //   }
  // }
  // 再合并到keep中
  // keep的第i个元素值为idx，这对应着y中的第idx个元素是第i个被保留的元素
  // bool *keep = (bool *)malloc(sizeof(bool) * nRows);
  printf("countZero_and_CheckConstraint_Keep_Wrapper合并开始\n");
  printf("keep_nnz: %lld\n", *keep_nnz);
  *keep = (long long *)malloc(sizeof(long long) * *keep_nnz);
  if (*keep == NULL)
  {
    printf("keep malloc失败\n");
  }
  else
  {
    printf("keep malloc成功\n");
  }
#pragma omp parallel num_threads(num_threads)
{
  int thread_id = omp_get_thread_num();
  if (thread_id == 0){
    for (int i = 0; i < keep_array[thread_id].size(); i++)
    {
      (*keep)[i] = keep_array[thread_id][i];
    }
  }
  else
  {
    for (int i = 0; i < keep_array[thread_id].size(); i++)
    {
      (*keep)[nnz_array[thread_id-1] + i] = keep_array[thread_id][i];
    }
  }
  
}
  // for (int i = 0; i < *keep_nnz; i++)
  // {
  // printf("keep[%d]: %lld\n", i, (*keep)[i]);
  // }
  printf("countZero_and_CheckConstraint_Keep_Wrapper结束\n");
}

extern "C" void countZero_Keep_Wrapper(long long **keep, long long *keep_nnz, double *y, int resolution_y, double thr, double violate_degree){
  printf("countZero_Keep_Wrapper开始\n");
  int num_threads = 128;
  if (num_threads > resolution_y)
  {
      num_threads = resolution_y;
  }
  std::vector<long long> nnz_array(num_threads);

  std::vector<std::vector<long long>> keep_array;
  for (int i = 0; i < num_threads; i++){
    std::vector<long long> keep_sub;
    keep_array.push_back(keep_sub);
  }
  
  int pow_resolution_2 = pow(resolution_y, 2);
  double scale = 2.0 * pow_resolution_2;
printf("countZero_Keep_Wrapper omp parallel开始\n");
#pragma omp parallel num_threads(num_threads)
{
  int thread_id = omp_get_thread_num();
  int start = thread_id * resolution_y / num_threads;
  int end = (thread_id + 1) * resolution_y / num_threads;
  nnz_array[thread_id] = 0;
  for (int i1 = start; i1 < end; i1++)
  {
    long long idx_temp = 0;
    long long idx_temp_min = 0;
    long long idx_1 = 0;
    long long idx_2 = 0;
    for (int j1 = 0; j1 < resolution_y; j1++)
    {
      idx_1 = i1 * resolution_y + j1;
      for (int i2 = 0; i2 < resolution_y; i2++)
      {
        for (int j2 = 0; j2 < resolution_y; j2++)
        {
          idx_2 = i2 * resolution_y + j2;
          idx_temp = idx_1 * pow_resolution_2 + idx_2;
          idx_temp_min = start * pow_resolution_2 + idx_2;
          idx_temp_min = idx_temp_min * resolution_y;
          // if (fabs(y[idx_temp]) < thr)
          // {
          //   if (x[idx_1] + x[pow_resolution_2 + idx_2] > (1+violate_degree) * ((i1 - i2) * (i1 - i2) + (j1 - j2) * (j1 - j2) + 1e-8) / scale)
          //   {
          //     keep_array[thread_id].push_back(idx_temp);
          //     nnz_array[thread_id] += 1;
          //   }
          // }
          // else
          // {
          //   keep_array[thread_id].push_back(idx_temp);
          //   nnz_array[thread_id] += 1;
          // }
          if (fabs(y[idx_temp] > thr))
          {
            keep_array[thread_id].push_back(idx_temp);
            nnz_array[thread_id] += 1;
          }
        }
      }
    }
  }
}
  // 先前的nnz_array是每个线程的非零元素个数，现在需要计算到每个线程为止有多少个非零元素
  for (int i = 1; i < num_threads; i++)
  {
    nnz_array[i] += nnz_array[i - 1];
  }
  *keep_nnz = nnz_array[num_threads - 1];
  // for (int i = 0; i < num_threads; i++){
  //   printf("nnz_array[%d]: %lld\n", i, nnz_array[i]);
  // }
  // for (int i = 0; i < num_threads; i++){
  //   for (int j = 0; j < keep_array[i].size(); j++){
  //     printf("keep_array[%d][%d]: %lld\n", i, j, keep_array[i][j]);
  //   }
  // }
  // 再合并到keep中
  // keep的第i个元素值为idx，这对应着y中的第idx个元素是第i个被保留的元素
  // bool *keep = (bool *)malloc(sizeof(bool) * nRows);
  printf("countZero_Keep_Wrapper合并开始\n");
  printf("keep_nnz: %lld\n", *keep_nnz);
  *keep = (long long *)malloc(sizeof(long long) * *keep_nnz);
  if (*keep == NULL)
  {
    printf("keep malloc失败\n");
  }
  else
  {
    printf("keep malloc成功\n");
  }
#pragma omp parallel num_threads(num_threads)
{
  int thread_id = omp_get_thread_num();
  if (thread_id == 0){
    for (int i = 0; i < keep_array[thread_id].size(); i++)
    {
      (*keep)[i] = keep_array[thread_id][i];
    }
  }
  else
  {
    for (int i = 0; i < keep_array[thread_id].size(); i++)
    {
      (*keep)[nnz_array[thread_id-1] + i] = keep_array[thread_id][i];
    }
  }
  
}
  // for (int i = 0; i < *keep_nnz; i++)
  // {
  // printf("keep[%d]: %lld\n", i, (*keep)[i]);
  // }
  printf("countZero_Keep_Wrapper结束\n");
}

extern "C" void countZero_and_checkConstraint_Keep_redundancy_Wrapper(long long **keep_fine_redundancy, long long *keep_fine_redundancy_len, double *y_solution_last, long long y_solution_last_len, double *x_init, int resolution_now, int resolution_last,double thr, double violate_degree){
  std::vector<long long> keep_nonzero_constraint_vec;
  int scale = resolution_now / resolution_last;
  long long pow_resolution_now_4 = pow(resolution_now, 4);
  long long pow_resolution_now_2 = pow(resolution_now, 2);
  long long pow_resolution_last_2 = pow(resolution_last, 2);
  double scale_constant = 2.0 * pow_resolution_now_2;
  *keep_fine_redundancy_len = 0;

  #pragma omp parallel
  {
    std::vector<long long> keep_nonzero_constraint_vec_local;
    long long local_count = 0;  // 定义局部变量来跟踪每个线程的计数
    #pragma omp for nowait
    for (long long i1 = 0; i1 < resolution_last; i1++) {
      for (long long i2 = 0; i2 < resolution_last; i2++){
        for (long long j1 = 0; j1 < resolution_last; j1++){
          for(long long j2 = 0; j2 < resolution_last; j2++){
            for(long long k1 = 0; k1 < scale; k1++){
              for (long long k2 = 0; k2 < scale; k2++){
                for (long long l1 = 0; l1 < scale; l1++){
                  for (long long l2 = 0; l2 < scale; l2++){
                    long long idx_coarse = (i1 * resolution_last + j1) * pow_resolution_last_2 + i2 * resolution_last + j2;
                    long long idx_i1 = i1 * scale + k1;
                    long long idx_i2 = i2 * scale + k2;
                    long long idx_j1 = j1 * scale + l1;
                    long long idx_j2 = j2 * scale + l2;
                    long long idx_1 = idx_i1 * resolution_now + idx_j1;
                    long long idx_2 = idx_i2 * resolution_now + idx_j2;
                    if (fabs(y_solution_last[idx_coarse]) >= thr)
                    {
                      long long idx_fine = idx_1 * pow_resolution_now_2 + idx_2;
                      keep_nonzero_constraint_vec_local.push_back(idx_fine);
                      local_count++;
                    }
                    else{
                      if(x_init[idx_1] + x_init[pow_resolution_now_2 + idx_2] > (1+violate_degree) * ((idx_i1 - idx_i2) * (idx_i1 - idx_i2) + (idx_j1 - idx_j2) * (idx_j1 - idx_j2)) / scale_constant){
                        long long idx_fine = idx_1 * pow_resolution_now_2 + idx_2;
                        keep_nonzero_constraint_vec_local.push_back(idx_fine);
                        local_count++;                        
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    
    // 合并局部向量到全局向量
    #pragma omp critical
    {
      keep_nonzero_constraint_vec.insert(keep_nonzero_constraint_vec.end(), keep_nonzero_constraint_vec_local.begin(), keep_nonzero_constraint_vec_local.end());
      *keep_fine_redundancy_len += local_count;  // 使用原子操作或临界区来更新全局长度
    }
  }
  *keep_fine_redundancy = (long long *)malloc(sizeof(long long) * *keep_fine_redundancy_len);
  #pragma omp parallel for
  for (long long i = 0; i < *keep_fine_redundancy_len; i++){
    (*keep_fine_redundancy)[i] = keep_nonzero_constraint_vec[i];
  }
  printf("countZero_and_checkConstraint_Keep_redundancy_Wrapper, keep_fine_redundancy_len: %lld\n", *keep_fine_redundancy_len);
}

