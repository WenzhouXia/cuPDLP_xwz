#include "wrapper_clp.h"

// #ifdef __cplusplus
#include <ClpPresolve.hpp>

#include "ClpModel.hpp"
#include "ClpSimplex.hpp"

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

  CUPDLP_INIT(constraint_type_clp, nRows_clp);
  CUPDLP_INIT(*constraint_new_idx, *nRows);

  // recalculate nRows and nnz for Ax - z = 0
  for (int i = 0; i < nRows_clp; i++) {
    has_lower = lhs_clp[i] > -1e20;
    has_upper = rhs_clp[i] < 1e20;

    // count number of equations and rows
    if (has_lower && has_upper && lhs_clp[i] == rhs_clp[i]) {
      constraint_type_clp[i] = EQ;
      (*nEqs)++;
    } else if (has_lower && !has_upper) {
      constraint_type_clp[i] = GEQ;
    } else if (!has_lower && has_upper) {
      constraint_type_clp[i] = LEQ;
    } else if (has_lower && has_upper) {
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

exit_cleanup:
  // free buffer memory
  if (constraint_type_clp != NULL) {
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

