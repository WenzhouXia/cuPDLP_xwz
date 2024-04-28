//
// Created by chuwen on 23-12-4.
//

#ifndef CUPDLP_WRAPPER_CLP_H
#define CUPDLP_WRAPPER_CLP_H

#ifdef __cplusplus
extern "C"
{
#endif

#define CUPDLP_INIT(var, size)                                  \
  {                                                             \
    (var) = (typeof(var))malloc((size) * sizeof(typeof(*var))); \
    if ((var) == NULL)                                          \
    {                                                           \
      retcode = 1;                                              \
      goto exit_cleanup;                                        \
    }                                                           \
  }

  typedef enum CONSTRAINT_TYPE
  {
    EQ = 0,
    LEQ,
    GEQ,
    BOUND
  } constraint_type;

  int formulateLP(void *model, double **cost, int *nCols, int *nRows, int *nnz,
                  int *nEqs, int **csc_beg, int **csc_idx, double **csc_val,
                  double **rhs, double **lower, double **upper, double *offset,
                  int *nCols_origin);

  int formulateLP_new(void *model, double **cost, int *nCols, int *nRows,
                      int *nnz, int *nEqs, int **csc_beg, int **csc_idx,
                      double **csc_val, double **rhs, double **lower,
                      double **upper, double *offset, double *sign_origin,
                      int *nCols_origin, int **constraint_new_idx);

  void loadMps(void *model, const char *filename);

  void loadProblemWrapper(void *model, int numCols, int numRows,
                          const int *start, const int *index,
                          const double *value, const double *colLower,
                          const double *colUpper, const double *obj,
                          const double *rowLower, const double *rowUpper);
  void loadProblem_byMatrix_Wrapper(void *model, int resolution,
                                    const double *colLower,
                                    const double *colUpper, const double *obj,
                                    const double *rowLower, const double *rowUpper);
  void loadProblem_delete_byMatrix_Wrapper(void *model, int resolution, int *zero_idx, int *zero_idx_len,
                                           const double *colLower,
                                           const double *colUpper, const double *obj,
                                           const double *rowLower, const double *rowUpper);
  void loadProblem_delete_byMatrix_Wrapper_longlong(void *model, int resolution, long long *zero_idx, long long *zero_idx_len,
                                                    const double *colLower,
                                                    const double *colUpper, const double *obj,
                                                    const double *rowLower, const double *rowUpper);
  void loadProblem_delete_byMatrix_Wrapper_longlong_parallel(void *model, int resolution, long long *zero_idx, long long *zero_idx_len,
                                                             const double *colLower,
                                                             const double *colUpper, const double *obj,
                                                             const double *rowLower, const double *rowUpper);
  void loadProblem_delete_byMatrix_byKeepIdx_Wrapper_longlong_parallel(void *model, int resolution, long long *keep_idx, long long *keep_nnz,
                                                                       const double *colLower,
                                                                       const double *colUpper, const double *obj,
                                                                       const double *rowLower, const double *rowUpper);
  void generate_ConstraintMatrix_byKeepIdx_Wrapper_longlong_parallel(void *model, int resolution, long long *keep_idx, long long *keep_nnz);
  void printSparseMatrix(void *mat_ptr);
  void loadProblem_delete_by_keep_byMatrix_Wrapper_longlong(void *model, int resolution, bool *keep, long long *keep_true_idx, long long *len_after_delete,
                                                            const double *colLower,
                                                            const double *colUpper, const double *obj,
                                                            const double *rowLower, const double *rowUpper);
  void writeLpWrapper(void *model, const char *filename);

  void deleteModel(void *model);
  void *createModel();

  void *createPresolve();
  void deletePresolve(void *presolve);
  void *presolvedModel(void *presolve, void *model);
  void Construct_dualOT_Matrix(void *matrix, int resolution);
  void countZero_and_CheckConstraint_Keep_Wrapper(long long **keep, long long *keep_nnz, double *y, double *x, int resolution_y, double thr, double violate_degree);
  void countZero_Keep_Wrapper(long long **keep, long long *keep_nnz, double *y, int resolution_y, double thr, double violate_degree);
  void countZero_and_checkConstraint_Keep_redundancy_Wrapper(long long **keep_fine_redundancy, long long *keep_fine_redundancy_len, double *y_solution_last, long long y_solution_last_len, double *x_init, int resolution_now, int resolution_last, double thr, double violate_degree);
  void countZero_Keep_redundancy_Wrapper(long long **keep_countZero, long long *keep_countZero_len, double *y_solution_last, long long y_solution_last_len, int resolution_now, int resolution_last, double thr);
  void countZero_Keep_redundancy_keepLast_Wrapper(long long **keep_countZero, long long *keep_countZero_len, double *y_solution_last, long long y_solution_last_len, long long *keep_last, long long keep_last_len, int resolution_now, int resolution_last, double thr);
#ifdef __cplusplus
}
#endif
#endif // CUPDLP_WRAPPER_CLP_H
