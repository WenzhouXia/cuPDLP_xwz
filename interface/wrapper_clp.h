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
  void writeLpWrapper(void *model, const char *filename);

  void deleteModel(void *model);
  void *createModel();

  void *createPresolve();
  void deletePresolve(void *presolve);
  void *presolvedModel(void *presolve, void *model);
  void Construct_dualOT_Matrix(void *matrix, int resolution);

#ifdef __cplusplus
}
#endif
#endif // CUPDLP_WRAPPER_CLP_H
