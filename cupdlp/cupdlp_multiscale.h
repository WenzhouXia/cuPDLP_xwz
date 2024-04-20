//
// created by Wenzhou Xia
//

#ifndef CUPDLP_CUPDLP_MULTISCALE_H
#define CUPDLP_CUPDLP_MULTISCALE_H

#include "cupdlp_defs.h"
#include "cupdlp_linalg.h"
#include "cupdlp_proj.h"
// #include "cupdlp_scaling.h"
#include "cupdlp_step.h"
#include "cupdlp_utils.h"
#include "glbopts.h"

void CUPDLP_multiscale_testprint();

void readCSVToFloatArray(cupdlp_float *data, const char *filename, const int numCols);
void normalizeArray1D(cupdlp_float *array, cupdlp_int array_len);
void coarsingArray1D(cupdlp_float *array_coarse, cupdlp_float *array, cupdlp_int resolution, const cupdlp_int coarse_degree);
void print_float_array2D(cupdlp_float **array, cupdlp_int numRows, cupdlp_int numCols);
void print_float_array1D(cupdlp_float *array, cupdlp_int num);
void print_int_array1D(cupdlp_int *array, cupdlp_int num);
cupdlp_float **mergeTwoArrays2D(cupdlp_float **a, cupdlp_float **b, cupdlp_int resolution);
cupdlp_float *mergeTwoArrays1D(const cupdlp_float *a, const cupdlp_float *b, int a_len, int b_len);
cupdlp_float *mergeTwoArrays1D_minus(const cupdlp_float *a, const cupdlp_float *b, cupdlp_int a_len, cupdlp_int b_len);
void mergeTwoArrays1D_minus_ptr(cupdlp_float *result, const cupdlp_float *a, const cupdlp_float *b, cupdlp_int a_len, cupdlp_int b_len);
void normalizedSquaredEuclideanDistance(cupdlp_float *c, cupdlp_int m, cupdlp_int n);
void normalizedSquaredEuclideanDistance_longlong(cupdlp_float *c, cupdlp_int m, cupdlp_int n);
void coarsing_normalizedSquaredEuclideanDistance(cupdlp_float *c_coarse, cupdlp_int m, cupdlp_int n, const cupdlp_int coarse_degree);
void coarsing_normalizedSquaredEuclideanDistance_longlong(cupdlp_float *c_coarse, cupdlp_int m, cupdlp_int n, const cupdlp_int coarse_degree);
void generate_c_coarse_delete_directly(cupdlp_float *c_coarse_delete, cupdlp_int coarse_degree, cupdlp_int m, cupdlp_int n, cupdlp_bool *keep, long long *keep_idx);
void generate_c_coarse_delete_directly_byKeepIdx_parallel(cupdlp_float *c_coarse_delete, cupdlp_int coarse_degree, cupdlp_int m, cupdlp_int n, long long *keep_idx);
void generate_minus_c_coarse_delete_directly_byKeepIdx_parallel(cupdlp_float *c_coarse_delete, cupdlp_int coarse_degree, cupdlp_int m, cupdlp_int n, long long *keep_idx);
void generate_c_coarse_delete_directly_byKeepIdx_parallel_pdbalance(cupdlp_float *c_coarse_delete, cupdlp_int coarse_degree, cupdlp_int m, cupdlp_int n, long long *keep_idx, long long *keep_nnz, cupdlp_float *balance_weight);
void generate_c_coarse_directly_parallel(cupdlp_float *c_coarse, cupdlp_int coarse_degree, cupdlp_int m, cupdlp_int n);
cupdlp_int *dualOT_startArray(cupdlp_int m, cupdlp_int n);
cupdlp_int *dualOT_indexArray(cupdlp_int m, cupdlp_int n);
cupdlp_float *dualOT_valueArray(cupdlp_int m, cupdlp_int n);
cupdlp_float *dualOT_rowLower(cupdlp_int m, cupdlp_int n);
cupdlp_float *dualOT_rowLower_delete(cupdlp_int len_after_delete);
// cupdlp_float *dualOT_rowUpper(cupdlp_int m, cupdlp_int n);
cupdlp_float *dualOT_colLower(cupdlp_int m, cupdlp_int n);
cupdlp_float *dualOT_colUpper(cupdlp_int m, cupdlp_int n);
cupdlp_float *dualOT_colLower_inf(cupdlp_int m, cupdlp_int n);
cupdlp_float *dualOT_colUpper_inf(cupdlp_int m, cupdlp_int n);
void generate_dualOT_model_from_distribution_and_cost(void *model, const cupdlp_float *a, const cupdlp_float *b, const cupdlp_int a_len, const cupdlp_int b_len, cupdlp_float *c);
void generate_dualOT_model_byMatrix_from_distribution_and_cost(void *model, const cupdlp_float *a, const cupdlp_float *b, const cupdlp_int a_len, const cupdlp_int b_len, cupdlp_float *c, cupdlp_int resolution);
void generate_dualOT_model_delete_byMatrix_from_distribution_and_cost(void *model, const cupdlp_float *a, const cupdlp_float *b, const cupdlp_int a_len, const cupdlp_int b_len, cupdlp_float *c_delete, int resolution, int *zero_idx, int *zero_idx_len);
void generate_dualOT_model_delete_byMatrix_from_distribution_and_cost_longlong(void *model, const cupdlp_float *a, const cupdlp_float *b, const long long a_len, const long long b_len, cupdlp_float *c_delete, int resolution, long long *zero_idx, long long *zero_idx_len);
void generate_dualOT_model_delete_byMatrix_from_distribution_and_cost_longlong_parallel(void *model, const cupdlp_float *a, const cupdlp_float *b, const long long a_len, const long long b_len, cupdlp_float *c_delete, int resolution, long long *zero_idx, long long *zero_idx_len);
void generate_dualOT_model_delete_byMatrix_byKeepIdx_from_distribution_and_cost_longlong_parallel(void *model, const cupdlp_float *a, const cupdlp_float *b, const long long a_len, const long long b_len, cupdlp_float *c_delete, int resolution, long long *keep_idx, long long *keep_nnz);
void generate_dualOT_model_delete_by_keep_byMatrix_from_distribution_and_cost_longlong(void *model, const cupdlp_float *a, const cupdlp_float *b, const long long a_len, const long long b_len, cupdlp_float *c_delete, int resolution, cupdlp_bool *keep, long long *keep_idx, long long *len_after_delete);
void generate_dualOT_model_from_csv(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution);
void generate_coarse_dualOT_model_from_csv(void *model_coarse, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, const cupdlp_int coarse_degree);
void generate_dualOT_model_delete_from_csv(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_int *zero_idx, cupdlp_int *zero_idx_len, cupdlp_int len_after_delete);
void generate_coarse_dualOT_model_delete_from_csv(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_int *zero_idx, cupdlp_int *zero_idx_len, cupdlp_int len_after_delete, cupdlp_int coarse_degree);
void generate_coarse_dualOT_model_delete_from_csv_longlong(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, long long *zero_idx, long long *zero_idx_len, long long len_after_delete, cupdlp_int coarse_degree);
void generate_coarse_dualOT_model_delete_from_csv_longlong_parallel(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, long long *zero_idx, long long *zero_idx_len, long long len_after_delete, cupdlp_int coarse_degree);
void generate_coarse_dualOT_model_delete_byKeepIdx_from_csv_longlong_parallel(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, long long *keep_idx, long long *keep_nnz, cupdlp_int coarse_degree);
void generate_coarse_dualOT_model_delete_byKeepIdx_from_csv_longlong_parallel_pdbalance(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, long long *keep_idx, long long *keep_nnz, cupdlp_int coarse_degree, cupdlp_float *balance_weight);
// void LP_Solve_Multiscale(w, w_coarse, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx);
void createCUPDLPwork(CUPDLPwork *w, void *model, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_int **constraint_new_idx, int *nCols_origin_ptr, int *nRows_ptr);
void createCUPDLPwork_clear(CUPDLPwork *w, void *model, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_int **constraint_new_idx);
void fine_dualOT_primal(cupdlp_float *x_init, cupdlp_float *x_coarse_solution, cupdlp_int x_len, cupdlp_int x_coarse_len, cupdlp_int resolution, cupdlp_int coarse_degree);
void fine_dualOT_dual(cupdlp_float *y_init, cupdlp_float *y_coarse_solution, long long y_len, long long y_coarse_len, cupdlp_int resolution, cupdlp_int coarse_degree);
void fine_dualOT_dual_parallel(cupdlp_float *y_init, cupdlp_float *y_coarse_solution, long long y_len, long long y_coarse_len, cupdlp_int resolution, cupdlp_int coarse_degree);
cupdlp_float *fine_and_delete_dualOT_dual_longlong(long long *len_after_delete, cupdlp_bool *keep, long long *keep_idx, cupdlp_float *y_coarse_solution, cupdlp_int resolution, cupdlp_int coarse_degree_now, cupdlp_int coarse_degree_last, cupdlp_float thr);
int countArray1D_Smaller_than_threshold(cupdlp_float *a, int a_len, cupdlp_float thr);
long long countArray1D_Smaller_than_threshold_longlong(cupdlp_float *a, long long a_len, cupdlp_float thr);
int *countArray1D_Smaller_than_threshold_with_Record(cupdlp_float *a, int a_len, int *a_record_len, cupdlp_float thr);
long long *countArray1D_Smaller_than_threshold_with_Record_longlong(cupdlp_float *a, long long a_len, long long *a_record_len, cupdlp_float thr);
void countArray1D_Smaller_than_threshold_with_KeepIdx_longlong_parallel(long long *keep_idx, cupdlp_float *a, long long a_len, long long *keep_nnz, cupdlp_float thr);
void countZero_and_checkConstraint_with_KeepIdx_longlong_parallel(long long *keep_idx, cupdlp_float *y, long long y_len, cupdlp_float *x, long long x_len, long long *keep_nnz, cupdlp_float thr, cupdlp_int resolution_now, cupdlp_float violate_degree);
void saveArray1D_to_csv(cupdlp_float *a, int a_len, const char *filename);
void analyseArray1D(cupdlp_float *a, long long a_len, cupdlp_float thr, const char *filename);
int countArray1D_same_elements(int *a, int a_len, int *b, int b_len);
void compareTwoArray1D(cupdlp_float *a, cupdlp_int a_len, cupdlp_float *b, cupdlp_int b_len, cupdlp_float thr);
void constructCoinPackedMatrix(void *mat, cupdlp_int resolution);
void deleteArray1DElements(cupdlp_float *a_dele, cupdlp_float *a, cupdlp_int a_len, cupdlp_int *zero_idx, cupdlp_int *zero_idx_len);
void deleteArray1DElements_longlong(cupdlp_float *a_dele, cupdlp_float *a, long long a_len, long long *zero_idx, long long *zero_idx_len);
void deleteArray1DElements_byKeepIdx_longlong_parallel(cupdlp_float *a_delete, cupdlp_float *a, long long a_len, long long *keep_idx, long long *keep_nnz);
void deleteArray1DElements_by_keep_longlong(cupdlp_float *a_delete, cupdlp_float *a, long long a_len, cupdlp_bool *keep);
void recoverArray1DElements(cupdlp_float *a_recover, cupdlp_float *a_delete, cupdlp_int a_recover_len, cupdlp_int a_delete_len, cupdlp_int *zero_idx, cupdlp_int *zero_idx_len);
void recoverArray1DElements_longlong(cupdlp_float *a_recover, cupdlp_float *a_delete, long long a_recover_len, long long a_delete_len, long long *zero_idx, long long *zero_idx_len);
void recoverArray1DElements_byKeepIdx_longlong_parallel(cupdlp_float *a_recover, cupdlp_float *a_delete, long long a_recover_len, long long *keep_nnz, long long *keep_idx);
void testPtr(cupdlp_float *a);
void construct_and_solve_Multiscale(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_bool *ifChangeIntParam, cupdlp_bool *ifChangeFloatParam, cupdlp_int *intParam, cupdlp_float *floatParam, cupdlp_bool ifSaveSol, cupdlp_int coarse_degree, cupdlp_int coarse_degree_last, cupdlp_float **x_solution, cupdlp_float **y_solution, cupdlp_float *x_solution_last, cupdlp_float *y_solution_last);
void construct_and_solve_Multiscale_longlong(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_bool *ifChangeIntParam, cupdlp_bool *ifChangeFloatParam, cupdlp_int *intParam, cupdlp_float *floatParam, cupdlp_bool ifSaveSol, cupdlp_int coarse_degree, cupdlp_int coarse_degree_last, cupdlp_float **x_solution, cupdlp_float **y_solution, cupdlp_float *x_solution_last, cupdlp_float *y_solution_last);
void directly_construct_and_solve_Multiscale_longlong(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_bool *ifChangeIntParam, cupdlp_bool *ifChangeFloatParam, cupdlp_int *intParam, cupdlp_float *floatParam, cupdlp_bool ifSaveSol, cupdlp_int coarse_degree, cupdlp_int coarse_degree_last, cupdlp_float **x_solution, cupdlp_float **y_solution, cupdlp_float *x_solution_last, cupdlp_float *y_solution_last);
void construct_and_solve_Multiscale_longlong_pdbalance(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_bool *ifChangeIntParam, cupdlp_bool *ifChangeFloatParam, cupdlp_int *intParam, cupdlp_float *floatParam, cupdlp_bool ifSaveSol, cupdlp_int coarse_degree, cupdlp_int coarse_degree_last, cupdlp_float **x_solution, cupdlp_float **y_solution, cupdlp_float *x_solution_last, cupdlp_float *y_solution_last);
// void construct_and_solve_Multiscale_by_keep_longlong(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_bool *ifChangeIntParam, cupdlp_bool *ifChangeFloatParam, cupdlp_int *intParam, cupdlp_float *floatParam, cupdlp_bool ifSaveSol, cupdlp_int coarse_degree, cupdlp_int coarse_degree_last, cupdlp_float **x_solution, cupdlp_float **y_solution, cupdlp_float *x_solution_last, cupdlp_float *y_solution_last);
// void construct_and_solve_Multiscale_withoutRecover(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_bool *ifChangeIntParam, cupdlp_bool *ifChangeFloatParam, cupdlp_int *intParam, cupdlp_float *floatParam, cupdlp_bool ifSaveSol, cupdlp_int coarse_degree, cupdlp_int coarse_degree_last, cupdlp_float **x_solution, cupdlp_float **y_solution, cupdlp_float *x_solution_last, cupdlp_float *y_solution_last);
void construct_and_solve_Multiscale_test(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_int coarse_degree);
int parallelTest(int N, int num_threads);
void computepPrimalFeas(cupdlp_float *x_solution, cupdlp_int resolution, cupdlp_int coarse_degree);
void compute_q_2norm(cupdlp_float *q_2norm, cupdlp_float *a, cupdlp_float *b, cupdlp_int vec_len);
void pdbalance_dualOT_primal_forward(cupdlp_float *x_init_balance, cupdlp_float *x_init, cupdlp_int x_len, cupdlp_float *balance_weight);
void pdbalance_dualOT_primal_backward(cupdlp_float *x_solution, cupdlp_float *x_solution_balance, cupdlp_int x_len, cupdlp_float *balance_weight);
void scale_floatArray1D(cupdlp_float *a_scaled, cupdlp_float *a, cupdlp_int a_len, cupdlp_float balance_weight);
void compute_2norm_floatArray1D(cupdlp_float *norm, cupdlp_float *a, long long a_len);
void checkTwoArray1D_whether_equal(long long *a, long long *b, long long vec_len);
void dualOT_formulateLP_directly(CUPDLPwork *w, cupdlp_float *a, cupdlp_float *b, cupdlp_int resolution, cupdlp_int coarse_degree, long long *keep_idx, long long *keep_nnz, int *keep_a_idx, int *keep_b_idx, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_int **constraint_new_idx);
void generate_keep_a_b_idx_from_keep_idx(int **keep_a_idx, int **keep_b_idx, long long *keep_idx, long long *keep_nnz, int resolution_now);
void generate_constraint_new_idx(int **constraint_new_idx, int nRows);
void printCUPDLPwork(CUPDLPwork *w);
#endif // CUPDLP_CUPDLP_MULTISCALE_H