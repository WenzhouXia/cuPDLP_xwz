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
void normalizedSquaredEuclideanDistance(cupdlp_float *c, cupdlp_int m, cupdlp_int n);
void coarsing_normalizedSquaredEuclideanDistance(cupdlp_float *c_coarse, cupdlp_int m, cupdlp_int n, const cupdlp_int coarse_degree);

cupdlp_int *dualOT_startArray(cupdlp_int m, cupdlp_int n);
cupdlp_int *dualOT_indexArray(cupdlp_int m, cupdlp_int n);
cupdlp_float *dualOT_valueArray(cupdlp_int m, cupdlp_int n);
cupdlp_float *dualOT_rowLower(cupdlp_int m, cupdlp_int n);
cupdlp_float *dualOT_rowLower_delete(cupdlp_int len_after_delete);
// cupdlp_float *dualOT_rowUpper(cupdlp_int m, cupdlp_int n);
cupdlp_float *dualOT_colLower(cupdlp_int m, cupdlp_int n);
cupdlp_float *dualOT_colUpper(cupdlp_int m, cupdlp_int n);
void generate_dualOT_model_from_distribution_and_cost(void *model, const cupdlp_float *a, const cupdlp_float *b, const cupdlp_int a_len, const cupdlp_int b_len, cupdlp_float *c);
void generate_dualOT_model_byMatrix_from_distribution_and_cost(void *model, const cupdlp_float *a, const cupdlp_float *b, const cupdlp_int a_len, const cupdlp_int b_len, cupdlp_float *c, cupdlp_int resolution);
void generate_dualOT_model_delete_byMatrix_from_distribution_and_cost(void *model, const cupdlp_float *a, const cupdlp_float *b, const cupdlp_int a_len, const cupdlp_int b_len, cupdlp_float *c_delete, int resolution, int *zero_idx, int *zero_idx_len);
void generate_dualOT_model_from_csv(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution);
void generate_coarse_dualOT_model_from_csv(void *model_coarse, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, const cupdlp_int coarse_degree);
void generate_dualOT_model_delete_from_csv(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_int *zero_idx, cupdlp_int *zero_idx_len, cupdlp_int len_after_delete);
void generate_coarse_dualOT_model_delete_from_csv(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_int *zero_idx, cupdlp_int *zero_idx_len, cupdlp_int len_after_delete, cupdlp_int coarse_degree);
void generate_coarse_dualOT_model_delete_from_csv_longlong(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, long long *zero_idx, long long *zero_idx_len, long long len_after_delete, cupdlp_int coarse_degree);
// void LP_Solve_Multiscale(w, w_coarse, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx);
void createCUPDLPwork(CUPDLPwork *w, void *model, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_int **constraint_new_idx, int *nCols_origin_ptr, int *nRows_ptr);
void createCUPDLPwork_clear(CUPDLPwork *w, void *model, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_int **constraint_new_idx);
void fine_dualOT_primal(cupdlp_float *x_init, cupdlp_float *x_coarse_solution, cupdlp_int x_len, cupdlp_int x_coarse_len, cupdlp_int resolution, cupdlp_int coarse_degree);
void fine_dualOT_dual(cupdlp_float *y_init, cupdlp_float *y_coarse_solution, long long y_len, long long y_coarse_len, cupdlp_int resolution, cupdlp_int coarse_degree);
int countArray1D_Smaller_than_threshold(cupdlp_float *a, int a_len, cupdlp_float thr);
long long countArray1D_Smaller_than_threshold_longlong(cupdlp_float *a, long long a_len, cupdlp_float thr);
int *countArray1D_Smaller_than_threshold_with_Record(cupdlp_float *a, int a_len, int *a_record_len, cupdlp_float thr);
long long *countArray1D_Smaller_than_threshold_with_Record_longlong(cupdlp_float *a, long long a_len, long long *a_record_len, cupdlp_float thr);
void saveArray1D_to_csv(cupdlp_float *a, int a_len, const char *filename);
void analyseArray1D(cupdlp_float *a, long long a_len, cupdlp_float thr, const char *filename);
int countArray1D_same_elements(int *a, int a_len, int *b, int b_len);
void compareTwoArray1D(cupdlp_float *a, cupdlp_int a_len, cupdlp_float *b, cupdlp_int b_len, cupdlp_float thr);
void constructCoinPackedMatrix(void *mat, cupdlp_int resolution);
void deleteArray1DElements(cupdlp_float *a_dele, cupdlp_float *a, cupdlp_int a_len, cupdlp_int *zero_idx, cupdlp_int *zero_idx_len);
void deleteArray1DElements_longlong(cupdlp_float *a_dele, cupdlp_float *a, long long a_len, long long *zero_idx, long long *zero_idx_len);
void recoverArray1DElements(cupdlp_float *a_recover, cupdlp_float *a_delete, cupdlp_int a_recover_len, cupdlp_int a_delete_len, cupdlp_int *zero_idx, cupdlp_int *zero_idx_len);
void recoverArray1DElements_longlong(cupdlp_float *a_recover, cupdlp_float *a_delete, long long a_recover_len, long long a_delete_len, long long *zero_idx, long long *zero_idx_len);
void testPtr(cupdlp_float *a);
void construct_and_solve_Multiscale(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_bool *ifChangeIntParam, cupdlp_bool *ifChangeFloatParam, cupdlp_int *intParam, cupdlp_float *floatParam, cupdlp_bool ifSaveSol, cupdlp_int coarse_degree, cupdlp_int coarse_degree_last, cupdlp_float **x_solution, cupdlp_float **y_solution, cupdlp_float *x_solution_last, cupdlp_float *y_solution_last);
void construct_and_solve_Multiscale_longlong(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_bool *ifChangeIntParam, cupdlp_bool *ifChangeFloatParam, cupdlp_int *intParam, cupdlp_float *floatParam, cupdlp_bool ifSaveSol, cupdlp_int coarse_degree, cupdlp_int coarse_degree_last, cupdlp_float **x_solution, cupdlp_float **y_solution, cupdlp_float *x_solution_last, cupdlp_float *y_solution_last);
void construct_and_solve_Multiscale_withoutRecover(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_bool *ifChangeIntParam, cupdlp_bool *ifChangeFloatParam, cupdlp_int *intParam, cupdlp_float *floatParam, cupdlp_bool ifSaveSol, cupdlp_int coarse_degree, cupdlp_int coarse_degree_last, cupdlp_float **x_solution, cupdlp_float **y_solution, cupdlp_float *x_solution_last, cupdlp_float *y_solution_last);
void construct_and_solve_Multiscale_test(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_int coarse_degree);
#endif // CUPDLP_CUPDLP_MULTISCALE_H