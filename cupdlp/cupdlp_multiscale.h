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

cupdlp_float *readCSVToFloatArray(const char *filename, cupdlp_int numRows, cupdlp_int numCols);
void normalizeArray1D(cupdlp_float *array, cupdlp_int array_len);
cupdlp_float *coarsingArray1D(cupdlp_float *array, cupdlp_int array_len, const cupdlp_int coarse_degree);

void print_float_array2D(cupdlp_float **array, cupdlp_int numRows, cupdlp_int numCols);
void print_float_array1D(cupdlp_float *array, cupdlp_int num);
void print_int_array1D(cupdlp_int *array, cupdlp_int num);
cupdlp_float **mergeTwoArrays2D(cupdlp_float **a, cupdlp_float **b, cupdlp_int resolution);
cupdlp_float *mergeTwoArrays1D(const cupdlp_float *a, const cupdlp_float *b, int a_len, int b_len);
cupdlp_float *mergeTwoArrays1D_minus(const cupdlp_float *a, const cupdlp_float *b, cupdlp_int a_len, cupdlp_int b_len);
cupdlp_float *normalizedSquaredEuclideanDistance(cupdlp_int m, cupdlp_int n);

cupdlp_int *dualOT_startArray(cupdlp_int m, cupdlp_int n);
cupdlp_int *dualOT_indexArray(cupdlp_int m, cupdlp_int n);
cupdlp_float *dualOT_valueArray(cupdlp_int m, cupdlp_int n);
cupdlp_float *dualOT_rowLower(cupdlp_int m, cupdlp_int n);
// cupdlp_float *dualOT_rowUpper(cupdlp_int m, cupdlp_int n);
cupdlp_float *dualOT_colLower(cupdlp_int m, cupdlp_int n);
cupdlp_float *dualOT_colUpper(cupdlp_int m, cupdlp_int n);
void generate_dualOT_model_from_distribution_and_cost(void *model, const cupdlp_float *a, const cupdlp_float *b, const cupdlp_int a_len, const cupdlp_int b_len, cupdlp_float *c);
void generate_dualOT_model_from_csv(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution);
void generate_coarse_dualOT_model(void *model_coarse, cupdlp_float *a, cupdlp_float *b, cupdlp_int a_len, cupdlp_int b_len, cupdlp_int resolution, const cupdlp_int coarse_degree);
void LP_Solve_Multiscale(w, w_coarse, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx);
CUPDLPwork *createCUPDLPwork(void *model, CUPDLPscaling *scaling, cupdlp_int ifScaling, cupdlp_bool ifPresolve);
#endif // CUPDLP_CUPDLP_MULTISCALE_H