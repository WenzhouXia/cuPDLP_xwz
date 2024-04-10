//
// Created by chuwen on 23-11-28.
//

#ifndef CUPDLP_CUPDLP_STEP_H
#define CUPDLP_CUPDLP_STEP_H

#include "cupdlp_defs.h"
// #include "cupdlp_scaling.h"
#include "glbopts.h"

void PDTEST_x_md_step(CUPDLPwork *work, cupdlp_float beta);
void PDTEST_x_ag_step(CUPDLPwork *work, cupdlp_float beta);
void PDTEST_y_ag_step(CUPDLPwork *work, cupdlp_float beta);
void PDTEST_x_bar_step(CUPDLPwork *work, cupdlp_float theta);
void PDTEST_primalGradientStep(CUPDLPwork *work, cupdlp_float dPrimalStepSize);
void PDTEST_dualGradientStep(CUPDLPwork *work, cupdlp_float dDualStepSize);

void PDTEST_x_ag_step_ReverseOrder(CUPDLPwork *work, cupdlp_float beta);
void PDTEST_y_ag_step_ReverseOrder(CUPDLPwork *work, cupdlp_float beta);
void PDTEST_x_bar_step_ReverseOrder(CUPDLPwork *work, cupdlp_float theta);
void PDTEST_primalGradientStep_ReverseOrder(CUPDLPwork *work, cupdlp_float dPrimalStepSize);
void PDTEST_dualGradientStep_ReverseOrder(CUPDLPwork *work, cupdlp_float dDualStepSize);

cupdlp_retcode PDHG_Power_Method(CUPDLPwork *work, double *lambda);
cupdlp_retcode PDTEST_Power_Method(CUPDLPwork *work, double *lambda);

void PDHG_Compute_Step_Size_Ratio(CUPDLPwork *pdhg);
void PDTEST_Compute_Step_Size_Ratio(CUPDLPwork *pdhg);

void PDHG_Update_Iterate_Constant_Step_Size(CUPDLPwork *pdhg);
void PDTEST_Update_Iterate_Constant_Step_Size(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
void PDTEST_Update_Iterate_Constant_Step_Size_ReverseOrder(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);

void PDHG_Update_Iterate_Malitsky_Pock(CUPDLPwork *pdhg);

cupdlp_retcode PDHG_Update_Iterate_Adaptive_Step_Size(CUPDLPwork *pdhg);
cupdlp_retcode PDHG_Update_Iterate_Adaptive_Step_Size_AdapTheta(CUPDLPwork *pdhg);
cupdlp_retcode PDTEST_Update_Iterate_Adaptive_Step_Size_ag(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
cupdlp_retcode PDTEST_Update_Iterate_Adaptive_Step_Size_ag_best(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
cupdlp_retcode PDTEST_Update_Iterate_Adaptive_Step_Size_ag_best2(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
cupdlp_retcode PDTEST_Update_Iterate_Adaptive_Step_Size_ag_best3(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
cupdlp_retcode PDTEST_Update_Iterate_Adaptive_Step_Size_ag_best4(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
cupdlp_retcode PDTEST_Update_Iterate_Adaptive_Step_Size_ag_best5(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
cupdlp_retcode PDTEST_Update_Iterate_Adaptive_Step_Size_ag_best6(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
cupdlp_retcode PDTEST_Update_Iterate_Adaptive_Step_Size(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);

cupdlp_retcode PDHG_Init_Step_Sizes(CUPDLPwork *pdhg);
cupdlp_retcode PDHG_Init_Step_Sizes_Multiscale(CUPDLPwork *pdhg, cupdlp_float *stepsize_init, cupdlp_float *weight_init, cupdlp_bool whether_first);
cupdlp_retcode PDHG_Init_Step_Sizes_AdapTheta(CUPDLPwork *pdhg);
cupdlp_retcode PDTEST_Init_Step_Sizes(CUPDLPwork *pdhg);
cupdlp_retcode PDTEST_Init_Step_Sizes_best(CUPDLPwork *pdhg);
cupdlp_retcode PDTEST_Init_Step_Sizes_best1(CUPDLPwork *pdhg);
cupdlp_retcode PDTEST_Init_Step_Sizes_best2(CUPDLPwork *pdhg);
cupdlp_retcode PDTEST_Init_Step_Sizes_best3(CUPDLPwork *pdhg);
cupdlp_retcode PDTEST_Init_Step_Sizes_best4(CUPDLPwork *pdhg);
cupdlp_retcode PDTEST_Init_Step_Sizes_best5(CUPDLPwork *pdhg);
void PDHG_Compute_Average_Iterate(CUPDLPwork *work);
void PDTEST_Compute_Average_Iterate(CUPDLPwork *work);

void PDHG_Update_Average(CUPDLPwork *work);
void PDTEST_Update_Average(CUPDLPwork *work);

cupdlp_retcode PDHG_Update_Iterate(CUPDLPwork *pdhg);
cupdlp_retcode PDHG_Update_Iterate_AdapTheta(CUPDLPwork *pdhg);
cupdlp_retcode PDTEST_Average_Update_Iterate(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
cupdlp_retcode PDTEST_Update_Iterate(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
cupdlp_retcode PDTEST_Update_Iterate_best(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
cupdlp_retcode PDTEST_Update_Iterate_best2(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);

void PDHG_primalGradientStep(CUPDLPwork *work, cupdlp_float dPrimalStepSize);
void PDHG_dualGradientStep(CUPDLPwork *work, cupdlp_float dDualStepSize);

#endif // CUPDLP_CUPDLP_STEP_H
