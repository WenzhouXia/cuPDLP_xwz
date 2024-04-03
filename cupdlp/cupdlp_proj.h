//
// Created by chuwen on 23-11-28.
//

#ifndef CUPDLP_CUPDLP_PROJ_H
#define CUPDLP_CUPDLP_PROJ_H

#include "cupdlp_defs.h"
#include "glbopts.h"

void PDHG_Project_Bounds(CUPDLPwork *work, double *r);

void PDHG_Project_Row_Duals(CUPDLPwork *work, double *r);

void PDHG_Restart_Iterate(CUPDLPwork *pdhg);
void PDTEST_Restart_Iterate(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
void PDTEST_Restart_Iterate_Only_Beta(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
void PDTEST_Restart_Iterate_best(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
void PDTEST_Restart_Iterate_best2(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
void PDTEST_Average_Restart_Iterate(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);

void PDHG_Restart_Iterate_GPU(CUPDLPwork *pdhg);
void PDTEST_Average_Restart_Iterate_GPU(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
void PDTEST_Restart_Iterate_GPU(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
void PDTEST_Restart_Iterate_GPU_Only_Beta(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
void PDTEST_Restart_Iterate_GPU_best(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
void PDTEST_Restart_Iterate_GPU_best1(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
void PDTEST_Restart_Iterate_GPU_best2(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
void PDTEST_Restart_Iterate_GPU_best3(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
void PDTEST_Restart_Iterate_GPU_best4(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
void PDTEST_Restart_Iterate_GPU_best5(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
void PDTEST_Restart_Iterate_GPU_best6(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);

void PDTEST_Restart_Iterate_GPU_Only_WeightUpdate(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);
void PDTEST_Restart_Iterate_GPU_Only_VariableUpdate(CUPDLPwork *pdhg, cupdlp_int *nIter_restart);

#endif // CUPDLP_CUPDLP_PROJ_H
