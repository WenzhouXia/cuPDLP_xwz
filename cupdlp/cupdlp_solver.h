//
// Created by chuwen on 23-11-27.
//

#ifndef CUPDLP_CUPDLP_SOLVER_H
#define CUPDLP_CUPDLP_SOLVER_H

#include "cupdlp_defs.h"
#include "glbopts.h"

#ifdef __cplusplus
extern "C"
{
#endif
#define CUPDLP_CHECK_TIMEOUT(pdhg)                             \
  {                                                            \
    PDHG_Compute_SolvingTime(pdhg);                            \
    if (pdhg->timers->dSolvingTime > pdhg->settings->dTimeLim) \
    {                                                          \
      retcode = RETCODE_FAILED;                                \
      goto exit_cleanup;                                       \
    }                                                          \
  }

  void PDHG_Compute_Primal_Feasibility(CUPDLPwork *work, double *primalResidual,
                                       const double *ax, const double *x,
                                       double *dPrimalFeasibility,
                                       double *dPrimalObj);

  void PDHG_Compute_Dual_Feasibility(CUPDLPwork *work, double *dualResidual,
                                     const double *aty, const double *x,
                                     const double *y, double *dDualFeasibility,
                                     double *dDualObj, double *dComplementarity);
  void PDTEST_printCudaDenseVecGPU(const CUPDLPvec *vec);
  void PDTEST_printCudaMatGPU(CUPDLPwork *work);
  void PDTEST_printCudafloatGPU(cupdlp_float *data, int size);

  void PDTEST_Compute_dDualObj(CUPDLPwork *work);

  void PDHG_Compute_Residuals(CUPDLPwork *work);
  void PDTEST_Average_Compute_Residuals(CUPDLPwork *work);
  void PDTEST_Compute_Residuals_best(CUPDLPwork *work);
  void PDTEST_Compute_Residuals(CUPDLPwork *work);
  void PDTEST_Compute_Residuals_CurrentandAverage(CUPDLPwork *work);
  void PDHG_Init_Variables_Multiscale(CUPDLPwork *work, cupdlp_float *x_init, cupdlp_float *y_init);
  void PDHG_Init_Variables(CUPDLPwork *work);
  void PDTEST_Average_Init_Variables(CUPDLPwork *work);
  void PDTEST_Init_Variables(CUPDLPwork *work);

  void PDHG_Check_Data(CUPDLPwork *work);

  cupdlp_bool PDTEST_Check_Termination_best(CUPDLPwork *pdhg, int bool_print);
  cupdlp_bool PDHG_Check_Termination(CUPDLPwork *pdhg, int bool_print);

  cupdlp_bool PDHG_Check_Termination_Average(CUPDLPwork *pdhg, int bool_print);

  void PDHG_Print_Header(CUPDLPwork *pdhg);

  void PDHG_Print_Iter(CUPDLPwork *pdhg);

  void PDHG_Print_Iter_Average(CUPDLPwork *pdhg);

  void PDHG_Compute_SolvingTime(CUPDLPwork *pdhg);
  cupdlp_retcode PDHG_Solve_Multiscale_withStepsize(CUPDLPwork *pdhg, cupdlp_float *x_init, cupdlp_float *y_init, cupdlp_float *stepsize_last, cupdlp_float *weight_last, cupdlp_float *stepsize_init, cupdlp_float *weight_init, cupdlp_bool whether_first);
  cupdlp_retcode PDHG_Solve(CUPDLPwork *pdhg);
  cupdlp_retcode PDHG_Solve_AdapTheta(CUPDLPwork *pdhg);
  cupdlp_retcode PDTEST_Average_Solve(CUPDLPwork *pdhg);
  cupdlp_retcode PDTEST_Solve(CUPDLPwork *pdhg);
  cupdlp_retcode PDTEST_Solve_best(CUPDLPwork *pdhg);
  cupdlp_retcode PDTEST_Solve_best2(CUPDLPwork *pdhg);

  void PDHG_PostSolve(CUPDLPwork *pdhg, cupdlp_int nCols_origin,
                      cupdlp_int *constraint_new_idx, cupdlp_float *x_origin,
                      cupdlp_float *y_origin);

  void PDTEST_PostSolve(CUPDLPwork *pdhg, cupdlp_int nCols_origin,
                        cupdlp_int *constraint_new_idx, cupdlp_float *x_origin,
                        cupdlp_float *y_origin);
  cupdlp_retcode LP_SolvePDHG_Multiscale_withStepsize(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_bool *ifChangeFloatParam, cupdlp_float *floatParam, char *fp, cupdlp_bool ifSaveSol, cupdlp_int *constraint_new_idx, cupdlp_float **x_solution, cupdlp_float **y_solution, cupdlp_float *x_init, cupdlp_float *y_init, cupdlp_float *stepsize_last, cupdlp_float *weight_last, cupdlp_float *stepsize_init, cupdlp_float *weight_init, cupdlp_bool whether_first);
  cupdlp_retcode LP_SolvePDHG_Multiscale(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_bool *ifChangeFloatParam, cupdlp_float *floatParam, char *fp, cupdlp_bool ifSaveSol, cupdlp_int *constraint_new_idx, cupdlp_float **x_solution, cupdlp_float **y_solution, cupdlp_float *x_init, cupdlp_float *y_init);
  cupdlp_retcode LP_SolvePDHG(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_bool *ifChangeFloatParam, cupdlp_float *floatParam, char *fp, cupdlp_bool ifSaveSol, cupdlp_int *constraint_new_idx, cupdlp_int resolution, cupdlp_float **x_solution, cupdlp_float **y_solution);

  cupdlp_retcode LP_SolvePDHG_AdapTheta(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam,
                                        cupdlp_int *intParam,
                                        cupdlp_bool *ifChangeFloatParam,
                                        cupdlp_float *floatParam, char *fp,
                                        cupdlp_float *x_origin, cupdlp_int nCols_origin,
                                        cupdlp_float *y_origin, cupdlp_bool ifSaveSol,
                                        cupdlp_int *constraint_new_idx);

  cupdlp_retcode LP_SolvePDTEST(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam,
                                cupdlp_int *intParam,
                                cupdlp_bool *ifChangeFloatParam,
                                cupdlp_float *floatParam, char *fp,
                                cupdlp_float *x_origin, cupdlp_int nCols_origin,
                                cupdlp_float *y_origin, cupdlp_bool ifSaveSol,
                                cupdlp_int *constraint_new_idx);
  cupdlp_retcode LP_SolvePDTEST_best(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam,
                                     cupdlp_int *intParam,
                                     cupdlp_bool *ifChangeFloatParam,
                                     cupdlp_float *floatParam, char *fp,
                                     cupdlp_float *x_origin, cupdlp_int nCols_origin,
                                     cupdlp_float *y_origin, cupdlp_bool ifSaveSol,
                                     cupdlp_int *constraint_new_idx);
  cupdlp_retcode LP_SolvePDTEST_best2(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam,
                                      cupdlp_int *intParam,
                                      cupdlp_bool *ifChangeFloatParam,
                                      cupdlp_float *floatParam, char *fp,
                                      cupdlp_float *x_origin, cupdlp_int nCols_origin,
                                      cupdlp_float *y_origin, cupdlp_bool ifSaveSol,
                                      cupdlp_int *constraint_new_idx);
  cupdlp_retcode LP_SolvePDTEST_Average(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam,
                                        cupdlp_int *intParam,
                                        cupdlp_bool *ifChangeFloatParam,
                                        cupdlp_float *floatParam, char *fp,
                                        cupdlp_float *x_origin, cupdlp_int nCols_origin,
                                        cupdlp_float *y_origin, cupdlp_bool ifSaveSol,
                                        cupdlp_int *constraint_new_idx);

  cupdlp_retcode LP_SolvePDTEST_min(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam,
                                    cupdlp_int *intParam,
                                    cupdlp_bool *ifChangeFloatParam,
                                    cupdlp_float *floatParam, char *fp,
                                    cupdlp_float *x_origin, cupdlp_int nCols_origin, cupdlp_float *y_origin, cupdlp_bool ifSaveSol, cupdlp_int *constraint_new_idx);

#ifdef __cplusplus
}
#endif
#endif // CUPDLP_CUPDLP_SOLVER_H
