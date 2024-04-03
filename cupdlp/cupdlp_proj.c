//
// Created by chuwen on 23-11-28.
//

#include "cupdlp_proj.h"

#include "cupdlp_defs.h"
#include "cupdlp_linalg.h"
#include "cupdlp_restart.h"
// #include "cupdlp_scaling.h"
#include "cupdlp_solver.h"
#include "cupdlp_step.h"
#include "cupdlp_utils.h"
#include "glbopts.h"

// primal projection: project x to [lower, upper]
void PDHG_Project_Bounds(CUPDLPwork *work, cupdlp_float *r)
{
  CUPDLPproblem *problem = work->problem;

  // cupdlp_projUpperBound(r, r, problem->upper, problem->nCols);
  // cupdlp_projLowerBound(r, r, problem->lower, problem->nCols);

  cupdlp_projub(r, problem->upper, problem->nCols);
  cupdlp_projlb(r, problem->lower, problem->nCols);
}

void PDHG_Project_Row_Duals(CUPDLPwork *work, cupdlp_float *r)
{
  CUPDLPproblem *problem = work->problem;

  // cupdlp_projPositive(r + problem->nEqs, r + problem->nEqs, problem->nRows -
  // problem->nEqs);
  cupdlp_projPos(r + problem->nEqs, problem->nRows - problem->nEqs);
}

// void PDHG_Restart_Iterate(CUPDLPwork *pdhg)
// {
//     CUPDLPproblem *problem = pdhg->problem;
//     CUPDLPiterates *iterates = pdhg->iterates;
//     CUPDLPstepsize *stepsize = pdhg->stepsize;
//     CUPDLPtimers *timers = pdhg->timers;

//     // PDHG_Compute_Average_Iterate(pdhg);
//     PDHG_restart_choice restart_choice = PDHG_Check_Restart(pdhg);

//     if (restart_choice == PDHG_NO_RESTART)
//         return;

//     PDHG_Compute_Step_Size_Ratio(pdhg);

//     stepsize->dSumPrimalStep = 0.0;
//     stepsize->dSumDualStep = 0.0;
//     cupdlp_zero(iterates->xSum, cupdlp_float, problem->nCols);
//     cupdlp_zero(iterates->ySum, cupdlp_float, problem->nRows);

//     if (restart_choice == PDHG_RESTART_TO_AVERAGE)
//     {
//         cupdlp_copy(iterates->x, iterates->xAverage, cupdlp_float,
//         problem->nCols); cupdlp_copy(iterates->y, iterates->yAverage,
//         cupdlp_float, problem->nRows); cupdlp_copy(iterates->ax,
//         iterates->axAverage, cupdlp_float, problem->nRows);
//         cupdlp_copy(iterates->aty, iterates->atyAverage, cupdlp_float,
//         problem->nCols);
//     }
//     cupdlp_copy(iterates->xLastRestart, iterates->x, cupdlp_float,
//     problem->nCols); cupdlp_copy(iterates->yLastRestart, iterates->y,
//     cupdlp_float, problem->nRows);

//     iterates->iLastRestartIter = timers->nIter;

//     PDHG_Compute_Residuals(pdhg);
//     // cupdlp_printf("Recomputed stepsize ratio: %e,  sqrt(ratio)=%e",
//     stepsize->dBeta, sqrt(stepsize->dBeta));
// }

void PDHG_Restart_Iterate(CUPDLPwork *pdhg)
{
  switch (pdhg->settings->eRestartMethod)
  {
  case PDHG_WITHOUT_RESTART:
    break;
  case PDHG_GPU_RESTART:
    PDHG_Restart_Iterate_GPU(pdhg);
    break;
  case PDHG_CPU_RESTART:
    // TODO: implement PDHG_Restart_Iterate_CPU(pdhg);
    break;
  }
}

void PDTEST_Restart_Iterate(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  switch (pdhg->settings->eRestartMethod)
  {
  case PDHG_WITHOUT_RESTART:
    // cupdlp_printf("PDHG_WITHOUT_RESTART\n");
    break;
  case PDHG_GPU_RESTART:
    // cupdlp_printf("PDHG_GPU_RESTART\n");
    PDTEST_Restart_Iterate_GPU(pdhg, nIter_restart);
    break;
    // case PDHG_CPU_RESTART:
    //   // TODO: implement PDHG_Restart_Iterate_CPU(pdhg);
    //   break;
  }
}

void PDTEST_Restart_Iterate_Only_Beta(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  switch (pdhg->settings->eRestartMethod)
  {
  case PDHG_WITHOUT_RESTART:
    // cupdlp_printf("PDHG_WITHOUT_RESTART\n");
    break;
  case PDHG_GPU_RESTART:
    // cupdlp_printf("PDHG_GPU_RESTART\n");
    PDTEST_Restart_Iterate_GPU_Only_Beta(pdhg, nIter_restart);
    break;
    // case PDHG_CPU_RESTART:
    //   // TODO: implement PDHG_Restart_Iterate_CPU(pdhg);
    //   break;
  }
}

void PDTEST_Restart_Iterate_best(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  switch (pdhg->settings->eRestartMethod)
  {
  case PDHG_WITHOUT_RESTART:
    // cupdlp_printf("PDHG_WITHOUT_RESTART\n");
    break;
  case PDHG_GPU_RESTART:
    // cupdlp_printf("PDHG_GPU_RESTART\n");
    PDTEST_Restart_Iterate_GPU_best(pdhg, nIter_restart);
    break;
    // case PDHG_CPU_RESTART:
    //   // TODO: implement PDHG_Restart_Iterate_CPU(pdhg);
    //   break;
  }
}
void PDTEST_Restart_Iterate_best2(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  switch (pdhg->settings->eRestartMethod)
  {
  case PDHG_WITHOUT_RESTART:
    // cupdlp_printf("PDHG_WITHOUT_RESTART\n");
    break;
  case PDHG_GPU_RESTART:
    // cupdlp_printf("PDHG_GPU_RESTART\n");
    PDTEST_Restart_Iterate_GPU_best2(pdhg, nIter_restart);
    break;
    // case PDHG_CPU_RESTART:
    //   // TODO: implement PDHG_Restart_Iterate_CPU(pdhg);
    //   break;
  }
}
void PDTEST_Average_Restart_Iterate(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  switch (pdhg->settings->eRestartMethod)
  {
  case PDHG_WITHOUT_RESTART:
    // cupdlp_printf("PDHG_WITHOUT_RESTART\n");
    break;
  case PDHG_GPU_RESTART:
    PDTEST_Average_Restart_Iterate_GPU(pdhg, nIter_restart);
    break;
    // case PDHG_CPU_RESTART:
    //   // TODO: implement PDHG_Restart_Iterate_CPU(pdhg);
    //   break;
  }
}

void PDHG_Restart_Iterate_GPU(CUPDLPwork *pdhg)
{
  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPtimers *timers = pdhg->timers;

  // PDHG_Compute_Average_Iterate(pdhg);
  PDHG_restart_choice restart_choice = PDHG_Check_Restart_GPU(pdhg);

  if (restart_choice == PDHG_NO_RESTART)
    return;

  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;
  CUPDLP_ZERO_VEC(iterates->xSum, cupdlp_float, problem->nCols);
  CUPDLP_ZERO_VEC(iterates->ySum, cupdlp_float, problem->nRows);

  if (restart_choice == PDHG_RESTART_TO_AVERAGE)
  {
    cupdlp_printf("Restart to average\n");
    resobj->dPrimalFeasLastRestart = resobj->dPrimalFeasAverage;
    resobj->dDualFeasLastRestart = resobj->dDualFeasAverage;
    resobj->dDualityGapLastRestart = resobj->dDualityGapAverage;

    // cupdlp_copy(iterates->x, iterates->xAverage, cupdlp_float,
    // problem->nCols); cupdlp_copy(iterates->y, iterates->yAverage,
    // cupdlp_float, problem->nRows); cupdlp_copy(iterates->ax,
    // iterates->axAverage, cupdlp_float, problem->nRows);
    // cupdlp_copy(iterates->aty, iterates->atyAverage, cupdlp_float,
    // problem->nCols);

    CUPDLP_COPY_VEC(iterates->x->data, iterates->xAverage->data, cupdlp_float,
                    problem->nCols);
    CUPDLP_COPY_VEC(iterates->y->data, iterates->yAverage->data, cupdlp_float,
                    problem->nRows);
    CUPDLP_COPY_VEC(iterates->ax->data, iterates->axAverage->data, cupdlp_float,
                    problem->nRows);
    CUPDLP_COPY_VEC(iterates->aty->data, iterates->atyAverage->data,
                    cupdlp_float, problem->nCols);
  }
  else
  {
    cupdlp_printf("Restart to current\n");
    resobj->dPrimalFeasLastRestart = resobj->dPrimalFeas;
    resobj->dDualFeasLastRestart = resobj->dDualFeas;
    resobj->dDualityGapLastRestart = resobj->dDualityGap;
  }

  PDHG_Compute_Step_Size_Ratio(pdhg);

  // cupdlp_copy(iterates->xLastRestart, iterates->x, cupdlp_float,
  // problem->nCols); cupdlp_copy(iterates->yLastRestart, iterates->y,
  // cupdlp_float, problem->nRows);
  CUPDLP_COPY_VEC(iterates->xLastRestart, iterates->x->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->yLastRestart, iterates->y->data, cupdlp_float,
                  problem->nRows);

  iterates->iLastRestartIter = timers->nIter;

  PDHG_Compute_Residuals(pdhg);
  // cupdlp_printf("Recomputed stepsize ratio: %e,  sqrt(ratio)=%e",
  // stepsize->dBeta, sqrt(stepsize->dBeta));
}
void PDTEST_Average_Restart_Iterate_GPU(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  CUPDLPproblem *problem = pdhg->problem;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPtimers *timers = pdhg->timers;

  // PDHG_Compute_Average_Iterate(pdhg);
  PDHG_restart_choice restart_choice = PDTEST_Check_Restart_GPU(pdhg);

  if (restart_choice == PDHG_NO_RESTART)
    return;

  // 如果restart了，就把nIter_restart置为0
  cupdlp_printf("Restart，重置nIter_restart为0\n");
  *nIter_restart = 0;

  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;
  CUPDLP_ZERO_VEC(iterates->x_agSum, cupdlp_float, problem->nCols);
  CUPDLP_ZERO_VEC(iterates->y_agSum, cupdlp_float, problem->nRows);

  if (restart_choice == PDHG_RESTART_TO_AVERAGE)
  {
    cupdlp_printf("Restart to average\n");
    resobj->dPrimalFeasLastRestart = resobj->dPrimalFeasAverage;
    resobj->dDualFeasLastRestart = resobj->dDualFeasAverage;
    resobj->dDualityGapLastRestart = resobj->dDualityGapAverage;

    // 需要思考：拷贝给x_ag还是x？
    CUPDLP_COPY_VEC(iterates->x->data, iterates->x_agAverage->data, cupdlp_float,
                    problem->nCols);
    CUPDLP_COPY_VEC(iterates->y->data, iterates->y_agAverage->data, cupdlp_float,
                    problem->nRows);
    CUPDLP_COPY_VEC(iterates->ax->data, iterates->ax_agAverage->data, cupdlp_float,
                    problem->nRows);
    CUPDLP_COPY_VEC(iterates->aty->data, iterates->aty_agAverage->data,
                    cupdlp_float, problem->nCols);

    CUPDLP_COPY_VEC(iterates->x_ag->data, iterates->x_agAverage->data, cupdlp_float,
                    problem->nCols);
    CUPDLP_COPY_VEC(iterates->y_ag->data, iterates->y_agAverage->data, cupdlp_float,
                    problem->nRows);
    CUPDLP_COPY_VEC(iterates->ax_ag->data, iterates->ax_agAverage->data, cupdlp_float,
                    problem->nRows);
    CUPDLP_COPY_VEC(iterates->aty_ag->data, iterates->aty_agAverage->data,
                    cupdlp_float, problem->nCols);

    CUPDLP_COPY_VEC(iterates->x_bar->data, iterates->x_agAverage->data, cupdlp_float,
                    problem->nCols);
    CUPDLP_COPY_VEC(iterates->ax_bar->data, iterates->ax_agAverage->data, cupdlp_float,
                    problem->nRows);
  }
  else
  {
    cupdlp_printf("Restart to current\n");
    resobj->dPrimalFeasLastRestart = resobj->dPrimalFeas;
    resobj->dDualFeasLastRestart = resobj->dDualFeas;
    resobj->dDualityGapLastRestart = resobj->dDualityGap;
    CUPDLP_COPY_VEC(iterates->x->data, iterates->x_ag->data, cupdlp_float,
                    problem->nCols);
    CUPDLP_COPY_VEC(iterates->y->data, iterates->y_ag->data, cupdlp_float,
                    problem->nRows);
    CUPDLP_COPY_VEC(iterates->ax->data, iterates->ax_ag->data, cupdlp_float,
                    problem->nRows);
    CUPDLP_COPY_VEC(iterates->aty->data, iterates->aty_ag->data,
                    cupdlp_float, problem->nCols);

    CUPDLP_COPY_VEC(iterates->x_bar->data, iterates->x_ag->data, cupdlp_float,
                    problem->nCols);
    CUPDLP_COPY_VEC(iterates->ax_bar->data, iterates->ax_ag->data, cupdlp_float,
                    problem->nRows);
  }
  //////////////////////////////////////
  PDTEST_Compute_Step_Size_Ratio(pdhg);
  //////////////////////////////////////

  CUPDLP_COPY_VEC(iterates->xLastRestart, iterates->x_ag->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->yLastRestart, iterates->y_ag->data, cupdlp_float,
                  problem->nRows);

  iterates->iLastRestartIter = timers->nIter;

  PDTEST_Average_Compute_Residuals(pdhg);
  // cupdlp_printf("Recomputed stepsize ratio: %e,  sqrt(ratio)=%e",
  // stepsize->dBeta, sqrt(stepsize->dBeta));
}
void PDTEST_Restart_Iterate_GPU(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  CUPDLPproblem *problem = pdhg->problem;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPtimers *timers = pdhg->timers;

  // PDHG_Compute_Average_Iterate(pdhg);
  PDHG_restart_choice restart_choice = PDTEST_Check_Restart_GPU(pdhg);

  if (restart_choice == PDHG_NO_RESTART)
    return;
  // 如果restart了，就把nIter_restart置为0
  *nIter_restart = 0;

  resobj->dPrimalFeasLastRestart = resobj->dPrimalFeas;
  resobj->dDualFeasLastRestart = resobj->dDualFeas;
  resobj->dDualityGapLastRestart = resobj->dDualityGap;

  // PDTEST_Compute_Step_Size_Ratio(pdhg);

  CUPDLP_COPY_VEC(iterates->xLastRestart, iterates->x_ag->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->yLastRestart, iterates->y_ag->data, cupdlp_float,
                  problem->nRows);

  iterates->iLastRestartIter = timers->nIter;

  PDTEST_Compute_Residuals(pdhg);
  // cupdlp_printf("Recomputed stepsize ratio: %e,  sqrt(ratio)=%e",
  // stepsize->dBeta, sqrt(stepsize->dBeta));
}
void PDTEST_Restart_Iterate_GPU_Only_Beta(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  CUPDLPproblem *problem = pdhg->problem;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPtimers *timers = pdhg->timers;

  PDHG_restart_choice restart_choice = PDTEST_Check_Restart_GPU_Only_Current(pdhg);

  if (restart_choice == PDHG_NO_RESTART)
    return;
  ///////////////////////////////////////////////////
  // 如果restart了，就把nIter_restart置为0
  // *nIter_restart = 0;
  // stepsize->dBeta_ag = 0.0;
  stepsize->dBeta_ag = 1.0;
  // CUPDLP_COPY_VEC(iterates->x->data, iterates->x_ag->data, cupdlp_float,
  //                 problem->nCols);
  // CUPDLP_COPY_VEC(iterates->y->data, iterates->y_ag->data, cupdlp_float,
  //                 problem->nRows);
  // CUPDLP_COPY_VEC(iterates->ax->data, iterates->ax_ag->data, cupdlp_float,
  //                 problem->nRows);
  // CUPDLP_COPY_VEC(iterates->aty->data, iterates->aty_ag->data,
  //                 cupdlp_float, problem->nCols);

  // CUPDLP_COPY_VEC(iterates->x_bar->data, iterates->x_ag->data, cupdlp_float,
  //                 problem->nCols);
  // CUPDLP_COPY_VEC(iterates->ax_bar->data, iterates->ax_ag->data, cupdlp_float,
  //                 problem->nRows);
  ///////////////////////////////////////////////////
  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;

  resobj->dPrimalFeasLastRestart = resobj->dPrimalFeas;
  resobj->dDualFeasLastRestart = resobj->dDualFeas;
  resobj->dDualityGapLastRestart = resobj->dDualityGap;

  ////////////////////////////////////////////////////////
  PDTEST_Compute_Step_Size_Ratio(pdhg);
  ////////////////////////////////////////////////////////

  CUPDLP_COPY_VEC(iterates->xLastRestart, iterates->x_ag->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->yLastRestart, iterates->y_ag->data, cupdlp_float,
                  problem->nRows);

  iterates->iLastRestartIter = timers->nIter;

  PDTEST_Compute_Residuals(pdhg);
  // cupdlp_printf("Recomputed stepsize ratio: %e,  sqrt(ratio)=%e",
  // stepsize->dBeta, sqrt(stepsize->dBeta));
}

void PDTEST_Restart_Iterate_GPU_best(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  CUPDLPsettings *settings = pdhg->settings;
  cupdlp_int bestID = settings->bestID;
  cupdlp_printf("bestID = %d\n", bestID);
  switch (bestID)
  {
  case 1:
    cupdlp_printf("case 1");
    PDTEST_Restart_Iterate_GPU_best1(pdhg, nIter_restart);
    break;
  case 2:
    cupdlp_printf("case 2");
    PDTEST_Restart_Iterate_GPU_best2(pdhg, nIter_restart);
    break;
  case 3:
    cupdlp_printf("case 3");
    PDTEST_Restart_Iterate_GPU_best3(pdhg, nIter_restart);
    break;
  case 4:
    cupdlp_printf("case 4");
    PDTEST_Restart_Iterate_GPU_best4(pdhg, nIter_restart);
    break;
  case 5:
    cupdlp_printf("case 5");
    PDTEST_Restart_Iterate_GPU_best5(pdhg, nIter_restart);
    break;
  default:
    cupdlp_printf("Error: bestID = %d, 不在取值范围内", bestID);
    break;
  }
}

void PDTEST_Restart_Iterate_GPU_best1(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  CUPDLPproblem *problem = pdhg->problem;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPtimers *timers = pdhg->timers;

  PDHG_restart_choice restart_choice = PDTEST_Check_Restart_GPU_Only_Current(pdhg);

  if (restart_choice == PDHG_NO_RESTART)
    return;
  ///////////////////////////////////////////////////
  // 如果restart了，就把nIter_restart置为0
  // *nIter_restart = 0;
  // stepsize->dBeta_ag = 0.0;
  stepsize->dBeta_ag = 1.0;
  // CUPDLP_COPY_VEC(iterates->x->data, iterates->x_ag->data, cupdlp_float,
  //                 problem->nCols);
  // CUPDLP_COPY_VEC(iterates->y->data, iterates->y_ag->data, cupdlp_float,
  //                 problem->nRows);
  // CUPDLP_COPY_VEC(iterates->ax->data, iterates->ax_ag->data, cupdlp_float,
  //                 problem->nRows);
  // CUPDLP_COPY_VEC(iterates->aty->data, iterates->aty_ag->data,
  //                 cupdlp_float, problem->nCols);

  // CUPDLP_COPY_VEC(iterates->x_bar->data, iterates->x_ag->data, cupdlp_float,
  //                 problem->nCols);
  // CUPDLP_COPY_VEC(iterates->ax_bar->data, iterates->ax_ag->data, cupdlp_float,
  //                 problem->nRows);
  ///////////////////////////////////////////////////
  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;

  resobj->dPrimalFeasLastRestart = resobj->dPrimalFeas;
  resobj->dDualFeasLastRestart = resobj->dDualFeas;
  resobj->dDualityGapLastRestart = resobj->dDualityGap;

  ////////////////////////////////////////////////////////
  PDTEST_Compute_Step_Size_Ratio(pdhg);
  ////////////////////////////////////////////////////////

  CUPDLP_COPY_VEC(iterates->xLastRestart, iterates->x_ag->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->yLastRestart, iterates->y_ag->data, cupdlp_float,
                  problem->nRows);

  iterates->iLastRestartIter = timers->nIter;

  PDTEST_Compute_Residuals(pdhg);
  // cupdlp_printf("Recomputed stepsize ratio: %e,  sqrt(ratio)=%e",
  // stepsize->dBeta, sqrt(stepsize->dBeta));
}

void PDTEST_Restart_Iterate_GPU_best2(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  CUPDLPproblem *problem = pdhg->problem;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPtimers *timers = pdhg->timers;

  PDHG_restart_choice restart_choice = PDTEST_Check_Restart_GPU_Only_Current(pdhg);

  if (restart_choice == PDHG_NO_RESTART)
    return;
  ///////////////////////////////////////////////////
  // 如果restart了，就把nIter_restart置为0
  // *nIter_restart = 0;
  // stepsize->dBeta_ag = 0.0;
  // stepsize->dBeta_ag = 1.0;
  // CUPDLP_COPY_VEC(iterates->x->data, iterates->x_ag->data, cupdlp_float,
  //                 problem->nCols);
  // CUPDLP_COPY_VEC(iterates->y->data, iterates->y_ag->data, cupdlp_float,
  //                 problem->nRows);
  // CUPDLP_COPY_VEC(iterates->ax->data, iterates->ax_ag->data, cupdlp_float,
  //                 problem->nRows);
  // CUPDLP_COPY_VEC(iterates->aty->data, iterates->aty_ag->data,
  //                 cupdlp_float, problem->nCols);

  // CUPDLP_COPY_VEC(iterates->x_bar->data, iterates->x_ag->data, cupdlp_float,
  //                 problem->nCols);
  // CUPDLP_COPY_VEC(iterates->ax_bar->data, iterates->ax_ag->data, cupdlp_float,
  //                 problem->nRows);
  ///////////////////////////////////////////////////
  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;

  resobj->dPrimalFeasLastRestart = resobj->dPrimalFeas;
  resobj->dDualFeasLastRestart = resobj->dDualFeas;
  resobj->dDualityGapLastRestart = resobj->dDualityGap;

  ////////////////////////////////////////////////////////
  PDTEST_Compute_Step_Size_Ratio(pdhg);
  ////////////////////////////////////////////////////////

  CUPDLP_COPY_VEC(iterates->xLastRestart, iterates->x_ag->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->yLastRestart, iterates->y_ag->data, cupdlp_float,
                  problem->nRows);

  iterates->iLastRestartIter = timers->nIter;

  PDTEST_Compute_Residuals(pdhg);
  // cupdlp_printf("Recomputed stepsize ratio: %e,  sqrt(ratio)=%e",
  // stepsize->dBeta, sqrt(stepsize->dBeta));
}

void PDTEST_Restart_Iterate_GPU_best3(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  CUPDLPproblem *problem = pdhg->problem;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPtimers *timers = pdhg->timers;

  PDHG_restart_choice restart_choice = PDTEST_Check_Restart_GPU_Only_Current(pdhg);

  if (restart_choice == PDHG_NO_RESTART)
    return;
  ///////////////////////////////////////////////////
  // 如果restart了，就把nIter_restart置为0
  // *nIter_restart = 0;
  // stepsize->dBeta_ag = 0.0;
  // stepsize->dBeta_ag = 1.0;
  CUPDLP_COPY_VEC(iterates->x->data, iterates->x_ag->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->y->data, iterates->y_ag->data, cupdlp_float,
                  problem->nRows);
  CUPDLP_COPY_VEC(iterates->ax->data, iterates->ax_ag->data, cupdlp_float,
                  problem->nRows);
  CUPDLP_COPY_VEC(iterates->aty->data, iterates->aty_ag->data,
                  cupdlp_float, problem->nCols);

  CUPDLP_COPY_VEC(iterates->x_bar->data, iterates->x_ag->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->ax_bar->data, iterates->ax_ag->data, cupdlp_float,
                  problem->nRows);
  ///////////////////////////////////////////////////
  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;

  resobj->dPrimalFeasLastRestart = resobj->dPrimalFeas;
  resobj->dDualFeasLastRestart = resobj->dDualFeas;
  resobj->dDualityGapLastRestart = resobj->dDualityGap;

  ////////////////////////////////////////////////////////
  PDTEST_Compute_Step_Size_Ratio(pdhg);
  ////////////////////////////////////////////////////////

  CUPDLP_COPY_VEC(iterates->xLastRestart, iterates->x_ag->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->yLastRestart, iterates->y_ag->data, cupdlp_float,
                  problem->nRows);

  iterates->iLastRestartIter = timers->nIter;

  PDTEST_Compute_Residuals(pdhg);
  // cupdlp_printf("Recomputed stepsize ratio: %e,  sqrt(ratio)=%e",
  // stepsize->dBeta, sqrt(stepsize->dBeta));
}

void PDTEST_Restart_Iterate_GPU_best4(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  CUPDLPproblem *problem = pdhg->problem;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPtimers *timers = pdhg->timers;

  PDHG_restart_choice restart_choice = PDTEST_Check_Restart_GPU_Only_Current(pdhg);

  if (restart_choice == PDHG_NO_RESTART)
    return;
  ///////////////////////////////////////////////////
  // 如果restart了，就把nIter_restart置为0
  // *nIter_restart = 0;
  // stepsize->dBeta_ag = 0.0;
  stepsize->dBeta_ag = 1.0;
  // CUPDLP_COPY_VEC(iterates->x->data, iterates->x_ag->data, cupdlp_float,
  //                 problem->nCols);
  // CUPDLP_COPY_VEC(iterates->y->data, iterates->y_ag->data, cupdlp_float,
  //                 problem->nRows);
  // CUPDLP_COPY_VEC(iterates->ax->data, iterates->ax_ag->data, cupdlp_float,
  //                 problem->nRows);
  // CUPDLP_COPY_VEC(iterates->aty->data, iterates->aty_ag->data,
  //                 cupdlp_float, problem->nCols);

  // CUPDLP_COPY_VEC(iterates->x_bar->data, iterates->x_ag->data, cupdlp_float,
  //                 problem->nCols);
  // CUPDLP_COPY_VEC(iterates->ax_bar->data, iterates->ax_ag->data, cupdlp_float,
  //                 problem->nRows);
  ///////////////////////////////////////////////////
  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;

  resobj->dPrimalFeasLastRestart = resobj->dPrimalFeas;
  resobj->dDualFeasLastRestart = resobj->dDualFeas;
  resobj->dDualityGapLastRestart = resobj->dDualityGap;

  ////////////////////////////////////////////////////////
  PDTEST_Compute_Step_Size_Ratio(pdhg);
  ////////////////////////////////////////////////////////

  CUPDLP_COPY_VEC(iterates->xLastRestart, iterates->x_ag->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->yLastRestart, iterates->y_ag->data, cupdlp_float,
                  problem->nRows);

  iterates->iLastRestartIter = timers->nIter;

  PDTEST_Compute_Residuals(pdhg);
  // cupdlp_printf("Recomputed stepsize ratio: %e,  sqrt(ratio)=%e",
  // stepsize->dBeta, sqrt(stepsize->dBeta));
}

void PDTEST_Restart_Iterate_GPU_best5(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  CUPDLPproblem *problem = pdhg->problem;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPtimers *timers = pdhg->timers;

  PDHG_restart_choice restart_choice = PDTEST_Check_Restart_GPU_Only_Current(pdhg);

  if (restart_choice == PDHG_NO_RESTART)
    return;
  ///////////////////////////////////////////////////
  // 如果restart了，就把nIter_restart置为0
  // *nIter_restart = 0;
  // stepsize->dBeta_ag = 0.0;
  // stepsize->dBeta_ag = 1.0;
  // CUPDLP_COPY_VEC(iterates->x->data, iterates->x_ag->data, cupdlp_float,
  //                 problem->nCols);
  // CUPDLP_COPY_VEC(iterates->y->data, iterates->y_ag->data, cupdlp_float,
  //                 problem->nRows);
  // CUPDLP_COPY_VEC(iterates->ax->data, iterates->ax_ag->data, cupdlp_float,
  //                 problem->nRows);
  // CUPDLP_COPY_VEC(iterates->aty->data, iterates->aty_ag->data,
  //                 cupdlp_float, problem->nCols);

  // CUPDLP_COPY_VEC(iterates->x_bar->data, iterates->x_ag->data, cupdlp_float,
  //                 problem->nCols);
  // CUPDLP_COPY_VEC(iterates->ax_bar->data, iterates->ax_ag->data, cupdlp_float,
  //                 problem->nRows);
  ///////////////////////////////////////////////////
  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;

  resobj->dPrimalFeasLastRestart = resobj->dPrimalFeas;
  resobj->dDualFeasLastRestart = resobj->dDualFeas;
  resobj->dDualityGapLastRestart = resobj->dDualityGap;

  ////////////////////////////////////////////////////////
  PDTEST_Compute_Step_Size_Ratio(pdhg);
  ////////////////////////////////////////////////////////

  CUPDLP_COPY_VEC(iterates->xLastRestart, iterates->x_ag->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->yLastRestart, iterates->y_ag->data, cupdlp_float,
                  problem->nRows);

  iterates->iLastRestartIter = timers->nIter;

  PDTEST_Compute_Residuals(pdhg);
  // cupdlp_printf("Recomputed stepsize ratio: %e,  sqrt(ratio)=%e",
  // stepsize->dBeta, sqrt(stepsize->dBeta));
}
void PDTEST_Restart_Iterate_GPU_Only_VariableUpdate(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  CUPDLPproblem *problem = pdhg->problem;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPtimers *timers = pdhg->timers;

  PDHG_restart_choice restart_choice = PDTEST_Check_Restart_GPU_Only_Current(pdhg);

  if (restart_choice == PDHG_NO_RESTART)
    return;
  ///////////////////////////////////////////////////
  // 如果restart了，就把nIter_restart置为0
  // *nIter_restart = 0;
  // stepsize->dBeta_ag = 0.0;
  stepsize->dBeta_ag = 1.0;
  CUPDLP_COPY_VEC(iterates->x->data, iterates->x_ag->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->y->data, iterates->y_ag->data, cupdlp_float,
                  problem->nRows);
  CUPDLP_COPY_VEC(iterates->ax->data, iterates->ax_ag->data, cupdlp_float,
                  problem->nRows);
  CUPDLP_COPY_VEC(iterates->aty->data, iterates->aty_ag->data,
                  cupdlp_float, problem->nCols);

  CUPDLP_COPY_VEC(iterates->x_bar->data, iterates->x_ag->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->ax_bar->data, iterates->ax_ag->data, cupdlp_float,
                  problem->nRows);
  ///////////////////////////////////////////////////
  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;

  resobj->dPrimalFeasLastRestart = resobj->dPrimalFeas;
  resobj->dDualFeasLastRestart = resobj->dDualFeas;
  resobj->dDualityGapLastRestart = resobj->dDualityGap;

  ////////////////////////////////////////////////////////
  PDTEST_Compute_Step_Size_Ratio(pdhg);
  ////////////////////////////////////////////////////////

  CUPDLP_COPY_VEC(iterates->xLastRestart, iterates->x_ag->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->yLastRestart, iterates->y_ag->data, cupdlp_float,
                  problem->nRows);

  iterates->iLastRestartIter = timers->nIter;

  PDTEST_Compute_Residuals(pdhg);
  // cupdlp_printf("Recomputed stepsize ratio: %e,  sqrt(ratio)=%e",
  // stepsize->dBeta, sqrt(stepsize->dBeta));
}

void PDTEST_Restart_Iterate_GPU_Only_WeightUpdate(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  CUPDLPproblem *problem = pdhg->problem;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPtimers *timers = pdhg->timers;

  PDHG_restart_choice restart_choice = PDTEST_Check_Restart_GPU_Only_Current(pdhg);

  if (restart_choice == PDHG_NO_RESTART)
    return;
  ///////////////////////////////////////////////////
  // 如果restart了，就把nIter_restart置为0
  // *nIter_restart = 0;
  // stepsize->dBeta_ag = 0.0;
  stepsize->dBeta_ag = 1.0;
  CUPDLP_COPY_VEC(iterates->x->data, iterates->x_ag->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->y->data, iterates->y_ag->data, cupdlp_float,
                  problem->nRows);
  CUPDLP_COPY_VEC(iterates->ax->data, iterates->ax_ag->data, cupdlp_float,
                  problem->nRows);
  CUPDLP_COPY_VEC(iterates->aty->data, iterates->aty_ag->data,
                  cupdlp_float, problem->nCols);

  CUPDLP_COPY_VEC(iterates->x_bar->data, iterates->x_ag->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->ax_bar->data, iterates->ax_ag->data, cupdlp_float,
                  problem->nRows);
  ///////////////////////////////////////////////////
  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;

  resobj->dPrimalFeasLastRestart = resobj->dPrimalFeas;
  resobj->dDualFeasLastRestart = resobj->dDualFeas;
  resobj->dDualityGapLastRestart = resobj->dDualityGap;

  ////////////////////////////////////////////////////////
  PDTEST_Compute_Step_Size_Ratio(pdhg);
  ////////////////////////////////////////////////////////

  CUPDLP_COPY_VEC(iterates->xLastRestart, iterates->x_ag->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->yLastRestart, iterates->y_ag->data, cupdlp_float,
                  problem->nRows);

  iterates->iLastRestartIter = timers->nIter;

  PDTEST_Compute_Residuals(pdhg);
  // cupdlp_printf("Recomputed stepsize ratio: %e,  sqrt(ratio)=%e",
  // stepsize->dBeta, sqrt(stepsize->dBeta));
}