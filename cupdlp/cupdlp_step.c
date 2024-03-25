//
// Created by chuwen on 23-11-28.
//

#include "cupdlp_step.h"

#include "cupdlp_defs.h"
#include "cupdlp_linalg.h"
#include "cupdlp_proj.h"
// #include "cupdlp_scaling.h"
#include "cupdlp_solver.h"
#include "cupdlp_utils.h"
#include "glbopts.h"

void PDTEST_x_md_step(CUPDLPwork *work, cupdlp_float beta)
{
  PDTESTiterates *iterates = work->PDTESTiterates;
  CUPDLPproblem *problem = work->problem;

#if !(CUPDLP_CPU) & USE_KERNELS
  pdtest_x_md_update_cuda(iterates->x_md->data, iterates->x_ag->data, iterates->x->data, beta, problem->nCols);
#else
  cupdlp_printf("PDTEST_x_md_step is not implemented\n");
#endif
}

void PDTEST_x_ag_step(CUPDLPwork *work, cupdlp_float beta)
{
  PDTESTiterates *iterates = work->PDTESTiterates;
  CUPDLPproblem *problem = work->problem;

#if !(CUPDLP_CPU) & USE_KERNELS
  pdtest_x_ag_update_cuda(iterates->x_agUpdate->data, iterates->x_ag->data, iterates->xUpdate->data, beta, problem->nCols);
#else
  cupdlp_printf("PDTEST_x_ag_step is not implemented\n");
#endif
}

void PDTEST_y_ag_step(CUPDLPwork *work, cupdlp_float beta)
{
  PDTESTiterates *iterates = work->PDTESTiterates;
  CUPDLPproblem *problem = work->problem;

#if !(CUPDLP_CPU) & USE_KERNELS
  pdtest_y_ag_update_cuda(iterates->y_agUpdate->data, iterates->y_ag->data, iterates->yUpdate->data, beta, problem->nRows);
#else
  cupdlp_printf("PDTEST_y_ag_step is not implemented\n");
#endif
}

void PDTEST_x_bar_step(CUPDLPwork *work, cupdlp_float theta)
{
  PDTESTiterates *iterates = work->PDTESTiterates;
  CUPDLPproblem *problem = work->problem;

#if !(CUPDLP_CPU) & USE_KERNELS
  pdtest_x_bar_update_cuda(iterates->x_barUpdate->data, iterates->xUpdate->data, iterates->x->data, theta, problem->nCols);
#else
  cupdlp_printf("PDTEST_x_bar_step is not implemented\n");
#endif
}

void PDTEST_primalGradientStep(CUPDLPwork *work, cupdlp_float dPrimalStepSize)
{
  PDTESTiterates *iterates = work->PDTESTiterates;
  CUPDLPproblem *problem = work->problem;

#if !(CUPDLP_CPU) & USE_KERNELS
  pdtest_pgrad_cuda(iterates->xUpdate->data, iterates->x->data, problem->cost,
                    iterates->atyUpdate->data, dPrimalStepSize, problem->nCols);
#else

  // cupdlp_copy(iterates->xUpdate, iterates->x, cupdlp_float, problem->nCols);
  CUPDLP_COPY_VEC(iterates->xUpdate->data, iterates->x->data, cupdlp_float,
                  problem->nCols);

  // AddToVector(iterates->xUpdate, -dPrimalStepSize, problem->cost,
  // problem->nCols); AddToVector(iterates->xUpdate, dPrimalStepSize,
  // iterates->aty, problem->nCols);

  cupdlp_float alpha = -dPrimalStepSize;
  cupdlp_axpy(work, problem->nCols, &alpha, problem->cost,
              iterates->xUpdate->data);
  alpha = dPrimalStepSize;
  cupdlp_axpy(work, problem->nCols, &alpha, iterates->aty->data,
              iterates->xUpdate->data);
#endif
}

void PDTEST_dualGradientStep(CUPDLPwork *work, cupdlp_float dDualStepSize)
{
  PDTESTiterates *iterates = work->PDTESTiterates;
  CUPDLPproblem *problem = work->problem;

#if !(CUPDLP_CPU) & USE_KERNELS
  pdtest_dgrad_cuda(iterates->yUpdate->data, iterates->y->data, problem->rhs,
                    iterates->ax_bar->data, dDualStepSize,
                    problem->nRows);
#else

  // cupdlp_copy(iterates->yUpdate, iterates->y, cupdlp_float, problem->nRows);
  CUPDLP_COPY_VEC(iterates->yUpdate->data, iterates->y->data, cupdlp_float,
                  problem->nRows);

  // AddToVector(iterates->yUpdate, dDualStepSize, problem->rhs,
  // problem->nRows); AddToVector(iterates->yUpdate, -2.0 * dDualStepSize,
  // iterates->axUpdate, problem->nRows); AddToVector(iterates->yUpdate,
  // dDualStepSize, iterates->ax, problem->nRows);

  cupdlp_float alpha = dDualStepSize;
  cupdlp_axpy(work, problem->nRows, &alpha, problem->rhs,
              iterates->yUpdate->data);
  alpha = -2.0 * dDualStepSize;
  cupdlp_axpy(work, problem->nRows, &alpha, iterates->axUpdate->data,
              iterates->yUpdate->data);
  alpha = dDualStepSize;
  cupdlp_axpy(work, problem->nRows, &alpha, iterates->ax->data,
              iterates->yUpdate->data);
#endif
}

// xUpdate = x^k - dPrimalStep * (c - A'y^k)
void PDHG_primalGradientStep(CUPDLPwork *work, cupdlp_float dPrimalStepSize)
{
  CUPDLPiterates *iterates = work->iterates;
  CUPDLPproblem *problem = work->problem;

#if !(CUPDLP_CPU) & USE_KERNELS
  cupdlp_pgrad_cuda(iterates->xUpdate->data, iterates->x->data, problem->cost,
                    iterates->aty->data, dPrimalStepSize, problem->nCols);
#else

  // cupdlp_copy(iterates->xUpdate, iterates->x, cupdlp_float, problem->nCols);
  CUPDLP_COPY_VEC(iterates->xUpdate->data, iterates->x->data, cupdlp_float,
                  problem->nCols);

  // AddToVector(iterates->xUpdate, -dPrimalStepSize, problem->cost,
  // problem->nCols); AddToVector(iterates->xUpdate, dPrimalStepSize,
  // iterates->aty, problem->nCols);

  cupdlp_float alpha = -dPrimalStepSize;
  cupdlp_axpy(work, problem->nCols, &alpha, problem->cost,
              iterates->xUpdate->data);
  alpha = dPrimalStepSize;
  cupdlp_axpy(work, problem->nCols, &alpha, iterates->aty->data,
              iterates->xUpdate->data);
#endif
}

// yUpdate = y^k + dDualStep * (b - A * (2x^{k+1} - x^{k})
void PDHG_dualGradientStep(CUPDLPwork *work, cupdlp_float dDualStepSize)
{
  CUPDLPiterates *iterates = work->iterates;
  CUPDLPproblem *problem = work->problem;

#if !(CUPDLP_CPU) & USE_KERNELS
  cupdlp_dgrad_cuda(iterates->yUpdate->data, iterates->y->data, problem->rhs,
                    iterates->ax->data, iterates->axUpdate->data, dDualStepSize,
                    problem->nRows);
#else

  // cupdlp_copy(iterates->yUpdate, iterates->y, cupdlp_float, problem->nRows);
  CUPDLP_COPY_VEC(iterates->yUpdate->data, iterates->y->data, cupdlp_float,
                  problem->nRows);

  // AddToVector(iterates->yUpdate, dDualStepSize, problem->rhs,
  // problem->nRows); AddToVector(iterates->yUpdate, -2.0 * dDualStepSize,
  // iterates->axUpdate, problem->nRows); AddToVector(iterates->yUpdate,
  // dDualStepSize, iterates->ax, problem->nRows);

  cupdlp_float alpha = dDualStepSize;
  cupdlp_axpy(work, problem->nRows, &alpha, problem->rhs,
              iterates->yUpdate->data);
  alpha = -2.0 * dDualStepSize;
  cupdlp_axpy(work, problem->nRows, &alpha, iterates->axUpdate->data,
              iterates->yUpdate->data);
  alpha = dDualStepSize;
  cupdlp_axpy(work, problem->nRows, &alpha, iterates->ax->data,
              iterates->yUpdate->data);
#endif
}

cupdlp_retcode PDHG_Power_Method(CUPDLPwork *work, cupdlp_float *lambda)
{
  cupdlp_retcode retcode = RETCODE_OK;
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPiterates *iterates = work->iterates;

  cupdlp_printf("Power Method:\n");

  cupdlp_float *q = work->buffer->data;

  cupdlp_initvec(q, 1.0, lp->nRows);

  double res = 0.0;
  for (cupdlp_int iter = 0; iter < 20; ++iter)
  {
    // z = A*A'*q
    ATy(work, iterates->aty, work->buffer);
    Ax(work, iterates->ax, iterates->aty);

    // q = z / norm(z)
    CUPDLP_COPY_VEC(q, iterates->ax->data, cupdlp_float, lp->nRows);
    cupdlp_float qNorm = 0.0;
    cupdlp_twoNorm(work, lp->nRows, q, &qNorm);
    cupdlp_scaleVector(work, 1.0 / qNorm, q, lp->nRows);

    ATy(work, iterates->aty, work->buffer);

    cupdlp_twoNormSquared(work, lp->nCols, iterates->aty->data, lambda);

    cupdlp_float alpha = -(*lambda);
    cupdlp_axpy(work, lp->nRows, &alpha, q, iterates->ax->data);

    cupdlp_twoNormSquared(work, lp->nCols, iterates->ax->data, &res);

    cupdlp_printf("% d  %e  %.3f\n", iter, *lambda, res);
  }

exit_cleanup:
  return retcode;
}
cupdlp_retcode PDTEST_Power_Method(CUPDLPwork *work, cupdlp_float *lambda)
{
  cupdlp_retcode retcode = RETCODE_OK;
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  PDTESTiterates *iterates = work->PDTESTiterates;

  cupdlp_printf("Power Method:\n");

  cupdlp_float *q = work->buffer->data;

  cupdlp_initvec(q, 1.0, lp->nRows);

  double res = 0.0;
  for (cupdlp_int iter = 0; iter < 20; ++iter)
  {
    // z = A*A'*q
    ATy(work, iterates->aty, work->buffer);
    Ax(work, iterates->ax, iterates->aty);

    // q = z / norm(z)
    CUPDLP_COPY_VEC(q, iterates->ax->data, cupdlp_float, lp->nRows);
    cupdlp_float qNorm = 0.0;
    cupdlp_twoNorm(work, lp->nRows, q, &qNorm);
    cupdlp_scaleVector(work, 1.0 / qNorm, q, lp->nRows);

    ATy(work, iterates->aty, work->buffer);

    cupdlp_twoNormSquared(work, lp->nCols, iterates->aty->data, lambda);

    cupdlp_float alpha = -(*lambda);
    cupdlp_axpy(work, lp->nRows, &alpha, q, iterates->ax->data);

    cupdlp_twoNormSquared(work, lp->nCols, iterates->ax->data, &res);

    cupdlp_printf("% d  %e  %.3f\n", iter, *lambda, res);
  }

exit_cleanup:
  return retcode;
}

// Primal Weight Update
void PDHG_Compute_Step_Size_Ratio(CUPDLPwork *pdhg)
{
  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  cupdlp_float dMeanStepSize =
      sqrt(stepsize->dPrimalStep * stepsize->dDualStep);

  // cupdlp_float dDiffPrimal = cupdlp_diffTwoNorm(iterates->x,
  // iterates->xLastRestart, problem->nCols); cupdlp_float dDiffDual =
  // cupdlp_diffTwoNorm(iterates->y, iterates->yLastRestart, problem->nRows);

  cupdlp_float dDiffPrimal = 0.0;
  cupdlp_diffTwoNorm(pdhg, iterates->x->data, iterates->xLastRestart,
                     problem->nCols, &dDiffPrimal);
  cupdlp_float dDiffDual = 0.0;
  cupdlp_diffTwoNorm(pdhg, iterates->y->data, iterates->yLastRestart,
                     problem->nRows, &dDiffDual);

  if (fmin(dDiffPrimal, dDiffDual) > 1e-10)
  {
    cupdlp_float dBetaUpdate = dDiffDual / dDiffPrimal;
    cupdlp_float dLogBetaUpdate =
        0.5 * log(dBetaUpdate) + 0.5 * log(sqrt(stepsize->dBeta));
    stepsize->dBeta = exp(dLogBetaUpdate) * exp(dLogBetaUpdate);
  }

  stepsize->dPrimalStep = dMeanStepSize / sqrt(stepsize->dBeta);
  stepsize->dDualStep = stepsize->dPrimalStep * stepsize->dBeta;
  stepsize->dTheta = 1.0;
}

void PDTEST_Compute_Step_Size_Ratio(CUPDLPwork *pdhg)
{
  CUPDLPproblem *problem = pdhg->problem;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  cupdlp_float dMeanStepSize =
      sqrt(stepsize->dPrimalStep * stepsize->dDualStep);

  // cupdlp_float dDiffPrimal = cupdlp_diffTwoNorm(iterates->x,
  // iterates->xLastRestart, problem->nCols); cupdlp_float dDiffDual =
  // cupdlp_diffTwoNorm(iterates->y, iterates->yLastRestart, problem->nRows);

  cupdlp_float dDiffPrimal = 0.0;
  cupdlp_diffTwoNorm(pdhg, iterates->x_ag->data, iterates->xLastRestart,
                     problem->nCols, &dDiffPrimal);
  cupdlp_float dDiffDual = 0.0;
  cupdlp_diffTwoNorm(pdhg, iterates->y_ag->data, iterates->yLastRestart,
                     problem->nRows, &dDiffDual);

  if (fmin(dDiffPrimal, dDiffDual) > 1e-10)
  {
    cupdlp_float dBetaUpdate = dDiffDual / dDiffPrimal;
    cupdlp_float dLogBetaUpdate =
        0.5 * log(dBetaUpdate) + 0.5 * log(sqrt(stepsize->dBeta));
    stepsize->dBeta = exp(dLogBetaUpdate) * exp(dLogBetaUpdate);
  }

  stepsize->dPrimalStep = dMeanStepSize / sqrt(stepsize->dBeta);
  stepsize->dDualStep = stepsize->dPrimalStep * stepsize->dBeta;
  stepsize->dTheta = 1.0;
}

void PDHG_Update_Iterate_Constant_Step_Size(CUPDLPwork *pdhg)
{
  //            CUPDLP_ASSERT(0);
  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  // Ax(pdhg, iterates->ax, iterates->x);
  // ATyCPU(pdhg, iterates->aty, iterates->y);
  Ax(pdhg, iterates->ax, iterates->x);
  ATy(pdhg, iterates->aty, iterates->y);

  // x^{k+1} = proj_{X}(x^k - dPrimalStep * (c - A'y^k))
  PDHG_primalGradientStep(pdhg, stepsize->dPrimalStep);

  PDHG_Project_Bounds(pdhg, iterates->xUpdate->data);
  // Ax(pdhg, iterates->axUpdate, iterates->xUpdate);
  Ax(pdhg, iterates->axUpdate, iterates->xUpdate);

  // y^{k+1} = y^k + dDualStep * (b - A * (2x^{k+1} - x^{k})
  PDHG_dualGradientStep(pdhg, stepsize->dDualStep);

  PDHG_Project_Row_Duals(pdhg, iterates->yUpdate->data);
  // ATyCPU(pdhg, iterates->atyUpdate, iterates->yUpdate);
  ATy(pdhg, iterates->atyUpdate, iterates->yUpdate);
}
void PDTEST_Update_Iterate_Constant_Step_Size(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  //            CUPDLP_ASSERT(0);
  CUPDLPproblem *problem = pdhg->problem;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPtimers *timers = pdhg->timers; // 各种耗时
  // dStepSizeUpdate是论文中的eta
  // cupdlp_float dStepSize =
  //     sqrt(stepsize->dPrimalStep * stepsize->dDualStep);

  // stepsize->dPrimalStep = dStepSize;
  // stepsize->dDualStep = dStepSize;
  cupdlp_float dMultiStartTime;
  cupdlp_float dAddStartTime;
  cupdlp_int t_count = *nIter_restart + 1;
  // cupdlp_int t_count = timers->nIter + 1;
  cupdlp_printf("t_count: %d\n", t_count);
  // cupdlp_float beta = (t_count + 1) / 2.0;
  cupdlp_float beta = (t_count + 3.8) / 4.8;
  // cupdlp_float beta = t_count + 0.0;
  cupdlp_float dIterTime = getTimeStamp();

  // 没有必要计算x_md, 因为梯度是常数，用不到x_md
  // x_md^{t} = (1 - 1 / beta^t)x_ag^{t} + (1 / beta^t)x^{t}
  // PDTEST_x_md_step(pdhg, beta);
  // 没有必要进行Project，因为都是已经投影过的x进行线性组合
  // PDHG_Project_Bounds(pdhg, iterates->x_md->data);

  // 计算Ax_bar^{t}, 后面计算y^{t+1}有用
  dMultiStartTime = getTimeStamp();
  Ax(pdhg, iterates->ax_bar, iterates->x_bar);
  timers->dMatVecMultiplyTime += getTimeStamp() - dMultiStartTime;

  // y^{t+1} = y^t + dDualStep * (b - A * (x_bar^{t})
  dAddStartTime = getTimeStamp();
  PDTEST_dualGradientStep(pdhg, stepsize->dDualStep);
  timers->dVecVecAddTime += getTimeStamp() - dAddStartTime;

  PDHG_Project_Row_Duals(pdhg, iterates->yUpdate->data);

  // 计算ATy^{t+1}, 后面计算x^{t+1}有用
  dMultiStartTime = getTimeStamp();
  ATy(pdhg, iterates->atyUpdate, iterates->yUpdate);
  timers->dMatVecMultiplyTime += getTimeStamp() - dMultiStartTime;

  // x^{t+1} = proj_{X}(x^t - dPrimalStep * (c - A'y^{t+1}))
  dAddStartTime = getTimeStamp();
  PDTEST_primalGradientStep(pdhg, stepsize->dPrimalStep);
  timers->dVecVecAddTime += getTimeStamp() - dAddStartTime;
  PDHG_Project_Bounds(pdhg, iterates->xUpdate->data);

  // x_ag^{t+1} = (1 - 1 / beta^t)x_ag^{t} + (1 / beta^t)x^{t}
  dAddStartTime = getTimeStamp();
  PDTEST_x_ag_step(pdhg, beta);
  timers->dVecVecAddTime += getTimeStamp() - dAddStartTime;
  // 没有必要进行Project，因为都是已经投影过的x进行线性组合
  // PDHG_Project_Bounds(pdhg, iterates->x_agUpdate->data);

  // y_ag^{t+1} = (1 - 1 / beta^t)y_ag^{t} + (1 / beta^t)y^{t}
  dAddStartTime = getTimeStamp();
  PDTEST_y_ag_step(pdhg, beta);
  timers->dVecVecAddTime += getTimeStamp() - dAddStartTime;

  // theta^{t+1}
  // cupdlp_float theta = t_count / (t_count + 1.0);
  cupdlp_float theta = 1.0;

  // x_bar^{t+1} = theta^{t+1}(x^{t+1} - x^{t}) + x^{t+1}
  dAddStartTime = getTimeStamp();
  PDTEST_x_bar_step(pdhg, theta);
  timers->dVecVecAddTime += getTimeStamp() - dAddStartTime;

  // 更新一下ax_ag和aty_ag
  dMultiStartTime = getTimeStamp();
  Ax(pdhg, iterates->ax_agUpdate, iterates->x_agUpdate);
  ATy(pdhg, iterates->aty_agUpdate, iterates->y_agUpdate);
  timers->dMatVecMultiplyTime += getTimeStamp() - dMultiStartTime;

  timers->dIterTime += getTimeStamp() - dIterTime;
}
void PDHG_Update_Iterate_Malitsky_Pock(CUPDLPwork *pdhg)
{
  cupdlp_printf("Malitsky-Pock is not implemented\n");
  cupdlp_printf(" - use %d and %d instead", PDHG_FIXED_LINESEARCH,
                PDHG_ADAPTIVE_LINESEARCH);
  exit(-1);
}

cupdlp_retcode PDHG_Update_Iterate_Adaptive_Step_Size(CUPDLPwork *pdhg)
{
  cupdlp_retcode retcode = RETCODE_OK;
  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;

  cupdlp_float dStepSizeUpdate =
      sqrt(stepsize->dPrimalStep * stepsize->dDualStep);

  cupdlp_bool isDone = false;
  // number of steps this round
  int stepIterThis = 0;
  while (!isDone)
  {
    ++stepsize->nStepSizeIter;
    ++stepIterThis;

    cupdlp_float dPrimalStepUpdate = dStepSizeUpdate / sqrt(stepsize->dBeta);
    cupdlp_float dDualStepUpdate = dStepSizeUpdate * sqrt(stepsize->dBeta);

#pragma region Update
    // x^{k+1} = proj_{X}(x^k - dPrimalStep * (cupdlp - A'y^k))
    PDHG_primalGradientStep(pdhg, dPrimalStepUpdate);

    PDHG_Project_Bounds(pdhg, iterates->xUpdate->data);
    Ax(pdhg, iterates->axUpdate, iterates->xUpdate);

    // y^{k+1} = proj_{Y}(y^k + dDualStep * (b - A * (2 * x^{k+1} - x^{k})))
    PDHG_dualGradientStep(pdhg, dDualStepUpdate);

    PDHG_Project_Row_Duals(pdhg, iterates->yUpdate->data);
    ATy(pdhg, iterates->atyUpdate, iterates->yUpdate);
#pragma endregion

    // dMovement是\|z_{t+1}-z_t\|_2^2/2, 分子
    // dInteraction是分母
    cupdlp_float dMovement = 0.0;
    cupdlp_float dInteraction = 0.0;

#if !(CUPDLP_CPU) & USE_KERNELS
    cupdlp_compute_interaction_and_movement(pdhg, &dMovement, &dInteraction);
#else
    cupdlp_float dX = 0.0;
    cupdlp_diffTwoNormSquared(pdhg, iterates->x->data, iterates->xUpdate->data,
                              problem->nCols, &dX);
    dX *= 0.5 * sqrt(stepsize->dBeta);

    cupdlp_float dY = 0.0;
    cupdlp_diffTwoNormSquared(pdhg, iterates->y->data, iterates->yUpdate->data,
                              problem->nRows, &dY);
    dY /= 2.0 * sqrt(stepsize->dBeta);
    dMovement = dX + dY;

    //      Δx' (AΔy)
    cupdlp_diffDotDiff(pdhg, iterates->x->data, iterates->xUpdate->data,
                       iterates->aty->data, iterates->atyUpdate->data,
                       problem->nCols, &dInteraction);
#endif

#if CUPDLP_DUMP_LINESEARCH_STATS & CUPDLP_DEBUG
    cupdlp_float dInteractiony = 0.0;
    //      Δy' (AΔx)
    cupdlp_diffDotDiff(pdhg, iterates->y->data, iterates->yUpdate->data,
                       iterates->ax->data, iterates->axUpdate->data,
                       problem->nRows, &dInteractiony);
#endif
    // 计算dStepSizeLimit
    cupdlp_float dStepSizeLimit;
    if (dInteraction != 0.0)
    {
      dStepSizeLimit = dMovement / fabs(dInteraction);
    }
    else
    {
      dStepSizeLimit = INFINITY;
    }

    if (dStepSizeUpdate <= dStepSizeLimit)
    {
      isDone = true;
      // break;
    }
    else
    {
      CUPDLP_CHECK_TIMEOUT(pdhg);
    }

    cupdlp_float dFirstTerm = (1.0 - pow(stepsize->nStepSizeIter + 1.0,
                                         -PDHG_STEPSIZE_REDUCTION_EXP)) *
                              dStepSizeLimit;
    cupdlp_float dSecondTerm =
        (1.0 + pow(stepsize->nStepSizeIter + 1.0, -PDHG_STEPSIZE_GROWTH_EXP)) *
        dStepSizeUpdate;
    dStepSizeUpdate = fmin(dFirstTerm, dSecondTerm);
#if CUPDLP_DUMP_LINESEARCH_STATS & CUPDLP_DEBUG
    cupdlp_printf(" -- stepsize iteration %d: %f %f\n", stepIterThis,
                  dStepSizeUpdate, dStepSizeLimit);

    cupdlp_printf(" -- PrimalStep DualStep: %f %f\n", stepsize->dPrimalStep,
                  stepsize->dDualStep);
    cupdlp_printf(" -- FirstTerm SecondTerm: %f %f\n", dFirstTerm, dSecondTerm);
    cupdlp_printf(" -- nStepSizeIter: %d\n", stepsize->nStepSizeIter);
    cupdlp_printf(" -- RED_EXP GRO_EXP: %f %f\n", PDHG_STEPSIZE_REDUCTION_EXP,
                  PDHG_STEPSIZE_GROWTH_EXP);

    cupdlp_printf("     -- iteraction(x) interaction(y): %f %f\n", dInteraction,
                  dInteractiony);
    cupdlp_printf("     -- movement (scaled norm)  : %f\n", dMovement);
    cupdlp_printf("     -- movement (scaled norm)  : %f\n", dMovement);
    if (stepIterThis > 200)
      break; // avoid unlimited runs due to bugs.
#endif
  }

  stepsize->dPrimalStep = dStepSizeUpdate / sqrt(stepsize->dBeta);
  stepsize->dDualStep = dStepSizeUpdate * sqrt(stepsize->dBeta);

exit_cleanup:
  return retcode;
}
cupdlp_retcode PDTEST_Update_Iterate_Adaptive_Step_Size_ag(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  cupdlp_retcode retcode = RETCODE_OK;
  CUPDLPproblem *problem = pdhg->problem;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPtimers *timers = pdhg->timers; // 各种耗时

  cupdlp_float dStepSizeUpdate =
      sqrt(stepsize->dPrimalStep * stepsize->dDualStep);

  cupdlp_bool isDone = false;
  // number of steps this round
  int stepIterThis = 0;
  while (!isDone)
  {
    ++stepsize->nStepSizeIter;
    ++stepIterThis;

    cupdlp_float dPrimalStepUpdate = dStepSizeUpdate / sqrt(stepsize->dBeta);
    cupdlp_float dDualStepUpdate = dStepSizeUpdate * sqrt(stepsize->dBeta);
#pragma region Update
    cupdlp_float dMultiStartTime;
    cupdlp_float dAddStartTime;
    cupdlp_int t_count = *nIter_restart + 1;
    cupdlp_printf("t_count: %d\n", t_count);
    cupdlp_float beta = (t_count + 1) / 2.0;
    cupdlp_float dIterTime = getTimeStamp();

    // 没有必要计算x_md, 因为梯度是常数，用不到x_md
    // x_md^{t} = (1 - 1 / beta^t)x_ag^{t} + (1 / beta^t)x^{t}
    // PDTEST_x_md_step(pdhg, beta);
    // 没有必要进行Project，因为都是已经投影过的x进行线性组合
    // PDHG_Project_Bounds(pdhg, iterates->x_md->data);

    // 计算Ax_bar^{t}, 后面计算y^{t+1}有用
    dMultiStartTime = getTimeStamp();
    Ax(pdhg, iterates->ax_bar, iterates->x_bar);
    timers->dMatVecMultiplyTime += getTimeStamp() - dMultiStartTime;

    // y^{t+1} = y^t + dDualStep * (b - A * (x_bar^{t})
    dAddStartTime = getTimeStamp();
    PDTEST_dualGradientStep(pdhg, stepsize->dDualStep);
    timers->dVecVecAddTime += getTimeStamp() - dAddStartTime;

    PDHG_Project_Row_Duals(pdhg, iterates->yUpdate->data);

    // 计算ATy^{t+1}, 后面计算x^{t+1}有用
    dMultiStartTime = getTimeStamp();
    ATy(pdhg, iterates->atyUpdate, iterates->yUpdate);
    timers->dMatVecMultiplyTime += getTimeStamp() - dMultiStartTime;

    // x^{t+1} = proj_{X}(x^t - dPrimalStep * (c - A'y^{t+1}))
    dAddStartTime = getTimeStamp();
    PDTEST_primalGradientStep(pdhg, stepsize->dPrimalStep);
    timers->dVecVecAddTime += getTimeStamp() - dAddStartTime;
    PDHG_Project_Bounds(pdhg, iterates->xUpdate->data);

    // x_ag^{t+1} = (1 - 1 / beta^t)x_ag^{t} + (1 / beta^t)x^{t}
    dAddStartTime = getTimeStamp();
    PDTEST_x_ag_step(pdhg, beta);
    timers->dVecVecAddTime += getTimeStamp() - dAddStartTime;
    // 没有必要进行Project，因为都是已经投影过的x进行线性组合
    // PDHG_Project_Bounds(pdhg, iterates->x_agUpdate->data);

    // y_ag^{t+1} = (1 - 1 / beta^t)y_ag^{t} + (1 / beta^t)y^{t}
    dAddStartTime = getTimeStamp();
    PDTEST_y_ag_step(pdhg, beta);
    timers->dVecVecAddTime += getTimeStamp() - dAddStartTime;

    // theta^{t+1}
    cupdlp_float theta = t_count / (t_count + 1.0);

    // x_bar^{t+1} = theta^{t+1}(x^{t+1} - x^{t}) + x^{t+1}
    dAddStartTime = getTimeStamp();
    PDTEST_x_bar_step(pdhg, theta);
    timers->dVecVecAddTime += getTimeStamp() - dAddStartTime;

    // 更新一下ax_ag和aty_ag
    dMultiStartTime = getTimeStamp();
    Ax(pdhg, iterates->ax_agUpdate, iterates->x_agUpdate);
    ATy(pdhg, iterates->aty_agUpdate, iterates->y_agUpdate);
    timers->dMatVecMultiplyTime += getTimeStamp() - dMultiStartTime;

    timers->dIterTime += getTimeStamp() - dIterTime;
#pragma endregion

    cupdlp_float dMovement = 0.0;
    cupdlp_float dInteraction = 0.0;

#if !(CUPDLP_CPU) & USE_KERNELS
    PDTEST_compute_interaction_and_movement(pdhg, &dMovement, &dInteraction);
#else
    cupdlp_float dX = 0.0;
    cupdlp_diffTwoNormSquared(pdhg, iterates->x->data, iterates->xUpdate->data,
                              problem->nCols, &dX);
    dX *= 0.5 * sqrt(stepsize->dBeta);

    cupdlp_float dY = 0.0;
    cupdlp_diffTwoNormSquared(pdhg, iterates->y->data, iterates->yUpdate->data,
                              problem->nRows, &dY);
    dY /= 2.0 * sqrt(stepsize->dBeta);
    dMovement = dX + dY;

    //      Δx' (AΔy)
    cupdlp_diffDotDiff(pdhg, iterates->x->data, iterates->xUpdate->data,
                       iterates->aty->data, iterates->atyUpdate->data,
                       problem->nCols, &dInteraction);
#endif

#if CUPDLP_DUMP_LINESEARCH_STATS & CUPDLP_DEBUG
    cupdlp_float dInteractiony = 0.0;
    //      Δy' (AΔx)
    cupdlp_diffDotDiff(pdhg, iterates->y->data, iterates->yUpdate->data,
                       iterates->ax->data, iterates->axUpdate->data,
                       problem->nRows, &dInteractiony);
#endif

    cupdlp_float dStepSizeLimit;
    if (dInteraction != 0.0)
    {
      dStepSizeLimit = dMovement / fabs(dInteraction);
    }
    else
    {
      dStepSizeLimit = INFINITY;
    }
    if (dStepSizeUpdate <= dStepSizeLimit)
    {
      isDone = true;
      // break;
    }
    else
    {
      CUPDLP_CHECK_TIMEOUT(pdhg);
    }

    cupdlp_float dFirstTerm = (1.0 - pow(stepsize->nStepSizeIter + 1.0,
                                         -PDHG_STEPSIZE_REDUCTION_EXP)) *
                              dStepSizeLimit;
    cupdlp_float dSecondTerm =
        (1.0 + pow(stepsize->nStepSizeIter + 1.0, -PDHG_STEPSIZE_GROWTH_EXP)) *
        dStepSizeUpdate;
    dStepSizeUpdate = fmin(dFirstTerm, dSecondTerm);
#if CUPDLP_DUMP_LINESEARCH_STATS & CUPDLP_DEBUG
    cupdlp_printf(" -- stepsize iteration %d: %f %f\n", stepIterThis,
                  dStepSizeUpdate, dStepSizeLimit);

    cupdlp_printf(" -- PrimalStep DualStep: %f %f\n", stepsize->dPrimalStep,
                  stepsize->dDualStep);
    cupdlp_printf(" -- FirstTerm SecondTerm: %f %f\n", dFirstTerm, dSecondTerm);
    cupdlp_printf(" -- nStepSizeIter: %d\n", stepsize->nStepSizeIter);
    cupdlp_printf(" -- RED_EXP GRO_EXP: %f %f\n", PDHG_STEPSIZE_REDUCTION_EXP,
                  PDHG_STEPSIZE_GROWTH_EXP);

    cupdlp_printf("     -- iteraction(x) interaction(y): %f %f\n", dInteraction,
                  dInteractiony);
    cupdlp_printf("     -- movement (scaled norm)  : %f\n", dMovement);
    cupdlp_printf("     -- movement (scaled norm)  : %f\n", dMovement);
    if (stepIterThis > 200)
      break; // avoid unlimited runs due to bugs.
#endif
  }

  stepsize->dPrimalStep = dStepSizeUpdate / sqrt(stepsize->dBeta);
  stepsize->dDualStep = dStepSizeUpdate * sqrt(stepsize->dBeta);

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDTEST_Update_Iterate_Adaptive_Step_Size(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  cupdlp_retcode retcode = RETCODE_OK;
  CUPDLPproblem *problem = pdhg->problem;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPtimers *timers = pdhg->timers; // 各种耗时

  cupdlp_float dStepSizeUpdate =
      sqrt(stepsize->dPrimalStep * stepsize->dDualStep);

  cupdlp_bool isDone = false;
  // number of steps this round
  int stepIterThis = 0;
  while (!isDone)
  {
    ++stepsize->nStepSizeIter;
    ++stepIterThis;

    cupdlp_float dPrimalStepUpdate = dStepSizeUpdate / sqrt(stepsize->dBeta);
    cupdlp_float dDualStepUpdate = dStepSizeUpdate * sqrt(stepsize->dBeta);
#pragma region Update
    cupdlp_float dMultiStartTime;
    cupdlp_float dAddStartTime;
    cupdlp_int t_count = *nIter_restart + 1;
    cupdlp_printf("t_count: %d\n", t_count);
    cupdlp_float beta = (t_count + 1) / 2.0;
    cupdlp_float dIterTime = getTimeStamp();

    // 没有必要计算x_md, 因为梯度是常数，用不到x_md
    // x_md^{t} = (1 - 1 / beta^t)x_ag^{t} + (1 / beta^t)x^{t}
    // PDTEST_x_md_step(pdhg, beta);
    // 没有必要进行Project，因为都是已经投影过的x进行线性组合
    // PDHG_Project_Bounds(pdhg, iterates->x_md->data);

    // 计算Ax_bar^{t}, 后面计算y^{t+1}有用
    dMultiStartTime = getTimeStamp();
    Ax(pdhg, iterates->ax_bar, iterates->x_bar);
    timers->dMatVecMultiplyTime += getTimeStamp() - dMultiStartTime;

    // y^{t+1} = y^t + dDualStep * (b - A * (x_bar^{t})
    dAddStartTime = getTimeStamp();
    PDTEST_dualGradientStep(pdhg, stepsize->dDualStep);
    timers->dVecVecAddTime += getTimeStamp() - dAddStartTime;

    PDHG_Project_Row_Duals(pdhg, iterates->yUpdate->data);

    // 计算ATy^{t+1}, 后面计算x^{t+1}有用
    dMultiStartTime = getTimeStamp();
    ATy(pdhg, iterates->atyUpdate, iterates->yUpdate);
    timers->dMatVecMultiplyTime += getTimeStamp() - dMultiStartTime;

    // x^{t+1} = proj_{X}(x^t - dPrimalStep * (c - A'y^{t+1}))
    dAddStartTime = getTimeStamp();
    PDTEST_primalGradientStep(pdhg, stepsize->dPrimalStep);
    timers->dVecVecAddTime += getTimeStamp() - dAddStartTime;
    PDHG_Project_Bounds(pdhg, iterates->xUpdate->data);

    // x_ag^{t+1} = (1 - 1 / beta^t)x_ag^{t} + (1 / beta^t)x^{t}
    dAddStartTime = getTimeStamp();
    PDTEST_x_ag_step(pdhg, beta);
    timers->dVecVecAddTime += getTimeStamp() - dAddStartTime;
    // 没有必要进行Project，因为都是已经投影过的x进行线性组合
    // PDHG_Project_Bounds(pdhg, iterates->x_agUpdate->data);

    // y_ag^{t+1} = (1 - 1 / beta^t)y_ag^{t} + (1 / beta^t)y^{t}
    dAddStartTime = getTimeStamp();
    PDTEST_y_ag_step(pdhg, beta);
    timers->dVecVecAddTime += getTimeStamp() - dAddStartTime;

    // theta^{t+1}
    cupdlp_float theta = t_count / (t_count + 1.0);

    // x_bar^{t+1} = theta^{t+1}(x^{t+1} - x^{t}) + x^{t+1}
    dAddStartTime = getTimeStamp();
    PDTEST_x_bar_step(pdhg, theta);
    timers->dVecVecAddTime += getTimeStamp() - dAddStartTime;

    // 更新一下ax_ag和aty_ag
    dMultiStartTime = getTimeStamp();
    Ax(pdhg, iterates->ax_agUpdate, iterates->x_agUpdate);
    ATy(pdhg, iterates->aty_agUpdate, iterates->y_agUpdate);
    timers->dMatVecMultiplyTime += getTimeStamp() - dMultiStartTime;

    timers->dIterTime += getTimeStamp() - dIterTime;
#pragma endregion

    cupdlp_float dMovement = 0.0;
    cupdlp_float dInteraction = 0.0;

#if !(CUPDLP_CPU) & USE_KERNELS
    PDTEST_compute_interaction_and_movement(pdhg, &dMovement, &dInteraction);
#else
    cupdlp_float dX = 0.0;
    cupdlp_diffTwoNormSquared(pdhg, iterates->x->data, iterates->xUpdate->data,
                              problem->nCols, &dX);
    dX *= 0.5 * sqrt(stepsize->dBeta);

    cupdlp_float dY = 0.0;
    cupdlp_diffTwoNormSquared(pdhg, iterates->y->data, iterates->yUpdate->data,
                              problem->nRows, &dY);
    dY /= 2.0 * sqrt(stepsize->dBeta);
    dMovement = dX + dY;

    //      Δx' (AΔy)
    cupdlp_diffDotDiff(pdhg, iterates->x->data, iterates->xUpdate->data,
                       iterates->aty->data, iterates->atyUpdate->data,
                       problem->nCols, &dInteraction);
#endif

#if CUPDLP_DUMP_LINESEARCH_STATS & CUPDLP_DEBUG
    cupdlp_float dInteractiony = 0.0;
    //      Δy' (AΔx)
    cupdlp_diffDotDiff(pdhg, iterates->y->data, iterates->yUpdate->data,
                       iterates->ax->data, iterates->axUpdate->data,
                       problem->nRows, &dInteractiony);
#endif

    cupdlp_float dStepSizeLimit;
    if (dInteraction != 0.0)
    {
      dStepSizeLimit = dMovement / fabs(dInteraction);
    }
    else
    {
      dStepSizeLimit = INFINITY;
    }
    if (dStepSizeUpdate <= dStepSizeLimit)
    {
      isDone = true;
      // break;
    }
    else
    {
      CUPDLP_CHECK_TIMEOUT(pdhg);
    }

    cupdlp_float dFirstTerm = (1.0 - pow(stepsize->nStepSizeIter + 1.0,
                                         -PDHG_STEPSIZE_REDUCTION_EXP)) *
                              dStepSizeLimit;
    cupdlp_float dSecondTerm =
        (1.0 + pow(stepsize->nStepSizeIter + 1.0, -PDHG_STEPSIZE_GROWTH_EXP)) *
        dStepSizeUpdate;
    dStepSizeUpdate = fmin(dFirstTerm, dSecondTerm);
#if CUPDLP_DUMP_LINESEARCH_STATS & CUPDLP_DEBUG
    cupdlp_printf(" -- stepsize iteration %d: %f %f\n", stepIterThis,
                  dStepSizeUpdate, dStepSizeLimit);

    cupdlp_printf(" -- PrimalStep DualStep: %f %f\n", stepsize->dPrimalStep,
                  stepsize->dDualStep);
    cupdlp_printf(" -- FirstTerm SecondTerm: %f %f\n", dFirstTerm, dSecondTerm);
    cupdlp_printf(" -- nStepSizeIter: %d\n", stepsize->nStepSizeIter);
    cupdlp_printf(" -- RED_EXP GRO_EXP: %f %f\n", PDHG_STEPSIZE_REDUCTION_EXP,
                  PDHG_STEPSIZE_GROWTH_EXP);

    cupdlp_printf("     -- iteraction(x) interaction(y): %f %f\n", dInteraction,
                  dInteractiony);
    cupdlp_printf("     -- movement (scaled norm)  : %f\n", dMovement);
    cupdlp_printf("     -- movement (scaled norm)  : %f\n", dMovement);
    if (stepIterThis > 200)
      break; // avoid unlimited runs due to bugs.
#endif
  }

  stepsize->dPrimalStep = dStepSizeUpdate / sqrt(stepsize->dBeta);
  stepsize->dDualStep = dStepSizeUpdate * sqrt(stepsize->dBeta);

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDHG_Init_Step_Sizes(CUPDLPwork *pdhg)
{
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;

  if (stepsize->eLineSearchMethod == PDHG_FIXED_LINESEARCH)
  {
    CUPDLP_CALL(PDHG_Power_Method(pdhg, &stepsize->dPrimalStep));
    // PDLP Intial primal weight = norm(cost) / norm(rhs) = sqrt(beta)
    // cupdlp_float a = twoNormSquared(problem->cost, problem->nCols);
    // cupdlp_float b = twoNormSquared(problem->rhs, problem->nRows);
    cupdlp_float a = 0.0;
    cupdlp_float b = 0.0;
    cupdlp_twoNormSquared(pdhg, problem->nCols, problem->cost, &a);
    cupdlp_twoNormSquared(pdhg, problem->nRows, problem->rhs, &b);

    if (fmin(a, b) > 1e-6)
    {
      stepsize->dBeta = a / b;
    }
    else
    {
      stepsize->dBeta = 1.0;
    }

    stepsize->dPrimalStep = 0.8 / sqrt(stepsize->dPrimalStep);
    stepsize->dDualStep = stepsize->dPrimalStep;
    stepsize->dPrimalStep /= sqrt(stepsize->dBeta);
    stepsize->dDualStep *= sqrt(stepsize->dBeta);
  }
  else
  {
    stepsize->dTheta = 1.0;

    // PDLP Intial primal weight = norm(cost) / norm(rhs) = sqrt(beta)
    // cupdlp_float a = twoNormSquared(problem->cost, problem->nCols);
    // cupdlp_float b = twoNormSquared(problem->rhs, problem->nRows);
    cupdlp_float a = 0.0;
    cupdlp_float b = 0.0;
    cupdlp_twoNormSquared(pdhg, problem->nCols, problem->cost, &a);
    cupdlp_twoNormSquared(pdhg, problem->nRows, problem->rhs, &b);

    if (fmin(a, b) > 1e-6)
    {
      stepsize->dBeta = a / b;
    }
    else
    {
      stepsize->dBeta = 1.0;
    }
    // infNorm can be avoid by previously calculated infNorm of csc matrix
    stepsize->dPrimalStep =
        // (1.0 / infNorm(problem->data->csc_matrix->colMatElem,
        // problem->data->csc_matrix->nMatElem)) /
        (1.0 / problem->data->csc_matrix->MatElemNormInf) /
        sqrt(stepsize->dBeta);
    stepsize->dDualStep = stepsize->dPrimalStep * stepsize->dBeta;
    iterates->dLastRestartBeta = stepsize->dBeta;
  }

  iterates->iLastRestartIter = 0;
  stepsize->dSumPrimalStep = 0;
  stepsize->dSumDualStep = 0;

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDTEST_Init_Step_Sizes(CUPDLPwork *pdhg)
{
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLPproblem *problem = pdhg->problem;

  PDTESTiterates *iterates = pdhg->PDTESTiterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;

  if (stepsize->eLineSearchMethod == PDHG_FIXED_LINESEARCH)
  {
    CUPDLP_CALL(PDTEST_Power_Method(pdhg, &stepsize->dPrimalStep));
    // PDLP Intial primal weight = norm(cost) / norm(rhs) = sqrt(beta)
    // cupdlp_float a = twoNormSquared(problem->cost, problem->nCols);
    // cupdlp_float b = twoNormSquared(problem->rhs, problem->nRows);

    // a, b记录的是cost和rhs的二范数的平方
    cupdlp_float a = 0.0;
    cupdlp_float b = 0.0;
    cupdlp_twoNormSquared(pdhg, problem->nCols, problem->cost, &a);
    cupdlp_twoNormSquared(pdhg, problem->nRows, problem->rhs, &b);

    if (fmin(a, b) > 1e-6)
    {
      stepsize->dBeta = a / b;
    }
    else
    {
      stepsize->dBeta = 1.0;
    }

    stepsize->dPrimalStep = 0.8 / sqrt(stepsize->dPrimalStep);
    stepsize->dDualStep = stepsize->dPrimalStep;
    // 算出dPrimalStep和dDualStep
    stepsize->dPrimalStep /= sqrt(stepsize->dBeta);
    stepsize->dDualStep *= sqrt(stepsize->dBeta);
  }
  else
  {
    stepsize->dTheta = 1.0;

    // PDLP Intial primal weight = norm(cost) / norm(rhs) = sqrt(beta)
    // cupdlp_float a = twoNormSquared(problem->cost, problem->nCols);
    // cupdlp_float b = twoNormSquared(problem->rhs, problem->nRows);
    cupdlp_float a = 0.0;
    cupdlp_float b = 0.0;
    cupdlp_twoNormSquared(pdhg, problem->nCols, problem->cost, &a);
    cupdlp_twoNormSquared(pdhg, problem->nRows, problem->rhs, &b);

    if (fmin(a, b) > 1e-6)
    {
      stepsize->dBeta = a / b;
    }
    else
    {
      stepsize->dBeta = 1.0;
    }
    // stepsize用无穷范数进行初始化，stepsize->dBeta是论文中的primal weight \omega的平方
    // infNorm can be avoid by previously calculated infNorm of csc matrix
    stepsize->dPrimalStep =
        // (1.0 / infNorm(problem->data->csc_matrix->colMatElem,
        // problem->data->csc_matrix->nMatElem)) /
        (1.0 / problem->data->csc_matrix->MatElemNormInf) /
        sqrt(stepsize->dBeta);
    stepsize->dDualStep = stepsize->dPrimalStep * stepsize->dBeta;
    iterates->dLastRestartBeta = stepsize->dBeta;
  }

  //////////////////////////////////////////////////////
  // CUPDLPdata *data = problem->data;
  // CUPDLP_MATRIX_FORMAT matrix_format = data->matrix_format;
  // cupdlp_printf("matrix_format: %d\n", matrix_format);
  // cupdlp_float matrix_2norm = problem->data->matrix_2norm;
  // stepsize->dPrimalStep = 1 / matrix_2norm;
  // stepsize->dDualStep = 1 / matrix_2norm;
  //////////////////////////////////////////////////////
  iterates->iLastRestartIter = 0;
  stepsize->dSumPrimalStep = 0;
  stepsize->dSumDualStep = 0;

exit_cleanup:
  return retcode;
}

void PDHG_Compute_Average_Iterate(CUPDLPwork *work)
{
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPstepsize *stepsize = work->stepsize;
  CUPDLPiterates *iterates = work->iterates;

  cupdlp_float dPrimalScale =
      stepsize->dSumPrimalStep > 0.0 ? 1.0 / stepsize->dSumPrimalStep : 1.0;
  cupdlp_float dDualScale =
      stepsize->dSumDualStep > 0.0 ? 1.0 / stepsize->dSumDualStep : 1.0;

  // cupdlp_scaleVector(iterates->xAverage, iterates->xSum, dPrimalScale,
  // lp->nCols); cupdlp_scaleVector(iterates->yAverage, iterates->ySum,
  // dDualScale, lp->nRows);

  CUPDLP_COPY_VEC(iterates->xAverage->data, iterates->xSum, cupdlp_float,
                  lp->nCols);
  CUPDLP_COPY_VEC(iterates->yAverage->data, iterates->ySum, cupdlp_float,
                  lp->nRows);
  cupdlp_scaleVector(work, dPrimalScale, iterates->xAverage->data, lp->nCols);
  cupdlp_scaleVector(work, dDualScale, iterates->yAverage->data, lp->nRows);

  // Ax(work, iterates->axAverage, iterates->xAverage);
  // ATyCPU(work, iterates->atyAverage, iterates->yAverage);
  Ax(work, iterates->axAverage, iterates->xAverage);
  ATy(work, iterates->atyAverage, iterates->yAverage);
}

void PDTEST_Compute_Average_Iterate(CUPDLPwork *work)
{
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPstepsize *stepsize = work->stepsize;
  PDTESTiterates *iterates = work->PDTESTiterates;

  cupdlp_float dPrimalScale =
      stepsize->dSumPrimalStep > 0.0 ? 1.0 / stepsize->dSumPrimalStep : 1.0;
  cupdlp_float dDualScale =
      stepsize->dSumDualStep > 0.0 ? 1.0 / stepsize->dSumDualStep : 1.0;

  // 把sum copy 给 average
  CUPDLP_COPY_VEC(iterates->x_agAverage->data, iterates->x_agSum, cupdlp_float,
                  lp->nCols);
  CUPDLP_COPY_VEC(iterates->y_agAverage->data, iterates->y_agSum, cupdlp_float,
                  lp->nRows);
  // 把 average 乘上 1/sumstep
  cupdlp_scaleVector(work, dPrimalScale, iterates->x_agAverage->data, lp->nCols);
  cupdlp_scaleVector(work, dDualScale, iterates->y_agAverage->data, lp->nRows);
  // cupdlp_printf("dPrimalScale: %f\n", dPrimalScale);
  // cupdlp_printf("x_agAverage:\n");
  // PDTEST_printCudaDenseVecGPU(iterates->x_agAverage);
  // Ax(work, iterates->axAverage, iterates->xAverage);
  // ATyCPU(work, iterates->atyAverage, iterates->yAverage);
  Ax(work, iterates->ax_agAverage, iterates->x_agAverage);
  ATy(work, iterates->aty_agAverage, iterates->y_agAverage);
}

void PDHG_Update_Average(CUPDLPwork *work)
{
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPstepsize *stepsize = work->stepsize;
  CUPDLPiterates *iterates = work->iterates;

  // PDLP weighs average iterates in this way
  cupdlp_float dMeanStepSize =
      sqrt(stepsize->dPrimalStep * stepsize->dDualStep);
  // AddToVector(iterates->xSum, dMeanStepSize, iterates->xUpdate,
  // lp->nCols); AddToVector(iterates->ySum, dMeanStepSize,
  // iterates->yUpdate, lp->nRows);
  cupdlp_axpy(work, lp->nCols, &dMeanStepSize, iterates->xUpdate->data,
              iterates->xSum);
  cupdlp_axpy(work, lp->nRows, &dMeanStepSize, iterates->yUpdate->data,
              iterates->ySum);

  stepsize->dSumPrimalStep += dMeanStepSize;
  stepsize->dSumDualStep += dMeanStepSize;
}

void PDTEST_Update_Average(CUPDLPwork *work)
{
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPstepsize *stepsize = work->stepsize;
  PDTESTiterates *iterates = work->PDTESTiterates;

  // PDLP weighs average iterates in this way
  cupdlp_float dMeanStepSize =
      sqrt(stepsize->dPrimalStep * stepsize->dDualStep);

  cupdlp_axpy(work, lp->nCols, &dMeanStepSize, iterates->x_agUpdate->data,
              iterates->x_agSum);
  cupdlp_axpy(work, lp->nRows, &dMeanStepSize, iterates->y_agUpdate->data,
              iterates->y_agSum);

  stepsize->dSumPrimalStep += dMeanStepSize;
  stepsize->dSumDualStep += dMeanStepSize;
}

cupdlp_retcode PDHG_Update_Iterate(CUPDLPwork *pdhg)
{
  cupdlp_retcode retcode = RETCODE_OK;

#if PDHG_USE_TIMERS
  CUPDLPtimers *timers = pdhg->timers;
  ++timers->nUpdateIterateCalls;
  cupdlp_float dStartTime = getTimeStamp();
#endif

  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPiterates *iterates = pdhg->iterates;

  switch (stepsize->eLineSearchMethod)
  {
  case PDHG_FIXED_LINESEARCH:
    PDHG_Update_Iterate_Constant_Step_Size(pdhg);
    break;
  case PDHG_MALITSKY_POCK_LINESEARCH:
    PDHG_Update_Iterate_Malitsky_Pock(pdhg);
    break;
  case PDHG_ADAPTIVE_LINESEARCH:
    CUPDLP_CALL(PDHG_Update_Iterate_Adaptive_Step_Size(pdhg));
    break;
  }

  PDHG_Update_Average(pdhg);

  CUPDLP_COPY_VEC(iterates->x->data, iterates->xUpdate->data, cupdlp_float,
                  problem->nCols);
  CUPDLP_COPY_VEC(iterates->y->data, iterates->yUpdate->data, cupdlp_float,
                  problem->nRows);
  CUPDLP_COPY_VEC(iterates->ax->data, iterates->axUpdate->data, cupdlp_float,
                  problem->nRows);
  CUPDLP_COPY_VEC(iterates->aty->data, iterates->atyUpdate->data, cupdlp_float,
                  problem->nCols);

#if PDHG_USE_TIMERS
  timers->dUpdateIterateTime += getTimeStamp() - dStartTime;
#endif

exit_cleanup:
  return RETCODE_OK;
}

cupdlp_retcode PDTEST_Average_Update_Iterate(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  cupdlp_retcode retcode = RETCODE_OK;

#if PDHG_USE_TIMERS
  CUPDLPtimers *timers = pdhg->timers;
  ++timers->nUpdateIterateCalls;
  cupdlp_float dStartTime = getTimeStamp();
#endif

  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;

  // PDTEST_Update_Iterate_Constant_Step_Size(pdhg, nIterL_restart);
  // 选择各种更新，看看效果
  switch (stepsize->eLineSearchMethod)
  {
  case PDHG_FIXED_LINESEARCH:
    PDTEST_Update_Iterate_Constant_Step_Size(pdhg, nIter_restart);
    break;
  case PDHG_ADAPTIVE_LINESEARCH:
    cupdlp_printf("PDTEST_ADAPTIVE_LINESEARCH\n");
    CUPDLP_CALL(PDTEST_Update_Iterate_Adaptive_Step_Size_ag(pdhg, nIter_restart));
    // CUPDLP_CALL(PDTEST_Update_Iterate_Adaptive_Step_Size(pdhg, nIter_restart));
    break;
  }

  // 其实是更新x_agSum和y_agSum
  PDTEST_Update_Average(pdhg);

  // 使用CUPDLP_COPY_VEC宏复制更新后的原始变量x和对偶变量y，以及它们对应的辅助变量ax和aty到更新变量xUpdate、yUpdate、axUpdate和atyUpdate。这一步确保了算法中使用的变量是最新的迭代结果
  CUPDLP_COPY_VEC(iterates->x->data, iterates->xUpdate->data, cupdlp_float, problem->nCols);
  CUPDLP_COPY_VEC(iterates->x_bar->data, iterates->x_barUpdate->data, cupdlp_float, problem->nCols);
  CUPDLP_COPY_VEC(iterates->x_ag->data, iterates->x_agUpdate->data, cupdlp_float, problem->nCols);
  CUPDLP_COPY_VEC(iterates->y->data, iterates->yUpdate->data, cupdlp_float, problem->nRows);
  CUPDLP_COPY_VEC(iterates->y_ag->data, iterates->y_agUpdate->data, cupdlp_float, problem->nRows);

  CUPDLP_COPY_VEC(iterates->ax_ag->data, iterates->ax_agUpdate->data, cupdlp_float,
                  problem->nRows);
  CUPDLP_COPY_VEC(iterates->aty_ag->data, iterates->aty_agUpdate->data, cupdlp_float,
                  problem->nCols);

#if PDHG_USE_TIMERS
  timers->dUpdateIterateTime += getTimeStamp() - dStartTime;
#endif

exit_cleanup:
  return RETCODE_OK;
}
cupdlp_retcode PDTEST_Update_Iterate(CUPDLPwork *pdhg, cupdlp_int *nIter_restart)
{
  cupdlp_retcode retcode = RETCODE_OK;

#if PDHG_USE_TIMERS
  CUPDLPtimers *timers = pdhg->timers;
  ++timers->nUpdateIterateCalls;
  cupdlp_float dStartTime = getTimeStamp();
#endif

  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;

  // PDTEST_Update_Iterate_Constant_Step_Size(pdhg, nIterL_restart);
  // 选择各种更新，看看效果
  switch (stepsize->eLineSearchMethod)
  {
  case PDHG_FIXED_LINESEARCH:
    PDTEST_Update_Iterate_Constant_Step_Size(pdhg, nIter_restart);
    break;
  case PDHG_ADAPTIVE_LINESEARCH:
    cupdlp_printf("PDTEST_ADAPTIVE_LINESEARCH\n");
    CUPDLP_CALL(PDTEST_Update_Iterate_Adaptive_Step_Size_ag(pdhg, nIter_restart));
    // CUPDLP_CALL(PDTEST_Update_Iterate_Adaptive_Step_Size(pdhg, nIter_restart));
    break;
  }

  PDTEST_Update_Average(pdhg);

  // 使用CUPDLP_COPY_VEC宏复制更新后的原始变量x和对偶变量y，以及它们对应的辅助变量ax和aty到更新变量xUpdate、yUpdate、axUpdate和atyUpdate。这一步确保了算法中使用的变量是最新的迭代结果
  CUPDLP_COPY_VEC(iterates->x->data, iterates->xUpdate->data, cupdlp_float, problem->nCols);
  CUPDLP_COPY_VEC(iterates->x_bar->data, iterates->x_barUpdate->data, cupdlp_float, problem->nCols);
  CUPDLP_COPY_VEC(iterates->x_ag->data, iterates->x_agUpdate->data, cupdlp_float, problem->nCols);
  CUPDLP_COPY_VEC(iterates->y->data, iterates->yUpdate->data, cupdlp_float, problem->nRows);
  CUPDLP_COPY_VEC(iterates->y_ag->data, iterates->y_agUpdate->data, cupdlp_float, problem->nRows);

  CUPDLP_COPY_VEC(iterates->ax_ag->data, iterates->ax_agUpdate->data, cupdlp_float,
                  problem->nRows);
  CUPDLP_COPY_VEC(iterates->aty_ag->data, iterates->aty_agUpdate->data, cupdlp_float,
                  problem->nCols);

#if PDHG_USE_TIMERS
  timers->dUpdateIterateTime += getTimeStamp() - dStartTime;
#endif

exit_cleanup:
  return RETCODE_OK;
}