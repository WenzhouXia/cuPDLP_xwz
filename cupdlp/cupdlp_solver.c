
#include "cupdlp_solver.h"

#include "cupdlp_defs.h"
#include "cupdlp_linalg.h"
#include "cupdlp_proj.h"
#include "cupdlp_restart.h"
// #include "cupdlp_scaling.h"
// #include "cupdlp_scaling_new.h"
#include "cupdlp_step.h"
#include "cupdlp_utils.h"
#include "cupdlp_multiscale.h"
#include "glbopts.h"

void PDHG_Compute_Primal_Feasibility(CUPDLPwork *work, double *primalResidual,
                                     const double *ax, const double *x,
                                     double *dPrimalFeasibility,
                                     double *dPrimalObj)
{
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPscaling *scaling = work->scaling;

  // primal variable violation

  // todo, add this
  //    *dPrimalObj = Dotprod_Neumaier(problem->cost, x, lp->nCols);
  cupdlp_dot(work, lp->nCols, x, problem->cost, dPrimalObj);
  *dPrimalObj = *dPrimalObj * problem->sign_origin - problem->offset;

  // cupdlp_copy(primalResidual, ax, cupdlp_float, lp->nRows);
  CUPDLP_COPY_VEC(primalResidual, ax, cupdlp_float, lp->nRows);

  // AddToVector(primalResidual, -1.0, problem->rhs, lp->nRows);
  cupdlp_float alpha = -1.0;
  cupdlp_axpy(work, lp->nRows, &alpha, problem->rhs, primalResidual);

  double dPrimalFeas = 0.0;

  // todo, check
  //  cupdlp_projNegative(primalResidual + problem->nEqs, primalResidual +
  //  problem->nEqs, lp->nRows - problem->nEqs);
  //

  cupdlp_projNeg(primalResidual + problem->nEqs, lp->nRows - problem->nEqs);

  if (scaling->ifScaled)
  {
    // cupdlp_edot(primalResidual, scaling->rowScale, lp->nRows);
    // cupdlp_edot(primalResidual, scaling->rowScale_gpu, lp->nRows);

    // cupdlp_copy_vec(work->buffer3, scaling->rowScale, cupdlp_float,
    // lp->nRows); cupdlp_edot(primalResidual, work->buffer3, lp->nRows);

    cupdlp_edot(primalResidual, work->rowScale, lp->nRows);
  }

  cupdlp_twoNorm(work, lp->nRows, primalResidual, dPrimalFeasibility);
}

void PDHG_Compute_Dual_Feasibility(CUPDLPwork *work, double *dualResidual,
                                   const double *aty, const double *x,
                                   const double *y, double *dDualFeasibility,
                                   double *dDualObj, double *dComplementarity)
{
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPresobj *resobj = work->resobj;
  CUPDLPscaling *scaling = work->scaling;
  // todo, compute Neumaier
  //    *dDualObj = Dotprod_Neumaier(problem->rhs, y, lp->nRows);
  cupdlp_dot(work, lp->nRows, y, problem->rhs, dDualObj);
  *dDualObj = *dDualObj * problem->sign_origin - problem->offset;
  *dComplementarity = 0.0;
  // @note:
  // original dual residual in pdlp:
  // they compute:
  //    violation   +  reduced cost
  //  |max(-y, 0)|  + |(I-Π)(c-Α'υ)|
  // compute c - A'y

  CUPDLP_COPY_VEC(dualResidual, aty, cupdlp_float, lp->nCols);
  cupdlp_float alpha = -1.0;
  cupdlp_scaleVector(work, alpha, dualResidual, lp->nCols);

  alpha = 1.0;
  cupdlp_axpy(work, lp->nCols, &alpha, problem->cost, dualResidual);

  //    julia version
  //        function compute_reduced_costs_from_primal_gradient_kernel!(
  //            primal_gradient::CuDeviceVector{Float64},
  //            isfinite_variable_lower_bound::CuDeviceVector{Bool},
  //            isfinite_variable_upper_bound::CuDeviceVector{Bool},
  //            num_variables::Int64,
  //            reduced_costs::CuDeviceVector{Float64},
  //            reduced_costs_violation::CuDeviceVector{Float64},
  //        )
  //            tx = threadIdx().x + (blockDim().x * (blockIdx().x - 0x1))
  //            if tx <= num_variables
  //                @inbounds begin
  //                    reduced_costs[tx] = max(primal_gradient[tx], 0.0) *
  //                    isfinite_variable_lower_bound[tx] +
  //                    min(primal_gradient[tx], 0.0) *
  //                    isfinite_variable_upper_bound[tx]
  //
  //                    reduced_costs_violation[tx] = primal_gradient[tx] -
  //                    reduced_costs[tx]
  //                end
  //            end
  //            return
  //        end

  // cupdlp_copy(resobj->dSlackPos, dualResidual, cupdlp_float, lp->nCols);
  CUPDLP_COPY_VEC(resobj->dSlackPos, dualResidual, cupdlp_float, lp->nCols);

  // cupdlp_projPositive(resobj->dSlackPos, resobj->dSlackPos, lp->nCols);
  cupdlp_projPos(resobj->dSlackPos, lp->nCols);

  // cupdlp_cdot_fb(resobj->dSlackPos, problem->hasLower, lp->nCols);
  cupdlp_edot(resobj->dSlackPos, problem->hasLower, lp->nCols);

  cupdlp_float temp = 0.0;
  cupdlp_dot(work, lp->nCols, x, resobj->dSlackPos, &temp);
  *dComplementarity += temp;
  cupdlp_dot(work, lp->nCols, resobj->dSlackPos, resobj->dLowerFiltered, &temp);
  *dComplementarity -= temp;
  cupdlp_dot(work, lp->nCols, resobj->dSlackPos, resobj->dLowerFiltered, &temp);
  *dDualObj += temp;

  CUPDLP_COPY_VEC(resobj->dSlackNeg, dualResidual, cupdlp_float, lp->nCols);

  cupdlp_projNeg(resobj->dSlackNeg, lp->nCols);

  // ScaleVector(-1.0, resobj->dSlackNeg, lp->nCols);
  cupdlp_scaleVector(work, -1.0, resobj->dSlackNeg, lp->nCols);

  // cupdlp_cdot_fb(resobj->dSlackNeg, problem->hasUpper, lp->nCols);
  cupdlp_edot(resobj->dSlackNeg, problem->hasUpper, lp->nCols);

  cupdlp_dot(work, lp->nCols, x, resobj->dSlackNeg, &temp);
  *dComplementarity -= temp;
  cupdlp_dot(work, lp->nCols, resobj->dSlackNeg, resobj->dUpperFiltered, &temp);
  *dComplementarity += temp;
  cupdlp_dot(work, lp->nCols, resobj->dSlackNeg, resobj->dUpperFiltered, &temp);
  *dDualObj -= temp;

  alpha = -1.0;
  cupdlp_axpy(work, lp->nCols, &alpha, resobj->dSlackPos, dualResidual);
  alpha = 1.0;
  cupdlp_axpy(work, lp->nCols, &alpha, resobj->dSlackNeg, dualResidual);

  if (scaling->ifScaled)
  {
    // cupdlp_edot(dualResidual, scaling->colScale, lp->nCols);
    // cupdlp_edot(dualResidual, scaling->colScale_gpu, lp->nCols);

    // cupdlp_copy_vec(work->buffer3, scaling->colScale, cupdlp_float,
    // lp->nCols); cupdlp_edot(dualResidual, work->buffer3, lp->nCols);

    cupdlp_edot(dualResidual, work->colScale, lp->nCols);
  }

  cupdlp_twoNorm(work, lp->nCols, dualResidual, dDualFeasibility);
}

void PDTEST_printCudaDenseVecGPU(const CUPDLPvec *vec)
{
  // 分配CPU内存
  cupdlp_float *hostData = (cupdlp_float *)malloc(vec->len * sizeof(cupdlp_float));
  // 从GPU到CPU复制数据
  cudaMemcpy(hostData, vec->data, vec->len * sizeof(cupdlp_float), cudaMemcpyDeviceToHost);
  cupdlp_printf("Vector length: %d\n", vec->len);
  cupdlp_printf("Vector elements:\n");
  for (int i = 0; i < vec->len; ++i)
  {
    cupdlp_printf("Element[%d]: %f\n", i, hostData[i]);
  }
  // 释放CPU内存
  free(hostData);
}
void PDTEST_printCudaMatGPU(CUPDLPwork *work)
{
  CUPDLPproblem *problem = work->problem;
  CUPDLPcsc *csc_matrix = problem->data->csc_matrix;
  cupdlp_int nMatElem = csc_matrix->nMatElem;
  cupdlp_float *hostData = (cupdlp_float *)malloc(problem->data->csc_matrix->nMatElem * sizeof(cupdlp_float));
  cupdlp_float *matelem = problem->data->csc_matrix->colMatElem;
  cudaMemcpy(hostData, matelem, problem->data->csc_matrix->nMatElem * sizeof(cupdlp_float), cudaMemcpyDeviceToHost);
  for (int i = 0; i < problem->data->csc_matrix->nMatElem; ++i)
  {
    cupdlp_printf("MatElem[%d]: %f\n", i, hostData[i]);
  }
}

void PDTEST_printCudafloatGPU(cupdlp_float *data, int size)
{
  // 在主机内存中分配空间来存储成本向量
  cupdlp_float *hostData = (cupdlp_float *)malloc(size * sizeof(cupdlp_float));

  // 从GPU复制到主机
  cudaMemcpy(hostData, data, size * sizeof(cupdlp_float), cudaMemcpyDeviceToHost);

  // 打印成本向量
  cupdlp_printf("Vector length: %d\n", size);
  cupdlp_printf("Vector elements:\n");
  for (int i = 0; i < size; ++i)
  {
    cupdlp_printf("Element[%d]: %f\n", i, hostData[i]);
  }

  // 释放主机内存
  free(hostData);
}
void PDTEST_Compute_dDualObj(CUPDLPwork *work)
{
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  PDTESTiterates *iterates = work->PDTESTiterates;
  CUPDLPresobj *resobj = work->resobj;

  // cupdlp_float *hostData = (cupdlp_float *)malloc(y_ag->len * sizeof(cupdlp_float));
  // // 从GPU复制到CPU
  // cudaMemcpy(hostData, y_ag->data, y_ag->len * sizeof(cupdlp_float), cudaMemcpyDeviceToHost);
  // cupdlp_printf("Vector length: %d\n", y_ag->len);
  // cupdlp_printf("Vector elements:\n");
  // for (cupdlp_int i = 0; i < y_ag->len; ++i)
  // {
  //   cupdlp_printf("Element[%d]: %f\n", i, hostData[i]);
  // }
  // *dDualObj = Dotprod_Neumaier(problem->rhs, y, lp->nRows);
  cupdlp_dot(work, lp->nRows, iterates->y_ag->data, problem->rhs, &resobj->dDualObj);
  cupdlp_printf("dDualObj: %f\n", resobj->dDualObj);
}

void PDHG_Compute_Residuals(CUPDLPwork *work)
{
#if problem_USE_TIMERS
  ++problem->nComputeResidualsCalls;
  double dStartTime = getTimeStamp();
#endif
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPresobj *resobj = work->resobj;
  CUPDLPiterates *iterates = work->iterates;
  CUPDLPscaling *scaling = work->scaling;
  CUPDLPsettings *settings = work->settings;

  PDHG_Compute_Primal_Feasibility(work, resobj->primalResidual,
                                  iterates->ax->data, iterates->x->data,
                                  &resobj->dPrimalFeas, &resobj->dPrimalObj);
  PDHG_Compute_Dual_Feasibility(work, resobj->dualResidual, iterates->aty->data,
                                iterates->x->data, iterates->y->data,
                                &resobj->dDualFeas, &resobj->dDualObj,
                                &resobj->dComplementarity);

  PDHG_Compute_Primal_Feasibility(
      work, resobj->primalResidualAverage, iterates->axAverage->data,
      iterates->xAverage->data, &resobj->dPrimalFeasAverage,
      &resobj->dPrimalObjAverage);
  PDHG_Compute_Dual_Feasibility(
      work, resobj->dualResidualAverage, iterates->atyAverage->data,
      iterates->xAverage->data, iterates->yAverage->data,
      &resobj->dDualFeasAverage, &resobj->dDualObjAverage,
      &resobj->dComplementarityAverage);

  // resobj->dPrimalObj /= (scaling->dObjScale * scaling->dObjScale);
  // resobj->dDualObj /= (scaling->dObjScale * scaling->dObjScale);
  resobj->dDualityGap = resobj->dPrimalObj - resobj->dDualObj;
  resobj->dRelObjGap =
      fabs(resobj->dPrimalObj - resobj->dDualObj) /
      (1.0 + fabs(resobj->dPrimalObj) + fabs(resobj->dDualObj));

  // resobj->dPrimalObjAverage /= scaling->dObjScale * scaling->dObjScale;
  // resobj->dDualObjAverage /= scaling->dObjScale * scaling->dObjScale;
  resobj->dDualityGapAverage =
      resobj->dPrimalObjAverage - resobj->dDualObjAverage;
  resobj->dRelObjGapAverage =
      fabs(resobj->dPrimalObjAverage - resobj->dDualObjAverage) /
      (1.0 + fabs(resobj->dPrimalObjAverage) + fabs(resobj->dDualObjAverage));

#if problem_USE_TIMERS
  problem->dComputeResidualsTime += getTimeStamp() - dStartTime;
#endif
}

void PDTEST_Average_Compute_Residuals(CUPDLPwork *work)
{
#if problem_USE_TIMERS
  ++problem->nComputeResidualsCalls;
  double dStartTime = getTimeStamp();
#endif
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPresobj *resobj = work->resobj;
  PDTESTiterates *iterates = work->PDTESTiterates;
  CUPDLPscaling *scaling = work->scaling;
  CUPDLPsettings *settings = work->settings;

  PDHG_Compute_Primal_Feasibility(work, resobj->primalResidual,
                                  iterates->ax_ag->data, iterates->x_ag->data,
                                  &resobj->dPrimalFeas, &resobj->dPrimalObj);
  PDHG_Compute_Dual_Feasibility(work, resobj->dualResidual, iterates->aty_ag->data, iterates->x_ag->data, iterates->y_ag->data, &resobj->dDualFeas, &resobj->dDualObj, &resobj->dComplementarity);

  PDHG_Compute_Primal_Feasibility(
      work, resobj->primalResidualAverage, iterates->ax_agAverage->data,
      iterates->x_agAverage->data, &resobj->dPrimalFeasAverage,
      &resobj->dPrimalObjAverage);
  PDHG_Compute_Dual_Feasibility(
      work, resobj->dualResidualAverage, iterates->aty_agAverage->data,
      iterates->x_agAverage->data, iterates->y_agAverage->data,
      &resobj->dDualFeasAverage, &resobj->dDualObjAverage,
      &resobj->dComplementarityAverage);

  resobj->dDualityGap = resobj->dPrimalObj - resobj->dDualObj;
  resobj->dRelObjGap =
      fabs(resobj->dPrimalObj - resobj->dDualObj) /
      (1.0 + fabs(resobj->dPrimalObj) + fabs(resobj->dDualObj));

  resobj->dDualityGapAverage =
      resobj->dPrimalObjAverage - resobj->dDualObjAverage;
  resobj->dRelObjGapAverage =
      fabs(resobj->dPrimalObjAverage - resobj->dDualObjAverage) /
      (1.0 + fabs(resobj->dPrimalObjAverage) + fabs(resobj->dDualObjAverage));

#if problem_USE_TIMERS
  problem->dComputeResidualsTime += getTimeStamp() - dStartTime;
#endif
}

void PDTEST_Compute_Residuals_best(CUPDLPwork *work)
{
#if problem_USE_TIMERS
  ++problem->nComputeResidualsCalls;
  double dStartTime = getTimeStamp();
#endif
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPresobj *resobj = work->resobj;
  PDTESTiterates *iterates = work->PDTESTiterates;
  CUPDLPscaling *scaling = work->scaling;
  CUPDLPsettings *settings = work->settings;
  cupdlp_int bestID = settings->bestID;

  switch (bestID)
  {
  case 0:
    PDTEST_Compute_Residuals(work);
    break;
  case 1:
    PDTEST_Compute_Residuals(work);
    break;
  case 2:
    PDTEST_Compute_Residuals(work);
    break;
  case 3:
    PDTEST_Compute_Residuals(work);
    break;
  case 4:
    PDTEST_Compute_Residuals(work);
    break;
  case 5:
    PDTEST_Compute_Residuals(work);
    break;
  case 6:
    PDTEST_Compute_Residuals_CurrentandAverage(work);
    break;
  default:
    cupdlp_printf("bestID not supported\n");
    break;
  }

#if problem_USE_TIMERS
  problem->dComputeResidualsTime += getTimeStamp() - dStartTime;
#endif
}

void PDTEST_Compute_Residuals(CUPDLPwork *work)
{
#if problem_USE_TIMERS
  ++problem->nComputeResidualsCalls;
  double dStartTime = getTimeStamp();
#endif
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPresobj *resobj = work->resobj;
  PDTESTiterates *iterates = work->PDTESTiterates;
  CUPDLPscaling *scaling = work->scaling;
  CUPDLPsettings *settings = work->settings;

  PDHG_Compute_Primal_Feasibility(work, resobj->primalResidual,
                                  iterates->ax_ag->data, iterates->x_ag->data,
                                  &resobj->dPrimalFeas, &resobj->dPrimalObj);
  PDHG_Compute_Dual_Feasibility(work, resobj->dualResidual, iterates->aty_ag->data, iterates->x_ag->data, iterates->y_ag->data, &resobj->dDualFeas, &resobj->dDualObj, &resobj->dComplementarity);

  // PDHG_Compute_Primal_Feasibility(
  //     work, resobj->primalResidualAverage, iterates->axAverage->data,
  //     iterates->xAverage->data, &resobj->dPrimalFeasAverage,
  //     &resobj->dPrimalObjAverage);
  // PDHG_Compute_Dual_Feasibility(
  //     work, resobj->dualResidualAverage, iterates->atyAverage->data,
  //     iterates->xAverage->data, iterates->yAverage->data,
  //     &resobj->dDualFeasAverage, &resobj->dDualObjAverage,
  //     &resobj->dComplementarityAverage);

  resobj->dDualityGap = resobj->dPrimalObj - resobj->dDualObj;
  resobj->dRelObjGap =
      fabs(resobj->dPrimalObj - resobj->dDualObj) /
      (1.0 + fabs(resobj->dPrimalObj) + fabs(resobj->dDualObj));

  // resobj->dDualityGapAverage =
  //     resobj->dPrimalObjAverage - resobj->dDualObjAverage;
  // resobj->dRelObjGapAverage =
  //     fabs(resobj->dPrimalObjAverage - resobj->dDualObjAverage) /
  //     (1.0 + fabs(resobj->dPrimalObjAverage) + fabs(resobj->dDualObjAverage));

#if problem_USE_TIMERS
  problem->dComputeResidualsTime += getTimeStamp() - dStartTime;
#endif
}

void PDTEST_Compute_Residuals_CurrentandAverage(CUPDLPwork *work)
{
#if problem_USE_TIMERS
  ++problem->nComputeResidualsCalls;
  double dStartTime = getTimeStamp();
#endif
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPresobj *resobj = work->resobj;
  PDTESTiterates *iterates = work->PDTESTiterates;
  CUPDLPscaling *scaling = work->scaling;
  CUPDLPsettings *settings = work->settings;

  PDHG_Compute_Primal_Feasibility(work, resobj->primalResidual,
                                  iterates->ax->data, iterates->x->data,
                                  &resobj->dPrimalFeas, &resobj->dPrimalObj);
  PDHG_Compute_Dual_Feasibility(work, resobj->dualResidual, iterates->aty->data, iterates->x->data, iterates->y->data, &resobj->dDualFeas, &resobj->dDualObj, &resobj->dComplementarity);

  PDHG_Compute_Primal_Feasibility(work, resobj->primalResidualAverage, iterates->ax_ag->data, iterates->x_ag->data, &resobj->dPrimalFeasAverage, &resobj->dPrimalObjAverage);
  PDHG_Compute_Dual_Feasibility(work, resobj->dualResidualAverage, iterates->aty_ag->data, iterates->x_ag->data, iterates->y_ag->data, &resobj->dDualFeasAverage, &resobj->dDualObjAverage, &resobj->dComplementarityAverage);

  resobj->dDualityGap = resobj->dPrimalObj - resobj->dDualObj;
  resobj->dRelObjGap =
      fabs(resobj->dPrimalObj - resobj->dDualObj) /
      (1.0 + fabs(resobj->dPrimalObj) + fabs(resobj->dDualObj));

  resobj->dDualityGapAverage = resobj->dPrimalObjAverage - resobj->dDualObjAverage;
  resobj->dRelObjGapAverage =
      fabs(resobj->dPrimalObjAverage - resobj->dDualObjAverage) /
      (1.0 + fabs(resobj->dPrimalObjAverage) + fabs(resobj->dDualObjAverage));

#if problem_USE_TIMERS
  problem->dComputeResidualsTime += getTimeStamp() - dStartTime;
#endif
}

void PDHG_Init_Variables_Multiscale(CUPDLPwork *work, cupdlp_float *x_init, cupdlp_float *y_init)
{
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPstepsize *stepsize = work->stepsize;
  CUPDLPiterates *iterates = work->iterates;

  // CUPDLP_ZERO_VEC(iterates->x->data, cupdlp_float, lp->nCols);
  CUPDLP_COPY_VEC(iterates->x->data, x_init, cupdlp_float, problem->nCols);
  // XXX: PDLP Does not project x0,  so we uncomment for 1-1 comparison

  PDHG_Project_Bounds(work, iterates->x->data);

  // CUPDLP_ZERO_VEC(iterates->y->data, cupdlp_float, lp->nRows);
  CUPDLP_COPY_VEC(iterates->y->data, y_init, cupdlp_float, problem->nRows);

  Ax(work, iterates->ax, iterates->x);
  ATy(work, iterates->aty, iterates->y);

  // CUPDLP_ZERO_VEC(iterates->xSum, cupdlp_float, lp->nCols);
  // CUPDLP_ZERO_VEC(iterates->ySum, cupdlp_float, lp->nRows);
  // CUPDLP_ZERO_VEC(iterates->xAverage->data, cupdlp_float, lp->nCols);
  // CUPDLP_ZERO_VEC(iterates->yAverage->data, cupdlp_float, lp->nRows);
  CUPDLP_COPY_VEC(iterates->xSum, iterates->x->data, cupdlp_float, problem->nCols);
  CUPDLP_COPY_VEC(iterates->ySum, iterates->y->data, cupdlp_float, problem->nRows);
  CUPDLP_COPY_VEC(iterates->xAverage->data, iterates->x->data, cupdlp_float, problem->nCols);
  CUPDLP_COPY_VEC(iterates->yAverage->data, iterates->y->data, cupdlp_float, problem->nRows);

  PDHG_Project_Bounds(work, iterates->xSum);
  PDHG_Project_Bounds(work, iterates->xAverage->data);

  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;

  CUPDLP_ZERO_VEC(iterates->xLastRestart, cupdlp_float, lp->nCols);
  CUPDLP_ZERO_VEC(iterates->yLastRestart, cupdlp_float, lp->nRows);
}

void PDHG_Init_Variables(CUPDLPwork *work)
{
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPstepsize *stepsize = work->stepsize;
  CUPDLPiterates *iterates = work->iterates;

  // cupdlp_zero(iterates->x, cupdlp_float, lp->nCols);
  CUPDLP_ZERO_VEC(iterates->x->data, cupdlp_float, lp->nCols);

  // XXX: PDLP Does not project x0,  so we uncomment for 1-1 comparison

  PDHG_Project_Bounds(work, iterates->x->data);

  // cupdlp_zero(iterates->y, cupdlp_float, lp->nRows);
  CUPDLP_ZERO_VEC(iterates->y->data, cupdlp_float, lp->nRows);

  // Ax(work, iterates->ax, iterates->x);
  // ATyCPU(work, iterates->aty, iterates->y);
  Ax(work, iterates->ax, iterates->x);
  ATy(work, iterates->aty, iterates->y);

  // cupdlp_zero(iterates->xSum, cupdlp_float, lp->nCols);
  // cupdlp_zero(iterates->ySum, cupdlp_float, lp->nRows);
  // cupdlp_zero(iterates->xAverage, cupdlp_float, lp->nCols);
  // cupdlp_zero(iterates->yAverage, cupdlp_float, lp->nRows);
  CUPDLP_ZERO_VEC(iterates->xSum, cupdlp_float, lp->nCols);
  CUPDLP_ZERO_VEC(iterates->ySum, cupdlp_float, lp->nRows);
  CUPDLP_ZERO_VEC(iterates->xAverage->data, cupdlp_float, lp->nCols);
  CUPDLP_ZERO_VEC(iterates->yAverage->data, cupdlp_float, lp->nRows);

  PDHG_Project_Bounds(work, iterates->xSum);
  PDHG_Project_Bounds(work, iterates->xAverage->data);

  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;

  CUPDLP_ZERO_VEC(iterates->xLastRestart, cupdlp_float, lp->nCols);
  CUPDLP_ZERO_VEC(iterates->yLastRestart, cupdlp_float, lp->nRows);
}

void PDTEST_Average_Init_Variables(CUPDLPwork *work)
{
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPstepsize *stepsize = work->stepsize;
  PDTESTiterates *iterates = work->PDTESTiterates;

  CUPDLP_ZERO_VEC(iterates->x->data, cupdlp_float, lp->nCols);
  PDHG_Project_Bounds(work, iterates->x->data);

  CUPDLP_ZERO_VEC(iterates->x_ag->data, cupdlp_float, lp->nCols);
  PDHG_Project_Bounds(work, iterates->x_ag->data);

  CUPDLP_ZERO_VEC(iterates->x_bar->data, cupdlp_float, lp->nCols);
  PDHG_Project_Bounds(work, iterates->x_bar->data);

  // CUPDLP_ZERO_VEC(iterates->x_md->data, cupdlp_float, lp->nCols);
  // PDHG_Project_Bounds(work, iterates->x_md->data);

  CUPDLP_ZERO_VEC(iterates->y->data, cupdlp_float, lp->nRows);
  CUPDLP_ZERO_VEC(iterates->y_ag->data, cupdlp_float, lp->nRows);

  Ax(work, iterates->ax, iterates->x);
  Ax(work, iterates->ax_bar, iterates->x_bar);
  Ax(work, iterates->ax_ag, iterates->x_ag);

  ATy(work, iterates->aty, iterates->y);
  ATy(work, iterates->aty_ag, iterates->y_ag);

  CUPDLP_ZERO_VEC(iterates->x_agSum, cupdlp_float, lp->nCols);
  CUPDLP_ZERO_VEC(iterates->y_agSum, cupdlp_float, lp->nRows);
  CUPDLP_ZERO_VEC(iterates->x_agAverage->data, cupdlp_float, lp->nCols);
  CUPDLP_ZERO_VEC(iterates->y_agAverage->data, cupdlp_float, lp->nRows);

  PDHG_Project_Bounds(work, iterates->x_agSum);
  PDHG_Project_Bounds(work, iterates->x_agAverage->data);

  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;

  CUPDLP_ZERO_VEC(iterates->xLastRestart, cupdlp_float, lp->nCols);
  CUPDLP_ZERO_VEC(iterates->yLastRestart, cupdlp_float, lp->nRows);
}

void PDTEST_Init_Variables(CUPDLPwork *work)
{
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPstepsize *stepsize = work->stepsize;
  PDTESTiterates *iterates = work->PDTESTiterates;

  CUPDLP_ZERO_VEC(iterates->x->data, cupdlp_float, lp->nCols);
  PDHG_Project_Bounds(work, iterates->x->data);

  CUPDLP_ZERO_VEC(iterates->x_ag->data, cupdlp_float, lp->nCols);
  PDHG_Project_Bounds(work, iterates->x_ag->data);

  CUPDLP_ZERO_VEC(iterates->x_bar->data, cupdlp_float, lp->nCols);
  PDHG_Project_Bounds(work, iterates->x_bar->data);

  // CUPDLP_ZERO_VEC(iterates->x_md->data, cupdlp_float, lp->nCols);
  // PDHG_Project_Bounds(work, iterates->x_md->data);

  CUPDLP_ZERO_VEC(iterates->y->data, cupdlp_float, lp->nRows);
  CUPDLP_ZERO_VEC(iterates->y_ag->data, cupdlp_float, lp->nRows);

  Ax(work, iterates->ax, iterates->x);
  Ax(work, iterates->ax_bar, iterates->x_bar);
  Ax(work, iterates->ax_ag, iterates->x_ag);

  ATy(work, iterates->aty, iterates->y);
  ATy(work, iterates->aty_ag, iterates->y_ag);

  // CUPDLP_ZERO_VEC(iterates->xSum, cupdlp_float, lp->nCols);
  // CUPDLP_ZERO_VEC(iterates->ySum, cupdlp_float, lp->nRows);
  // CUPDLP_ZERO_VEC(iterates->xAverage->data, cupdlp_float, lp->nCols);
  // CUPDLP_ZERO_VEC(iterates->yAverage->data, cupdlp_float, lp->nRows);

  // PDHG_Project_Bounds(work, iterates->xSum);
  // PDHG_Project_Bounds(work, iterates->xAverage->data);

  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;
  stepsize->dStepSizeSum = 0.0;

  CUPDLP_ZERO_VEC(iterates->xLastRestart, cupdlp_float, lp->nCols);
  CUPDLP_ZERO_VEC(iterates->yLastRestart, cupdlp_float, lp->nRows);
}
/* TODO: this function seems considering
 *       l1 <= Ax <= u1
 *       l2 <=  x <= u2
 *       needs rewritten for current formulation
 */
void PDHG_Check_Data(CUPDLPwork *work)
{
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPstepsize *stepsize = work->stepsize;
  CUPDLPiterates *iterates = work->iterates;
  cupdlp_int nFreeCol = 0;
  cupdlp_int nFixedCol = 0;
  cupdlp_int nUpperCol = 0;
  cupdlp_int nLowerCol = 0;
  cupdlp_int nRangedCol = 0;
  cupdlp_int nFreeRow = 0;
  cupdlp_int nFixedRow = 0;
  cupdlp_int nUpperRow = 0;
  cupdlp_int nLowerRow = 0;
  cupdlp_int nRangedRow = 0;

  for (cupdlp_int iSeq = 0; iSeq < lp->nCols; ++iSeq)
  {
    cupdlp_bool hasLower = problem->lower[iSeq] > -INFINITY;
    cupdlp_bool hasUpper = problem->upper[iSeq] < +INFINITY;

    if (!hasLower && !hasUpper)
    {
      ++nFreeCol;
      cupdlp_printf("Warning: variable %d is free.", iSeq);
    }

    if (hasLower && hasUpper)
    {
      if (problem->lower[iSeq] == problem->upper[iSeq])
      {
        ++nFixedCol;
        // cupdlp_printf( "Warning: variable %d is fixed.", iSeq);
      }
      else
        ++nRangedCol;
    }

    if (hasLower)
    {
      // XXX: uncommented for PDLP comparison
      // CUPDLP_ASSERT(iterates->x[iSeq] >= problem->lower[iSeq]);
      nLowerCol += !hasUpper;
    }

    if (hasUpper)
    {
      // XXX: uncommented for PDLP comparison
      // CUPDLP_ASSERT(iterates->x[iSeq] <= problem->upper[iSeq]);
      nUpperCol += !hasLower;
    }
  }

  for (cupdlp_int iSeq = lp->nCols; iSeq < lp->nCols; ++iSeq)
  {
    cupdlp_bool hasLower = problem->lower[iSeq] > -INFINITY;
    cupdlp_bool hasUpper = problem->upper[iSeq] < +INFINITY;

    if (!hasLower && !hasUpper)
    {
      ++nFreeRow;
      cupdlp_printf("Warning: row %d is free.", iSeq - lp->nCols);
    }

    if (hasLower && hasUpper)
    {
      if (problem->lower[iSeq] == problem->upper[iSeq])
        ++nFixedRow;
      else
        ++nRangedRow;
    }

    if (hasLower)
    {
      // CUPDLP_ASSERT(iterates->x[iSeq] >= problem->lower[iSeq]);
      nLowerRow += !hasUpper;
    }

    if (hasUpper)
    {
      // CUPDLP_ASSERT(iterates->x[iSeq] <= problem->upper[iSeq]);
      nUpperRow += !hasLower;
    }
  }

  for (cupdlp_int iRow = 0; iRow < lp->nRows; ++iRow)
  {
    CUPDLP_ASSERT(iterates->y->data[iRow] < +INFINITY);
    CUPDLP_ASSERT(iterates->y->data[iRow] > -INFINITY);
  }

  for (cupdlp_int iRow = 0; iRow < lp->nRows; ++iRow)
  {
    if (problem->data->csr_matrix->rowMatBeg[iRow + 1] -
            problem->data->csr_matrix->rowMatBeg[iRow] ==
        1)
    {
      cupdlp_printf("Warning: row %d is a singleton row.", iRow);
    }
  }

  CUPDLP_ASSERT(nRangedRow == 0);
  cupdlp_printf("nFreeCol  : %d\n", nFreeCol);
  cupdlp_printf("nFixedCol : %d\n", nFixedCol);
  cupdlp_printf("nRangedCol: %d\n", nRangedCol);
  cupdlp_printf("nLowerCol : %d\n", nLowerCol);
  cupdlp_printf("nUpperCol : %d\n", nUpperCol);
  cupdlp_printf("nFreeRow  : %d\n", nFreeRow);
  cupdlp_printf("nFixedRow : %d\n", nFixedRow);
  cupdlp_printf("nRangedRow: %d\n", nRangedRow);
  cupdlp_printf("nLowerRow : %d\n", nLowerRow);
  cupdlp_printf("nUpperRow : %d\n", nUpperRow);

  // We need to test problems ranged row-bounds more carefully.
  CUPDLP_ASSERT(nRangedRow == 0);
}

cupdlp_bool PDTEST_Check_Termination_best(CUPDLPwork *pdhg, int bool_print)
{
  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPsettings *settings = pdhg->settings;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPscaling *scaling = pdhg->scaling;
  cupdlp_int bestID = settings->bestID;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;

  switch (bestID)
  {
  case 0:
    if (PDHG_Check_Termination(pdhg, bool_print))
    {
      cupdlp_printf("Optimal current solution.\n");
      resobj->termIterate = LAST_ITERATE;
      resobj->termCode = OPTIMAL;
      break;
    }
    break;
  case 1:
    if (PDHG_Check_Termination(pdhg, bool_print))
    {
      cupdlp_printf("Optimal current solution.\n");
      resobj->termIterate = LAST_ITERATE;
      resobj->termCode = OPTIMAL;
      break;
    }
    break;
  case 2:
    if (PDHG_Check_Termination(pdhg, bool_print))
    {
      cupdlp_printf("Optimal current solution.\n");
      resobj->termIterate = LAST_ITERATE;
      resobj->termCode = OPTIMAL;
      break;
    }
    break;
  case 3:
    if (PDHG_Check_Termination(pdhg, bool_print))
    {
      cupdlp_printf("Optimal current solution.\n");
      resobj->termIterate = LAST_ITERATE;
      resobj->termCode = OPTIMAL;
      break;
    }
    break;
  case 4:
    if (PDHG_Check_Termination(pdhg, bool_print))
    {
      cupdlp_printf("Optimal current solution.\n");
      resobj->termIterate = LAST_ITERATE;
      resobj->termCode = OPTIMAL;
      break;
    }
    break;
  case 5:
    if (PDHG_Check_Termination(pdhg, bool_print))
    {
      cupdlp_printf("Optimal current solution.\n");
      resobj->termIterate = LAST_ITERATE;
      resobj->termCode = OPTIMAL;
      break;
    }
    break;
  case 6:
    if (PDHG_Check_Termination(pdhg, bool_print))
    {
      cupdlp_printf("Optimal current solution.\n");
      resobj->termIterate = LAST_ITERATE;
      resobj->termCode = OPTIMAL;
      break;
    }

    if (PDHG_Check_Termination_Average(pdhg, bool_print))
    {
      cupdlp_printf("Optimal average solution.\n");

      CUPDLP_COPY_VEC(iterates->x->data, iterates->x_ag->data,
                      cupdlp_float, problem->nCols);
      CUPDLP_COPY_VEC(iterates->y->data, iterates->y_ag->data,
                      cupdlp_float, problem->nRows);

      resobj->termIterate = AVERAGE_ITERATE;
      resobj->termCode = OPTIMAL;
      break;
    }
    break;
  }
}

cupdlp_bool PDHG_Check_Termination(CUPDLPwork *pdhg, int bool_print)
{
  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPsettings *settings = pdhg->settings;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPscaling *scaling = pdhg->scaling;
#if PDHG_DISPLAY_TERMINATION_CHECK
  // todo, check, is it correct
  if (bool_print)
  {
    cupdlp_printf(
        "Termination check: %e|%e  %e|%e  %e|%e\n", resobj->dPrimalFeas,
        settings->dPrimalTol * (1.0 + scaling->dNormRhs), resobj->dDualFeas,
        settings->dDualTol * (1.0 + scaling->dNormCost), resobj->dRelObjGap,
        settings->dGapTol);
  }

#endif
  int bool_pass =
      ((resobj->dPrimalFeas <
        settings->dPrimalTol * (1.0 + scaling->dNormRhs)) &&
       (resobj->dDualFeas < settings->dDualTol * (1.0 + scaling->dNormCost)) &&
       (resobj->dRelObjGap < settings->dGapTol));
  return bool_pass;
}

cupdlp_bool PDHG_Check_Termination_Average(CUPDLPwork *pdhg, int bool_print)
{
  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPsettings *settings = pdhg->settings;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPscaling *scaling = pdhg->scaling;
#if PDHG_DISPLAY_TERMINATION_CHECK
  if (bool_print)
  {
    cupdlp_printf("Termination check: %e|%e  %e|%e  %e|%e\n",
                  resobj->dPrimalFeasAverage,
                  settings->dPrimalTol * (1.0 + scaling->dNormRhs),
                  resobj->dDualFeasAverage,
                  settings->dDualTol * (1.0 + scaling->dNormCost),
                  resobj->dRelObjGapAverage, settings->dGapTol);
  }
#endif
  int bool_pass = ((resobj->dPrimalFeasAverage <
                    settings->dPrimalTol * (1.0 + scaling->dNormRhs)) &&
                   (resobj->dDualFeasAverage <
                    settings->dDualTol * (1.0 + scaling->dNormCost)) &&
                   (resobj->dRelObjGapAverage < settings->dGapTol));
  return bool_pass;
}

void PDHG_Print_Header(CUPDLPwork *pdhg)
{
  cupdlp_printf("%5s  %15s  %15s   %8s  %8s  %10s  %8s %7s\n", "Iter",
                "Primal.Obj", "Dual.Obj", "Gap", "Compl", "Primal.Inf",
                "Dual.Inf", "Time");
}

void PDHG_Print_Iter(CUPDLPwork *pdhg)
{
  /* Format time as xxx.yy for < 1000s and as integer afterwards. */
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPtimers *timers = pdhg->timers;
  char timeString[8];
  if (timers->dSolvingTime < 100.0)
    cupdlp_snprintf(timeString, 8, "%6.2fs", timers->dSolvingTime);
  else
    cupdlp_snprintf(timeString, 8, "%6ds", (cupdlp_int)timers->dSolvingTime);

  cupdlp_printf("%5d  %+15.8e  %+15.8e  %+8.2e  %8.2e  %10.2e  %8.2e %7s [L]\n",
                timers->nIter, resobj->dPrimalObj, resobj->dDualObj,
                resobj->dDualityGap, resobj->dComplementarity,
                resobj->dPrimalFeas, resobj->dDualFeas, timeString);
}

void PDHG_Print_Iter_Average(CUPDLPwork *pdhg)
{
  /* Format time as xxx.yy for < 1000s and as integer afterwards. */
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPtimers *timers = pdhg->timers;
  char timeString[8];
  if (timers->dSolvingTime < 100.0)
    cupdlp_snprintf(timeString, 8, "%6.2fs", timers->dSolvingTime);
  else
    cupdlp_snprintf(timeString, 8, "%6ds", (cupdlp_int)timers->dSolvingTime);

  cupdlp_printf("%5d  %+15.8e  %+15.8e  %+8.2e  %8.2e  %10.2e  %8.2e %7s [A]\n",
                timers->nIter, resobj->dPrimalObjAverage,
                resobj->dDualObjAverage, resobj->dDualityGapAverage,
                resobj->dComplementarityAverage, resobj->dPrimalFeasAverage,
                resobj->dDualFeasAverage, timeString);
}

void PDHG_Compute_SolvingTime(CUPDLPwork *pdhg)
{
  CUPDLPtimers *timers = pdhg->timers;
  timers->dSolvingTime = getTimeStamp() - timers->dSolvingBeg;
}
cupdlp_retcode PDHG_Solve_Multiscale_withStepsize(CUPDLPwork *pdhg, cupdlp_float *x_init, cupdlp_float *y_init, cupdlp_float *stepsize_last, cupdlp_float *weight_last, cupdlp_float *stepsize_init, cupdlp_float *weight_init, cupdlp_bool whether_first)
{
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPsettings *settings = pdhg->settings;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPtimers *timers = pdhg->timers;

  cupdlp_float dCheckTime = 0.0;
  cupdlp_float dCheckTimetemp = 0.0;
  cupdlp_float dUpdateTime = 0.0;
  cupdlp_float dUpdateTimetemp = 0.0;
  cupdlp_float dCheckTerminationTime = 0.0;
  cupdlp_float dCheckTerminationTimetemp = 0.0;
  cupdlp_float dRestartTime = 0.0;
  cupdlp_float dRestartTimetemp = 0.0;

  timers->dSolvingBeg = getTimeStamp();

  PDHG_Init_Data(pdhg);

  CUPDLP_CALL(PDHG_Init_Step_Sizes_Multiscale_withStepsize(pdhg, stepsize_init, weight_init, whether_first));

  PDHG_Init_Variables_Multiscale(pdhg, x_init, y_init);

  // todo: translate check_data into cuda or do it on cpu
  // PDHG_Check_Data(pdhg);

  // PDHG_Print_Header(pdhg);

  for (timers->nIter = 0; timers->nIter < settings->nIterLim; ++timers->nIter)
  {
    dCheckTimetemp = getTimeStamp();
    PDHG_Compute_SolvingTime(pdhg);
#if CUPDLP_DUMP_ITERATES_STATS & CUPDLP_DEBUG
    PDHG_Dump_Stats(pdhg);
#endif
    int bool_checking = (timers->nIter < 10) ||
                        (timers->nIter == (settings->nIterLim - 1)) ||
                        (timers->dSolvingTime > settings->dTimeLim);
    int bool_print = 0;
#if CUPDLP_DEBUG
    bool_checking = (bool_checking || !(timers->nIter % CUPDLP_DEBUG_INTERVAL));
    bool_print = bool_checking;
#else
    bool_checking =
        (bool_checking || !(timers->nIter % CUPDLP_RELEASE_INTERVAL));
    bool_print =
        (bool_checking && !(timers->nIter % (CUPDLP_RELEASE_INTERVAL *
                                             settings->nLogInterval))) ||
        (timers->nIter == (settings->nIterLim - 1)) ||
        (timers->dSolvingTime > settings->dTimeLim);
#endif
    if (bool_checking)
    {
      PDHG_Compute_Average_Iterate(pdhg);
      PDHG_Compute_Residuals(pdhg);
      if (bool_print)
      {
        PDHG_Print_Header(pdhg);
        PDHG_Print_Iter(pdhg);
        PDHG_Print_Iter_Average(pdhg);
      }
      dCheckTerminationTimetemp = getTimeStamp();
      if (PDHG_Check_Termination(pdhg, bool_print))
      {
        cupdlp_printf("Optimal current solution.\n");
        resobj->termIterate = LAST_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }

      if (PDHG_Check_Termination_Average(pdhg, bool_print))
      {
        cupdlp_printf("Optimal average solution.\n");

        CUPDLP_COPY_VEC(iterates->x->data, iterates->xAverage->data,
                        cupdlp_float, problem->nCols);
        CUPDLP_COPY_VEC(iterates->y->data, iterates->yAverage->data,
                        cupdlp_float, problem->nRows);

        resobj->termIterate = AVERAGE_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }

      if (timers->dSolvingTime > settings->dTimeLim)
      {
        cupdlp_printf("Time limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }

      if (timers->nIter == (settings->nIterLim - 1))
      {
        cupdlp_printf("Iteration limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }
      dCheckTerminationTime += getTimeStamp() - dCheckTerminationTimetemp;
      dRestartTimetemp = getTimeStamp();
      PDHG_Restart_Iterate(pdhg);
      dRestartTime += getTimeStamp() - dRestartTimetemp;
      dCheckTime += getTimeStamp() - dCheckTimetemp;
    }
    dUpdateTimetemp = getTimeStamp();
    CUPDLP_CALL(PDHG_Update_Iterate(pdhg));
    dUpdateTime += getTimeStamp() - dUpdateTimetemp;
  }
  cupdlp_printf("Check time: %e\n", dCheckTime);
  cupdlp_printf("Update time: %e\n", dUpdateTime);
  cupdlp_printf("Check termination time: %e\n", dCheckTerminationTime);
  cupdlp_printf("Restart time: %e\n", dRestartTime);
  cupdlp_printf("Compute Update time: %e\n", timers->dIterTime);
  // print at last
  PDHG_Print_Header(pdhg);
  PDHG_Print_Iter(pdhg);
  PDHG_Print_Iter_Average(pdhg);

  // 保存最后一步迭代的步长和权重
  *stepsize_last = sqrt(stepsize->dPrimalStep * stepsize->dDualStep);
  *weight_last = sqrt(stepsize->dBeta);

#if PDHG_USE_TIMERS
  cupdlp_printf("Timing information:\n");
  // cupdlp_printf("%20s %e in %d iterations\n", "Total solver time",
  //               timers->dSolvingTime, timers->nIter);
  cupdlp_printf(
      "%20s %e in %d iterations\n", "Total solver time",
      timers->dSolvingTime + timers->dScalingTime + timers->dPresolveTime,
      timers->nIter);
  cupdlp_printf("%20s %e in %d iterations\n", "Solve time",
                timers->dSolvingTime, timers->nIter);
  cupdlp_printf("%20s %e \n", "Iters per sec",
                timers->nIter / timers->dSolvingTime);
  cupdlp_printf("%20s %e\n", "Scaling time", timers->dScalingTime);
  cupdlp_printf("%20s %e\n", "Presolve time", timers->dPresolveTime);
  cupdlp_printf("%20s %e in %d calls\n", "Ax", timers->dAxTime,
                timers->nAxCalls);
  cupdlp_printf("%20s %e in %d calls\n", "Aty", timers->dAtyTime,
                timers->nAtyCalls);
  cupdlp_printf("%20s %e in %d calls\n", "ComputeResiduals",
                timers->dComputeResidualsTime, timers->nComputeResidualsCalls);
  cupdlp_printf("%20s %e in %d calls\n", "UpdateIterates",
                timers->dUpdateIterateTime, timers->nUpdateIterateCalls);
#endif

#if !(CUPDLP_CPU)
  cupdlp_printf("GPU Timing information:\n");
  cupdlp_printf("%20s %e\n", "CudaPrepare", timers->CudaPrepareTime);
  cupdlp_printf("%20s %e\n", "Alloc&CopyMatToDevice",
                timers->AllocMem_CopyMatToDeviceTime);
  cupdlp_printf("%20s %e\n", "CopyVecToDevice", timers->CopyVecToDeviceTime);
  cupdlp_printf("%20s %e\n", "DeviceMatVecProd", timers->DeviceMatVecProdTime);
  cupdlp_printf("%20s %e\n", "CopyVecToHost", timers->CopyVecToHostTime);
#endif

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDHG_Solve_Multiscale(CUPDLPwork *pdhg, cupdlp_float *x_init, cupdlp_float *y_init)
{
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPsettings *settings = pdhg->settings;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPtimers *timers = pdhg->timers;

  cupdlp_float dCheckTime = 0.0;
  cupdlp_float dCheckTimetemp = 0.0;
  cupdlp_float dUpdateTime = 0.0;
  cupdlp_float dUpdateTimetemp = 0.0;
  cupdlp_float dCheckTerminationTime = 0.0;
  cupdlp_float dCheckTerminationTimetemp = 0.0;
  cupdlp_float dRestartTime = 0.0;
  cupdlp_float dRestartTimetemp = 0.0;

  timers->dSolvingBeg = getTimeStamp();

  PDHG_Init_Data(pdhg);

  CUPDLP_CALL(PDHG_Init_Step_Sizes_Multiscale(pdhg));

  PDHG_Init_Variables_Multiscale(pdhg, x_init, y_init);

  // todo: translate check_data into cuda or do it on cpu
  // PDHG_Check_Data(pdhg);

  // PDHG_Print_Header(pdhg);

  for (timers->nIter = 0; timers->nIter < settings->nIterLim; ++timers->nIter)
  {
    dCheckTimetemp = getTimeStamp();
    PDHG_Compute_SolvingTime(pdhg);
#if CUPDLP_DUMP_ITERATES_STATS & CUPDLP_DEBUG
    PDHG_Dump_Stats(pdhg);
#endif
    int bool_checking = (timers->nIter < 10) ||
                        (timers->nIter == (settings->nIterLim - 1)) ||
                        (timers->dSolvingTime > settings->dTimeLim);
    int bool_print = 0;
#if CUPDLP_DEBUG
    bool_checking = (bool_checking || !(timers->nIter % CUPDLP_DEBUG_INTERVAL));
    bool_print = bool_checking;
#else
    bool_checking =
        (bool_checking || !(timers->nIter % CUPDLP_RELEASE_INTERVAL));
    bool_print =
        (bool_checking && !(timers->nIter % (CUPDLP_RELEASE_INTERVAL *
                                             settings->nLogInterval))) ||
        (timers->nIter == (settings->nIterLim - 1)) ||
        (timers->dSolvingTime > settings->dTimeLim);
#endif
    if (bool_checking)
    {
      PDHG_Compute_Average_Iterate(pdhg);
      PDHG_Compute_Residuals(pdhg);
      if (bool_print)
      {
        PDHG_Print_Header(pdhg);
        PDHG_Print_Iter(pdhg);
        PDHG_Print_Iter_Average(pdhg);
      }
      dCheckTerminationTimetemp = getTimeStamp();
      if (PDHG_Check_Termination(pdhg, bool_print))
      {
        cupdlp_printf("Optimal current solution.\n");
        resobj->termIterate = LAST_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }

      if (PDHG_Check_Termination_Average(pdhg, bool_print))
      {
        cupdlp_printf("Optimal average solution.\n");

        CUPDLP_COPY_VEC(iterates->x->data, iterates->xAverage->data,
                        cupdlp_float, problem->nCols);
        CUPDLP_COPY_VEC(iterates->y->data, iterates->yAverage->data,
                        cupdlp_float, problem->nRows);

        resobj->termIterate = AVERAGE_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }

      if (timers->dSolvingTime > settings->dTimeLim)
      {
        cupdlp_printf("Time limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }

      if (timers->nIter == (settings->nIterLim - 1))
      {
        cupdlp_printf("Iteration limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }
      dCheckTerminationTime += getTimeStamp() - dCheckTerminationTimetemp;
      dRestartTimetemp = getTimeStamp();
      PDHG_Restart_Iterate(pdhg);
      dRestartTime += getTimeStamp() - dRestartTimetemp;
      dCheckTime += getTimeStamp() - dCheckTimetemp;
    }
    dUpdateTimetemp = getTimeStamp();
    CUPDLP_CALL(PDHG_Update_Iterate(pdhg));
    dUpdateTime += getTimeStamp() - dUpdateTimetemp;
  }
  cupdlp_printf("Check time: %e\n", dCheckTime);
  cupdlp_printf("Update time: %e\n", dUpdateTime);
  cupdlp_printf("Check termination time: %e\n", dCheckTerminationTime);
  cupdlp_printf("Restart time: %e\n", dRestartTime);
  cupdlp_printf("Compute Update time: %e\n", timers->dIterTime);
  // print at last
  PDHG_Print_Header(pdhg);
  PDHG_Print_Iter(pdhg);
  PDHG_Print_Iter_Average(pdhg);

#if PDHG_USE_TIMERS
  cupdlp_printf("Timing information:\n");
  // cupdlp_printf("%20s %e in %d iterations\n", "Total solver time",
  //               timers->dSolvingTime, timers->nIter);
  cupdlp_printf(
      "%20s %e in %d iterations\n", "Total solver time",
      timers->dSolvingTime + timers->dScalingTime + timers->dPresolveTime,
      timers->nIter);
  cupdlp_printf("%20s %e in %d iterations\n", "Solve time",
                timers->dSolvingTime, timers->nIter);
  cupdlp_printf("%20s %e \n", "Iters per sec",
                timers->nIter / timers->dSolvingTime);
  cupdlp_printf("%20s %e\n", "Scaling time", timers->dScalingTime);
  cupdlp_printf("%20s %e\n", "Presolve time", timers->dPresolveTime);
  cupdlp_printf("%20s %e in %d calls\n", "Ax", timers->dAxTime,
                timers->nAxCalls);
  cupdlp_printf("%20s %e in %d calls\n", "Aty", timers->dAtyTime,
                timers->nAtyCalls);
  cupdlp_printf("%20s %e in %d calls\n", "ComputeResiduals",
                timers->dComputeResidualsTime, timers->nComputeResidualsCalls);
  cupdlp_printf("%20s %e in %d calls\n", "UpdateIterates",
                timers->dUpdateIterateTime, timers->nUpdateIterateCalls);
#endif

#if !(CUPDLP_CPU)
  cupdlp_printf("GPU Timing information:\n");
  cupdlp_printf("%20s %e\n", "CudaPrepare", timers->CudaPrepareTime);
  cupdlp_printf("%20s %e\n", "Alloc&CopyMatToDevice",
                timers->AllocMem_CopyMatToDeviceTime);
  cupdlp_printf("%20s %e\n", "CopyVecToDevice", timers->CopyVecToDeviceTime);
  cupdlp_printf("%20s %e\n", "DeviceMatVecProd", timers->DeviceMatVecProdTime);
  cupdlp_printf("%20s %e\n", "CopyVecToHost", timers->CopyVecToHostTime);
#endif

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDHG_Solve(CUPDLPwork *pdhg)
{
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPsettings *settings = pdhg->settings;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPtimers *timers = pdhg->timers;

  cupdlp_float dCheckTime = 0.0;
  cupdlp_float dCheckTimetemp = 0.0;
  cupdlp_float dUpdateTime = 0.0;
  cupdlp_float dUpdateTimetemp = 0.0;
  cupdlp_float dCheckTerminationTime = 0.0;
  cupdlp_float dCheckTerminationTimetemp = 0.0;
  cupdlp_float dRestartTime = 0.0;
  cupdlp_float dRestartTimetemp = 0.0;

  timers->dSolvingBeg = getTimeStamp();

  PDHG_Init_Data(pdhg);

  CUPDLP_CALL(PDHG_Init_Step_Sizes(pdhg));

  PDHG_Init_Variables(pdhg);

  // todo: translate check_data into cuda or do it on cpu
  // PDHG_Check_Data(pdhg);

  // PDHG_Print_Header(pdhg);

  for (timers->nIter = 0; timers->nIter < settings->nIterLim; ++timers->nIter)
  {
    dCheckTimetemp = getTimeStamp();
    PDHG_Compute_SolvingTime(pdhg);
#if CUPDLP_DUMP_ITERATES_STATS & CUPDLP_DEBUG
    PDHG_Dump_Stats(pdhg);
#endif
    int bool_checking = (timers->nIter < 10) ||
                        (timers->nIter == (settings->nIterLim - 1)) ||
                        (timers->dSolvingTime > settings->dTimeLim);
    int bool_print = 0;
#if CUPDLP_DEBUG
    bool_checking = (bool_checking || !(timers->nIter % CUPDLP_DEBUG_INTERVAL));
    bool_print = bool_checking;
#else
    bool_checking =
        (bool_checking || !(timers->nIter % CUPDLP_RELEASE_INTERVAL));
    bool_print =
        (bool_checking && !(timers->nIter % (CUPDLP_RELEASE_INTERVAL *
                                             settings->nLogInterval))) ||
        (timers->nIter == (settings->nIterLim - 1)) ||
        (timers->dSolvingTime > settings->dTimeLim);
#endif
    if (bool_checking)
    {
      PDHG_Compute_Average_Iterate(pdhg);
      PDHG_Compute_Residuals(pdhg);
      if (bool_print)
      {
        PDHG_Print_Header(pdhg);
        PDHG_Print_Iter(pdhg);
        PDHG_Print_Iter_Average(pdhg);
      }
      dCheckTerminationTimetemp = getTimeStamp();
      if (PDHG_Check_Termination(pdhg, bool_print))
      {
        cupdlp_printf("Optimal current solution.\n");
        resobj->termIterate = LAST_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }

      if (PDHG_Check_Termination_Average(pdhg, bool_print))
      {
        cupdlp_printf("Optimal average solution.\n");

        CUPDLP_COPY_VEC(iterates->x->data, iterates->xAverage->data,
                        cupdlp_float, problem->nCols);
        CUPDLP_COPY_VEC(iterates->y->data, iterates->yAverage->data,
                        cupdlp_float, problem->nRows);

        resobj->termIterate = AVERAGE_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }

      if (timers->dSolvingTime > settings->dTimeLim)
      {
        cupdlp_printf("Time limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }

      if (timers->nIter == (settings->nIterLim - 1))
      {
        cupdlp_printf("Iteration limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }
      dCheckTerminationTime += getTimeStamp() - dCheckTerminationTimetemp;
      dRestartTimetemp = getTimeStamp();
      PDHG_Restart_Iterate(pdhg);
      dRestartTime += getTimeStamp() - dRestartTimetemp;
      dCheckTime += getTimeStamp() - dCheckTimetemp;
    }
    dUpdateTimetemp = getTimeStamp();
    CUPDLP_CALL(PDHG_Update_Iterate(pdhg));
    dUpdateTime += getTimeStamp() - dUpdateTimetemp;
  }
  cupdlp_printf("Check time: %e\n", dCheckTime);
  cupdlp_printf("Update time: %e\n", dUpdateTime);
  cupdlp_printf("Check termination time: %e\n", dCheckTerminationTime);
  cupdlp_printf("Restart time: %e\n", dRestartTime);
  cupdlp_printf("Compute Update time: %e\n", timers->dIterTime);
  // print at last
  PDHG_Print_Header(pdhg);
  PDHG_Print_Iter(pdhg);
  PDHG_Print_Iter_Average(pdhg);

#if PDHG_USE_TIMERS
  cupdlp_printf("Timing information:\n");
  // cupdlp_printf("%20s %e in %d iterations\n", "Total solver time",
  //               timers->dSolvingTime, timers->nIter);
  cupdlp_printf(
      "%20s %e in %d iterations\n", "Total solver time",
      timers->dSolvingTime + timers->dScalingTime + timers->dPresolveTime,
      timers->nIter);
  cupdlp_printf("%20s %e in %d iterations\n", "Solve time",
                timers->dSolvingTime, timers->nIter);
  cupdlp_printf("%20s %e \n", "Iters per sec",
                timers->nIter / timers->dSolvingTime);
  cupdlp_printf("%20s %e\n", "Scaling time", timers->dScalingTime);
  cupdlp_printf("%20s %e\n", "Presolve time", timers->dPresolveTime);
  cupdlp_printf("%20s %e in %d calls\n", "Ax", timers->dAxTime,
                timers->nAxCalls);
  cupdlp_printf("%20s %e in %d calls\n", "Aty", timers->dAtyTime,
                timers->nAtyCalls);
  cupdlp_printf("%20s %e in %d calls\n", "ComputeResiduals",
                timers->dComputeResidualsTime, timers->nComputeResidualsCalls);
  cupdlp_printf("%20s %e in %d calls\n", "UpdateIterates",
                timers->dUpdateIterateTime, timers->nUpdateIterateCalls);
#endif

#if !(CUPDLP_CPU)
  cupdlp_printf("GPU Timing information:\n");
  cupdlp_printf("%20s %e\n", "CudaPrepare", timers->CudaPrepareTime);
  cupdlp_printf("%20s %e\n", "Alloc&CopyMatToDevice",
                timers->AllocMem_CopyMatToDeviceTime);
  cupdlp_printf("%20s %e\n", "CopyVecToDevice", timers->CopyVecToDeviceTime);
  cupdlp_printf("%20s %e\n", "DeviceMatVecProd", timers->DeviceMatVecProdTime);
  cupdlp_printf("%20s %e\n", "CopyVecToHost", timers->CopyVecToHostTime);
#endif

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDHG_Solve_AdapTheta(CUPDLPwork *pdhg)
{
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPsettings *settings = pdhg->settings;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPtimers *timers = pdhg->timers;

  cupdlp_float dCheckTime = 0.0;
  cupdlp_float dCheckTimetemp = 0.0;
  cupdlp_float dUpdateTime = 0.0;
  cupdlp_float dUpdateTimetemp = 0.0;
  cupdlp_float dCheckTerminationTime = 0.0;
  cupdlp_float dCheckTerminationTimetemp = 0.0;
  cupdlp_float dRestartTime = 0.0;
  cupdlp_float dRestartTimetemp = 0.0;

  timers->dSolvingBeg = getTimeStamp();

  PDHG_Init_Data(pdhg);

  CUPDLP_CALL(PDHG_Init_Step_Sizes_AdapTheta(pdhg));

  PDHG_Init_Variables(pdhg);

  // todo: translate check_data into cuda or do it on cpu
  // PDHG_Check_Data(pdhg);

  // PDHG_Print_Header(pdhg);

  for (timers->nIter = 0; timers->nIter < settings->nIterLim; ++timers->nIter)
  {
    dCheckTimetemp = getTimeStamp();
    PDHG_Compute_SolvingTime(pdhg);
#if CUPDLP_DUMP_ITERATES_STATS & CUPDLP_DEBUG
    PDHG_Dump_Stats(pdhg);
#endif
    int bool_checking = (timers->nIter < 10) ||
                        (timers->nIter == (settings->nIterLim - 1)) ||
                        (timers->dSolvingTime > settings->dTimeLim);
    int bool_print = 0;
#if CUPDLP_DEBUG
    bool_checking = (bool_checking || !(timers->nIter % CUPDLP_DEBUG_INTERVAL));
    bool_print = bool_checking;
#else
    bool_checking =
        (bool_checking || !(timers->nIter % CUPDLP_RELEASE_INTERVAL));
    bool_print =
        (bool_checking && !(timers->nIter % (CUPDLP_RELEASE_INTERVAL *
                                             settings->nLogInterval))) ||
        (timers->nIter == (settings->nIterLim - 1)) ||
        (timers->dSolvingTime > settings->dTimeLim);
#endif
    if (bool_checking)
    {
      PDHG_Compute_Average_Iterate(pdhg);
      PDHG_Compute_Residuals(pdhg);
      if (bool_print)
      {
        PDHG_Print_Header(pdhg);
        PDHG_Print_Iter(pdhg);
        PDHG_Print_Iter_Average(pdhg);
      }
      dCheckTerminationTimetemp = getTimeStamp();
      if (PDHG_Check_Termination(pdhg, bool_print))
      {
        cupdlp_printf("Optimal current solution.\n");
        resobj->termIterate = LAST_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }

      if (PDHG_Check_Termination_Average(pdhg, bool_print))
      {
        cupdlp_printf("Optimal average solution.\n");

        CUPDLP_COPY_VEC(iterates->x->data, iterates->xAverage->data,
                        cupdlp_float, problem->nCols);
        CUPDLP_COPY_VEC(iterates->y->data, iterates->yAverage->data,
                        cupdlp_float, problem->nRows);

        resobj->termIterate = AVERAGE_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }

      if (timers->dSolvingTime > settings->dTimeLim)
      {
        cupdlp_printf("Time limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }

      if (timers->nIter == (settings->nIterLim - 1))
      {
        cupdlp_printf("Iteration limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }
      dCheckTerminationTime += getTimeStamp() - dCheckTerminationTimetemp;
      dRestartTimetemp = getTimeStamp();
      PDHG_Restart_Iterate(pdhg);
      dRestartTime += getTimeStamp() - dRestartTimetemp;
      dCheckTime += getTimeStamp() - dCheckTimetemp;
    }
    dUpdateTimetemp = getTimeStamp();
    CUPDLP_CALL(PDHG_Update_Iterate_AdapTheta(pdhg));
    dUpdateTime += getTimeStamp() - dUpdateTimetemp;
  }
  cupdlp_printf("Check time: %e\n", dCheckTime);
  cupdlp_printf("Update time: %e\n", dUpdateTime);
  cupdlp_printf("Check termination time: %e\n", dCheckTerminationTime);
  cupdlp_printf("Restart time: %e\n", dRestartTime);
  cupdlp_printf("Compute Update time: %e\n", timers->dIterTime);
  // print at last
  PDHG_Print_Header(pdhg);
  PDHG_Print_Iter(pdhg);
  PDHG_Print_Iter_Average(pdhg);

#if PDHG_USE_TIMERS
  cupdlp_printf("Timing information:\n");
  // cupdlp_printf("%20s %e in %d iterations\n", "Total solver time",
  //               timers->dSolvingTime, timers->nIter);
  cupdlp_printf(
      "%20s %e in %d iterations\n", "Total solver time",
      timers->dSolvingTime + timers->dScalingTime + timers->dPresolveTime,
      timers->nIter);
  cupdlp_printf("%20s %e in %d iterations\n", "Solve time",
                timers->dSolvingTime, timers->nIter);
  cupdlp_printf("%20s %e \n", "Iters per sec",
                timers->nIter / timers->dSolvingTime);
  cupdlp_printf("%20s %e\n", "Scaling time", timers->dScalingTime);
  cupdlp_printf("%20s %e\n", "Presolve time", timers->dPresolveTime);
  cupdlp_printf("%20s %e in %d calls\n", "Ax", timers->dAxTime,
                timers->nAxCalls);
  cupdlp_printf("%20s %e in %d calls\n", "Aty", timers->dAtyTime,
                timers->nAtyCalls);
  cupdlp_printf("%20s %e in %d calls\n", "ComputeResiduals",
                timers->dComputeResidualsTime, timers->nComputeResidualsCalls);
  cupdlp_printf("%20s %e in %d calls\n", "UpdateIterates",
                timers->dUpdateIterateTime, timers->nUpdateIterateCalls);
#endif

#if !(CUPDLP_CPU)
  cupdlp_printf("GPU Timing information:\n");
  cupdlp_printf("%20s %e\n", "CudaPrepare", timers->CudaPrepareTime);
  cupdlp_printf("%20s %e\n", "Alloc&CopyMatToDevice",
                timers->AllocMem_CopyMatToDeviceTime);
  cupdlp_printf("%20s %e\n", "CopyVecToDevice", timers->CopyVecToDeviceTime);
  cupdlp_printf("%20s %e\n", "DeviceMatVecProd", timers->DeviceMatVecProdTime);
  cupdlp_printf("%20s %e\n", "CopyVecToHost", timers->CopyVecToHostTime);
#endif

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDTEST_Average_Solve(CUPDLPwork *pdhg)
{
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLPproblem *problem = pdhg->problem;    // 需要求解的问题
  CUPDLPstepsize *stepsize = pdhg->stepsize; // 步长
  CUPDLPsettings *settings = pdhg->settings; // 设置，包括scaling, termination criteria, max iter and time, restart

  CUPDLPresobj *resobj = pdhg->resobj;             /* residuals and objectives 残差和目标函数值*/
  PDTESTiterates *iterates = pdhg->PDTESTiterates; // 迭代变量
  CUPDLPtimers *timers = pdhg->timers;             // 各种耗时
  cupdlp_float dCheckTime = 0.0;
  cupdlp_float dCheckTimetemp = 0.0;
  cupdlp_float dUpdateTime = 0.0;
  cupdlp_float dUpdateTimetemp = 0.0;

  timers->dSolvingBeg = getTimeStamp();

  // PDHG_Init_Data(pdhg);

  CUPDLP_CALL(PDTEST_Init_Step_Sizes(pdhg));

  PDTEST_Average_Init_Variables(pdhg);

  // todo: translate check_data into cuda or do it on cpu
  // PDHG_Check_Data(pdhg);

  // PDHG_Print_Header(pdhg);
  // 主要的迭代！！
  cupdlp_int nIter_restart_value = 0;
  cupdlp_int *nIter_restart = &nIter_restart_value;
  for (timers->nIter = 0; timers->nIter < settings->nIterLim; ++timers->nIter)
  {
    PDHG_Compute_SolvingTime(pdhg);
#if CUPDLP_DUMP_ITERATES_STATS & CUPDLP_DEBUG
    PDTEST_Dump_Stats(pdhg);
#endif
    int bool_checking = (timers->nIter < 10) ||
                        (timers->nIter == (settings->nIterLim - 1)) ||
                        (timers->dSolvingTime > settings->dTimeLim);
    int bool_print = 0;
#if CUPDLP_DEBUG
    bool_checking = (bool_checking || !(timers->nIter % CUPDLP_DEBUG_INTERVAL));
    bool_print = bool_checking;
#else
    // 检查：迭代次数<10/迭代次数到了/迭代时间到了/迭代次数=0(mod40)
    bool_checking =
        (bool_checking || !(timers->nIter % CUPDLP_RELEASE_INTERVAL));
    // 打印：检查且迭代次数=0(mod特定值)/迭代次数到了/迭代时间到了
    bool_print =
        (bool_checking && !(timers->nIter % (CUPDLP_RELEASE_INTERVAL *
                                             settings->nLogInterval))) ||
        (timers->nIter == (settings->nIterLim - 1)) ||
        (timers->dSolvingTime > settings->dTimeLim);
#endif
    if (bool_checking)
    {
      dCheckTimetemp = getTimeStamp();
      PDTEST_Compute_Average_Iterate(pdhg);
      PDTEST_Average_Compute_Residuals(pdhg);
      if (bool_print)
      {
        PDHG_Print_Header(pdhg);
        PDHG_Print_Iter(pdhg);
        PDHG_Print_Iter_Average(pdhg);
      }

      if (PDHG_Check_Termination(pdhg, bool_print))
      {
        cupdlp_printf("Optimal current solution.\n");
        resobj->termIterate = LAST_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }

      if (PDHG_Check_Termination_Average(pdhg, bool_print))
      {
        cupdlp_printf("Optimal average solution.\n");

        CUPDLP_COPY_VEC(iterates->x_ag->data, iterates->x_agAverage->data,
                        cupdlp_float, problem->nCols);
        CUPDLP_COPY_VEC(iterates->y_ag->data, iterates->y_agAverage->data,
                        cupdlp_float, problem->nRows);

        resobj->termIterate = AVERAGE_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }

      if (timers->dSolvingTime > settings->dTimeLim)
      {
        cupdlp_printf("Time limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }

      if (timers->nIter == (settings->nIterLim - 1))
      {
        cupdlp_printf("Iteration limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }

      PDTEST_Average_Restart_Iterate(pdhg, nIter_restart); // restart策略
      dCheckTime += getTimeStamp() - dCheckTimetemp;
    }
    dUpdateTimetemp = getTimeStamp();
    CUPDLP_CALL(PDTEST_Average_Update_Iterate(pdhg, nIter_restart)); // 迭代更新
    *nIter_restart += 1;
    dUpdateTime += getTimeStamp() - dUpdateTimetemp;
  }
  cupdlp_printf("Check time: %e\n", dCheckTime);
  cupdlp_printf("Update time: %e\n", dUpdateTime);
  // print at last
  PDHG_Print_Header(pdhg);
  PDHG_Print_Iter(pdhg);
  PDHG_Print_Iter_Average(pdhg);

#if PDHG_USE_TIMERS
  cupdlp_printf("Timing information:\n");
  // cupdlp_printf("%20s %e in %d iterations\n", "Total solver time",
  //               timers->dSolvingTime, timers->nIter);
  cupdlp_printf(
      "%20s %e in %d iterations\n", "Total solver time",
      timers->dSolvingTime + timers->dScalingTime + timers->dPresolveTime,
      timers->nIter);
  cupdlp_printf("%20s %e in %d iterations\n", "Solve time",
                timers->dSolvingTime, timers->nIter);
  cupdlp_printf("%20s %e \n", "Iters per sec",
                timers->nIter / timers->dSolvingTime);
  cupdlp_printf("%20s %e\n", "Scaling time", timers->dScalingTime);
  cupdlp_printf("%20s %e\n", "Presolve time", timers->dPresolveTime);
  cupdlp_printf("%20s %e in %d calls\n", "Ax", timers->dAxTime,
                timers->nAxCalls);
  cupdlp_printf("%20s %e in %d calls\n", "Aty", timers->dAtyTime,
                timers->nAtyCalls);
  cupdlp_printf("%20s %e in %d calls\n", "ComputeResiduals",
                timers->dComputeResidualsTime, timers->nComputeResidualsCalls);
  cupdlp_printf("%20s %e in %d calls\n", "UpdateIterates",
                timers->dUpdateIterateTime, timers->nUpdateIterateCalls);
  cupdlp_printf("dMatVecMultiplyTime: %e\n", timers->dMatVecMultiplyTime);
  cupdlp_printf("dVecVecAddTime: %e\n", timers->dVecVecAddTime);
  cupdlp_printf("dIterTime: %e\n", timers->dIterTime);
#endif

#if !(CUPDLP_CPU)
  cupdlp_printf("GPU Timing information:\n");
  cupdlp_printf("%20s %e\n", "CudaPrepare", timers->CudaPrepareTime);
  cupdlp_printf("%20s %e\n", "Alloc&CopyMatToDevice",
                timers->AllocMem_CopyMatToDeviceTime);
  cupdlp_printf("%20s %e\n", "CopyVecToDevice", timers->CopyVecToDeviceTime);
  cupdlp_printf("%20s %e\n", "DeviceMatVecProd", timers->DeviceMatVecProdTime);
  cupdlp_printf("%20s %e\n", "CopyVecToHost", timers->CopyVecToHostTime);
#endif

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDTEST_Solve(CUPDLPwork *pdhg)
{
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLPproblem *problem = pdhg->problem;    // 需要求解的问题
  CUPDLPstepsize *stepsize = pdhg->stepsize; // 步长
  CUPDLPsettings *settings = pdhg->settings; // 设置，包括scaling, termination criteria, max iter and time, restart

  CUPDLPresobj *resobj = pdhg->resobj;             /* residuals and objectives 残差和目标函数值*/
  PDTESTiterates *iterates = pdhg->PDTESTiterates; // 迭代变量
  CUPDLPtimers *timers = pdhg->timers;             // 各种耗时

  timers->dSolvingBeg = getTimeStamp();

  // PDHG_Init_Data(pdhg);

  CUPDLP_CALL(PDTEST_Init_Step_Sizes(pdhg));

  PDTEST_Init_Variables(pdhg);

  // todo: translate check_data into cuda or do it on cpu
  // PDHG_Check_Data(pdhg);

  // PDHG_Print_Header(pdhg);
  // 主要的迭代！！
  cupdlp_int nIter_resetart_value = 0;
  cupdlp_int *nIter_restart = &nIter_resetart_value;
  for (timers->nIter = 0; timers->nIter < settings->nIterLim; ++timers->nIter)
  {
    PDHG_Compute_SolvingTime(pdhg);
#if CUPDLP_DUMP_ITERATES_STATS & CUPDLP_DEBUG
    PDTEST_Dump_Stats(pdhg);
#endif
    int bool_checking = (timers->nIter < 10) ||
                        (timers->nIter == (settings->nIterLim - 1)) ||
                        (timers->dSolvingTime > settings->dTimeLim);
    int bool_print = 0;
#if CUPDLP_DEBUG
    bool_checking = (bool_checking || !(timers->nIter % CUPDLP_DEBUG_INTERVAL));
    bool_print = bool_checking;
#else
    // 检查：迭代次数<10/迭代次数到了/迭代时间到了/迭代次数=0(mod40)
    bool_checking =
        (bool_checking || !(timers->nIter % CUPDLP_RELEASE_INTERVAL));
    // 打印：检查且迭代次数=0(mod特定值)/迭代次数到了/迭代时间到了
    bool_print =
        (bool_checking && !(timers->nIter % (CUPDLP_RELEASE_INTERVAL *
                                             settings->nLogInterval))) ||
        (timers->nIter == (settings->nIterLim - 1)) ||
        (timers->dSolvingTime > settings->dTimeLim);
#endif
    if (bool_checking)
    {
      PDTEST_Compute_Residuals(pdhg);
      if (bool_print)
      {
        PDHG_Print_Header(pdhg);
        PDHG_Print_Iter(pdhg);
        // PDHG_Print_Iter_Average(pdhg);
      }

      if (PDHG_Check_Termination(pdhg, bool_print))
      {
        cupdlp_printf("Optimal current solution.\n");
        resobj->termIterate = LAST_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }

      // if (PDHG_Check_Termination_Average(pdhg, bool_print))
      // {
      //   cupdlp_printf("Optimal average solution.\n");

      //   CUPDLP_COPY_VEC(iterates->x->data, iterates->xAverage->data,
      //                   cupdlp_float, problem->nCols);
      //   CUPDLP_COPY_VEC(iterates->y->data, iterates->yAverage->data,
      //                   cupdlp_float, problem->nRows);

      //   resobj->termIterate = AVERAGE_ITERATE;
      //   resobj->termCode = OPTIMAL;
      //   break;
      // }

      if (timers->dSolvingTime > settings->dTimeLim)
      {
        cupdlp_printf("Time limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }

      if (timers->nIter == (settings->nIterLim - 1))
      {
        cupdlp_printf("Iteration limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }

      PDTEST_Restart_Iterate_Only_Beta(pdhg, nIter_restart); // restart策略
    }
    CUPDLP_CALL(PDTEST_Update_Iterate(pdhg, nIter_restart)); // 迭代更新
    *nIter_restart += 1;
  }
  // print at last
  PDHG_Print_Header(pdhg);
  PDHG_Print_Iter(pdhg);

  // PDHG_Print_Iter_Average(pdhg);

#if PDHG_USE_TIMERS
  cupdlp_printf("Timing information:\n");
  // cupdlp_printf("%20s %e in %d iterations\n", "Total solver time",
  //               timers->dSolvingTime, timers->nIter);
  cupdlp_printf(
      "%20s %e in %d iterations\n", "Total solver time",
      timers->dSolvingTime + timers->dScalingTime + timers->dPresolveTime,
      timers->nIter);
  cupdlp_printf("%20s %e in %d iterations\n", "Solve time",
                timers->dSolvingTime, timers->nIter);
  cupdlp_printf("%20s %e \n", "Iters per sec",
                timers->nIter / timers->dSolvingTime);
  cupdlp_printf("%20s %e\n", "Scaling time", timers->dScalingTime);
  cupdlp_printf("%20s %e\n", "Presolve time", timers->dPresolveTime);
  cupdlp_printf("%20s %e in %d calls\n", "Ax", timers->dAxTime,
                timers->nAxCalls);
  cupdlp_printf("%20s %e in %d calls\n", "Aty", timers->dAtyTime,
                timers->nAtyCalls);
  cupdlp_printf("%20s %e in %d calls\n", "ComputeResiduals",
                timers->dComputeResidualsTime, timers->nComputeResidualsCalls);
  cupdlp_printf("%20s %e in %d calls\n", "UpdateIterates",
                timers->dUpdateIterateTime, timers->nUpdateIterateCalls);
  cupdlp_printf("dMatVecMultiplyTime: %e\n", timers->dMatVecMultiplyTime);
  cupdlp_printf("dVecVecAddTime: %e\n", timers->dVecVecAddTime);
  cupdlp_printf("dIterTime: %e\n", timers->dIterTime);
#endif

#if !(CUPDLP_CPU)
  cupdlp_printf("GPU Timing information:\n");
  cupdlp_printf("%20s %e\n", "CudaPrepare", timers->CudaPrepareTime);
  cupdlp_printf("%20s %e\n", "Alloc&CopyMatToDevice",
                timers->AllocMem_CopyMatToDeviceTime);
  cupdlp_printf("%20s %e\n", "CopyVecToDevice", timers->CopyVecToDeviceTime);
  cupdlp_printf("%20s %e\n", "DeviceMatVecProd", timers->DeviceMatVecProdTime);
  cupdlp_printf("%20s %e\n", "CopyVecToHost", timers->CopyVecToHostTime);
#endif

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDTEST_Solve_best(CUPDLPwork *pdhg)
{
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLPproblem *problem = pdhg->problem;    // 需要求解的问题
  CUPDLPstepsize *stepsize = pdhg->stepsize; // 步长
  CUPDLPsettings *settings = pdhg->settings; // 设置，包括scaling, termination criteria, max iter and time, restart

  CUPDLPresobj *resobj = pdhg->resobj;             /* residuals and objectives 残差和目标函数值*/
  PDTESTiterates *iterates = pdhg->PDTESTiterates; // 迭代变量
  CUPDLPtimers *timers = pdhg->timers;             // 各种耗时

  timers->dSolvingBeg = getTimeStamp();

  // PDHG_Init_Data(pdhg);

  CUPDLP_CALL(PDTEST_Init_Step_Sizes_best(pdhg));

  PDTEST_Init_Variables(pdhg);

  // todo: translate check_data into cuda or do it on cpu
  // PDHG_Check_Data(pdhg);

  // PDHG_Print_Header(pdhg);
  // 主要的迭代！！
  cupdlp_int nIter_resetart_value = 0;
  cupdlp_int *nIter_restart = &nIter_resetart_value;
  for (timers->nIter = 0; timers->nIter < settings->nIterLim; ++timers->nIter)
  {
    PDHG_Compute_SolvingTime(pdhg);
#if CUPDLP_DUMP_ITERATES_STATS & CUPDLP_DEBUG
    PDTEST_Dump_Stats(pdhg);
#endif
    int bool_checking = (timers->nIter < 10) ||
                        (timers->nIter == (settings->nIterLim - 1)) ||
                        (timers->dSolvingTime > settings->dTimeLim);
    int bool_print = 0;
#if CUPDLP_DEBUG
    bool_checking = (bool_checking || !(timers->nIter % CUPDLP_DEBUG_INTERVAL));
    bool_print = bool_checking;
#else
    // 检查：迭代次数<10/迭代次数到了/迭代时间到了/迭代次数=0(mod40)
    bool_checking =
        (bool_checking || !(timers->nIter % CUPDLP_RELEASE_INTERVAL));
    // 打印：检查且迭代次数=0(mod特定值)/迭代次数到了/迭代时间到了
    bool_print =
        (bool_checking && !(timers->nIter % (CUPDLP_RELEASE_INTERVAL *
                                             settings->nLogInterval))) ||
        (timers->nIter == (settings->nIterLim - 1)) ||
        (timers->dSolvingTime > settings->dTimeLim);
#endif
    if (bool_checking)
    {
      PDTEST_Compute_Residuals_best(pdhg);
      if (bool_print)
      {
        PDHG_Print_Header(pdhg);
        PDHG_Print_Iter(pdhg);
        PDHG_Print_Iter_Average(pdhg);
      }
      // PDTEST_Check_Termination_best(pdhg, bool_print);

      if (PDHG_Check_Termination(pdhg, bool_print))
      {
        cupdlp_printf("Optimal current solution.\n");
        resobj->termIterate = LAST_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }

      if (PDHG_Check_Termination_Average(pdhg, bool_print))
      {
        cupdlp_printf("Optimal average solution.\n");

        CUPDLP_COPY_VEC(iterates->x->data, iterates->xAverage->data,
                        cupdlp_float, problem->nCols);
        CUPDLP_COPY_VEC(iterates->y->data, iterates->yAverage->data,
                        cupdlp_float, problem->nRows);

        resobj->termIterate = AVERAGE_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }

      if (timers->dSolvingTime > settings->dTimeLim)
      {
        cupdlp_printf("Time limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }

      if (timers->nIter == (settings->nIterLim - 1))
      {
        cupdlp_printf("Iteration limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }

      PDTEST_Restart_Iterate_best(pdhg, nIter_restart); // restart策略
    }
    CUPDLP_CALL(PDTEST_Update_Iterate_best(pdhg, nIter_restart)); // 迭代更新
    *nIter_restart += 1;
  }
  // print at last
  PDHG_Print_Header(pdhg);
  PDHG_Print_Iter(pdhg);

  // PDHG_Print_Iter_Average(pdhg);

#if PDHG_USE_TIMERS
  cupdlp_printf("Timing information:\n");
  // cupdlp_printf("%20s %e in %d iterations\n", "Total solver time",
  //               timers->dSolvingTime, timers->nIter);
  cupdlp_printf(
      "%20s %e in %d iterations\n", "Total solver time",
      timers->dSolvingTime + timers->dScalingTime + timers->dPresolveTime,
      timers->nIter);
  cupdlp_printf("%20s %e in %d iterations\n", "Solve time",
                timers->dSolvingTime, timers->nIter);
  cupdlp_printf("%20s %e \n", "Iters per sec",
                timers->nIter / timers->dSolvingTime);
  cupdlp_printf("%20s %e\n", "Scaling time", timers->dScalingTime);
  cupdlp_printf("%20s %e\n", "Presolve time", timers->dPresolveTime);
  cupdlp_printf("%20s %e in %d calls\n", "Ax", timers->dAxTime,
                timers->nAxCalls);
  cupdlp_printf("%20s %e in %d calls\n", "Aty", timers->dAtyTime,
                timers->nAtyCalls);
  cupdlp_printf("%20s %e in %d calls\n", "ComputeResiduals",
                timers->dComputeResidualsTime, timers->nComputeResidualsCalls);
  cupdlp_printf("%20s %e in %d calls\n", "UpdateIterates",
                timers->dUpdateIterateTime, timers->nUpdateIterateCalls);
  cupdlp_printf("dMatVecMultiplyTime: %e\n", timers->dMatVecMultiplyTime);
  cupdlp_printf("dVecVecAddTime: %e\n", timers->dVecVecAddTime);
  cupdlp_printf("dIterTime: %e\n", timers->dIterTime);
#endif

#if !(CUPDLP_CPU)
  cupdlp_printf("GPU Timing information:\n");
  cupdlp_printf("%20s %e\n", "CudaPrepare", timers->CudaPrepareTime);
  cupdlp_printf("%20s %e\n", "Alloc&CopyMatToDevice",
                timers->AllocMem_CopyMatToDeviceTime);
  cupdlp_printf("%20s %e\n", "CopyVecToDevice", timers->CopyVecToDeviceTime);
  cupdlp_printf("%20s %e\n", "DeviceMatVecProd", timers->DeviceMatVecProdTime);
  cupdlp_printf("%20s %e\n", "CopyVecToHost", timers->CopyVecToHostTime);
#endif

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDTEST_Solve_best2(CUPDLPwork *pdhg)
{
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLPproblem *problem = pdhg->problem;    // 需要求解的问题
  CUPDLPstepsize *stepsize = pdhg->stepsize; // 步长
  CUPDLPsettings *settings = pdhg->settings; // 设置，包括scaling, termination criteria, max iter and time, restart

  CUPDLPresobj *resobj = pdhg->resobj;             /* residuals and objectives 残差和目标函数值*/
  PDTESTiterates *iterates = pdhg->PDTESTiterates; // 迭代变量
  CUPDLPtimers *timers = pdhg->timers;             // 各种耗时

  timers->dSolvingBeg = getTimeStamp();

  // PDHG_Init_Data(pdhg);

  CUPDLP_CALL(PDTEST_Init_Step_Sizes(pdhg));

  PDTEST_Init_Variables(pdhg);

  // todo: translate check_data into cuda or do it on cpu
  // PDHG_Check_Data(pdhg);

  // PDHG_Print_Header(pdhg);
  // 主要的迭代！！
  cupdlp_int nIter_resetart_value = 0;
  cupdlp_int *nIter_restart = &nIter_resetart_value;
  for (timers->nIter = 0; timers->nIter < settings->nIterLim; ++timers->nIter)
  {
    PDHG_Compute_SolvingTime(pdhg);
#if CUPDLP_DUMP_ITERATES_STATS & CUPDLP_DEBUG
    PDTEST_Dump_Stats(pdhg);
#endif
    int bool_checking = (timers->nIter < 10) ||
                        (timers->nIter == (settings->nIterLim - 1)) ||
                        (timers->dSolvingTime > settings->dTimeLim);
    int bool_print = 0;
#if CUPDLP_DEBUG
    bool_checking = (bool_checking || !(timers->nIter % CUPDLP_DEBUG_INTERVAL));
    bool_print = bool_checking;
#else
    // 检查：迭代次数<10/迭代次数到了/迭代时间到了/迭代次数=0(mod40)
    bool_checking =
        (bool_checking || !(timers->nIter % CUPDLP_RELEASE_INTERVAL));
    // 打印：检查且迭代次数=0(mod特定值)/迭代次数到了/迭代时间到了
    bool_print =
        (bool_checking && !(timers->nIter % (CUPDLP_RELEASE_INTERVAL *
                                             settings->nLogInterval))) ||
        (timers->nIter == (settings->nIterLim - 1)) ||
        (timers->dSolvingTime > settings->dTimeLim);
#endif
    if (bool_checking)
    {
      PDTEST_Compute_Residuals(pdhg);
      if (bool_print)
      {
        PDHG_Print_Header(pdhg);
        PDHG_Print_Iter(pdhg);
        // PDHG_Print_Iter_Average(pdhg);
      }

      if (PDHG_Check_Termination(pdhg, bool_print))
      {
        cupdlp_printf("Optimal current solution.\n");
        resobj->termIterate = LAST_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }

      // if (PDHG_Check_Termination_Average(pdhg, bool_print))
      // {
      //   cupdlp_printf("Optimal average solution.\n");

      //   CUPDLP_COPY_VEC(iterates->x->data, iterates->xAverage->data,
      //                   cupdlp_float, problem->nCols);
      //   CUPDLP_COPY_VEC(iterates->y->data, iterates->yAverage->data,
      //                   cupdlp_float, problem->nRows);

      //   resobj->termIterate = AVERAGE_ITERATE;
      //   resobj->termCode = OPTIMAL;
      //   break;
      // }

      if (timers->dSolvingTime > settings->dTimeLim)
      {
        cupdlp_printf("Time limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }

      if (timers->nIter == (settings->nIterLim - 1))
      {
        cupdlp_printf("Iteration limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }

      PDTEST_Restart_Iterate_best2(pdhg, nIter_restart); // restart策略
    }
    CUPDLP_CALL(PDTEST_Update_Iterate_best2(pdhg, nIter_restart)); // 迭代更新
    *nIter_restart += 1;
  }
  // print at last
  PDHG_Print_Header(pdhg);
  PDHG_Print_Iter(pdhg);

  // PDHG_Print_Iter_Average(pdhg);

#if PDHG_USE_TIMERS
  cupdlp_printf("Timing information:\n");
  // cupdlp_printf("%20s %e in %d iterations\n", "Total solver time",
  //               timers->dSolvingTime, timers->nIter);
  cupdlp_printf(
      "%20s %e in %d iterations\n", "Total solver time",
      timers->dSolvingTime + timers->dScalingTime + timers->dPresolveTime,
      timers->nIter);
  cupdlp_printf("%20s %e in %d iterations\n", "Solve time",
                timers->dSolvingTime, timers->nIter);
  cupdlp_printf("%20s %e \n", "Iters per sec",
                timers->nIter / timers->dSolvingTime);
  cupdlp_printf("%20s %e\n", "Scaling time", timers->dScalingTime);
  cupdlp_printf("%20s %e\n", "Presolve time", timers->dPresolveTime);
  cupdlp_printf("%20s %e in %d calls\n", "Ax", timers->dAxTime,
                timers->nAxCalls);
  cupdlp_printf("%20s %e in %d calls\n", "Aty", timers->dAtyTime,
                timers->nAtyCalls);
  cupdlp_printf("%20s %e in %d calls\n", "ComputeResiduals",
                timers->dComputeResidualsTime, timers->nComputeResidualsCalls);
  cupdlp_printf("%20s %e in %d calls\n", "UpdateIterates",
                timers->dUpdateIterateTime, timers->nUpdateIterateCalls);
  cupdlp_printf("dMatVecMultiplyTime: %e\n", timers->dMatVecMultiplyTime);
  cupdlp_printf("dVecVecAddTime: %e\n", timers->dVecVecAddTime);
  cupdlp_printf("dIterTime: %e\n", timers->dIterTime);
#endif

#if !(CUPDLP_CPU)
  cupdlp_printf("GPU Timing information:\n");
  cupdlp_printf("%20s %e\n", "CudaPrepare", timers->CudaPrepareTime);
  cupdlp_printf("%20s %e\n", "Alloc&CopyMatToDevice",
                timers->AllocMem_CopyMatToDeviceTime);
  cupdlp_printf("%20s %e\n", "CopyVecToDevice", timers->CopyVecToDeviceTime);
  cupdlp_printf("%20s %e\n", "DeviceMatVecProd", timers->DeviceMatVecProdTime);
  cupdlp_printf("%20s %e\n", "CopyVecToHost", timers->CopyVecToHostTime);
#endif

exit_cleanup:
  return retcode;
}

void PDHG_PostSolve(CUPDLPwork *pdhg, cupdlp_int nCols_origin,
                    cupdlp_int *constraint_new_idx, cupdlp_float *x_origin,
                    cupdlp_float *y_origin)
{
  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPscaling *scaling = pdhg->scaling;

  // unscale
  if (scaling->ifScaled)
  {
    cupdlp_ediv(iterates->x->data, pdhg->colScale, problem->nCols);
    cupdlp_ediv(iterates->y->data, pdhg->rowScale, problem->nRows);

    // cupdlp_ediv(iterates->x->data, scaling->colScale_gpu, problem->nCols);
    // cupdlp_ediv(iterates->y->data, scaling->rowScale_gpu, problem->nRows);
  }

  // extract x from (x, z)
  CUPDLP_COPY_VEC(x_origin, iterates->x->data, cupdlp_float, nCols_origin);

  cupdlp_float *ytmp =
      (cupdlp_float *)cupdlp_malloc(problem->nRows * sizeof(cupdlp_float));
  CUPDLP_COPY_VEC(ytmp, iterates->y->data, cupdlp_float, problem->nRows);
  // un-permute y
  for (int i = 0; i < problem->nRows; i++)
  {
    y_origin[i] = ytmp[constraint_new_idx[i]];
  }
  cupdlp_free(ytmp);
}
void PDTEST_PostSolve(CUPDLPwork *pdhg, cupdlp_int nCols_origin,
                      cupdlp_int *constraint_new_idx, cupdlp_float *x_origin,
                      cupdlp_float *y_origin)
{
  CUPDLPproblem *problem = pdhg->problem;
  PDTESTiterates *iterates = pdhg->PDTESTiterates;
  CUPDLPscaling *scaling = pdhg->scaling;

  // unscale
  if (scaling->ifScaled)
  {
    cupdlp_ediv(iterates->x_ag->data, pdhg->colScale, problem->nCols);
    cupdlp_ediv(iterates->y_ag->data, pdhg->rowScale, problem->nRows);

    // cupdlp_ediv(iterates->x->data, scaling->colScale_gpu, problem->nCols);
    // cupdlp_ediv(iterates->y->data, scaling->rowScale_gpu, problem->nRows);
  }

  // extract x from (x, z)
  CUPDLP_COPY_VEC(x_origin, iterates->x_ag->data, cupdlp_float, nCols_origin);

  cupdlp_float *ytmp =
      (cupdlp_float *)cupdlp_malloc(problem->nRows * sizeof(cupdlp_float));
  CUPDLP_COPY_VEC(ytmp, iterates->y_ag->data, cupdlp_float, problem->nRows);
  // un-permute y
  for (int i = 0; i < problem->nRows; i++)
  {
    y_origin[i] = ytmp[constraint_new_idx[i]];
  }
  cupdlp_free(ytmp);
}
cupdlp_retcode LP_SolvePDHG_Multiscale_withStepsize(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_bool *ifChangeFloatParam, cupdlp_float *floatParam, char *fp, cupdlp_bool ifSaveSol, cupdlp_int *constraint_new_idx, cupdlp_float **x_solution, cupdlp_float **y_solution, cupdlp_float *x_init, cupdlp_float *y_init, cupdlp_float *stepsize_last, cupdlp_float *weight_last, cupdlp_float *stepsize_init, cupdlp_float *weight_init, cupdlp_bool whether_first)
{
  cupdlp_retcode retcode = RETCODE_OK;

  PDHG_PrintHugeCUPDHG();

  CUPDLP_CALL(PDHG_SetUserParam(pdhg, ifChangeIntParam, intParam,
                                ifChangeFloatParam, floatParam));

  CUPDLP_CALL(PDHG_Solve_Multiscale_withStepsize(pdhg, x_init, y_init, stepsize_last, weight_last, stepsize_init, weight_init, whether_first));
  cupdlp_printf("fout: %s\n", fp);

  // cupdlp_int nCols_origin = 2 * pow(resolution, 2);
  // cupdlp_int nRows = pow(resolution, 4);
  cupdlp_int nCols_origin = pdhg->problem->nCols;
  cupdlp_int nRows = pdhg->problem->nRows;
  cupdlp_float *x_origin = cupdlp_NULL;
  cupdlp_float *y_origin = cupdlp_NULL;
  CUPDLP_INIT_ZERO(x_origin, nCols_origin);
  CUPDLP_INIT_ZERO(y_origin, nRows);

  PDHG_PostSolve(pdhg, nCols_origin, constraint_new_idx, x_origin, y_origin);

  *x_solution = x_origin;
  *y_solution = y_origin;
  cupdlp_printf("Begin writing!\n");
  writeJson(fp, pdhg, x_origin, nCols_origin, y_origin, pdhg->problem->nRows, ifSaveSol);

exit_cleanup:
  if (retcode != RETCODE_OK)
  {
    cupdlp_printf("Error in LP_SolvePDHG\n");
  }
  // PDHG_Destroy(&pdhg);
  return retcode;
}

cupdlp_retcode LP_SolvePDHG_Multiscale(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_bool *ifChangeFloatParam, cupdlp_float *floatParam, char *fp, cupdlp_bool ifSaveSol, cupdlp_int *constraint_new_idx, cupdlp_float **x_solution, cupdlp_float **y_solution, cupdlp_float *x_init, cupdlp_float *y_init)
{
  cupdlp_retcode retcode = RETCODE_OK;

  PDHG_PrintHugeCUPDHG();

  CUPDLP_CALL(PDHG_SetUserParam(pdhg, ifChangeIntParam, intParam,
                                ifChangeFloatParam, floatParam));

  CUPDLP_CALL(PDHG_Solve_Multiscale(pdhg, x_init, y_init));
  cupdlp_printf("fout: %s\n", fp);

  // cupdlp_int nCols_origin = 2 * pow(resolution, 2);
  // cupdlp_int nRows = pow(resolution, 4);
  cupdlp_int nCols_origin = pdhg->problem->nCols;
  cupdlp_int nRows = pdhg->problem->nRows;
  cupdlp_float *x_origin = cupdlp_NULL;
  cupdlp_float *y_origin = cupdlp_NULL;
  CUPDLP_INIT_ZERO(x_origin, nCols_origin);
  CUPDLP_INIT_ZERO(y_origin, nRows);

  PDHG_PostSolve(pdhg, nCols_origin, constraint_new_idx, x_origin, y_origin);

  *x_solution = x_origin;
  *y_solution = y_origin;
  cupdlp_printf("Begin writing!\n");
  cupdlp_printf("约束矩阵尺寸为nRows：%d, nCols：%d\n", pdhg->problem->nRows, pdhg->problem->nCols);
  writeJson(fp, pdhg, x_origin, nCols_origin, y_origin, pdhg->problem->nRows, ifSaveSol);

exit_cleanup:
  if (retcode != RETCODE_OK)
  {
    cupdlp_printf("Error in LP_SolvePDHG\n");
  }
  // PDHG_Destroy(&pdhg);
  return retcode;
}

cupdlp_retcode LP_SolvePDHG(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_bool *ifChangeFloatParam, cupdlp_float *floatParam, char *fp, cupdlp_bool ifSaveSol, cupdlp_int *constraint_new_idx, cupdlp_int resolution, cupdlp_float **x_solution, cupdlp_float **y_solution)
{
  cupdlp_retcode retcode = RETCODE_OK;

  PDHG_PrintHugeCUPDHG();

  CUPDLP_CALL(PDHG_SetUserParam(pdhg, ifChangeIntParam, intParam,
                                ifChangeFloatParam, floatParam));

  CUPDLP_CALL(PDHG_Solve(pdhg));
  cupdlp_printf("fout: %s\n", fp);

  cupdlp_int nCols_origin = 2 * pow(resolution, 2);
  cupdlp_int nRows = pow(resolution, 4);
  cupdlp_float *x_origin = cupdlp_NULL;
  cupdlp_float *y_origin = cupdlp_NULL;
  CUPDLP_INIT_ZERO(x_origin, nCols_origin);
  CUPDLP_INIT_ZERO(y_origin, nRows);

  PDHG_PostSolve(pdhg, nCols_origin, constraint_new_idx, x_origin, y_origin);

  *x_solution = x_origin;
  *y_solution = y_origin;
  // x_solution = &x_origin;
  // y_solution = &y_origin;
  cupdlp_printf("Begin writing!\n");
  writeJson(fp, pdhg, x_origin, nCols_origin, y_origin, pdhg->problem->nRows, ifSaveSol);

exit_cleanup:
  if (retcode != RETCODE_OK)
  {
    cupdlp_printf("Error in LP_SolvePDHG\n");
  }
  PDHG_Destroy(&pdhg);
  return retcode;
}

cupdlp_retcode LP_SolvePDHG_AdapTheta(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_bool *ifChangeFloatParam, cupdlp_float *floatParam, char *fp, cupdlp_float *x_origin, cupdlp_int nCols_origin, cupdlp_float *y_origin, cupdlp_bool ifSaveSol, cupdlp_int *constraint_new_idx)
{
  cupdlp_retcode retcode = RETCODE_OK;

  PDHG_PrintHugeCUPDHG();

  CUPDLP_CALL(PDHG_SetUserParam(pdhg, ifChangeIntParam, intParam,
                                ifChangeFloatParam, floatParam));

  CUPDLP_CALL(PDHG_Solve_AdapTheta(pdhg));

  PDHG_PostSolve(pdhg, nCols_origin, constraint_new_idx, x_origin, y_origin);

  writeJson(fp, pdhg, x_origin, nCols_origin, y_origin, pdhg->problem->nRows,
            ifSaveSol);

exit_cleanup:
  PDHG_Destroy(&pdhg);
  return retcode;
}

cupdlp_retcode LP_SolvePDTEST(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_bool *ifChangeFloatParam, cupdlp_float *floatParam, char *fp, cupdlp_float *x_origin, cupdlp_int nCols_origin, cupdlp_float *y_origin, cupdlp_bool ifSaveSol, cupdlp_int *constraint_new_idx)
{
  cupdlp_retcode retcode = RETCODE_OK;

  PDHG_PrintHugeCUPDHG();

  CUPDLP_CALL(PDTEST_SetUserParam(pdhg, ifChangeIntParam, intParam,
                                  ifChangeFloatParam, floatParam));

  CUPDLP_CALL(PDTEST_Solve(pdhg));

  // 后处理，应该是把原来的scaling之类的变回去
  PDTEST_PostSolve(pdhg, nCols_origin, constraint_new_idx, x_origin, y_origin);

  writeJson(fp, pdhg, x_origin, nCols_origin, y_origin, pdhg->problem->nRows,
            ifSaveSol);

exit_cleanup:
  PDHG_Destroy(&pdhg);
  return retcode;
}

cupdlp_retcode LP_SolvePDTEST_best(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_bool *ifChangeFloatParam, cupdlp_float *floatParam, char *fp, cupdlp_float *x_origin, cupdlp_int nCols_origin, cupdlp_float *y_origin, cupdlp_bool ifSaveSol, cupdlp_int *constraint_new_idx)
{
  cupdlp_retcode retcode = RETCODE_OK;

  PDHG_PrintHugeCUPDHG();

  CUPDLP_CALL(PDTEST_SetUserParam(pdhg, ifChangeIntParam, intParam,
                                  ifChangeFloatParam, floatParam));

  CUPDLP_CALL(PDTEST_Solve_best(pdhg));

  // 后处理，应该是把原来的scaling之类的变回去
  PDTEST_PostSolve(pdhg, nCols_origin, constraint_new_idx, x_origin, y_origin);

  writeJson(fp, pdhg, x_origin, nCols_origin, y_origin, pdhg->problem->nRows,
            ifSaveSol);

exit_cleanup:
  PDHG_Destroy(&pdhg);
  return retcode;
}
cupdlp_retcode LP_SolvePDTEST_best2(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_bool *ifChangeFloatParam, cupdlp_float *floatParam, char *fp, cupdlp_float *x_origin, cupdlp_int nCols_origin, cupdlp_float *y_origin, cupdlp_bool ifSaveSol, cupdlp_int *constraint_new_idx)
{
  cupdlp_retcode retcode = RETCODE_OK;

  PDHG_PrintHugeCUPDHG();

  CUPDLP_CALL(PDTEST_SetUserParam(pdhg, ifChangeIntParam, intParam,
                                  ifChangeFloatParam, floatParam));

  CUPDLP_CALL(PDTEST_Solve_best2(pdhg));

  // 后处理，应该是把原来的scaling之类的变回去
  PDTEST_PostSolve(pdhg, nCols_origin, constraint_new_idx, x_origin, y_origin);

  writeJson(fp, pdhg, x_origin, nCols_origin, y_origin, pdhg->problem->nRows,
            ifSaveSol);

exit_cleanup:
  PDHG_Destroy(&pdhg);
  return retcode;
}

cupdlp_retcode LP_SolvePDTEST_Average(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_bool *ifChangeFloatParam, cupdlp_float *floatParam, char *fp, cupdlp_float *x_origin, cupdlp_int nCols_origin, cupdlp_float *y_origin, cupdlp_bool ifSaveSol, cupdlp_int *constraint_new_idx)
{
  cupdlp_retcode retcode = RETCODE_OK;

  PDHG_PrintHugeCUPDHG();

  CUPDLP_CALL(PDTEST_SetUserParam(pdhg, ifChangeIntParam, intParam,
                                  ifChangeFloatParam, floatParam));

  CUPDLP_CALL(PDTEST_Average_Solve(pdhg));

  // 后处理，应该是把原来的scaling之类的变回去
  PDTEST_PostSolve(pdhg, nCols_origin, constraint_new_idx, x_origin, y_origin);

  writeJson(fp, pdhg, x_origin, nCols_origin, y_origin, pdhg->problem->nRows,
            ifSaveSol);

exit_cleanup:
  PDHG_Destroy(&pdhg);
  return retcode;
}

cupdlp_retcode LP_SolvePDTEST_min(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_bool *ifChangeFloatParam, cupdlp_float *floatParam, char *fp, cupdlp_float *x_origin, cupdlp_int nCols_origin, cupdlp_float *y_origin, cupdlp_bool ifSaveSol, cupdlp_int *constraint_new_idx)
{
  cupdlp_retcode retcode = RETCODE_OK;

  PDHG_PrintHugeCUPDHG();

  CUPDLP_CALL(PDTEST_SetUserParam(pdhg, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam));
  PDHG_Init_Data(pdhg);
  CUPDLP_CALL(PDTEST_Init_Step_Sizes(pdhg));
  PDTEST_Init_Variables(pdhg);
  cupdlp_printf("--------------------------------------------------\n");
  // PDTEST_printCudaMatGPU(pdhg);

  cupdlp_printf("--------------------------------------------------\n");
  CUPDLPtimers *timers = pdhg->timers;             // 各种耗时
  CUPDLPsettings *settings = pdhg->settings;       // 设置，包括scaling, termination criteria, max iter and time, restart
  PDTESTiterates *iterates = pdhg->PDTESTiterates; // 迭代变量
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPproblem *problem = pdhg->problem;
  cupdlp_int restartFreq = 1000;
  CUPDLPresobj *resobj = pdhg->resobj;
  cupdlp_float stopThr = 1.0E-5;
  cupdlp_printf("stopThr = %f\n", stopThr);
  timers->dSolvingTime = getTimeStamp();
  for (timers->nIter = 0; timers->nIter < settings->nIterLim; ++timers->nIter)
  {
    cupdlp_int t = timers->nIter + 1;
    cupdlp_float beta = (t + 3.8) / 4.8;
    cupdlp_float theta = 1.0;

    // 计算平均解
    // PDHG_Compute_Average_Iterate(pdhg);

#pragma region 迭代
    // PDTEST_printCudaDenseVecGPU(iterates->y_ag);
    // PDTEST_printCudaDenseVecGPU(iterates->x_ag);
    // 计算Ax_bar^{t}, 后面计算y^{t+1}有用
    // Ax(pdhg, iterates->ax_bar, iterates->x_bar);
    // // y^{t+1} = y^t + dDualStep * (b - A * (x_bar^{t})
    // PDTEST_dualGradientStep(pdhg, stepsize->dDualStep);
    // cupdlp_printf("y^{t+1} = y^t + dDualStep * (b - A * (x_bar^{t})\n");
    // PDHG_Project_Row_Duals(pdhg, iterates->yUpdate->data);
    // // PDTEST_printCudaDenseVecGPU(iterates->yUpdate);
    // // 计算ATy^{t+1}, 后面计算x^{t+1}有用
    // ATy(pdhg, iterates->atyUpdate, iterates->yUpdate);
    // // x^{t+1} = proj_{X}(x^t - dPrimalStep * (c - A'y^{t+1}))
    // PDTEST_primalGradientStep(pdhg, stepsize->dPrimalStep);
    // cupdlp_printf("x^{t+1} = proj_{X}(x^t - dPrimalStep * (c - A'y^{t+1}))\n");
    // PDHG_Project_Bounds(pdhg, iterates->xUpdate->data);
    // // PDTEST_printCudaDenseVecGPU(iterates->xUpdate);
    // // x_ag^{t+1} = (1 - 1 / beta^t)x_ag^{t} + (1 / beta^t)x^{t+1}
    // cupdlp_printf("beta = %f\n", beta);
    // PDTEST_x_ag_step(pdhg, beta);
    // cupdlp_printf("x_ag^{t+1} = (1 - 1 / beta^t)x_ag^{t} + (1 / beta^t)x^{t+1}\n");
    // // PDTEST_printCudaDenseVecGPU(iterates->x_agUpdate);
    // // PDTEST_printCudaDenseVecGPU(iterates->x_ag);
    // // PDTEST_printCudaDenseVecGPU(iterates->xUpdate);
    // // 没有必要进行Project，因为都是已经投影过的x进行线性组合
    // // PDHG_Project_Bounds(pdhg, iterates->x_agUpdate->data);
    // // y_ag^{t+1} = (1 - 1 / beta^t)y_ag^{t} + (1 / beta^t)y^{t}
    // PDTEST_y_ag_step(pdhg, beta);
    // // x_bar^{t+1} = theta^{t+1}(x^{t+1} - x^{t}) + x^{t+1}
    // PDTEST_x_bar_step(pdhg, theta);
    // // 更新一下ax_ag和aty_ag
    // Ax(pdhg, iterates->ax_ag, iterates->x_agUpdate);
    // ATy(pdhg, iterates->aty_ag, iterates->y_agUpdate);

    Ax(pdhg, iterates->ax_bar, iterates->x_bar);
    PDTEST_dualGradientStep(pdhg, stepsize->dDualStep);
    PDHG_Project_Row_Duals(pdhg, iterates->yUpdate->data);
    ATy(pdhg, iterates->atyUpdate, iterates->yUpdate);
    PDTEST_primalGradientStep(pdhg, stepsize->dPrimalStep);
    PDHG_Project_Bounds(pdhg, iterates->xUpdate->data);
    PDTEST_x_ag_step(pdhg, beta);
    PDTEST_y_ag_step(pdhg, beta);
    PDTEST_x_bar_step(pdhg, theta);
    Ax(pdhg, iterates->ax_agUpdate, iterates->x_agUpdate);
    ATy(pdhg, iterates->aty_agUpdate, iterates->y_agUpdate);

#pragma endregion

#pragma region 检查是否收敛
    int bool_print = 1;
    if (timers->nIter % 40 == 39)
    {
      PDTEST_Compute_Residuals(pdhg);
      PDHG_Print_Header(pdhg);
      PDHG_Print_Iter(pdhg);
      if (fabs(resobj->dDualityGap) < stopThr)
      {
        cupdlp_printf("Gap = %f\n", resobj->dDualityGap);

        cupdlp_printf("Optimal current solution.\n");
        resobj->termIterate = LAST_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }
      //   if (PDHG_Check_Termination_Average(pdhg, bool_print))
      //   {
      //     cupdlp_printf("Optimal average solution.\n");

      //     CUPDLP_COPY_VEC(iterates->x_ag->data, iterates->x_agAverage->data,
      //                     cupdlp_float, problem->nCols);
      //     CUPDLP_COPY_VEC(iterates->y_ag->data, iterates->y_agAverage->data,
      //                     cupdlp_float, problem->nRows);

      //     resobj->termIterate = AVERAGE_ITERATE;
      //     resobj->termCode = OPTIMAL;
      //     break;
      //   }
    }
#pragma endregion

#pragma region 更新
    resobj->dPrimalFeasLastRestart = resobj->dPrimalFeas;
    resobj->dDualFeasLastRestart = resobj->dDualFeas;
    resobj->dDualityGapLastRestart = resobj->dDualityGap;
    CUPDLP_COPY_VEC(iterates->x->data, iterates->xUpdate->data, cupdlp_float, problem->nCols);
    CUPDLP_COPY_VEC(iterates->x_bar->data, iterates->x_barUpdate->data, cupdlp_float, problem->nCols);
    CUPDLP_COPY_VEC(iterates->x_ag->data, iterates->x_agUpdate->data, cupdlp_float, problem->nCols);
    CUPDLP_COPY_VEC(iterates->y->data, iterates->yUpdate->data, cupdlp_float, problem->nRows);
    CUPDLP_COPY_VEC(iterates->y_ag->data, iterates->y_agUpdate->data, cupdlp_float, problem->nRows);
#pragma endregion

#pragma region 重启
    if (timers->nIter % restartFreq == restartFreq - 1)
    {
      // 重启
      cupdlp_printf("Restart to current\n");
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
#pragma endregion
  }
  timers->dSolvingTime = getTimeStamp() - timers->dSolvingTime;
  PDTEST_Compute_Residuals(pdhg);
  PDHG_Print_Header(pdhg);
  PDHG_Print_Iter(pdhg);
  // 后处理，应该是把原来的scaling之类的变回去
  PDTEST_PostSolve(pdhg, nCols_origin, constraint_new_idx, x_origin, y_origin);

  writeJson(fp, pdhg, x_origin, nCols_origin, y_origin, pdhg->problem->nRows,
            ifSaveSol);

exit_cleanup:
  PDHG_Destroy(&pdhg);
  return retcode;
}