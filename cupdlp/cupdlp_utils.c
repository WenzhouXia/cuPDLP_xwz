//
// Created by chuwen on 23-11-26.
//

#include "cupdlp_utils.h"

#include <stdio.h>

#include "cupdlp_cs.h"
#include "cupdlp_linalg.h"
#include "glbopts.h"

#if !(CUPDLP_CPU)

#include "cuda/cupdlp_cudalinalg.cuh"

#endif

void dense_clear(CUPDLPdense *dense)
{
  if (dense)
  {
    if (dense->data)
    {
      cupdlp_free(dense->data);
    }
    cupdlp_free(dense);
  }
}

cupdlp_int csc_clear(CUPDLPcsc *csc)
{
  if (csc)
  {
#if !(CUPDLP_CPU)
    if (csc->cuda_csc != NULL)
    {
      CHECK_CUSPARSE(cusparseDestroySpMat(csc->cuda_csc))
    }
#endif
    if (csc->colMatBeg)
    {
      CUPDLP_FREE_VEC(csc->colMatBeg);
    }
    if (csc->colMatIdx)
    {
      CUPDLP_FREE_VEC(csc->colMatIdx);
    }
    if (csc->colMatElem)
    {
      CUPDLP_FREE_VEC(csc->colMatElem);
    }
    CUPDLP_FREE_VEC(csc);
  }
  return 0;
}

cupdlp_int csr_clear(CUPDLPcsr *csr)
{
  if (csr)
  {
#if !(CUPDLP_CPU)
    if (csr->cuda_csr != NULL)
    {
      CHECK_CUSPARSE(cusparseDestroySpMat(csr->cuda_csr))
    }
#endif
    if (csr->rowMatBeg)
    {
      CUPDLP_FREE_VEC(csr->rowMatBeg);
    }
    if (csr->rowMatIdx)
    {
      CUPDLP_FREE_VEC(csr->rowMatIdx);
    }
    if (csr->rowMatElem)
    {
      CUPDLP_FREE_VEC(csr->rowMatElem);
    }
    CUPDLP_FREE_VEC(csr);
  }
  return 0;
}

void data_clear(CUPDLPdata *data)
{
  if (data)
  {
    switch (data->matrix_format)
    {
    case DENSE:
      dense_clear(data->dense_matrix);
      break;
    case CSR:
      csr_clear(data->csr_matrix);
      break;
    case CSC:
      csc_clear(data->csc_matrix);
      break;
    case CSR_CSC:
      csr_clear(data->csr_matrix);
      csc_clear(data->csc_matrix);
      break;
    }
    cupdlp_free(data);
  }
}

void problem_clear(CUPDLPproblem *problem)
{
  if (problem)
  {
    if (problem->data)
    {
      data_clear(problem->data);
    }
    //        if (problem->colMatBeg) {
    //            cupdlp_free(problem->colMatBeg);
    //        }
    //        if (problem->colMatIdx) {
    //            cupdlp_free(problem->colMatIdx);
    //        }
    //        if (problem->colMatElem) {
    //            cupdlp_free(problem->colMatElem);
    //        }
    //        if (problem->rowMatBeg) {
    //            cupdlp_free(problem->rowMatBeg);
    //        }
    //        if (problem->rowMatIdx) {
    //            cupdlp_free(problem->rowMatIdx);
    //        }
    //        if (problem->rowMatElem) {
    //            cupdlp_free(problem->rowMatElem);
    //        }
    if (problem->lower)
    {
      CUPDLP_FREE_VEC(problem->lower);
    }
    if (problem->upper)
    {
      CUPDLP_FREE_VEC(problem->upper);
    }
    if (problem->cost)
    {
      CUPDLP_FREE_VEC(problem->cost);
    }
    if (problem->rhs)
    {
      CUPDLP_FREE_VEC(problem->rhs);
    }
    if (problem->hasLower)
    {
      CUPDLP_FREE_VEC(problem->hasLower);
    }
    if (problem->hasUpper)
    {
      CUPDLP_FREE_VEC(problem->hasUpper);
    }
    CUPDLP_FREE_VEC(problem);
  }
}

void settings_clear(CUPDLPsettings *settings)
{
  if (settings)
  {
    cupdlp_free(settings);
  }
}

cupdlp_int vec_clear(CUPDLPvec *vec)
{
  if (vec)
  {
    if (vec->data)
    {
      CUPDLP_FREE_VEC(vec->data);
    }
#if !(CUPDLP_CPU)
    CHECK_CUSPARSE(cusparseDestroyDnVec(vec->cuda_vec))
#endif
    cupdlp_free(vec);
  }

  return 0;
}

void iterates_clear(CUPDLPiterates *iterates)
{
  if (iterates)
  {
    if (iterates->x)
    {
      // CUPDLP_FREE_VEC(iterates->x);
      vec_clear(iterates->x);
    }
    if (iterates->y)
    {
      // CUPDLP_FREE_VEC(iterates->y);
      vec_clear(iterates->y);
    }
    if (iterates->xUpdate)
    {
      // CUPDLP_FREE_VEC(iterates->xUpdate);
      vec_clear(iterates->xUpdate);
    }
    if (iterates->yUpdate)
    {
      // CUPDLP_FREE_VEC(iterates->yUpdate);
      vec_clear(iterates->yUpdate);
    }
    if (iterates->xSum)
    {
      CUPDLP_FREE_VEC(iterates->xSum);
    }
    if (iterates->ySum)
    {
      CUPDLP_FREE_VEC(iterates->ySum);
    }
    if (iterates->xAverage)
    {
      // CUPDLP_FREE_VEC(iterates->xAverage);
      vec_clear(iterates->xAverage);
    }
    if (iterates->yAverage)
    {
      // CUPDLP_FREE_VEC(iterates->yAverage);
      vec_clear(iterates->yAverage);
    }
    if (iterates->xLastRestart)
    {
      CUPDLP_FREE_VEC(iterates->xLastRestart);
    }
    if (iterates->yLastRestart)
    {
      CUPDLP_FREE_VEC(iterates->yLastRestart);
    }
    if (iterates->ax)
    {
      // CUPDLP_FREE_VEC(iterates->ax);
      vec_clear(iterates->ax);
    }
    if (iterates->axUpdate)
    {
      // CUPDLP_FREE_VEC(iterates->axUpdate);
      vec_clear(iterates->axUpdate);
    }
    if (iterates->axAverage)
    {
      // CUPDLP_FREE_VEC(iterates->axAverage);
      vec_clear(iterates->axAverage);
    }
    if (iterates->aty)
    {
      // CUPDLP_FREE_VEC(iterates->aty);
      vec_clear(iterates->aty);
    }
    if (iterates->atyUpdate)
    {
      // CUPDLP_FREE_VEC(iterates->atyUpdate);
      vec_clear(iterates->atyUpdate);
    }
    if (iterates->atyAverage)
    {
      // CUPDLP_FREE_VEC(iterates->atyAverage);
      vec_clear(iterates->atyAverage);
    }
    CUPDLP_FREE_VEC(iterates);
  }
}

void resobj_clear(CUPDLPresobj *resobj)
{
  if (resobj)
  {
    if (resobj->primalResidual)
    {
      CUPDLP_FREE_VEC(resobj->primalResidual);
    }
    if (resobj->dualResidual)
    {
      CUPDLP_FREE_VEC(resobj->dualResidual);
    }
    if (resobj->primalResidualAverage)
    {
      CUPDLP_FREE_VEC(resobj->primalResidualAverage);
    }
    if (resobj->dualResidualAverage)
    {
      CUPDLP_FREE_VEC(resobj->dualResidualAverage);
    }
    if (resobj->dSlackPos)
    {
      CUPDLP_FREE_VEC(resobj->dSlackPos);
    }
    if (resobj->dSlackNeg)
    {
      CUPDLP_FREE_VEC(resobj->dSlackNeg);
    }
    if (resobj->dLowerFiltered)
    {
      CUPDLP_FREE_VEC(resobj->dLowerFiltered);
    }
    if (resobj->dUpperFiltered)
    {
      CUPDLP_FREE_VEC(resobj->dUpperFiltered);
    }
    CUPDLP_FREE_VEC(resobj);
  }
}

void stepsize_clear(CUPDLPstepsize *stepsize)
{
  if (stepsize)
  {
    cupdlp_free(stepsize);
  }
}

void timers_clear(CUPDLPtimers *timers)
{
#if !(CUPDLP_CPU)
  cupdlp_printf("%20s %e\n", "Free Device memory", timers->FreeDeviceMemTime);
#endif

  if (timers)
  {
    cupdlp_free(timers);
  }
}

void scaling_clear(CUPDLPscaling *scaling)
{
  if (scaling)
  {
    if (scaling->colScale)
    {
      // cupdlp_free(scaling->colScale);
      CUPDLP_FREE_VEC(scaling->colScale); // now on gpu
    }
    if (scaling->rowScale)
    {
      // cupdlp_free(scaling->rowScale);
      CUPDLP_FREE_VEC(scaling->rowScale); // now on gpu
    }
    cupdlp_free(scaling);
  }
}

cupdlp_int PDHG_Clear(CUPDLPwork *w)
{
  CUPDLPproblem *problem = w->problem;
  CUPDLPsettings *settings = w->settings;
  CUPDLPiterates *iterates = w->iterates;
  CUPDLPresobj *resobj = w->resobj;
  CUPDLPstepsize *stepsize = w->stepsize;
  CUPDLPtimers *timers = w->timers;
  CUPDLPscaling *scaling = w->scaling;

  if (w)
  {
    cupdlp_float begin = getTimeStamp();
#if !(CUPDLP_CPU)

    // CUDAmv *MV = w->MV;
    // if (MV)
    // {
    //     cupdlp_float begin = getTimeStamp();
    //     cuda_free_mv(MV);
    //     timers->FreeDeviceMemTime += getTimeStamp() - begin;
    // }
    CHECK_CUBLAS(cublasDestroy(w->cublashandle))
    CHECK_CUSPARSE(cusparseDestroy(w->cusparsehandle))
    CHECK_CUDA(cudaFree(w->dBuffer))
    if (w->buffer2)
      CUPDLP_FREE_VEC(w->buffer2);
    if (w->buffer3)
      CUPDLP_FREE_VEC(w->buffer3);
#endif
    if (w->colScale)
      CUPDLP_FREE_VEC(w->colScale);
    if (w->rowScale)
      CUPDLP_FREE_VEC(w->rowScale);

    if (w->buffer)
    {
      // CUPDLP_FREE_VEC(w->buffer);
      vec_clear(w->buffer);
    }

    if (problem)
    {
      // problem_clear(problem);
      problem = cupdlp_NULL;
    }

    if (iterates)
    {
      iterates_clear(iterates);
    }

    if (resobj)
    {
      resobj_clear(resobj);
    }

#if !(CUPDLP_CPU)
    timers->FreeDeviceMemTime += getTimeStamp() - begin;
#endif

    if (settings)
    {
      settings_clear(settings);
    }
    if (stepsize)
    {
      stepsize_clear(stepsize);
    }
    if (timers)
    {
      timers_clear(timers);
    }
    if (scaling)
    {
      // scaling_clear(scaling);
      scaling = cupdlp_NULL;
    }
    cupdlp_free(w);
  }

  return 0;
}

void PDHG_PrintPDHGParam(CUPDLPwork *w)
{
  CUPDLPsettings *settings = w->settings;
  CUPDLPstepsize *stepsize = w->stepsize;
  CUPDLPresobj *resobj = w->resobj;
  CUPDLPiterates *iterates = w->iterates;
  CUPDLPscaling *scaling = w->scaling;
  CUPDLPtimers *timers = w->timers;

  cupdlp_printf("\n");

  cupdlp_printf("\n");
  cupdlp_printf("--------------------------------------------------\n");
  cupdlp_printf("CUPDHG Parameters:\n");
  cupdlp_printf("--------------------------------------------------\n");
  cupdlp_printf("\n");

  cupdlp_printf("    nIterLim:          %d\n", settings->nIterLim);
  cupdlp_printf("    dTimeLim (sec):    %.2f\n", settings->dTimeLim);
  cupdlp_printf("    ifScaling:         %d\n", settings->ifScaling);
  // cupdlp_printf("    iScalingMethod:    %d\n", settings->iScalingMethod);)
  cupdlp_printf("    ifRuizScaling:     %d\n", scaling->ifRuizScaling);
  cupdlp_printf("    ifL2Scaling:       %d\n", scaling->ifL2Scaling);
  cupdlp_printf("    ifPcScaling:       %d\n", scaling->ifPcScaling);
  cupdlp_printf("    eLineSearchMethod: %d\n", stepsize->eLineSearchMethod);
  // cupdlp_printf("    dScalingLimit:     %.4e\n", settings->dScalingLimit);
  cupdlp_printf("    dPrimalTol:        %.4e\n", settings->dPrimalTol);
  cupdlp_printf("    dDualTol:          %.4e\n", settings->dDualTol);
  cupdlp_printf("    dGapTol:           %.4e\n", settings->dGapTol);
  cupdlp_printf("    dFeasTol:          %.4e\n", resobj->dFeasTol);
  cupdlp_printf("    eRestartMethod:    %d\n", settings->eRestartMethod);
  cupdlp_printf("    ifPDTEST:          %d\n", settings->ifPDTEST);
  cupdlp_printf("\n");
  cupdlp_printf("--------------------------------------------------\n");
  cupdlp_printf("\n");
}

void PDHG_PrintHugeCUPDHG()
{
  cupdlp_printf("\n");
  cupdlp_printf("  ____ _   _ ____  ____  _     ____\n");
  cupdlp_printf(" / ___| | | |  _ \\|  _ \\| |   |  _ \\\n");
  cupdlp_printf("| |   | | | | |_) | | | | |   | |_) |\n");
  cupdlp_printf("| |___| |_| |  __/| |_| | |___|  __/\n");
  cupdlp_printf(" \\____|\\___/|_|   |____/|_____|_|\n");
  cupdlp_printf("\n");
}

void PDHG_PrintUserParamHelper()
{
  PDHG_PrintHugeCUPDHG();

  cupdlp_printf("CUPDHG User Parameters:\n");
  cupdlp_printf("\n");

  cupdlp_printf("    -h: print this helper\n");
  cupdlp_printf("\n");

  cupdlp_printf("    -nIterLim:  maximum iteration number\n");
  cupdlp_printf("                type:    int\n");
  cupdlp_printf("                default: 10000000\n");
  cupdlp_printf("                range:   >= 0\n");
  cupdlp_printf("\n");

  cupdlp_printf("    -ifScaling: whether to use scaling\n");
  cupdlp_printf("                type:    bool\n");
  cupdlp_printf("                default: true\n");
  cupdlp_printf("                range:   true or false\n");
  cupdlp_printf("\n");

  cupdlp_printf("    -eLineSearchMethod: which line search method to use\n");
  cupdlp_printf(
      "                        0-Fixed, 1-Malitsky-Pock, 2-Adaptive\n");
  cupdlp_printf("                type:    int\n");
  cupdlp_printf("                default: 2\n");
  cupdlp_printf("                range:   0 to 2\n");
  cupdlp_printf("\n");

  cupdlp_printf("    -dPrimalTol: primal tolerance\n");
  cupdlp_printf("                type:    double\n");
  cupdlp_printf("                default: 1e-4\n");
  cupdlp_printf("                range:   >= 0\n");
  cupdlp_printf("\n");

  cupdlp_printf("    -dDualTol: dual tolerance\n");
  cupdlp_printf("                type:    double\n");
  cupdlp_printf("                default: 1e-4\n");
  cupdlp_printf("                range:   >= 0\n");
  cupdlp_printf("\n");

  cupdlp_printf("    -dGapTol: gap tolerance\n");
  cupdlp_printf("                type:    double\n");
  cupdlp_printf("                default: 1e-4\n");
  cupdlp_printf("                range:   >= 0\n");
  cupdlp_printf("\n");

  cupdlp_printf("    -dFeasTol: feasibility tolerance\n");
  cupdlp_printf("                type:    double\n");
  cupdlp_printf("                default: 1e-8\n");
  cupdlp_printf("                range:   >= 0\n");
  cupdlp_printf("\n");

  cupdlp_printf("    -dTimeLim: time limit (in seconds)\n");
  cupdlp_printf("                type:    double\n");
  cupdlp_printf("                default: 3600\n");
  cupdlp_printf("                range:   > 0\n");
  cupdlp_printf("\n");

  cupdlp_printf("    -eRestartMethod: which restart method to use\n");
  cupdlp_printf("                     0-None, 1-GPU\n");
  cupdlp_printf("                type:    int\n");
  cupdlp_printf("                default: 1\n");
  cupdlp_printf("                range:   0 to 1\n");
  cupdlp_printf("\n");

  cupdlp_printf("    -ifRuizScaling: whether to use Ruiz scaling\n");
  cupdlp_printf("                type:    bool\n");
  cupdlp_printf("                default: true\n");
  cupdlp_printf("                range:   true or false\n");
  cupdlp_printf("\n");

  cupdlp_printf("    -ifL2Scaling: whether to use L2 scaling\n");
  cupdlp_printf("                type:    bool\n");
  cupdlp_printf("                default: false\n");
  cupdlp_printf("                range:   true or false\n");
  cupdlp_printf("\n");

  cupdlp_printf("    -ifPcScaling: whether to use Pock-Chambolle scaling\n");
  cupdlp_printf("                type:    bool\n");
  cupdlp_printf("                default: true\n");
  cupdlp_printf("                range:   true or false\n");
  cupdlp_printf("\n");

  cupdlp_printf("    -ifPresolve: whether to presolve problem\n");
  cupdlp_printf("                type:    bool\n");
  cupdlp_printf("                default: true\n");
  cupdlp_printf("                range:   true or false\n");
  cupdlp_printf("\n");

  cupdlp_printf("    -ifPDTEST: whether to use PDTEST\n");
  cupdlp_printf("                type:    int\n");
  cupdlp_printf("                default: 0, use PDHG\n");
  cupdlp_printf("                range:   1 or 0\n");
  cupdlp_printf("\n");
}

cupdlp_retcode getUserParam(int argc, char **argv,
                            cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam,
                            cupdlp_bool *ifChangeFloatParam,
                            cupdlp_float *floatParam)
{
  cupdlp_retcode retcode = RETCODE_OK;

  // set ifChangeIntParam and ifChangeFloatParam to false
  for (cupdlp_int i = 0; i < N_INT_USER_PARAM; ++i)
  {
    ifChangeIntParam[i] = false;
  }

  for (cupdlp_int i = 0; i < N_FLOAT_USER_PARAM; ++i)
  {
    ifChangeFloatParam[i] = false;
  }

  // parse command line arguments
  for (cupdlp_int i = 0; i < argc - 1; ++i)
  {
    if (strcmp(argv[i], "-h") == 0)
    {
      PDHG_PrintUserParamHelper();

      retcode = RETCODE_FAILED;
      goto exit_cleanup;
    }

    if (strcmp(argv[i], "-nIterLim") == 0)
    {
      ifChangeIntParam[N_ITER_LIM] = true;
      intParam[N_ITER_LIM] = atoi(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-ifScaling") == 0)
    {
      ifChangeIntParam[IF_SCALING] = true;
      intParam[IF_SCALING] = atoi(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-iScalingMethod") == 0)
    {
      ifChangeIntParam[I_SCALING_METHOD] = true;
      intParam[I_SCALING_METHOD] = atoi(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-eLineSearchMethod") == 0)
    {
      ifChangeIntParam[E_LINE_SEARCH_METHOD] = true;
      intParam[E_LINE_SEARCH_METHOD] = atoi(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-dScalingLimit") == 0)
    {
      ifChangeFloatParam[D_SCALING_LIMIT] = true;
      floatParam[D_SCALING_LIMIT] = atof(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-dPrimalTol") == 0)
    {
      ifChangeFloatParam[D_PRIMAL_TOL] = true;
      floatParam[D_PRIMAL_TOL] = atof(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-dDualTol") == 0)
    {
      ifChangeFloatParam[D_DUAL_TOL] = true;
      floatParam[D_DUAL_TOL] = atof(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-dGapTol") == 0)
    {
      ifChangeFloatParam[D_GAP_TOL] = true;
      floatParam[D_GAP_TOL] = atof(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-dFeasTol") == 0)
    {
      ifChangeFloatParam[D_FEAS_TOL] = true;
      floatParam[D_FEAS_TOL] = atof(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-dTimeLim") == 0)
    {
      ifChangeFloatParam[D_TIME_LIM] = true;
      floatParam[D_TIME_LIM] = atof(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-eRestartMethod") == 0)
    {
      ifChangeIntParam[E_RESTART_METHOD] = true;
      intParam[E_RESTART_METHOD] = atoi(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-ifRuizScaling") == 0)
    {
      ifChangeIntParam[IF_RUIZ_SCALING] = true;
      intParam[IF_RUIZ_SCALING] = atoi(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-ifL2Scaling") == 0)
    {
      ifChangeIntParam[IF_L2_SCALING] = true;
      intParam[IF_L2_SCALING] = atoi(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-ifPcScaling") == 0)
    {
      ifChangeIntParam[IF_PC_SCALING] = true;
      intParam[IF_PC_SCALING] = atoi(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-nLogInt") == 0)
    {
      ifChangeIntParam[N_LOG_INTERVAL] = true;
      intParam[N_LOG_INTERVAL] = atoi(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-ifPresolve") == 0)
    {
      ifChangeIntParam[IF_PRESOLVE] = true;
      intParam[IF_PRESOLVE] = atoi(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-ifPDTEST") == 0)
    {
      ifChangeIntParam[IF_PDTEST] = true;
      intParam[IF_PDTEST] = atoi(argv[i + 1]);
      cupdlp_printf("提取的ifPDTEST: %d\n", intParam[IF_PDTEST]);
    }
    else if (strcmp(argv[i], "-bestID") == 0)
    {
      ifChangeIntParam[BEST_ID] = true;
      intParam[BEST_ID] = atoi(argv[i + 1]);
      cupdlp_printf("bestID: %d\n", intParam[BEST_ID]);
    }
  }

  if (strcmp(argv[argc - 1], "-h") == 0)
  {
    PDHG_PrintUserParamHelper();

    retcode = RETCODE_FAILED;
    goto exit_cleanup;
  }

exit_cleanup:
  return retcode;
}

cupdlp_retcode settings_SetUserParam(CUPDLPsettings *settings,
                                     cupdlp_bool *ifChangeIntParam,
                                     cupdlp_int *intParam,
                                     cupdlp_bool *ifChangeFloatParam,
                                     cupdlp_float *floatParam)
{
  cupdlp_retcode retcode = RETCODE_OK;

  if (ifChangeIntParam[N_ITER_LIM])
  {
    settings->nIterLim = intParam[N_ITER_LIM];
  }

  if (ifChangeIntParam[N_LOG_INTERVAL])
  {
    settings->nLogInterval = intParam[N_LOG_INTERVAL];
  }

  if (ifChangeIntParam[IF_SCALING])
  {
    settings->ifScaling = intParam[IF_SCALING];
  }

  if (ifChangeIntParam[I_SCALING_METHOD])
  {
    settings->iScalingMethod = intParam[I_SCALING_METHOD];
  }

  if (ifChangeFloatParam[D_SCALING_LIMIT])
  {
    settings->dScalingLimit = floatParam[D_SCALING_LIMIT];
  }

  if (ifChangeFloatParam[D_PRIMAL_TOL])
  {
    settings->dPrimalTol = floatParam[D_PRIMAL_TOL];
  }

  if (ifChangeFloatParam[D_DUAL_TOL])
  {
    settings->dDualTol = floatParam[D_DUAL_TOL];
  }

  if (ifChangeFloatParam[D_GAP_TOL])
  {
    settings->dGapTol = floatParam[D_GAP_TOL];
  }

  if (ifChangeFloatParam[D_TIME_LIM])
  {
    settings->dTimeLim = floatParam[D_TIME_LIM];
  }

  if (ifChangeIntParam[E_RESTART_METHOD])
  {
    settings->eRestartMethod = intParam[E_RESTART_METHOD];
  }
  if (ifChangeIntParam[IF_PDTEST])
  {
    settings->ifPDTEST = intParam[IF_PDTEST];
  }
  if (ifChangeIntParam[BEST_ID])
  {
    settings->bestID = intParam[BEST_ID];
  }

exit_cleanup:
  return retcode;
}

cupdlp_retcode resobj_SetUserParam(CUPDLPresobj *resobj,
                                   cupdlp_bool *ifChangeIntParam,
                                   cupdlp_int *intParam,
                                   cupdlp_bool *ifChangeFloatParam,
                                   cupdlp_float *floatParam)
{
  cupdlp_retcode retcode = RETCODE_OK;

  if (ifChangeFloatParam[D_FEAS_TOL])
  {
    resobj->dFeasTol = floatParam[D_FEAS_TOL];
  }

exit_cleanup:
  return retcode;
}

cupdlp_retcode iterates_SetUserParam(CUPDLPiterates *iterates,
                                     cupdlp_bool *ifChangeIntParam,
                                     cupdlp_int *intParam,
                                     cupdlp_bool *ifChangeFloatParam,
                                     cupdlp_float *floatParam)
{
  cupdlp_retcode retcode = RETCODE_OK;

exit_cleanup:
  return retcode;
}

cupdlp_retcode iterates_SetUserParam_PDTEST(PDTESTiterates *iterates,
                                            cupdlp_bool *ifChangeIntParam,
                                            cupdlp_int *intParam,
                                            cupdlp_bool *ifChangeFloatParam,
                                            cupdlp_float *floatParam)
{
  cupdlp_retcode retcode = RETCODE_OK;

exit_cleanup:
  return retcode;
}
cupdlp_retcode stepsize_SetUserParam(CUPDLPstepsize *stepsize,
                                     cupdlp_bool *ifChangeIntParam,
                                     cupdlp_int *intParam,
                                     cupdlp_bool *ifChangeFloatParam,
                                     cupdlp_float *floatParam)
{
  cupdlp_retcode retcode = RETCODE_OK;

  if (ifChangeIntParam[E_LINE_SEARCH_METHOD])
  {
    stepsize->eLineSearchMethod = intParam[E_LINE_SEARCH_METHOD];
  }

exit_cleanup:
  return retcode;
}

cupdlp_retcode scaling_SetUserParam(CUPDLPscaling *scaling,
                                    cupdlp_bool *ifChangeIntParam,
                                    cupdlp_int *intParam,
                                    cupdlp_bool *ifChangeFloatParam,
                                    cupdlp_float *floatParam)
{
  cupdlp_retcode retcode = RETCODE_OK;

  if (ifChangeIntParam[IF_RUIZ_SCALING])
  {
    scaling->ifRuizScaling = intParam[IF_RUIZ_SCALING];
  }

  if (ifChangeIntParam[IF_L2_SCALING])
  {
    scaling->ifL2Scaling = intParam[IF_L2_SCALING];
  }

  if (ifChangeIntParam[IF_PC_SCALING])
  {
    scaling->ifPcScaling = intParam[IF_PC_SCALING];
  }

exit_cleanup:
  return retcode;
}

cupdlp_retcode timers_SetUserParam(CUPDLPtimers *timers,
                                   cupdlp_bool *ifChangeIntParam,
                                   cupdlp_int *intParam,
                                   cupdlp_bool *ifChangeFloatParam,
                                   cupdlp_float *floatParam)
{
  cupdlp_retcode retcode = RETCODE_OK;

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDHG_SetUserParam(CUPDLPwork *w, cupdlp_bool *ifChangeIntParam,
                                 cupdlp_int *intParam,
                                 cupdlp_bool *ifChangeFloatParam,
                                 cupdlp_float *floatParam)
{
  cupdlp_retcode retcode = RETCODE_OK;
  printf("PDHG_SetUserParam\n");
  cupdlp_printf("PDHG_SetUserParam, w->settings->dTimeLim: %f\n", w->settings->dTimeLim);
  CUPDLPsettings *settings = w->settings;
  printf("settings point!");
  CUPDLPstepsize *stepsize = w->stepsize;
  printf("stepsize point!");
  CUPDLPresobj *resobj = w->resobj;
  printf("resobj point!");
  CUPDLPiterates *iterates = w->iterates;
  CUPDLPscaling *scaling = w->scaling;
  CUPDLPtimers *timers = w->timers;

  CUPDLP_CALL(settings_SetUserParam(settings, ifChangeIntParam, intParam,
                                    ifChangeFloatParam, floatParam));
  CUPDLP_CALL(stepsize_SetUserParam(stepsize, ifChangeIntParam, intParam,
                                    ifChangeFloatParam, floatParam));
  CUPDLP_CALL(resobj_SetUserParam(resobj, ifChangeIntParam, intParam,
                                  ifChangeFloatParam, floatParam));
  CUPDLP_CALL(iterates_SetUserParam(iterates, ifChangeIntParam, intParam,
                                    ifChangeFloatParam, floatParam));
  CUPDLP_CALL(scaling_SetUserParam(scaling, ifChangeIntParam, intParam,
                                   ifChangeFloatParam, floatParam));
  CUPDLP_CALL(timers_SetUserParam(timers, ifChangeIntParam, intParam,
                                  ifChangeFloatParam, floatParam));

  PDHG_PrintPDHGParam(w);

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDTEST_SetUserParam(CUPDLPwork *w, cupdlp_bool *ifChangeIntParam,
                                   cupdlp_int *intParam,
                                   cupdlp_bool *ifChangeFloatParam,
                                   cupdlp_float *floatParam)
{
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLPsettings *settings = w->settings;
  CUPDLPstepsize *stepsize = w->stepsize;
  CUPDLPresobj *resobj = w->resobj;
  PDTESTiterates *iterates = w->PDTESTiterates;
  CUPDLPscaling *scaling = w->scaling;
  CUPDLPtimers *timers = w->timers;

  CUPDLP_CALL(settings_SetUserParam(settings, ifChangeIntParam, intParam,
                                    ifChangeFloatParam, floatParam));
  CUPDLP_CALL(stepsize_SetUserParam(stepsize, ifChangeIntParam, intParam,
                                    ifChangeFloatParam, floatParam));
  CUPDLP_CALL(resobj_SetUserParam(resobj, ifChangeIntParam, intParam,
                                  ifChangeFloatParam, floatParam));
  CUPDLP_CALL(iterates_SetUserParam_PDTEST(iterates, ifChangeIntParam, intParam,
                                           ifChangeFloatParam, floatParam));
  CUPDLP_CALL(scaling_SetUserParam(scaling, ifChangeIntParam, intParam,
                                   ifChangeFloatParam, floatParam));
  CUPDLP_CALL(timers_SetUserParam(timers, ifChangeIntParam, intParam,
                                  ifChangeFloatParam, floatParam));

  PDHG_PrintPDHGParam(w);

exit_cleanup:
  return retcode;
}

cupdlp_retcode settings_Alloc(CUPDLPsettings *settings)
{
  cupdlp_retcode retcode = RETCODE_OK;
  settings->nIterLim = INFINITY;
  settings->nLogInterval = 100;
  // settings->dTimeLim = INFINITY;
  settings->dTimeLim = 3600;
  settings->ifScaling = true;
  settings->iScalingMethod = 3; // no use
  settings->dScalingLimit = 5;  // no use
  settings->eRestartMethod = PDHG_GPU_RESTART;

  // termination criteria
  settings->dPrimalTol = 1e-4;
  settings->dDualTol = 1e-4;
  settings->dGapTol = 1e-4;

  return retcode;
}

cupdlp_retcode resobj_Alloc(CUPDLPresobj *resobj, CUPDLPproblem *problem,
                            cupdlp_int ncols, cupdlp_int nrows)
{
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLP_INIT_ZERO_VEC(resobj->primalResidual, nrows);
  CUPDLP_INIT_ZERO_VEC(resobj->dualResidual, ncols);
  CUPDLP_INIT_ZERO_VEC(resobj->primalResidualAverage, nrows);
  CUPDLP_INIT_ZERO_VEC(resobj->dualResidualAverage, ncols);
  CUPDLP_INIT_ZERO_VEC(resobj->dSlackPos, ncols);
  CUPDLP_INIT_ZERO_VEC(resobj->dSlackNeg, ncols);
  CUPDLP_INIT_ZERO_VEC(resobj->dLowerFiltered, ncols);
  CUPDLP_INIT_ZERO_VEC(resobj->dUpperFiltered, ncols);

  // need to translate to cuda type
  // for (int i = 0; i < ncols; i++)
  // {
  //     resobj->dLowerFiltered[i] = problem->lower[i] > -INFINITY ?
  //     problem->lower[i] : 0.0; resobj->dUpperFiltered[i] = problem->upper[i]
  //     < +INFINITY ? problem->upper[i] : 0.0;
  // }

  cupdlp_filterlb(resobj->dLowerFiltered, problem->lower, -INFINITY, ncols);
  cupdlp_filterub(resobj->dUpperFiltered, problem->upper, +INFINITY, ncols);

  // initialization
  resobj->dFeasTol = 1e-8;
  resobj->dPrimalObj = 0.0;
  resobj->dDualObj = 0.0;
  resobj->dDualityGap = 0.0;
  resobj->dComplementarity = 0.0;
  resobj->dPrimalFeas = 0.0;
  resobj->dDualFeas = 0.0;
  resobj->dRelObjGap = 0.0;
  resobj->dPrimalObjAverage = 0.0;
  resobj->dDualObjAverage = 0.0;
  resobj->dDualityGapAverage = 0.0;
  resobj->dComplementarityAverage = 0.0;
  resobj->dPrimalFeasAverage = 0.0;
  resobj->dDualFeasAverage = 0.0;
  resobj->dRelObjGapAverage = 0.0;
  resobj->dPrimalFeasLastRestart = 0.0;
  resobj->dDualFeasLastRestart = 0.0;
  resobj->dDualityGapLastRestart = 0.0;
  resobj->dPrimalFeasLastCandidate = 0.0;
  resobj->dDualFeasLastCandidate = 0.0;
  resobj->dDualityGapLastCandidate = 0.0;

  resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
  resobj->termIterate = LAST_ITERATE;

  // todo, pass work
  //   cupdlp_twoNorm(problem->cost, ncols, &resobj->dNormCost);
  //   twoNorm(problem->rhs, nrows);

exit_cleanup:
  return retcode;
}

cupdlp_retcode iterates_Alloc(CUPDLPiterates *iterates, cupdlp_int ncols,
                              cupdlp_int nrows)
{
  cupdlp_retcode retcode = RETCODE_OK;

  iterates->nCols = ncols;
  iterates->nRows = nrows;

  CUPDLP_INIT_ZERO_VEC(iterates->xSum, ncols);
  CUPDLP_INIT_ZERO_VEC(iterates->ySum, nrows);
  CUPDLP_INIT_ZERO_VEC(iterates->xLastRestart, ncols);
  CUPDLP_INIT_ZERO_VEC(iterates->yLastRestart, nrows);

  CUPDLP_INIT(iterates->x, 1);
  CUPDLP_INIT(iterates->xUpdate, 1);
  CUPDLP_INIT(iterates->xAverage, 1);
  CUPDLP_INIT(iterates->y, 1);
  CUPDLP_INIT(iterates->yUpdate, 1);
  CUPDLP_INIT(iterates->yAverage, 1);
  CUPDLP_INIT(iterates->ax, 1);
  CUPDLP_INIT(iterates->axUpdate, 1);
  CUPDLP_INIT(iterates->axAverage, 1);
  CUPDLP_INIT(iterates->aty, 1);
  CUPDLP_INIT(iterates->atyUpdate, 1);
  CUPDLP_INIT(iterates->atyAverage, 1);

  CUPDLP_CALL(vec_Alloc(iterates->x, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->xUpdate, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->xAverage, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->y, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->yUpdate, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->yAverage, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->ax, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->axUpdate, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->axAverage, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->aty, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->atyUpdate, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->atyAverage, ncols));

  // initialization
  iterates->iLastRestartIter = 0;
  iterates->dLastRestartDualityGap = 0.0;
  iterates->dLastRestartBeta = 0.0;

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDTESTiterates_Alloc(PDTESTiterates *iterates, cupdlp_int ncols,
                                    cupdlp_int nrows)
{
  cupdlp_retcode retcode = RETCODE_OK;

  iterates->nCols = ncols;
  iterates->nRows = nrows;

  CUPDLP_INIT_ZERO_VEC(iterates->xSum, ncols);
  CUPDLP_INIT_ZERO_VEC(iterates->ySum, nrows);
  // 我的修改///////////////////////////////////////////////////
  CUPDLP_INIT_ZERO_VEC(iterates->x_agSum, ncols);
  CUPDLP_INIT_ZERO_VEC(iterates->y_agSum, nrows);
  /////////////////////////////////////////////////////
  CUPDLP_INIT_ZERO_VEC(iterates->xLastRestart, ncols);
  CUPDLP_INIT_ZERO_VEC(iterates->yLastRestart, nrows);

  CUPDLP_INIT(iterates->x, 1);
  CUPDLP_INIT(iterates->xUpdate, 1);
  CUPDLP_INIT(iterates->xAverage, 1);
  CUPDLP_INIT(iterates->y, 1);
  CUPDLP_INIT(iterates->yUpdate, 1);
  CUPDLP_INIT(iterates->yAverage, 1);
  CUPDLP_INIT(iterates->ax, 1);
  CUPDLP_INIT(iterates->axUpdate, 1);
  CUPDLP_INIT(iterates->axAverage, 1);
  CUPDLP_INIT(iterates->aty, 1);
  CUPDLP_INIT(iterates->atyUpdate, 1);
  CUPDLP_INIT(iterates->atyAverage, 1);
  // 我的修改///////////////////////////////////////////////////
  CUPDLP_INIT(iterates->x_ag, 1);
  CUPDLP_INIT(iterates->y_ag, 1);
  CUPDLP_INIT(iterates->x_bar, 1);
  CUPDLP_INIT(iterates->x_md, 1);
  CUPDLP_INIT(iterates->ax_bar, 1);
  CUPDLP_INIT(iterates->ax_ag, 1);
  CUPDLP_INIT(iterates->aty_ag, 1);
  CUPDLP_INIT(iterates->x_agUpdate, 1);
  CUPDLP_INIT(iterates->y_agUpdate, 1);
  CUPDLP_INIT(iterates->x_barUpdate, 1);
  CUPDLP_INIT(iterates->ax_agUpdate, 1);
  CUPDLP_INIT(iterates->aty_agUpdate, 1);
  CUPDLP_INIT(iterates->x_agAverage, 1);
  CUPDLP_INIT(iterates->y_agAverage, 1);
  CUPDLP_INIT(iterates->ax_agAverage, 1);
  CUPDLP_INIT(iterates->aty_agAverage, 1);
  /////////////////////////////////////////////////////

  CUPDLP_CALL(vec_Alloc(iterates->x, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->xUpdate, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->xAverage, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->y, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->yUpdate, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->yAverage, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->ax, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->axUpdate, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->axAverage, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->aty, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->atyUpdate, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->atyAverage, ncols));
  // 我的修改///////////////////////////////////////////////////
  CUPDLP_CALL(vec_Alloc(iterates->x_ag, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->x_bar, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->x_md, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->x_agUpdate, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->x_barUpdate, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->aty_ag, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->y_ag, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->y_agUpdate, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->ax_bar, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->ax_ag, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->ax_agUpdate, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->aty_agUpdate, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->x_agAverage, ncols));
  CUPDLP_CALL(vec_Alloc(iterates->y_agAverage, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->ax_agAverage, nrows));
  CUPDLP_CALL(vec_Alloc(iterates->aty_agAverage, ncols));
  /////////////////////////////////////////////////////

  // initialization
  iterates->iLastRestartIter = 0;
  iterates->dLastRestartDualityGap = 0.0;
  iterates->dLastRestartBeta = 0.0;

exit_cleanup:
  return retcode;
}

cupdlp_retcode stepsize_Alloc(CUPDLPstepsize *stepsize)
{
  cupdlp_retcode retcode = RETCODE_OK;

  stepsize->eLineSearchMethod = PDHG_ADAPTIVE_LINESEARCH;

  // initialization
  stepsize->nStepSizeIter = 0;
  stepsize->dPrimalStep = 0.0;
  stepsize->dDualStep = 0.0;
  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;
  stepsize->dBeta = 0.0;
  stepsize->dTheta = 0.0;

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDTEST_stepsize_Alloc(CUPDLPstepsize *stepsize)
{
  cupdlp_retcode retcode = RETCODE_OK;

  stepsize->eLineSearchMethod = PDHG_ADAPTIVE_LINESEARCH;

  // initialization
  stepsize->nStepSizeIter = 0;
  stepsize->dPrimalStep = 0.0;
  stepsize->dDualStep = 0.0;
  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;
  stepsize->dBeta = 0.0;
  stepsize->dTheta = 0.0;
  stepsize->dStepSizeLastIter = 0.0;
  stepsize->dStepSizeSum = 0.0;

exit_cleanup:
  return retcode;
}

cupdlp_retcode scaling_Alloc(CUPDLPscaling *scaling, CUPDLPproblem *problem,
                             cupdlp_int ncols, cupdlp_int nrows)
{
  cupdlp_retcode retcode = RETCODE_OK;
  scaling->ifScaled = 0;

  CUPDLP_INIT(scaling->colScale, ncols);
  CUPDLP_INIT(scaling->rowScale, nrows);

  scaling->ifRuizScaling = 1;
  scaling->ifL2Scaling = 0;
  scaling->ifPcScaling = 1;

  scaling->dNormCost = twoNorm(problem->cost, problem->nCols);
  scaling->dNormRhs = twoNorm(problem->rhs, problem->nRows);

exit_cleanup:
  return retcode;
}

cupdlp_retcode timers_Alloc(CUPDLPtimers *timers)
{
  cupdlp_retcode retcode = RETCODE_OK;

  timers->nIter = 0;
  timers->dSolvingTime = 0.0;
  timers->dSolvingBeg = 0.0;
  timers->dScalingTime = 0.0;
  timers->dPresolveTime = 0.0;

#if PDHG_USE_TIMERS
  timers->dAtyTime = 0.0;
  timers->dAxTime = 0.0;
  timers->dComputeResidualsTime = 0.0;
  timers->dUpdateIterateTime = 0.0;
  timers->nAtyCalls = 0;
  timers->nAxCalls = 0;
  timers->nComputeResidualsCalls = 0;
  timers->nUpdateIterateCalls = 0;
  timers->dMatVecMultiplyTime = 0.0;
  timers->dVecVecAddTime = 0.0;
  timers->dIterTime = 0.0;
  timers->dAdapStepTime = 0.0;
#endif
#if !(CUPDLP_CPU)
  // GPU timers
  timers->AllocMem_CopyMatToDeviceTime = 0.0;
  timers->CopyVecToDeviceTime = 0.0;
  timers->DeviceMatVecProdTime = 0.0;
  timers->CopyVecToHostTime = 0.0;
  timers->FreeDeviceMemTime = 0.0;
  timers->CudaPrepareTime = 0.0;
#endif

exit_cleanup:
  return retcode;
}

cupdlp_retcode vec_Alloc(CUPDLPvec *vec, cupdlp_int n)
{
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLP_INIT_ZERO_VEC(vec->data, n);
  vec->len = n;
#if !(CUPDLP_CPU)
  CHECK_CUSPARSE(
      cusparseCreateDnVec(&vec->cuda_vec, n, vec->data, CudaComputeType));
#endif

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDHG_Alloc(CUPDLPwork *w)
{
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLP_INIT(w->settings, 1);
  CUPDLP_INIT(w->resobj, 1);
  CUPDLP_INIT(w->iterates, 1);
  CUPDLP_INIT(w->stepsize, 1);

  CUPDLP_INIT(w->timers, 1);
  CUPDLP_CALL(timers_Alloc(w->timers));

  cupdlp_float begin = getTimeStamp();
  // buffer
  CUPDLP_INIT(w->buffer, 1);
  CUPDLP_CALL(vec_Alloc(w->buffer, w->problem->data->nRows));
  CUPDLP_INIT_ZERO_VEC(w->buffer2,
                       MAX(w->problem->data->nCols, w->problem->data->nRows));
  CUPDLP_INIT_ZERO_VEC(w->buffer3,
                       MAX(w->problem->data->nCols, w->problem->data->nRows));

  // for scaling
  CUPDLP_INIT_ZERO_VEC(w->colScale, w->problem->data->nCols);
  CUPDLP_INIT_ZERO_VEC(w->rowScale, w->problem->data->nRows);

  CUPDLP_CALL(settings_Alloc(w->settings));
  CUPDLP_CALL(resobj_Alloc(w->resobj, w->problem, w->problem->data->nCols,
                           w->problem->data->nRows));
  CUPDLP_CALL(iterates_Alloc(w->iterates, w->problem->data->nCols,
                             w->problem->data->nRows));
  CUPDLP_CALL(stepsize_Alloc(w->stepsize));

#if !(CUPDLP_CPU)
  //   CHECK_CUSPARSE(cusparseCreate(&w->cusparsehandle));
  //   CHECK_CUBLAS(cublasCreate(&w->cublashandle));
  cuda_alloc_MVbuffer(
      //   w->problem->data->matrix_format,
      w->cusparsehandle, w->problem->data->csc_matrix->cuda_csc,
      w->iterates->x->cuda_vec, w->iterates->ax->cuda_vec,
      w->problem->data->csr_matrix->cuda_csr, w->iterates->y->cuda_vec,
      w->iterates->aty->cuda_vec, &w->dBuffer);
  w->timers->AllocMem_CopyMatToDeviceTime += getTimeStamp() - begin;
#endif

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDTEST_Alloc(CUPDLPwork *w)
{
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLP_INIT(w->settings, 1);
  CUPDLP_INIT(w->resobj, 1);
  /////////////////////////////////////////////
  // PDTESTiterates or iterates
  // CUPDLP_INIT(w->iterates, 1);
  CUPDLP_INIT(w->PDTESTiterates, 1);

  CUPDLP_INIT(w->stepsize, 1);

  CUPDLP_INIT(w->timers, 1);
  CUPDLP_CALL(timers_Alloc(w->timers));

  cupdlp_float begin = getTimeStamp();
  // buffer
  CUPDLP_INIT(w->buffer, 1);
  CUPDLP_CALL(vec_Alloc(w->buffer, w->problem->data->nRows));
  CUPDLP_INIT_ZERO_VEC(w->buffer2,
                       MAX(w->problem->data->nCols, w->problem->data->nRows));
  CUPDLP_INIT_ZERO_VEC(w->buffer3,
                       MAX(w->problem->data->nCols, w->problem->data->nRows));

  // for scaling
  CUPDLP_INIT_ZERO_VEC(w->colScale, w->problem->data->nCols);
  CUPDLP_INIT_ZERO_VEC(w->rowScale, w->problem->data->nRows);

  CUPDLP_CALL(settings_Alloc(w->settings));
  CUPDLP_CALL(resobj_Alloc(w->resobj, w->problem, w->problem->data->nCols,
                           w->problem->data->nRows));
  /////////////////////////////////////////////
  // PDTESTiterates or iterates
  // CUPDLP_CALL(iterates_Alloc(w->iterates, w->problem->data->nCols, w->problem->data->nRows));
  CUPDLP_CALL(PDTESTiterates_Alloc(w->PDTESTiterates, w->problem->data->nCols, w->problem->data->nRows));

  CUPDLP_CALL(PDTEST_stepsize_Alloc(w->stepsize));

#if !(CUPDLP_CPU)
  //   CHECK_CUSPARSE(cusparseCreate(&w->cusparsehandle));
  //   CHECK_CUBLAS(cublasCreate(&w->cublashandle));
  cuda_alloc_MVbuffer(
      //   w->problem->data->matrix_format,
      w->cusparsehandle, w->problem->data->csc_matrix->cuda_csc,
      w->PDTESTiterates->x->cuda_vec, w->PDTESTiterates->ax->cuda_vec,
      w->problem->data->csr_matrix->cuda_csr, w->PDTESTiterates->y->cuda_vec,
      w->PDTESTiterates->aty->cuda_vec, &w->dBuffer);
  w->timers->AllocMem_CopyMatToDeviceTime += getTimeStamp() - begin;
#endif

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDHG_Create(CUPDLPwork **ww, CUPDLPproblem *lp,
                           CUPDLPscaling *scaling)
{
  cupdlp_retcode retcode = RETCODE_OK;
  CUPDLP_INIT_ZERO(*ww, 1);

  CUPDLPwork *w = *ww;
  w->problem = lp;
  w->scaling = scaling;

exit_cleanup:
  return retcode;
}

void PDHG_Destroy(CUPDLPwork **w)
{
  if (w && *w)
  {
    PDHG_Clear(*w);
#if !(CUPDLP_CPU)
    cudaDeviceReset();
#endif
  }
}

void PDHG_Init_Data(CUPDLPwork *work) {}

double my_clock(void)
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return (1e-06 * t.tv_usec + t.tv_sec);
}

double getTimeStamp(void) { return my_clock(); }

void dense_copy(CUPDLPdense *dst, CUPDLPdense *src)
{
  dst->nRows = src->nRows;
  dst->nCols = src->nCols;
  cupdlp_copy(dst->data, src->data, cupdlp_float, src->nRows * src->nCols);

  return;
}
// void csr_copy(CUPDLPcsr *dst, CUPDLPcsr *src)
// {
//   dst->nRows = src->nRows;
//   dst->nCols = src->nCols;
//   dst->nMatElem = src->nMatElem;
//   cupdlp_copy(dst->rowMatBeg, src->rowMatBeg, cupdlp_int, src->nRows + 1);
//   cupdlp_copy(dst->rowMatIdx, src->rowMatIdx, cupdlp_int, src->nMatElem);
//   cupdlp_copy(dst->rowMatElem, src->rowMatElem, cupdlp_float, src->nMatElem);

//   return;
// }
void csr_copy(CUPDLPcsr *dst, CUPDLPcsr *src)
{
  dst->nRows = src->nRows;
  dst->nCols = src->nCols;
  dst->nMatElem = src->nMatElem;
  CUPDLP_COPY_VEC(dst->rowMatBeg, src->rowMatBeg, cupdlp_int, src->nRows + 1);
  CUPDLP_COPY_VEC(dst->rowMatIdx, src->rowMatIdx, cupdlp_int, src->nMatElem);
  CUPDLP_COPY_VEC(dst->rowMatElem, src->rowMatElem, cupdlp_float, src->nMatElem);

#if !(CUPDLP_CPU)
  // Pointer to GPU csc matrix
  CHECK_CUSPARSE(cusparseCreateCsr(
      &dst->cuda_csr, src->nRows, src->nCols, src->nMatElem, dst->rowMatBeg,
      dst->rowMatIdx, dst->rowMatElem, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
      CUSPARSE_INDEX_BASE_ZERO, CudaComputeType));
#endif

  return;
}

cupdlp_int csc_copy(CUPDLPcsc *dst, CUPDLPcsc *src)
{
  dst->nRows = src->nRows;
  dst->nCols = src->nCols;
  dst->nMatElem = src->nMatElem;
  CUPDLP_COPY_VEC(dst->colMatBeg, src->colMatBeg, cupdlp_int, src->nCols + 1);
  CUPDLP_COPY_VEC(dst->colMatIdx, src->colMatIdx, cupdlp_int, src->nMatElem);
  CUPDLP_COPY_VEC(dst->colMatElem, src->colMatElem, cupdlp_float,
                  src->nMatElem);

#if !(CUPDLP_CPU)
  // Pointer to GPU csc matrix
  CHECK_CUSPARSE(cusparseCreateCsc(
      &dst->cuda_csc, src->nRows, src->nCols, src->nMatElem, dst->colMatBeg,
      dst->colMatIdx, dst->colMatElem, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
      CUSPARSE_INDEX_BASE_ZERO, CudaComputeType));
#endif

  return 0;
}
// void csr2csc(CUPDLPcsc *csc, CUPDLPcsr *csr)
// {
//   cupdlp_dcs *cs_csr =
//       cupdlp_dcs_spalloc(csr->nCols, csc->nRows, csc->nMatElem, 1, 0);
//   cupdlp_copy(cs_csr->p, csr->rowMatBeg, cupdlp_int, csr->nRows + 1);
//   cupdlp_copy(cs_csr->i, csr->rowMatIdx, cupdlp_int, csr->nMatElem);
//   cupdlp_copy(cs_csr->x, csr->rowMatElem, cupdlp_float, csr->nMatElem);

//   cupdlp_dcs *cs_csc = cupdlp_dcs_transpose(cs_csr, 1);
//   csc->nCols = cs_csc->m;
//   csc->nRows = cs_csc->n;
//   csc->nMatElem = cs_csc->nzmax;
//   cupdlp_copy(csc->colMatBeg, cs_csc->p, cupdlp_int, cs_csc->n + 1);
//   cupdlp_copy(csc->colMatIdx, cs_csc->i, cupdlp_int, cs_csc->nzmax);
//   cupdlp_copy(csc->colMatElem, cs_csc->x, cupdlp_float, cs_csc->nzmax);

//   cupdlp_dcs_free(cs_csc);
//   cupdlp_dcs_free(cs_csr);

//   return;
// }
void csr2csc(CUPDLPcsc *csc, CUPDLPcsr *csr)
{
  cupdlp_printf("csr->nRows = %d, csr->nCols = %d, csr->nMatElem = %d\n", csr->nRows, csr->nCols, csr->nMatElem);
  // cupdlp_dcs *cs_csr = cupdlp_dcs_spalloc(csr->nCols, csr->nRows, csr->nMatElem, 1, 0);
  cupdlp_dcs *cs_csr = cupdlp_dcs_sp_csr_alloc(csr->nRows, csr->nCols, csr->nMatElem, 1, 0);
  cupdlp_printf("cs_csr->m = %d, cs_csr->n = %d, cs_csr->nzmax = %d\n", cs_csr->m, cs_csr->n, cs_csr->nzmax);
  cupdlp_copy(cs_csr->p, csr->rowMatBeg, cupdlp_int, csr->nRows + 1);
  cupdlp_copy(cs_csr->i, csr->rowMatIdx, cupdlp_int, csr->nMatElem);
  cupdlp_copy(cs_csr->x, csr->rowMatElem, cupdlp_float, csr->nMatElem);

  cupdlp_dcs *cs_csc = cupdlp_dcs_csr_transpose(cs_csr, 1);
  cupdlp_printf("cs_csc->m = %d, cs_csc->n = %d, cs_csc->nzmax = %d\n", cs_csc->m, cs_csc->n, cs_csc->nzmax);
  csc->nCols = cs_csc->m;
  csc->nRows = cs_csc->n;
  csc->nMatElem = cs_csc->nzmax;
  cupdlp_printf("csc->nRows = %d, csc->nCols = %d, csc->nMatElem = %d\n", csc->nRows, csc->nCols, csc->nMatElem);
  cupdlp_printf("cs_csc->nzmax = %d\n", cs_csc->nzmax);
  // for (int i = 0; i < cs_csc->nzmax; i++)
  // {
  //   printf("colMatIdx,cs_csc->i[%d] = %d\n", i, cs_csc->i[i]);
  // }
  // for (int i = 0; i < cs_csc->nzmax; i++)
  // {
  //   printf("colMatElem,cs_csc->x[%d] = %f\n", i, cs_csc->x[i]);
  // }
  // for (int i = 0; i < cs_csc->m + 1; i++)
  // {
  //   printf("colMatBeg,cs_csc->p[%d] = %d\n", i, cs_csc->p[i]);
  // }
  // for (int i = 0; i < cs_csc->m + 1; i++)
  // {
  //   printf("cs_csc->p[%d] = %d\n", i, cs_csc->p[i]);
  // }
  // for (int i = 0; i < cs_csc->nzmax; i++)
  // {
  //   printf("cs_csc->x[%d] = %f\n", i, cs_csc->x[i]);
  // }
  CUPDLP_COPY_VEC(csc->colMatBeg, cs_csc->p, cupdlp_int, cs_csc->m + 1);
  CUPDLP_COPY_VEC(csc->colMatIdx, cs_csc->i, cupdlp_int, cs_csc->nzmax);
  CUPDLP_COPY_VEC(csc->colMatElem, cs_csc->x, cupdlp_float, cs_csc->nzmax);

  cupdlp_dcs_free(cs_csc);
  cupdlp_dcs_free(cs_csr);

#if !(CUPDLP_CPU)
  // Pointer to GPU csc matrix
  CHECK_CUSPARSE(cusparseCreateCsr(
      &csc->cuda_csc, csc->nRows, csc->nCols, csc->nMatElem, csc->colMatBeg,
      csc->colMatIdx, csc->colMatElem, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
      CUSPARSE_INDEX_BASE_ZERO, CudaComputeType));
#endif

  return;
}

cupdlp_int csc2csr(CUPDLPcsr *csr, CUPDLPcsc *csc)
{
  // The transpose may need to be done on the GPU
  // Currently, it is done on the CPU
  cupdlp_printf("csc->nRows = %d, csc->nCols = %d, csc->nMatElem = %d\n", csc->nRows, csc->nCols, csc->nMatElem);
  cupdlp_dcs *cs_csc =
      cupdlp_dcs_spalloc(csc->nRows, csc->nCols, csc->nMatElem, 1, 0);
  cupdlp_copy(cs_csc->p, csc->colMatBeg, cupdlp_int, csc->nCols + 1);
  cupdlp_copy(cs_csc->i, csc->colMatIdx, cupdlp_int, csc->nMatElem);
  cupdlp_copy(cs_csc->x, csc->colMatElem, cupdlp_float, csc->nMatElem);
  cupdlp_printf("cs_csc->m = %d, cs_csc->n = %d, cs_csc->nzmax = %d\n", cs_csc->m, cs_csc->n, cs_csc->nzmax);
  cupdlp_dcs *cs_csr = cupdlp_dcs_transpose(cs_csc, 1);
  cupdlp_printf("cs_csr->m = %d, cs_csr->n = %d, cs_csr->nzmax = %d\n", cs_csr->m, cs_csr->n, cs_csr->nzmax);
  csr->nCols = cs_csr->m;
  csr->nRows = cs_csr->n;
  csr->nMatElem = cs_csr->nzmax;
  cupdlp_printf("csr->nRows = %d, csr->nCols = %d, csr->nMatElem = %d\n", csr->nRows, csr->nCols, csr->nMatElem);
  CUPDLP_COPY_VEC(csr->rowMatBeg, cs_csr->p, cupdlp_int, cs_csr->n + 1);
  CUPDLP_COPY_VEC(csr->rowMatIdx, cs_csr->i, cupdlp_int, cs_csr->nzmax);
  CUPDLP_COPY_VEC(csr->rowMatElem, cs_csr->x, cupdlp_float, cs_csr->nzmax);
  // printf("cs_csr->p[0]: %d\n", cs_csr->p[0]);
  // int *host_rowMatBeg = (int *)malloc((cs_csr->n + 1) * sizeof(int)); // 分配 CPU 内存
  // if (host_rowMatBeg != NULL)
  // {
  //   cudaMemcpy(host_rowMatBeg, csr->rowMatBeg, (cs_csr->n + 1) * sizeof(int), cudaMemcpyDeviceToHost); // 从 GPU 复制到 CPU
  //   printf("csr->rowMatBeg[0] = %d\n", host_rowMatBeg[0]);                                             // 在 CPU 上打印
  //   free(host_rowMatBeg);                                                                              // 释放 CPU 内存
  // }
  // else
  // {
  //   printf("Failed to allocate host memory\n");
  // }
  // printf("csr->rowMatBeg[0] = %d\n", csr->rowMatBeg[0]);
  cupdlp_dcs_free(cs_csc);
  cupdlp_dcs_free(cs_csr);

#if !(CUPDLP_CPU)
  // Pointer to GPU csc matrix
  CHECK_CUSPARSE(cusparseCreateCsr(
      &csr->cuda_csr, csr->nRows, csr->nCols, csr->nMatElem, csr->rowMatBeg,
      csr->rowMatIdx, csr->rowMatElem, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
      CUSPARSE_INDEX_BASE_ZERO, CudaComputeType));
#endif

  return 0;
}

void dense2csr(CUPDLPcsr *csr, CUPDLPdense *dense)
{
  csr->nRows = dense->nRows;
  csr->nCols = dense->nCols;

  cupdlp_int nnz = 0;
  cupdlp_int iCol = 0;
  cupdlp_int iRow = 0;
  csr->rowMatBeg[0] = 0;
  for (iRow = 0; iRow < csr->nRows; ++iRow)
  {
    for (iCol = 0; iCol < csr->nCols; ++iCol)
    {
      if (dense->data[iCol * csr->nRows + iRow] != 0)
      {
        csr->rowMatIdx[nnz] = iCol;
        csr->rowMatElem[nnz] = dense->data[iCol * csr->nRows + iRow];
        ++nnz;
      }
    }
    csr->rowMatBeg[iRow + 1] = nnz;
  }
  csr->nMatElem = nnz;

  return;
}

void dense2csc(CUPDLPcsc *csc, CUPDLPdense *dense)
{
  csc->nRows = dense->nRows;
  csc->nCols = dense->nCols;

  cupdlp_int nnz = 0;
  cupdlp_int iCol = 0;
  cupdlp_int iRow = 0;
  csc->colMatBeg[0] = 0;
  for (iCol = 0; iCol < csc->nCols; ++iCol)
  {
    for (iRow = 0; iRow < csc->nRows; ++iRow)
    {
      if (dense->data[iCol * csc->nRows + iRow] != 0)
      {
        csc->colMatIdx[nnz] = iRow;
        csc->colMatElem[nnz] = dense->data[iCol * csc->nRows + iRow];
        ++nnz;
      }
    }
    csc->colMatBeg[iCol + 1] = nnz;
  }

  csc->nMatElem = nnz;

  return;
}

void csr2dense(CUPDLPdense *dense, CUPDLPcsr *csr)
{
  dense->nRows = csr->nRows;
  dense->nCols = csr->nCols;

  cupdlp_int iRow = 0;
  cupdlp_int iCol = 0;
  cupdlp_int iMatElem = 0;
  for (iRow = 0; iRow < dense->nRows; ++iRow)
    for (iCol = 0; iCol < dense->nCols; ++iCol)
    {
      if (iCol == csr->rowMatIdx[iMatElem])
      {
        dense->data[iRow * dense->nCols + iCol] = csr->rowMatElem[iMatElem];
        ++iMatElem;
      }
      else
      {
        dense->data[iRow * dense->nCols + iCol] = 0;
      }
    }

  return;
}

void csc2dense(CUPDLPdense *dense, CUPDLPcsc *csc)
{
  dense->nRows = csc->nRows;
  dense->nCols = csc->nCols;

  cupdlp_int iRow = 0;
  cupdlp_int iCol = 0;
  cupdlp_int iMatElem = 0;
  for (iCol = 0; iCol < dense->nCols; ++iCol)
    for (iRow = 0; iRow < dense->nRows; ++iRow)
    {
      if (iRow == csc->colMatIdx[iMatElem])
      {
        dense->data[iRow * dense->nCols + iCol] = csc->colMatElem[iMatElem];
        ++iMatElem;
      }
      else
      {
        dense->data[iRow * dense->nCols + iCol] = 0;
      }
    }

  return;
}

cupdlp_retcode dense_create(CUPDLPdense **dense)
{
  cupdlp_retcode retcode = RETCODE_OK;
  CUPDLP_INIT(*dense, 1);

exit_cleanup:
  return retcode;
}

cupdlp_retcode csr_create(CUPDLPcsr **csr)
{
  cupdlp_retcode retcode = RETCODE_OK;
  CUPDLP_INIT(*csr, 1);

exit_cleanup:
  return retcode;
}

cupdlp_retcode csc_create(CUPDLPcsc **csc)
{
  cupdlp_retcode retcode = RETCODE_OK;
  CUPDLP_INIT(*csc, 1);

exit_cleanup:
  return retcode;
}

cupdlp_retcode dense_alloc_matrix(CUPDLPdense *dense, cupdlp_int nRows,
                                  cupdlp_int nCols, void *src,
                                  CUPDLP_MATRIX_FORMAT src_matrix_format)
{
  cupdlp_retcode retcode = RETCODE_OK;
  CUPDLP_INIT_ZERO_VEC(dense->data, nRows * nCols);

  switch (src_matrix_format)
  {
  case DENSE:
    dense_copy(dense, (CUPDLPdense *)src);
    break;
  case CSR:
    csr2dense(dense, (CUPDLPcsr *)src);
    break;
  case CSC:
    csc2dense(dense, (CUPDLPcsc *)src);
    break;
  default:
    break;
  }
exit_cleanup:
  return retcode;
}

cupdlp_retcode csr_alloc_matrix(CUPDLPcsr *csr, cupdlp_int nRows,
                                cupdlp_int nCols, void *src,
                                CUPDLP_MATRIX_FORMAT src_matrix_format)
{
  cupdlp_retcode retcode = RETCODE_OK;
  cupdlp_int nnz = 0;
  switch (src_matrix_format)
  {
  case DENSE:
    nnz = nRows * nCols;
    break;
  case CSR:
    nnz = ((CUPDLPcsr *)src)->nMatElem;
    break;
  case CSC:
    nnz = ((CUPDLPcsc *)src)->nMatElem;
    break;
  default:
    break;
  }
  // todo make sure this is right
  CUPDLP_INIT_ZERO_VEC(csr->rowMatBeg, nRows + 1);
  CUPDLP_INIT_ZERO_VEC(csr->rowMatIdx, nnz);
  CUPDLP_INIT_ZERO_VEC(csr->rowMatElem, nnz);
  // printf("csr->rowMatBeg[0] = %d\n", csr->rowMatBeg[0]);
  switch (src_matrix_format)
  {
  case DENSE:
    dense2csr(csr, (CUPDLPdense *)src);
    break;
  case CSR:
    csr_copy(csr, (CUPDLPcsr *)src);
    break;
  case CSC:
    csc2csr(csr, (CUPDLPcsc *)src);
    break;
  default:
    break;
  }
exit_cleanup:
  return retcode;
}

cupdlp_retcode csc_alloc_matrix(CUPDLPcsc *csc, cupdlp_int nRows,
                                cupdlp_int nCols, void *src,
                                CUPDLP_MATRIX_FORMAT src_matrix_format)
{
  cupdlp_retcode retcode = RETCODE_OK;
  cupdlp_int nnz = 0;
  switch (src_matrix_format)
  {
  case DENSE:
    nnz = nRows * nCols;
    break;
  case CSR:
    nnz = ((CUPDLPcsr *)src)->nMatElem;
    break;
  case CSC:
    nnz = ((CUPDLPcsc *)src)->nMatElem;
    break;
  default:
    break;
  }

  CUPDLP_INIT_ZERO_VEC(csc->colMatBeg, nCols + 1);
  CUPDLP_INIT_ZERO_VEC(csc->colMatIdx, nnz);
  CUPDLP_INIT_ZERO_VEC(csc->colMatElem, nnz);
  // cupdlp_printf("csc_alloc_matrix断点2\n");
  switch (src_matrix_format)
  {
  case DENSE:
    dense2csc(csc, (CUPDLPdense *)src);
    break;
  case CSR:
    csr2csc(csc, (CUPDLPcsr *)src);
    break;
  case CSC:
    csc_copy(csc, (CUPDLPcsc *)src);
    break;
  default:
    break;
  }
  // cupdlp_printf("csc_alloc_matrix断点3\n");

exit_cleanup:
  return retcode;
}

cupdlp_retcode dense_alloc(CUPDLPdense *dense, cupdlp_int nRows,
                           cupdlp_int nCols, cupdlp_float *val)
{
  cupdlp_retcode retcode = RETCODE_OK;
  dense->nRows = nRows;
  dense->nCols = nCols;
  dense->data = cupdlp_NULL;
  CUPDLP_INIT_ZERO_VEC(dense->data, nRows * nCols);

  CUPDLP_COPY_VEC(dense->data, val, cupdlp_float, nRows * nCols);
exit_cleanup:
  return retcode;
}

cupdlp_retcode csr_alloc(CUPDLPcsr *csr, cupdlp_int nRows, cupdlp_int nCols,
                         cupdlp_int nnz, cupdlp_int *row_ptr,
                         cupdlp_int *col_ind, cupdlp_float *val)
{
  cupdlp_retcode retcode = RETCODE_OK;
  csr->nRows = nRows;
  csr->nCols = nCols;
  csr->nMatElem = nnz;
  csr->rowMatBeg = cupdlp_NULL;
  csr->rowMatIdx = cupdlp_NULL;
  csr->rowMatElem = cupdlp_NULL;
  CUPDLP_INIT_ZERO_VEC(csr->rowMatBeg, nRows + 1);
  CUPDLP_INIT_ZERO_VEC(csr->rowMatIdx, nnz);
  CUPDLP_INIT_ZERO_VEC(csr->rowMatElem, nnz);

  CUPDLP_COPY_VEC(csr->rowMatBeg, row_ptr, cupdlp_int, nRows + 1);
  CUPDLP_COPY_VEC(csr->rowMatIdx, col_ind, cupdlp_int, nnz);
  CUPDLP_COPY_VEC(csr->rowMatElem, val, cupdlp_float, nnz);
exit_cleanup:
  return retcode;
}

cupdlp_retcode csc_alloc(CUPDLPcsc *csc, cupdlp_int nRows, cupdlp_int nCols,
                         cupdlp_int nnz, cupdlp_int *col_ptr,
                         cupdlp_int *row_ind, cupdlp_float *val)
{
  cupdlp_retcode retcode = RETCODE_OK;
  csc->nRows = nRows;
  csc->nCols = nCols;
  csc->nMatElem = nnz;
  csc->colMatBeg = cupdlp_NULL;
  csc->colMatIdx = cupdlp_NULL;
  csc->colMatElem = cupdlp_NULL;
  CUPDLP_INIT_ZERO_VEC(csc->colMatBeg, nCols + 1);
  CUPDLP_INIT_ZERO_VEC(csc->colMatIdx, nnz);
  CUPDLP_INIT_ZERO_VEC(csc->colMatElem, nnz);

  CUPDLP_COPY_VEC(csc->colMatBeg, col_ptr, cupdlp_int, nCols + 1);
  CUPDLP_COPY_VEC(csc->colMatIdx, row_ind, cupdlp_int, nnz);
  CUPDLP_COPY_VEC(csc->colMatElem, val, cupdlp_float, nnz);
exit_cleanup:
  return retcode;
}

void vecPrint(const char *s, const cupdlp_float *a, cupdlp_int n)
{
  cupdlp_printf("%s: ", s);
  for (cupdlp_int i = 0; i < n; ++i)
  {
    cupdlp_printf("%.3f ", a[i]);
  }
  cupdlp_printf("\n");
}

void vecIntPrint(const char *s, const cupdlp_int *a, cupdlp_int n)
{
  cupdlp_printf("%s: ", s);
  for (cupdlp_int i = 0; i < n; ++i)
  {
    cupdlp_printf("%d ", a[i]);
  }
  cupdlp_printf("\n");
}

void PDHG_Dump_Stats(CUPDLPwork *w)
{
  cupdlp_int nCols = w->iterates->nCols;
  cupdlp_int nRows = w->iterates->nRows;
  CUPDLPiterates *iterates = w->iterates;
  CUPDLPstepsize *stepsize = w->stepsize;

  cupdlp_printf("------------------------------------------------\n");
  cupdlp_printf("Iteration % 3d\n", w->timers->nIter);
#if CUPDLP_DUMP_ITERATES
  vecPrint("x", iterates->x->data, nCols);
  vecPrint("y", iterates->y->data, nRows);
  vecPrint("xSum", iterates->xSum, nCols);
  vecPrint("ySum", iterates->ySum, nRows);
  vecPrint("Ax ", iterates->ax->data, nRows);
  vecPrint("A'y", iterates->aty->data, nCols);
  vecPrint("xLastRestart", iterates->xLastRestart, nCols);
  vecPrint("yLastRestart", iterates->yLastRestart, nRows);
#endif
  cupdlp_printf(
      "PrimalStep: %e, SumPrimalStep: %e, DualStep: %e, SumDualStep: %e\n",
      stepsize->dPrimalStep, stepsize->dSumPrimalStep, stepsize->dDualStep,
      stepsize->dSumDualStep);
  cupdlp_printf("Stepsize: %e, Primal weight: %e Ratio: %e\n",
                sqrt(stepsize->dPrimalStep * stepsize->dDualStep),
                sqrt(stepsize->dBeta), stepsize->dTheta);
}

void PDTEST_Dump_Stats(CUPDLPwork *w)
{
  cupdlp_int nCols = w->PDTESTiterates->nCols;
  cupdlp_int nRows = w->PDTESTiterates->nRows;
  PDTESTiterates *iterates = w->PDTESTiterates;
  CUPDLPstepsize *stepsize = w->stepsize;

  cupdlp_printf("------------------------------------------------\n");
  cupdlp_printf("Iteration % 3d\n", w->timers->nIter);
#if CUPDLP_DUMP_ITERATES
  vecPrint("x", iterates->x->data, nCols);
  vecPrint("y", iterates->y->data, nRows);
  vecPrint("xSum", iterates->xSum, nCols);
  vecPrint("ySum", iterates->ySum, nRows);
  vecPrint("Ax ", iterates->ax->data, nRows);
  vecPrint("A'y", iterates->aty->data, nCols);
  vecPrint("xLastRestart", iterates->xLastRestart, nCols);
  vecPrint("yLastRestart", iterates->yLastRestart, nRows);
#endif
  cupdlp_printf(
      "PrimalStep: %e, SumPrimalStep: %e, DualStep: %e, SumDualStep: %e\n",
      stepsize->dPrimalStep, stepsize->dSumPrimalStep, stepsize->dDualStep,
      stepsize->dSumDualStep);
  cupdlp_printf("Stepsize: %e, Primal weight: %e Ratio: %e\n",
                sqrt(stepsize->dPrimalStep * stepsize->dDualStep),
                sqrt(stepsize->dBeta), stepsize->dTheta);
}

void csrPrintDense(const char *s, CUPDLPcsr *csr)
{
  cupdlp_printf("------------------------------------------------\n");
  cupdlp_printf("%s:\n", s);
  cupdlp_int deltaCol = 0;
  for (cupdlp_int iRow = 0; iRow < csr->nRows; ++iRow)
  {
    for (cupdlp_int iElem = csr->rowMatBeg[iRow];
         iElem < csr->rowMatBeg[iRow + 1]; ++iElem)
    {
      if (iElem == csr->rowMatBeg[iRow])
        deltaCol = csr->rowMatIdx[iElem];
      else
        deltaCol = csr->rowMatIdx[iElem] - csr->rowMatIdx[iElem - 1] - 1;
      for (cupdlp_int i = 0; i < deltaCol; ++i)
      {
        cupdlp_printf("       ");
      }
      cupdlp_printf("%6.3f ", csr->rowMatElem[iElem]);
    }
    cupdlp_printf("\n");
  }
  cupdlp_printf("------------------------------------------------\n");
}

void cscPrintDense(const char *s, CUPDLPcsc *csc)
{
  cupdlp_printf("------------------------------------------------\n");
  cupdlp_printf("%s (Trans):\n", s);
  cupdlp_int deltaRow = 0;
  for (cupdlp_int iCol = 0; iCol < csc->nCols; ++iCol)
  {
    for (cupdlp_int iElem = csc->colMatBeg[iCol];
         iElem < csc->colMatBeg[iCol + 1]; ++iElem)
    {
      if (iElem == csc->colMatBeg[iCol])
        deltaRow = csc->colMatIdx[iElem];
      else
        deltaRow = csc->colMatIdx[iElem] - csc->colMatIdx[iElem - 1] - 1;
      for (cupdlp_int i = 0; i < deltaRow; ++i)
      {
        cupdlp_printf("       ");
      }
      cupdlp_printf("%6.3f ", csc->colMatElem[iElem]);
    }
    cupdlp_printf("\n");
  }
  cupdlp_printf("------------------------------------------------\n");
}

#ifndef CUPDLP_OUTPUT_NAMES
const char *termCodeNames[] = {
    "OPTIMAL",
    "INFEASIBLE",
    "UNBOUNDED",
    "INFEASIBLE_OR_UNBOUNDED",
    "TIMELIMIT_OR_ITERLIMIT",
};
const char *termIterateNames[] = {
    "LAST_ITERATE",
    "AVERAGE_ITERATE",
};
#endif

void writeJson(const char *fout, CUPDLPwork *work, cupdlp_float *x,
               cupdlp_int nx, cupdlp_float *y, cupdlp_int ny,
               cupdlp_bool ifSaveSol)
{
  FILE *fptr;

  cupdlp_printf("--------------------------------\n");
  cupdlp_printf("--- saving to %s\n", fout);
  cupdlp_printf("--------------------------------\n");
  // Open a file in writing mode
  fptr = fopen(fout, "w");

  fprintf(fptr, "{\n");

  // timers
  fprintf(fptr, "%d, ", work->timers->nIter);
  fprintf(fptr, "%.2fs\n", work->timers->dSolvingTime);
  fprintf(fptr, "\"dPrimalObj\":%.10f,", work->resobj->dPrimalObj);
  fprintf(fptr, "\"nIter\":%d,", work->timers->nIter);
  fprintf(fptr, "\"nAtyCalls\":%d,", work->timers->nAtyCalls);
  fprintf(fptr, "\"nAxCalls\":%d,", work->timers->nAxCalls);
  fprintf(fptr, "\"dSolvingBeg\":%f,", work->timers->dSolvingBeg);
  fprintf(fptr, "\"dSolvingTime\":%f,", work->timers->dSolvingTime);
  fprintf(fptr, "\"dPresolveTime\":%f,", work->timers->dPresolveTime);
  fprintf(fptr, "\"dScalingTime\":%f,", work->timers->dScalingTime);
#if !(CUPDLP_CPU)
  fprintf(fptr, "\"AllocMem_CopyMatToDeviceTime\":%f,",
          work->timers->AllocMem_CopyMatToDeviceTime);
  fprintf(fptr, "\"CopyVecToDeviceTime\":%f,",
          work->timers->CopyVecToDeviceTime);
  fprintf(fptr, "\"CopyVecToHostTime\":%f,", work->timers->CopyVecToHostTime);
  fprintf(fptr, "\"DeviceMatVecProdTime\":%f,",
          work->timers->DeviceMatVecProdTime);
#endif
  // residuals
  fprintf(fptr, "\"dPrimalObj\":%f,", work->resobj->dPrimalObj);
  fprintf(fptr, "\"dDualObj\":%f,", work->resobj->dDualObj);
  fprintf(fptr, "\"dPrimalFeas\":%f,", work->resobj->dPrimalFeas);
  fprintf(fptr, "\"dDualFeas\":%f,", work->resobj->dDualFeas);
  fprintf(fptr, "\"dPrimalObjAverage\":%f,", work->resobj->dPrimalObjAverage);
  fprintf(fptr, "\"dDualObjAverage\":%f,", work->resobj->dDualObjAverage);
  fprintf(fptr, "\"dPrimalFeasAverage\":%f,", work->resobj->dPrimalFeasAverage);
  fprintf(fptr, "\"dDualFeasAverage\":%f,", work->resobj->dDualFeasAverage);
  fprintf(fptr, "\"dDualityGap\":%f,", work->resobj->dDualityGap);
  fprintf(fptr, "\"dDualityGapAverage\":%f,", work->resobj->dDualityGapAverage);
  fprintf(fptr, "\"dComplementarity\":%f,", work->resobj->dComplementarity);
  fprintf(fptr, "\"dComplementarityAverage\":%f,",
          work->resobj->dComplementarityAverage);

  //  todo should this be added to postsolve?
  // todo, fix dNormCost and this
  if (work->resobj->termIterate == AVERAGE_ITERATE)
  {
    fprintf(fptr, "\"dRelPrimalFeas\":%f,",
            work->resobj->dPrimalFeasAverage / (1.0 + work->scaling->dNormRhs));
    fprintf(fptr, "\"dRelDualFeas\":%f,",
            work->resobj->dDualFeasAverage / (1.0 + work->scaling->dNormCost));
    fprintf(fptr, "\"dRelDualityGap\":%f,", work->resobj->dRelObjGapAverage);
  }
  else
  {
    fprintf(fptr, "\"dRelPrimalFeas\":%f,",
            work->resobj->dPrimalFeas / (1.0 + work->scaling->dNormRhs));
    fprintf(fptr, "\"dRelDualFeas\":%f,",
            work->resobj->dDualFeas / (1.0 + work->scaling->dNormCost));
    fprintf(fptr, "\"dRelDualityGap\":%f,", work->resobj->dRelObjGap);
  }
  fprintf(fptr, "\"terminationCode\":\"%s\",",
          termCodeNames[work->resobj->termCode]);
  fprintf(fptr, "\"terminationIterate\":\"%s\"",
          termIterateNames[work->resobj->termIterate]);

  // print solutions
  if (ifSaveSol)
  {
    // primal solution
    fprintf(fptr, ",\"x\":[");
    if (x && nx > 0)
    {
      for (int i = 0; i < nx - 1; ++i)
      {
        fprintf(fptr, "%f,", x[i]);
      }
      fprintf(fptr, "%f", x[nx - 1]);
    }
    fprintf(fptr, "]");

    // dual solution
    fprintf(fptr, ",\"y\":[");
    if (y && ny > 0)
    {
      for (int i = 0; i < ny - 1; ++i)
      {
        fprintf(fptr, "%f,", y[i]);
      }
      fprintf(fptr, "%f", y[ny - 1]);
    }
    fprintf(fptr, "]");
  }

  fprintf(fptr, "}");
  // Close the file
  fclose(fptr);
}
