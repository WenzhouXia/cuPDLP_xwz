#include "mps_lp.h"
#include "wrapper_clp.h"

/*
    min cTx
    s.t. Aeq x = b
         Aineq x <= bineq
         ub >= x >= 0
         colUbIdx: index of columns with upper bound (not all columns have upper
   bound)
*/
cupdlp_retcode main(int argc, char **argv)
{
  cupdlp_retcode retcode = RETCODE_OK;

  char *fname = "./example/afiro.mps.gz";
  char *fout = "./solu·on.json";

  int nCols;
  int nRows;
  int nEqs;
  int nCols_origin;
  cupdlp_bool ifSaveSol = false;
  cupdlp_bool ifPresolve = false;
  cupdlp_int ifPDTEST = false;
  int nnz = 0;
  double *rhs = NULL;
  double *cost = NULL;

  cupdlp_float *lower = NULL;
  cupdlp_float *upper = NULL;

  // -------------------------
  int *csc_beg = NULL, *csc_idx = NULL;
  double *csc_val = NULL;
  double offset =
      0.0;                // true objVal = sig * c'x - offset, sig = 1 (min) or -1 (max)
  double sign_origin = 1; // 1 (min) or -1 (max)
  int *constraint_new_idx = NULL;
  cupdlp_float *x_origin = cupdlp_NULL;
  cupdlp_float *y_origin = cupdlp_NULL;

  void *model = NULL;
  void *presolvedmodel = NULL;
  void *presolveinfo = NULL;

  CUPDLPscaling *scaling =
      (CUPDLPscaling *)cupdlp_malloc(sizeof(CUPDLPscaling));

  // claim solvers variables
  // prepare pointers
  CUPDLP_MATRIX_FORMAT src_matrix_format = CSC;     // 原矩阵格式
  CUPDLP_MATRIX_FORMAT dst_matrix_format = CSR_CSC; // 目标矩阵格式
  CUPDLPcsc *csc_cpu = cupdlp_NULL;
  CUPDLPproblem *prob = cupdlp_NULL;

  // load parameters
  for (cupdlp_int i = 0; i < argc - 1; i++)
  {
    // if (strcmp(argv[i], "-niter") == 0) {
    //   niters = atof(argv[i + 1]);
    // } else
    if (strcmp(argv[i], "-fname") == 0)
    {
      fname = argv[i + 1];
    }
    else if (strcmp(argv[i], "-out") == 0)
    {
      fout = argv[i + 1];
    }
    else if (strcmp(argv[i], "-h") == 0)
    {
      print_script_usage();
    }
    else if (strcmp(argv[i], "-savesol") == 0)
    {
      ifSaveSol = atoi(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-ifPre") == 0)
    {
      ifPresolve = atoi(argv[i + 1]);
    }
    else if (strcmp(argv[i], "-ifPDTEST") == 0)
    {
      ifPDTEST = atoi(argv[i + 1]);
    }
  }
  if (strcmp(argv[argc - 1], "-h") == 0)
  {
    print_script_usage();
  }

  // set solver parameters
  cupdlp_bool ifChangeIntParam[N_INT_USER_PARAM] = {false};
  cupdlp_int intParam[N_INT_USER_PARAM] = {0};
  cupdlp_bool ifChangeFloatParam[N_FLOAT_USER_PARAM] = {false};
  cupdlp_float floatParam[N_FLOAT_USER_PARAM] = {0.0};
  CUPDLP_CALL(getUserParam(argc, argv, ifChangeIntParam, intParam,
                           ifChangeFloatParam, floatParam));

  model = createModel();
  loadMps(model, fname);

  void *model2solve = model;

  if (ifChangeIntParam[IF_PRESOLVE])
  {
    ifPresolve = intParam[IF_PRESOLVE];
  }

  cupdlp_float presolve_time = getTimeStamp();
  if (ifPresolve)
  {
    presolveinfo = createPresolve();
    presolvedmodel = presolvedModel(presolveinfo, model);
    model2solve = presolvedmodel;
  }
  presolve_time = getTimeStamp() - presolve_time;

  // CUPDLP_CALL(formulateLP(model, &cost, &nCols, &nRows, &nnz, &nEqs,
  // &csc_beg,
  //                         &csc_idx, &csc_val, &rhs, &lower, &upper,
  //                         &offset));

  // CUPDLP_CALL(formulateLP_new(
  //     model, &cost, &nCols, &nRows, &nnz, &nEqs, &csc_beg, &csc_idx,
  //     &csc_val, &rhs, &lower, &upper, &offset, &nCols_origin,
  //     &constraint_new_idx));

  CUPDLP_CALL(formulateLP_new(model2solve, &cost, &nCols, &nRows, &nnz, &nEqs,
                              &csc_beg, &csc_idx, &csc_val, &rhs, &lower,
                              &upper, &offset, &sign_origin, &nCols_origin,
                              &constraint_new_idx));

  /*
      min cTx
      s.t. Aeq x = b
           Aineq x <= bineq
           ub >= x >= 0
           colUbIdx: index of columns with upper bound (not all columns have
     upper bound)
  */

  if (retcode != RETCODE_OK)
  {
    cupdlp_printf("Error reading MPS file\n");
    retcode = RETCODE_FAILED;
    goto exit_cleanup;
  }

  CUPDLP_CALL(Init_Scaling(scaling, nCols, nRows, cost, rhs));
  cupdlp_int ifScaling = 1;

  if (ifChangeIntParam[IF_SCALING])
  {
    ifScaling = intParam[IF_SCALING];
  }

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

  // these two handles need to be established first
  CUPDLPwork *w = cupdlp_NULL;
  CUPDLP_INIT_ZERO(w, 1);
#if !(CUPDLP_CPU)
  cupdlp_float cuda_prepare_time = getTimeStamp();
  CHECK_CUSPARSE(cusparseCreate(&w->cusparsehandle));
  CHECK_CUBLAS(cublasCreate(&w->cublashandle));
  cuda_prepare_time = getTimeStamp() - cuda_prepare_time;
#endif

  CUPDLP_CALL(problem_create(&prob));

  // currently, only supprot that input matrix is CSC, and store both CSC and
  // CSR
  CUPDLP_CALL(csc_create(&csc_cpu));
  csc_cpu->nRows = nRows;
  csc_cpu->nCols = nCols;
  csc_cpu->nMatElem = nnz;
  csc_cpu->colMatBeg = (int *)malloc((1 + nCols) * sizeof(int));
  csc_cpu->colMatIdx = (int *)malloc(nnz * sizeof(int));
  csc_cpu->colMatElem = (double *)malloc(nnz * sizeof(double));
  memcpy(csc_cpu->colMatBeg, csc_beg, (nCols + 1) * sizeof(int));
  memcpy(csc_cpu->colMatIdx, csc_idx, nnz * sizeof(int));
  memcpy(csc_cpu->colMatElem, csc_val, nnz * sizeof(double));
#if !(CUPDLP_CPU)
  csc_cpu->cuda_csc = NULL;
#endif

  cupdlp_float scaling_time = getTimeStamp();
  CUPDLP_CALL(PDHG_Scale_Data_cuda(csc_cpu, ifScaling, scaling, cost, lower,
                                   upper, rhs));
  scaling_time = getTimeStamp() - scaling_time;

  cupdlp_float alloc_matrix_time = 0.0;
  cupdlp_float copy_vec_time = 0.0;

  if (ifChangeIntParam[IF_PDTEST])
  {
    ifPDTEST = intParam[IF_PDTEST]; // 默认用PDHG不用PDTEST，ifPDTEST默认值为0
  }
  if (!ifPDTEST)
  {
    CUPDLP_CALL(problem_alloc(prob, nRows, nCols, nEqs, cost, offset, sign_origin,
                              csc_cpu, src_matrix_format, dst_matrix_format, rhs,
                              lower, upper, &alloc_matrix_time, &copy_vec_time));
  }
  else
  {
    /////////////////////////////////////////////////////////////////////
    cupdlp_printf("--------------------------------------------------\n");
    cupdlp_float sum = 0.0;
    cupdlp_printf("csc_cpu->nMatElem = %d\n", csc_cpu->nMatElem);
    for (cupdlp_int i = 0; i < csc_cpu->nMatElem; ++i)
    {
      cupdlp_float value = csc_cpu->colMatElem[i]; // 获取矩阵非零元素的值
      sum += value * value;                        // 将元素的平方累加到总和中
    }
    cupdlp_float matrix_2norm = sqrt(sum);
    cupdlp_printf("matrix_2norm = %f\n", matrix_2norm);
    /////////////////////////////////////////////////////////////////////
    CUPDLP_CALL(PDTEST_problem_alloc(prob, nRows, nCols, nEqs, cost, offset, sign_origin,
                                     csc_cpu, src_matrix_format, dst_matrix_format, rhs,
                                     lower, upper, &alloc_matrix_time, &copy_vec_time, matrix_2norm));
  }

  // solve
  // cupdlp_printf("Enter main solve loop\n");

  w->problem = prob;
  w->scaling = scaling;
  if (!ifPDTEST)
  {
    PDHG_Alloc(w);
  }
  else
  {
    PDTEST_Alloc(w);
  }
  // // PDHG_Alloc(w);
  // PDTEST_Alloc(w);
  w->timers->dScalingTime = scaling_time;
  w->timers->dPresolveTime = presolve_time;
  CUPDLP_COPY_VEC(w->rowScale, scaling->rowScale, cupdlp_float, nRows);
  CUPDLP_COPY_VEC(w->colScale, scaling->colScale, cupdlp_float, nCols);

#if !(CUPDLP_CPU)
  w->timers->AllocMem_CopyMatToDeviceTime += alloc_matrix_time;
  w->timers->CopyVecToDeviceTime += copy_vec_time;
  w->timers->CudaPrepareTime = cuda_prepare_time;
#endif

  // CUPDLP_CALL(LP_SolvePDHG(prob, cupdlp_NULL, cupdlp_NULL, cupdlp_NULL,
  // cupdlp_NULL));
  //   CUPDLP_CALL(LP_SolvePDHG(prob, ifChangeIntParam, intParam,
  //                               ifChangeFloatParam, floatParam, fout));

  CUPDLP_INIT(x_origin, nCols_origin);
  CUPDLP_INIT(y_origin, nRows);
  switch (ifPDTEST)
  {
  case 0:
    cupdlp_printf("--------------------------------------------------\n");
    cupdlp_printf("enter main solve loop, PDHG\n");
    cupdlp_printf("--------------------------------------------------\n");
    CUPDLP_CALL(LP_SolvePDHG(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx));
    break;
  case 1:
    cupdlp_printf("--------------------------------------------------\n");
    cupdlp_printf("enter main solve loop, PDTEST\n");
    cupdlp_printf("--------------------------------------------------\n");
    CUPDLP_CALL(LP_SolvePDTEST(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx));
    break;
  case 2:
    cupdlp_printf("--------------------------------------------------\n");
    cupdlp_printf("enter main solve loop, PDTESTm_min\n");
    cupdlp_printf("--------------------------------------------------\n");
    CUPDLP_CALL(LP_SolvePDTEST_min(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx));
    break;
  default:
    cupdlp_printf("Error: ifPDTEST = %d\n, 不在取值范围内", ifPDTEST);
    break;
  }
  //

  // print result
  // TODO: implement after adding IO

exit_cleanup:
  deleteModel(model);
  if (ifPresolve)
  {
    deletePresolve(presolveinfo);
    deleteModel(presolvedmodel);
  }

  if (scaling)
  {
    scaling_clear(scaling);
  }
  if (cost != NULL)
    cupdlp_free(cost);
  if (csc_beg != NULL)
    cupdlp_free(csc_beg);
  if (csc_idx != NULL)
    cupdlp_free(csc_idx);
  if (csc_val != NULL)
    cupdlp_free(csc_val);
  if (rhs != NULL)
    cupdlp_free(rhs);
  if (lower != NULL)
    cupdlp_free(lower);
  if (upper != NULL)
    cupdlp_free(upper);
  if (constraint_new_idx != NULL)
    cupdlp_free(constraint_new_idx);
  if (x_origin != NULL)
    cupdlp_free(x_origin);
  if (y_origin != NULL)
    cupdlp_free(y_origin);
  // free memory
  csc_clear(csc_cpu);
  problem_clear(prob);

  // freealldata(Aeqp, Aeqi, Aeqx, Aineqp, Aineqi, Aineqx, colUbIdx, colUbElem,
  //             rhs, cost, x, s, t, sx, ss, st, y, lower, upper);

  return retcode;
}
