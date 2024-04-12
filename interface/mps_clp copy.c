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
  char *fout = "./solution.json";

  double *cost = NULL;
  int nCols;
  int nRows;
  int nnz = 0;
  int nEqs;
  int *csc_beg = NULL, *csc_idx = NULL;
  double *csc_val = NULL;
  double *rhs = NULL;
  cupdlp_float *lower = NULL;
  cupdlp_float *upper = NULL;
  double offset = 0.0;    // true objVal = sig * c'x - offset, sig = 1 (min) or -1 (max)
  double sign_origin = 1; // 1 (min) or -1 (max)
  int nCols_origin;

  cupdlp_bool ifSaveSol = false;
  cupdlp_bool ifPresolve = false;
  cupdlp_int ifPDTEST = 0;
  cupdlp_int bestID = 1;

  // -------------------------

  // cupdlp_float *x_origin = cupdlp_NULL;
  // cupdlp_float *y_origin = cupdlp_NULL;

  void *presolvedmodel = NULL;
  void *presolveinfo = NULL;

  // claim solvers variables
  // prepare pointers
  CUPDLP_MATRIX_FORMAT src_matrix_format = CSC;     // 原矩阵格式
  CUPDLP_MATRIX_FORMAT dst_matrix_format = CSR_CSC; // 目标矩阵格式
  CUPDLPcsc *csc_cpu = cupdlp_NULL;
  CUPDLPproblem *prob = cupdlp_NULL;

#pragma region load parameters
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
    else if (strcmp(argv[i], "-bestID") == 0)
    {
      bestID = atoi(argv[i + 1]);
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
#pragma endregion

  char basePath[] = "/home/xiawenzhou/Documents/XiaWenzhou/OptimalTransport/OT_instances/DOTmark/Data";
  char type[] = "CauchyDensity";
  int resolution = 32;
  int fileNumber_1 = 1001;
  int fileNumber_2 = 1002;
  char csvpath_1[256];
  char csvpath_2[256];
  sprintf(csvpath_1, "%s/%s/data%d_%d.csv", basePath, type, resolution, fileNumber_1);
  sprintf(csvpath_2, "%s/%s/data%d_%d.csv", basePath, type, resolution, fileNumber_2);
  // void *model = NULL;
  // model = createModel();
  // void *model_coarse = NULL;
  // model_coarse = createModel();
  // loadMps(model, fname);
  // writeLpWrapper(model, "test_ot_load");
  // // coarse_degree = 0, 精确问题
  // // 粗糙问题
  // cupdlp_int coarse_degree = 1;
  // generate_coarse_dualOT_model_from_csv(model_coarse, csvpath_1, csvpath_2, resolution, coarse_degree);
  ///////////////////////////////////////////////////////////////////////////
#pragma region 创建CUPDLPwork，但是被我注释掉了

  // void *model2solve = model;

  // if (ifChangeIntParam[IF_PRESOLVE])
  // {
  //   ifPresolve = intParam[IF_PRESOLVE];
  // }

  // cupdlp_float presolve_time = getTimeStamp();
  // if (ifPresolve)
  // {
  //   presolveinfo = createPresolve();
  //   presolvedmodel = presolvedModel(presolveinfo, model);
  //   model2solve = presolvedmodel;
  // }
  // presolve_time = getTimeStamp() - presolve_time;

  // CUPDLP_CALL(formulateLP_new(model2solve, &cost, &nCols, &nRows, &nnz, &nEqs,
  //                             &csc_beg, &csc_idx, &csc_val, &rhs, &lower,
  //                             &upper, &offset, &sign_origin, &nCols_origin,
  //                             &constraint_new_idx));

  /*
      min cTx
      s.t. Aeq x = b
           Aineq x <= bineq
           ub >= x >= 0
           colUbIdx: index of columns with upper bound (not all columns have
     upper bound)
  */

  // if (retcode != RETCODE_OK)
  // {
  //   cupdlp_printf("Error reading MPS file\n");
  //   retcode = RETCODE_FAILED;
  //   goto exit_cleanup;
  // }

  // CUPDLP_CALL(Init_Scaling(scaling, nCols, nRows, cost, rhs));
  // cupdlp_int ifScaling = 1;

  // if (ifChangeIntParam[IF_SCALING])
  // {
  //   ifScaling = intParam[IF_SCALING];
  // }

  // if (ifChangeIntParam[IF_RUIZ_SCALING])
  // {
  //   scaling->ifRuizScaling = intParam[IF_RUIZ_SCALING];
  // }

  // if (ifChangeIntParam[IF_L2_SCALING])
  // {
  //   scaling->ifL2Scaling = intParam[IF_L2_SCALING];
  // }

  // if (ifChangeIntParam[IF_PC_SCALING])
  // {
  //   scaling->ifPcScaling = intParam[IF_PC_SCALING];
  // }`

  ////////////////////////////////////////////////////////////////////////////////////////
  //   // these two handles need to be established first
  //   CUPDLPwork *w = cupdlp_NULL;
  //   CUPDLP_INIT_ZERO(w, 1);
  // #if !(CUPDLP_CPU)
  //   cupdlp_float cuda_prepare_time = getTimeStamp();
  //   CHECK_CUSPARSE(cusparseCreate(&w->cusparsehandle));
  //   CHECK_CUBLAS(cublasCreate(&w->cublashandle));
  //   cuda_prepare_time = getTimeStamp() - cuda_prepare_time;
  // #endif

  //   CUPDLP_CALL(problem_create(&prob));

  //   // currently, only supprot that input matrix is CSC, and store both CSC and
  //   // CSR
  //   CUPDLP_CALL(csc_create(&csc_cpu));
  //   csc_cpu->nRows = nRows;
  //   csc_cpu->nCols = nCols;
  //   csc_cpu->nMatElem = nnz;
  //   csc_cpu->colMatBeg = (int *)malloc((1 + nCols) * sizeof(int));
  //   csc_cpu->colMatIdx = (int *)malloc(nnz * sizeof(int));
  //   csc_cpu->colMatElem = (double *)malloc(nnz * sizeof(double));
  //   memcpy(csc_cpu->colMatBeg, csc_beg, (nCols + 1) * sizeof(int));
  //   memcpy(csc_cpu->colMatIdx, csc_idx, nnz * sizeof(int));
  //   memcpy(csc_cpu->colMatElem, csc_val, nnz * sizeof(double));
  // #if !(CUPDLP_CPU)
  //   csc_cpu->cuda_csc = NULL;
  // #endif

  //   cupdlp_float scaling_time = getTimeStamp();
  //   CUPDLP_CALL(PDHG_Scale_Data_cuda(csc_cpu, ifScaling, scaling, cost, lower,
  //                                    upper, rhs));
  //   scaling_time = getTimeStamp() - scaling_time;

  //   cupdlp_float alloc_matrix_time = 0.0;
  //   cupdlp_float copy_vec_time = 0.0;

  //     switch (ifPDTEST)
  //   {
  //   case 0:
  //   case 5:
  //     CUPDLP_CALL(problem_alloc(prob, nRows, nCols, nEqs, cost, offset, sign_origin, csc_cpu, src_matrix_format, dst_matrix_format, rhs, lower, upper, &alloc_matrix_time, &copy_vec_time));
  //     break;
  //   case 1:
  //   case 2:
  //   case 3:
  //   case 4:
  //     ///////////////////////////////////////////////////////////////////
  //     cupdlp_printf("--------------------------------------------------\n");
  //     cupdlp_float sum = 0.0;
  //     cupdlp_printf("csc_cpu->nMatElem = %d\n", csc_cpu->nMatElem);
  //     for (cupdlp_int i = 0; i < csc_cpu->nMatElem; ++i)
  //     {
  //       cupdlp_float value = csc_cpu->colMatElem[i]; // 获取矩阵非零元素的值
  //       sum += value * value;                        // 将元素的平方累加到总和中
  //     }
  //     cupdlp_float matrix_2norm = sqrt(sum);
  //     cupdlp_printf("matrix_2norm = %f\n", matrix_2norm);
  //     /////////////////////////////////////////////////////////////////////
  //     CUPDLP_CALL(PDTEST_problem_alloc(prob, nRows, nCols, nEqs, cost, offset, sign_origin,
  //                                      csc_cpu, src_matrix_format, dst_matrix_format, rhs,
  //                                      lower, upper, &alloc_matrix_time, &copy_vec_time, matrix_2norm));
  //     break;
  //   default:
  //     cupdlp_printf("Error: ifPDTEST = %d, 不在取值范围内", ifPDTEST);
  //     break;
  //   }

  //   // solve
  //   // cupdlp_printf("Enter main solve loop\n");

  //   w->problem = prob;
  //   w->scaling = scaling;
  //   switch (ifPDTEST)
  //   {
  //   case 0:
  //   case 5:
  //     PDHG_Alloc(w);
  //     break;
  //   case 1:
  //   case 2:
  //   case 3:
  //   case 4:
  //     PDTEST_Alloc(w);
  //     break;
  //   default:
  //     cupdlp_printf("Error: ifPDTEST = %d, 不在取值范围内", ifPDTEST);
  //     break;
  //   }
  //   // // PDHG_Alloc(w);
  //   // PDTEST_Alloc(w);
  //   w->timers->dScalingTime = scaling_time;
  //   w->timers->dPresolveTime = presolve_time;
  //   CUPDLP_COPY_VEC(w->rowScale, scaling->rowScale, cupdlp_float, nRows);
  //   CUPDLP_COPY_VEC(w->colScale, scaling->colScale, cupdlp_float, nCols);

  // #if !(CUPDLP_CPU)
  //   w->timers->AllocMem_CopyMatToDeviceTime += alloc_matrix_time;
  //   w->timers->CopyVecToDeviceTime += copy_vec_time;
  //   w->timers->CudaPrepareTime = cuda_prepare_time;
  // #endif

  //   // CUPDLP_CALL(LP_SolvePDHG(prob, cupdlp_NULL, cupdlp_NULL, cupdlp_NULL,
  //   // cupdlp_NULL));
  //   //   CUPDLP_CALL(LP_SolvePDHG(prob, ifChangeIntParam, intParam,
  //   //                               ifChangeFloatParam, floatParam, fout));
  // CUPDLP_INIT(x_origin, nCols_origin);
  // CUPDLP_INIT(y_origin, nRows);
  // nCols_origin = 2 * pow(resolution, 2);
  // nRows = pow(resolution, 4);
  // CUPDLP_INIT(x_origin, nCols_origin);
  // CUPDLP_INIT(y_origin, nRows);
////////////////////////////////////////////////////////////////////////////////////////
#pragma endregion

  if (ifChangeIntParam[IF_PDTEST])
  {
    ifPDTEST = intParam[IF_PDTEST]; // 默认用PDHG不用PDTEST，ifPDTEST默认值为0
  }

  // 第2层级
  void *model_2 = NULL;
  model_2 = createModel();
  cupdlp_int coarse_degree_3 = -1;
  cupdlp_int coarse_degree_2 = 2;
  cupdlp_int resolution_coarse_2 = resolution / pow(2, coarse_degree_2);
  cupdlp_float *x_solution_2 = cupdlp_NULL;
  CUPDLP_INIT(x_solution_2, 2 * pow(resolution_coarse_2, 2));
  cupdlp_float *y_solution_2 = cupdlp_NULL;
  CUPDLP_INIT(y_solution_2, pow(resolution_coarse_2, 4));
  cupdlp_float *x_init_2 = cupdlp_NULL;
  CUPDLP_INIT(x_init_2, 2 * pow(resolution_coarse_2, 2));
  cupdlp_float *y_init_2 = cupdlp_NULL;
  CUPDLP_INIT(y_init_2, pow(resolution_coarse_2, 4));
  construct_and_solve_Multiscale(model_2, csvpath_1, csvpath_2, resolution, ifChangeIntParam, ifChangeFloatParam, intParam, floatParam, ifSaveSol, coarse_degree_2, coarse_degree_3, &x_solution_2, &y_solution_2, x_init_2, y_init_2);
  analyseArray1D(y_solution_2, pow(resolution_coarse_2, 4), 1e-20, "./y_solution.txt");

  // 第1层级
  void *model_1 = NULL;
  model_1 = createModel();
  cupdlp_int coarse_degree_1 = 1;
  cupdlp_int resolution_coarse_1 = resolution / pow(2, coarse_degree_1);
  cupdlp_float *x_solution_1 = cupdlp_NULL;
  CUPDLP_INIT(x_solution_1, 2 * pow(resolution_coarse_1, 2));
  cupdlp_float *y_solution_1 = cupdlp_NULL;
  CUPDLP_INIT(y_solution_1, pow(resolution_coarse_1, 4));
  construct_and_solve_Multiscale(model_1, csvpath_1, csvpath_2, resolution, ifChangeIntParam, ifChangeFloatParam, intParam, floatParam, ifSaveSol, coarse_degree_1, coarse_degree_2, &x_solution_1, &y_solution_1, x_solution_2, y_solution_2);
  analyseArray1D(y_solution_2, pow(resolution_coarse_2, 4), 1e-20, "./y_solution.txt");
  analyseArray1D(y_solution_1, pow(resolution_coarse_1, 4), 1e-20, "./y_solution.txt");

  // 第0层级
  // void *model_0 = NULL;
  // model_0 = createModel();
  // cupdlp_int coarse_degree_0 = 0;
  // cupdlp_int resolution_coarse_0 = resolution / pow(2, coarse_degree_0);
  // cupdlp_float *x_solution_0 = cupdlp_NULL;
  // CUPDLP_INIT(x_solution_0, 2 * pow(resolution_coarse_0, 2));
  // cupdlp_float *y_solution_0 = cupdlp_NULL;
  // CUPDLP_INIT(y_solution_0, pow(resolution_coarse_0, 4));
  // construct_and_solve_Multiscale(model_0, csvpath_1, csvpath_2, resolution, ifChangeIntParam, ifChangeFloatParam, intParam, floatParam, ifSaveSol, coarse_degree_0, coarse_degree_1, &x_solution_0, &y_solution_0, x_solution_1, y_solution_1);
  // analyseArray1D(y_solution_2, pow(resolution_coarse_2, 4), 1e-20, "./y_solution.txt");
  // analyseArray1D(y_solution_1, pow(resolution_coarse_1, 4), 1e-20, "./y_solution.txt");
  // analyseArray1D(y_solution_0, pow(resolution_coarse_0, 4), 1e-20, "./y_solution.txt");

  // // 释放内存
  // deleteModel(model_0);
  // cupdlp_free(x_solution_0);
  // cupdlp_free(y_solution_0);
  // deleteModel(model_1);
  // cupdlp_free(x_solution_1);
  // cupdlp_free(y_solution_1);
  // // deleteModel(model_2);
  // cupdlp_free(x_solution_2);
  // cupdlp_free(y_solution_2);

#pragma region 手动Multiscale
//   int *nCols_origin_ptr = cupdlp_NULL;
//   int *nRows_ptr = cupdlp_NULL;
//   int *nCols_origin_ptr_coarse = cupdlp_NULL;
//   int *nRows_ptr_coarse = cupdlp_NULL;
//   CUPDLP_INIT_ZERO(nCols_origin_ptr, 1);
//   CUPDLP_INIT_ZERO(nRows_ptr, 1);
//   CUPDLP_INIT_ZERO(nCols_origin_ptr_coarse, 1);
//   CUPDLP_INIT_ZERO(nRows_ptr_coarse, 1);

//   // 定义了一个指针，想要在createCUPDLPwork中改变这个指针，就要传入指针的地址

//   int *constraint_new_idx_coarse = NULL;

//   CUPDLPwork *w_coarse = cupdlp_NULL;
//   CUPDLP_INIT_ZERO(w_coarse, 1);
//   createCUPDLPwork(w_coarse, model_coarse, ifChangeIntParam, intParam, &constraint_new_idx_coarse, nCols_origin_ptr_coarse, nRows_ptr_coarse);

//   // 接收最优解
//   cupdlp_float *x_coarse_solution = cupdlp_NULL;
//   cupdlp_float *y_coarse_solution = cupdlp_NULL;
//   CUPDLP_INIT_ZERO(x_coarse_solution, *nCols_origin_ptr_coarse);
//   CUPDLP_INIT_ZERO(y_coarse_solution, *nRows_ptr_coarse);
//   // 接受粗糙part最后的步长和权重
//   cupdlp_float *stepsize_coarse_last = cupdlp_NULL;
//   cupdlp_float *weight_coarse_last = cupdlp_NULL;
//   CUPDLP_INIT_ZERO(stepsize_coarse_last, 1);
//   CUPDLP_INIT_ZERO(weight_coarse_last, 1);
//   cupdlp_float *stepsize_coarse_init = cupdlp_NULL;
//   cupdlp_float *weight_coarse_init = cupdlp_NULL;
//   CUPDLP_INIT_ZERO(stepsize_coarse_init, 1);
//   CUPDLP_INIT_ZERO(weight_coarse_init, 1);
//   // 作为初始解
//   cupdlp_float *x_coarse_init = cupdlp_NULL;
//   cupdlp_float *y_coarse_init = cupdlp_NULL;
//   cupdlp_int nCols_coarse = w_coarse->problem->nCols;
//   cupdlp_int nRows_coarse = w_coarse->problem->nRows;
//   CUPDLP_INIT_ZERO(x_coarse_init, w_coarse->problem->nCols);
//   CUPDLP_INIT_ZERO(y_coarse_init, w_coarse->problem->nRows);

//   cupdlp_printf("--------------------------------------------------\n");
//   cupdlp_printf("enter main solve loop, PDHG_Multiscale\n");
//   cupdlp_printf("--------------------------------------------------\n");
//   char *fout_coarse = "./coarse_solution.txt";
//   cupdlp_bool whether_first_coarse = true;
//   CUPDLP_CALL(LP_SolvePDHG_Multiscale(w_coarse, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout_coarse, ifSaveSol, constraint_new_idx_coarse, &x_coarse_solution, &y_coarse_solution, x_coarse_init, y_coarse_init, stepsize_coarse_last, weight_coarse_last, stepsize_coarse_init, weight_coarse_init, whether_first_coarse));

//   cupdlp_free(x_coarse_init);
//   cupdlp_free(y_coarse_init);
//   cupdlp_free(stepsize_coarse_init);
//   cupdlp_free(weight_coarse_init);
// #pragma region 考虑稀疏性的版本
//   cupdlp_printf("--------------------------------------------------\n");
//   cupdlp_printf("--------------------------------------------------\n");

//   int x_init_len = 2 * pow(resolution, 2);
//   int y_init_len = pow(resolution, 4);
//   // coarse的最优解作为fine的初始解
//   cupdlp_float *x_init = cupdlp_NULL;
//   cupdlp_float *y_init = cupdlp_NULL;
//   CUPDLP_INIT_ZERO(x_init, x_init_len);
//   CUPDLP_INIT_ZERO(y_init, y_init_len);
//   cupdlp_printf("x_coarse_len = %d, x_len = %d\n", nCols_coarse, x_init_len);
//   cupdlp_printf("y_coarse_len = %d, y_len = %d\n", nRows_coarse, y_init_len);
//   fine_dualOT_primal(x_init, x_coarse_solution, x_init_len, nCols_coarse, resolution, coarse_degree);
//   fine_dualOT_dual(y_init, y_coarse_solution, y_init_len, nRows_coarse, resolution, coarse_degree);

//   // 找出y中零元素的index，方便后续对A进行修剪
//   cupdlp_int *y_init_zero_idx_len = cupdlp_NULL;
//   CUPDLP_INIT_ZERO(y_init_zero_idx_len, 1);
//   cupdlp_int *y_init_zero_idx = countArray1D_Smaller_than_threshold_with_Record(y_init, y_init_len, y_init_zero_idx_len, 1e-20);
//   cupdlp_int y_init_delete_len = y_init_len - *y_init_zero_idx_len;
//   cupdlp_printf("y_int_len = %d, y_init_delete_len = %d\n", y_init_len, y_init_delete_len);
//   cupdlp_float *y_init_delete = cupdlp_NULL;
//   CUPDLP_INIT_ZERO(y_init_delete, y_init_delete_len);
//   deleteArray1DElements(y_init_delete, y_init, y_init_len, y_init_zero_idx, y_init_zero_idx_len);

//   // generate_coarse_dualOT_model_from_csv(model, csvpath_1, csvpath_2, resolution, 0);
//   generate_dualOT_model_delete_from_csv(model, csvpath_1, csvpath_2, resolution, y_init_zero_idx, y_init_zero_idx_len, y_init_delete_len);
//   cupdlp_int *constraint_new_idx = NULL;
//   CUPDLPwork *w_delete = cupdlp_NULL;
//   CUPDLP_INIT_ZERO(w_delete, 1);
//   createCUPDLPwork(w_delete, model, ifChangeIntParam, intParam, &constraint_new_idx, nCols_origin_ptr, nRows_ptr);

//   cupdlp_float *x_solution = cupdlp_NULL;
//   cupdlp_float *y_solution = cupdlp_NULL;
//   CUPDLP_INIT_ZERO(x_solution, *nCols_origin_ptr);
//   CUPDLP_INIT_ZERO(y_solution, *nRows_ptr);

//   cupdlp_float *stepsize_last = cupdlp_NULL;
//   cupdlp_float *weight_last = cupdlp_NULL;
//   CUPDLP_INIT_ZERO(stepsize_last, 1);
//   CUPDLP_INIT_ZERO(weight_last, 1);
//   cupdlp_float *stepsize_init = cupdlp_NULL;
//   cupdlp_float *weight_init = cupdlp_NULL;
//   CUPDLP_INIT_ZERO(stepsize_init, 1);
//   CUPDLP_INIT_ZERO(weight_init, 1);
//   cupdlp_printf("stepsize_coarse_last = %f\n", *stepsize_coarse_last);
//   cupdlp_printf("weight_coarse_last = %f\n", *weight_coarse_last);
//   *stepsize_init = *stepsize_coarse_last;
//   *weight_init = *weight_coarse_last;
//   cupdlp_printf("stepsize_init = %f\n", *stepsize_init);
//   cupdlp_printf("weight_init = %f\n", *weight_init);

//   cupdlp_float *x_zero = cupdlp_NULL;
//   cupdlp_float *y_zero = cupdlp_NULL;
//   CUPDLP_INIT_ZERO(x_zero, x_init_len);
//   CUPDLP_INIT_ZERO(y_zero, y_init_len);

//   fout = "./solution.txt";
//   cupdlp_bool whether_first_fine = true;
//   CUPDLP_CALL(LP_SolvePDHG_Multiscale(w_delete, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, ifSaveSol, constraint_new_idx, &x_solution, &y_solution, x_init, y_init_delete, stepsize_last, weight_last, stepsize_init, weight_init, whether_first_fine));
//   // CUPDLP_CALL(LP_SolvePDHG_Multiscale(w_delete, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, ifSaveSol, constraint_new_idx, &x_solution, &y_solution, x_zero, y_zero, stepsize_last, weight_last, stepsize_init, weight_init, whether_first_fine));
//   // CUPDLP_CALL(LP_SolvePDHG_Multiscale(w_delete, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, ifSaveSol, constraint_new_idx, &x_solution, &y_solution, x_init, y_init, stepsize_last, weight_last, stepsize_init, weight_init, whether_first_fine));

//   analyseArray1D(x_coarse_solution, nCols_coarse, 1e-20, "./solution_x_coarse_optimal.csv");
//   cupdlp_printf("--------------------------------------------------\n");
//   analyseArray1D(x_init, w_delete->problem->nCols, 1e-20, "./solution_x_fine_init.csv");
//   cupdlp_printf("--------------------------------------------------\n");
//   analyseArray1D(x_solution, w_delete->problem->nCols, 1e-20, "./solution_x_fine_optimal.csv");
//   cupdlp_printf("--------------------------------------------------\n");
//   cupdlp_printf("--------------------------------------------------\n");
//   cupdlp_printf("y_int_len = %d, y_init_delete_len = %d\n", y_init_len, y_init_delete_len);
//   analyseArray1D(y_init, y_init_len, 1e-20, "./solution_y_fine_init.csv");
//   cupdlp_printf("--------------------------------------------------\n");
//   analyseArray1D(y_init_delete, w_delete->problem->nRows, 1e-20, "./solution_y_fine_delete_init.csv");
//   cupdlp_printf("--------------------------------------------------\n");
//   analyseArray1D(y_solution, w_delete->problem->nRows, 1e-20, "./solution_y_fine_optimal.csv");
//   cupdlp_printf("--------------------------------------------------\n");
//   compareTwoArray1D(y_init, w_delete->problem->nRows, y_solution, w_delete->problem->nRows, 1e-20);
#pragma endregion

#pragma endregion

#pragma region 未考虑稀疏性的版本
  // cupdlp_printf("--------------------------------------------------\n");
  // cupdlp_printf("--------------------------------------------------\n");
  // generate_coarse_dualOT_model_from_csv(model, csvpath_1, csvpath_2, resolution, 0);
  // int *constraint_new_idx = NULL;
  // CUPDLPwork *w = cupdlp_NULL;
  // CUPDLP_INIT_ZERO(w, 1);
  // createCUPDLPwork(w, model, ifChangeIntParam, intParam, &constraint_new_idx, nCols_origin_ptr, nRows_ptr);
  // // coarse的最优解作为fine的初始解
  // cupdlp_float *x_init = cupdlp_NULL;
  // cupdlp_float *y_init = cupdlp_NULL;
  // CUPDLP_INIT_ZERO(x_init, w->problem->nCols);
  // CUPDLP_INIT_ZERO(y_init, w->problem->nRows);
  // cupdlp_printf("x_coarse_len = %d, x_len = %d\n", nCols_coarse, w->problem->nCols);
  // cupdlp_printf("y_coarse_len = %d, y_len = %d\n", nRows_coarse, w->problem->nRows);
  // fine_dualOT_primal(x_init, x_coarse_solution, w->problem->nCols, nCols_coarse, resolution, coarse_degree);
  // fine_dualOT_dual(y_init, y_coarse_solution, w->problem->nRows, nRows_coarse, resolution, coarse_degree);

  // cupdlp_float *x_solution = cupdlp_NULL;
  // cupdlp_float *y_solution = cupdlp_NULL;
  // CUPDLP_INIT_ZERO(x_solution, *nCols_origin_ptr);
  // CUPDLP_INIT_ZERO(y_solution, *nRows_ptr);

  // cupdlp_float *stepsize_last = cupdlp_NULL;
  // cupdlp_float *weight_last = cupdlp_NULL;
  // CUPDLP_INIT_ZERO(stepsize_last, 1);
  // CUPDLP_INIT_ZERO(weight_last, 1);
  // cupdlp_float *stepsize_init = cupdlp_NULL;
  // cupdlp_float *weight_init = cupdlp_NULL;
  // CUPDLP_INIT_ZERO(stepsize_init, 1);
  // CUPDLP_INIT_ZERO(weight_init, 1);
  // cupdlp_printf("stepsize_coarse_last = %f\n", *stepsize_coarse_last);
  // cupdlp_printf("weight_coarse_last = %f\n", *weight_coarse_last);
  // *stepsize_init = *stepsize_coarse_last;
  // *weight_init = *weight_coarse_last;
  // cupdlp_printf("stepsize_init = %f\n", *stepsize_init);
  // cupdlp_printf("weight_init = %f\n", *weight_init);

  // fout = "./solution.txt";
  // CUPDLP_CALL(LP_SolvePDHG_Multiscale(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, ifSaveSol, constraint_new_idx, &x_solution, &y_solution, x_init, y_init, stepsize_last, weight_last, stepsize_init, weight_init, false));

  // analyseArray1D(x_coarse_solution, nCols_coarse, 1e-20, "./solution_x_coarse_optimal.csv");
  // cupdlp_printf("--------------------------------------------------\n");
  // analyseArray1D(x_init, w->problem->nCols, 1e-20, "./solution_x_fine_init.csv");
  // cupdlp_printf("--------------------------------------------------\n");
  // analyseArray1D(x_solution, w->problem->nCols, 1e-20, "./solution_x_fine_optimal.csv");
  // cupdlp_printf("--------------------------------------------------\n");
  // cupdlp_printf("--------------------------------------------------\n");
  // analyseArray1D(y_coarse_solution, nRows_coarse, 1e-20, "./solution_y_coarse_optimal.csv");
  // cupdlp_printf("--------------------------------------------------\n");
  // analyseArray1D(y_init, w->problem->nRows, 1e-20, "./solution_y_fine_init.csv");
  // cupdlp_printf("--------------------------------------------------\n");
  // analyseArray1D(y_solution, w->problem->nRows, 1e-20, "./solution_y_fine_optimal.csv");
  // cupdlp_printf("--------------------------------------------------\n");
  // compareTwoArray1D(y_init, w->problem->nRows, y_solution, w->problem->nRows, 1e-20);
#pragma endregion

  // cupdlp_printf("测试构造矩阵\n");
  // void *mat = cupdlp_NULL;
  // CUPDLP_INIT_ZERO(mat, 1);
  // constructCoinPackedMatrix(mat, 2);

  // switch (ifPDTEST)
  // {
  // case 0:
  //   cupdlp_printf("--------------------------------------------------\n");
  //   cupdlp_printf("enter main solve loop, PDHG\n");
  //   cupdlp_printf("--------------------------------------------------\n");
  //   CUPDLP_CALL(LP_SolvePDHG_Multiscale(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, ifSaveSol, constraint_new_idx, resolution, &x_solution, &y_solution, x_init, y_init));
  //   break;
  // case 1:
  //   cupdlp_printf("--------------------------------------------------\n");
  //   cupdlp_printf("enter main solve loop, PDTEST\n");
  //   cupdlp_printf("--------------------------------------------------\n");
  //   CUPDLP_CALL(LP_SolvePDTEST(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx));
  //   break;
  // case 2:
  //   cupdlp_printf("--------------------------------------------------\n");
  //   cupdlp_printf("enter main solve loop, PDTEST_Average\n");
  //   cupdlp_printf("--------------------------------------------------\n");
  //   CUPDLP_CALL(LP_SolvePDTEST_Average(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx));
  //   break;
  // case 3:
  //   cupdlp_printf("--------------------------------------------------\n");
  //   cupdlp_printf("enter main solve loop, PDTEST_min\n");
  //   cupdlp_printf("--------------------------------------------------\n");
  //   CUPDLP_CALL(LP_SolvePDTEST_min(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx));
  //   break;
  // case 4:
  //   cupdlp_printf("--------------------------------------------------\n");
  //   cupdlp_printf("enter main solve loop, PDTEST_best\n");
  //   cupdlp_printf("--------------------------------------------------\n");
  //   CUPDLP_CALL(LP_SolvePDTEST_best(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx));
  //   break;
  // case 5:
  //   cupdlp_printf("--------------------------------------------------\n");
  //   cupdlp_printf("enter main solve loop, PDHG_AdapTheta\n");
  //   cupdlp_printf("--------------------------------------------------\n");
  //   CUPDLP_CALL(LP_SolvePDHG_AdapTheta(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx));
  //   break;
  // default:
  //   cupdlp_printf("Error: ifPDTEST = %d, 不在取值范围内", ifPDTEST);
  //   break;
  // }
  //

  // print result
  // TODO: implement after adding IO

exit_cleanup:
  if (retcode != RETCODE_OK)
  {
    cupdlp_printf("main Exit cleanup\n");
  }
  // deleteModel(model);
  // PDHG_Destroy(&w_coarse);
  // PDHG_Destroy(&w);
  // cupdlp_free(x_origin);
  // cupdlp_free(y_origin);
  // if (ifPresolve)
  // {
  //   deletePresolve(presolveinfo);
  //   deleteModel(presolvedmodel);
  // }

  // if (scaling)
  // {
  //   scaling_clear(scaling);
  // }
  // if (cost != NULL)
  //   cupdlp_free(cost);
  // if (csc_beg != NULL)
  //   cupdlp_free(csc_beg);
  // if (csc_idx != NULL)
  //   cupdlp_free(csc_idx);
  // if (csc_val != NULL)
  //   cupdlp_free(csc_val);
  // if (rhs != NULL)
  //   cupdlp_free(rhs);
  // if (lower != NULL)
  //   cupdlp_free(lower);
  // if (upper != NULL)
  //   cupdlp_free(upper);
  // if (constraint_new_idx != NULL)
  //   cupdlp_free(constraint_new_idx);
  // if (x_origin != NULL)
  //   cupdlp_free(x_origin);
  // if (y_origin != NULL)
  //   cupdlp_free(y_origin);
  // // free memory
  // csc_clear(csc_cpu);
  // problem_clear(prob);

  // freealldata(Aeqp, Aeqi, Aeqx, Aineqp, Aineqi, Aineqx, colUbIdx, colUbElem,
  //             rhs, cost, x, s, t, sx, ss, st, y, lower, upper);

  return retcode;
}
