#include "cupdlp_multiscale.h"

void CUPDLP_multiscale_testprint()
{
    cupdlp_printf("CUPDLP_multiscale_testprint\n");
    return 0;
}

void readCSVToFloatArray(cupdlp_float *data, const char *filename, const int numCols)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        printf("无法打开文件 %s\n", filename);
        fclose(file);
        return NULL;
    }

    cupdlp_int buffer_length = 20 * numCols;
    char buffer[buffer_length]; // 假设每行数据不会超过buffer_length字节

    cupdlp_int row = 0;
    while (fgets(buffer, sizeof(buffer), file))
    {
        char *token = strtok(buffer, ",");
        cupdlp_int col = 0;
        while (token)
        {
            // 计算当前元素在一维数组中的索引，并填充数据
            cupdlp_int idx = row * numCols + col;
            data[idx] = atof(token);
            token = strtok(NULL, ",");
            col++;
        }
        row++;
    }

    fclose(file);
}

void normalizeArray1D(cupdlp_float *array, cupdlp_int array_len)
{
    cupdlp_float sum = 0.0;
    for (cupdlp_int i = 0; i < array_len; i++)
    {
        sum += array[i];
    }
    for (cupdlp_int i = 0; i < array_len; i++)
    {
        array[i] = array[i] / sum;
    }
}

void coarsingArray1D(cupdlp_float *array_coarse, cupdlp_float *array, cupdlp_int resolution, const cupdlp_int coarse_degree)
{
    cupdlp_float sum = 0.0;
    cupdlp_int scale = pow(2, coarse_degree);
    cupdlp_int resolution_coarse = resolution / scale;
    // for (cupdlp_int i = 0; i < array_coarse_len; i++)
    // {
    //     array_coarse[i] = 0.0;
    //     for (cupdlp_int j = 0; j < scale; j++)
    //     {
    //         array_coarse[i] += array[i * scale + j];
    //     }
    // }
    for (cupdlp_int i = 0; i < resolution_coarse; i++)
    {
        for (cupdlp_int j = 0; j < resolution_coarse; j++)
        {
            array_coarse[i * resolution_coarse + j] = 0.0;
            for (cupdlp_int k = 0; k < scale; k++)
            {
                for (cupdlp_int l = 0; l < scale; l++)
                {
                    array_coarse[i * resolution_coarse + j] += array[(i * scale + k) * resolution + j * scale + l];
                }
            }
        }
    }
}

void print_float_array2D(cupdlp_float **array, cupdlp_int numRows, cupdlp_int numCols)
{
    if (!array)
    {
        printf("数组为空\n");
        return;
    }
    cupdlp_printf("读取的数据：\n");
    for (int i = 0; i < numRows; i++)
    {
        for (int j = 0; j < numCols; j++)
        {
            cupdlp_printf("%f ", array[i][j]);
            cupdlp_printf("\n");
        }
    }
}

void print_float_array1D(cupdlp_float *array, cupdlp_int num)
{
    if (!array)
    {
        printf("数组为空\n");
        return;
    }
    cupdlp_printf("读取的数据：\n");
    for (int i = 0; i < num; i++)
    {
        cupdlp_printf("%f ", array[i]);
        cupdlp_printf("\n");
    }
}
void print_int_array1D(cupdlp_int *array, cupdlp_int num)
{
    if (!array)
    {
        printf("数组为空\n");
        return;
    }
    cupdlp_printf("读取的数据：\n");
    for (int i = 0; i < num; i++)
    {
        // cupdlp_printf("%f ", array[i]);
        cupdlp_printf("%d ", array[i]);
        cupdlp_printf("\n");
    }
}

cupdlp_float **mergeTwoArrays2D(cupdlp_float **a, cupdlp_float **b, cupdlp_int resolution)
{
    // 创建新的二维数组
    cupdlp_float **mergedArray = (cupdlp_float **)malloc(resolution * sizeof(cupdlp_float *));
    for (cupdlp_int i = 0; i < resolution; i++)
    {
        mergedArray[i] = (cupdlp_float *)malloc(2 * resolution * sizeof(cupdlp_float));

        // 复制a数组的数据到新数组的前半部分
        for (cupdlp_int j = 0; j < resolution; j++)
        {
            mergedArray[i][j] = a[i][j];
        }

        // 复制b数组的数据到新数组的后半部分
        for (cupdlp_int j = 0; j < resolution; j++)
        {
            mergedArray[i][j + resolution] = b[i][j];
        }
    }

    return mergedArray;
}
cupdlp_float *mergeTwoArrays1D(const cupdlp_float *a, const cupdlp_float *b, int a_len, int b_len)
{
    cupdlp_float *result = (cupdlp_float *)malloc((a_len + b_len) * sizeof(cupdlp_float));

    if (!result)
    {
        printf("内存分配失败\n");
        return NULL;
    }

    // 先复制数组a到结果数组
    for (int i = 0; i < a_len; i++)
    {
        result[i] = a[i];
    }

    // 接着复制数组b到结果数组
    for (int i = 0; i < b_len; i++)
    {
        result[a_len + i] = b[i];
    }

    return result;
}

cupdlp_float *mergeTwoArrays1D_minus(const cupdlp_float *a, const cupdlp_float *b, cupdlp_int a_len, cupdlp_int b_len)
{
    cupdlp_float *result = (cupdlp_float *)malloc((a_len + b_len) * sizeof(cupdlp_float));

    if (!result)
    {
        printf("内存分配失败\n");
        return NULL;
    }

    // 先复制数组a到结果数组
    for (int i = 0; i < a_len; i++)
    {
        result[i] = -a[i];
    }

    // 接着复制数组b到结果数组
    for (int i = 0; i < b_len; i++)
    {
        result[a_len + i] = -b[i];
    }

    return result;
}

// m是图片的行数，n是图片的列数，m*n是分布的长度，m*m*n*n是c的长度
void normalizedSquaredEuclideanDistance(cupdlp_float *c, cupdlp_int m, cupdlp_int n)
{
    for (cupdlp_int i = 0; i < m; i++)
    {
        for (cupdlp_int j = 0; j < n; j++)
        {
            for (cupdlp_int k = 0; k < m; k++)
            {
                for (cupdlp_int l = 0; l < n; l++)
                {
                    // 从第一张图第i行第j列，到第二张图第k行第l列，需要的cost
                    c[(i * n + j) * m * n + (k * n + l)] = ((i - k) * (i - k) + (j - l) * (j - l)) / (m * m + n * n + 0.0);
                }
            }
        }
    }
}

void coarsing_normalizedSquaredEuclideanDistance(cupdlp_float *c_coarse, cupdlp_int m, cupdlp_int n, const cupdlp_int coarse_degree)
{
    normalizedSquaredEuclideanDistance(c_coarse, m / pow(2, coarse_degree), n / pow(2, coarse_degree));
    return c_coarse;
}

// start数组定义了每一列的起始元素在 index[] 和 value[] 数组中的位置，或者说每一列的起始元素是第几个非零元素；start数组的最后一个为非零元素个数
cupdlp_int *dualOT_startArray(cupdlp_int m, cupdlp_int n)
{
    cupdlp_int *startArray = (cupdlp_int *)malloc((n + m + 1) * sizeof(cupdlp_int));
    if (!startArray)
    {
        printf("内存分配失败\n");
        return NULL;
    }
    for (cupdlp_int i = 0; i < n; i++)
    {
        startArray[i] = i * m;
    }
    for (cupdlp_int j = 0; j < m; j++)
    {
        startArray[j + n] = n * m + n * j;
    }
    startArray[n + m] = 2 * m * n;
    return startArray;
}

// index数组定义了每一个非零元素所在的行数
cupdlp_int *dualOT_indexArray(cupdlp_int m, cupdlp_int n)
{
    cupdlp_int *indexArray = (cupdlp_int *)malloc(2 * m * n * sizeof(cupdlp_int));
    if (!indexArray)
    {
        printf("内存分配失败\n");
        return NULL;
    }
    for (cupdlp_int i = 0; i < n; i++)
    {
        for (cupdlp_int j = 0; j < m; j++)
        {
            indexArray[i * m + j] = i * m + j;
        }
    }
    // 第二部分总共m列，每列有n个元素
    for (cupdlp_int i = 0; i < m; i++)
    {
        for (cupdlp_int j = 0; j < n; j++)
        {
            indexArray[n * m + i * n + j] = i + j * m;
        }
    }
    return indexArray;
}

cupdlp_float *dualOT_valueArray(cupdlp_int m, cupdlp_int n)
{
    cupdlp_float *valueArray = (cupdlp_float *)malloc(2 * m * n * sizeof(cupdlp_float));
    if (!valueArray)
    {
        printf("内存分配失败\n");
        return NULL;
    }
    for (cupdlp_int i = 0; i < 2 * m * n; i++)
    {
        valueArray[i] = 1.0;
    }
    return valueArray;
}

cupdlp_float *dualOT_rowLower(cupdlp_int m, cupdlp_int n)
{
    cupdlp_float *rowLower = (cupdlp_float *)malloc((n * m) * sizeof(cupdlp_float));
    if (!rowLower)
    {
        printf("内存分配失败\n");
        return NULL;
    }
    for (cupdlp_int i = 0; i < n * m; i++)
    {
        rowLower[i] = -1e30;
    }
    return rowLower;
}

cupdlp_float *dualOT_colLower(cupdlp_int m, cupdlp_int n)
{
    cupdlp_float *colLower = (cupdlp_float *)malloc((m + n) * sizeof(cupdlp_float));
    if (!colLower)
    {
        printf("内存分配失败\n");
        return NULL;
    }
    for (cupdlp_int i = 0; i < m + n; i++)
    {
        colLower[i] = -1e30;
    }
    return colLower;
}

cupdlp_float *dualOT_colUpper(cupdlp_int m, cupdlp_int n)
{
    cupdlp_float *colUpper = (cupdlp_float *)malloc((m + n) * sizeof(cupdlp_float));
    if (!colUpper)
    {
        printf("内存分配失败\n");
        return NULL;
    }
    for (cupdlp_int i = 0; i < m + n; i++)
    {
        colUpper[i] = 1e30;
    }
    return colUpper;
}

void generate_dualOT_model_from_distribution_and_cost(void *model, const cupdlp_float *a, const cupdlp_float *b, const cupdlp_int a_len, const cupdlp_int b_len, cupdlp_float *c)
{
    // obj = -q
    cupdlp_float *obj = mergeTwoArrays1D_minus(a, b, a_len, b_len);
    cupdlp_int numCols = a_len + b_len;
    cupdlp_int numRows = a_len * b_len;
    cupdlp_int *start = dualOT_startArray(a_len, b_len);
    cupdlp_int *index = dualOT_indexArray(a_len, b_len);
    cupdlp_float *value = dualOT_valueArray(a_len, b_len);
    // 对x的约束(没有约束)
    cupdlp_float *rowLower = dualOT_rowLower(a_len, b_len);
    cupdlp_float *colLower = dualOT_colLower(a_len, b_len);
    cupdlp_float *colUpper = dualOT_colUpper(a_len, b_len);
    loadProblemWrapper(model, numCols, numRows, start, index, value, colLower, colUpper, obj, rowLower, c);
    free(start);
    free(index);
    free(value);
    free(rowLower);
    free(colLower);
    free(colUpper);
    free(obj);
}

void generate_dualOT_model_byMatrix_from_distribution_and_cost(void *model, const cupdlp_float *a, const cupdlp_float *b, const cupdlp_int a_len, const cupdlp_int b_len, cupdlp_float *c, int resolution)
{
    // obj = -q
    cupdlp_float *obj = mergeTwoArrays1D_minus(a, b, a_len, b_len);
    cupdlp_int numCols = a_len + b_len;
    cupdlp_int numRows = a_len * b_len;
    // cupdlp_int *start = dualOT_startArray(a_len, b_len);
    // cupdlp_int *index = dualOT_indexArray(a_len, b_len);
    // cupdlp_float *value = dualOT_valueArray(a_len, b_len);

    // 对x的约束(没有约束)
    cupdlp_float *rowLower = dualOT_rowLower(a_len, b_len);
    cupdlp_float *colLower = dualOT_colLower(a_len, b_len);
    cupdlp_float *colUpper = dualOT_colUpper(a_len, b_len);
    loadProblem_byMatrix_Wrapper(model, resolution, colLower, colUpper, obj, rowLower, c);
    // loadProblemWrapper(model, numCols, numRows, start, index, value, colLower, colUpper, obj, rowLower, c);
    // free(start);
    // free(index);
    // free(value);
    free(rowLower);
    free(colLower);
    free(colUpper);
    free(obj);
}

void generate_dualOT_model_from_csv(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution)
{
    cupdlp_retcode retcode = RETCODE_OK;

    // void *model = NULL;
    // model = createModel();
    cupdlp_printf("开始读取数据\n");
    cupdlp_int a_len = resolution * resolution;
    cupdlp_int b_len = resolution * resolution;

    cupdlp_float *a = cupdlp_NULL;
    CUPDLP_INIT(a, a_len);
    readCSVToFloatArray(a, csvpath_1, resolution);

    cupdlp_float *b = cupdlp_NULL;
    CUPDLP_INIT(b, b_len);
    readCSVToFloatArray(b, csvpath_2, resolution);

    normalizeArray1D(a, a_len);
    normalizeArray1D(b, b_len);

    cupdlp_float *c = cupdlp_NULL;
    CUPDLP_INIT(c, a_len * b_len);
    normalizedSquaredEuclideanDistance(c, resolution, resolution);

    cupdlp_printf("开始构建模型\n");
    // generate_dualOT_model_from_distribution_and_cost(model, a, b, a_len, b_len, c);
    generate_dualOT_model_byMatrix_from_distribution_and_cost(model, a, b, a_len, b_len, c, resolution);
    cupdlp_printf("模型构建完成\n");

exit_cleanup:
{
    if (retcode != RETCODE_OK)
    {
        cupdlp_printf("generate_dualOT_model_from_csv, exit_cleanup\n");
    }
    cupdlp_free(a);
    cupdlp_free(b);
    cupdlp_free(c);
}
}

void generate_coarse_dualOT_model_from_csv(void *model_coarse, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, const cupdlp_int coarse_degree)
{
    cupdlp_retcode retcode = RETCODE_OK;

    cupdlp_printf("开始读取数据\n");
    cupdlp_int a_len = resolution * resolution;
    cupdlp_int b_len = resolution * resolution;

    cupdlp_float *a = cupdlp_NULL;
    CUPDLP_INIT(a, a_len);
    readCSVToFloatArray(a, csvpath_1, resolution);

    cupdlp_float *b = cupdlp_NULL;
    CUPDLP_INIT(b, b_len);
    readCSVToFloatArray(b, csvpath_2, resolution);

    normalizeArray1D(a, a_len);
    normalizeArray1D(b, b_len);

    cupdlp_int a_coarse_len = a_len / pow(pow(2, coarse_degree), 2);
    cupdlp_int b_coarse_len = b_len / pow(pow(2, coarse_degree), 2);

    cupdlp_float *a_coarse = cupdlp_NULL;
    CUPDLP_INIT(a_coarse, a_coarse_len);
    coarsingArray1D(a_coarse, a, resolution, coarse_degree);

    cupdlp_float *b_coarse = cupdlp_NULL;
    CUPDLP_INIT(b_coarse, b_coarse_len);
    coarsingArray1D(b_coarse, b, resolution, coarse_degree);

    cupdlp_float *c_coarse = cupdlp_NULL;
    CUPDLP_INIT(c_coarse, a_coarse_len * b_coarse_len);
    coarsing_normalizedSquaredEuclideanDistance(c_coarse, resolution, resolution, coarse_degree);

    cupdlp_printf("开始构建模型\n");
    // generate_dualOT_model_from_distribution_and_cost(model_coarse, a_coarse, b_coarse, a_coarse_len, b_coarse_len, c_coarse);
    generate_dualOT_model_byMatrix_from_distribution_and_cost(model_coarse, a_coarse, b_coarse, a_coarse_len, b_coarse_len, c_coarse, resolution / pow(2, coarse_degree));
    cupdlp_printf("模型构建完成\n");
exit_cleanup:
{
    if (retcode != RETCODE_OK)
    {
        cupdlp_printf("generate_coarse_dualOT_model_from_csv, exit_cleanup\n");
    }
    cupdlp_free(a);
    cupdlp_free(b);
    cupdlp_free(a_coarse);
    cupdlp_free(b_coarse);
    cupdlp_free(c_coarse);
}
}

void LP_Solve_Multiscale(w, w_coarse, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx)
{
    LP_SolvePDHG(w_coarse, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx);
    // 提取x_coarse, y_coarse, 转化为下一层的x_origin, y_origin
    LP_SolvePDHG(w_coarse, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx);
}

void createCUPDLPwork(CUPDLPwork *w, void *model, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_int **constraint_new_idx, int *nCols_origin_ptr, int *nRows_ptr)
{
    cupdlp_retcode retcode = RETCODE_OK;

    void *model2solve = model;
    int nCols;
    int nRows;
    CUPDLPcsc *csc_cpu = cupdlp_NULL;
    int nnz = 0;
    int *csc_beg = NULL, *csc_idx = NULL;
    double *csc_val = NULL;
    double *cost = NULL;
    int nEqs;
    double *rhs = NULL;
    cupdlp_float *lower = NULL;
    cupdlp_float *upper = NULL;
    double offset = 0.0;    // true objVal = sig * c'x - offset, sig = 1 (min) or -1 (max)
    double sign_origin = 1; // 1 (min) or -1 (max)
    int nCols_origin;

    CUPDLP_MATRIX_FORMAT src_matrix_format = CSC;     // 原矩阵格式
    CUPDLP_MATRIX_FORMAT dst_matrix_format = CSR_CSC; // 目标矩阵格式
    cupdlp_float alloc_matrix_time = 0.0;
    cupdlp_float copy_vec_time = 0.0;

    void *presolveinfo = NULL;
    void *presolvedmodel = NULL;

    // cupdlp_float *x_origin = cupdlp_NULL;
    // cupdlp_float *y_origin = cupdlp_NULL;

    cupdlp_float presolve_time = getTimeStamp();

    cupdlp_bool ifPresolve = false;
    if (ifChangeIntParam[IF_PRESOLVE])
    {
        ifPresolve = intParam[IF_PRESOLVE];
    }
    if (ifPresolve)
    {
        presolveinfo = createPresolve();
        presolvedmodel = presolvedModel(presolveinfo, model);
        model2solve = presolvedmodel;
    }
    presolve_time = getTimeStamp() - presolve_time;

    CUPDLP_CALL(formulateLP_new(model, &cost, &nCols, &nRows, &nnz, &nEqs,
                                &csc_beg, &csc_idx, &csc_val, &rhs, &lower,
                                &upper, &offset, &sign_origin, &nCols_origin,
                                constraint_new_idx));
    cupdlp_printf("nCols_origin: %d\n", nCols_origin);
    cupdlp_printf("nCols: %d\n", nCols);
    cupdlp_printf("nRows: %d\n", nRows);

    cupdlp_float cuda_prepare_time = getTimeStamp();
    CHECK_CUSPARSE(cusparseCreate(&w->cusparsehandle));
    CHECK_CUBLAS(cublasCreate(&w->cublashandle));
    cuda_prepare_time = getTimeStamp() - cuda_prepare_time;

    CUPDLPproblem *prob = cupdlp_NULL;
    CUPDLP_CALL(problem_create(&prob));

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
    csc_cpu->cuda_csc = NULL;

    CUPDLPscaling *scaling =
        (CUPDLPscaling *)cupdlp_malloc(sizeof(CUPDLPscaling));
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
    cupdlp_float scaling_time = getTimeStamp();
    CUPDLP_CALL(PDHG_Scale_Data_cuda(csc_cpu, ifScaling, scaling, cost, lower,
                                     upper, rhs));
    scaling_time = getTimeStamp() - scaling_time;

    CUPDLP_CALL(problem_alloc(prob, nRows, nCols, nEqs, cost, offset, sign_origin, csc_cpu, src_matrix_format, dst_matrix_format, rhs, lower, upper, &alloc_matrix_time, &copy_vec_time));

    w->problem = prob;
    w->scaling = scaling;
    PDHG_Alloc(w);
    w->timers->dScalingTime = scaling_time;
    w->timers->dPresolveTime = presolve_time;
    CUPDLP_COPY_VEC(w->rowScale, scaling->rowScale, cupdlp_float, nRows);
    CUPDLP_COPY_VEC(w->colScale, scaling->colScale, cupdlp_float, nCols);
    w->timers->AllocMem_CopyMatToDeviceTime += alloc_matrix_time;
    w->timers->CopyVecToDeviceTime += copy_vec_time;
    w->timers->CudaPrepareTime = cuda_prepare_time;

    // CUPDLP_INIT(x_origin, nCols_origin);
    // CUPDLP_INIT(y_origin, nRows);
    nCols_origin_ptr = &nCols_origin;
    nRows_ptr = &nRows;
exit_cleanup:
    if (retcode != RETCODE_OK)
    {
        cupdlp_printf("createCUPDLPwork exit_cleanup\n");
    }
}

void fine_dualOT_primal(cupdlp_float *x_init, cupdlp_float *x_coarse_solution, cupdlp_int x_len, cupdlp_int x_coarse_len, cupdlp_int resolution, cupdlp_int coarse_degree)
{
    cupdlp_int scale = pow(2, coarse_degree);
    cupdlp_int resolution_coarse = resolution / scale;
    if (x_len != x_coarse_len * pow(scale, 2))
    {
        cupdlp_printf("x_len != x_coarse_len * pow(scale, 2)\n");
        return NULL;
    }
    // for (cupdlp_int i = 0; i < x_coarse_len; i++)
    // {
    //     for (cupdlp_int j = 0; j < scale; j++)
    //     {
    //         x_init[i * scale + j] = x_coarse_solution[i];
    //     }
    // }
    cupdlp_int idx_temp = 0;
    cupdlp_int idx_coarse_temp = 0;
    for (cupdlp_int i = 0; i < resolution_coarse; i++)
    {
        for (cupdlp_int j = 0; j < resolution_coarse; j++)
        {
            for (cupdlp_int k = 0; k < scale; k++)
            {
                for (cupdlp_int l = 0; l < scale; l++)
                {
                    x_init[(i * scale + k) * resolution + j * scale + l] = x_coarse_solution[i * resolution_coarse + j];
                    idx_temp = (i * scale + k) * resolution + j * scale + l + pow(resolution, 2);

                    idx_coarse_temp = i * resolution_coarse + j + pow(resolution_coarse, 2);
                    x_init[idx_temp] = x_coarse_solution[idx_coarse_temp];
                }
            }
        }
    }
    cupdlp_printf("fine_dualOT_primal完成\n");
}

void fine_dualOT_dual(cupdlp_float *y_init, cupdlp_float *y_coarse_solution, cupdlp_int y_len, cupdlp_int y_coarse_len, cupdlp_int resolution, cupdlp_int coarse_degree)
{
    cupdlp_int scale = pow(2, coarse_degree);
    if (y_len != y_coarse_len * pow(scale, 4))
    {
        cupdlp_printf("y_len != y_coarse_len * pow(scale, 4)\n");
        return NULL;
    }
    cupdlp_printf("fine_dualOT_dual开始\n");
    cupdlp_float y_temp = 0.0;
    cupdlp_int resolution_coarse = resolution / scale;
    cupdlp_int idx_coarse_temp = 0;
    cupdlp_int idx_1 = 0;
    cupdlp_int idx_2 = 0;
    cupdlp_int idx_temp = 0;
    cupdlp_int idx_max = 0;
    for (cupdlp_int i1 = 0; i1 < resolution_coarse; i1++)
    {
        for (cupdlp_int j1 = 0; j1 < resolution_coarse; j1++)
        {
            for (cupdlp_int i2 = 0; i2 < resolution_coarse; i2++)
            {
                for (cupdlp_int j2 = 0; j2 < resolution_coarse; j2++)
                {
                    idx_coarse_temp = (i1 * resolution_coarse + j1) * pow(resolution_coarse, 2) + (i2 * resolution_coarse + j2);
                    y_temp = y_coarse_solution[idx_coarse_temp];
                    // 以coarse_degree=1为例，一个粗点对之间的输运变成了4*4个细点对之间的输运
                    y_temp = y_temp / (pow(scale, 4) + 0.0);
                    for (cupdlp_int k1 = 0; k1 < scale; k1++)
                    {
                        for (cupdlp_int l1 = 0; l1 < scale; l1++)
                        {
                            for (cupdlp_int k2 = 0; k2 < scale; k2++)
                            {
                                for (cupdlp_int l2 = 0; l2 < scale; l2++)
                                {
                                    idx_1 = (i1 * scale + k1) * resolution + j1 * scale + l1;
                                    idx_2 = (i2 * scale + k2) * resolution + j2 * scale + l2;
                                    idx_temp = idx_1 * pow(resolution, 2) + idx_2;
                                    // if (idx_max < idx_temp)
                                    // {
                                    //     idx_max = idx_temp;
                                    //     printf("idx_init: %d\n", idx_max);
                                    // }

                                    y_init[idx_temp] = y_temp;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    cupdlp_printf("fine_dualOT_dual完成\n");
}

int countArray1D_Smaller_than_threshold(cupdlp_float *a, int a_len, cupdlp_float thr)
{
    int count = 0;
    for (int i = 0; i < a_len; i++)
    {
        if (fabs(a[i]) < thr)
            count++;
    }
    return count;
}

int *countArray1D_Smaller_than_threshold_with_Record(cupdlp_float *a, int a_len, int *a_record_len, cupdlp_float thr)
{

    int count = 0;
    for (int i = 0; i < a_len; i++)
    {
        if (fabs(a[i]) < thr)
            count++;
    }
    int *record = (int *)malloc(count * sizeof(int));
    *a_record_len = count;
    count = 0;
    for (int i = 0; i < a_len; i++)
    {
        if (fabs(a[i]) < thr)
        {
            record[count] = i;
            count++;
        }
    }
    return record;
}

void saveArray1D_to_csv(cupdlp_float *a, int a_len, const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        printf("无法创建文件\n");
        return;
    }

    for (int i = 0; i < a_len; i++)
    {
        fprintf(file, "%f\n", a[i]);
    }
    cupdlp_printf("保存完毕\n");

    fclose(file);
}

void analyseArray1D(cupdlp_float *a, cupdlp_int a_len, cupdlp_float thr, const char *filename)
{
    cupdlp_int count = countArray1D_Smaller_than_threshold(a, a_len, thr);
    cupdlp_printf("向量总长度为：%d，小于阈值%.10f的元素个数为：%d, 稀疏元素占比为：%.2f%%\n", a_len, thr, count, 100 * count / (a_len + 0.0));
    saveArray1D_to_csv(a, a_len, filename);
}

int countArray1D_same_elements(int *a, int a_len, int *b, int b_len)
{
    int count = 0;
    for (int i = 0; i < a_len; i++)
    {
        for (int j = 0; j < b_len; j++)
        {
            if (a[i] == b[j])
            {
                count++;
                break; // 找到相同元素后不需要继续比较当前元素
            }
        }
    }
    return count;
}

void compareTwoArray1D(cupdlp_float *a, cupdlp_int a_len, cupdlp_float *b, cupdlp_int b_len, cupdlp_float thr)
{
    int *a_record_len = cupdlp_NULL;
    a_record_len = (int *)malloc(1 * sizeof(int));
    int *b_record_len = cupdlp_NULL;
    b_record_len = (int *)malloc(1 * sizeof(int));

    int *a_record = countArray1D_Smaller_than_threshold_with_Record(a, a_len, a_record_len, thr);
    int *b_record = countArray1D_Smaller_than_threshold_with_Record(b, b_len, b_record_len, thr);

    int same_count = countArray1D_same_elements(a_record, *a_record_len, b_record, *b_record_len);
    cupdlp_printf("两个数组中小于阈值%.10f的元素个数分别为：%d, %d, 其中相同位置的元素个数为：%d, 分别占比：%f%%, %f%%\n", thr, *a_record_len, *b_record_len, same_count, 100 * same_count / (*a_record_len + 0.0), 100 * same_count / (*b_record_len + 0.0));
}

void constructCoinPackedMatrix(void *mat, cupdlp_int resolution)
{
    printf("开始构建矩阵\n");
    Construct_dualOT_Matrix(mat, resolution);
}
