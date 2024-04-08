#include "cupdlp_multiscale.h"

void CUPDLP_multiscale_testprint()
{
    cupdlp_printf("CUPDLP_multiscale_testprint\n");
    return 0;
}

cupdlp_float *readCSVToFloatArray(const char *filename, cupdlp_int numRows, cupdlp_int numCols)
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

    // 现在我们知道了总的行数和列数，可以分配足够的内存
    cupdlp_float *data = (cupdlp_float *)malloc(numRows * numCols * sizeof(cupdlp_float));
    if (!data)
        return NULL; // 内存分配失败

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
    return data;
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

cupdlp_float *coarsingArray1D(cupdlp_float *array, cupdlp_int array_len, const cupdlp_int coarse_degree)
{
    cupdlp_float array_coarse_len = array_len / pow(2, coarse_degree);
    cupdlp_float *array_coarse = (cupdlp_float *)malloc(array_coarse_len * sizeof(cupdlp_float));
    cupdlp_float sum = 0.0;

    for (cupdlp_int i = 0; i < array_coarse_len; i++)
    {
        array_coarse[i] = 0.0;
        for (cupdlp_int j = 0; j < array_coarse_len; j++)
        {
            array_coarse[i] += array[i * coarse_degree + j];
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
cupdlp_float *normalizedSquaredEuclideanDistance(cupdlp_int m, cupdlp_int n)
{
    cupdlp_float *c = (cupdlp_float *)malloc(m * m * n * n * sizeof(cupdlp_float));
    if (!c)
    {
        printf("内存分配失败\n");
        return NULL;
    }
    for (cupdlp_int i = 0; i < m; i++)
    {
        for (cupdlp_int j = 0; j < n; j++)
        {
            for (cupdlp_int k = 0; k < m; k++)
            {
                for (cupdlp_int l = 0; l < n; l++)
                {
                    c[(i * n + j) * m * n + (k * n + l)] = ((i - k) * (i - k) + (j - l) * (j - l)) / (m * m + n * n + 0.0);
                    // c[(i * n + j) * m * n + (k * n + l)] = 1.0;
                }
            }
        }
    }
    return c;
}

cupdlp_float *coarsing_normalizedSquaredEuclideanDistance(cupdlp_int m, cupdlp_int n, const cupdlp_int coarse_degree)
{
    cupdlp_float *c_coarse = normalizedSquaredEuclideanDistance(m / pow(2, coarse_degree), n / pow(2, coarse_degree));
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

void generate_dualOT_model_from_csv(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution)
{
    // void *model = NULL;
    // model = createModel();
    cupdlp_printf("开始读取数据\n");
    cupdlp_int a_len = resolution * resolution;
    cupdlp_int b_len = resolution * resolution;
    cupdlp_float *a = readCSVToFloatArray(csvpath_1, resolution, resolution);
    cupdlp_float *b = readCSVToFloatArray(csvpath_2, resolution, resolution);
    normalizeArray1D(a, a_len);
    normalizeArray1D(b, b_len);
    cupdlp_float *c = normalizedSquaredEuclideanDistance(resolution, resolution);
    cupdlp_printf("开始构建模型\n");
    generate_dualOT_model_from_distribution_and_cost(model, a, b, a_len, b_len, c);
    cupdlp_printf("模型构建完成\n");
    free(c);
    free(a);
    free(b);
}

void generate_coarse_dualOT_model(void *model_coarse, cupdlp_float *a, cupdlp_float *b, cupdlp_int a_len, cupdlp_int b_len, cupdlp_int resolution, const cupdlp_int coarse_degree)
{
    cupdlp_int a_coarse_len = a_len / pow(2, coarse_degree);
    cupdlp_int b_coarse_len = b_len / pow(2, coarse_degree);
    cupdlp_float *a_coarse = coarsingArray1D(a, a_len, coarse_degree);
    cupdlp_float *b_coarse = coarsingArray1D(b, b_len, coarse_degree);
    cupdlp_float *c_coarse = coarsing_normalizedSquaredEuclideanDistance(resolution, resolution, coarse_degree);
    generate_dualOT_model_from_distribution_and_cost(model_coarse, a_coarse, b_coarse, a_coarse_len, b_coarse_len, c_coarse);
    free(a_coarse);
    free(b_coarse);
    free(c_coarse);
}

void LP_Solve_Multiscale(w, w_coarse, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx)
{
    LP_SolvePDHG(w_coarse, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx);
    // 提取x_coarse, y_coarse, 转化为下一层的x_origin, y_origin
    LP_SolvePDHG(w_coarse, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx);
}

CUPDLPwork *createCUPDLPwork(void *model, CUPDLPscaling *scaling, cupdlp_int ifScaling, cupdlp_bool ifPresolve)
{
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
    int *constraint_new_idx = NULL;
    int nCols_origin;
    cupdlp_retcode retcode = RETCODE_OK;
    CUPDLP_MATRIX_FORMAT src_matrix_format = CSC;     // 原矩阵格式
    CUPDLP_MATRIX_FORMAT dst_matrix_format = CSR_CSC; // 目标矩阵格式
    cupdlp_float alloc_matrix_time = 0.0;
    cupdlp_float copy_vec_time = 0.0;

    void *presolveinfo = NULL;
    void *presolvedmodel = NULL;

    cupdlp_float *x_origin = cupdlp_NULL;
    cupdlp_float *y_origin = cupdlp_NULL;

    cupdlp_float presolve_time = getTimeStamp();
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
                                &constraint_new_idx));

    CUPDLPwork *w = cupdlp_NULL;
    CUPDLP_INIT_ZERO(w, 1);
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

    CUPDLP_INIT(x_origin, nCols_origin);
    CUPDLP_INIT(y_origin, nRows);
exit_cleanup:
    cupdlp_printf("exit_cleanup\n");

    return w;
}