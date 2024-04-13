#include "cupdlp_multiscale.h"

void CUPDLP_multiscale_testprint()
{
    cupdlp_printf("CUPDLP_multiscale_testprint\n");
    return;
}

void readCSVToFloatArray(cupdlp_float *data, const char *filename, const int numCols)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        printf("无法打开文件 %s\n", filename);
        fclose(file);
        return;
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
    cupdlp_printf("数组长度为：%d, 读取的数据：\n", num);
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
    cupdlp_printf("c_coarse长度为: %d\n", m * n * m * n);
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
void normalizedSquaredEuclideanDistance_longlong(cupdlp_float *c, cupdlp_int m, cupdlp_int n)
{
    // long long c_coarse_len = m * m * n * n;
    // cupdlp_printf("c_coarse长度为: %lld\n", c_coarse_len);
    long long idx = 0;
    long long idx_1 = 0;
    long long idx_2 = 0;
    int vec_len = m * n;
    for (cupdlp_int i = 0; i < m; i++)
    {
        for (cupdlp_int j = 0; j < n; j++)
        {
            for (cupdlp_int k = 0; k < m; k++)
            {

                printf("i: %d, j: %d, k: %d\n", i, j, k);
                for (cupdlp_int l = 0; l < n; l++)
                {
                    // 从第一张图第i行第j列，到第二张图第k行第l列，需要的cost
                    idx_1 = i * n + j;
                    idx_1 = idx_1 * vec_len;
                    idx_2 = k * n + l;
                    idx = idx_1 + idx_2;
                    c[idx] = ((i - k) * (i - k) + (j - l) * (j - l)) / (m * m + n * n + 0.0);
                }
            }
        }
    }
}

void coarsing_normalizedSquaredEuclideanDistance(cupdlp_float *c_coarse, cupdlp_int m, cupdlp_int n, const cupdlp_int coarse_degree)
{
    cupdlp_int m_coarse = m / pow(2, coarse_degree);
    cupdlp_int n_coarse = n / pow(2, coarse_degree);
    normalizedSquaredEuclideanDistance(c_coarse, m_coarse, n_coarse);
}

void coarsing_normalizedSquaredEuclideanDistance_longlong(cupdlp_float *c_coarse, cupdlp_int m, cupdlp_int n, const cupdlp_int coarse_degree)
{
    cupdlp_int m_coarse = m / pow(2, coarse_degree);
    cupdlp_int n_coarse = n / pow(2, coarse_degree);
    normalizedSquaredEuclideanDistance_longlong(c_coarse, m_coarse, n_coarse);
}

void generate_c_coarse_delete_directly(cupdlp_float *c_coarse_delete, cupdlp_int coarse_degree, cupdlp_int m, cupdlp_int n, cupdlp_bool *keep, long long *keep_idx)
{
    // m_coarse和n_coarse是当前分辨率下的图片长宽
    // keep记录了该元素是否被保留
    // keep[i]若为true，则是第keep_idx[i]个被保留的元素
    cupdlp_int m_coarse = m / pow(2, coarse_degree);
    cupdlp_int n_coarse = n / pow(2, coarse_degree);
    long long idx = 0;
    long long idx_1 = 0;
    long long idx_2 = 0;
    long long c_coarse_delete_idx = 0;
    int vec_len = m_coarse * n_coarse;
    float scale_const = pow(m_coarse, 2) + pow(n_coarse, 2) + 0.0;

    for (cupdlp_int i = 0; i < m_coarse; i++)
    {
        cupdlp_printf("i: %d\n", i);
        for (cupdlp_int j = 0; j < n_coarse; j++)
        {
            for (cupdlp_int k = 0; k < m_coarse; k++)
            {
                // printf("i: %d, j: %d, k: %d\n", i, j, k);
                for (cupdlp_int l = 0; l < n_coarse; l++)
                {

                    // 从第一张图第i行第j列，到第二张图第k行第l列，需要的cost
                    // idx = (i * n_coarse + j) * m_coarse * n_coarse + (k * n_coarse + l);
                    idx_1 = i * n_coarse + j;
                    idx_1 = idx_1 * vec_len;
                    idx_2 = k * n_coarse + l;
                    idx = idx_1 + idx_2;
                    if (keep[idx])
                    {
                        // printf("keep_idx[%lld]: %lld\n", idx, keep_idx[idx]);
                        // c_coarse_delete_idx = keep_idx[idx];
                        // printf("c_coarse_delete_idx: %lld\n", c_coarse_delete_idx);
                        c_coarse_delete[c_coarse_delete_idx] = ((i - k) * (i - k) + (j - l) * (j - l)) / scale_const;
                        c_coarse_delete_idx += 1;
                    }
                    // c_coarse_delete[idx] = ((i - k) * (i - k) + (j - l) * (j - l)) / scale_const;
                }
            }
        }
    }
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
    cupdlp_printf("rowLower长度为: %d\n", n * m);
    return rowLower;
}

cupdlp_float *dualOT_rowLower_delete(cupdlp_int len_after_delete)
{
    cupdlp_float *rowLower = (cupdlp_float *)malloc((len_after_delete) * sizeof(cupdlp_float));
    if (!rowLower)
    {
        printf("内存分配失败\n");
        return NULL;
    }
    for (cupdlp_int i = 0; i < len_after_delete; i++)
    {
        rowLower[i] = -1e30;
    }
    cupdlp_printf("rowLower长度为: %d\n", len_after_delete);
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
    cupdlp_printf("colLower长度为: %d\n", m + n);
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
    cupdlp_printf("colUpper长度为: %d\n", m + n);
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

void generate_dualOT_model_byMatrix_from_distribution_and_cost(void *model, const cupdlp_float *a, const cupdlp_float *b, const cupdlp_int a_len, const cupdlp_int b_len, cupdlp_float *c, cupdlp_int resolution)
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
    cupdlp_printf("begin loadProblem_byMatrix_Wrapper\n");
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

void generate_dualOT_model_delete_byMatrix_from_distribution_and_cost(void *model, const cupdlp_float *a, const cupdlp_float *b, const cupdlp_int a_len, const cupdlp_int b_len, cupdlp_float *c_delete, int resolution, int *zero_idx, int *zero_idx_len)
{
    // obj = -q
    cupdlp_float *obj = mergeTwoArrays1D_minus(a, b, a_len, b_len);
    cupdlp_int numCols = a_len + b_len;
    cupdlp_int numRows = a_len * b_len;
    // cupdlp_int *start = dualOT_startArray(a_len, b_len);
    // cupdlp_int *index = dualOT_indexArray(a_len, b_len);
    // cupdlp_float *value = dualOT_valueArray(a_len, b_len);

    // 对x的约束(没有约束)
    // cupdlp_float *rowLower = dualOT_rowLower(a_len, b_len);
    cupdlp_float *rowLower = dualOT_rowLower_delete(a_len * b_len - *zero_idx_len);
    cupdlp_float *colLower = dualOT_colLower(a_len, b_len);
    cupdlp_float *colUpper = dualOT_colUpper(a_len, b_len);
    loadProblem_delete_byMatrix_Wrapper(model, resolution, zero_idx, zero_idx_len, colLower, colUpper, obj, rowLower, c_delete);
    // loadProblemWrapper(model, numCols, numRows, start, index, value, colLower, colUpper, obj, rowLower, c);
    // free(start);
    // free(index);
    // free(value);
    free(rowLower);
    free(colLower);
    free(colUpper);
    free(obj);
}

void generate_dualOT_model_delete_byMatrix_from_distribution_and_cost_longlong(void *model, const cupdlp_float *a, const cupdlp_float *b, const long long a_len, const long long b_len, cupdlp_float *c_delete, int resolution, long long *zero_idx, long long *zero_idx_len)
{
    // obj = -q
    cupdlp_float *obj = mergeTwoArrays1D_minus(a, b, a_len, b_len);
    cupdlp_int numCols = a_len + b_len;
    long long numRows = a_len * b_len;
    // cupdlp_int *start = dualOT_startArray(a_len, b_len);
    // cupdlp_int *index = dualOT_indexArray(a_len, b_len);
    // cupdlp_float *value = dualOT_valueArray(a_len, b_len);

    // 对x的约束(没有约束)
    // cupdlp_float *rowLower = dualOT_rowLower(a_len, b_len);
    cupdlp_float *rowLower = dualOT_rowLower_delete(numRows - *zero_idx_len);
    cupdlp_float *colLower = dualOT_colLower(a_len, b_len);
    cupdlp_float *colUpper = dualOT_colUpper(a_len, b_len);
    loadProblem_delete_byMatrix_Wrapper_longlong(model, resolution, zero_idx, zero_idx_len, colLower, colUpper, obj, rowLower, c_delete);
    // loadProblemWrapper(model, numCols, numRows, start, index, value, colLower, colUpper, obj, rowLower, c);
    // free(start);
    // free(index);
    // free(value);
    free(rowLower);
    free(colLower);
    free(colUpper);
    free(obj);
}

void generate_dualOT_model_delete_by_keep_byMatrix_from_distribution_and_cost_longlong(void *model, const cupdlp_float *a, const cupdlp_float *b, const long long a_len, const long long b_len, cupdlp_float *c_delete, int resolution, cupdlp_bool *keep, long long *keep_idx, long long *len_after_delete)
{
    cupdlp_bool retcode = RETCODE_OK;
    // obj = -q
    cupdlp_float *obj = mergeTwoArrays1D_minus(a, b, a_len, b_len);
    cupdlp_int numCols = a_len + b_len;
    long long numRows = a_len * b_len;
    // cupdlp_int *start = dualOT_startArray(a_len, b_len);
    // cupdlp_int *index = dualOT_indexArray(a_len, b_len);
    // cupdlp_float *value = dualOT_valueArray(a_len, b_len);

    // 对x的约束(没有约束)
    // cupdlp_float *rowLower = dualOT_rowLower(a_len, b_len);
    int len_after_delete_int = *len_after_delete;
    cupdlp_float *rowLower = dualOT_rowLower_delete(len_after_delete_int);
    cupdlp_float *colLower = dualOT_colLower(a_len, b_len);
    cupdlp_float *colUpper = dualOT_colUpper(a_len, b_len);
    bool *keep_Cpp = cupdlp_NULL;
    long long keep_len = pow(resolution, 4);
    CUPDLP_INIT(keep_Cpp, keep_len);
    for (long long i = 0; i < keep_len; i++)
    {
        keep_Cpp[i] = keep[i];
    }
    // for (long long i = 0; i < *len_after_delete; i++)
    // {
    //     if (keep[i])
    //     {
    //         printf("keep_true_idx: %lld\n", keep_idx[i]);
    //     }
    // }
    loadProblem_delete_by_keep_byMatrix_Wrapper_longlong(model, resolution, keep_Cpp, keep_idx, len_after_delete, colLower, colUpper, obj, rowLower, c_delete);
exit_cleanup:
    if (retcode != RETCODE_OK)
    {
        cupdlp_printf("generate_dualOT_model_delete_by_keep_byMatrix_from_distribution_and_cost_longlong, exit_cleanup\n");
    }
    cupdlp_free(rowLower);
    cupdlp_free(colLower);
    cupdlp_free(colUpper);
    cupdlp_free(obj);
    cupdlp_free(keep_Cpp)
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

    cupdlp_int resolution_coarse = resolution / pow(2, coarse_degree);
    cupdlp_printf("开始读取csv数据\n");
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
    cupdlp_int c_coarse_len = a_coarse_len * b_coarse_len;
    CUPDLP_INIT(c_coarse, c_coarse_len);
    coarsing_normalizedSquaredEuclideanDistance(c_coarse, resolution, resolution, coarse_degree);
    cupdlp_printf("csv数据读取完成，开始构建模型\n");
    // generate_dualOT_model_from_distribution_and_cost(model_coarse, a_coarse, b_coarse, a_coarse_len, b_coarse_len, c_coarse);
    cupdlp_printf("coarse_degree: %d\n", coarse_degree);
    cupdlp_printf("a_coarse_len: %d\n", a_coarse_len);
    cupdlp_printf("b_coarse_len: %d\n", b_coarse_len);
    generate_dualOT_model_byMatrix_from_distribution_and_cost(model_coarse, a_coarse, b_coarse, a_coarse_len, b_coarse_len, c_coarse, resolution_coarse);
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

void generate_dualOT_model_delete_from_csv(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_int *zero_idx, cupdlp_int *zero_idx_len, cupdlp_int len_after_delete)
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

    cupdlp_float *c_delete = cupdlp_NULL;
    CUPDLP_INIT(c_delete, len_after_delete);
    deleteArray1DElements(c_delete, c, a_len * b_len, zero_idx, zero_idx_len);

    cupdlp_printf("开始构建模型\n");
    // 统计耗时
    cupdlp_float time_generate = getTimeStamp();
    // generate_dualOT_model_from_distribution_and_cost(model, a, b, a_len, b_len, c);
    // generate_dualOT_model_byMatrix_from_distribution_and_cost(model, a, b, a_len, b_len, c_delete, resolution);
    generate_dualOT_model_delete_byMatrix_from_distribution_and_cost(model, a, b, a_len, b_len, c_delete, resolution, zero_idx, zero_idx_len);
    time_generate = getTimeStamp() - time_generate;
    cupdlp_printf("模型构建完成, 耗时: %.3f 秒\n", time_generate);

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

void generate_coarse_dualOT_model_delete_from_csv(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_int *zero_idx, cupdlp_int *zero_idx_len, cupdlp_int len_after_delete, cupdlp_int coarse_degree)
{
    cupdlp_retcode retcode = RETCODE_OK;

    cupdlp_int resolution_coarse = resolution / pow(2, coarse_degree);
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

    cupdlp_float *c_delete = cupdlp_NULL;
    CUPDLP_INIT(c_delete, len_after_delete);
    cupdlp_printf("c_delete_len: %d\n", len_after_delete);
    deleteArray1DElements(c_delete, c_coarse, a_coarse_len * b_coarse_len, zero_idx, zero_idx_len);

    cupdlp_printf("开始构建模型\n");
    // 统计耗时
    cupdlp_float time_generate = getTimeStamp();
    // generate_dualOT_model_from_distribution_and_cost(model, a, b, a_len, b_len, c);
    // generate_dualOT_model_byMatrix_from_distribution_and_cost(model, a, b, a_len, b_len, c_delete, resolution);
    generate_dualOT_model_delete_byMatrix_from_distribution_and_cost(model, a_coarse, b_coarse, a_coarse_len, b_coarse_len, c_delete, resolution_coarse, zero_idx, zero_idx_len);
    time_generate = getTimeStamp() - time_generate;
    cupdlp_printf("模型构建完成, 耗时: %.3f 秒\n", time_generate);

exit_cleanup:
{
    if (retcode != RETCODE_OK)
    {
        cupdlp_printf("generate_dualOT_model_from_csv, exit_cleanup\n");
    }
    cupdlp_free(a);
    cupdlp_free(b);
    cupdlp_free(a_coarse);
    cupdlp_free(b_coarse);
    cupdlp_free(c_coarse);
    cupdlp_free(c_delete)
}
}

void generate_coarse_dualOT_model_delete_from_csv_longlong(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, long long *zero_idx, long long *zero_idx_len, long long len_after_delete, cupdlp_int coarse_degree)
{
    cupdlp_retcode retcode = RETCODE_OK;

    cupdlp_int resolution_coarse = resolution / pow(2, coarse_degree);
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

    cupdlp_int a_coarse_len = a_len / pow(pow(2, coarse_degree), 2);
    cupdlp_int b_coarse_len = b_len / pow(pow(2, coarse_degree), 2);
    cupdlp_float *a_coarse = cupdlp_NULL;
    CUPDLP_INIT(a_coarse, a_coarse_len);
    coarsingArray1D(a_coarse, a, resolution, coarse_degree);
    cupdlp_float *b_coarse = cupdlp_NULL;
    CUPDLP_INIT(b_coarse, b_coarse_len);
    coarsingArray1D(b_coarse, b, resolution, coarse_degree);

    // 构造keep数组，用于方便判断一个数是否需要保留
    cupdlp_bool *keep = cupdlp_NULL;
    long long *keep_idx = cupdlp_NULL;
    long long keep_len = pow(resolution_coarse, 4);
    CUPDLP_INIT(keep, keep_len);
    CUPDLP_INIT(keep_idx, keep_len);
    for (long long i = 0; i < keep_len; i++)
    {
        keep[i] = true;
        keep_idx[i] = 0;
    }
    for (long long i = 0; i < *zero_idx_len; i++)
    {
        keep[zero_idx[i]] = false;
    }

    cupdlp_float *c_coarse_delete = cupdlp_NULL;
    CUPDLP_INIT(c_coarse_delete, len_after_delete);
    cupdlp_printf("c_delete_len: %lld\n", len_after_delete);
    generate_c_coarse_delete_directly(c_coarse_delete, coarse_degree, resolution, resolution, keep, keep_idx);
    cupdlp_printf("开始构建模型\n");
    // 统计耗时
    cupdlp_float time_generate = getTimeStamp();
    // generate_dualOT_model_from_distribution_and_cost(model, a, b, a_len, b_len, c);
    // generate_dualOT_model_byMatrix_from_distribution_and_cost(model, a, b, a_len, b_len, c_delete, resolution);
    generate_dualOT_model_delete_byMatrix_from_distribution_and_cost_longlong(model, a_coarse, b_coarse, a_coarse_len, b_coarse_len, c_coarse_delete, resolution_coarse, zero_idx, zero_idx_len);
    time_generate = getTimeStamp() - time_generate;
    cupdlp_printf("模型构建完成, 耗时: %.3f 秒\n", time_generate);

exit_cleanup:

    if (retcode != RETCODE_OK)
    {
        cupdlp_printf("generate_dualOT_model_from_csv, exit_cleanup\n");
    }
    cupdlp_free(a);
    cupdlp_free(b);
    cupdlp_free(a_coarse);
    cupdlp_free(b_coarse);
    cupdlp_free(c_coarse_delete);
    cupdlp_free(keep);
    cupdlp_free(keep_idx);

    //     cupdlp_float *c_coarse = cupdlp_NULL;
    //     long long c_coarse_len = a_coarse_len * b_coarse_len;
    //     CUPDLP_INIT(c_coarse, c_coarse_len);
    //     coarsing_normalizedSquaredEuclideanDistance_longlong(c_coarse, resolution, resolution, coarse_degree);

    //     cupdlp_float *c_delete = cupdlp_NULL;
    //     CUPDLP_INIT(c_delete, len_after_delete);
    //     cupdlp_printf("c_delete_len: %lld\n", len_after_delete);
    //     deleteArray1DElements_longlong(c_delete, c_coarse, c_coarse_len, zero_idx, zero_idx_len);

    //     cupdlp_printf("开始构建模型\n");
    //     // 统计耗时
    //     cupdlp_float time_generate = getTimeStamp();
    //     // generate_dualOT_model_from_distribution_and_cost(model, a, b, a_len, b_len, c);
    //     // generate_dualOT_model_byMatrix_from_distribution_and_cost(model, a, b, a_len, b_len, c_delete, resolution);
    //     generate_dualOT_model_delete_byMatrix_from_distribution_and_cost_longlong(model, a_coarse, b_coarse, a_coarse_len, b_coarse_len, c_delete, resolution_coarse, zero_idx, zero_idx_len);
    //     time_generate = getTimeStamp() - time_generate;
    //     cupdlp_printf("模型构建完成, 耗时: %.3f 秒\n", time_generate);

    // exit_cleanup:
    // {
    //     if (retcode != RETCODE_OK)
    //     {
    //         cupdlp_printf("generate_dualOT_model_from_csv, exit_cleanup\n");
    //     }
    //     cupdlp_free(a);
    //     cupdlp_free(b);
    //     cupdlp_free(a_coarse);
    //     cupdlp_free(b_coarse);
    //     cupdlp_free(c_coarse);
    //     cupdlp_free(c_delete)
    // }
}

void generate_coarse_dualOT_model_delete_by_keep_from_csv_longlong(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_bool *keep, long long *keep_idx, long long len_after_delete, cupdlp_int coarse_degree)
{
    cupdlp_retcode retcode = RETCODE_OK;

    cupdlp_int resolution_coarse = resolution / pow(2, coarse_degree);
    // void *model = NULL;
    // model = createModel();
    cupdlp_printf("开始csv读取数据\n");
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

    cupdlp_float *c_coarse_delete = cupdlp_NULL;
    CUPDLP_INIT(c_coarse_delete, len_after_delete);
    cupdlp_printf("c_delete_len: %lld\n", len_after_delete);

    generate_c_coarse_delete_directly(c_coarse_delete, coarse_degree, resolution, resolution, keep, keep_idx);
    cupdlp_printf("开始构建模型\n");
    // 统计耗时
    cupdlp_float time_generate = getTimeStamp();
    // generate_dualOT_model_from_distribution_and_cost(model, a, b, a_len, b_len, c);
    // generate_dualOT_model_byMatrix_from_distribution_and_cost(model, a, b, a_len, b_len, c_delete, resolution);
    generate_dualOT_model_delete_by_keep_byMatrix_from_distribution_and_cost_longlong(model, a_coarse, b_coarse, a_coarse_len, b_coarse_len, c_coarse_delete, resolution_coarse, keep, keep_idx, &len_after_delete);
    time_generate = getTimeStamp() - time_generate;
    cupdlp_printf("模型构建完成, 耗时: %.3f 秒\n", time_generate);

exit_cleanup:

    if (retcode != RETCODE_OK)
    {
        cupdlp_printf("generate_dualOT_model_from_csv, exit_cleanup\n");
    }
    cupdlp_free(a);
    cupdlp_free(b);
    cupdlp_free(a_coarse);
    cupdlp_free(b_coarse);
    cupdlp_free(c_coarse_delete);

    //     cupdlp_float *c_coarse = cupdlp_NULL;
    //     long long c_coarse_len = a_coarse_len * b_coarse_len;
    //     CUPDLP_INIT(c_coarse, c_coarse_len);
    //     coarsing_normalizedSquaredEuclideanDistance_longlong(c_coarse, resolution, resolution, coarse_degree);

    //     cupdlp_float *c_delete = cupdlp_NULL;
    //     CUPDLP_INIT(c_delete, len_after_delete);
    //     cupdlp_printf("c_delete_len: %lld\n", len_after_delete);
    //     deleteArray1DElements_longlong(c_delete, c_coarse, c_coarse_len, zero_idx, zero_idx_len);

    //     cupdlp_printf("开始构建模型\n");
    //     // 统计耗时
    //     cupdlp_float time_generate = getTimeStamp();
    //     // generate_dualOT_model_from_distribution_and_cost(model, a, b, a_len, b_len, c);
    //     // generate_dualOT_model_byMatrix_from_distribution_and_cost(model, a, b, a_len, b_len, c_delete, resolution);
    //     generate_dualOT_model_delete_byMatrix_from_distribution_and_cost_longlong(model, a_coarse, b_coarse, a_coarse_len, b_coarse_len, c_delete, resolution_coarse, zero_idx, zero_idx_len);
    //     time_generate = getTimeStamp() - time_generate;
    //     cupdlp_printf("模型构建完成, 耗时: %.3f 秒\n", time_generate);

    // exit_cleanup:
    // {
    //     if (retcode != RETCODE_OK)
    //     {
    //         cupdlp_printf("generate_dualOT_model_from_csv, exit_cleanup\n");
    //     }
    //     cupdlp_free(a);
    //     cupdlp_free(b);
    //     cupdlp_free(a_coarse);
    //     cupdlp_free(b_coarse);
    //     cupdlp_free(c_coarse);
    //     cupdlp_free(c_delete)
    // }
}

// void LP_Solve_Multiscale(w, w_coarse, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx)
// {
//     LP_SolvePDHG(w_coarse, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx);
//     // 提取x_coarse, y_coarse, 转化为下一层的x_origin, y_origin
//     LP_SolvePDHG(w_coarse, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, x_origin, nCols_origin, y_origin, ifSaveSol, constraint_new_idx);
// }

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
        cupdlp_printf("Presolve\n");
        presolveinfo = createPresolve();
        presolvedmodel = presolvedModel(presolveinfo, model);
        model2solve = presolvedmodel;
    }
    presolve_time = getTimeStamp() - presolve_time;

    CUPDLP_CALL(formulateLP_new(model2solve, &cost, &nCols, &nRows, &nnz, &nEqs,
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

void createCUPDLPwork_clear(CUPDLPwork *w, void *model, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam, cupdlp_int **constraint_new_idx)
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
        cupdlp_printf("Presolve\n");
        presolveinfo = createPresolve();
        presolvedmodel = presolvedModel(presolveinfo, model);
        model2solve = presolvedmodel;
    }
    presolve_time = getTimeStamp() - presolve_time;

    CUPDLP_CALL(formulateLP_new(model2solve, &cost, &nCols, &nRows, &nnz, &nEqs,
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
exit_cleanup:
    if (retcode != RETCODE_OK)
    {
        cupdlp_printf("createCUPDLPwork exit_cleanup\n");
    }
}

void fine_dualOT_primal(cupdlp_float *x_init, cupdlp_float *x_coarse_solution, cupdlp_int x_len, cupdlp_int x_coarse_len, cupdlp_int resolution, cupdlp_int coarse_degree)
{
    // resolution是当前的分辨率， resolution_coarse是上一层次的分辨率
    cupdlp_int scale = pow(2, coarse_degree);
    cupdlp_int resolution_coarse = resolution / scale;
    if (x_len != x_coarse_len * pow(scale, 2))
    {
        cupdlp_printf("x_len != x_coarse_len * pow(scale, 2)\n");
        return;
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

void fine_dualOT_dual(cupdlp_float *y_init, cupdlp_float *y_coarse_solution, long long y_len, long long y_coarse_len, cupdlp_int resolution, cupdlp_int coarse_degree)
{
    long long scale = pow(2, coarse_degree);
    long long y_len_temp = y_len;
    long long y_coarse_len_temp = y_coarse_len;
    long long multiple_temp = pow(scale, 4);
    long long y_len_temp_compare = y_coarse_len_temp * multiple_temp;
    if (y_len_temp != y_len_temp_compare)
    {
        cupdlp_printf("y_len: %lld\n", y_len_temp);
        cupdlp_printf("y_coarse_len: %lld\n", y_coarse_len_temp);
        cupdlp_printf("scale: %lld\n", scale);
        cupdlp_printf("y_len != y_coarse_len * pow(scale, 4)\n");
        return;
    }
    cupdlp_printf("fine_dualOT_dual开始\n");
    cupdlp_float y_temp = 0.0;
    cupdlp_int resolution_coarse = resolution / scale;
    long long idx_coarse_temp = 0;
    long long idx_1 = 0;
    long long idx_2 = 0;
    long long idx_temp = 0;
    for (cupdlp_int i1 = 0; i1 < resolution_coarse; i1++)
    {
        cupdlp_printf("i1: %d\n", i1);
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

cupdlp_float *fine_and_delete_dualOT_dual_longlong(long long *len_after_delete, cupdlp_bool *keep, long long *keep_idx, cupdlp_float *y_coarse_solution, cupdlp_int resolution, cupdlp_int coarse_degree_now, cupdlp_int coarse_degree_last, cupdlp_float thr)
{
    cupdlp_bool retcode = RETCODE_OK;
    // 输出的keep代表的是当前的分辨率下哪些约束被delete了，哪些被保留了
    cupdlp_int coarse_degree_diff = coarse_degree_last - coarse_degree_now;
    long long scale = pow(2, coarse_degree_diff);
    cupdlp_printf("fine_and_delete_dualOT_dual_longlong开始\n");
    cupdlp_int resolution_coarse = resolution / pow(2, coarse_degree_last);
    long long idx_coarse_temp = 0;
    long long idx_coarse_col_temp = 0;
    cupdlp_bool *keep_coarse = cupdlp_NULL;
    long long keep_coarse_len = pow(resolution_coarse, 4);
    CUPDLP_INIT(keep_coarse, keep_coarse_len);
    long long *keep_coarse_idx = cupdlp_NULL;
    CUPDLP_INIT_ZERO(keep_coarse_idx, keep_coarse_len);
    // 用来记录非零的粗元素，在粗分辨率图像中是所在列的第几个非零元素
    long long *keep_coarse_col_idx = cupdlp_NULL;
    CUPDLP_INIT_ZERO(keep_coarse_col_idx, keep_coarse_len);
    cupdlp_float y_temp = 0.0;

    for (long long i = 0; i < keep_coarse_len; i++)
    {
        keep_coarse[i] = false;
    }
    long long true_idx = 0;
    // 考虑粗分辨率下(i1, j1)到(i2, j2)的输运，对应到细分辨率下的输运
    // 先跑出keep中有几个true
    for (cupdlp_int i1 = 0; i1 < resolution_coarse; i1++)
    {
        cupdlp_printf("i1: %d\n", i1);
        for (cupdlp_int j1 = 0; j1 < resolution_coarse; j1++)
        {
            for (cupdlp_int i2 = 0; i2 < resolution_coarse; i2++)
            {
                for (cupdlp_int j2 = 0; j2 < resolution_coarse; j2++)
                {
                    idx_coarse_temp = (i1 * resolution_coarse + j1) * pow(resolution_coarse, 2) + (i2 * resolution_coarse + j2);
                    y_temp = y_coarse_solution[idx_coarse_temp];
                    if (fabs(y_temp) >= thr)
                    {
                        keep_coarse[idx_coarse_temp] = true;
                        keep_coarse_idx[true_idx] = idx_coarse_temp;
                        keep_coarse_col_idx[true_idx] = idx_coarse_col_temp;
                        true_idx += 1;
                    }
                }
            }
        }
    }
    cupdlp_printf("keep_coarse完成, 其中有%lld个true\n", true_idx);

    long long y_coarse_delete_len = true_idx;
    long long y_init_delete_len = pow(scale, 4);
    y_init_delete_len = y_init_delete_len * y_coarse_delete_len;
    cupdlp_float *y_init_delete = cupdlp_NULL;
    CUPDLP_INIT_ZERO(y_init_delete, y_init_delete_len);
    *len_after_delete = y_init_delete_len;

    long long y_coarse_len = keep_coarse_len;

    // 根据keep_coarse修改y_init_delete
    cupdlp_printf("修改y_init_delete\n");
    long long idx_temp = 0;
    for (long long i = 0; i < y_coarse_len; i++)
    {
        if (keep_coarse[i])
        {
            y_temp = y_coarse_solution[i] / pow(scale, 4);
            for (int j = 0; j < pow(scale, 4); j++)
            {
                // idx_temp表示这是y_init_delete中的第几个元素
                idx_temp = keep_coarse_idx[i] * pow(scale, 4) + j;
                y_init_delete[idx_temp] = y_temp;
            }
        }
    }
    // 修改keep和keep_idx
    cupdlp_int resolution_now = resolution / pow(2, coarse_degree_now);
    long long idx_1 = 0;
    long long idx_2 = 0;
    true_idx = 0;
    for (cupdlp_int i1 = 0; i1 < resolution_coarse; i1++)
    {
        cupdlp_printf("i1: %d\n", i1);
        for (cupdlp_int j1 = 0; j1 < resolution_coarse; j1++)
        {
            for (cupdlp_int i2 = 0; i2 < resolution_coarse; i2++)
            {
                for (cupdlp_int j2 = 0; j2 < resolution_coarse; j2++)
                {
                    idx_coarse_temp = (i1 * resolution_coarse + j1) * pow(resolution_coarse, 2) + (i2 * resolution_coarse + j2);
                    if (keep_coarse[idx_coarse_temp])
                    {
                        for (cupdlp_int k1 = 0; k1 < scale; k1++)
                        {
                            for (cupdlp_int l1 = 0; l1 < scale; l1++)
                            {
                                for (cupdlp_int k2 = 0; k2 < scale; k2++)
                                {
                                    for (cupdlp_int l2 = 0; l2 < scale; l2++)
                                    {
                                        idx_1 = (i1 * scale + k1) * resolution_now + j1 * scale + l1;
                                        idx_2 = (i2 * scale + k2) * resolution_now + j2 * scale + l2;
                                        idx_temp = idx_1 * pow(resolution_now, 2) + idx_2;

                                        keep[idx_temp] = true;
                                        keep_idx[idx_temp] = true_idx;
                                        true_idx += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    cupdlp_printf("keep中有%lld个true\n", true_idx);
    // for (long long i = 0; i < *len_after_delete; i++)
    // {
    //     if (keep[i])
    //     {
    //         printf("keep_true_idx: %lld\n", keep_idx[i]);
    //     }
    // }
    cupdlp_printf("fine_and_delete_dualOT_dual_longlong完成\n");
    // print_float_array1D(y_init_delete, *len_after_delete);
exit_cleanup:
    if (retcode != RETCODE_OK)
    {
        cupdlp_printf("fine_and_delete_dualOT_dual exit_cleanup\n");
    }
    cupdlp_free(keep_coarse);
    cupdlp_free(keep_coarse_idx);
    cupdlp_free(keep_coarse_col_idx);
    return y_init_delete;
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
long long countArray1D_Smaller_than_threshold_longlong(cupdlp_float *a, long long a_len, cupdlp_float thr)
{
    long long count = 0;
    for (long long i = 0; i < a_len; i++)
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
    cupdlp_printf("countArray1D_Smaller_than_threshold_with_Record完成\n");
    return record;
}

long long *countArray1D_Smaller_than_threshold_with_Record_longlong(cupdlp_float *a, long long a_len, long long *a_record_len, cupdlp_float thr)
{

    long long count = 0;
    for (long long i = 0; i < a_len; i++)
    {
        if (fabs(a[i]) < thr)
            count += 1;
        if (count % 1000000 == 0)
            cupdlp_printf("count: %lld\n", count);
    }
    cupdlp_printf("零元素个数为：%lld\n", count);
    long long *record = (long long *)malloc(count * sizeof(long long));
    *a_record_len = count;
    count = 0;
    for (long long i = 0; i < a_len; i++)
    {
        if (fabs(a[i]) < thr)
        {
            record[count] = i;
            count++;
        }
    }
    cupdlp_printf("countArray1D_Smaller_than_threshold_with_Record完成, 零元素个数为：%lld\n", count);
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
    cupdlp_printf("保存完毕, 文件名：%s\n", filename);

    fclose(file);
}

void analyseArray1D(cupdlp_float *a, long long a_len, cupdlp_float thr, const char *filename)
{
    long long count = countArray1D_Smaller_than_threshold_longlong(a, a_len, thr);
    cupdlp_float a_len_float = a_len + 0.0;
    cupdlp_float propotion = 100 * count / a_len_float;
    cupdlp_printf("向量总长度为：%lld，小于阈值%.10f的元素个数为：%lld, 稀疏元素占比为：%.2f%%\n", a_len, thr, count, propotion);
    // saveArray1D_to_csv(a, a_len, filename);
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

// a为原数组，zero_idx为需要删除的元素的idx，a_dele为删除元素后的数组
void deleteArray1DElements(cupdlp_float *a_dele, cupdlp_float *a, cupdlp_int a_len, cupdlp_int *zero_idx, cupdlp_int *zero_idx_len)
{
    // 创建一个布尔数组，用于标记哪些元素需要被保留
    bool *keep = (bool *)malloc(a_len * sizeof(bool));
    for (int i = 0; i < a_len; i++)
    {
        keep[i] = true; // 初始时假设所有元素都保留
    }

    // 标记需要删除的元素
    for (int i = 0; i < *zero_idx_len; i++)
    {
        keep[zero_idx[i]] = false;
        // cupdlp_printf("需要删除的元素的idx：%d\n", zero_idx[i]);
    }

    int count = 0;
    for (int i = 0; i < a_len; i++)
    {
        if (keep[i])
        {
            a_dele[count] = a[i];
            count += 1;
            // cupdlp_printf("保留的元素：%f\n", a[i]);
        }
    }
    // 释放内存
    free(keep);
}

void deleteArray1DElements_longlong(cupdlp_float *a_dele, cupdlp_float *a, long long a_len, long long *zero_idx, long long *zero_idx_len)
{
    // 创建一个布尔数组，用于标记哪些元素需要被保留
    bool *keep = (bool *)malloc(a_len * sizeof(bool));
    for (long long i = 0; i < a_len; i++)
    {
        keep[i] = true; // 初始时假设所有元素都保留
    }

    // 标记需要删除的元素
    for (long long i = 0; i < *zero_idx_len; i++)
    {
        keep[zero_idx[i]] = false;
        // cupdlp_printf("需要删除的元素的idx：%d\n", zero_idx[i]);
    }

    long long count = 0;
    for (long long i = 0; i < a_len; i++)
    {
        if (keep[i])
        {
            a_dele[count] = a[i];
            count += 1;
            // cupdlp_printf("保留的元素：%f\n", a[i]);
        }
    }
    // 释放内存
    free(keep);
    keep = NULL;
}

void deleteArray1DElements_by_keep_longlong(cupdlp_float *a_delete, cupdlp_float *a, long long a_len, cupdlp_bool *keep)
{
    long long count = 0;
    for (long long i = 0; i < a_len; i++)
    {
        if (keep[i])
        {
            a_delete[count] = a[i];
            count += 1;
            // cupdlp_printf("保留的元素：%f\n", a[i]);
        }
    }
}
void recoverArray1DElements(cupdlp_float *a_recover, cupdlp_float *a_delete, cupdlp_int a_recover_len, cupdlp_int a_delete_len, cupdlp_int *zero_idx, cupdlp_int *zero_idx_len)
{
    // a_delete为删除元素后的数组，a_recover为恢复元素后的数组
    // zero_idx为被删除的元素的idx
    int count1 = 0;
    int count2 = 0;
    for (int i = 0; i < a_recover_len; i++)
    {
        if (i == zero_idx[count1])
        {
            a_recover[i] = 0.0;
            count1 += 1;
        }
        else
        {
            a_recover[i] = a_delete[count2];
            count2 += 1;
        }
    }
}

void recoverArray1DElements_longlong(cupdlp_float *a_recover, cupdlp_float *a_delete, long long a_recover_len, long long a_delete_len, long long *zero_idx, long long *zero_idx_len)
{
    // a_delete为删除元素后的数组，a_recover为恢复元素后的数组
    // zero_idx为被删除的元素的idx
    long long count1 = 0;
    long long count2 = 0;
    for (long long i = 0; i < a_recover_len; i++)
    {
        if (i == zero_idx[count1])
        {
            a_recover[i] = 0.0;
            count1 += 1;
        }
        else
        {
            a_recover[i] = a_delete[count2];
            count2 += 1;
        }
    }
}

void recoverArray1DElements_by_keep_longlong(cupdlp_float *a_recover, cupdlp_float *a_delete, long long a_recover_len, long long a_delete_len, cupdlp_bool *keep, long long *keep_idx)
{
    // a_delete为删除元素后的数组，a_recover为恢复元素后的数组
    // zero_idx为被删除的元素的idx
    long long count1 = 0;
    long long count2 = 0;
    long long count = 0;
    for (long long i = 0; i < a_delete_len; i++)
    {
        if (keep[i])
        {
            a_recover[i] = a_delete[keep_idx[i]];
        }
        else
        {
            a_recover[i] = 0.0;
        }
    }
    // for (long long i = 0; i < a_recover_len; i++)
    // {
    //     if (i == zero_idx[count1])
    //     {
    //         a_recover[i] = 0.0;
    //         count1 += 1;
    //     }
    //     else
    //     {
    //         a_recover[i] = a_delete[count2];
    //         count2 += 1;
    //     }
    // }
}

void testPtr(cupdlp_float *a)
{
    a[0] = 99999.9;
    cupdlp_printf("a[0] = %f\n", a[0]);
}

void construct_and_solve_Multiscale(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_bool *ifChangeIntParam, cupdlp_bool *ifChangeFloatParam, cupdlp_int *intParam, cupdlp_float *floatParam, cupdlp_bool ifSaveSol, cupdlp_int coarse_degree, cupdlp_int coarse_degree_last, cupdlp_float **x_solution, cupdlp_float **y_solution, cupdlp_float *x_solution_last, cupdlp_float *y_solution_last)
{
    cupdlp_bool retcode = RETCODE_OK;
    if (coarse_degree_last != -1)
    {
        cupdlp_printf("coarse_degree_last != -1, 利用稀疏性构造模型\n");
        // 对上一步的解进行细化
        int resolution_now = resolution / pow(2, coarse_degree);
        int resolution_last = resolution / pow(2, coarse_degree_last);
        int x_init_len = 2 * pow(resolution_now, 2);
        // int y_init_len = pow(resolution_now, 4);
        long long y_init_len = pow(resolution_now, 4);
        int x_solution_last_len = 2 * pow(resolution_last, 2);
        long long y_solution_last_len = pow(resolution_last, 4);
        int coarse_degree_diff = coarse_degree_last - coarse_degree;
        cupdlp_float *x_init = cupdlp_NULL;
        cupdlp_float *y_init = cupdlp_NULL;
        CUPDLP_INIT_ZERO(x_init, x_init_len);
        CUPDLP_INIT_ZERO(y_init, y_init_len);
        fine_dualOT_primal(x_init, x_solution_last, x_init_len, x_solution_last_len, resolution_now, coarse_degree_diff);
        fine_dualOT_dual(y_init, y_solution_last, y_init_len, y_solution_last_len, resolution_now, coarse_degree_diff);
        // 找出y_init中小于阈值的元素
        long long *y_init_zero_idx_len = cupdlp_NULL;
        CUPDLP_INIT_ZERO(y_init_zero_idx_len, 1);
        long long *y_init_zero_idx = countArray1D_Smaller_than_threshold_with_Record_longlong(y_init, y_init_len, y_init_zero_idx_len, 1e-20);
        long long y_init_delete_len = y_init_len - *y_init_zero_idx_len;
        cupdlp_printf("y_init_len: %d\n", y_init_len);
        cupdlp_printf("y_init_zero_idx_len: %d\n", *y_init_zero_idx_len);
        cupdlp_printf("y_init_delete_len: %d\n", y_init_delete_len);
        cupdlp_float *y_init_delete = cupdlp_NULL;
        CUPDLP_INIT_ZERO(y_init_delete, y_init_delete_len);
        deleteArray1DElements_longlong(y_init_delete, y_init, y_init_len, y_init_zero_idx, y_init_zero_idx_len);
        // 创建y_solution_delete，用于暂时储存delete过后的解
        cupdlp_float *y_solution_delete = cupdlp_NULL;
        CUPDLP_INIT_ZERO(y_solution_delete, y_init_delete_len);
        // 创建模型
        // void *model = NULL;
        // model = createModel();
        generate_coarse_dualOT_model_delete_from_csv_longlong(model, csvpath_1, csvpath_2, resolution, y_init_zero_idx, y_init_zero_idx_len, y_init_delete_len, coarse_degree);
        int *nCols_origin_ptr = cupdlp_NULL;
        int *nRows_ptr = cupdlp_NULL;
        CUPDLP_INIT_ZERO(nCols_origin_ptr, 1);
        CUPDLP_INIT_ZERO(nRows_ptr, 1);
        cupdlp_int *constraint_new_idx = NULL;
        CUPDLPwork *w = cupdlp_NULL;
        CUPDLP_INIT_ZERO(w, 1);
        createCUPDLPwork(w, model, ifChangeIntParam, intParam, &constraint_new_idx, nCols_origin_ptr, nRows_ptr);
        // 求解
        char fout[256];
        cupdlp_printf("开始求解\n");
        sprintf(fout, "./solution_Resolution%d_CoarseDegree%d.txt", resolution, coarse_degree);
        cupdlp_bool whether_first_fine = true; // 默认使用cuPDLP自带的初始步长和权重

        CUPDLP_CALL(LP_SolvePDHG_Multiscale(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, ifSaveSol, constraint_new_idx, x_solution, &y_solution_delete, x_init, y_init_delete));
        cupdlp_printf("利用稀疏性构造模型求解完毕\n");
        // 求解完后，把y_delete_solution恢复到y_solution中
        recoverArray1DElements_longlong(*y_solution, y_solution_delete, y_init_len, y_init_delete_len, y_init_zero_idx, y_init_zero_idx_len);
        // cupdlp_free(nCols_origin_ptr);
        // cupdlp_free(nRows_ptr);
        // cupdlp_free(constraint_new_idx);
        // // PDHG_Destroy(w);
        // cupdlp_free(x_init);
        // cupdlp_free(y_init);
        // cupdlp_free(y_init_zero_idx_len);
        // cupdlp_free(y_init_zero_idx);
        // cupdlp_free(y_init_delete);
        // cupdlp_free(stepsize_last);
        // cupdlp_free(weight_last);
        // cupdlp_free(stepsize_init);
        // cupdlp_free(weight_init);
        // cupdlp_free(y_solution_delete);
    }
    else
    {
        cupdlp_printf("coarse_degree_last == -1, 直接构造模型\n");

        // 创建模型
        // void *model = NULL;
        // model = createModel();
        generate_coarse_dualOT_model_from_csv(model, csvpath_1, csvpath_2, resolution, coarse_degree);
        int *nCols_origin_ptr = cupdlp_NULL;
        int *nRows_ptr = cupdlp_NULL;
        CUPDLP_INIT_ZERO(nCols_origin_ptr, 1);
        CUPDLP_INIT_ZERO(nRows_ptr, 1);
        cupdlp_int *constraint_new_idx = NULL;
        CUPDLPwork *w = cupdlp_NULL;
        CUPDLP_INIT_ZERO(w, 1);
        printf("开始构建CUPDLPwork...\n");
        createCUPDLPwork(w, model, ifChangeIntParam, intParam, &constraint_new_idx, nCols_origin_ptr, nRows_ptr);
        // 直接用初始值
        int resolution_now = resolution / pow(2, coarse_degree);
        int x_init_len = 2 * pow(resolution_now, 2);
        int y_init_len = pow(resolution_now, 4);
        int coarse_degree_diff = coarse_degree_last - coarse_degree;
        cupdlp_float *x_init = cupdlp_NULL;
        cupdlp_float *y_init = cupdlp_NULL;
        CUPDLP_INIT_ZERO(x_init, x_init_len);
        CUPDLP_INIT_ZERO(y_init, y_init_len);
        *x_init = *x_solution_last;
        *y_init = *y_solution_last;
        // 求解
        char fout[256];
        sprintf(fout, "./solution_Resolution%d_CoarseDegree%d.txt", resolution, coarse_degree);
        cupdlp_bool whether_first_fine = true; // 默认使用cuPDLP自带的初始步长和权重
        // 传入的x_solution和y_solution是二级指针，所以不用再写成&x_solution, &y_solution
        CUPDLP_CALL(LP_SolvePDHG_Multiscale(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, ifSaveSol, constraint_new_idx, x_solution, y_solution, x_init, y_init));

        // cupdlp_free(nCols_origin_ptr);
        // cupdlp_free(nRows_ptr);
        // cupdlp_free(constraint_new_idx);
        // // PDHG_Destroy(w);
        // cupdlp_free(x_init);
        // cupdlp_free(y_init);
        // cupdlp_free(stepsize_last);
        // cupdlp_free(weight_last);
        // cupdlp_free(stepsize_init);
        // cupdlp_free(weight_init);
    }

exit_cleanup:
{
    if (retcode != RETCODE_OK)
    {
        cupdlp_printf("construct_and_solve_Multiscale, exit_cleanup\n");
    }
}
}

void construct_and_solve_Multiscale_longlong(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_bool *ifChangeIntParam, cupdlp_bool *ifChangeFloatParam, cupdlp_int *intParam, cupdlp_float *floatParam, cupdlp_bool ifSaveSol, cupdlp_int coarse_degree, cupdlp_int coarse_degree_last, cupdlp_float **x_solution, cupdlp_float **y_solution, cupdlp_float *x_solution_last, cupdlp_float *y_solution_last)
{
    cupdlp_bool retcode = RETCODE_OK;
    if (coarse_degree_last != -1)
    {
        cupdlp_printf("coarse_degree_last != -1, 利用稀疏性构造模型\n");
        // 对上一步的解进行细化
        int resolution_now = resolution / pow(2, coarse_degree);
        int resolution_last = resolution / pow(2, coarse_degree_last);
        int x_init_len = 2 * pow(resolution_now, 2);
        // int y_init_len = pow(resolution_now, 4);
        long long y_init_len = pow(resolution_now, 4);
        int x_solution_last_len = 2 * pow(resolution_last, 2);
        long long y_solution_last_len = pow(resolution_last, 4);
        int coarse_degree_diff = coarse_degree_last - coarse_degree;
        cupdlp_float *x_init = cupdlp_NULL;
        CUPDLP_INIT_ZERO(x_init, x_init_len);
        fine_dualOT_primal(x_init, x_solution_last, x_init_len, x_solution_last_len, resolution_now, coarse_degree_diff);

        // cupdlp_bool *keep = cupdlp_NULL;
        // long long *keep_idx = cupdlp_NULL;
        // CUPDLP_INIT(keep, y_init_len);
        // CUPDLP_INIT(keep_idx, y_init_len);
        // for (long long i = 0; i < y_init_len; i++)
        // {
        //     keep[i] = false;
        //     keep_idx[i] = 0;
        // }
        // long long *len_after_delete = cupdlp_NULL;
        // CUPDLP_INIT_ZERO(len_after_delete, 1);
        // cupdlp_float *y_init_delete = fine_and_delete_dualOT_dual_longlong(len_after_delete, keep, keep_idx, y_solution_last, resolution, coarse_degree, coarse_degree_last, 1e-20);
        // generate_coarse_dualOT_model_delete_by_keep_from_csv_longlong(model, csvpath_1, csvpath_2, resolution, keep, keep_idx, *len_after_delete, coarse_degree);
        // cupdlp_int *constraint_new_idx = NULL;
        // CUPDLPwork *w = cupdlp_NULL;
        // CUPDLP_INIT_ZERO(w, 1);
        // createCUPDLPwork_clear(w, model, ifChangeIntParam, intParam, &constraint_new_idx);
        // // 求解
        // char fout[256];
        // cupdlp_printf("开始求解\n");
        // sprintf(fout, "./solution_Resolution%d_CoarseDegree%d.txt", resolution, coarse_degree);
        // cupdlp_bool whether_first_fine = true; // 默认使用cuPDLP自带的初始步长和权重

        // // 创建y_solution_delete，用于暂时储存delete过后的解
        // cupdlp_float *y_solution_delete = cupdlp_NULL;
        // CUPDLP_INIT_ZERO(y_solution_delete, *len_after_delete);
        // CUPDLP_CALL(LP_SolvePDHG_Multiscale(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, ifSaveSol, constraint_new_idx, x_solution, &y_solution_delete, x_init, y_init_delete));
        // cupdlp_printf("利用稀疏性构造模型求解完毕\n");
        // // 求解完后，把y_delete_solution恢复到y_solution中
        // recoverArray1DElements_by_keep_longlong(*y_solution, y_solution_delete, y_init_len, *len_after_delete, keep, keep_idx);

        // cupdlp_free(constraint_new_idx);
        // cupdlp_free(w);
        // // // PDHG_Destroy(w);
        // cupdlp_free(keep);
        // cupdlp_free(keep_idx);
        // cupdlp_free(len_after_delete);
        // cupdlp_free(x_init);
        // cupdlp_free(y_init_delete);

        cupdlp_float *y_init = cupdlp_NULL;
        CUPDLP_INIT_ZERO(y_init, y_init_len);
        fine_dualOT_dual(y_init, y_solution_last, y_init_len, y_solution_last_len, resolution_now, coarse_degree_diff);
        // 找出y_init中小于阈值的元素
        long long *y_init_zero_idx_len = cupdlp_NULL;
        CUPDLP_INIT_ZERO(y_init_zero_idx_len, 1);
        long long *y_init_zero_idx = countArray1D_Smaller_than_threshold_with_Record_longlong(y_init, y_init_len, y_init_zero_idx_len, 1e-20);
        long long y_init_delete_len = y_init_len - *y_init_zero_idx_len;
        cupdlp_printf("y_init_len: %lld\n", y_init_len);
        cupdlp_printf("y_init_zero_idx_len: %lld\n", *y_init_zero_idx_len);
        cupdlp_printf("y_init_delete_len: %lld\n", y_init_delete_len);
        cupdlp_float *y_init_delete = cupdlp_NULL;
        CUPDLP_INIT_ZERO(y_init_delete, y_init_delete_len);
        deleteArray1DElements_longlong(y_init_delete, y_init, y_init_len, y_init_zero_idx, y_init_zero_idx_len);
        // 创建y_solution_delete，用于暂时储存delete过后的解
        cupdlp_float *y_solution_delete = cupdlp_NULL;
        CUPDLP_INIT_ZERO(y_solution_delete, y_init_delete_len);

        // 创建模型
        // void *model = NULL;
        // model = createModel();
        generate_coarse_dualOT_model_delete_from_csv_longlong(model, csvpath_1, csvpath_2, resolution, y_init_zero_idx, y_init_zero_idx_len, y_init_delete_len, coarse_degree);
        cupdlp_int *constraint_new_idx = NULL;
        CUPDLPwork *w = cupdlp_NULL;
        CUPDLP_INIT_ZERO(w, 1);
        createCUPDLPwork_clear(w, model, ifChangeIntParam, intParam, &constraint_new_idx);
        // 求解
        char fout[256];
        cupdlp_printf("开始求解\n");
        sprintf(fout, "./solution_Resolution%d_CoarseDegree%d.txt", resolution, coarse_degree);
        cupdlp_bool whether_first_fine = true; // 默认使用cuPDLP自带的初始步长和权重

        CUPDLP_CALL(LP_SolvePDHG_Multiscale(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, ifSaveSol, constraint_new_idx, x_solution, &y_solution_delete, x_init, y_init_delete));
        cupdlp_printf("利用稀疏性构造模型求解完毕\n");
        // 求解完后，把y_delete_solution恢复到y_solution中
        recoverArray1DElements_longlong(*y_solution, y_solution_delete, y_init_len, y_init_delete_len, y_init_zero_idx, y_init_zero_idx_len);

        cupdlp_free(constraint_new_idx);
        cupdlp_free(w);
        // // PDHG_Destroy(w);
        cupdlp_free(x_init);
        cupdlp_free(y_init);
        cupdlp_free(y_init_zero_idx_len);
        cupdlp_free(y_init_zero_idx);
        cupdlp_free(y_init_delete);
    }
    else
    {
        cupdlp_printf("coarse_degree_last == -1, 直接构造模型\n");

        // 创建模型
        // void *model = NULL;
        // model = createModel();
        generate_coarse_dualOT_model_from_csv(model, csvpath_1, csvpath_2, resolution, coarse_degree);
        cupdlp_int *constraint_new_idx = NULL;
        CUPDLPwork *w = cupdlp_NULL;
        CUPDLP_INIT_ZERO(w, 1);
        printf("开始构建CUPDLPwork...\n");
        createCUPDLPwork_clear(w, model, ifChangeIntParam, intParam, &constraint_new_idx);
        // 直接用初始值
        int resolution_now = resolution / pow(2, coarse_degree);
        int x_init_len = 2 * pow(resolution_now, 2);
        int y_init_len = pow(resolution_now, 4);
        int coarse_degree_diff = coarse_degree_last - coarse_degree;
        cupdlp_float *x_init = cupdlp_NULL;
        cupdlp_float *y_init = cupdlp_NULL;
        CUPDLP_INIT_ZERO(x_init, x_init_len);
        CUPDLP_INIT_ZERO(y_init, y_init_len);
        *x_init = *x_solution_last;
        *y_init = *y_solution_last;
        // 求解
        char fout[256];
        sprintf(fout, "./solution_Resolution%d_CoarseDegree%d.txt", resolution, coarse_degree);
        cupdlp_bool whether_first_fine = true; // 默认使用cuPDLP自带的初始步长和权重
        // 传入的x_solution和y_solution是二级指针，所以不用再写成&x_solution, &y_solution
        CUPDLP_CALL(LP_SolvePDHG_Multiscale(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, ifSaveSol, constraint_new_idx, x_solution, y_solution, x_init, y_init));

        cupdlp_free(constraint_new_idx);
        cupdlp_free(w);
        // PDHG_Destroy(w);
        cupdlp_free(x_init);
        cupdlp_free(y_init);
    }

exit_cleanup:
{
    if (retcode != RETCODE_OK)
    {
        cupdlp_printf("construct_and_solve_Multiscale, exit_cleanup\n");
    }
}
}

void construct_and_solve_Multiscale_by_keep_longlong(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_bool *ifChangeIntParam, cupdlp_bool *ifChangeFloatParam, cupdlp_int *intParam, cupdlp_float *floatParam, cupdlp_bool ifSaveSol, cupdlp_int coarse_degree, cupdlp_int coarse_degree_last, cupdlp_float **x_solution, cupdlp_float **y_solution, cupdlp_float *x_solution_last, cupdlp_float *y_solution_last)
{
    cupdlp_bool retcode = RETCODE_OK;
    if (coarse_degree_last != -1)
    {
        cupdlp_printf("coarse_degree_last != -1, 利用稀疏性构造模型\n");
        // 对上一步的解进行细化
        int resolution_now = resolution / pow(2, coarse_degree);
        int resolution_last = resolution / pow(2, coarse_degree_last);
        int x_init_len = 2 * pow(resolution_now, 2);
        // int y_init_len = pow(resolution_now, 4);
        long long y_init_len = pow(resolution_now, 4);
        int x_solution_last_len = 2 * pow(resolution_last, 2);
        long long y_solution_last_len = pow(resolution_last, 4);
        int coarse_degree_diff = coarse_degree_last - coarse_degree;
        cupdlp_float *x_init = cupdlp_NULL;
        CUPDLP_INIT_ZERO(x_init, x_init_len);
        fine_dualOT_primal(x_init, x_solution_last, x_init_len, x_solution_last_len, resolution_now, coarse_degree_diff);

        cupdlp_bool *keep = cupdlp_NULL;
        long long *keep_idx = cupdlp_NULL;
        CUPDLP_INIT(keep, y_init_len);
        CUPDLP_INIT(keep_idx, y_init_len);
        for (long long i = 0; i < y_init_len; i++)
        {
            keep[i] = false;
            keep_idx[i] = 0;
        }
        long long *len_after_delete = cupdlp_NULL;
        CUPDLP_INIT_ZERO(len_after_delete, 1);
        cupdlp_float *y_init_delete = fine_and_delete_dualOT_dual_longlong(len_after_delete, keep, keep_idx, y_solution_last, resolution, coarse_degree, coarse_degree_last, 1e-20);
        generate_coarse_dualOT_model_delete_by_keep_from_csv_longlong(model, csvpath_1, csvpath_2, resolution, keep, keep_idx, *len_after_delete, coarse_degree);
        cupdlp_int *constraint_new_idx = NULL;
        CUPDLPwork *w = cupdlp_NULL;
        CUPDLP_INIT_ZERO(w, 1);
        createCUPDLPwork_clear(w, model, ifChangeIntParam, intParam, &constraint_new_idx);
        // 求解
        char fout[256];
        cupdlp_printf("开始求解\n");
        sprintf(fout, "./solution_Resolution%d_CoarseDegree%d.txt", resolution, coarse_degree);
        cupdlp_bool whether_first_fine = true; // 默认使用cuPDLP自带的初始步长和权重

        // 创建y_solution_delete，用于暂时储存delete过后的解
        cupdlp_float *y_solution_delete = cupdlp_NULL;
        CUPDLP_INIT_ZERO(y_solution_delete, *len_after_delete);
        CUPDLP_CALL(LP_SolvePDHG_Multiscale(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, ifSaveSol, constraint_new_idx, x_solution, &y_solution_delete, x_init, y_init_delete));
        cupdlp_printf("利用稀疏性构造模型求解完毕\n");
        // 求解完后，把y_delete_solution恢复到y_solution中
        recoverArray1DElements_by_keep_longlong(*y_solution, y_solution_delete, y_init_len, *len_after_delete, keep, keep_idx);

        cupdlp_free(constraint_new_idx);
        cupdlp_free(w);
        // // PDHG_Destroy(w);
        cupdlp_free(keep);
        cupdlp_free(keep_idx);
        cupdlp_free(len_after_delete);
        cupdlp_free(x_init);
        cupdlp_free(y_init_delete);

        // fine_dualOT_dual(y_init, y_solution_last, y_init_len, y_solution_last_len, resolution_now, coarse_degree_diff);
        // // 找出y_init中小于阈值的元素
        // long long *y_init_zero_idx_len = cupdlp_NULL;
        // CUPDLP_INIT_ZERO(y_init_zero_idx_len, 1);
        // long long *y_init_zero_idx = countArray1D_Smaller_than_threshold_with_Record_longlong(y_init, y_init_len, y_init_zero_idx_len, 1e-20);
        // long long y_init_delete_len = y_init_len - *y_init_zero_idx_len;
        // cupdlp_printf("y_init_len: %lld\n", y_init_len);
        // cupdlp_printf("y_init_zero_idx_len: %lld\n", *y_init_zero_idx_len);
        // cupdlp_printf("y_init_delete_len: %lld\n", y_init_delete_len);
        // cupdlp_float *y_init_delete = cupdlp_NULL;
        // CUPDLP_INIT_ZERO(y_init_delete, y_init_delete_len);
        // deleteArray1DElements_longlong(y_init_delete, y_init, y_init_len, y_init_zero_idx, y_init_zero_idx_len);
        // // 创建y_solution_delete，用于暂时储存delete过后的解
        // cupdlp_float *y_solution_delete = cupdlp_NULL;
        // CUPDLP_INIT_ZERO(y_solution_delete, y_init_delete_len);

        // // 创建模型
        // // void *model = NULL;
        // // model = createModel();
        // generate_coarse_dualOT_model_delete_from_csv_longlong(model, csvpath_1, csvpath_2, resolution, y_init_zero_idx, y_init_zero_idx_len, y_init_delete_len, coarse_degree);
        // cupdlp_int *constraint_new_idx = NULL;
        // CUPDLPwork *w = cupdlp_NULL;
        // CUPDLP_INIT_ZERO(w, 1);
        // createCUPDLPwork_clear(w, model, ifChangeIntParam, intParam, &constraint_new_idx);
        // // 求解
        // char fout[256];
        // cupdlp_printf("开始求解\n");
        // sprintf(fout, "./solution_Resolution%d_CoarseDegree%d.txt", resolution, coarse_degree);
        // cupdlp_bool whether_first_fine = true; // 默认使用cuPDLP自带的初始步长和权重

        // CUPDLP_CALL(LP_SolvePDHG_Multiscale(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, ifSaveSol, constraint_new_idx, x_solution, &y_solution_delete, x_init, y_init_delete));
        // cupdlp_printf("利用稀疏性构造模型求解完毕\n");
        // // 求解完后，把y_delete_solution恢复到y_solution中
        // recoverArray1DElements_longlong(*y_solution, y_solution_delete, y_init_len, y_init_delete_len, y_init_zero_idx, y_init_zero_idx_len);

        // cupdlp_free(constraint_new_idx);
        // cupdlp_free(w);
        // // // PDHG_Destroy(w);
        // cupdlp_free(x_init);
        // cupdlp_free(y_init);
        // cupdlp_free(y_init_zero_idx_len);
        // cupdlp_free(y_init_zero_idx);
        // cupdlp_free(y_init_delete);
    }
    else
    {
        cupdlp_printf("coarse_degree_last == -1, 直接构造模型\n");

        // 创建模型
        // void *model = NULL;
        // model = createModel();
        generate_coarse_dualOT_model_from_csv(model, csvpath_1, csvpath_2, resolution, coarse_degree);
        cupdlp_int *constraint_new_idx = NULL;
        CUPDLPwork *w = cupdlp_NULL;
        CUPDLP_INIT_ZERO(w, 1);
        printf("开始构建CUPDLPwork...\n");
        createCUPDLPwork_clear(w, model, ifChangeIntParam, intParam, &constraint_new_idx);
        // 直接用初始值
        int resolution_now = resolution / pow(2, coarse_degree);
        int x_init_len = 2 * pow(resolution_now, 2);
        int y_init_len = pow(resolution_now, 4);
        int coarse_degree_diff = coarse_degree_last - coarse_degree;
        cupdlp_float *x_init = cupdlp_NULL;
        cupdlp_float *y_init = cupdlp_NULL;
        CUPDLP_INIT_ZERO(x_init, x_init_len);
        CUPDLP_INIT_ZERO(y_init, y_init_len);
        *x_init = *x_solution_last;
        *y_init = *y_solution_last;
        // 求解
        char fout[256];
        sprintf(fout, "./solution_Resolution%d_CoarseDegree%d.txt", resolution, coarse_degree);
        cupdlp_bool whether_first_fine = true; // 默认使用cuPDLP自带的初始步长和权重
        // 传入的x_solution和y_solution是二级指针，所以不用再写成&x_solution, &y_solution
        CUPDLP_CALL(LP_SolvePDHG_Multiscale(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, ifSaveSol, constraint_new_idx, x_solution, y_solution, x_init, y_init));

        cupdlp_free(constraint_new_idx);
        cupdlp_free(w);
        // PDHG_Destroy(w);
        cupdlp_free(x_init);
        cupdlp_free(y_init);
    }

exit_cleanup:
{
    if (retcode != RETCODE_OK)
    {
        cupdlp_printf("construct_and_solve_Multiscale, exit_cleanup\n");
    }
}
}

// void construct_and_solve_Multiscale_withoutRecover(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_bool *ifChangeIntParam, cupdlp_bool *ifChangeFloatParam, cupdlp_int *intParam, cupdlp_float *floatParam, cupdlp_bool ifSaveSol, cupdlp_int coarse_degree, cupdlp_int coarse_degree_last, cupdlp_float **x_solution, cupdlp_float **y_solution, cupdlp_float *x_solution_last, cupdlp_float *y_solution_last)
// {
//     cupdlp_bool retcode = RETCODE_OK;
//     if (coarse_degree_last != -1)
//     {
//         cupdlp_printf("coarse_degree_last != -1, 利用稀疏性构造模型\n");
//         // 对上一步的解进行细化
//         int resolution_now = resolution / pow(2, coarse_degree);
//         int resolution_last = resolution / pow(2, coarse_degree_last);
//         int x_init_len = 2 * pow(resolution_now, 2);
//         // int y_init_len = pow(resolution_now, 4);
//         long long y_init_len = pow(resolution_now, 4);
//         int x_solution_last_len = 2 * pow(resolution_last, 2);
//         long long y_solution_last_len = pow(resolution_last, 4);
//         int coarse_degree_diff = coarse_degree_last - coarse_degree;
//         cupdlp_float *x_init = cupdlp_NULL;
//         cupdlp_float *y_init = cupdlp_NULL;
//         CUPDLP_INIT_ZERO(x_init, x_init_len);
//         CUPDLP_INIT_ZERO(y_init, y_init_len);
//         fine_dualOT_primal(x_init, x_solution_last, x_init_len, x_solution_last_len, resolution_now, coarse_degree_diff);
//         fine_dualOT_dual(y_init, y_solution_last, y_init_len, y_solution_last_len, resolution_now, coarse_degree_diff);
//         // 找出y_init中小于阈值的元素
//         long long *y_init_zero_idx_len = cupdlp_NULL;
//         CUPDLP_INIT_ZERO(y_init_zero_idx_len, 1);
//         long long *y_init_zero_idx = countArray1D_Smaller_than_threshold_with_Record(y_init, y_init_len, y_init_zero_idx_len, 1e-20);
//         long long y_init_delete_len = y_init_len - *y_init_zero_idx_len;
//         cupdlp_printf("y_init_len: %lld\n", y_init_len);
//         cupdlp_printf("y_init_zero_idx_len: %lld\n", *y_init_zero_idx_len);
//         cupdlp_printf("y_init_delete_len: %lld\n", y_init_delete_len);
//         cupdlp_float *y_init_delete = cupdlp_NULL;
//         CUPDLP_INIT_ZERO(y_init_delete, y_init_delete_len);
//         deleteArray1DElements(y_init_delete, y_init, y_init_len, y_init_zero_idx, y_init_zero_idx_len);
//         // 创建y_solution_delete，用于暂时储存delete过后的解
//         cupdlp_float *y_solution_delete = cupdlp_NULL;
//         CUPDLP_INIT_ZERO(y_solution_delete, y_init_delete_len);
//         // 创建模型
//         // void *model = NULL;
//         // model = createModel();
//         generate_coarse_dualOT_model_delete_from_csv(model, csvpath_1, csvpath_2, resolution, y_init_zero_idx, y_init_zero_idx_len, y_init_delete_len, coarse_degree);
//         int *nCols_origin_ptr = cupdlp_NULL;
//         int *nRows_ptr = cupdlp_NULL;
//         CUPDLP_INIT_ZERO(nCols_origin_ptr, 1);
//         CUPDLP_INIT_ZERO(nRows_ptr, 1);
//         cupdlp_int *constraint_new_idx = NULL;
//         CUPDLPwork *w = cupdlp_NULL;
//         CUPDLP_INIT_ZERO(w, 1);
//         createCUPDLPwork(w, model, ifChangeIntParam, intParam, &constraint_new_idx, nCols_origin_ptr, nRows_ptr);
//         // 求解
//         char fout[256];
//         cupdlp_printf("开始求解\n");
//         sprintf(fout, "./solution_Resolution%d_CoarseDegree%d.txt", resolution, coarse_degree);
//         cupdlp_bool whether_first_fine = true; // 默认使用cuPDLP自带的初始步长和权重

//         CUPDLP_CALL(LP_SolvePDHG_Multiscale(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, ifSaveSol, constraint_new_idx, x_solution, &y_solution_delete, x_init, y_init_delete));
//         cupdlp_printf("利用稀疏性构造模型求解完毕\n");
//         // // 求解完后，把y_delete_solution恢复到y_solution中
//         // recoverArray1DElements(*y_solution, y_solution_delete, y_init_len, y_init_delete_len, y_init_zero_idx, y_init_zero_idx_len);

//         // cupdlp_free(nCols_origin_ptr);
//         // cupdlp_free(nRows_ptr);
//         // cupdlp_free(constraint_new_idx);
//         // // PDHG_Destroy(w);
//         // cupdlp_free(x_init);
//         // cupdlp_free(y_init);
//         // cupdlp_free(y_init_zero_idx_len);
//         // cupdlp_free(y_init_zero_idx);
//         // cupdlp_free(y_init_delete);
//         // cupdlp_free(stepsize_last);
//         // cupdlp_free(weight_last);
//         // cupdlp_free(stepsize_init);
//         // cupdlp_free(weight_init);
//         // cupdlp_free(y_solution_delete);
//     }
//     else
//     {
//         cupdlp_printf("coarse_degree_last == -1, 直接构造模型\n");

//         // 创建模型
//         // void *model = NULL;
//         // model = createModel();
//         generate_coarse_dualOT_model_from_csv(model, csvpath_1, csvpath_2, resolution, coarse_degree);
//         int *nCols_origin_ptr = cupdlp_NULL;
//         int *nRows_ptr = cupdlp_NULL;
//         CUPDLP_INIT_ZERO(nCols_origin_ptr, 1);
//         CUPDLP_INIT_ZERO(nRows_ptr, 1);
//         cupdlp_int *constraint_new_idx = NULL;
//         CUPDLPwork *w = cupdlp_NULL;
//         CUPDLP_INIT_ZERO(w, 1);
//         printf("开始构建CUPDLPwork...\n");
//         createCUPDLPwork(w, model, ifChangeIntParam, intParam, &constraint_new_idx, nCols_origin_ptr, nRows_ptr);
//         // 直接用初始值
//         int resolution_now = resolution / pow(2, coarse_degree);
//         int x_init_len = 2 * pow(resolution_now, 2);
//         int y_init_len = pow(resolution_now, 4);
//         int coarse_degree_diff = coarse_degree_last - coarse_degree;
//         cupdlp_float *x_init = cupdlp_NULL;
//         cupdlp_float *y_init = cupdlp_NULL;
//         CUPDLP_INIT_ZERO(x_init, x_init_len);
//         CUPDLP_INIT_ZERO(y_init, y_init_len);
//         *x_init = *x_solution_last;
//         *y_init = *y_solution_last;
//         // 求解
//         char fout[256];
//         sprintf(fout, "./solution_Resolution%d_CoarseDegree%d.txt", resolution, coarse_degree);
//         cupdlp_bool whether_first_fine = true; // 默认使用cuPDLP自带的初始步长和权重
//         // 传入的x_solution和y_solution是二级指针，所以不用再写成&x_solution, &y_solution
//         CUPDLP_CALL(LP_SolvePDHG_Multiscale(w, ifChangeIntParam, intParam, ifChangeFloatParam, floatParam, fout, ifSaveSol, constraint_new_idx, x_solution, y_solution, x_init, y_init));

//         // cupdlp_free(nCols_origin_ptr);
//         // cupdlp_free(nRows_ptr);
//         // cupdlp_free(constraint_new_idx);
//         // // PDHG_Destroy(w);
//         // cupdlp_free(x_init);
//         // cupdlp_free(y_init);
//         // cupdlp_free(stepsize_last);
//         // cupdlp_free(weight_last);
//         // cupdlp_free(stepsize_init);
//         // cupdlp_free(weight_init);
//     }

// exit_cleanup:
// {
//     if (retcode != RETCODE_OK)
//     {
//         cupdlp_printf("construct_and_solve_Multiscale, exit_cleanup\n");
//     }
// }
// }

void construct_and_solve_Multiscale_test(void *model, const char *csvpath_1, const char *csvpath_2, cupdlp_int resolution, cupdlp_int coarse_degree)
{
    // void *model = NULL;
    // model = createModel();
    generate_coarse_dualOT_model_from_csv(model, csvpath_1, csvpath_2, resolution, coarse_degree);
}
