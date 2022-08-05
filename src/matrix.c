#include "matrix.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Include SSE intrinsics
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#include <x86intrin.h>
#endif

/* Below are some intel intrinsics that might be useful
 * void _mm256_storeu_pd (double * mem_addr, __m256d a)
 * __m256d _mm256_set1_pd (double a)
 * __m256d _mm256_set_pd (double e3, double e2, double e1, double e0)
 * __m256d _mm256_loadu_pd (double const * mem_addr)
 * __m256d _mm256_add_pd (__m256d a, __m256d b)
 * __m256d _mm256_sub_pd (__m256d a, __m256d b)
 * __m256d _mm256_fmadd_pd (__m256d a, __m256d b, __m256d c)
 * __m256d _mm256_mul_pd (__m256d a, __m256d b)
 * __m256d _mm256_cmp_pd (__m256d a, __m256d b, const int imm8)
 * __m256d _mm256_and_pd (__m256d a, __m256d b)
 * __m256d _mm256_max_pd (__m256d a, __m256d b)
*/

/* Generates a random double between low and high */
double rand_double(double low, double high) {
    double range = (high - low);
    double div = RAND_MAX / range;
    return low + (rand() / div);
}

/* Generates a random matrix */
void rand_matrix(matrix *result, unsigned int seed, double low, double high) {
    srand(seed);
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            set(result, i, j, rand_double(low, high));
        }
    }
}

/*
 * Returns the double value of the matrix at the given row and column.
 * You may assume `row` and `col` are valid. Note that the matrix is in row-major order.
 */
double get(matrix *mat, int row, int col) {
    // Task 1.1 TODO
    return mat->data[(mat->cols * row) + col];
}

/*
 * Sets the value at the given row and column to val. You may assume `row` and
 * `col` are valid. Note that the matrix is in row-major order.
 */
void set(matrix *mat, int row, int col, double val) {
    // Task 1.1 TODO
    mat->data[(mat->cols * row) + col] = val;
}

/*
 * Allocates space for a matrix struct pointed to by the double pointer mat with
 * `rows` rows and `cols` columns. You should also allocate memory for the data array
 * and initialize all entries to be zeros. `parent` should be set to NULL to indicate that
 * this matrix is not a slice. You should also set `ref_cnt` to 1.
 * You should return -1 if either `rows` or `cols` or both have invalid values. Return -2 if any
 * call to allocate memory in this function fails.
 * Return 0 upon success.
 */
int allocate_matrix(matrix **mat, int rows, int cols) {
    // Task 1.2 TODO
    // HINTS: Follow these steps.
    // 1. Check if the dimensions are valid. Return -1 if either dimension is not positive.
    // 2. Allocate space for the new matrix struct. Return -2 if allocating memory failed.
    // 3. Allocate space for the matrix data, initializing all entries to be 0. Return -2 if allocating memory failed.
    // 4. Set the number of rows and columns in the matrix struct according to the arguments provided.
    // 5. Set the `parent` field to NULL, since this matrix was not created from a slice.
    // 6. Set the `ref_cnt` field to 1.
    // 7. Store the address of the allocated matrix struct at the location `mat` is pointing at.
    // 8. Return 0 upon success.
    int mat_size = rows * cols;
    if (mat_size <= 0) {
        return -1;
    } 

    // Allocation and throwing error if allocation fails
    matrix* matrix_ptr = (matrix*) calloc(1, sizeof(matrix));
    double* data_ptr = (double*) calloc((size_t)mat_size, sizeof(double));
    if (matrix_ptr == NULL || data_ptr == NULL){
        return -2;
    }

    matrix_ptr -> rows = rows;
    matrix_ptr -> cols = cols;
    matrix_ptr -> data = data_ptr;
    matrix_ptr -> parent = NULL;
    matrix_ptr -> ref_cnt = 1;
    
    * mat = matrix_ptr;

    //Success
    return 0;
}

/*
 * You need to make sure that you only free `mat->data` if `mat` is not a slice and has no existing slices,
 * or that you free `mat->parent->data` if `mat` is the last existing slice of its parent matrix and its parent
 * matrix has no other references (including itself).
 */
void deallocate_matrix(matrix *mat) {
    // Task 1.3 TODO
    // HINTS: Follow these steps.
    // 1. If the matrix pointer `mat` is NULL, return.
    // 2. If `mat` has no parent: decrement its `ref_cnt` field by 1. If the `ref_cnt` field becomes 0, then free `mat` and its `data` field.
    // 3. Otherwise, recursively call `deallocate_matrix` on `mat`'s parent, then free `mat`.
    if (mat == NULL) {
        return;
    }
    if (mat->parent == NULL) {
        mat->ref_cnt = mat->ref_cnt - 1;
        if (mat->ref_cnt == 0) {
            free(mat->data);
            free(mat);
        }
    }
    else {
        deallocate_matrix(mat->parent);
        free(mat);
    }
    
}

/*
 * Allocates space for a matrix struct pointed to by `mat` with `rows` rows and `cols` columns.
 * Its data should point to the `offset`th entry of `from`'s data (you do not need to allocate memory)
 * for the data field. `parent` should be set to `from` to indicate this matrix is a slice of `from`
 * and the reference counter for `from` should be incremented. Lastly, do not forget to set the
 * matrix's row and column values as well.
 * You should return -1 if either `rows` or `cols` or both have invalid values. Return -2 if any
 * call to allocate memory in this function fails.
 * Return 0 upon success.
 * NOTE: Here we're allocating a matrix struct that refers to already allocated data, so
 * there is no need to allocate space for matrix data.
 */
int allocate_matrix_ref(matrix **mat, matrix *from, int offset, int rows, int cols) {
    // Task 1.4 TODO
    // HINTS: Follow these steps.
    // 1. Check if the dimensions are valid. Return -1 if either dimension is not positive.
    // 2. Allocate space for the new matrix struct. Return -2 if allocating memory failed.
    // 3. Set the `data` field of the new struct to be the `data` field of the `from` struct plus `offset`.
    // 4. Set the number of rows and columns in the new struct according to the arguments provided.
    // 5. Set the `parent` field of the new struct to the `from` struct pointer.
    // 6. Increment the `ref_cnt` field of the `from` struct by 1.
    // 7. Store the address of the allocated matrix struct at the location `mat` is pointing at.
    // 8. Return 0 upon success.
    int mat_size = rows * cols;
    if (mat_size <= 0) {
        return -1;
    }

    matrix* new_mat = (matrix*) calloc(1, sizeof(matrix));
    if (new_mat == NULL){
        return -2;
    }

    // #3
    new_mat->data = from->data + offset;
    
    // #4
    new_mat->rows = rows;
    new_mat->cols = cols;
    new_mat->parent = from;
    from->ref_cnt ++;

    *mat = new_mat;

    return 0;
}

/*
 * Sets all entries in mat to val. Note that the matrix is in row-major order.
 */
void fill_matrix(matrix *mat, double val) {
    // Task 1.5 TODO
    double* data = mat->data;
    
    __m256d vals = _mm256_set1_pd (val);
    int size = mat->rows * mat->cols;
    int i;
    for (i = 0; i < size / 4 * 4; i+=4 ){
        _mm256_storeu_pd(data + i, vals);
    }
    
    //tail case
    for(; i < size; i++){
        data[i] = val;
    }

    
}

/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success.
 * Note that the matrix is in row-major order.
 */
int abs_matrix(matrix *result, matrix *mat) {
    // Task 1.5 TODO
    __m256d data_vec, mask, masked_data;
    double* data = mat->data;
    int i;
    int size = mat->rows * mat-> cols;
    __m256d zeroes = _mm256_set1_pd (0);
    __m256d neg_one = _mm256_set1_pd (-1);
    for (i = 0; i < size / 4 * 4; i+=4) {
        data_vec = _mm256_loadu_pd(data+i);
        mask = _mm256_cmp_pd (data_vec, zeroes, 1);
        masked_data = _mm256_and_pd(mask, data_vec);

        // masked_data * -2  + data
        data_vec = _mm256_fmadd_pd (masked_data, neg_one, data_vec);
        data_vec = _mm256_fmadd_pd (masked_data, neg_one, data_vec);

        //storing back to memory
        _mm256_storeu_pd(result->data + i, data_vec);
    }
    
    //tail case
    for(; i < size; i++) {
        result->data[i] = fabs(mat->data[i]);
    }
    
    return 0;
}

/*
 * (OPTIONAL)
 * Store the result of element-wise negating mat's entries to `result`.
 * Return 0 upon success.
 * Note that the matrix is in row-major order.
 */
int neg_matrix(matrix *result, matrix *mat) {
    // Task 1.5 TODO
    int mat_size = mat->rows * mat-> cols;
    double* r_data = result->data;
    double* mat_data = mat->data;
    for (int i = 0; i < mat_size; i++) {
        r_data[i] = -(mat_data[i]);
    }
    return 0;
}

/*
 * Store the result of adding mat1 and mat2 to `result`.
 * Return 0 upon success.
 * You may assume `mat1` and `mat2` have the same dimensions.
 * Note that the matrix is in row-major order.
 */
int add_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    // Task 1.5 TODO

    double* r_data = result->data;
    double* mat1_data = mat1->data;
    double* mat2_data = mat2->data;


    __m256d data1_vec, data2_vec, sum;
    int i;
    int size = mat1->rows * mat1-> cols;
    for (i = 0; i < size / 4 * 4; i+=4) {
        data1_vec = _mm256_loadu_pd(mat1_data + i);
        data2_vec = _mm256_loadu_pd(mat2_data + i);
        sum = _mm256_add_pd(data1_vec, data2_vec);
        _mm256_storeu_pd(result->data + i, sum);
    }

    //tail case
    for (; i < size; i++) {
        r_data[i] = mat1_data[i] + mat2_data[i];
    }
    return 0;
}

/*
 * (OPTIONAL)
 * Store the result of subtracting mat2 from mat1 to `result`.
 * Return 0 upon success.
 * You may assume `mat1` and `mat2` have the same dimensions.
 * Note that the matrix is in row-major order.
 */
int sub_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    // Task 1.5 TODO
    int mat_size = mat1->rows * mat1-> cols;
    double* r_data = result->data;
    double* mat1_data = mat1->data;
    double* mat2_data = mat2->data;
    for (int i = 0; i < mat_size; i++) {
        r_data[i] = mat1_data[i] - mat2_data[i];
    }
    return 0;
}

/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success.
 * Remember that matrix multiplication is not the same as multiplying individual elements.
 * You may assume `mat1`'s number of columns is equal to `mat2`'s number of rows.
 * Note that the matrix is in row-major order.
 */
int mul_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    // Task 1.6 TODO
    __m256d sum, data1_vec, data2_vec;
    double data2_0, data2_1, data2_2, data2_3, data_1, data_2, summ;
    int i = 0;
    int j = 0;
    int k = 0;

    for(i = 0; i < mat1->rows / 4 * 4; i+= 4) {
        for(j = 0; j < mat2->cols / 4 * 4; j+= 4) {
            sum = _mm256_set1_pd(0.0);
            for(k = 0; k < mat2->rows /4 * 4; k+=4) {
                data1_vec = _mm256_loadu_pd(mat1->data + (k + (i * mat1->cols)));
                data2_0 = get(mat2, k, j);
                data2_1 = get(mat2, k+1, j);
                data2_2 = get(mat2, k+2, j);
                data2_3 = get(mat2, k+3, j);
                data2_vec = _mm256_set_pd(data2_3, data2_2, data2_1, data2_0);
                
                sum = _mm256_fmadd_pd(data1_vec, data2_vec, sum);
            }
            // Accumulating the sums
            double dim_sum[4];
            _mm256_storeu_pd(dim_sum, sum);
            double sum_acc = dim_sum[0] + dim_sum[1] + dim_sum[2] + dim_sum[3];

            // Tail Case for K
            for(; k < mat2->rows; k++) {
                data_1 = get(mat1, i, k);
                data_2 = get(mat2, k, j);
                sum_acc += data_1 * data_2;
            }
            set(result, i, j, sum_acc);
        }
    }

    // Tail Case for i and j

    for(; i < mat1->rows; i++) {
        for(j = 0; j < mat2->cols; j++) {
            sum = _mm256_set1_pd(0.0);
            for(k = 0; k < mat2->rows /4 * 4; k+=4) {
                data1_vec = _mm256_loadu_pd(mat1->data + (k + (i * mat1->cols)));
                data2_0 = get(mat2, k, j);
                data2_1 = get(mat2, k+1, j);
                data2_2 = get(mat2, k+2, j);
                data2_3 = get(mat2, k+3, j);
                data2_vec = _mm256_set_pd(data2_3, data2_2, data2_1, data2_0);
                
                sum = _mm256_fmadd_pd(data1_vec, data2_vec, sum);
            }
            // Accumulating the sums
            double dim_sum[4];
            _mm256_storeu_pd(dim_sum, sum);
            double sum_acc = dim_sum[0] + dim_sum[1] + dim_sum[2] + dim_sum[3];

            // Tail Case for K
            for(; k < mat2->rows; k++) {
                data_1 = get(mat1, i, k);
                data_2 = get(mat2, k, j);
                sum_acc += data_1 * data_2;
            }
            set(result, i, j, sum_acc);
        }
    }
    return 0;


}

/*
 * Store the result of raising mat to the (pow)th power to `result`.
 * Return 0 upon success.
 * Remember that pow is defined with matrix multiplication, not element-wise multiplication.
 * You may assume `mat` is a square matrix and `pow` is a non-negative integer.
 * Note that the matrix is in row-major order.
 */
int pow_matrix(matrix *result, matrix *mat, int pow) {

    // Task 1.6 TODO
    matrix** temp_array = (matrix**) calloc(1, sizeof(matrix*));
    allocate_matrix(temp_array,result->rows, result->cols);
    matrix * temp = *temp_array;
    fill_matrix(result, 0.0);
    for (int i = 0; i < mat -> rows; i ++) {
        set(result, i, i, 1.0);
        set(temp, i, i, 1.0);
    }

    int size = mat->rows*mat->cols;
    if (pow != 0) {
        for (int i = 0; i < pow; i++) {
            mul_matrix(result, temp, mat);
            for (int j = 0; j < size; j++){
                temp->data[j] = result->data[j];
            }
        }
    }

    //Assigning results
    deallocate_matrix(temp);
    free(temp_array);
    return 0;
}

