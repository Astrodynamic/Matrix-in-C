#ifndef MATRIC_INCLUDE_MATRIX_H_
#define MATRIC_INCLUDE_MATRIX_H_

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define SUCCESS 1
#define FAILURE 0

#define EPS 1e-7

typedef enum { F_OK, F_EM, F_EC } e_error_t;

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

int castom_create_matrix(int rows, int columns, matrix_t *result);
void castom_remove_matrix(matrix_t *A);
int castom_eq_matrix(matrix_t *A, matrix_t *B);
int castom_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int castom_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int castom_mult_number(matrix_t *A, double number, matrix_t *result);
int castom_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int castom_transpose(matrix_t *A, matrix_t *result);
int castom_calc_complements(matrix_t *A, matrix_t *result);
matrix_t castom_minor_of_matrix(matrix_t *A, int row, int collumn);
int castom_determinant(matrix_t *A, double *result);
int castom_inverse_matrix(matrix_t *A, matrix_t *result);

#endif  // MATRIC_INCLUDE_MATRIX_H_
