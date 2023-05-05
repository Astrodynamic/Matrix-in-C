#include "matrix.h"

int castom_create_matrix(int rows, int columns, matrix_t *result) {
  int flag = F_EM;
  memset(result, 0, sizeof(matrix_t));
  if (rows > 0 && columns > 0) {
    flag = F_OK;
    result->rows = rows;
    result->columns = columns;
    result->matrix = (double **)calloc(rows, sizeof(double *));
    for (int i = 0; i < rows; ++i)
      result->matrix[i] = (double *)calloc(columns, sizeof(double));
  }
  return flag;
}

void castom_remove_matrix(matrix_t *A) {
  if (A) {
    if (A->matrix && A->rows && A->columns) {
      for (int i = 0; i < A->rows; ++i)
        if (A->matrix[i]) free(A->matrix[i]);
      free(A->matrix);
    }
    memset(A, 0, sizeof(matrix_t));
  }
}

int castom_eq_matrix(matrix_t *A, matrix_t *B) {
  int result = FAILURE;
  int A_r = A->rows, A_c = A->columns;
  int B_r = B->rows, B_c = B->columns;
  if (A_r == B_r && A_c == B_c && A_r && A_c) {
    result = SUCCESS;
    for (int i = 0; i < A->rows; ++i) {
      for (int j = 0; j < B->columns; ++j) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) >= EPS) result = FAILURE;
      }
    }
  }
  return result;
}

int castom_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int flag = F_OK;
  int A_r = A->rows, A_c = A->columns;
  int B_r = B->rows, B_c = B->columns;
  memset(result, 0, sizeof(matrix_t));
  if (A_r <= 0 || A_c <= 0 || B_r <= 0 || B_c <= 0) {
    flag = F_EM;
  } else if (A_r != B_r || A_c != B_c) {
    flag = F_EC;
  } else {
    castom_create_matrix(A_r, B_c, result);
    for (int i = 0; i < result->rows; ++i) {
      for (int j = 0; j < result->columns; ++j) {
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
      }
    }
  }
  return flag;
}

int castom_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int flag = F_OK;
  int A_r = A->rows, A_c = A->columns;
  int B_r = B->rows, B_c = B->columns;
  memset(result, 0, sizeof(matrix_t));
  if (A_r <= 0 || A_c <= 0 || B_r <= 0 || B_c <= 0) {
    flag = F_EM;
  } else if (A_r != B_r || A_c != B_c) {
    flag = F_EC;
  } else {
    castom_create_matrix(A_r, B_c, result);
    for (int i = 0; i < result->rows; ++i) {
      for (int j = 0; j < result->columns; ++j) {
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
      }
    }
  }
  return flag;
}

int castom_mult_number(matrix_t *A, double number, matrix_t *result) {
  int flag = F_EM;
  memset(result, 0, sizeof(matrix_t));
  if (A->rows > 0 && A->columns > 0) {
    flag = F_OK;
    castom_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < result->rows; ++i) {
      for (int j = 0; j < result->columns; ++j) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  }
  return flag;
}

int castom_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int flag = F_OK;
  int A_r = A->rows, A_c = A->columns;
  int B_r = B->rows, B_c = B->columns;
  memset(result, 0, sizeof(matrix_t));
  if (A_r <= 0 || A_c <= 0 || B_r <= 0 || B_c <= 0) {
    flag = F_EM;
  } else if (A_c != B_r) {
    flag = F_EC;
  } else {
    castom_create_matrix(A_r, B_c, result);
    for (int i = 0; i < result->rows; ++i) {
      for (int j = 0; j < result->columns; ++j) {
        result->matrix[i][j] = 0.0;
        for (int k = 0; k < A->columns; ++k) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    }
  }
  return flag;
}

int castom_transpose(matrix_t *A, matrix_t *result) {
  int flag = F_EM;
  memset(result, 0, sizeof(matrix_t));
  if (A->rows > 0 && A->columns > 0) {
    flag = F_OK;
    castom_create_matrix(A->columns, A->rows, result);
    for (int i = 0; i < A->rows; ++i) {
      for (int j = 0; j < A->columns; ++j) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  }
  return flag;
}

int castom_calc_complements(matrix_t *A, matrix_t *result) {
  int flag = F_OK;
  memset(result, 0, sizeof(matrix_t));
  int A_r = A->rows, A_c = A->columns;
  if (A_r <= 0 || A_c <= 0) {
    flag = F_EM;
  } else if (A_r != A_c) {
    flag = F_EC;
  } else {
    double det = 0;
    matrix_t minor;
    castom_create_matrix(A_r, A_c, result);
    if (A_r == 1) {
      result->matrix[0][0] = 1.0;
    } else {
      for (int i = 0; i < result->rows; ++i) {
        for (int j = 0; j < result->columns; ++j) {
          minor = castom_minor_of_matrix(A, i, j);
          castom_determinant(&minor, &det);
          result->matrix[i][j] = pow(-1.0, i + j) * det;
          castom_remove_matrix(&minor);
        }
      }
    }
  }
  return flag;
}

matrix_t castom_minor_of_matrix(matrix_t *A, int row, int collumn) {
  matrix_t t_matrix;
  castom_create_matrix(A->rows - 1, A->columns - 1, &t_matrix);
  for (int i = 0, t_i = 0; i < A->rows; ++i) {
    if (i != row) {
      for (int j = 0, t_j = 0; j < A->columns; ++j) {
        if (j != collumn) t_matrix.matrix[t_i][t_j++] = A->matrix[i][j];
      }
      ++t_i;
    }
  }
  return t_matrix;
}

int castom_determinant(matrix_t *A, double *result) {
  int flag = F_OK;
  if (A->rows <= 0 || A->columns <= 0) {
    flag = F_EM;
  } else if (A->rows != A->columns) {
    flag = F_EC;
  } else {
    *result = 0.0;
    if (A->rows == 1) {
      *result = A->matrix[0][0];
    } else if (A->rows == 2) {
      *result =
          A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
    } else {
      double det = 0;
      for (int i = 0; i < A->rows; ++i) {
        matrix_t minor = castom_minor_of_matrix(A, i, 0);
        castom_determinant(&minor, &det);
        *result += A->matrix[i][0] * pow(-1.0, i) * det;
        castom_remove_matrix(&minor);
      }
    }
  }
  return flag;
}

int castom_inverse_matrix(matrix_t *A, matrix_t *result) {
  int flag = F_OK;
  double determinant = 0;
  memset(result, 0, sizeof(matrix_t));
  flag = castom_determinant(A, &determinant);
  if (fabs(determinant) < EPS && flag == F_OK) {
    flag = F_EC;
  } else if (flag == F_OK) {
    matrix_t comp_temp;
    castom_calc_complements(A, &comp_temp);
    matrix_t tran_temp;
    castom_transpose(&comp_temp, &tran_temp);
    castom_remove_matrix(&comp_temp);
    castom_mult_number(&tran_temp, 1.0 / determinant, result);
    castom_remove_matrix(&tran_temp);
  }
  return flag;
}
