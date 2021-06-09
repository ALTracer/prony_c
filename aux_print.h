#ifndef AUX_PRINT_H
#define AUX_PRINT_H

#include "types.h"

void print_vector(const double*, size_t);
void print_matrix(myArray_t, size_t, size_t);
void print_vector_c(const _Dcomplex *, size_t);
void print_matrix_c(myZArray_t M, size_t mla, size_t mlb);
void print_alpha_freq(const _Dcomplex * z, size_t len);

#endif // AUX_PRINT_H
