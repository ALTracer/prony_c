#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "config.h"
#include "types.h"

void print_vector(const double * vector, size_t len){
    size_t k;
    for (k = 0; k<len; k++){
        printf("%.2f\t",vector[k]);
        if (k%10 == 9)
            printf("\n");
    }
    printf("\n");
    return;
}

void print_matrix(myArray_t M, size_t mla, size_t mlb){
    size_t a,b;
    for (a = 0; a < mla; a++){
        for (b = 0; b < mlb; b++){
            printf("%.2f\t", M[a][b]);
        }
        printf("\n");
    }
    printf("\n");
    return;
}

void print_vector_c(const _Dcomplex * vector, size_t len){
    size_t k;
    for (k = 0; k<len; k++){
        printf("%.4f%+.4fi\n", creal(vector[k]), cimag(vector[k]));
//        if (k%10 == 9) printf("\n");
    }
    printf("\n");
    return;
}

void print_matrix_c(myZArray_t M, size_t mla, size_t mlb){
    size_t a,b;
    for (a = 0; a < mla; a++){
        for (b = 0; b < mlb; b++){
            printf("%.3f%+.3fi\t", creal(M[a][b]), cimag(M[a][b]));
        }
        printf("\n");
    }
    printf("\n");
    return;
}

void print_alpha_freq(const _Dcomplex * z, size_t len){
    size_t k;
    double alpha, freq;
    for (k=0; k<len; k++){
        alpha = log(cabs(z[k]))/Ts;
        freq = atan(cimag(z[k])/creal(z[k]))/(2*M_PI*Ts);
        printf("[%lld]: alpha = %.3f, freq = %.3f\n", k, alpha, freq);
    }
    return;
}
