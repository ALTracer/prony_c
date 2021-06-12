#include <stdio.h> //printf
#include <stdlib.h> //malloc,calloc,free
#include <assert.h> //assert

#define _USE_MATH_DEFINES
#include <math.h> //sin, M_PI

#define __STDC_WANT_LIB_EXT1__ 1
#include <string.h> //memcpy_s

#if __has_include(<mkl_lapacke.h>)
    #include <mkl_lapacke.h>
    #define USE_LAPACKE
#endif

// Signal storage
double * signal;
size_t len;

#include "config.h"
#include "types.h"
#include "aux_print.h"

void compose_test_signal(double* signal, size_t len);

int lsolve(const myArray_t * a, const double * b, double * x);
int polyroots(const double * polycoeffs, size_t order, _Dcomplex * roots);
int lsolvez(const myZArray_t * a, const double * b, _Dcomplex * x);

void compose_matrix1(double * signal, size_t len, myArray_t M, size_t mla, size_t mlb){
    size_t a, b, idx;
    for (a = 0; a < mla; a++) {
        for (b = 0; b < mlb; b++) {
            idx = M_ORDER - b + a - 1;
            assert(idx < len);
            assert(idx >= 0);
            M[a][b] = signal[idx];
        }
    }
    printf("Matrix 1 (%lld x %lld):\n", mla, mlb);
    print_matrix(M, mla, mlb);
    return;
}
void compose_vector1(double * signal, size_t len, double * rem, size_t rem_len) {
    size_t k, idx;
    for (k = 0; k<rem_len; k++) {
        idx = M_ORDER + k;
        assert(idx < len);
        rem[k] = -signal[idx];
    }
    printf("Vector 1 (-remainder):\n");
    print_vector(rem, rem_len);
    return;
}

void compose_matrix2(const _Dcomplex * z, myZArray_t Z) {
    size_t a, b;
    _Dcomplex temp, e;
    for (b = 0; b < M_ORDER; b++){
        Z[0][b] = _Cbuild(1.0F, 0.0F);
    }
    for (a = 1; a < M_ORDER; a++){
        for (b = 0; b < M_ORDER; b++){
            e = _Cbuild(a*1.0F, 0.0F);
            temp = cpow(z[b], e);
            Z[a][b] = temp;
        }
    }
    printf("Matrix 2 (%lld x %lld):\n", M_ORDER, M_ORDER);
    print_matrix_c(Z, M_ORDER, M_ORDER);
}

int main()
{
    // Create input signal
    len = 40;
    signal = calloc(len,sizeof(double));

    compose_test_signal(signal, len);
    printf("Input signal (%lld points):\n", len);
    print_vector(signal, len);

/*    for (size_t i = 0; i < len-M_ORDER; i++) { */
    // Create first matrix
    assert(len >= MLA + MLB);
    myArray_t matrix1;
    compose_matrix1(signal, len, matrix1, MLA, MLB);

    // Create first remainder vector
    size_t rem_len = MLA;
    double * rem = calloc(rem_len, sizeof(double));
    compose_vector1(signal, len, rem, rem_len);

    // Solve matrix equation for polynomial
    double * a_solved = calloc(M_ORDER, sizeof(double));
    lsolve(&matrix1, rem, a_solved);

    double * a = calloc(M_ORDER+1, sizeof(double));
    a[0] = 1.0f;
    memcpy(&a[1], a_solved, M_ORDER*sizeof(double));
    printf("Vector 2*: a (%lld order polynomial):\n", M_ORDER);
    print_vector(a, M_ORDER+1);

    // Solve polynomial for complex roots
    _Dcomplex * z = calloc(M_ORDER, sizeof(_Dcomplex));
    double a1[M_ORDER+1] = {1.0, 1.353, 1.028, -2.808, 3.942};
    polyroots(a1, M_ORDER, z);
    // Optionally print the alpha and Freq.
    print_alpha_freq(z, M_ORDER);

    // Create second matrix (complex)
    myZArray_t Z;
    compose_matrix2(z, Z);
    double * rem2 = calloc(M_ORDER, sizeof(double));
    memcpy(rem2, signal, M_ORDER*sizeof(double));
    //Solve complex matrix equation for complex polynomial
    _Dcomplex * h = calloc(M_ORDER-1, sizeof(_Dcomplex));
    lsolvez(&Z, rem2, h);
    printf("Vector 4: h\n");
    print_vector_c(h, M_ORDER-1);
    //Optionally print the A and Psi

    // Clean up dynamic storage
    free(h);
    free(rem2);
    free(z);
    free(a);
    free(a_solved);
    free(rem);
/*    memcpy_s(signal, len-1, signal, len-1);
    realloc(signal, len-1);
    len--;
    }*/
    free(signal);
    return 0;
}

void compose_test_signal(double* signal, size_t len){
    size_t k;
    double A = 1.0f;
    double w = 2*M_PI*50;
    double Phi = M_PI/2;
    for(k = 0; k<len; k++){
        signal[k] = A*sin(w*k*Ts + Phi);
//        signal[k] = k;
    }
    return;
}

int lsolve(const myArray_t * a, const double * b, double * x){
#ifdef USE_LAPACKE
//#include "mkl_lapacke.h"
    lapack_int info,m,n,lda,ldb,nrhs;
    nrhs = 1;
    m = MLA; n = MLB;
    lda = MLB; ldb = 1;

    double * b_mod = calloc(MLA, sizeof(double));
    const size_t alen = MLA*sizeof(double);
    memcpy_s(b_mod, alen, b, alen);
    info = LAPACKE_dgels(LAPACK_ROW_MAJOR,'N',m,n,nrhs,
                         (double*)a,lda,
                         b_mod,ldb);
    if (info)
        printf("LAPACK solver error at info = %d\n", info);

    const size_t blen = MLB*sizeof(double);
    memcpy_s(x, blen, b_mod, blen);
    free(b_mod);
    return info;
#else
    for (int k = 0; k < M_ORDER; k++)
        x[k] = pow(-1.0f,k+1);
    printf("Lapacke unavailable, skipping lsolve!\n");
    return -1;
#endif
}

int polyroots(const double * polycoeffs, size_t order, _Dcomplex * roots) {
#include "tclague.h"
    int status;
    double eps = 1e-8, eps2 = 1e-6;
    Complex * A = calloc(order+1, sizeof(Complex));
    for (size_t k = 0; k<order; k++) {
        A[k][0] = polycoeffs[k];
        A[k][1] = 0.0f;
    }
    Complex * wtmp1 = calloc(order+1, sizeof(Complex));
    Complex * wtmp2 = calloc(order+1, sizeof(Complex));
    Complex * wtmp3 = calloc(order+1, sizeof(Complex));
    Complex * Z = calloc(order+1, sizeof(Complex));
    Z[0][0] = 0.8; Z[0][1] = 0.6;
    CLAGUE((int)order, A, 10, eps, eps2, &status, Z, wtmp1, wtmp2, wtmp3);
    if (status)
        printf("CLAGUE did not converge under 1e-6.\n");
    free(A);
    free(wtmp1);
    free(wtmp2);
    free(wtmp3);
    for (size_t k = 0; k < order; k++) {
//        roots[k] = Z[k][0] + Z[k][1]*I;
        roots[k] = _Cbuild(Z[k][0], Z[k][1]);
    }
    free(Z);
    printf("Vector 3: z (%lld complex roots):\n", order);
    print_vector_c(roots, order);
    return status;
}

int lsolvez(const myZArray_t * Z, const double * b, _Dcomplex * x){
#ifdef USE_LAPACKE
//#include "mkl_lapacke.h"
    lapack_int info,n,lda,ldb,nrhs;
    nrhs = 1;
    n = M_ORDER - 1;
    lda = n; ldb = 1;
    lapack_int ipiv[M_ORDER];

    MKL_Complex16 * a = calloc(M_ORDER*M_ORDER, sizeof(MKL_Complex16));
    const size_t arrlen = M_ORDER*M_ORDER * sizeof(_Dcomplex);
    memcpy_s(a, arrlen, Z, arrlen);

    MKL_Complex16 * b_mod = calloc(M_ORDER, sizeof(MKL_Complex16));
    for (size_t k = 0; k < M_ORDER; k++){
        b_mod[k].real = b[k];
        b_mod[k].imag = 0.0f;
    }

    info = LAPACKE_zgesv(LAPACK_ROW_MAJOR, n, nrhs,
                         a, lda, ipiv,
                         b_mod, ldb);
    const size_t veclen = M_ORDER * sizeof(_Dcomplex);
    memcpy_s(x, veclen, b_mod, veclen);
    free(b_mod);
    return info;
#else
    printf("Lapacke unavailable, skipping lsolvez!\n");
    return -1;
#endif
}
