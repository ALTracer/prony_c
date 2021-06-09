#ifndef TCLAGUE_H
#define TCLAGUE_H

typedef double Complex[2];

/**
 * @brief CLAGUE Complex Laguerre polynomial roots solver
 * @param N Order of polynomial
 * @param A Table of size (0:N) of complex coef. in decr. order of X powers
 * @param ITMAX Maximum number of iterations for each root
 * @param EPS Minimal relative error
 * @param EPS2 Maximal relative error
 * @param IMP Flag of no convergence better than eps2
 * @param X X(0) = 0 First root (approx.); Table (0:N) of found roots
 * @param B Working table, part 1
 * @param C Working table, part 2
 * @param D Working table, part 3
 */
void CLAGUE(int N, Complex *A, int ITMAX, double EPS, double EPS2,
            int *IMP, Complex *X, Complex *B, Complex *C, Complex *D);
#endif // TCLAGUE_H
