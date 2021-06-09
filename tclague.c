/****************************************************
*   Find all roots of a complex polynomial using    *
*   Laguerre formulation in complex domain.         *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
* (Find roots of complex polynomial:                *
*  (2.5 + I) X3 + (7.5 - 12 I) X2 + (-3.75 + 0.4 I  *
*  ) X + (4.25 - I)  )                              *
*                                                   *
*  Error code = 0                                   *
*                                                   *
*  Complex roots are:                               *
*                                                   *
*  (-1.079156,4.904068)                             *
*  (-0.111846,0.668007)                             *
*  (0.259967,-0.399661)                             *
*                                                   *
* ------------------------------------------------- *
* Ref.: From Numath Library By Tuan Dang Trong in   *
*       Fortran 77 [BIBLI 18].                      *
*                                                   *
*                C++ Release By J-P Moreau, Paris.  *
*                       (www.jpmoreau.fr)           *
****************************************************/
#include <stdio.h>
#include <math.h>

#define  NMAX  25      

typedef double Complex[2];

Complex A[NMAX], RAC[NMAX], W1[NMAX], W2[NMAX], W3[NMAX]; //tables of complex numbers
int I, IMP, ITMAX, N;
double EPS,EPS2;
Complex CZERO;


//*** Utility functions for complex numbers ***

// absolute value of z
double CABS(Complex z) {
    double XX,YY,X,Y,W;
    XX=z[0]; YY=z[1];
    X=fabs(XX);
    Y=fabs(YY);
    if (X == 0.0)
        W=Y;
    else {
        if (Y == 0.0)
            W=X;
        else {
            if (X > Y)
                W=X*sqrt(1.0+(Y/X)*(Y/X));
            else
                W=Y*sqrt(1.0+(X/Y)*(X/Y));
        }
    }
    return W;
}

// z3=z1+z2
void CADD(Complex z1, Complex z2, Complex z3) {
    z3[0]=z1[0]+z2[0];
    z3[1]=z1[1]+z2[1];
}

// let z1=z
void CAssign(Complex z, Complex z1) {
    z1[0]=z[0]; z1[1]=z[1];
}

// Z=Z1/Z2
void CDIV(Complex Z1, Complex Z2, Complex Z) {
    double D;
    D=Z2[0]*Z2[0]+Z2[1]*Z2[1];
    if (D<1e-12) return;
    Z[0]=(Z1[0]*Z2[0]+Z1[1]*Z2[1])/D;
    Z[1]=(Z1[1]*Z2[0]-Z1[0]*Z2[1])/D;
}

// Z=Z1*Z2
void CMUL(Complex Z1, Complex Z2, Complex Z) {
    Z[0]=Z1[0]*Z2[0] - Z1[1]*Z2[1];
    Z[1]=Z1[0]*Z2[1] + Z1[1]*Z2[0];
}

// z=x+iy
void CMPLX(Complex z, double x, double y) {
    z[0]=x; z[1]=y;
}

//print a complex number
void CPrint(Complex z) {
    printf(" (%f,%f)\n", z[0], z[1]);
}

void CSQRT(Complex z, Complex z1) {
    //  SQUARE ROOT OF A COMPLEX NUMBER  A+I*B = SQRT(X+I*Y)
    double X,Y,A,B;
    X=z[0]; Y=z[1];
    if (X == 0.0 && Y == 0.0) {
        A=0.0;
        B=0.0;
    }
    else {
        A=sqrt(fabs(X)+CABS(z)*0.5);
        if (X >= 0.0)
            B=Y/(A+A);
        else
            if (Y < 0.0)
                B=-A;
            else {
                B=A;
                A=Y/(B+B);
            }
    }
    CMPLX(z1,A,B);
}

// z3=z1-z2
void CSUB(Complex z1, Complex z2, Complex z3) {
    z3[0]=z1[0]-z2[0];
    z3[1]=z1[1]-z2[1];
}


void CLAGUE(int N, Complex *A, int ITMAX, double EPS, double EPS2,
            int *IMP, Complex *X, Complex *B, Complex *C, Complex *D)  {
/*================================================================
!     ROOTS OF A COMPLEX COEFFICIENTS POLYNOMIAL
!     BY LAGUERRE FORMULA IN COMPLEX DOMAIN
!=================================================================
!     CALLING MODE:
!       CLAGUE(N,A,ITMAX,EPS,EPS2,IMP,X,B,C,D);
!     INPUTS:
!     N : ORDER OF POLYNOMIAL
!     A : TABLE OF SIZE (0:N) OF COMPLEX COEFFICIENTS STORED IN
!         DECREASING ORDER OF X POWERS
!     ITMAX: MAXIMUN NUMBER OF ITERATIONS FOR EACH ROOT
!     EPS:   MINIMAL RELATIVE ERROR
!     EPS2:  MAXIMAL RELATIVE ERROR
!     X(0):  APPROXIMATE VALUE OF FIRST ROOT (=0. GENERALLY)
!     OUTPUTS:
!     IMP:  FLAG = 0  CONVERGENCE WITH AT LEAST EPS2 PRECISION
!                = 1  NO CONVERGENCE
!     X:  TABLE OF SIZE (0:N) OF FOUND ROOTS
!         X(I),I=1,N
!     WORKING ZONE:
!     W:  TABLE OF SIZE (0:3*N), here divided into B, C, D.
! ----------------------------------------------------------------
!     REFERENCE:
!     E.DURAND. SOLUTIONS NUMERIQUES DES EQUATIONS ALGEBRIQUES
!               TOME I, MASSON & CIE PAGES 269-270
!================================================================*/
    //Labels: e1,e2,e3,e4
    Complex XK,XR,F,FP,FS,H,DEN,SQ,D2,TMP1,TMP2,TMP3;
    double EPS1,TEST;
    int I,IK,IT;

    *IMP=0;
    CAssign(A[0],B[0]);
    CAssign(B[0],C[0]);
    CAssign(C[0],D[0]);

    if (N == 1) {
        //X[1]=-A[1]/A[0]
        CDIV(A[1],A[0],X[1]);
        X[1][0]=-X[1][0];  //change sign of X[1]
        X[1][1]=-X[1][1];
        return;
    }

    CAssign(X[0],XK);
e1:   IK=0;
    EPS1=EPS;
    IT=0;

e2:   for (I=1; I<=N; I++) {
        //B[I]=A[I]+XK*B[I-1];
        CMUL(XK,B[I-1],TMP1);
        CADD(A[I],TMP1,B[I]);
        if (I <= N-1) {
            //C[I]=B[I)+XK*C[I-1)
            CMUL(XK,C[I-1],TMP1);
            CADD(B[I],TMP1,C[I]);
        }
        if (I <= N-2) {
            //D[I]=C[I)+XK*D[I-1)
            CMUL(XK,D[I-1],TMP1);
            CADD(C[I],TMP1,D[I]);
        }
    }

    CAssign(B[N],F);
    CAssign(C[N-1],FP);
    CAssign(D[N-2],FS);
    //H=((N-1)*FP)^2-N*(N-1)*F*FS
    CMPLX(TMP1,1.0*(N-1),0.0); CMUL(TMP1,FP,TMP2);
    CMUL(TMP2,TMP2,TMP1); //TMP1=((N-1)*FP)^2
    CMPLX(TMP2,1.0*N*(N-1),0.0); CMUL(TMP2,F,TMP3);
    CMUL(TMP3,FS,TMP2);   //TMP2=N*(N-1)*F*FS
    CSUB(TMP1,TMP2,H);
    if (CABS(H) == 0.0) {
        for (I=N; I>0; I--) {
            //X[I]=-A[1]/(N*A[0])
            CMPLX(TMP1,1.0*N,0.0); CMUL(TMP1,A[0],TMP2);
            CDIV(A[1],TMP2,X[I]);
            X[I][0]=-X[I][0];
            X[I][1]=-X[I][1];
        }
        return;
    }
    if (CABS(FP) == 0.0)
        CSQRT(H,DEN);
    else {
        CSQRT(H,SQ);
        CADD(FP,SQ,DEN);
        CSUB(FP,SQ,D2);
        if (CABS(DEN) < CABS(D2)) CAssign(D2,DEN);
    }
    if (CABS(DEN) == 0.0) {
        XK[0]=XK[0] + 0.1;
        XK[1]=XK[1] + 0.1;
        goto e2;
    }
    IK=IK+1;
    //XR=XK-N*F/DEN
    CMPLX(TMP1,1.0*N,0.0); CMUL(TMP1,F,TMP2);
    CDIV(TMP2,DEN,TMP1);
    CSUB(XK,TMP1,XR);
    CSUB(XR,XK,TMP1);
    TEST=CABS(TMP1);
    CAssign(XR,XK);
    if (TEST < EPS1*CABS(XR)) goto e3;
    //    WRITELN(IK,H[1],H[2],DEN[1],DEN[2],XK[1],XK[2],TEST);
    if (IK <= ITMAX) goto e2;
    EPS1=EPS1*10.0;
    IK=0;
    if (EPS1 < EPS2) goto e2;
    *IMP=1;
    return;

e3:   CAssign(XR,X[N]);
    N=N-1;
    if (N == 1) goto e4;
    if (N <= 0) return;
    for (I=0; I<=N; I++)  CAssign(B[I],A[I]);
    goto e1;

e4:   //X[N]=-B[1]/B[0]
    CDIV(B[1],B[0],X[N]);
    X[N][0]=-X[N][0];
    X[N][1]=-X[N][1];

}

#ifdef TCLAGUE_MAIN
int main() {

    CMPLX(CZERO,0.0,0.0);

    // Example #1
    N=3;
    CMPLX(A[0],2.5,1.0);
    CMPLX(A[1],7.5,-12.0);
    CMPLX(A[2],-3.75,0.4);
    CMPLX(A[3],4.25,-1.0);

    /* Example #2 (real roots are -3, 1, 2)
  CMPLX(A[0],1.0,0.0);
  CMPLX(A[1],0d0,0d0);
  CMPLX(A[2],-7.0,0.0);
  CMPLX(A[3],6.0,0.0);
*/

    ITMAX=10;                 //Maximum number of iterations
    EPS=1e-8; EPS2=1e-6;      //Minimum, maximum relative error
    CAssign(CZERO,RAC[0]);    //Approximate value of 1st root

    // call complex Laguerre procedure
    CLAGUE(N,A,ITMAX,EPS,EPS2,&IMP,RAC,W1,W2,W3);

    // print results
    printf("\n Error code = %d\n\n", IMP);
    printf(" Complex roots are:\n\n");
    for (I=1; I<=N; I++) CPrint(RAC[I]);
    printf("\n");

}
#endif
// end of file tclague.cpp
