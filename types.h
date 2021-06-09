#ifndef TYPES_H
#define TYPES_H

// Prony model constants
#define M_ORDER 4ll
#define MLB M_ORDER
#define MLA 5ll
//const size_t m_order = 4;
//const size_t mla = m_order, mlb = m_order;

#include <complex.h>
#if ! defined _MSC_VER
#define _Dcomplex double _Complex
static _Dcomplex _Cbuild(const double A, const double B){
    return A + B*I;
}
#endif

// complex double imaginary unit
//#define I _Cbuild(0.0F,1.0F)

typedef double myArray_t[MLA][MLB];
typedef _Dcomplex myZArray_t[M_ORDER][M_ORDER];

#endif // TYPES_H
