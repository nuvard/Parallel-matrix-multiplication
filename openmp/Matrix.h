//
// Created by nastya on 11/23/19.
//

#ifndef OPENMP_MULTIPLICATION_MATRIX_H
#define OPENMP_MULTIPLICATION_MATRIX_H
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>

using namespace std;
class Matrix {
public:
    int ncol;
    int nrow;
    double ** data;
    Matrix(int ncol, int nrow);
    void randomInit();
    void toZeros();
    void print(string iden);
    ~Matrix();

};


#endif //OPENMP_MULTIPLICATION_MATRIX_H

