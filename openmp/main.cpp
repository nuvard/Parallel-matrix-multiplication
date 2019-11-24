#include <omp.h>
#include <iostream>
#include "cannon.h"
#include "fox.h"
#include "utils.h"
#include "Matrix.h"
#include "band.h"
#define CHUNK 16
using namespace std;

int main(){
    int nThreads = omp_get_max_threads();
    int Size = 10;
    Matrix * A = new Matrix(Size, Size);
    Matrix * B = new Matrix(Size, Size);
    Matrix * C = new Matrix(Size, Size);


    cout << "==> Initialization" << endl;
    initialize(A, B, C);
    A->print("A");
    B->print("B");
    C->print("C");
    
    cout << "==> Multiplication" << endl;
    cout << "=====> Naive multiplication" << endl;
    double begin = omp_get_wtime();
    matrix_product(C, A, B);
    double end = omp_get_wtime();
    double el_time = end-begin;
    log(C, el_time, "C (naive)");
    cout << "=====> Band multiplication" << endl;
    BandMultiplication(A,B,C);
    cout << "=====> Cannon multiplication(openmp)" << endl;
    CannonMultiplex(A,B,C);
    cout << "=====> Fox multiplication (openmp)" << endl;
    FoxMultiplex(A,B,C);
    cout << "==> Finished" << endl;
    cout << "==> Chunk " << CHUNK  << " * Time " << el_time<< endl;
    return 0;
}

