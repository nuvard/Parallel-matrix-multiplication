
#ifndef OPENMP_MULTIPLICATION_CANNON_H
#define OPENMP_MULTIPLICATION_CANNON_H



#define P_SQRT 2
#define P (P_SQRT * P_SQRT) //number of processes
#define N 10
#define BLOCK_SZ (N / P_SQRT) //block size
#include "Matrix.h"
#include "utils.h"
#include "omp.h"
/**
 * multiply the corresponding submatrices of A and B.
 */

void process_mult(Matrix *A, Matrix *B, Matrix *C) {
    int r, c, id, k,
            rbegin, rend, cbegin, cend, // block delimiters
            l, m;
    //Matrix *sa = new Matrix(BLOCK_SZ, BLOCK_SZ);
    //Matrix *sb = new Matrix(BLOCK_SZ, BLOCK_SZ);
    //Matrix *sc = new Matrix(BLOCK_SZ, BLOCK_SZ);


#pragma omp parallel default(none) private(l, m, r, c, k, rbegin, rend, cbegin, cend, id) shared(A, B, C) num_threads(P)
    {
        id = omp_get_thread_num();
        rbegin = (id / P_SQRT) * BLOCK_SZ;
        rend = rbegin + BLOCK_SZ;

        cbegin = (id % P_SQRT) * BLOCK_SZ;
        cend = cbegin + BLOCK_SZ;
        Matrix *sa = new Matrix(BLOCK_SZ, BLOCK_SZ);
	Matrix *sb = new Matrix(BLOCK_SZ, BLOCK_SZ);
	Matrix *sc = new Matrix(BLOCK_SZ, BLOCK_SZ);
        sa->toZeros();
        sb->toZeros();
        sc->toZeros();
        //copy the blocks for this process
        for(r = rbegin, l = 0; r < rend; r++, l++){
            for(c = cbegin, m = 0; c < cend; c++, m++){
                sa->data[l][m] = A->data[r][c];
                sb->data[l][m] = B->data[r][c];
                sc->data[l][m] = C->data[r][c];
            }
        }

        matrix_product(sc, sa, sb);

        //put results back to C
        for(r = rbegin, l = 0; r < rend; r++, l++){
            for(c = cbegin, m = 0; c < cend; c++, m++){
                C->data[r][c] = sc->data[l][m];
            }
        }
    }
}


void CannonMultiplex(Matrix *A, Matrix *B, Matrix *C){
    shift_matrix_left(A, BLOCK_SZ, 1);
    shift_matrix_up(B, BLOCK_SZ, 1);
    double t1, t2;
    t1 = omp_get_wtime();
    for(int i = 0; i < P_SQRT; i++){
        process_mult(A, B, C);
        shift_matrix_left(A, BLOCK_SZ, 0);
        shift_matrix_up(B, BLOCK_SZ, 0);
    }
    t2 = omp_get_wtime();
    log(C, t2-t1, "Cannon");
}

#endif //OPENMP_MULTIPLICATION_CANNON_H
