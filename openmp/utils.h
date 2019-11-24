//
// Created by nastya on 11/23/19.
//

#ifndef OPENMP_MULTIPLICATION_UTILS_H
#define OPENMP_MULTIPLICATION_UTILS_H


#include "Matrix.h"

void Flip(double **&B, double**&Bt, int &Size){
    for(int i = 0; i < Size; ++i)
        for(int j = 0; j < Size; ++j)
            Bt[i][j] = B[j][i];
}



void initialize(Matrix*&A, Matrix*&B, Matrix*&C) {
    A->randomInit();
    B->randomInit();
    C->toZeros();
}

void shift_matrix_left(Matrix *&m, int block_sz, int initial) {
    int i, j, k, s, step = block_sz;

    Matrix * aux = new Matrix( m->ncol,1);
    //printf("aux: %d*%d\n", aux->nrow, aux->ncol);
    aux->toZeros();
    //printf("**Shifting\n");
    for(k = 0, s = 0; k < m->ncol; k += block_sz, s++){
        for(i = k; i < (k + block_sz); i++){
            if(initial > 0){
		//printf("Step is: %d\n", s*block_sz);
                step = s * block_sz;
            }
            for(j = 0; j < m->ncol; j++){
	        //printf("%d: %d", j, (j+step)%m->ncol );
                aux->data[0][j] = m->data[i][(j + step) % m->ncol];
            }
            for(j = 0; j < m->ncol; j++){
		//printf("Before: %f, after: %f", m->data[i][j], aux->data[0][j]);
                m->data[i][j] = aux->data[0][j];
            }
        }
    }

}


/**
 * shift the given matrix up.
 *
 * @param m: the matrix to shift.
 * @param initial: a value > 0 indicates that it is a first shift, otherwise is
 *                 normal shift.
 */
void shift_matrix_up(Matrix *&m, int block_sz, int initial) {
    int i, j, k, s, step = block_sz;
    Matrix *aux = new Matrix(m->nrow,1);
    //printf("**Shifting\n");
    for(k = 0, s = 0; k < m->nrow; k += block_sz, s++){
        for(i = k; i < (k + block_sz); ++i){
            if(initial > 0){
                step = s * block_sz;
            }
            for(j = 0; j < m->nrow; j++){
                aux->data[0][j] = m->data[(j + step) % m->nrow][i];
            }
            for(j = 0; j < m->nrow; j++){
                m->data[j][i] = aux->data[0][j];
            }
        }
    }
}

/**
 * Matrix multiplication
 **/
void matrix_product(Matrix *&c, Matrix *&a, Matrix *&b) {
    int r, s, k;
    //cout << "A rows: " << a->nrow << " A cols: " << a->ncol << " B rows: " << b->ncol << " B cols: " << b->ncol<<endl;
#pragma omp parallel for shared(a,b,c) private(r,s,k)
    for(r = 0; r < a -> nrow; r++){
        for(s = 0; s < b->ncol; s++){
            for(k = 0; k < a->ncol; k++){
#pragma omp critical
                c->data[r][s] += a->data[r][k] * b->data[k][s];
	//	cout << "C"<<"["<< r<< "]["<< s <<"]"<< " = " << c->data[r][s]<< endl;
            }
        }
    }
    
}

void log(Matrix *&C, double t, string iden){
    C->print("C");
    FILE * f = fopen("log.txt", "a");
    fprintf(f, "%s Elapsed time: %f\n", iden, t);
    fclose(f);
    cout << "=====>" << iden << " * Time " << t << endl;
}

#endif //OPENMP_MULTIPLICATION_UTILS_H

