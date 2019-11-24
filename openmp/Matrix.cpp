//
// Created by nastya on 11/23/19.
//

#include "Matrix.h"
Matrix::Matrix(int ncol, int nrow){
    this-> ncol = ncol;
    this -> nrow = nrow;
    this -> data = new double * [nrow];
    for (int i=0; i< ncol; ++i){
        this -> data[i] = new double [ncol];
    }

};


void Matrix::randomInit(){
    for (int i = 0; i < this-> nrow; ++i){
        for (int j = 0; j< this-> ncol; ++j){
            this-> data[i][j]=0.1 *( rand() % 10);
        }
    }
};

void Matrix::toZeros(){
    for (int i = 0; i < this-> nrow; ++i){
        for (int j = 0; j< this-> ncol; ++j){
            this-> data[i][j]=0;
        }
    }
};

void Matrix::print(string iden){
    if(this->ncol>10 || this->nrow > 10){
        cout <<"Matrix " << iden << " is too big!"<< endl;
    }
    else{
        cout << "Matrix:" << iden << endl;
        for(int i = 0; i < this->nrow; i++){
            for(int j = 0; j < this->ncol; j++){
                cout << this -> data[i][j] << " ";
            }
            cout << endl;
        }
    }
};

Matrix::~Matrix(){
    for(int i=0; i<this-> nrow;++i){
        delete data[i];
    }
};
