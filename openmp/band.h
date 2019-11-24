#ifndef OPENMP_MULTIPLICATION_BAND_H
#define OPENMP_MULTIPLICATION_BAND_H

void BandMultiplication(Matrix *&A, Matrix *&B, Matrix*&C)
{
    int i, j,k;
    double sum;
    double t1 = omp_get_wtime();
    #pragma omp parallel for private(j,k,sum)
    for(i=0;i<A->nrow;i++)
    {
        for(k=0;k<B->ncol;k++)
	{
            sum=0;
            for(j=0;j<A->ncol;j++)
	    {
                sum = sum + (A->data[i][j]*B->data[j][k]);
                
	    }

            #pragma omp critical
            C->data[i][k]=sum;
        }
    }
    double el_time = omp_get_wtime() - t1;
    log(C, el_time, " C (band)");
}



#endif //OPENMP_MULTIPLICATION_BAND_H

