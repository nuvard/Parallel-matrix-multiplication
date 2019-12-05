#include<stdlib.h>
#include <ctime>
#include<math.h>
#include"mpi.h"
#include <iostream>
#include <cstdio>
#define MATRIX_SIZE 100
using namespace std;
int first_matrix[MATRIX_SIZE][MATRIX_SIZE];
int second_matrix[MATRIX_SIZE][MATRIX_SIZE];

typedef struct {
    int proc_count;
    int dim;
    int row;
    int col;
    int rank;
    MPI_Comm grid_comm;
    MPI_Comm row_comm;
    MPI_Comm col_comm;
} GridStructure;

void MultiplyLocal(int **a, int **b, int **c, int size) {
    int temp = 0;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++) {
            temp = 0;
            for (int k = 0; k < size; k++)
                temp += (a[i][k] * b[k][j]);
            c[i][j] += temp;
        }
}
void UnpackMatrix(int *buff, int **a, int size) {
    int k = 0;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++) {
            a[i][j] = buff[k];
            k++;
        }
}
void PackMatrix(int *buff, int **a, int size) {
    int k = 0;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++) {
            buff[k] = a[i][j];
            k++;
        }
}
void PrintMatrix(int **matrix, int size) {
    for (int i = 0; i < size; i++) {
        std::cout << "|";
        for (int j = 0; j < size; j++) {
            int el = matrix[i][j];
            if (el < 10)
                cout << " ";
            cout << el;
            cout << "|";
        }
        cout <<endl;
    }
}
void PrintPackedMatrix(int *matrix, int size) {
    for (int i = 0; i < size; i++) {
        cout << "|";
        for (int j = 0; j < size; j++) {
            int el = matrix[i * size + j];
            if (el < 10)
               cout << " ";
            cout << el;
            cout << "|";
        }
        cout << endl;
    }
}
void GenerateMatrices(int **a, int **b) {
	printf("Generating matrix!\n");
     printf("proc num: ");
     for (int i = 0; i < MATRIX_SIZE; i++)
        for (int j = 0; j < MATRIX_SIZE; j++) {
            a[i][j] =0.1*(rand() % 100)+1;//(i==j) ? 1 : 0;
            b[i][j] =0.1*( rand() % 100)+1;
        }
    printf("Matrix generated!\n");
}


void FoxMultiply(int n, GridStructure *grid, int **a, int **b, int **c) {
    int **temp_a, *buff, stage, root, submat_dim, src, dst;
    MPI_Status status;

    submat_dim = n / grid->dim;
    //printf("submat dim: %i\n", submat_dim);
    temp_a = new int*[submat_dim];
    for(int i = 0; i < submat_dim; ++i)
        temp_a[i] = new int[submat_dim];
    for (int i = 0; i < submat_dim; i++)
        for (int j = 0; j < submat_dim; j++)
            temp_a[i][j] = 0;
    
    buff = new int[submat_dim*submat_dim];
    for (int i = 0; i < submat_dim * submat_dim; i++)
        buff[i] = 0;

    src = (grid->row + 1) % grid->dim;
    dst = (grid->row + grid->dim - 1) % grid->dim;

    for (stage = 0; stage < grid->dim; stage++) {
        root = (grid->row + stage) % grid->dim;
        if (root == grid->col) {
            PackMatrix(buff, a, submat_dim);
            MPI_Bcast(buff, submat_dim * submat_dim, MPI_INT, root, grid->row_comm);
            UnpackMatrix(buff, a, submat_dim);
            MultiplyLocal(a, b, c, submat_dim);
        } else {
            PackMatrix(buff, temp_a, submat_dim);
            MPI_Bcast(buff, submat_dim * submat_dim, MPI_INT, root, grid->row_comm);
            UnpackMatrix(buff, temp_a, submat_dim);
            MultiplyLocal(temp_a, b, c, submat_dim);
        }
        PackMatrix(buff, b, submat_dim);
        MPI_Sendrecv_replace(buff, submat_dim * submat_dim, MPI_INT, dst, 0, src, 0, grid->col_comm, &status);
        UnpackMatrix(buff, b, submat_dim);
    }
}

int** TestSingleThread()
{
    int **test_result, **test_a, **test_b;
    test_a = new int*[MATRIX_SIZE];
    test_b = new int*[MATRIX_SIZE];
    test_result = new int*[MATRIX_SIZE];
    for(int i = 0; i < MATRIX_SIZE; ++i) {
        test_a[i] = new int[MATRIX_SIZE];
        test_b[i] = new int[MATRIX_SIZE];
        test_result[i] = new int[MATRIX_SIZE];
    }
    for (int i = 0; i < MATRIX_SIZE; i++)
        for (int j = 0; j < MATRIX_SIZE; j++) {
            test_a[i][j] = first_matrix[i][j];
            test_b[i][j] = second_matrix[i][j];
            test_result[i][j] = 0;
        }
    MultiplyLocal(test_a, test_b, test_result, MATRIX_SIZE);

    //std::cout << "Local result:" << std::endl;
    //PrintMatrix(test_result, MATRIX_SIZE);
    return test_result;
}


int main(int argc, char *argv[]) {
    int block_size;
    int **a = new int*[MATRIX_SIZE];
    int **b = new int*[MATRIX_SIZE];
    for (int i=0;i<MATRIX_SIZE; ++i)
    {
	    a[i]=new int [MATRIX_SIZE];
	    b[i]=new int [MATRIX_SIZE];
    }

    int **local_a, **local_b, **local_c;
    clock_t start,finish;
    MPI_Init(&argc, &argv);
    //printf("==> Setup grid\n");
    GridStructure grid;
    int dimensions[2];
    int wrap_around[2];
    int coordinates[2];
    int free_coords[2];
    int world_rank;
    
    MPI_Comm_size(MPI_COMM_WORLD, &(grid.proc_count));
    //printf("Count of processors: %i\n", grid->proc_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    grid.dim = (int) sqrt((double) grid.proc_count);
    //printf("Dim: %i\n", grid->dim);
    dimensions[0] = dimensions[1] = grid.dim;
    //is grid periodic or not?
    wrap_around[0] = 0;
    wrap_around[1] = 1;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, wrap_around, 1, &(grid.grid_comm));
    if (grid.grid_comm==MPI_COMM_NULL){
	    MPI_Finalize();
	    return 0;
    }
    MPI_Comm_rank(grid.grid_comm, &(grid.rank));
    MPI_Cart_coords(grid.grid_comm, grid.rank, 2, coordinates);
    grid.row = coordinates[0];
    grid.col = coordinates[1];
    //printf("New coords of proc %i = %i %i\n", grid->rank, coordinates[0], coordinates[1]);
    free_coords[0] = 0;
    free_coords[1] = 1;
               //split grid to comms rowwise (communicator for rows)
    MPI_Cart_sub(grid.grid_comm, free_coords, &(grid.row_comm));
    free_coords[0] = 1;
    free_coords[1] = 0;
    MPI_Cart_sub(grid.grid_comm, free_coords, &(grid.col_comm));
       printf("==> Generating matrices\n");
    int * buff = new int[MATRIX_SIZE*MATRIX_SIZE];
    if(grid.rank==0)
    {
         GenerateMatrices(a,b);
    }
    PackMatrix(buff, a, MATRIX_SIZE);
    MPI_Bcast(buff, MATRIX_SIZE*MATRIX_SIZE, MPI_INT, 0, grid.grid_comm);
    UnpackMatrix(buff, a, MATRIX_SIZE);
    
    PackMatrix(buff, b, MATRIX_SIZE);
    MPI_Bcast(buff, MATRIX_SIZE*MATRIX_SIZE, MPI_INT,0, grid.grid_comm);
    UnpackMatrix(buff, b, MATRIX_SIZE);

    block_size = MATRIX_SIZE / grid.dim;
    int base_row = grid.row * block_size;
    int base_col = grid.col * block_size;

    local_a = new int*[MATRIX_SIZE];
    local_b = new int*[MATRIX_SIZE];
    local_c = new int*[MATRIX_SIZE];
    for(int i = 0; i < MATRIX_SIZE; ++i)
    {
        local_a[i] = new int[MATRIX_SIZE];
        local_b[i] = new int[MATRIX_SIZE];
        local_c[i] = new int[MATRIX_SIZE];
    }
    for (int i = base_row; i < base_row + block_size; i++)
        for (int j = base_col; j < base_col + block_size; j++) {
            local_a[i - (base_row)][j - (base_col)] = first_matrix[i][j];
            local_b[i - (base_row)][j - (base_col)] = second_matrix[i][j];
            local_c[i - (base_row)][j - (base_col)] = 0;
        }

    if (grid.rank == 0)
    {
        cout<<"Ready..." <<endl;
    }

    MPI_Barrier(grid.grid_comm);
    if (grid.rank == 0)
    {
        start = clock();
    }
    FoxMultiply(MATRIX_SIZE, &grid, local_a, local_b, local_c);
    MPI_Barrier(grid.grid_comm);
    if (grid.rank == 0)
    {
        finish = clock();
        clock_t result_time = finish - start;
        cout << "Time: " << double(result_time)/CLOCKS_PER_SEC <<endl;
    }
    printf("Multiplied! Collecting result...\n");

    int *result_buff = new int[MATRIX_SIZE * MATRIX_SIZE];
    int *local_buff = new int[block_size * block_size];
    PackMatrix(local_buff, local_c, block_size);

    MPI_Gather(local_buff, block_size * block_size, MPI_INT, result_buff, block_size * block_size, MPI_INT, 0, grid.grid_comm);
   // MPI_Barrier(grid.grid_comm);
    if (grid.rank == 0) {
        int *data = new int[MATRIX_SIZE * MATRIX_SIZE];
        int k = 0;
        for (int bi = 0; bi < grid.dim; bi++)
            for (int bj = 0; bj < grid.dim; bj++)
                for (int i = bi * block_size; i < bi * block_size + block_size; i++)
                    for (int j = bj * block_size; j < bj * block_size + block_size; j++) {
                        data[i * MATRIX_SIZE + j] = result_buff[k];
                        k++;
                    }

        cout << "Fox result check:";

        //PrintPackedMatrix(data, MATRIX_SIZE);
        int ** res = TestSingleThread();
        for(int i=0;i<MATRIX_SIZE;i++)
            for(int j=0;j<MATRIX_SIZE;j++)
            {
                if(res[i][j]!=data[i*MATRIX_SIZE+j])
                {
                    cout << "ERROR!" << endl;
                    exit(1);
                }
            }
       cout << "OK" << endl;
    }

    MPI_Finalize();
    exit(0);
}		
