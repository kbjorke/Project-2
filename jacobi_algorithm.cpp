#include <cmath>
#include "jacobi_algorithm.h"

using namespace std;

void jacobi_algorithm ()
{

}

void similarity_transformation(int n, double **A, int *k, int *l)
{
    int i;
    int j;
}











double max_nondiag(int n, double **A, int *k, int *l)
{
    int i;
    int j;

    double max_element = 0;

    for( i = 0; i < n; i++ ){
        for( j = 0; j < n; j ++ ){
            if( i != j ){
                if( fabs(A[i][j]) > max_element ){
                    *k = i;
                    *l = j;
                    max_element = A[i][j];
                }
            }
        }
    }
    return max_element;
}
