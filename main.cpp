#include <unittest++/UnitTest++.h>
#include <iostream>
#include "jacobi_algorithm.h"

using namespace std;

int main()
{
    return UnitTest::RunAllTests();
}



TEST(max_nondiag_func){

    int i;
    int size = 2;
    double max_element;

    double **matrix = new double*[size];
    for( i = 0; i < size; i++ ){
        matrix[i] = new double[size];
    }

    matrix[0][0] = 1;
    matrix[0][1] = 2;
    matrix[1][0] = 4;
    matrix[1][1] = 3;

    int k;
    int l;

    max_element = max_nondiag(size, matrix, &k, &l);

    CHECK( ( max_element == 4 ) && k == 1 && l == 0);

    for( i = 0; i < size; i++ ){
         delete[] matrix[i];
    }
    delete[] matrix;
}
