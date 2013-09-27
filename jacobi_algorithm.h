#ifndef JACOBI_ALGORITHM_H
#define JACOBI_ALGORITHM_H

void jacobi_algorithm (int n, double **A, double epsilon=1e-8);
void similarity_transformation(int n, double **A, int *k, int *l);
double max_nondiag(int n, double **A, int *k, int *l);

#endif // JACOBI_ALGORITHM_H
