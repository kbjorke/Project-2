#ifndef OUTPUT_FUNCTIONS_H
#define OUTPUT_FUNCTIONS_H

void eigenstates_index(int m, double *index_list, int n, double *eigenvalues);
void output_eigenvalues(int n, double omega_r, double rho_max,
                        int m, const char *method = "jacobi");
void output_stability(int n_start, int n_stop, double resolution, double rho_max,
                      int m = 3, const char *method = "jacobi");
void output_eigenstates(int m, int n_steps, double rho_max,
                        double omega_r = 0, const char *method = "jacobi");
#endif // OUTPUT_FUNCTIONS_H
