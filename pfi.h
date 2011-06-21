#ifndef PFI_H
#define PFI_H

#define TRUE 1
#define FALSE 0

#ifdef _MPI_
inline static void gather_matrix(double *temp_rbuf, double *temp_sbuf, double **matrix, int i_length, int j_length);
#endif

inline static void save_omega_psi(int Nx, int My, double **omega, double **omega_0, double **psi, double **psi_0);

inline static void omega_calculation(int Nx, int My, double Cx2, double Cy2, double alpha, double alphaX, double alphaY, double **omega, double **omega_0, double **u, double **v, double *x, double *y, double **d_sq);

inline static void psi_calculation(int Nx, int My, double kappa2, double kappaA, double dxx, double **psi, double **psi_calc_0, double **omega, double k_psi, double *x, double *y, double tolerance_level, double **d_sq);

inline static void boundary_conditions(int Nx, int My, double **omega, double **psi, double **u, double **v, double *x, double *y, double dyy, double A, int ibl);

inline static void velocity_calculation(int Nx, int My, double **psi, double **u, double **v, double *x, double *y, double dx2, double dy2);

inline static void find_max_values(int Nx, int My, double **omega, double **omega_0, double **psi, double **psi_0, double *OmTol, double *PsiTol);

inline static void initialize_u_v_psi_omega(int i, int j, double **u, double **v, double **psi, double **omega, double *x, double *y, double A);

inline static void write_checkpoint(char *filename, int iteration, double **omega, double **psi, double **u, double **v, int length1, int length2);
inline static short read_checkpoint(char *filename, int *iteration, double **omega, double **psi, double **u, double **v, int length1, int length2);


void print_flow_field_initialization(FILE* out, double *x, double *y, double **u, double **v, double **omega, double **psi, double Nx, double My);

void print_1d(const char *name, double *x, int length, unsigned short only_total, FILE *file);

void print_2d(const char *name, double **x, int length1, int length2, unsigned short only_total, FILE *file);


#endif
