#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pfi.h"

#ifdef _WIN32
    #include "boinc_win.h"
    #include "str_util.h"
#endif

#ifdef _BOINC_
    #include "diagnostics.h"
    #include "util.h"
    #include "filesys.h"
    #include "boinc_api.h"
    #include "mfile.h"
#endif

#ifdef _MPI_
    #include "mpi.h"

    int process_id;
    int number_processes;

    int hostname_length;
    char hostname[MPI_MAX_PROCESSOR_NAME];

    int left_column, right_column;
    int *gather_recv_counts, *gather_displacements;
#endif

int i_min, i_max;
int i_min_inner, i_max_inner;

void linspace(double min, double max, int length, double **x, double *dx) {
    int i;

    (*dx) = fabs(max - min) / (length - 1);
    (*x)[0] = min;
    (*x)[length - 1] = max;

    for (i = 1; i < (length - 1); i++) {
        (*x)[i] = min + ((*dx) * i);
    }
}

int main(int number_arguments, char** arguments) {
    int i, j, k;
    int k_psi;
    int retval;

//    double t;

    double **v;                             // the y velocity component matrix
    double **u;                             // the x velocity component matrix
    double **psi;                           //the streamline matrix
    double **omega;                         //the vorticity matrix
    double **psi_0;
    double **psi_calc_0;
    double **omega_0;
    double **d_sq;

    double reynolds_number;                 // Re
    double cell_reynolds_number;            // Rc
    double A;
    double Ot;                              //number of points
    double dt;                              //spacing in time
    double tolerance_level;

    int initial_iteration;
    int report_frequency;
    int ibl;                                //number of boundary level BCs
    int incremental_file_frequency;
    int output_flow_field_initialization;

    char *name;
    char* output_file;
    char checkpoint_filename[512], report_filename[512];

    double dx, dx2, dxx;
    double dy, dy2, dyy;
    double KappaA, Kappa2;
    double C, Cx, Cx2, Cy, Cy2;
    double alpha, alphaX, alphaY;
    double d6, sqrt_xy_squares;

    double *x;
    double *y;

    double Xmin = -20;
    double Xmax = 20;
    double Ymin = 1;
    double Ymax = 11;
    double Nx = 199;
    double My = 399;
    double Umax = 1;

    double OmTol, PsiTol;

    double *temp_rbuf, *temp_sbuf;

#ifdef _MPI_
    double slice_size;
    int displacement;

    MPI_Init(&number_arguments, &arguments);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
    MPI_Comm_size(MPI_COMM_WORLD, &number_processes);
    MPI_Get_processor_name(hostname, &hostname_length);

    if (number_processes > (Nx + 2)) {
        fprintf(stderr, "MPI ERROR: number_processes [%d] > matrix columns (Nx + 2) [%d]\n", number_processes, (int)(Nx + 2));
        exit(0);
    }

    slice_size = ((double)Nx + 2.0) / number_processes;

    left_column = (int)floor(slice_size * (double)process_id);
    if (process_id > 0) left_column++;

    right_column = (int)floor(slice_size * (double)(process_id + 1));

    i_min = left_column;

    i_max = right_column;
    if (process_id != (number_processes - 1)) i_max++;

    i_min_inner = left_column;
    i_max_inner = right_column;
    if (process_id == (number_processes - 1)) {
        i_max_inner--;
    } else if (process_id == 0) {
        i_min_inner++;
        i_max_inner++;
    } else {
        i_max_inner++;
    }

    displacement = 0;
    gather_recv_counts = (int*)malloc(sizeof(int) * number_processes);
    gather_displacements = (int*)malloc(sizeof(int) * number_processes);
    for (i = 0; i < number_processes; i++) {
        gather_displacements[i] = displacement;

        gather_recv_counts[i] = ((int)floor(slice_size * (double)(i + 1))) - ((int)floor(slice_size * (double)i)) + 1;
        if (i > 0) gather_recv_counts[i]--;
        if (i == (number_processes - 1)) gather_recv_counts[i]--;
        gather_recv_counts[i] *= (My + 2);

        displacement += gather_recv_counts[i];
    }
    
    temp_rbuf = (double*)malloc(sizeof(double) * (Nx + 2) * (My + 2));
    temp_sbuf = (double*)malloc(sizeof(double) * (Nx + 2) * (My + 2));

//    fprintf(stderr, "[rank: %d] slice size: %lf, left_column: %d, right_column: %d, i_min: %d, i_min_inner: %d, i_max: %d, i_max_inner: %d\n", process_id, slice_size, left_column, right_column, i_min, i_min_inner, i_max, i_max_inner);
#else
    i_min = 0;
    i_min_inner = 1;
    i_max = Nx + 2;
    i_max_inner = Nx + 1;
#endif

    name = NULL;

    x = (double*)malloc(sizeof(double) * (Nx + 2));
    y = (double*)malloc(sizeof(double) * (My + 2));

    /**
     *  The following initializes BOINC client interaction
     */
    retval = 0;
    #ifdef _BOINC_
        #ifdef BOINC_APP_GRAPHICS
            #if defined(_WIN32) || defined(__APPLE)
                retval = boinc_init_graphics(worker);
            #else
                retval = boinc_init_graphics(worker, argv[0]);
            #endif
        #else
            retval = boinc_init();
        #endif

        if (retval) exit(retval);
    #endif

    output_flow_field_initialization = 0;

    for (i = 0; i < number_arguments; i++) {
#ifdef _MPI_
        if (process_id == 0) fprintf(stderr, "argument[%d]: %s\n", i, arguments[i]);
#else
        fprintf(stderr, "argument[%d]: %s\n", i, arguments[i]);
#endif

        if (!strcmp(arguments[i], "--name")) {
            name = arguments[++i];
        } else if (!strcmp(arguments[i], "--reynolds_number")) {
            reynolds_number = atof(arguments[++i]);
        } else if (!strcmp(arguments[i], "--circulation_parameter")) {
            A = atof(arguments[++i]);
        } else if (!strcmp(arguments[i], "--time_steps")) {
            Ot = atoi(arguments[++i]);
        } else if (!strcmp(arguments[i], "--report_frequency")) {
            report_frequency = atoi(arguments[++i]);
        } else if (!strcmp(arguments[i], "--dt")) {
            dt = atof(arguments[++i]);
        } else if (!strcmp(arguments[i], "--tolerance_level")) {
            tolerance_level = atof(arguments[++i]);
        } else if (!strcmp(arguments[i], "--number_boundary_level_bcs")) {
            ibl = atoi(arguments[++i]);
        } else if (!strcmp(arguments[i], "--output_file")) {
            output_file = arguments[++i];
        } else if (!strcmp(arguments[i], "--incremental_file_frequency")) {
            incremental_file_frequency = atoi(arguments[++i]);
        } else if (!strcmp(arguments[i], "--print_flow_field_initialization")) {
            output_flow_field_initialization = 1;
        }
    }

    if (name == NULL) name = "pfi";
    sprintf(checkpoint_filename, "%s_checkpoint", name);

#ifdef _MPI_
    if (process_id == 0) {
#endif
        fprintf(stderr, "reynolds number: %lf\n", reynolds_number);
        fprintf(stderr, "cirulation parameter (A): %lf\n", A);
        fprintf(stderr, "time_steps: %lf\n", Ot);
        fprintf(stderr, "dt: %lf\n", dt);
        fprintf(stderr, "tolerance_level: %lf\n", tolerance_level);
        fprintf(stderr, "number_boundary_level_bcs: %d\n", ibl);
        fprintf(stderr, "report_frequency: %d\n", report_frequency);
        fprintf(stderr, "Ot: %lf\n", Ot);
#ifdef _MPI_
    }
#endif

    linspace(Xmin, Xmax, Nx + 2, &x, &dx);
    dx2 = 2 * dx;
    dxx = dx * dx;

    linspace(Ymin, Ymax, My + 2, &y, &dy);
    dy2 = 2 * dy;
    dyy = dy * dy;

    Kappa2 = (dx / dy);
    Kappa2 = Kappa2 * Kappa2;
    KappaA = 1.0 / (2.0 * (1.0 + Kappa2));

    cell_reynolds_number = reynolds_number * dx;

    Cx = dt / dx;
    Cx2 = 0.5 * Cx;
    Cy = dt / dy;
    Cy2 = 0.5 * Cy;

    if (Cx > Cy) {
        C = Cx;
    } else {
        C = Cy;
    }

    alphaX = dt / (dxx * reynolds_number);
    alphaY = dt / (dyy * reynolds_number);
    alpha = (2 * alphaX) + (2 * alphaY);

    /**
     *  Check for invalid parameters
     */
    if (C >= 1) {
        fprintf(stderr, "ERROR: The Courant Number C (%10.5lf) must be less than 1. [FILE: %s, LINE: %d]\n", C, __FILE__, __LINE__);
        exit(0);
    }

    if (cell_reynolds_number >= (4.0 / C)) {
        fprintf(stderr, "ERROR: The Cell Reynolds Number Rc (%10.5lf) must be less than 4.0 / C (%10.5lf). [FILE: %s, LINE: %d]\n", cell_reynolds_number, (4.0 / C), __FILE__, __LINE__);
        exit(0);
    }

#ifdef _MPI_
    if (process_id == 0) {
#endif
        fprintf(stderr, "Grid Spacing: (dx: %10.5lf, dy: %10.5lf, dt: %10.5lf), [FILE: %s, LINE: %d]\n", dx, dy, dt, __FILE__, __LINE__);
#ifdef _MPI_
    }
#endif

    v = (double**)malloc(sizeof(double*) * (Nx + 2));
    u = (double**)malloc(sizeof(double*) * (Nx + 2));
    psi = (double**)malloc(sizeof(double*) * (Nx + 2));
    psi_0 = (double**)malloc(sizeof(double*) * (Nx + 2));
    psi_calc_0 = (double**)malloc(sizeof(double*) * (Nx + 2));
    omega = (double**)malloc(sizeof(double*) * (Nx + 2));
    omega_0 = (double**)malloc(sizeof(double*) * (Nx + 2));
    d_sq = (double**)malloc(sizeof(double*) * (Nx + 2));

    for (i = 0; i < (Nx + 2); i++) {
        v[i] = (double*)malloc(sizeof(double) * (My + 2));
        u[i] = (double*)malloc(sizeof(double) * (My + 2));
        psi[i] = (double*)malloc(sizeof(double) * (My + 2));
        psi_0[i] = (double*)malloc(sizeof(double) * (My + 2));
        psi_calc_0[i] = (double*)malloc(sizeof(double) * (My + 2));
        omega[i] = (double*)malloc(sizeof(double) * (My + 2));
        omega_0[i] = (double*)malloc(sizeof(double) * (My + 2));
        d_sq[i] = (double*)malloc(sizeof(double) * (My + 2));
    }

    if (read_checkpoint(checkpoint_filename, &initial_iteration, omega, psi, u, v, Nx + 2, My + 2)) {    // returns true if the checkpoint was read
#ifdef _MPI_
        if (process_id == 0)
#endif
            fprintf(stderr, "Starting from checkpoint at %d, initial iteration: %d\n", time(NULL), initial_iteration);
    } else {
        initial_iteration = 0;

        for (i = 0; i < (Nx + 2); i++) {
            d6 = (x[i] * x[i]) + (y[0] * y[0]);
            sqrt_xy_squares = sqrt( d6 );
            v[i][0] = -(y[0] - 1) / sqrt_xy_squares;
            u[i][0] = 0;

            psi[i][0] = (x[i] + A) * (y[0] - 1);

            for (j = 1; j < (My + 2); j++) {
                initialize_u_v_psi_omega(i, j, u, v, psi, omega, x, y, A);
            }

            omega[i][0] = ((7.0 * psi[i][0]) - (8.0 * psi[i][1]) + psi[i][2]) / (2.0 * dyy) / d6;
        }

        if (output_flow_field_initialization) print_flow_field_initialization(stderr, x, y, u, v, omega, psi, Nx, My);

#ifdef _MPI_
        if (process_id == 0) {
#endif
            fprintf(stderr, "Flow-field finished initializing at %d\n", time(NULL));
#ifdef _MPI_
        }
#endif
    }

    for (k = 0; k < Ot; k++) {
        save_omega_psi(Nx, My, omega, omega_0, psi, psi_0);
        omega_calculation(Nx, My, Cx2, Cy2, alpha, alphaX, alphaY, omega, omega_0, u, v, x, y, d_sq);
        psi_calculation(Nx, My, Kappa2, KappaA, dxx, psi, psi_calc_0, omega, k_psi, x, y, tolerance_level, d_sq);
        boundary_conditions(Nx, My, omega, psi, u, v, x, y, dyy, A, ibl);
        velocity_calculation(Nx, My, psi, u, v, x, y, dx2, dy2);
        find_max_values(Nx, My, omega, omega_0, psi, psi_0, &OmTol, &PsiTol);

#ifdef _MPI_
//      print_flow_field_initialization(stderr, x, y, u, v, omega, psi, Nx, My);
        if (process_id == 0 && (k % report_frequency) == 0) fprintf(stderr, "[rank: %d] k: %d, Ot: %lf, OmTol: %.10lf, PsiTol: %.10lf\n\n", process_id, k + initial_iteration, Ot, OmTol, PsiTol);
#else
//      print_flow_field_initialization(stderr, x, y, u, v, omega, psi, Nx, My);
        if ((k % report_frequency) == 0) fprintf(stderr, "k: %d, Ot: %lf, OmTol: %.10lf, PsiTol: %.10lf\n\n", k + initial_iteration, Ot, OmTol, PsiTol);
#endif

        if (k > 0 && incremental_file_frequency > 0 && (k % incremental_file_frequency) == 0) {
#ifdef _MPI_
            /**
             *  Need to collect omega, psi, u, v
             */
            gather_matrix(temp_rbuf, temp_sbuf, omega, Nx + 2, My + 2);
            gather_matrix(temp_rbuf, temp_sbuf, psi, Nx + 2, My + 2);
            gather_matrix(temp_rbuf, temp_sbuf, u, Nx + 2, My + 2);
            gather_matrix(temp_rbuf, temp_sbuf, v, Nx + 2, My + 2);
#endif

#ifdef _MPI_
            if (process_id == 0) {
#endif
                write_checkpoint(checkpoint_filename, initial_iteration + k, omega, psi, u, v, Nx + 2, My + 2);

                sprintf(report_filename, "%s_report_%d", name, initial_iteration + k);
                write_checkpoint(report_filename, initial_iteration + k, omega, psi, u, v, Nx + 2, My + 2);
#ifdef _MPI_
            }
#endif
        }
    }

    /**
     * Free the memory
     */
    for (i = 0; i < (Nx + 2); i++) {
        free(v[i]);
        free(u[i]);
        free(psi[i]);
        free(psi_0[i]);
        free(psi_calc_0[i]);
        free(omega[i]);
        free(omega_0[i]);
        free(d_sq[i]);
    }
    free(v);
    free(u);
    free(psi);
    free(psi_0);
    free(psi_calc_0);
    free(omega);
    free(omega_0);
    free(d_sq);

#ifdef _MPI_
    free(temp_rbuf);
    free(temp_sbuf);
#endif

    free(x);
    free(y);

    #ifdef _BOINC_
        boinc_finish(0);
    #endif

    return 0;
}

#ifdef _MPI_
inline static
void gather_matrix(double *temp_rbuf, double *temp_sbuf, double **matrix, int i_length, int j_length) {
    int i, j;

    for (i = i_min; i < i_max; i++) {
        for (j = 0; j < j_length; j++) {
            temp_sbuf[((i - i_min) * j_length) + j] = matrix[i][j];
        }
    }

    MPI_Gatherv(temp_sbuf, (i_max - i_min) * j_length, MPI_DOUBLE, temp_rbuf, gather_recv_counts, gather_displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (process_id == 0) {
        for (i = 0; i < i_length; i++) {
            for (j = 0; j < j_length; j++) {
                matrix[i][j] = temp_rbuf[(i * j_length) + j];
            }
        }
    }
}

inline static
void send_adjacent_columns(double **matrix, int j_length) {
    MPI_Status status;

    if (process_id == 0) {
        MPI_Sendrecv(matrix[right_column], j_length, MPI_DOUBLE, process_id + 1, 0, matrix[right_column + 1], j_length, MPI_DOUBLE, process_id + 1, 0, MPI_COMM_WORLD, &status);
    } else if (process_id == (number_processes - 1)) {
        MPI_Sendrecv(matrix[left_column], j_length, MPI_DOUBLE, process_id - 1, 0, matrix[left_column - 1], j_length, MPI_DOUBLE, process_id - 1, 0, MPI_COMM_WORLD, &status);
    } else if (process_id % 2 == 1) {
        MPI_Sendrecv(matrix[left_column], j_length, MPI_DOUBLE, process_id - 1, 0, matrix[left_column - 1], j_length, MPI_DOUBLE, process_id - 1, 0, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(matrix[right_column], j_length, MPI_DOUBLE, process_id + 1, 0, matrix[right_column + 1], j_length, MPI_DOUBLE, process_id + 1, 0, MPI_COMM_WORLD, &status);
    } else {
        MPI_Sendrecv(matrix[right_column], j_length, MPI_DOUBLE, process_id + 1, 0, matrix[right_column + 1], j_length, MPI_DOUBLE, process_id + 1, 0, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(matrix[left_column], j_length, MPI_DOUBLE, process_id - 1, 0, matrix[left_column - 1], j_length, MPI_DOUBLE, process_id - 1, 0, MPI_COMM_WORLD, &status);
    }
}
#endif

inline static
void initialize_u_v_psi_omega(int i, int j, double **u, double **v, double **psi, double **omega, double *x, double *y, double A) {
    double sqrt_xy_squares, xi_A, yj_1;

    xi_A = x[i] + A;
    yj_1 = y[j] - 1;

//    sqrt_xy_squares = sqrt( pow(x[i], 2.0) + pow(y[j], 2.0) );
    sqrt_xy_squares = sqrt( (x[i] * x[i]) + (y[j] * y[j]) );

    omega[i][j] = 0.0;
    psi[i][j] = xi_A * yj_1;
    u[i][j] = xi_A / sqrt_xy_squares;
    v[i][j] = -yj_1 / sqrt_xy_squares;
}

inline static
void save_omega_psi(int Nx, int My, double **omega, double **omega_0, double **psi, double **psi_0) {
    int i, j;
    /**
     *  Save the current values of omega and psi
     */
    for (i = i_min; i < i_max; i++) {
        for (j = 0; j < (My + 2); j++) {
            omega_0[i][j] = omega[i][j];
            psi_0[i][j] = psi[i][j];
        }
    }
}


/**
 *  Finite difference approximation for vorticity
 */
inline static
void omega_calculation(int Nx, int My, double Cx2, double Cy2, double alpha, double alphaX, double alphaY, double **omega, double **omega_0, double **u, double **v, double *x, double *y, double **d_sq) {
    int i, j;
    double d2, dip1, dim1, djp1, djm1;

    for (i = i_min; i < i_max; i++) {
        for (j = 0; j < (My + 2); j++) {
            d_sq[i][j] = sqrt( (x[i] * x[i]) + (y[j] * y[j]) );
        }
    }
#ifdef _MPI_
    send_adjacent_columns(d_sq, My + 2);
    send_adjacent_columns(u, My + 2);
    send_adjacent_columns(v, My + 2);
    send_adjacent_columns(omega_0, My + 2);
#endif

    for (i = i_min_inner; i < i_max_inner; i++) {
        for (j = 1; j < (My + 1); j++) {
            d2 = d_sq[i][j] * d_sq[i][j];
            dip1 = d_sq[i + 1][j];
            dim1 = d_sq[i - 1][j];
            djp1 = d_sq[i][j + 1];
            djm1 = d_sq[i][j - 1];

            omega[i][j] =   omega_0[i][j] * (1 - alpha / d2);
            omega[i][j] += (omega_0[i + 1][j] * (-Cx2 * u[i + 1][j] * dip1 + alphaX) +
                            omega_0[i - 1][j] * ( Cx2 * u[i - 1][j] * dim1 + alphaX) +
                            omega_0[i][j + 1] * (-Cy2 * v[i][j + 1] * djp1 + alphaY) +
                            omega_0[i][j - 1] * ( Cy2 * v[i][j - 1] * djm1 + alphaY)) / d2;
        }
    }

//    fprintf(stderr, "did first loop\n");
/*
    for (j = My; j >= 1; j--) {
        i = 0;
        d = sqrt( (x[i] * x[i]) + (y[j] * y[j]) );
        omega[i][j] =   omega_0[i][j] * (1 + 3.0 * Cx2 * u[i][j] / d + 2.0 * alphaX / (d * d) - 2.0 * alphaY / (d * d)) +
                        omega_0[i + 1][j] * (-4.0 * Cx2 * u[i][j] / d - 5.0 * alphaX / (d * d)) +
                        omega_0[i + 2][j] * (Cx2 * u[i][j] / d + 4.0 * alphaX / (d * d)) +
                        omega_0[i + 3][j] * (-alphaX / (d * d)) +
                        omega_0[i][j + 1] * (-Cy2 * v[i][j] / d + alphaY / (d * d)) +
                        omega_0[i][j - 1] * ( Cy2 * v[i][j] / d + alphaY / (d * d));

        fprintf(stderr,"omega[%d][%d] = %lf\n", i, j, omega[i][j]);

        i = Nx + 1;
        d = sqrt( (x[i] * x[i]) + (y[j] * y[j]) );
        omega[i][j] =   omega_0[i][j] * (1 - 3.0 * Cx2 * u[i][j] / d + 2.0 * alphaX / (d * d) - 2.0 * alphaY / (d * d)) +
                        omega_0[i - 1][j] * (4.0 * Cx2 * u[i][j] / d - 5.0 * alphaX / (d * d)) +
                        omega_0[i - 2][j] * (-Cx2 * u[i][j] / d + 4.0 * alphaX / (d * d)) +                 // should Cx2 be negative?
                        omega_0[i - 3][j] * (-alphaX / (d * d)) +
                        omega_0[i][j + 1] * (-Cy2 * v[i][j] / d + alphaY / (d * d)) +
                        omega_0[i][j - 1] * ( Cy2 * v[i][j] / d + alphaY / (d * d));
    }
*/
    /*
     *
    fprintf(stderr, "omega:\n");
    for (i = 0; i < Nx + 2; i++) {
        for (j = 0; j < My + 2; j++) {
            fprintf(stderr, " %.9lf\n", omega[i][j]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n\n");
    */

}


inline static
void psi_calculation(int Nx, int My, double kappa2, double kappaA, double dxx, double **psi, double **psi_calc_0, double **omega, double k_psi, double *x, double *y, double tolerance_level, double** d_sq) {
    int i, j, kPsi;
    double PsiTol, psi_diff;
    double PsiMax;

    for (i = i_min_inner; i < i_max_inner; i++) {
        for (j = 1; j < (My + 1); j++) {
            d_sq[i][j] = dxx * omega[i][j] * ((x[i] * x[i]) + (y[j] * y[j]));
        }
    }
#ifdef _MPI_
    send_adjacent_columns(d_sq, My + 2);
#endif

    kPsi = 0;
    do {
        kPsi += 1;

//        fprintf(stderr, "%d -- PsiTol: %.10lf\n", kPsi, PsiTol);

        for (i = i_min; i < i_max; i++) {
            for (j = 0; j < (My + 2); j++) {
                psi_calc_0[i][j] = psi[i][j];
            }
        }
#ifdef _MPI_
        send_adjacent_columns(psi_calc_0, My + 2);
#endif

        PsiTol = 0;
        for (i = i_min_inner; i < i_max_inner; i++) {
            for (j = 1; j < (My + 1); j++) {
                psi[i][j] = kappaA * (d_sq[i][j] + psi_calc_0[i + 1][j] + psi_calc_0[i - 1][j] + kappa2 * (psi_calc_0[i][j + 1] + psi_calc_0[i][j - 1]));

                psi_diff = fabs(psi[i][j] - psi_calc_0[i][j]);

                if (psi_diff > PsiTol) PsiTol = psi_diff; 
            }
        }

#ifdef _MPI_
        MPI_Allreduce(&PsiTol, &PsiMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        PsiTol = PsiMax;
#endif
    } while (PsiTol > tolerance_level);

#ifdef _MPI_
//    if (process_id == 0) fprintf(stderr, "[rank: %d] kPsi: %d, PsiTol: %.10lf\n", process_id, kPsi, PsiTol);
#else
    fprintf(stderr, "kPsi: %d, PsiTol: %.10lf\n", kPsi, PsiTol);
#endif
}

inline static
void boundary_conditions(int Nx, int My, double **omega, double **psi, double **u, double **v, double *x, double *y, double dyy, double A, int ibl) {
    int i, j;
    double d6;

    /**
     *  Lower & Upper Boundary Conditions
     */
    for (i = i_min; i < i_max; i++) {
        /**
         *  Lower
         */
        j = 0;
        d6 = (x[i] * x[i]) + (y[0] * y[0]);
        psi[i][0] = 0.0;
        omega[i][0] = ((7.0 * psi[i][0]) - (8.0 * psi[i][1]) + psi[i][2]) / (2.0 * dyy) / d6;
        u[i][0] = 0.0;
        v[i][0] = 0.0;

        /**
         *  Upper
         */
        initialize_u_v_psi_omega(i, My + 1, u, v, psi, omega, x, y, A);
    }

    /**
     *  Side Boundary Conditions
     */
    for (j = 1; j < (My + 1); j++) {
        if (j > ibl) {
            initialize_u_v_psi_omega(0, j, u, v, psi, omega, x, y, A);
            initialize_u_v_psi_omega(Nx + 1, j, u, v, psi, omega, x, y, A);
        } else {
            i = 0;
            omega[i][j] = omega[i + 1][j];
            psi[i][j] = psi[i + 1][j];
            u[i][j] = u[i + 1][j];
            v[i][j] = v[i + 1][j];

            i = Nx + 1;
            omega[i][j] = omega[i - 1][j];
            psi[i][j] = psi[i - 1][j];
            u[i][j] = u[i - 1][j];
            v[i][j] = v[i - 1][j];
        }
    }

}

inline static
void velocity_calculation(int Nx, int My, double **psi, double **u, double **v, double *x, double *y, double dx2, double dy2) {
    int i, j;
    double sqrt_xy_squares;

#ifdef _MPI_
    send_adjacent_columns(psi, My + 2);
#endif

    for (i = i_min_inner; i < i_max_inner; i++) {
        for (j = 1; j < (My + 1); j++) {
            sqrt_xy_squares = sqrt( (x[i] * x[i]) + (y[j] * y[j]) );

            u[i][j] =  ((psi[i][j + 1] - psi[i][j - 1]) / dy2) / sqrt_xy_squares;
            v[i][j] = -((psi[i + 1][j] - psi[i - 1][j]) / dx2) / sqrt_xy_squares;
        }
    }

/*        for (j = 1; j < (My + 1); j++) {
        i = 0;
        sqrt_xy_squares = sqrt( (x[i] * x[i]) + (y[i] * y[i]) );
        u[i][j] = (psi[i][j+1] - psi[i][j-1]) / dx2 / sqrt_xy_squares;
        v[i][j] = ((-3.0 * psi[i][j]) + (4.0 * psi[i+1][j]) - psi[i+2][j]) / dx2 / sqrt_xy_squares;

        i = Nx + 1;
        sqrt_xy_squares = sqrt( (x[i] * x[i]) + (y[i] * y[i]) );
        u[i][j] = (psi[i][j+1] - psi[i][j-1]) / dx2 / sqrt_xy_squares;
        v[i][j] = ((-3.0 * psi[i][j]) + (4.0 * psi[i-1][j]) - psi[i-2][j]) / dx2 / sqrt_xy_squares;
    }
*/

}

inline static
void find_max_values(int Nx, int My, double **omega, double **omega_0, double **psi, double **psi_0, double *OmTol, double *PsiTol) {
    int i, j;
    double om_current, psi_current, psi_max, om_max;

    (*OmTol) = 0;
    (*PsiTol) = 0;

    /**
     *  Check max value change
     */
    for (i = i_min; i < i_max; i++) {
        for (j = 0; j < (My + 2); j++) {
            om_current = fabs(omega[i][j] - omega_0[i][j]);
            if (om_current > (*OmTol)) {
                (*OmTol) = om_current;
            }

            psi_current = fabs(psi[i][j] - psi_0[i][j]);
            if (psi_current > (*PsiTol)) {
                (*PsiTol) = psi_current;
            }
        }
    }

#ifdef _MPI_
        MPI_Allreduce(PsiTol, &psi_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(OmTol, &om_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        (*OmTol) = om_max;
        (*PsiTol) = psi_max;
#endif
}

void print_1d(const char *name, double *x, int length, unsigned short only_total, FILE *file) {
    int i;
    double sum = 0.0;

    if (!only_total) fprintf(file, "[%s]\n", name);
    for (i = 0; i < length; i++) {
        if (!only_total) fprintf(file, "%s[%d] = %5.15lf\n", name, i, x[i]);
        sum += x[i];
    }
    fprintf(file, "[%s] total: %5.15lf\n", name, sum);
}

void print_2d(const char *name, double **x, int length1, int length2, unsigned short only_total, FILE *file) {
    int i, j;
    double sum = 0.0;

    for (i = i_min; i < i_max; i++) {
        for (j = 0; j < length2; j++) {
            if (!only_total) fprintf(file, "%s[%d][%d] = %5.15lf\n", name, i, j, x[i][j]);
            sum += x[i][j];
        }
    }
    fprintf(file, "[%s] total: %5.15lf\n", name, sum);
}

void print_flow_field_initialization(FILE* out, double *x, double *y, double **u, double **v, double **omega, double **psi, double Nx, double My) {
    int i, j;

#ifdef _MPI_
    if (process_id != 0) return;
#endif

    print_1d("x", x, Nx + 2, TRUE, stderr);
    print_1d("y", y, My + 2, TRUE, stderr);
    print_2d("u", u, Nx + 2, My + 2, TRUE, stderr);
    print_2d("v", v, Nx + 2, My + 2, TRUE, stderr);
    print_2d("psi", psi, Nx + 2, My + 2, TRUE, stderr);
    print_2d("omega", omega, Nx + 2, My + 2, TRUE, stderr);
}

inline static
short read_checkpoint(char *filename, int *iteration, double **omega, double **psi, double **u, double **v, int length1, int length2) {
    int i;
    FILE *file;

    file = fopen(filename, "r");
    if (file == NULL) return FALSE;

    fread(iteration, sizeof(int), 1, file);
    for (i = 0; i < length1; i++) {
        fread(omega[i], sizeof(double), length2, file);
        fread(psi[i], sizeof(double), length2, file);
        fread(u[i], sizeof(double), length2, file);
        fread(v[i], sizeof(double), length2, file);
    }
    fclose(file);

    return TRUE;
}

inline static 
void write_checkpoint(char *filename, int iteration, double **omega, double **psi, double **u, double **v, int length1, int length2) {
    int i;
    FILE *file;

    file = fopen(filename, "w+");

    fwrite(&iteration, sizeof(int), 1, file);
    for (i = 0; i < length1; i++) {
        fwrite(omega[i], sizeof(double), length2, file);            //  vorticity how much things are swirling around
        fwrite(psi[i], sizeof(double), length2, file);              //  stream function (tangent to the velocity)
        fwrite(u[i], sizeof(double), length2, file);                //
        fwrite(v[i], sizeof(double), length2, file);                //
    }
    fflush(file);
    fclose(file);
}

#ifdef _WIN32
int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR Args, int WinMode){
    LPSTR command_line;
    char* argv[100];
    int argc;

    command_line = GetCommandLine();
    argc = parse_command_line( command_line, argv );
    return main(argc, argv);
}
#endif
