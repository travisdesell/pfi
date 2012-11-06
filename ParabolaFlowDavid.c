#include<omp.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

/*  Parabola Flow
 *
 *  This program solves the Vorticity-Streamline equations for the flow of  
 *  an incompressible fluid around a canonic parabola at various modified
 *  Reynolds Numbers and Circulation parameters with or without a synthetic
 *  jet modification.
 *
 *  Authors:
 *
 *  Wallace J. Morris II - Original author of ParabolaFlow in Fortran77,
 *    added Omega,Psi, U, and V, boundary conditions, calculations and output,
 *    as well as the overarching program structure
 *
 *  Jean-Paul Reddinger - Added synthetic jet functionality, user defined
 *    parameters, variable mesh sizing, and output file reading
 *  
 *  Craig Hoffstein - Added synthetic jet functionality, user defined
 *    parameters, and variable mesh sizing
 *
 *  David McWilliams - Added synthetic jet functionality, converted code into C
 *
 *
 *  OUTPUT FILE INDEX
 *
 *  Omega		(Matrix of vorticity)
 *  Psi			(Matrix of streamlines)
 *  U			(Matrix of X velocities)
 *  V			(Matrix of Y velocities)
 *
 *  Integer Variables (OutIN)
 *    1.	Nx+2	(Total X direction points)
 *    2.	My+2	(Total Y direction points)
 *    3.	k	(Number of timesteps)
 *    4.	no longer used (previously kerr)
 *    5.	IBL	(Rows treated as within boundary layer (5% of total Y))
 *    6.	ia	(Jet start location)
 *    7.	ib	(Jet end location)
 *
 *  Double Precision Variables (OutDP)
 *    1.	Re	(Reynolds number)
 *    2.	Omtol	(Max change in omega per iteration)
 *    3.	PsiTol	(Max change in psi per iteration)
 *    4.	dx	(Grid spacing in X)
 *    5.	dy	(Grid spacing in Y)
 *    6.	dt	(Time-step size)
 *    7.	Tol	(Tolerance)
 *    8.	xmin	(Lower X boundary)
 *    9.	xmax	(Upper X boundary)
 *   10.	ymin	(Lower Y boundary)
 *   11.	ymax	(Upper Y boundary)
 *   12.	A	(Circulation parameter)
 *   13.	freq	(Jet frequency)
 *   14.	c0	(Jet strength)
 *   15.	xskew	(Placeholder for x direction mesh skewing)
 *   15.	yskew	(Placeholder for y direction mesh skewing)
 *
 *  Assign all variables
 *    i,j,k		Counting Variables
 *    Nx,My		Number of points in the X and Y directions
 *    dx,dy		Grid spacing in X and Y directions
 *    Ot,dt		Number of Points and spacing in time
 *    Xmin, Xmax 	X Boundaries
 *    Ymin, Ymax 	Y Boundaries
 *    x,y		Point Vectors
 *    Omega		Matrix of Vortices
 *    Omega0		Previous Omega Matrix
 *    Psi		Matrix of Streamlines
 *    Psi0		Previous Psi Matrix
 *    u,v		Velocity Component Matrices
 *    Re		Reynolds Number
 *    Rc		Cell Reynolds Number Rc = Re*dx/H
 *    C			Courant Number
 *    Kappa		Ratio of x and y spacing
 *    OmTol		Max change in Omega per iteration
 *    PsiTol		Max change in Psi per iteration
 *
 *
 *    VERSION INFO
 *    Version 3.3
 *
 *    Changes: 
 *    3.3               Added Parallelization using OpenMP
 *
 *    3.2               Uses Serial Optimization via DM and DM2 matrices
 *                      Uses pointers for Omega and Psi iterations
 *
 *    3.1               Original C conversion with jets
 */

// Function Prototypes
void linspace(const double min, const double max, const int N,
	      double* x, double* dx);
void OmegaCalc(const int Nx, int My, const double Cx2, const double Cy2,
	       const double alpha, const double alphaX, const double alphaY,
	       double** Omega, double** Omega0, double** u, double** v,
	       double** DM, double** DM2);
void PsiCalc(const int Nx, int My, const double Kappa2,
	     const double KappaA, const double dxx, double*** Psi,
	     double*** Psi0i, double** Omega, double** DM2, int* KPsi,
	     const double Tol);
int writeFile(double** Omega, double** Psi, double** u, double** v, int* OutIN,
	      double* OutDP, const int Nx, const int My, char* outfile);
void BCs (double** Omega, double** Psi, double** u, double** v, double* x,
	  double* y, double** DM, double** DM2, const double A, const int IBL,
	  const double dyy, const int Nx, const int My, const double dx,
	  const double t, const int ia, const int ib, const double c0,
	  const double freq);


int main() {

  const double Xmin=-20, Xmax=20, Ymin=1, Ymax=11;

  int Nx,My,Ot,report,IBL,psave,i,j,k,Omtol,ia,ib;

  int *OutIN;

  double Re,A,dt,Tol,Cx2,Cy2,alpha,alphaX,alphaY,Kappa2,KappaA,dx,dx2,dxx,
    dy,dy2,dyy,c0,freq;

  double *x,*y,*OutDP;

  double **Omega,**Psi,**u,**v,**Omega0,**Psi0,**Psi0i,**DM,**DM2;

  char filename[45],outfile[45];

  time_t now;

  //   User Input
  printf("\nWelcome to Parabola Flow Interactive\n\n");
  printf("What would you like to do?\n");
  printf("  1) Start new Simulation:\n");
  printf("  2) Continue Previous simulation\n");
  printf("  3) Exit\n");
  int index;
  scanf("%d",&index);

  // New Simulation
  if(index==1) {
    printf("Nx=?\n");
    scanf("%d",&Nx);
    printf("My=?\n");
    scanf("%d",&My);
    
    // General user inputs
    printf("What Reynolds Number would you like to run?\n");
    scanf("%lf",&Re);
    printf("What value of circulation parameter (A-Tilde),\n");
    printf("would you like to use?\n");
    scanf("%lf",&A);
    printf("How many time steps would you like to run?\n");
    scanf("%d",&Ot);
    printf("How many time steps between reports?\n");
    scanf("%d",&report);
    printf("dt = ?\n");
    scanf("%lf",&dt);
    printf("To what tolerance level would you like to iterate?\n");
    scanf("%lf",&Tol);
    printf("There are %12d grid lines in the vertical direction.\n",My);
    printf("How many do you want to be treated with,\n");
    printf("Boundary Layer BCs?\n");
    scanf("%d",&IBL);

    // User input for synthetic jets
    printf("Would you like to include a jet in the simulation? (y/n)\n");
    char index2[45];
    scanf("%s",index2);
    if(index2[0]=='y') {
      printf("What is the start location of the jet?\n");
      printf("Range of 1 to %d\n",Nx);
      scanf("%d",&ia);
      printf("What is the end location?\n");
      scanf("%d",&ib);
      printf("What is the amplitude of the jet?\n");
      scanf("%lf",&c0);
      printf("What is the frequency of the jet?\n");
      scanf("%lf",&freq);
    }
    else {
      freq=0;
      c0=0;
      ia=Nx+4;
      ib=Nx+4;
    }

    // Additional user input
    printf("What would you like to call the output file?\n");
    scanf("%s",filename);
    printf("How many iterations between incremental file writes?\n");
    printf("(use -1 to turn off incremental save)\n");
    scanf("%d",&psave);

    // Allocates arrays based on user defined mesh sizes
    x=malloc((Nx+2)*sizeof(double));
    y=malloc((My+2)*sizeof(double));
    Omega=malloc((My+2)*sizeof(double*));
    Psi=malloc((My+2)*sizeof(double*));
    u=malloc((My+2)*sizeof(double*));
    v=malloc((My+2)*sizeof(double*));
    Omega0=malloc((My+2)*sizeof(double*));
    Psi0=malloc((My+2)*sizeof(double*));
    Psi0i=malloc((My+2)*sizeof(double*));
    DM=malloc((My+2)*sizeof(double*));
    DM2=malloc((My+2)*sizeof(double*));
    OutDP=malloc((Nx+2)*sizeof(double));
    OutIN=malloc((Nx+2)*sizeof(int));
    for(i=0; i<My+2; ++i) {
      Omega[i]=malloc((Nx+2)*sizeof(double));
      Psi[i]=malloc((Nx+2)*sizeof(double));
      u[i]=malloc((Nx+2)*sizeof(double));
      v[i]=malloc((Nx+2)*sizeof(double));
      Omega0[i]=malloc((Nx+2)*sizeof(double));
      DM[i]=malloc((Nx+2)*sizeof(double));
      DM2[i]=malloc((Nx+2)*sizeof(double));
      Psi0[i]=malloc((Nx+2)*sizeof(double));
      Psi0i[i]=malloc((Nx+2)*sizeof(double));
    }

    // Generate grid vectors, Calculate dx and dy
    linspace(Xmin,Xmax,Nx+2,x,&dx);
    dx2=2*dx;
    dxx=dx*dx;

    linspace(Ymin,Ymax,My+2,y,&dy);
    dy2=2*dy;
    dyy=dy*dy;

    for(i=0; i<My+2; ++i)
      for(j=0; j<Nx+2; ++j) {
	DM2[i][j]=pow(x[j],2)+pow(y[i],2);
	DM[i][j]=sqrt(DM2[i][j]);
      }

    Kappa2=pow((dx/dy),2);
    KappaA=1/(2*(1+Kappa2));
    const double Rc=Re*dx;

    // Calculate Corant number
    const double Cx=dt/dx;
    Cx2=.5*Cx;
    const double Cy=dt/dy;
    Cy2=.5*Cy;

    double C;
    if(Cx>Cy)
      C=Cx;
    else
      C=Cy;

    alphaX=dt/(dxx*Re);
    alphaY=dt/(dyy*Re);
    alpha=2*alphaX+2*alphaY;

    // Check stability
    printf("The Courant Number C = %22.14le Must be less than 1\n\n",C);
    printf("The Cell Reynolds Number Rc = %20.14lf\n",Rc);
    printf("MUST BE LESS THAN 4/C = %20.14lf\n\n",4/C);
    printf("Grid Spacing dx = %20.14le, dy = %20.14le\n",dx,dy);
    printf("Continue? (y/n)\n");
    scanf("%s",index2);

    if(index2[0]!='y')
      return 0;

    // Initialize v, u, Psi, and Omega
    for(i=1; i<My+2; ++i)
      for(j=0; j<Nx+2; ++j) {
	v[i][j]=-(y[i]-1)/DM[i][j];
	u[i][j]=(x[j]+A)/DM[i][j];
	Psi[i][j]=(x[j]+A)*(y[i]-1);
	Omega[i][j]=0;
      }
    for(j=0; j<Nx+2; ++j) {
      v[0][j]=-(y[0]-1)/DM[0][j];
      u[0][j]=0;
      Psi[0][j]=(x[j]+A)*(y[0]-1);
      Omega[0][j]=((7*Psi[0][j]-8*Psi[1][j]+Psi[2][j])/(2*dyy))/DM2[0][j];
    }
    time(&now);
    printf("\nFlow-field finished initializing at %s\n",ctime(&now));
    k=0;
  }

  // Continuing a Simulation
  else if(index==2) {
    printf("Enter the Previous Simulation File Name:\n");
    char data[45];
    scanf("%s",data);
    FILE* f_in=fopen(data,"r");
    if(f_in==NULL) {
      printf("Error, file %s could not be opened for read.\n",data);
      return 0;
    }
    // Determine mesh size by looking back two lines from EOF
    int linelength=1;
    char c=getc(f_in);
    while(c!='\n') {
      c=getc(f_in);
      ++linelength;
    }
    fseek(f_in,-2*linelength,SEEK_END);
    int NxR,MyR;
    fscanf(f_in,"%d %d",&NxR,&MyR);
    Nx=NxR-2;
    My=MyR-2;

    // Allocates arrays
    x=malloc((Nx+2)*sizeof(double));
    y=malloc((My+2)*sizeof(double));
    Omega=malloc((My+2)*sizeof(double*));
    Psi=malloc((My+2)*sizeof(double*));
    u=malloc((My+2)*sizeof(double*));
    v=malloc((My+2)*sizeof(double*));
    Omega0=malloc((My+2)*sizeof(double*));
    Psi0=malloc((My+2)*sizeof(double*));
    Psi0i=malloc((My+2)*sizeof(double*));
    DM=malloc((My+2)*sizeof(double*));
    DM2=malloc((My+2)*sizeof(double*));
    OutDP=malloc((Nx+2)*sizeof(double));
    OutIN=malloc((Nx+2)*sizeof(int));
    for(i=0; i<My+2; ++i) {
      Omega[i]=malloc((Nx+2)*sizeof(double));
      Psi[i]=malloc((Nx+2)*sizeof(double));
      u[i]=malloc((Nx+2)*sizeof(double));
      v[i]=malloc((Nx+2)*sizeof(double));
      Omega0[i]=malloc((Nx+2)*sizeof(double));
      DM[i]=malloc((Nx+2)*sizeof(double));
      DM2[i]=malloc((Nx+2)*sizeof(double));
      Psi0[i]=malloc((Nx+2)*sizeof(double));
      Psi0i[i]=malloc((Nx+2)*sizeof(double));
    }

    // Input arrays and parameters from previous file
    rewind(f_in);
    for(i=0; i<My+2; ++i)
      for(j=0; j<Nx+2; ++j)
	fscanf(f_in,"%lg",&Omega[i][j]);
    for(i=0; i<My+2; ++i)
      for(j=0; j<Nx+2; ++j)
	fscanf(f_in,"%lg",&Psi[i][j]);
    for(i=0; i<My+2; ++i)
      for(j=0; j<Nx+2; ++j)
	fscanf(f_in,"%lg",&u[i][j]);
    for(i=0; i<My+2; ++i)
      for(j=0; j<Nx+2; ++j)
	fscanf(f_in,"%lg",&v[i][j]);
    for(i=0; i<Nx+2; ++i)
      fscanf(f_in,"%d",&OutIN[i]);
    for(i=0; i<Nx+2; ++i)
      fscanf(f_in,"%lg",&OutDP[i]);
    if(NxR==OutIN[0] && MyR==OutIN[1]) {
      k=OutIN[2];
      IBL=OutIN[4];
      Re=OutDP[0];
      dt=OutDP[5];
      Tol=OutDP[6];
      A=OutDP[11];

      // Print previous simulation parameters
      printf("\nThis file is compatible.\n\n");
      printf("The number of points in the x direction, Nx = %d\n",Nx);
      printf("The number of points in the y direction, My = %d\n",My);
      printf("The Reynolds number, Re = %22.14le\n",Re);
      printf("The circulation parameter, A = %22.14le\n",A);
      printf("The number of grid lines treated with boundary layer ");
      printf("BCs, IBL = %d\n",IBL);
      printf("The time step size, dt = %22.14le\n",dt);
      printf("The tolerance value, Tol = %22.14le\n\n",Tol);
      printf("Would you like to change any of the simulation ");
      printf("parameters? (y/n)\n");
      char index2[45];
      scanf("%s",index2);

      if(index2[0]=='y') {
	printf("What value of circulation parameter, ");
	printf("(A-Tilde) would you like to use?\n");
	scanf("%lf",&A);
	printf("How many grid lines do you want to ");
	printf("be treated with Boundary Layer BCs?\n");
	scanf("%d",&IBL);
	printf("\nInput a new Tolerance value equal or less than %22.14le\n",
	       Tol);
	scanf("%lf",&Tol);
	printf("\nCurrently dt = %22.14le\n",dt);
	printf("What value of time step (dt) would you like to use?\n");
	scanf("%lf",&dt);
      }

      // Print previous jet parameters
      if(OutDP[13]>0) {
	printf("\nThere is a jet in the simulation.\n");
	printf("The start location, ia = %d\n",OutIN[5]);
	printf("The end location, ib = %d\n",OutIN[6]);
	printf("The amplitude, c0 = %lf\n",OutDP[13]);
	printf("The frequency, freq = %lf\n\n",OutDP[12]);
      }
      else printf("There is no jet in the simulation\n");
      printf("Would you like to change the jet parameters? (y/n)\n");
      scanf("%s",index2);
      if(index2[0]=='y') {
	printf("What is the new start location of the jet?\n");
	printf("Max of %d\n",Nx);
	scanf("%d",&ia);
	printf("What is the new end location?\n");
	scanf("%d",&ib);
	printf("What is the new amplitude of the jet?\n");
	scanf("%lf",&c0);
	printf("What is the new frequency of the jet?\n");
	scanf("%lf",&freq);
      }
      else {
	ia=OutIN[5];
	ib=OutIN[6];
	c0=OutDP[13];
	freq=OutDP[12];
      }

      // Generate grid vectors, Calculate dx and dy
      linspace(Xmin,Xmax,Nx+2,x,&dx);
      dx2=2*dx;
      dxx=dx*dx;

      linspace(Ymin,Ymax,My+2,y,&dy);
      dy2=2*dy;
      dyy=dy*dy;

      for(i=0; i<My+2; ++i)
	for(j=0; j<Nx+2; ++j) {
	  DM2[i][j]=pow(x[j],2)+pow(y[i],2);
	  DM[i][j]=sqrt(DM2[i][j]);
	}

      Kappa2=pow((dx/dy),2);
      KappaA=1/(2*(1+Kappa2));
      const double Rc=Re*dx;

      // Calculate Corant number
      const double Cx=dt/dx;
      Cx2=.5*Cx;
      const double Cy=dt/dy;
      Cy2=.5*Cy;

      double C;
      if(Cx>Cy)
	C=Cx;
      else
	C=Cy;

      alphaX=dt/(dxx*Re);
      alphaY=dt/(dyy*Re);
      alpha=2*alphaX+2*alphaY;

      // Additional user input
      printf("There have been %d time steps for this model.\n",k);
      printf("How many more would you like to perform?\n");
      int Iter;
      scanf("%d",&Iter);
      Ot=Iter+k;
      printf("How many time steps between reports?\n");
      scanf("%d",&report);
      printf("What would you like to call the output file?\n");
      scanf("%s",filename);
      printf("How many iterations between incremental file writes?\n");
      printf("(use -1 to turn off incremental save)\n");
      scanf("%d",&psave);

      // Check stability
      printf("The Courant Number C = %22.14le Must be less than 1\n\n",C);
      printf("The Cell Reynolds Number Rc = %20.14lf\n",Rc);
      printf("MUST BE LESS THAN 4/C = %20.14lf\n\n",4/C);
      printf("Grid Spacing dx = %20.14le, dy = %20.14le\n",dx,dy);
      printf("Continue to %d time steps? (y/n)\n",Ot);
      scanf("%s",index2);
      if(index2[0]!='y')
	return 0;
    }
    else {
      printf("This file is not compatible.\n");
      printf("Exit to check file format\n\n");
      return 0;
    }
    fclose(f_in);
  }

  else
    return 0;

  // Begin Calculations for each time step, k
  int KPsi,tic1=0,tic2=0,ct=0;
  double OmTol,PsiTol,t;
  double** temp;
  while(k++<Ot) {
    OmTol=PsiTol=0;

#pragma omp parallel for
    for(i=1; i<My+1; ++i) {
      Omega0[i][0]=Omega[i][0];
      Omega0[i][Nx+1]=Omega[i][Nx+1];
      Psi0i[i][0]=Psi[i][0];
      Psi0i[i][Nx+1]=Psi[i][Nx+1];
    }

    memcpy(Omega0[0],Omega[0],(Nx+2)*sizeof(double));
    memcpy(Omega0[My+1],Omega[My+1],(Nx+2)*sizeof(double));
    memcpy(Psi0i[0],Psi[0],(Nx+2)*sizeof(double));
    memcpy(Psi0i[My+1],Psi[My+1],(Nx+2)*sizeof(double));

    temp=Omega0;
    Omega0=Omega;
    Omega=temp;

#pragma omp parallel for
    for(i=0; i<My+2; ++i)
      memcpy(Psi0[i],Psi[i],(Nx+2)*sizeof(double));

    // Omega and Psi Calculations
    OmegaCalc(Nx,My,Cx2,Cy2,alpha,alphaX,alphaY,Omega,Omega0,u,v,DM,DM2);
    PsiCalc(Nx,My,Kappa2,KappaA,dxx,&Psi,&Psi0i,Omega,DM2,&KPsi,Tol);

    // Boundary Conditions
    t=k*dt;
    BCs(Omega,Psi,u,v,x,y,DM,DM2,A,IBL,dyy,Nx,My,dx,t,ia,ib,c0,freq);

    // Calculate velocities
#pragma omp parallel default(shared)
    {
#pragma omp for
      for(i=1; i<My+1; ++i)
#pragma omp parallel shared(i, My)
	{
#pragma omp for
	  for(j=1; j<Nx+1; ++j) {
	    u[i][j]=(Psi[i+1][j]-Psi[i-1][j])/dy2/DM[i][j];
	    v[i][j]=-(Psi[i][j+1]-Psi[i][j-1])/dx2/DM[i][j];
	  }
	}
    }

    // Check max value change

#pragma omp parallel default(shared)
    {
#pragma omp for
      for(i=0; i<My+2; ++i)
#pragma omp parallel shared(i, My, OmTol, PsiTol)
	{
#pragma omp for
	  for(j=0; j<Nx+2; ++j) {
	    if(fabs(Omega[i][j]-Omega0[i][j])>OmTol)
	      OmTol=fabs(Omega[i][j]-Omega0[i][j]);
	    if(fabs(Psi[i][j]-Psi0[i][j])>PsiTol)
	      PsiTol=fabs(Psi[i][j]-Psi0[i][j]);
	  }
	}
    }

    if(k==Ot)
      break;

    // Write Incremental File      

    // Modifies filename
    if(++tic1==psave) {
      tic1=0;
      char outfile[45];
      sprintf(outfile,"%sP%03d",filename,++ct);
      printf("Writing Incremental File %s ",outfile);

      for(j=0; j<Nx+2; ++j) {
	OutIN[j]=0;
	OutDP[j]=0;
      }

      OutIN[0]=Nx+2;
      OutIN[1]=My+2;
      OutIN[2]=k;

      OutIN[4]=IBL;
      OutIN[5]=ia;
      OutIN[6]=ib;

      OutDP[0]=Re;
      OutDP[1]=OmTol;
      OutDP[2]=PsiTol;
      OutDP[3]=dx;
      OutDP[4]=dy;
      OutDP[5]=dt;
      OutDP[6]=Tol;
      OutDP[7]=Xmin;
      OutDP[8]=Xmax;
      OutDP[9]=Ymin;
      OutDP[10]=Ymax;
      OutDP[11]=A;
      OutDP[12]=freq;
      OutDP[13]=c0;

      if(writeFile(Omega,Psi,u,v,OutIN,OutDP,Nx,My,outfile)==-1) {
	printf("Error, file %s could not be opened for write.\n",outfile);
	return 0;
      }

      time(&now);
      printf("%s\n",ctime(&now));
    }

    // Print Report
    if(++tic2==report) {
      tic2=0;
      time(&now);
      printf("%d %d %22.14le %22.14le %s\n",
	     k,KPsi,OmTol,PsiTol,ctime(&now));
    }
  }

  // Writes Final Output File
  strcpy(outfile,filename);
  printf("Writing output file ");
  for(j=0; j<Nx+2; ++j) {
    OutIN[j]=0;
    OutDP[j]=0;
  }

  OutIN[0]=Nx+2;
  OutIN[1]=My+2;
  OutIN[2]=k;

  OutIN[4]=IBL;
  OutIN[5]=ia;
  OutIN[6]=ib;

  OutDP[0]=Re;
  OutDP[1]=OmTol;
  OutDP[2]=PsiTol;
  OutDP[3]=dx;
  OutDP[4]=dy;
  OutDP[5]=dt;
  OutDP[6]=Tol;
  OutDP[7]=Xmin;
  OutDP[8]=Xmax;
  OutDP[9]=Ymin;
  OutDP[10]=Ymax;
  OutDP[11]=A;
  OutDP[12]=freq;
  OutDP[13]=c0;

  if(writeFile(Omega,Psi,u,v,OutIN,OutDP,Nx,My,outfile)==-1) {
    printf("Error, file %s could not be opened for write.\n",outfile);
    return 0;
  }

  time(&now);
  printf("%s\n",ctime(&now));

  for(i=0; i<My+2; ++i) {
    free(Omega[i]);
    free(Psi[i]);
    free(u[i]);
    free(v[i]);
    free(Omega0[i]);
    free(Psi0[i]);
    free(Psi0i[i]);
    free(DM[i]);
    free(DM2[i]);
  }
  free(x);
  free(y);
  free(Omega);
  free(Psi);
  free(u);
  free(v);
  free(Omega0);
  free(Psi0);
  free(Psi0i);
  free(DM);
  free(DM2);
  free(OutDP);
  free(OutIN);

  return 0;
}


void linspace(const double min, const double max, const int N,
	      double* x, double* dx) {
  int i;
  *dx=(max-min)/(N-1);
  x[0]=min;
  x[N-1]=max;
  for(i=1; i<N-1; ++i)
    x[i]=min+(*dx)*i;
}

// Finite difference approximation for vorticity
void OmegaCalc(const int Nx, int My, const double Cx2, const double Cy2,
	       const double alpha, const double alphaX, const double alphaY,
	       double** Omega, double** Omega0, double** u, double** v,
	       double** DM, double** DM2) {
  int i,j;

#pragma omp parallel default(shared)
  {
#pragma omp for
    for(i=1; i<My+1; ++i)
#pragma omp parallel shared(i, My)
      {
#pragma omp for
	for(j=1; j<Nx+1; ++j) {
	  Omega[i][j]=Omega0[i][j]*(1-alpha/DM2[i][j])+
	    Omega0[i][j+1]*(-Cx2*u[i][j+1]*DM[i][j+1]+alphaX)/DM2[i][j]+
	    Omega0[i][j-1]*(Cx2*u[i][j-1]*DM[i][j-1]+alphaX)/DM2[i][j]+
	    Omega0[i+1][j]*(-Cy2*v[i+1][j]*DM[i+1][j]+alphaY)/DM2[i][j]+
	    Omega0[i-1][j]*(Cy2*v[i-1][j]*DM[i-1][j]+alphaY)/DM2[i][j];
	}
      }
  }
}

// Iterative Stream Function Routine
void PsiCalc(const int Nx, int My, const double Kappa2,
	     const double KappaA, const double dxx, double*** Psi,
	     double*** Psi0i, double** Omega, double** DM2, int* KPsi,
	     const double Tol) {
  int i,j;
  double PsiTol=1;
  double** temp;

  *KPsi=0;

  while(PsiTol>Tol) {
    ++*KPsi;
    PsiTol=0;

    temp=*Psi0i;
    *Psi0i=*Psi;
    *Psi=temp;

#pragma omp parallel default(shared)
    {
#pragma omp for
      for(i=1; i<My+1; ++i)
#pragma omp parallel shared(i, My, PsiTol)
	{
#pragma omp for
	  for(j=1; j<Nx+1; ++j) {
	    (*Psi)[i][j]=KappaA*(dxx*Omega[i][j]*DM2[i][j]+(*Psi0i)[i][j+1]+
				 (*Psi0i)[i][j-1]+Kappa2*((*Psi0i)[i+1][j]+
							  (*Psi0i)[i-1][j]));
	    if(fabs((*Psi)[i][j]-(*Psi0i)[i][j])>PsiTol)
	      PsiTol=fabs((*Psi)[i][j]-(*Psi0i)[i][j]);
	  }
	}
    }
  }
}

// Program output to file
int writeFile(double** Omega, double** Psi, double** u,double** v, int* OutIN,
	      double* OutDP, const int Nx, const int My, char* outfile) {
  int i,j;

  FILE* f_out=fopen(outfile,"w");

  if(f_out==NULL)
    return -1;

  for(i=0; i<My+2; ++i) {
    for(j=0; j<Nx+2; ++j)
      fprintf(f_out,"%30.13le",Omega[i][j]);
    fprintf(f_out,"\n");
  }

  for(i=0; i<My+2; ++i) {
    for(j=0; j<Nx+2; ++j)
      fprintf(f_out,"%30.13le",Psi[i][j]);
    fprintf(f_out,"\n");
  }

  for(i=0; i<My+2; ++i) {
    for(j=0; j<Nx+2; ++j)
      fprintf(f_out,"%30.13le",u[i][j]);
    fprintf(f_out,"\n");
  }

  for(i=0; i<My+2; ++i) {
    for(j=0; j<Nx+2; ++j)
      fprintf(f_out,"%30.13le",v[i][j]);
    fprintf(f_out,"\n");
  }

  for(j=0; j<Nx+2; ++j)
    fprintf(f_out,"%30d",OutIN[j]);
  fprintf(f_out,"\n");

  for(j=0; j<Nx+2; ++j)
    fprintf(f_out,"%30.13le",OutDP[j]);

  fclose(f_out);
  return 0;
}

// Boundary Conditions
void BCs (double** Omega, double** Psi, double** u, double** v, double* x,
	  double* y, double** DM, double** DM2, const double A, const int IBL,
	  const double dyy, const int Nx, const int My, const double dx,
	  const double t, const int ia, const int ib, const double c0,
	  const double freq) {
  int i,j;

  const double pi=3.14159265358979;
  const double amewa=(ia-((Nx+1)/2+1))*dx;
  const double f=sin(2*pi*freq*t);

  // Upper and Lower BCs
  for(j=0; j<Nx+2; ++j) {
    Psi[0][j]=0;
    if(j>=ia-1 && j<=ib-1)
      Psi[0][j]=(-c0*(-0.5*amewa*sqrt(pow(amewa,2)+1)-0.5*sinh(amewa)
		      +0.5*x[j]*sqrt(pow(x[j],2)+1)+0.5*sinh(x[j])))*f;
    if(j>ib-1)
      Psi[0][j]=Psi[0][ib-1];
    Omega[0][j]=(7*Psi[0][j]-8*Psi[1][j]+Psi[2][j])/(2*dyy)/DM2[0][j];
    u[0][j]=0;
    v[0][j]=0;
    if(j>ia-1 && j<ib-1)
      v[0][j]=c0*f;
    if(j>=ia-2 && j<=ib)
      Omega[0][j]+=(v[0][j+1]*sqrt(pow(x[j+1],2)+1)
		    -v[0][j-1]*sqrt(pow(x[j-1],2)+1))/(2*dx)/DM2[0][j];
    Omega[My+1][j]=0;
    Psi[My+1][j]=(x[j]+A)*(y[My+1]-1);
    u[My+1][j]=(x[j]+A)/DM[My+1][j];
    v[My+1][j]=-(y[My+1]-1)/DM[My+1][j];
  }

  // Side BCs
  for(i=1; i<My+1; ++i) {
    if(i<IBL) {
      Omega[i][0]=Omega[i][1];
      Psi[i][0]=Psi[i][1];
      u[i][0]=u[i][1];
      v[i][0]=v[i][1];
      Omega[i][Nx+1]=Omega[i][Nx];
      Psi[i][Nx+1]=Psi[i][Nx];
      u[i][Nx+1]=u[i][Nx];
      v[i][Nx+1]=v[i][Nx];
    }
    else {
      Omega[i][0]=0;
      Psi[i][0]=(x[0]+A)*(y[i]-1);
      u[i][0]=(x[0]+A)/DM[i][0];
      v[i][0]=-(y[i]-1)/DM[i][0];
      Omega[i][Nx+1]=0;
      Psi[i][Nx+1]=(x[Nx+1]+A)*(y[i]-1);
      u[i][Nx+1]=(x[Nx+1]+A)/DM[i][Nx+1];
      v[i][Nx+1]=-(y[i]-1)/DM[i][Nx+1];
    }
  }
}
