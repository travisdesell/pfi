	PROGRAM ParabolaFlow
*******************************************************************************
*******************************************************************************
*	This program solves the Vorticity-Streamline equations for the flow of  
*	an incompressible fluid around a canonic parabola at various modified
*	Reynolds Numbers and Circulation paramters
*
*	By Wallace J. Morris II
*******************************************************************************
*******************************************************************************

*	Assign all variables
*******************************************************************************
*	i,j,k,l,m,n	Counting Variables
*	Nx,My		Number of points in the X and Y directions
*	dx,dy		Grid spacing in X and Y directions
*	Ot,dt		Number of Points and spacing in time
* 	Xmin, Xmax	X Boundaries
*	Ymin, Ymax	Y Boundaries
*	x,y		Point Vectors
*	Omega		Matrix of Vorticies
*	Omega0		Previous Omega Matrix
*	Psi		Matrix of Streamlines
*	Psi0		Previous Psi Matrix
*	u,v		Velocity Component Matricies
*	Re		Reynolds Number
*	Rc		Cell Reynolds Number Rc = Re*dx/H
*	C		Courant Number
*	Kappa		Ratio of x and y spacing
*	OmTol		Max change in Omega per iteration
*	PsiTol		Max change in Psi per iteration

*	1-9		format designations
*	10-999		Loop designations
*	1000-9999	File Designations
*******************************************************************************

	character*45 data, stat, title, filename, index2, line, pfile, outfile,STRNG,idx
	character*3 string
	integer i,j,k,l,m,n,Nx,My,Ot,KPsi,index,kerr, kosc, report, LL, length, 
     c	TT, yy,psave,tic,ct,rcc,flag,IBL
	double precision Xmin,Xmax,Ymin,Ymax,xi,Umax,Re,C,
     c	dx,dx2,dxx,dy,dy2,dyy,dt,Kappa2,KappaA,Rc,Cx,Cx2,Cy,Cy2,
     c	alphaX,alphaY,alpha,Tol,OmTol,PsiTol,
     c	Amp,pi,Lambda,freq,yhalf,t,A, d3

*	Parameters to be adjusted
*******************************************************************************
	parameter (Xmin=-20, Xmax=20, Ymin=1, Ymax=11)
	parameter (Nx=199, My=399, pi = 3.14159265358979, Umax = 1)
*	parameter (Umax = 1, Re, A, Ot)


*	REMEMBER TO CHANGE THE SIZE OF THE FORMATTED DATA AS WELL



1	format(0X,(201D30.14))
2	format(0X,(201I30.1))
3	format(0X,(1D30.14))
4	format(0X,(1I30.1,2D30.14))
*******************************************************************************



	double precision x(Nx+2),y(My+2),Omega(Nx+2,My+2),Psi(Nx+2,My+2),
     c	u(Nx+2,My+2),v(Nx+2,My+2),Omega0(Nx+2,My+2),Psi0(Nx+2,My+2),
     c	OutDP(Nx+2),d6,PsiCalc0(Nx+2,My+2)
	integer OutIN(Nx+2)

*	open (unit = 3, file = 'flowerr500_osc1', status = 'old')



	


*	User Input Parameters
*******************************************************************************

1000	print *,' '
	print *,'Welcome to Parabola Flow Interactive'
	print *,''
	tic=0	
	ct=0
	flag=0
	report=0

	print *,'What would you like to do?'
	print *,'  1) Start new Simulation: ',
     c	' Mesh Dim =',Nx,' x',My
	print *,'  2) Continue Previous simulation '
	print *,'  3) Exit'
	read *,index

	if (index.eq.1) then
		print *,'What Reynolds Number would you like to run?'
		read *,Re
		print *,'What value of circulation parameter (A-Tilde) would you like to use?'
		read *,A

		print *,'A is',A
*		A = 1.8
		print *,'A is',A

		print *,'How many time steps would you like to run?'
		read *,Ot
		print *,'How many time steps between reports?'
		read *,report
		print *,'dt = ?'
		read *,dt
		print *,'To what tolerance level would you like to iterate?'
		read *,Tol
		print *,''
		print *,'There are ',My,' grid lines in the vertical direction.'
		print *,'How many do you want to be treated with Boundary Layer BCs?'
		read *,IBL
		print *,''
		print *,'And what would you like to call the output file?'
		read *,filename
		LL=length(filename,45)
*		print *,'Filename is ',LL,' characters long'
		pfile=filename
		TT=LL+3
		pfile(LL+1:)='P'

*		print *,pfile
		print *,'How many iterations between incremental file writes?'
		print *,'(Use negative number to turn off incremtnal save)'
		read *,psave



		print *,"A is: ",A, "Re is: ",Re



*		Generate the grid vectors and calculate dx and dy


			call linspace(Xmin,Xmax,Nx+2,x,dx)
			dx2 = 2*dx
			dxx = dx*dx

			call linspace(Ymin,Ymax,My+2,y,dy)
			dy2 = 2*dy
			dyy = dy*dy

			Kappa2 = (dx/dy)**2.0
			KappaA = 1.0/(2.0*(1+Kappa2))
			Rc = Re*dx

*			dt = .5*(Rc*dx)/4.0
*		*	dt = 0.5*(4.0/Re)
*		*	dt = 0.0006*(Re/2)/(1/(dxx)+1/(dyy))
*			dt = 0.0005
*			IBL=(My+1)/20

			Cx = dt/dx
			CX2 = .5*Cx
			Cy = dt/dy
			Cy2 = .5*Cy

			if (Cx.gt.Cy) then
				C = Cx
			else
				C = Cy
			endif

			alphaX = dt/(dxx*Re)
			alphaY = dt/(dyy*Re)
			alpha = 2*alphaX + 2*alphaY

			print *,'The Courant Number C =',C,' Must be less than 1'
			print *,''
			print *,'The Cell Reynolds Number Rc =',Rc,' MUST BE LESS THAN 4/C =', 4/C
		    	print *,''
			print *,'Grid Spacing dx, dy, dt:	',dx,dy,dt
			print *,''
			print *,'Continue? (y/n)'
			read *,index2

			if (index2.eq.'y') then
				print *,'Do you want to see screen output of the flow-field initialization? (y/n)'
				read *,idx
				goto 2000
			else
				goto 5000
			endif

	elseif (index.eq.2) then
		print *,' Enter the Previous Simulation File Name: '
		read *, data
		open (unit = 1, file = data, status = 'old')
		read (1,1) Omega
		read (1,1) Psi
		read (1,1) u
		read (1,1) v
		read (1,2) OutIN
		read (1,1) OutDP
		if (Nx.eq.OutIN(1)-2) then
			k = OutIN(3)
			kosc = 0
			kerr = k + 1
*			kerr = 1
			Re = OutDP(1) 
			Omtol = OutDP(2)
			PsiTol = OutDP(3) 
			dx = OutDP(4)
			dy = OutDP(5)
*			k  = OutIN(3)
			dt = OutDP(6)
			print *,' This file is compatible.'
			print *,''
			print *,'The Reynolds Number Re =',Re
			print *,''
			print *,'Currently A =',OutDP(12)
			print *,'What value of circulation parameter (A-Tilde) would you like to use?'
			read *,A
			print *,''
			print *,'IBL is =',OutIN(5)
			print *,'How many grid lines do you want to be treated with Boundary Layer BCs?'
			read *,IBL
			print *,''
			print *,'Input a new Tolerance value equal or less than ', OutDP(7)
			read *,Tol
			print *,''
			print *,'Currently dt =',OutDP(6)
1020			print *,'What value of time step (dt) would you like to use?'
			read *,dt



*			Generate the grid vectors and calculate dx and dy


				call linspace(Xmin,Xmax,Nx+2,x,dx)
				dx2 = 2*dx
				dxx = dx*dx

				call linspace(Ymin,Ymax,My+2,y,dy)
				dy2 = 2*dy
				dyy = dy*dy

				Kappa2 = (dx/dy)**2.0
				KappaA = 1.0/(2.0*(1+Kappa2))
				Rc = Re*dx

*			*	dt = .5*(Rc*dx)/4.0
*			*	dt = 0.5*(4.0/Re)
*			*	dt = 0.0006*(Re/2)/(1/(dxx)+1/(dyy))
*				dt = 0.0005
*				IBL=(My+1)/20

				Cx = dt/dx
				CX2 = .5*Cx
				Cy = dt/dy
				Cy2 = .5*Cy

				if (Cx.gt.Cy) then
					C = Cx
				else
					C = Cy
				endif

				alphaX = dt/(dxx*Re)
				alphaY = dt/(dyy*Re)
				alpha = 2*alphaX + 2*alphaY




			if (Rc.gt.4/C) then
				print *,'The Courant Number C =',C,' Must be less than 1 and'
				print *,'The Cell Reynolds Number Rc =',Rc,' MUST BE LESS THAN 4/C =', 4/C
				goto 1020
			endif	

			print *,''
			print *,'There have been ',k,' time steps for this model.'
			print *,'How many more would you like to perfom?'
			read *,Iter
			Ot = Iter + k
			print *,'How many time steps between reports?'
			read *,report
			print *,''
			print *,'And what would you like to call the output file?'
			read *,filename
			LL=length(filename,45)
*			print *,'Filename is ',LL,' characters long'
			pfile=filename
			TT=LL+3
			pfile(LL+1:)='P'
*			print *,pfile
			print *,'How many iterations between incremental file writes?'
			print *,'(Use negative number to turn off incremtnal save)'
			read *,psave	

				print *,'The Courant Number C =',C,' Must be less than 1'
				print *,''
				print *,'The Reynolds Number Re =',Re
				print *,'The Cell Reynolds Number Rc =',Rc,' MUST BE LESS THAN 4/C =', 4/C
     				print *,''
				print *,'Grid Spacing dx, dy, dt:	',dx,dy,dt
				print *,''
				print *,'Continue to ',Ot,'? (y/n)'
				read *,index2

				if (index2.eq.'y') then
					goto 3000
				else
					goto 5000
				endif



		else
			print *,' This file is not compatible.'
			print *,' Exit to change program parameters'
			print *,''
			goto 1000
		endif
	else
		goto 5000
	endif



*******************************************************************************
*******************************************************************************
*Set up the Initial Conditions
*******************************************************************************
*******************************************************************************

*Generate the y velocity component matrix v
**************************************************
2000	call ZEROS(Nx+2,My+2,1,v)
	do 5 i = 1,Nx+2
		do 6 j = 1,My+2
		v(i,j) = -(y(j)-1)/sqrt(x(i)**2+y(j)**2)
6		enddo
5	enddo		


*Generate the x velocity component matrix u
**************************************************

	print*,"A is: ",A
	call ZEROS(Nx+2,My+2,1,u)
	do 10 i = 1,Nx+2
		do 20 j = 2,My+2
		u(i,j) = (x(i)+A)/sqrt(x(i)**2+y(j)**2)
*		print*,"u[",i,"][",j,"] = (",x(i)," + ",A," / ",sqrt(x(i)**2+y(j)**2)," = ",u(i,j)
20		enddo
10	enddo


*	Generate the Streamline matrix Psi
******************************************
	call ZEROS(Nx+2,My+2,1,Psi)
	do 30 i = 1,Nx+2
		do 40 j = 1,My+2
		Psi(i,j) = (x(i)+A)*(y(j)-1)

40		enddo
30	enddo


*	Generate the Vorticity matrix Omega
*******************************************
	call ZEROS(Nx+2,My+2,1,Omega)
	do 45 i = 1,Nx+2
		d6=x(i)**2+y(1)**2
		Omega(i,1) = (7.0*Psi(i,1)-8.0*Psi(i,2)+Psi(i,3))/(2.0*dyy)/d6
45      enddo

	if (idx.eq.'y') then
        print*,"x:"
		call PRINTVECTOR(Nx+2,x)
        print*, ""

        print*,"y:"
		call PRINTVECTOR(My+2,y)
        print*,""

        print*,"u:"
		call PRINTMATRIX(Nx+2, My+2, 1, u)
        print*,""

        print*,"v:"
		call PRINTMATRIX(Nx+2, My+2, 1, v)
        print*,""

        print*,"psi:"
		call PRINTMATRIX(Nx+2, My+2, 1, Psi)
        print*,""

        print*,"omega:"
		call PRINTMATRIX(Nx+2, My+2, 1, Omega)
        print*,""
	endif
		call FDATE(STRNG)
		print*,'Flow-field finished initializing at ',STRNG

*****************************************************************************
*****************************************************************************
*	Perform the Numerical Calculations for each time step k
*****************************************************************************
*****************************************************************************
	k = 0
	kerr = 1
6000	continue
*	write (3,4) k,OmTol,PsiTol
3000	k = k+1
*	print*, 'k= ',k

	OmTol = 0
	PsiTol = 0

	do 72 i = 1,Nx+2
		do 71 j = 1,My+2
			Omega0(i,j) = Omega(i,j)
			Psi0(i,j) = Psi(i,j)
71		continue
72	continue

*	Omega,Psi, and Velocity Calculations
****************************************************
	
*	print*,"before OmegaCalc"
*	print*,"omega:"
*	call PRINTMATRIX(Nx+2, My+2, 1, Omega)

	call OmegaCalc(Nx,My,Cx2,Cy2,alpha,alphaX,alphaY,Omega,
     c	Omega0,u,v,x,y)

*	print*,"after OmegaCalc"
*	print*,"omega:"
*	call PRINTMATRIX(Nx+2, My+2, 1, Omega)

	call PsiCalc(Nx,My,Kappa2,KappaA,dxx,Psi,PsiCalc0,Omega,kPsi,x,y,Tol)

*	print*,"after PsiCalc"
*	print*,"psi:"
*	call PRINTMATRIX(Nx+2, My+2, 1, Psi)

	t = k*dt
*

*	Lower & Upper BCs
***************************************************
  
	do 87 i = 1,Nx+2
*	Lower 
		j=1
		d6=x(i)**2+y(1)**2
		Psi(i,1) = 0.0
		Omega(i,1) = (7.0*Psi(i,1)-8.0*Psi(i,2)+Psi(i,3))/(2.0*dyy)/d6
		u(i,1) = 0.0
		v(i,1) = 0.0

*	Upper 
		j=My+2
		Omega(i,My+2) = 0.0
		Psi(i,My+2) = (x(i)+A)*(y(j)-1)
		u(i,My+2) = (x(i)+A)/sqrt(x(i)**2+y(j)**2)
		v(i,My+2) =-(y(j)-1)/sqrt(x(i)**2+y(j)**2)

*		u(i,My+2) = (4.0*u(i,My+1)-u(i,My))/3.0
87	continue

*	Side BCs
*******************************************************

	i=1
	do 88 j = 2,My+1
		If(j.gt.IBL) go to 90
                Omega(i,j) = Omega(i+1,j)
		Psi(i,j) = Psi(i+1,j) 
		u(i,j) = u(i+1,j)
		v(i,j) = v(i+1,j)
		go to 88
		
90		Omega(i,j) = 0.0
		Psi(i,j) = (x(i)+A)*(y(j)-1)
		u(i,j) = (x(i)+A)/sqrt(x(i)**2+y(j)**2)
		v(i,j) =-(y(j)-1)/sqrt(x(i)**2+y(j)**2)
88	continue

	i=Nx+2
	do 89 j = 2,My+1
		If(j.gt.IBL) go to 91
                Omega(i,j) = Omega(i-1,j)
		Psi(i,j) = Psi(i-1,j) 
		u(i,j) = u(i-1,j)
		v(i,j) = v(i-1,j)
		go to 89
		
91		Omega(i,j) = 0.0
		Psi(i,j) = (x(i)+A)*(y(j)-1)
		u(i,j) = (x(i)+A)/sqrt(x(i)**2+y(j)**2)
		v(i,j) =-(y(j)-1)/sqrt(x(i)**2+y(j)**2)
89	continue

*	print*,'after lower / upper BC'

*		print*,Omega
*		print*,Psi

*	print*,""
*	print*,"after lower / upper /side BC"
*	print*,"psi:"
*	call PRINTMATRIX(Nx+2, My+2, 1, Psi)

*	print*,""
*	print*,"omega:"
*	call PRINTMATRIX(Nx+2, My+2, 1, Omega)

*	print*,""
*	print*,"u:"
*	call PRINTMATRIX(Nx+2, My+2, 1, u)

*	print*,""
*	print*,"v:"
*	call PRINTMATRIX(Nx+2, My+2, 1, v)

*	print*,"before velocities calculation: dy2 [",dy2,"] dx2 [",dx2,"]"
*	print*,"u:"
*	call PRINTMATRIX(Nx+2, My+2, 1, u)
*	print*,"v:"
*	call PRINTMATRIX(Nx+2, My+2, 1, v)


*	Calculate velocities
********************************************************
	do 100 i = 2,Nx+1
	do 110 j = 2,My+1
*		call UCalc(Nx,My,i,j,dy2,Psi,u,x,y)
*		call VCalc(Nx,My,i,j,dx2,Psi,v,x,y)
		d3=sqrt(x(i)**2+y(j)**2)
		u(i,j) = (Psi(i,j+1)-Psi(i,j-1))/(dy2)/d3
		v(i,j) = -(Psi(i+1,j)-Psi(i-1,j))/(dx2)/d3
110	continue
100	continue

*	print*,"after velocities calculation:"
*	print*,"u:"
*	call PRINTMATRIX(Nx+2, My+2, 1, u)
*	print*,"v:"
*	call PRINTMATRIX(Nx+2, My+2, 1, v)


*	i=1
*	do 111 j=2,My+1
*	d8=sqrt(x(i)**2+y(j)**2)
*	u(i,j)=(Psi(i,j+1)-Psi(i,j-1))/dx2/d8
*	v(i,j)=-(-3.0*Psi(i,j)+4.0*Psi(i+1,j)-Psi(i+2,j))/dx2/d8
*111     continue
*
*	i=Nx+2
*	do 112 j=2,My+1
*	d9=sqrt(x(i)**2+y(j)**2)
*	u(i,j)=(Psi(i,j+1)-Psi(i,j-1))/dx2/d9
*	v(i,j)=-(3.0*Psi(i,j)-4.0*Psi(i-1,j)+Psi(i-2,j))/dx2/d9
*112     continue

*	print*,'u:'
*	call PRINTMATRIX(Nx+2, My+2, 1, u)
*	print*,''

*	print*,'v:'
*	call PRINTMATRIX(Nx+2, My+2, 1, v)
*	print*,''

*	print*,'Calc velocities'

*	Check max value change
********************************************************

	do 120 i = 1,Nx+2
	do 130 j = 1,My+2
		if (abs(Omega(i,j)-Omega0(i,j)).gt.OmTol) then
			OmTol = abs(Omega(i,j)-Omega0(i,j))
		endif
		if (abs(Psi(i,j)-Psi0(i,j)).gt.PsiTol) then
			PsiTol = abs(Psi(i,j)-Psi0(i,j))
		endif

130 	continue
120	continue

	print*,'k: ', k, ' OmTol: ',OmTol,' PsiTol: ',PsiTol


*        print*,"x:"
*		call PRINTVECTOR(Nx+2,x)

*        print*,"y:"
*		call PRINTVECTOR(My+2,y)

*        print*,"u:"
*		call PRINTMATRIX(Nx+2, My+2, 1, u)

*        print*,"v:"
*		call PRINTMATRIX(Nx+2, My+2, 1, v)

*        print*,"psi:"
*		call PRINTMATRIX(Nx+2, My+2, 1, Psi)

*        print*,"omega:"
*		call PRINTMATRIX(Nx+2, My+2, 1, Omega)


***********************************************************************
*	Output iterations and tolerances
***********************************************************************

	
*	Modifies filename for IDW
***********************************************************************
	tic=tic+1
	if (tic.eq.psave) then
		if(k.ne.Ot) then
		tic=0
		ct=ct+1
		
		call btd(ct,string,3,rcc)

		do 19 yy=2,4
			pfile(LL+yy:)=string(yy-1:)
19 			continue
		print *,'Writing Incrimental File=',pfile
		outfile=pfile

*		goto 4000
*690
		flag=1
		
		
*		print *,'Incremental file saved.'
		print *,''
*					
		endif
	endif

	if (flag .eq. 1) goto 4000
	
690	flag=0


		
*	print*,'OmTol=',OmTol,'PsiTol =', PsiTol
***********************************************************************	
	

	if (k.lt.Ot) then	                                       
	    if (k.eq.kerr) then
	    	call FDATE(STRNG)
			print *,k,Kpsi,OmTol,PsiTol,kerr,' ',STRNG
			kerr = kerr + report
			goto 6000
		endif
		goto 3000
	endif
	outfile = filename
	print *,'Writing output file'
	
	
*	Write results to file
***********************************************************************
*	title = 'Re'

    

4000	open (unit = 2, file = outfile, status = 'new')

	write (2,1) Omega
	write (2,1) Psi
	write (2,1) u
	write (2,1) v

	do 155 i=1,Nx+2
	OutIn(i)=0
	OutDP(i)=0
155	CONTINUE

	OutIN(1) = Nx+2
	OutIN(2) = My+2
	OutIN(3) = k
	OutIN(4) = Kerr
	OutIN(5) = IBL
	OutDP(1) = Re
	OutDP(2) = Omtol
	OutDP(3) = PsiTol
	OutDP(4) = dx
	OutDP(5) = dy
	OutDP(6) = dt
	OutDP(7) = Tol
	OutDP(8) = xmin
	OutDP(9) = xmax
	OutDP(10) = ymin
	OutDP(11) = ymax
	OutDP(12) = A

	write (2,2) OutIN
	write (2,1) OutDP
	
	call FDATE(STRNG)
	print*,STRNG
	
*	if (k.eq.OTN) then
*		print *, 'To continue to iterate press 1 to exit press 2'
*		print *, 'Note: Remember to change the name of your output file'
*		read *,index
*		if (index.eq.1) then
*			print *, 'How many more iterations would you like to complete?'
*			read *, Iter
*			OTN = OTN + Iter
*			goto 3000
*		endif
*	endif

*	if (k .lt. Ot) goto 5000

	print*,'k: ',k,', Ot: ',Ot

	if (flag .eq. 1) goto 690





*	End of Program
***********************************************************************

5000	END



***********************************************************************
*	Finite difference approximation for vorticity
***********************************************************************
*	
*
*
*
***********************************************************************
	subroutine OmegaCalc(Nx,My,Cx2,Cy2,alpha,alphaX,alphaY,Omega,
     c	Omega0,u,v,x,y)
	integer Nx,My,i,j
	double precision Omega(Nx+2,My+2),U(Nx+2,My+2),v(Nx+2,My+2),x(Nx+2),y(My+2)
	double precision Omega0(Nx+2,My+2),Cx2,Cy2,alpha,alphaX,alphaY,d,d2,d3
	double precision dip1, dim1, djp1, djm1

	do 80 i = 2,Nx+1
		do 90 j = 2,My+1
	d=sqrt(x(i)**2+y(j)**2)
	dip1=sqrt(x(i+1)**2+y(j)**2)
	dim1=sqrt(x(i-1)**2+y(j)**2)
	djp1=sqrt(x(i)**2+y(j+1)**2)
	djm1=sqrt(x(i)**2+y(j-1)**2)
	Omega(i,j) = Omega0(i,j)*(1-alpha/(d*d)) + 
     c	Omega0(i+1,j)*(-Cx2*u(i+1,j)*dip1 + alphaX)/(d*d) +
     c	Omega0(i-1,j)*( Cx2*u(i-1,j)*dim1 + alphaX)/(d*d) +
     c	Omega0(i,j+1)*(-Cy2*v(i,j+1)*djp1 + alphaY)/(d*d) +
     c	Omega0(i,j-1)*( Cy2*v(i,j-1)*djm1 + alphaY)/(d*d)

90		continue
80	continue
	goto 75
	i=1
	do 60 j = My+1,2,-1
	d2=sqrt(x(i)**2+y(j)**2)
	Omega(i,j) = Omega0(i,j)*(1+3.0*Cx2*u(i,j)/d2+2.0*alphaX/(d2*d2)-2.0*alphaY/(d2*d2)) + 
     c	Omega0(i+1,j)*(-4.0*Cx2*u(i,j)/d2 - 5.0*alphaX/(d2*d2)) +
     c	Omega0(i+2,j)*(Cx2*u(i,j)/d2 + 4.0*alphaX/(d2*d2)) +
     c	Omega0(i+3,j)*(-alphaX/(d2*d2)) +
     c	Omega0(i,j+1)*(-Cy2*v(i,j)/d2 + alphaY/(d2*d2)) +
     c	Omega0(i,j-1)*(Cy2*v(i,j)/d2 + alphaY/(d2*d2))
	 
*		print *,"omega[",i,"][",j,"] = ",omega(i,j)
60	enddo

	i=Nx+2
	do 70 j = My+1,2,-1
	d3=sqrt(x(i)**2+y(j)**2)
	Omega(i,j) = Omega0(i,j)*(1-3.0*Cx2*u(i,j)/d3+2.0*alphaX/(d3*d3)-2.0*alphaY/(d3*d3)) + 
     c	Omega0(i-1,j)*(4.0*Cx2*u(i,j)/d3 - 5.0*alphaX/(d3*d3)) +
     c	Omega0(i-2,j)*(-Cx2*u(i,j)/d3 + 4.0*alphaX/(d3*d3)) +
     c	Omega0(i-3,j)*(-alphaX/(d3*d3)) +
     c	Omega0(i,j+1)*(-Cy2*v(i,j)/d3 + alphaY/(d3*d3)) +
     c	Omega0(i,j-1)*(Cy2*v(i,j)/d3 + alphaY/(d3*d3))
	
70	enddo
75	continue
	end


************************************************************************
*	Iterative Stream Function Routine
************************************************************************
*	Tol 		Tolerance level to cease iterations
*	PsiTol		Calculated difference in Psi values per iteration
*	Psi1m,Psi0m	Maximum single value of the previous Psi matrix
*	Psi0		Previous Psi Matrix Values
************************************************************************
	subroutine PsiCalc(Nx,My,Kappa2,KappaA,dxx,Psi,PsiCalc0,Omega,kPsi,
     c  x,y,Tol)
	integer Nx,My,Ot,i,j,kPsi,n,m
	double precision Omega(Nx+2,My+2),Psi(Nx+2,My+2),
     c  PsiCalc0(Nx+2,My+2),
     c	Kappa2,KappaA,dxx,Tol,PsiTol,Psi1m,Psi0m,x(Nx+2),y(My+2),d2,
     c  d12,d20

*	Tol = 1E-7
	PsiTol = 1
	Psi1m = 0
	kPsi = 0

10	if (abs(PsiTol) .gt. Tol) then
	kPsi = kPsi+1
*	print*,kPsi,' -- PsiTol: ',PsiTol
*	print*,'PsiTol: ',PsiTol
	Psi0m = Psi1m
	Psi1m = 0

	do 20 n = 1,Nx+2
	do 30 m = 1,My+2
		PsiCalc0(n,m) = Psi(n,m)
30	enddo
20	enddo


	do 40 i = 2,Nx+1
	do 50 j = My+1,2,-1
	d2=x(i)**2+y(j)**2
	Psi(i,j) = KappaA*(dxx*Omega(i,j)*d2 + PsiCalc0(i+1,j) +
     c	PsiCalc0(i-1,j) + Kappa2*(PsiCalc0(i,j+1) + PsiCalc0(i,j-1)))
*	print*,'psi[',i,'][',j,'] = ',Psi(i,j),', omega[',i,'][',j,'] = ',Omega(i,j)
	
	if (abs(Psi(i,j)-PsiCalc0(i,j)) .gt. Psi1m) then
		Psi1m = abs(Psi(i,j)-PsiCalc0(i,j))
	endif

50	enddo
40	enddo
	goto 90

	i=1
	do 60 j = My-1,2,-1
	d12=x(i)**2+y(j)**2
	Psi(i,j) = 1/(2*(1.0+Kappa2))*(-dxx*Omega(i,j)*d12 + 5.0*PsiCalc0(i+1,j) - 4.0*PsiCalc0(i+2,j) + PsiCalc0(i+3,j) + Kappa2*(5.0*PsiCalc0(i,j+1) - 4.0*PsiCalc0(i,j+2)) + PsiCalc0(i,j+3))

	if (abs(Psi(i,j)-PsiCalc0(i,j)) .gt. Psi1m) then
		Psi1m = abs(Psi(i,j)-PsiCalc0(i,j))
	endif
	 
60	enddo

	i=Nx+2
	do 70 j = My-1,2,-1
	d20=x(i)**2+y(j)**2
	Psi(i,j) = 1/(2*(1.0+Kappa2))*(-dxx*Omega(i,j)*d20 + 5.0*PsiCalc0(i-1,j) - 4.0*PsiCalc0(i-2,j) + PsiCalc0(i-3,j)+ Kappa2*(5.0*PsiCalc0(i,j+1) - 4.0*PsiCalc0(i,j+2)) + PsiCalc0(i,j+3))

	if (abs(Psi(i,j)-PsiCalc0(i,j)) .gt. Psi1m) then
		Psi1m = abs(Psi(i,j)-PsiCalc0(i,j))
	endif

70	enddo
	
90	PsiTol = abs(Psi1m)
	goto 10

	endif

	print*,''
	print*,'kPsi: ',kPsi,', PsiTol: ',PsiTol

	end




***********************************************************************
*	Finite difference approximation for u velocities
***********************************************************************
	subroutine UCalc(Nx,My,i,j,dy2,Psi,u,x,y)
	integer Nx,My,i,j
	double precision Psi(Nx+2,My+2),u(Nx+2,My+2),x(Nx+2),y(My+2)
	double precision dy2,d3
	d3=sqrt(x(i)**2+y(j)**2)
	u(i,j) = (Psi(i,j+1)-Psi(i,j-1))/(dy2)/d3

	end


***********************************************************************
*	Finite difference approximation for v velocities
***********************************************************************
	subroutine VCalc(Nx,My,i,j,dx2,Psi,v,x,y)
	integer Nx,My,i,j
	double precision Psi(Nx+2,My+2),v(Nx+2,My+2),x(Nx+2),y(My+2)
	double precision dx2,d4
	d4=sqrt(x(i)**2+y(j)**2)
	v(i,j) = -(Psi(i+1,j)-Psi(i-1,j))/(dx2)/d4

	end


***********************************************************************
*	Linspace Subroutine
***********************************************************************

	subroutine linspace(min,max,N,x,dx)
	integer N, i
	double precision min, max, dX, x(N)

	
	dX = abs(max-min)/(N-1)
	x(1) = min
	x(N) = max

	do 10 i = 2,N-1
		x(i) = min + dX*(i-1)
10	continue

	end



***********************************************************************
* A Subroutine for printing vectors
***********************************************************************

	subroutine PRINTVECTOR(N,x)
	integer i,j,N
	double precision x(N),s


	s = 0
	do 30 i = 1,N
*	print*,i, " -- ",x(i)
	s = s + x(i)
30	enddo
	print*,"total: ",s
	end



***********************************************************************
* A Subroutine for printing matrices
***********************************************************************

	subroutine PRINTMATRIX(N,M,O,x)
	integer i,j,N,M,O
	double precision x(N,M,O), s

	s = 0
	do 30 i = 1,N
		do 20 j = 1,M
			do 10 k = 1,O
*			print*,x(i,j,k)
			s = s + x(i,j,k)
10			enddo
20		enddo
30	enddo
	print*,"total: ",s
	end


***********************************************************************
* Generate the Zeros Matrix
***********************************************************************

	subroutine ZEROS(N,M,O,x)
	integer i,j,N,M,O
	double precision x(N,M,O)


	do 30 i = 1,N
		do 20 j = 1,M
			do 10 k = 1,O
			x(i,j,k) = 0
10			enddo
20		enddo
30	enddo
	end


*	String Length Reading Subroutine for Incremental-Data-Writer (IDW)
*	Code adapted from Michael Kupferschmid
***********************************************************************

	FUNCTION length(line,q)
	integer q, aa
	character*1 line(q)
	do 1 aa=q,1,-1
		length=aa
		if(line(aa).NE.' ') return
1	continue
	length=0
	return
	
	end
	
	
*	Converts count to string for IDW 
*	Code adapted from Michael Kupferschmid
***********************************************************************

	subroutine btd(ct, string, ls, rc)
	character*1 string(LS)
	integer*4 rc,ct
*	character*1 numerl(10)/'0','1','2','3','4','5','6','7','8','9'/
	character*1 numerl(10)
	
	numerl(1) = '0'
	numerl(2) = '1'
	numerl(3) = '2'
	numerl(4) = '3'
	numerl(5) = '4'
	numerl(6) = '5'
	numerl(7) = '6'
	numerl(8) = '7'
	numerl(9) = '8'
	numerl(10) = '9'
	
	rc=1
	if (ls .le. 0) return


*	convert the abs val of the number and place it in string	

	m=iabs(ct)
	do 1 k=ls,1,-1
		string(k)=numerl(1+M-10*(m/10))
		m=m/10
1	continue

*	does the number fit in the space allowed?

	if (m .ne. 0) return

*	blank out extra zeros (from th e left)
	
	do 2 k=1,ls
		kfnb=k
		if (string(k) .ne. '0') goto 3
		string(k)='0'
2	continue
	
*	if th string is alll 0s make put the last one back

	string(ls)='0'
	
*	put in a negative sign if the number is <0

3	if (n .lt. 0 .and. kfnb .eq. 1) return
	if (n .lt. 0) string(kfnb-1)='-'
	rc=0
	
	return
	
	end
	
	


