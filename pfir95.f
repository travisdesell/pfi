	PROGRAM ParabolaFlow
*******************************************************************************
*******************************************************************************
*	This program solves the Vorticity-Streamline equations for the flow of  
*	an incompressible fluid around a canonic parabola at various modified
*	Reynolds Numbers and Circulation paramters with or without a synthetic
*	jet modification.
*
*	Authors:
*	Wallace J. Morris II - Original author of ParabolaFlow, added Omega,
*		Psi, U, and V, boundary conditions, calculations and output,
*		as well as the overarching program structure
*	Jean-Paul Reddinger - Added synthetic jet functionality, user defined
*		parameters, variable mesh sizing, and output file reading
*	Craig Hoffstein - Added synthetic jet functionality, user defined
*		parameters, and variable mesh sizing
*       David McWilliams - Added synthetic jet functionality
*******************************************************************************
*******************************************************************************

*******************************************************************************
*	OUTPUT FILE INDEX
*******************************************************************************

*	Omega		(Matrix of vorticity)

*	Psi		(Matrix of streamlines)

*	U		(Matrix of X velocities)

*	V		(Matrix of Y velocities)

*	Integer Variables (OutIN)
*******************************************************************************
*   1.	Nx+2		(Total X direction points)
*   2.	My+2		(Total Y direction points)
*   3.	k		(Number of timesteps)
*   4.	Kerr		(Incrimental output file counter)
*   5.	IBL		(Rows treated as within boundary layer (5% of total Y))
*   6.	ia		(Jet start location)
*   7.	ib		(Jet end location)

	
*	Double Precision Variables (OutDP)
*******************************************************************************
*   1.	Re		(Reynolds number)
*   2.	Omtol		(Max change in omega per iteration)
*   3.	PsiTol		(Max change in psi per iteration)
*   4.	dx		(Grid spacing in X)
*   5.	dy		(Grid spacing in Y)
*   6.	dt		(Time-step size)
*   7.	Tol		(Tolerance)
*   8. 	xmin		(Lower X boundary)
*   9.	xmax		(Upper X boundary)
*  10.	ymin		(Lower Y boundary)
*  11.	ymax		(Upper Y boundary)
*  12.	A		(Circulation parameter)
*  13.	freq		(Jet frequency)
*  14.	c0		(Jet strength)
*  15.  xskew		(Placeholder for x direction mesh skewing)
*  15.  yskew		(Placeholder for y direction mesh skewing)


*	Assign all variables
*******************************************************************************
*	i,j,k,l,m,n,o	Counting Variables
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


*	10-999		Loop designations
*	1000-9999	File Designations
*******************************************************************************
	implicit none

	character(len = 20) form, form1

	character(len=45) data,stat,title,filename,index2,line,pfile,
     c	outfile,STRNG,index3,string,getlength
	integer i,j,k,l,m,n,o,Nx,My,Ot,KPsi,index,kerr,kosc,report, Iter,
     c	LL, length, TT, yy, psave,tic, ct, rcc, flag, IBL,ia,ib, totrow
	double precision Xmin, Xmax, Ymin, Ymax, xi, Umax, Re, C,
     c	dx,dx2,dxx,dy,dy2,dyy,dt,Kappa2,KappaA,Rc,Cx,Cx2,Cy,Cy2,
     c	alphaX,alphaY,alpha,Tol,OmTol,PsiTol,d6,
     c	Amp,pi,Lambda,freq,yhalf,t,c0,amewa,f,A

*	Parameters to be adjusted
*******************************************************************************
	parameter (Xmin=-20, Xmax=20, Ymin=1, Ymax=11)
	parameter (pi = 3.14159265358979d0, Umax = 1)

*	Preallocated arrays	
******************************************************************************* 
	double precision, allocatable :: x(:),y(:),Omega(:,:),Psi(:,:),
     c	u(:,:),v(:,:),Omega0(:,:),Psi0(:,:),
     c	OutDP(:)
	integer, allocatable :: OutIN(:)

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
	print *,'  1) Start new Simulation'
	print *,'  2) Continue Previous simulation '
	print *,'  3) Exit'
	read *,index

	if (index.eq.1) then
		
		print *,'What size mesh would you like to run?'
		print *,'N ='
		read *,Nx
		print *,'M ='
		read *,My
		
		freq = 0
		c0 = 0
		ia = Nx + 4
		ib = Nx + 4
		
		form = '(1X,E30.14E3)'
		form1 = '(1X,I30.1)'
		write (form(2:4),'(I3)')(Nx+2)
		write (form1(2:4),'(I3)')(Nx+2)
		
		
		allocate (x(Nx+2))
		allocate (y(My+2))
		allocate (Omega(Nx+2,My+2))
		allocate (Psi(Nx+2,My+2))
		allocate (u(Nx+2,My+2))
		allocate (v(Nx+2,My+2))
		allocate (Omega0(Nx+2,My+2))
		allocate (Psi0(Nx+2,My+2))
		allocate (OutDP(Nx+2))
		allocate (OutIN(Nx+2))
		
		
		
1010		print *,'Would you like to include a jet in the simulation?
     c		(y/n)'
		read *,index3
			
		if(index3.eq.'y') then
			print *,'What is the start location of the jet?'
			print *,'Range of 0 to', Nx
			read *,ia
			print *,'What is the end location?'
			read *,ib
			print *,'What is the amplitude of the jet?'
			read *,c0
			print *,'What is the frequency of the jet?'
			read *,freq
		
		elseif (index3.eq.'n') then
			freq = 0
			c0 = 0
			ia = Nx+4
			ib = Nx+4
		else 
		  	goto 1010
		endif		
		
		print *,'What Reynolds Number would you like to run?'
		read *,Re
		print *,'What value of circulation parameter (A-Tilde),'
     		print*,'would you like to use?'
		read *,A
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
		print *,'How many do you want to be treated with,'
                print *,'Boundary Layer BCs?'
		read *,IBL
		print *,''
		print *,'And what would you like to call the output file?'
		read *,filename
		LL=length(filename,45)
		
		pfile=filename
		TT=LL+3
		pfile(LL+1:)='P'

		print *,'How many iterations between incremental file writes?'
		print *,'(Use negative number to turn off incremtnal save)'
		read *,psave
		
		
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
			print *,'The Cell Reynolds Number Rc =',Rc
     			print *,'MUST BE LESS THAN 4/C =', 4/C
		    	print *,''
			print *,'Grid Spacing dx, dy, dt:	',dx,dy,dt
			print *,''
			print *,'Continue? (y/n)'
			read *,index2
			
			if (index2.eq.'y') then
				goto 2000
			else
				goto 5000
			endif
		
	elseif (index.eq.2) then
		print *,' Enter the Previous Simulation File Name: '
		read *, data
		open (unit = 1, file = data, status = 'old')

		form = '(A10)'
		
		o = 0
1012		if (getlength.NE.'          ') then
			read (1,form) getlength
			o = o + 1
			goto 1012
		endif
		
		My = (o - 1)/4 - 2
		close(1)
		
		open (unit = 1, file = data, status = 'old')
		
		form = '(A30)'
		
		o = 0
1013		if (o.lt.4*(My +2)+1) then
			read (1,form) getlength
			o = o + 1
			goto 1013
		endif

		form = '(    E30.14E3)'
		form1 = '(    I30.1)'
		
		read(getlength,form1) Nx
		Nx = Nx - 2
		close(1)		
		
		open (unit = 1, file = data, status = 'old')		
	
		allocate (x(Nx+2))
		allocate (y(My+2))
		allocate (Omega(Nx+2,My+2))
		allocate (Psi(Nx+2,My+2))
		allocate (u(Nx+2,My+2))
		allocate (v(Nx+2,My+2))
		allocate (Omega0(Nx+2,My+2))
		allocate (Psi0(Nx+2,My+2))
		allocate (OutDP(Nx+2))
		allocate (OutIN(Nx+2))
		
		write (form(2:5),'(I3)')(Nx+2)
		write (form1(2:5),'(I3)')(Nx+2)
		
		read (1,form) Omega
		read (1,form) Psi
		read (1,form) u
		read (1,form) v
		read (1,form1) OutIN
		read (1,form) OutDP
		
		IBL = OutIN(5)
		ia = OutIN(6)
		ib = OutIN(7)
		k = OutIN(3)
		kosc = 0
		kerr = k + 1
		Re = OutDP(1) 
		Omtol = OutDP(2)
		PsiTol = OutDP(3) 
		dx = OutDP(4)
		dy = OutDP(5)
		dt = OutDP(6)
		Tol = OutDP(7)
		A = OutDP(12)
		freq = OutDP(13)
		c0 = OutDP(14)
		
		print *,'The number of points along the airfoil (Nx) is ',Nx
		print *,'The number of points in y direction (My) is ',My
		print *,'The Reynolds Number (Re) is ',Re
		print *,'The circulation parameter (A~) is ',OutDP(12)
		print *,'The number of grid lines treated with boundary'
		print *,'layer boundary conditions (IBL) is ',OutIN(5)
		print *,'The Tolerance (Tol) value is ',OutDP(7)
		print *,'The time step (dt) is ',OutDP(6)
		print *,'The PsiTol is ', OutDP(3)
		print *,'The Omtol is ', OutDP(2)
		print *,'The k is ', OutIN(3)
		print *,'The kerr is ',kerr
		print *,'dx is ',dx
		print *,'dy is ',dy
		print *,'dt is ',dt

1015		print *,'Would you like to change any of the simulation'
     		print *,'parameters? (y/n)'
     		read *,index3
     			
     		if(index3.eq.'y') then
     			
		print *,'What value of circulation parameter,
     c		(A-Tilde) would you like to use?'
		read *,A			
		print *,'How many grid lines do you want to,
     c		be treated with Boundary Layer BCs?'
		read *,IBL
		print *,'Input a new Tolerance value equal or 
     c		less than ',OutDP(7)
		read *,Tol			
1020		print *,'What value of time step (dt) would you like to use?'
		read *,dt	
		
		elseif(index3.eq.'n') then
		else
		goto 1015
		endif
		
		
		if(OutDP(14).gt.kosc) then
		print *,'There is a jet in the simulation.'
		print *,''
		
		print *,'The start location is ',OutIN(6),''
		print *,'The end location is ',OutIN(7),''
		print *,'The amplitude is ',OutDP(14)
		print *,'The frequency is ',OutDP(13)
		
		else
		print *,'There is no jet in the simulation'
		print *,''
		endif
		
1030		print *,'Would you like to change the jet parameters? (y/n)'
		read *,index3
		
		if(index3.eq.'y') then
			
			print *,'What is the new start location of the jet?'
			print *,'Max of',Nx
			read *,ia
			print *,'What is the new end location?'
			read *,ib
			print *,'What is the new amplitude of the jet?'
			read *,c0				
			print *,'What is the new frequency of the jet?'
			read *,freq
	
		elseif(index3.eq.'n') then
			freq = OutDP(13)
			c0 = OutDP(14)
			ia = OutIN(6)
			ib = OutIN(7)
		else
			goto 1030
		endif
			
*		Generate the grid vectors and calculate dx and dy

			
1035			call linspace(Xmin,Xmax,Nx+2,x,dx)
			dx2 = 2*dx
			dxx = dx*dx
		
			call linspace(Ymin,Ymax,My+2,y,dy)
			dy2 = 2*dy
			dyy = dy*dy
		
			Kappa2 = (dx/dy)**2.0
			KappaA = 1.0/(2.0*(1+Kappa2))
			Rc = Re*dx
				
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
			print *,'The Courant Number (C) is 
     c			',C,' Must be less than 1 and'
			print *,'The Cell Reynolds Number (Rc) is ',Rc,',
     c			MUST BE LESS THAN (4/C) ', 4/C
     			print *,'What value of time step (dt) would you like to use?'
			read *,dt
			goto 1035
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
		pfile=filename
		TT=LL+3
		pfile(LL+1:)='P'
		print *,'How many iterations between incremental file writes?'
		print *,'(Use negative number to turn off incremtnal save)'
		read *,psave	
		
			print *,'The Courant Number (C) is ',C,' Must be less than 1'
			print *,''
			print *,'The Reynolds Number (Re) is ',Re
			print *,'The Cell Reynolds Number (Rc) ',Rc,',
     c			MUST BE LESS THAN (4/C) ', 4/C
     			print *,''
			print *,'Grid Spacing dx, dy, dt:	',dx,dy,dt
			print *,''
1038			print *,'Continue to ',Ot,'? (y/n)'
			read *,index2
			
			if (index2.eq.'y') then
				goto 3000
			elseif (index2.eq.'n') then
				goto 1000
			else
			goto 1038
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
		v(i,j) = -(y(j)-1)/Dsqrt(x(i)**2+y(j)**2)
6		enddo
5	enddo		


*Generate the x velocity component matrix u
**************************************************
	call ZEROS(Nx+2,My+2,1,u)
	do 10 i = 1,Nx+2
		do 20 j = 2,My+2
		u(i,j) = (x(i)+A)/Dsqrt(x(i)**2+y(j)**2)
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
		Omega(i,1) = (7.0d0*Psi(i,1)-8.0*Psi(i,2)+Psi(i,3))/(2.0*dyy)/d6
45      enddo

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

3000	k = k+1

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
	
	call OmegaCalc(Nx,My,Cx2,Cy2,alpha,alphaX,alphaY,Omega,
     c	Omega0,u,v,x,y)
	call PsiCalc(Nx,My,Kappa2,KappaA,dxx,Psi,Omega,kPsi,x,y,Tol)

	t = k*dt
*

*	Lower & Upper BC's
*******************************************************


	amewa=(ia-((Nx+1)/2 +1))*dx
	f=Dsin(2*pi*freq*t)
	
	
	do 87 i = 1,Nx+2
*	Lower 
		j=1
		d6=x(i)**2+y(1)**2
		Psi(i,1) = 0.0d0
		
		if ((i.ge.ia) .and. (i.le.ib)) then

	Psi(i,1)=(-c0*(0.0-.5*amewa*Dsqrt(amewa**2+1.0d0)-.5*Dsinh(amewa)
     c  +.5*x(i)*Dsqrt(x(i)**2+1)+.5*Dsinh(x(i))))*f

		else
			goto 910
		endif		
910		continue
		
		if (i.gt.ib) then
			Psi(i,1)=Psi(ib,1)
		else
			goto 920
		endif					
920		continue		
		
		
		Omega(i,1) = (7.0d0*Psi(i,1)-8.0*Psi(i,2)+Psi(i,3))/(2.0*dyy)/d6
		u(i,1) = 0.0
		v(i,1) = 0.0
		
		if ((i.gt.ia) .and. (i.lt.ib)) then
			v(i,1)=c0*f
		else 
			go to 930
		endif		
930		continue
			
		if ((i.ge.(ia-1)) .and. (i.le.(ib+1))) then
		
	Omega(i,1) = Omega(i,1)
     c  +(v(i+1,1)*Dsqrt(x(i+1)**2+1)-v(i-1,1)*Dsqrt(x(i-1)**2+1))
     c  /(2.0*dx)/d6

		
		else
			go to 940
		endif
		
940		continue		

*	Upper 
		j=My+2
		Omega(i,My+2) = 0.0
		Psi(i,My+2) = (x(i)+A)*(y(j)-1)
		u(i,My+2) = (x(i)+A)/Dsqrt(x(i)**2+y(j)**2)
		v(i,My+2) =-(y(j)-1)/Dsqrt(x(i)**2+y(j)**2)

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
		u(i,j) = (x(i)+A)/Dsqrt(x(i)**2+y(j)**2)
		v(i,j) =-(y(j)-1)/Dsqrt(x(i)**2+y(j)**2)
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
		u(i,j) = (x(i)+A)/Dsqrt(x(i)**2+y(j)**2)
		v(i,j) =-(y(j)-1)/Dsqrt(x(i)**2+y(j)**2)
89	continue

*	Calculate velocities
********************************************************
	do 100 i = 2,Nx+1
	do 110 j = 2,My+1
		call UCalc(Nx,My,i,j,dy2,Psi,u,x,y)
		call VCalc(Nx,My,i,j,dx2,Psi,v,x,y)
110	continue
100	continue

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
	
* NOTE:	following line clears the mess currently stored in the string

		string = '                                             '
		
		call btd(ct,string,3,rcc)
		
		do 19 yy=2,4
			pfile(LL+yy:)=string(yy-1:)
19 			continue
		print *,'Writing Incrimental File=',pfile
		outfile=pfile

		flag=1
		print *,''		
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
	outfile=filename
	print *,'Writing output file'
	
	
*	Write results to file
***********************************************************************    

4000	open(unit=2,file=outfile,status='new')

	write (2,form) Omega
	write (2,form) Psi
	write (2,form) u
	write (2,form) v

	do 155 i=1,Nx+2
	OutIn(i)=0
	OutDP(i)=0
155	CONTINUE

	OutIN(1) = Nx+2
	OutIN(2) = My+2
	OutIN(3) = k
	OutIN(4) = Kerr
	OutIN(5) = IBL

*NOTE: changes Jet beginning and end point to be saved as integers.

	OutIN(6) = ia
	OutIN(7) = ib
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
	OutDP(13) = freq
	OutDP(14) = c0


	write (2,form1) OutIN
	write (2,form) OutDP
	
	call FDATE(STRNG)
	print*,STRNG
	
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
	double precision Omega(Nx+2,My+2),U(Nx+2,My+2),v(Nx+2,My+2),
     c	x(Nx+2),y(My+2)
	double precision Omega0(Nx+2,My+2),Cx2,Cy2,alpha,alphaX,alphaY,
     c	d,d2,d3,dip1,dim1,djp1,djm1


	do 80 i = 2,Nx+1
		do 90 j = 2,My+1
	d=Dsqrt(x(i)**2+y(j)**2)
	dip1=Dsqrt(x(i+1)**2+y(j)**2)
	dim1=Dsqrt(x(i-1)**2+y(j)**2)
	djp1=Dsqrt(x(i)**2+y(j+1)**2)
	djm1=Dsqrt(x(i)**2+y(j-1)**2)
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
	d2=Dsqrt(x(i)**2+y(j)**2)
	Omega(i,j) = Omega0(i,j)*(1+3.0*Cx2*u(i,j)/d2+2.0*alphaX/
     c	(d2*d2)-2.0*alphaY/(d2*d2)) + 
     c	Omega0(i+1,j)*(-4.0*Cx2*u(i,j)/d2 - 5.0*alphaX/(d2*d2)) +
     c	Omega0(i+2,j)*(Cx2*u(i,j)/d2 + 4.0*alphaX/(d2*d2)) +
     c	Omega0(i+3,j)*(-alphaX/(d2*d2)) +
     c	Omega0(i,j+1)*(-Cy2*v(i,j)/d2 + alphaY/(d2*d2)) +
     c	Omega0(i,j-1)*(Cy2*v(i,j)/d2 + alphaY/(d2*d2))
	 
60	enddo

	i=Nx+2
	do 70 j = My+1,2,-1
	d3=Dsqrt(x(i)**2+y(j)**2)
	Omega(i,j) = Omega0(i,j)*(1-3.0*Cx2*u(i,j)/d3+2.0*alphaX/
     c	(d3*d3)-2.0*alphaY/(d3*d3)) + 
     c	Omega0(i-1,j)*(4.0*Cx2*u(i,j)/d3 - 5.0*alphaX/(d3*d3)) +
     c	Omega0(i-2,j)*(-Cx2*u(i,j)/d3 + 4.0*alphaX/(d3*d3)) +
     c	Omega0(i-3,j)*(-alphaX/(d3*d3)) +
     c	Omega0(i,j+1)*(-Cy2*v(i,j)/d3 + alphaY/(d3*d3)) +
     c	Omega0(i,j-1)*(Cy2*v(i,j)/d3 + alphaY/(d3*d3))
	
70	enddo
75	continue
	end


*************************************************************************
*	Iterative Stream Function Routine
*************************************************************************
*	Tol 		Tolerance level to cease iterations
*	PsiTol		Calculated difference in Psi values per iteration
*	Psi1m,Psi0m	Maximum single value of the previous Psi matrix
*	Psi0		Previous Psi Matrix Values
*************************************************************************
	subroutine PsiCalc(Nx,My,Kappa2,KappaA,dxx,Psi,Omega,kPsi,x,y,Tol)
	integer Nx,My,Ot,i,j,kPsi,n,m
	double precision Omega(Nx+2,My+2),Psi(Nx+2,My+2),Psi0(Nx+2,My+2),
     c	Kappa2,KappaA,dxx,Tol,PsiTol,Psi1m,Psi0m,x(Nx+2),y(My+2),d2

	PsiTol = 1
	Psi1m = 0
	kPsi = 0

10	if (abs(PsiTol) .gt. Tol) then
	kPsi = kPsi+1
	Psi0m = Psi1m
	Psi1m = 0

	do 20 n = 1,Nx+2
	do 30 m = 1,My+2
		Psi0(n,m) = Psi(n,m)
30	enddo
20	enddo


	do 40 i = 2,Nx+1
	do 50 j = My+1,2,-1
	d2=x(i)**2+y(j)**2
	Psi(i,j) = KappaA*(dxx*Omega(i,j)*d2 + Psi0(i+1,j) +
     c	Psi0(i-1,j) + Kappa2*(Psi0(i,j+1) + Psi0(i,j-1)))

	
	if (abs(Psi(i,j)-Psi0(i,j)) .gt. Psi1m) then
		Psi1m = abs(Psi(i,j)-Psi0(i,j))
	endif

50	enddo
40	enddo
	
90	PsiTol = abs(Psi1m)
	goto 10

	endif

	end

***********************************************************************
*	Finite difference approximation for u velocities
***********************************************************************
	subroutine UCalc(Nx,My,i,j,dy2,Psi,u,x,y)
	integer Nx,My,i,j
	double precision Psi(Nx+2,My+2),u(Nx+2,My+2),x(Nx+2),y(My+2)
	double precision dy2,d3
	d3=Dsqrt(x(i)**2+y(j)**2)
	u(i,j) = (Psi(i,j+1)-Psi(i,j-1))/(dy2)/d3

	end

***********************************************************************
*	Finite difference approximation for v velocities
***********************************************************************
	subroutine VCalc(Nx,My,i,j,dx2,Psi,v,x,y)
	integer Nx,My,i,j
	double precision Psi(Nx+2,My+2),v(Nx+2,My+2),x(Nx+2),y(My+2)
	double precision dx2,d4
	d4=Dsqrt(x(i)**2+y(j)**2)
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
	character(len=1) string(LS)
	integer rc,ct
	
	character(len=10), save ::
     c	numerl(10) = ['0','1','2','3','4','5','6','7','8','9']
	
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
	
	


