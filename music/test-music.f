*       This version has been written by V.A.Kudryavtsev, University of Sheffield,
*       v.kudryavtsev@sheffield.ac.uk, April 2008.
*
*	This is a simple test file which 
*       1. Calculates muon cross-sections and energy losses (CALL MUCRSEC and CALL MULOS), 
*       and write them on the disk if necessary (if INIT = 1);
*       2. Read the cross-sections and energy losses from the disk (CALL INITIALIZE_MUSIC);
*       3. Propagate muons down to a predefined depth/distance in a material 
*       (CALL MUON_TRANSPORT). 
*
*       It uses the following files:
*       - music-crosssections.f - file with all subroutines to calculate muon interaction
*       cross-sections and energy losses;
*	- music.f - file with all subroutines for muon transport including 
*		initialize_music, muon_transport, music etc.;
*	- music-cross-sections.dat - formatted file with tables with 
*		integral cross-sections; the file is created if the parameter INIT = 1;
*               once it is created the parameter INIT can be set to 0.
*	- music-eloss.dat - formatted file with continuous muon energy
*		losses (v<v_cut=10^(-3)); the file is created if the parameter INIT = 1;
*               once it is created the parameter INIT can be set to 0.
*	- music-double-diff-rock.dat - formatted file with inelastic
*		scattering cross-section to sample angular deviation; this cross-section
*               was calculated by me for rock and supplied with the code; the cross-section
*               is very similar for all types of rock; if you need it for a different material,
*               please, contact me.
*               
*	Initial parameters (all real*8):
*	x,y,z - initial muon coordinates in cm;
*	cx,cy,cz - initial direction cosines;
*	emu - muon energy in GeV;
*	depth - distance to the next output from the code (depth or distance to propagate
*	muons) in cm;
*	ttime - time in nsec.
*       idim and idim1 - parameters that determine 1D (both equal to 0) or 3D (both equal to 1)
*       simulations; if you are interested in angular deviation and lateral displacement
*       of muons you should set both of them to 1.
*	Final parameters are the same as initial parameters but contain
*	the values after muon was propagated down to 'depth'. The value of
*	depth is the same as before.
*	Note: the subroutines use standard CERN library routines and 
*	random number generators RANLUX. The code uses also
*	the generator of 2 numbers sampled according to Gaussians with a
*	predefined correlation coefficient - CORSET and CORGEN
*	from CERN library. The sequence of random numbers should be
*	initialized at the beginning of each run using RLUXGO and
*	RMARIN (see CERN library documentation)
*
*	As an example, this program propagates 10000 muons to a distance of
*	1000 m (100000 cm) in a predefined rock.
*	Muon energy is sampled according to the power-law spectrum at surface
*	with index 2.7 starting from 900 GeV.
*
*! For any problems please contact:
*! 
*! v.kudryavtsev@sheffield.ac.uk
*!
*! References:
*! 1. P. Antonioli, C. Ghetti, E.V. Korolkova, V.A. Kudryavtsev, G. Sartorelli,
*! Astroparticle Physics, 7 (1997) 357.
*! 2. V.A. Kudryavtsev, E. V. Korolkova, N. J. C. Spooner,
*! Physics Letters B, 471 (1999) 251.
*! 3. V.A.Kudryavtsev, Computer Physics Communications, 180 (2009) 339; 
*! http://dx.doi.org/10.1016/j.cpc.2008.10.013; arXiv:0810.4635 [physics.comp-ph]
*! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
*
*       The definitions below should be present in your 'main' program.
	implicit real*8 (a-h,o-z)
	parameter (pi=3.141592654)
	real*4 yfl
	real*4 zz0(20),a0(20),fr0(20)     !max number of elements for any material - 20
	real*4 par_ion(6)

*       Define the rock or other material here; max number of elements = 20

*       Example for standard rock
	data zz0/11.,19*0./    ! standard rock
	data a0/22.,19*0./     ! standard rock
	data fr0/1.,19*0./     ! standard rock
	data par_ion/136.4,-3.774,0.083,3.412,3.055,0.049/ !density correction parameters (st. rock)

*       Example for another type of rock
c	data zz0/1.,6.,7.,8.,11.,
c	1    12.,13.,19.,20.,26.,10*0./  !array of atomic numbers for all elements
c	data a0/1.00794,12.0107,14.067,15.9994,22.98977,
c	1    24.3050,26.98154,39.0983,40.078,55.845,10*0./ !array of atomic weights
c	data fr0/0.01,0.05,0.15,0.40,0.05,
c	1    0.06,0.06,0.01,0.20,0.01,10*0./ !fraction by mass for all elements
c	data par_ion/136.4,-3.774,0.083,3.412,3.055,0.049/ !density correction parameters (st. rock)

*       Example for pure water
c	data zz0/1.,8.,18*0./
c	data a0/1.008,15.999,18*0./
c	data fr0/0.1119,0.8881,18*0./
c	data par_ion/75.0,-3.502,0.2065,3.007,2.5,0.24/ !density correction parameters (water)
c	data par_ion/75.0,-3.5017,0.09116,3.4773,2.8004,0.24/ !density correction parameters 
c       (pure water - Sternheimer)

*       If you need density correction parameters (for ionisation energy loss) for other materials,
*       please, contact the code developer

	rad = 26.48d0     ! radiation length (for standard rock)
	rho = 2.65d0      ! density (for standard rock)
c	rad = 36.08d0     ! radiation length (for water)
c	rho = 1.00d0      ! density (for water)
	idim = 1        ! to switch off mult.scatt. set idim=0
	idim1 = 1       ! to switch off scatt. due to other processes set idim1=0
	
	minv = -30 ! corresponds to the cut between stichastic and continuous 
	         ! energy losses of 10^(-3), do not change this without expert advice
	init = 1 !if cross-sections have not been calculated yet set init = 1, otherwise =0

*       Calculate muon interaction cross-sections and continuous energy losses
*       and write them on the disk
	if(init.eq.1) then
	   call mucrsec(minv,zz0,a0,fr0)
	   call mulos(minv,zz0,a0,fr0,par_ion)
	end if

*       Read cross-sections and muon energy losses from the disk; do not change
	call initialize_music(minv,rho,rad)
c	print *,minv,rho,rad

	nmumax=1000  ! number of muons to transport

*       Initial muon parameters
	theta0=0.   ! zenith angle in radians
	phi0=0.     ! azimutal angle in radians
	x0=0.       ! x-ccordinate
	y0=0.       ! y-ccordinate
	z0=0.       ! z-ccordinate; muons are transported along z-axis
	
*       fixed muon energy or lowest energy if sampled according to the power-law spectrum
	emu00=900.

	depth=300000./rho  ! depth or distance to transport muons
	ttime0=0.          ! initial time
	iranlux=1          ! initial seed for random number generator RLUXGO
	iranlux1=1         ! initial seed for random number generator RMARIN
	call rluxgo(3,iranlux,0,0)  ! initialization of random number generator RLUXGO
	call rmarin(iranlux1,0,0)   ! initialization of random number generator RMARIN

*       Initialization of parameters to calculate mean muon energy, coordinates etc.
*       Needed only for test runs
	kmu=0
	eemu=0.
	rmu=0.
	xmu=0.
	ymu=0.
	zmu=0.
	themu=0.
	phimu=0.

*       Start muon transport (the number of muons to be transported is NMUMAX)
	do i=1,nmumax
	call ranlux(yfl,1)

*       sample muons according to the power-law spectrum; a user should 
*       change this to the appropriate formula or use fixed muon energy
	emu0=(emu00**(-2.7)*yfl)**(-1./2.7)

*       uncomment this if you want to set energy to emu00
c	emu0=emu00
	x=x0
	y=y0
	z=z0
	theta=theta0
	phi=phi0
	cz0=dcos(theta0)
	cx0=dsin(theta0)*dcos(phi0)
	cy0=dsin(theta0)*dsin(phi0)
	cx=cx0
	cy=cy0
	cz=cz0
	emu=emu0
	ttime=ttime0

*       Call muon transport subroutine
	call muon_transport(x,y,z,cx,cy,cz,emu,depth,ttime,idim,idim1)

*       Convert direction cosines to theta and phi (optional)
	theta=dacos(cz)
	IF(cx.NE.0.) THEN
	   phi=datan(cy/cx)
	   IF(cx.LT.0.) phi=pi+phi
	END IF
	IF(cx.EQ.0.) THEN
	   phi=pi/2.
	   IF(cy.LT.0.) phi=pi*3./2.
	END IF
	IF(phi.GT.pi*2.) phi=-pi*2.+phi
	IF(phi.LT.0.) phi=phi+pi*2.

*	Add final muon parameters to the sum of all previous muons but only if
*       muon energy exceeds muon rest mass; this is needed only if you want
*       to calculate mean muon parameters after transport
	if(emu.gt.0.106) then
	   kmu=kmu+1
	   eemu=eemu+emu
	   xmu=xmu+x
	   ymu=ymu+y
	   zmu=zmu+z
	   themu=themu+theta
	   phimu=phimu+phi
	   rmu=sqrt(x*x+y*y)+rmu
	end if

c         print *,'theta0,phi0  ',theta0,phi0
c         print *,'theta,phi,  ',theta,phi
c         print *,'x0,y0,z0  ',x0,y0,z0
c         print *,'x,y,z  ',x,y,z
c         print *,'emu0  ',emu0
c         print *,'emu  ',emu		!if emu<0.106 muon was stopped

	end do

*       Calculate mean muon parameters after transport of NMUMAX muons
	if(kmu.gt.0) then
	   spmu=1.*kmu/nmumax
	   eemu=eemu/kmu
	   xmu=xmu/kmu
	   ymu=ymu/kmu
	   zmu=zmu/kmu
	   themu=themu/kmu
	   phimu=phimu/kmu
	   rmu=rmu/kmu
	end if
	print *,'surv. prob.',spmu
	print *,'mean energy',eemu
	print *,'x,y,z',xmu,ymu,zmu
	print *,'lateral displ.',rmu
	print *,'theta',themu
	print *,'phi',phimu
	
	stop
	end

