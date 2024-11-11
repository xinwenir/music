C	CALCULATION OF MUON INTERACTION CROSS-SECTIONS
	subroutine mucrsec(minv,z0,a0,fr0)
	implicit real*4 (a-h,o-z)
	real*4 CS(100,71,4)
	DATA CS/28400*0./
	real*4 z0(20),a0(20),fr0(20)     !max number of elements for any material - 20
	real*4 fr1(20)
	data fr1/20*0./

	s = 0.
	n = 0
	do i = 1,20
	   if(z0(i).gt.0.5) then
	      s = s + fr0(i)
	      n = n + 1
	   end if
	end do
	do i = 1,20
	   if(z0(i).gt.0.5) then
	      fr1(i) = fr0(i) / s
	   end if
	end do

c	print *,z0
c	print *,a0
c	print *,fr1
	   
	DO 1 K=1,71
	EM2=0.0+0.1*(K-1)
	DV=0.01
	SB=0.
	SP=0.
	SN=0.
	SE=0.
	EM=10.**EM2
	V2=0.
	I=0
	I1=1
2	V1=V2-DV/2.
	V=10.**V1
	E=V*EM
	I=I+1


	SB1=0.
	SP1=0.
	SN1=0.
	SE1=0.

	do j = 1, 20
	   if(z0(j).gt.0.5) then

	fb=0
	fp=0.
	fn=0.
	fe=0.
	z = z0(j)
	a = a0(j)
	fr = fr1(j)

	CALL SIGBK(E,EM,Z,A,FB)
cc        CALL SIGB(E,EM,Z,A,FB)
	CALL SIGP(E,EM,Z,A,FP)
cc	CALL SIGNU(E,EM,A,FN)
	CALL SIGNU_BS(E,EM,A,FN)
cc	CALL SIGNU_F2(E,EM,A,FN)
	CALL SIGE(E,EM,Z,A,FE)

	SB1=SB1+FB*fr
	SP1=SP1+FP*fr
	SN1=SN1+FN*fr
	SE1=SE1+FE*fr

	end if
	end do

	sb = sb + sb1*DV*2.303
	sp = sp + sp1*DV*2.303
	sn = sn + sn1*DV*2.303
	se = se + se1*DV*2.303

	V2=V2-DV
	IF(I.EQ.5*I1) THEN
	I1=I1+1
	CS(I1,K,1)=SB
	CS(I1,K,2)=SP
	CS(I1,K,3)=SN
	CS(I1,K,4)=SE
	END IF
	IF(V2.GE.-4.94999) GO TO 2
c	write(6,3)EM2,SB,SP,SN,SE
  3	FORMAT(1X,'LOG(E)=',F4.1,1X,'B=',E10.3,
     *  1X,'P=',E10.3,1X,'N=',E10.3,1X,'E=',E10.3)
  1	CONTINUE
	OPEN(UNIT=3,NAME=
     *  'music-cross-sections.dat',
     *  FORM='FORMATTED',STATUS='unknown')
	WRITE(3,4)CS
4	FORMAT(4(71(20(5E14.6/)/)/)/)
	CLOSE(3)
	return
	END

C	CALCULATION OF ENERGY LOSSES OF MUONS
C
	subroutine mulos(minv,z0,a0,fr0,par_ion)	
	implicit real*4 (a-h,o-z)
	real*4 z0(20),a0(20),fr0(20)     !max number of elements for any material - 20
	real*4 fr1(20)
	real*4 par_ion(6)
	DIMENSION EMULO(81,6)
	DATA EMULO/486*0./,fr1/20*0./

	ame=0.000511
	ammu=0.105658

	ai = par_ion(1)
	cc = par_ion(2)
	aa = par_ion(3)
	bb = par_ion(4)
	y1 = par_ion(5)
	y0 = par_ion(6)

	s = 0.
	n = 0
	s1 = 0.
	s2 = 0.
	do i = 1,20
	   if(z0(i).gt.0.5) then
	      s = s + fr0(i)
	      n = n + 1
	   end if
	end do
	do i = 1,20
	   if(z0(i).gt.0.5) then
	      fr1(i) = fr0(i) / s
	   end if
	end do
	do i = 1,20
	   if(z0(i).gt.0.5) then
	      s1 = s1 + z0(i)/a0(i)*fr1(i)
	      s2 = s2 + z0(i)*z0(i)/a0(i)*fr1(i)
	   end if
	end do
c	print *,s1,s2,fr1
	zmean = s2/s1
	amean = zmean/s1
	
	DO 1 K=1,81	!81
	EM2=-1.0+0.1*(K-1)
	DV=0.01
	SB=0.
	SP=0.
	SN=0.
	se=0.
	b0=0.
	EM=10.**EM2

	gam=em/ammu
	if(gam.gt.1.) then
	beta=sqrt(gam*gam-1.)/gam
	ems=2.*beta*beta*gam*gam*ame/
     *	(1.+2.*gam*ame/ammu+(ame/ammu)**2)
        A1=0.1535*zmean/amean*1.E-3/beta/beta
	a2=alog(ame*ammu*1.e18/(ai*ai))
	B11=A2+0.693-2.*beta*beta+0.25*(EMS/EM)**2-cc+
     *	ALOG(EMS/ammu)
cc	e_cut=0.01
	v_cut=10.**(0.1*minv)	!0.001
	e_cut=em*v_cut
cc	e_cut=ems
	b12=a2+0.693-beta*beta*(1.+e_cut/ems)-cc+
	1    alog(e_cut/ammu)+0.25*(e_cut/EM)**2
	b1=b12		!b12
	if(em.lt.1.) b1=b11
	IF(ALOG10(beta*gam).GE.y1) b0=a1*b1
        IF(ALOG10(beta*gam).lt.y1.and.ALOG10(beta*gam).ge.y0) 
     *	B0=A1*(B1-aa*(y1-ALOG10(beta*gam))**bb)
        IF(ALOG10(beta*gam).lt.y0)
     *	b0=a1*(b1+4.6052*alog10(beta*gam))
	end if
cc	print *,'gam,beta,ems,a1,a2,b1,b0',gam,beta,ems,a1,a2,b1,b0
	EMULO(K,6)=B0

	V2=0.1*minv		!-3.0
	if(em2.ge.0) then
	I=0
2	V1=V2-DV/2.
	V=10.**V1
	E=V*EM
	I=I+1

	SB1=0.
	SP1=0.
	SN1=0.
	se1=0.

	do i = 1, 20
	   if(z0(i).gt.0.5) then

	fb=0.
	fp=0.
	fn=0.
	fe=0.
	z = z0(i)
	a = a0(i)
	fr = fr1(i)

cc	CALL SIGB(E,EM,Z,A,FB)
        CALL SIGBK(E,EM,Z,A,FB)
        CALL SIGP(E,EM,Z,A,FP)
cc	CALL SIGNU(E,EM,A,FN)
	CALL SIGNU_BS(E,EM,A,FN)
cc	CALL SIGNU_f2(E,EM,A,FN)
	call sige1(e,em,z,a,fe)

	SB1=SB1+FB*fr
	SP1=SP1+FP*fr
	SN1=SN1+FN*fr
	SE1=SE1+FE*fr

	end if
	end do

	sb = sb + sb1*DV*2.303*V
	sp = sp + sp1*DV*2.303*V
	sn = sn + sn1*DV*2.303*V
	se = se + se1*DV*2.303*V

	V2=V2-DV
	IF(V2.GE.-9.0) GO TO 2
	end if
	EMULO(K,1)=SB
	EMULO(K,2)=SP
	EMULO(K,3)=SN
	emulo(k,4)=se
	emulo(k,5)=sb+sn+sp+se
***	EMULO(K,5)=SB+SN+SP
C	WRITE(6,11)EM,EMULO(K,4)
11	FORMAT(2X,'ENERGY',F12.2,2X,F14.5)
c	fac=10./6.023*A*em
c	WRITE(6,3)EM2,SB,SP,SN,EMULO(K,5)*em,B0
3	FORMAT(1X,F4.1,1X,'B=',E10.3,
     *  1X,'P=',E10.3,1X,'N=',E10.3,1X,
     *	E11.4,1X,'I=',E11.4)
1	CONTINUE
	OPEN(UNIT=3,NAME=
     *  'music-eloss.dat',
	2    FORM='FORMATTED',STATUS='unknown')
c        WRITE(6,3)EM2,SB,SP,SN,se,EMULO(71,5),B0
	WRITE(3,4)EMULO
4	FORMAT(16(5E12.5/),E12.5/)
	WRITE(3,5)zmean,amean
 5	format(2f12.3)
	CLOSE(3)
	return
	END


C	CALCULATION OF THE MUON BREMSSTRAHLUNG CROSS-SECTION
c	using the formula by Kelner, Kokoulin, Petrukhin (1995)
C	INPUT PARAMETERS: ENERGY OF PHOTON E, MUON ENERGY EM,
C	NUCLEAR CHARGE Z, NUCLEAR MASS A
C	OUTPUT PARAMETER: VALUE OF V*DSIGMA/DV, WHERE V=E/EM
C
	SUBROUTINE SIGBK(E,EM,Z,A,FB)	
	V=E/EM
	fmmu=0.105658
	fme=0.000511
	IF(V.GE.1.0) GO TO 1
	VMAX=1.-3.*fmmu/4./EM*sqrt(2.718)*Z**(1./3.)
	IF(V.GE.VMAX) GO TO 1
	IF(EM.LE.fmmu) GO TO 1
	B1=182.7
	B2=1440.
	IF(Z.LT.1.5) THEN
	   B1=202.4
	   B2=446.
	END IF
	del=fmmu**2*V/(2.*EM*(1.-V))
	dn=1.54*A**0.27
	del_el=alog(dn/(1.+del*(dn*sqrt(2.718)-2.)/fmmu))
	FI=alog(183.*fmmu/fme*Z**(-1./3.)/(1.+del*sqrt(2.718)*183.*
     *	Z**(-1./3.)/fme))-del_el
	FI_in=alog(fmmu/del/(del*fmmu/fme/fme+sqrt(2.718)))-
     *	alog(1.+1./(del*sqrt(2.718)*1429.*Z**(-2./3)/fme))+del_el
	sig0=(2.*fme/fmmu)**2*2.81794e-13*6.023e23*2.81794e-13/137.*
	1    (4./3.-4./3.*v+v*v)
	fb=sig0*(Z*Z/A*FI+Z/A*FI_in)
	IF(FB.LE.0.) FB=0.
	GO TO 2
1	FB=0.
2	CONTINUE
	RETURN
	END


C	CALCULATION OF THE MUON BREMSSTRAHLUNG CROSS-SECTION
C	INPUT PARAMETERS: ENERGY OF PHOTON E, MUON ENERGY EM,
C	NUCLEAR CHARGE Z, NUCLEAR MASS A
C	OUTPUT PARAMETER: VALUE OF V*DSIGMA/DV, WHERE V=E/EM
C
	SUBROUTINE SIGB(E,EM,Z,A,FB)	
	V=E/EM
	IF(V.GE.1.0) GO TO 1
	fmmu=0.105658
	fme=0.000511
	VMAX=1.-3.*fmmu/4./EM*sqrt(2.718)*Z**(1./3.)
	IF(V.GE.VMAX) GO TO 1
	IF(EM.LE.fmmu) GO TO 1
	B1=184.15
	B2=1194.
	IF(Z.EQ.1) THEN
	   B1=202.4
	   B2=446.
	END IF
	T1=B1/sqrt(2.718)/fme*Z**(-1./3.)
	T2=B2/sqrt(2.718)/fme*Z**(-2./3.)
	T3=0.105658
	T4=T3**2*V/(2.*EM*(1.-V))
	T5=1.9*T3*Z**(-1./3.)
	T6=(1.+4.*T3**2/T5**2)**0.5
	DEL1=ALOG(T3/T5)+T6/2.*ALOG((T6+1.)/(T6-1.))
	DEL2=ALOG(T3/T5)+0.25*(3.*T6-T6**3)*ALOG((T6+1.)/(T6-1.))+
     *  2.*T3**2/T5**2
	FI01=0.5*(ALOG(T3**2*T1**2/(1.+T4**2*T1**2))+1.)-T4*T1*
     * 	ATAN(1./(T4*T1))+0.5/Z*(ALOG(T3**2*T2**2/
     *  (1.+T4**2*T2**2))+1.)-T4*T2/Z*ATAN(1./(T4*T2))
	FI02=0.5*(ALOG(T3**2*T1**2/(1.+T4**2*T1**2))+2./3.)+
     *  2.*T4**2*T1**2*(1.-T4*T1*ATAN(1./(T4*T1))-0.75*ALOG(1.+
     *  1./(T4**2*T1**2)))+0.5/Z*(ALOG(T3**2*T2**2/(1.+T4**2*T2**2))+
     *  2./3.)+2./Z*T4**2*T2**2*(1.-T4*T2*ATAN(1./(T4*T2))-
     *  0.75*ALOG(1.+1./(T4**2*T2**2)))
	FI1=FI01-DEL1
	FI2=FI02-DEL2
	FB=(2.*fme/fmmu)**2*2.81794e-13*6.023e23*2.81794e-13/137.*
	1    Z**2/A*((2.-2.*V+V*V)*FI1-2./3.*(1.-V)*FI2)
	IF(FB.LE.0.) FB=0.
	GO TO 2
1	FB=0.
2	CONTINUE
	RETURN
	END


C	CALCULATION OF THE PAIR PRODUCTION CROSS-SECTION OF MUONS
C	INPUT PARAMETERS: ENERGY OF PAIR E, MUON ENERGY EM,
C	NUCLEAR CHARGE Z, NUCLEAR MASS A
C	OUTPUT PARAMETER: VALUE OF V*DSIGMA/DV, WHERE V=E/EM
C
	SUBROUTINE SIGP(E,EM,Z,A,FP)	
	REAL*8 B,S,YE,YM,AE,AM,FE,FM,csi2(8),em3(8),csi3(8)
**	data em3/0.,1.,2.,3.,4.,5.,6.,7./
**	data csi2/0.,0.3,0.8,1.1,1.23,1.28,1.3,1.3/
	data em3/0.,1.,1.5,2.,2.5,3.,3.5,7./
        data csi2/0.,0.60,0.82,0.91,0.97,1.03,1.07,1.07/	! rock
        data csi3/0.,0.62,0.85,0.94,1.00,1.06,1.10,1.10/	! hydrogen
	V=E/EM
	PE=0.000511
	PM=0.105655
	VMIN=4.*PE/EM
	VMAX=1.-3.*SQRT(2.718)*PM/4./EM*Z**(1./3.)
	IF(V.GE.VMAX.OR.V.LE.VMIN) GO TO 1
	AL=1./137.036
	AN=6.023E23
	AR=3.8616E-11
	PI=3.1416
	RMAX=SQRT(1.-4.*PE/EM/V)*(1.-6.*PM**2/EM**2/(1.-V))
	RMIN=0.
	B=V*V/2./(1.-V)
	S2=0.
	X1=RMIN
	DX=0.01
3	X=X1+DX/2.
	S=(PM*V/2./PE)**2*(1.-X*X)/(1.-V)
	YM=(4.+X*X+3.*B*(1.+X*X))/((1.+X*X)*(3./2.+2.*B)*DLOG(3.+S)+
     *  1.-3./2.*X*X)
	YE=(5.-X*X+4.*B*(1.+X*X))/(2.*(1.+3.*B)*DLOG(3.+1./S)-
     *  X*X-2.*B*(2.-X*X))
	AM=DLOG(2./3.*PM/PE*183.*Z**(-2./3.)/(1.+2.*PE*SQRT(2.718)*
     *  183.*Z**(-1./3.)*(1.+S)*(1.+YM)/EM/V/(1.-X*X)))
	AE=DLOG(183.*Z**(-1./3.)*DSQRT((1.+S)*(1.+YE))/(1.+2.*PE*
     *  SQRT(2.718)*183.*Z**(-1./3.)*(1.+S)*(1.+YE)/EM/V/(1.-X*X)))-
     *  0.5*DLOG(1.+(1.5*PE/PM*Z**(1./3.))**2*(1.+S)*(1.+YE))
	FM=(((1.+X*X)*(1.+1.5*B)-1./S*(1.+2.*B)*(1.-X*X))*DLOG(1.+S)+
     *  S*(1.-X*X-B)/(1.+S)+(1.+2.*B)*(1.-X*X))*AM
	FE=(((2.+X*X)*(1.+B)+S*(3.+X*X))*DLOG(1.+1./S)+
     *  (1.-X*X-B)/(1.+S)-(3.+X*X))*AE
	S1=FE+(PE/PM)**2*FM
	S2=S2+S1*DX
	X1=X1+DX
	IF(X1.GE.RMAX-0.06) DX=0.01
	IF(X1.GE.RMAX-0.011) DX=0.001
	IF(X1.GE.RMAX-0.0011) DX=0.0001
	IF(X1.GE.RMAX-0.00011) DX=0.00001
	IF(X1.LT.RMAX-DX/2.) GO TO 3
	i=1
10	i=i+1
	if(i.gt.8) go to 111
	if(alog10(em).gt.em3(i)) go to 10
	go to 112
 111	i=8
 112	continue
	csi=(csi2(i)-csi2(i-1))/(em3(i)-em3(i-1))*(alog10(em)-
     *	em3(i-1))+csi2(i-1)
	if(z.lt.1.5) 
     *	csi=(csi3(i)-csi3(i-1))/(em3(i)-em3(i-1))*(alog10(em)-
     *	em3(i-1))+csi3(i-1)
	FP=2.*AL**4*2./3./PI*AR*AN*AR*Z*(Z+csi)/A*(1.-V)*S2
	IF(FP.LE.0.) FP=0.
	GO TO 2
1	FP=0.
2	CONTINUE
	RETURN
	END


C	CALCULATION OF MUON INELASTIC SCATTERING CROSS-SECTION
C	INPUT PARAMETERS: MUON ENERGY LOSS E, MUON ENERGY EM,
C	ATOMIC MASS A
C	OUTPUT PARAMETER: VALUE OF V*DSIGMA/DV, WHERE V=E/EM
C
	SUBROUTINE SIGNU(aE,aEM,aA,aFN)	
	implicit double precision (b-z)
	e=ae
	em=aem
	a=aa
	V=E/EM
        prm=0.938
        pm=0.105658*0.105658
        if(em.le.dsqrt(pm)) go to 1
        vmin=0.14/em
c	vmin = 0.2/em
        vmax=(2.*prm*em*em-2.*prm*pm)/(2.*prm*em+pm+prm*prm)/em
        if(v.lt.vmin.or.v.gt.vmax) go to 1
        em1=em-e
        if(em1.le.dsqrt(pm)) go to 1
	SIGMA=114.3+1.647*(dLOG(0.0213*E))**2
	IF(E.LE.17.) SIGMA=96.1+82./dSQRT(E)
	X5=0.00282*SIGMA*A**(1./3.)
	GX=3./X5**3*(X5**2/2.-1.+EXP(-X5)*(1.+X5))
	if(a.le.1.5) gx=1.
	T1=PM*V**2/(1.-V)
	T2=1.-2./V+2./V/V
	T3=0.54
	T4=1.8
	Y1=7.E-10*SIGMA*V*V
	Y2=T2*dLOG(1.+T3/T1)-T2*T3/(T3+T1)-2.*PM/T1
	Y3=T2*dLOG(1.+T4/T1)-2.*PM/T1
	Y4=T3/(T3+T1)
	Y5=T4/T1*dLOG(1.+T1/T4)
	FN=Y1*(0.75*GX*Y2+0.25*Y3+0.5*PM/T1*(0.75*GX*Y4+0.25*Y5))
	IF(FN.LE.0.) FN=0.
	GO TO 2
1	FN=0.
2	CONTINUE
	afn=fn
	RETURN
	END


C	CALCULATION OF MUON INELASTIC SCATTERING CROSS-SECTION
C	INPUT PARAMETERS: MUON ENERGY LOSS E, MUON ENERGY EM,
C	ATOMIC MASS A
C	OUTPUT PARAMETER: VALUE OF V*DSIGMA/DV, WHERE V=E/EM
C
	SUBROUTINE SIGNU_BS(aE,aEM,aA,aFN)	
	implicit double precision (b-z)
	real*8 atau(8,7),emul(7)
	data emul/3.d0,4.d0,5.d0,6.d0,7.d0,8.d0,9.d0/
	data atau/7.174409d-4,-0.2436045d0,-0.2942209d0,-0.1658391d0,
	1    -0.05227727d0,-9.328318d-3,-8.751909d-4,-3.343145d-5,
	2    1.7132d-3,-0.5756682d0,-0.68615d0,-0.3825223d0,
	3    -0.1196482d0,-0.02124577d0,-1.987841d-3,-7.584046d-5,
	4    4.082304d-3,-1.553973d0,-2.004218d0,-1.207777d0, 
	5    -0.4033373d0,-0.07555636d0,-7.399682d-3,-2.943396d-4,
	6    8.628455d-3,-3.251305d0,-3.999623d0,-2.33175d0,
	7    -0.7614046d0,-0.1402496d0,-0.01354059d0,-5.3155d-4,
	8    0.01244159d0,-5.976818d0,-6.855045d0,-3.88775d0,
	9    -1.270677d0,-0.2370768d0,-0.02325118d0,-9.265136d-4,
	1    0.02204591d0,-9.495636d0,-10.05705d0,-5.636636d0,
	2    -1.883845d0,-0.3614146d0,-0.03629659d0,-1.473118d-3,
	3    0.03228755d0,-13.92918d0,-14.37232d0,-8.418409d0,
	4    -2.948277d0,-0.5819409d0,-0.059275d0,-2.419946d-3/
	e=ae
	em=aem
	a=aa
	V=E/EM
        prm=0.938
        pm=0.105658*0.105658
        if(em.le.dsqrt(pm)) go to 1
        vmin=0.14/em
        vmax=(2.*prm*em*em-2.*prm*pm)/(2.*prm*em+pm+prm*prm)/em
        if(v.lt.vmin.or.v.gt.vmax) go to 1
        em1=em-e
        if(em1.le.dsqrt(pm)) go to 1
	SIGMA=114.3+1.647*(dLOG(0.0213*E))**2
	IF(E.LE.17.) SIGMA=96.1+82./dSQRT(E)
	X5=0.00282*SIGMA*A**(1./3.)
	GX=3./X5**3*(X5**2/2.-1.+EXP(-X5)*(1.+X5))
	if(a.le.1.5) gx=1.
	T1=PM*V**2/(1.-V)
	T2=1.-2./V+2./V/V
	T3=0.54
	T4=1.8
	Y1=7.E-10*SIGMA*V*V
	Y2=T2*dLOG(1.+T3/T1)-T2*T3/(T3+T1)-2.*PM/T1
	1    *(T3+2.*T1)/(T3+T1)+4.*PM/T3*dlog(1.+T3/T1) 
	Y3=(T2+2.*PM/T4)*dLOG(1.+T4/T1)-2.*PM/T1
	Y4=T3/(T3+T1)
	Y5=T4/T1*dLOG(1.+T1/T4)
	FN1=Y1*(0.75*GX*Y2+0.25*Y3+0.5*PM/T1*(0.75*GX*Y4+0.25*Y5))

	eml = dlog10(em)
	i = 1
 10	i = i + 1
	if(i.ge.7) go to 11
	if(eml.gt.emul(i)) go to 10
 11	ii = i - 1
	s1 = 0.d0
	s2 = 0.d0
	s = 0.d0
	do k = 1, 8
	   s1 = s1 + atau(k,ii)*(dlog10(v))**(1.d0*(k-1))
	   s2 = s2 + atau(k,ii+1)*(dlog10(v))**(1.d0*(k-1))
	end do
	if(s1.gt.0.d0.and.s2.gt.0.d0) then
c	s = s1 + (s2 - s1) / (emul(ii+1) - emul(ii)) * 
c	1    (eml - emul(ii))
	   s = dlog10(s1) + (dlog10(s2) - dlog10(s1)) / 
	1	(emul(ii+1) - emul(ii)) * (eml - emul(ii))
	   s = 10.d0**s
	   if(s.le.0.d0) s = 0.d0
	end if
	FN2 = s * 6.023d-7

	FN = FN1 + FN2

	IF(FN.LE.0.) FN=0.
	GO TO 2
1	FN=0.
2	CONTINUE
	afn=fn
	RETURN
	END


C	CALCULATION OF MUON INELASTIC SCATTERING CROSS-SECTION
C       FOLLOWING S.I DUTTA ET AL. (PRD 63 (2001) 094020)
C	INPUT PARAMETERS: MUON ENERGY LOSS E, MUON ENERGY EM,
C	ATOMIC MASS A
C	OUTPUT PARAMETER: VALUE OF V*DSIGMA/DV, WHERE V=E/EM
C
	SUBROUTINE SIGNU_f2(ae,aem,aa,afn)	
	implicit none
	real*4 ae,aem,aa,afn
	real*8 prm,pm,ppi,pi,pal,pm0,pmp,pmr,q02,q2,x,v,e,em,fn
	real*8 vmin,vmax,qmin,qmax,em1,qa2,dqa
	real*8 plam,cp1,cp2,cp3,pap1,pap2,pap3,bp1,bp2,bp3
	real*8 cr1,cr2,cr3,par1,par2,par3,br1,br2,br3
	real*8 W2,F2,F2p,F21,F22,px,A,s,The,qa1,sha,R,t
	real*8 xp,xr,par,pap,br,cr,bp,cp,sigma
	parameter (prm=0.93827)
        parameter (pm=0.105658)
	parameter (ppi=0.13957)
	parameter (pi=3.14159265)
	parameter (pal=1./137.)
	parameter (pm0=0.31985)
	parameter (pmp=49.457)
	parameter (pmr=0.15052)
	parameter (q02=0.46017)
	parameter (plam=0.06527)
	parameter (cp1=0.28067)
	parameter (cp2=0.22291)
	parameter (cp3=2.1979)
	parameter (pap1=-0.0808)
	parameter (pap2=-0.44812)
	parameter (pap3=1.1709)
	parameter (bp1=0.60243)
	parameter (bp2=1.3754)
	parameter (bp3=1.8439)
	parameter (cr1=0.80107)
	parameter (cr2=0.97307)
	parameter (cr3=3.4942)
	parameter (par1=0.58400)
	parameter (par2=0.37888)
	parameter (par3=2.6063)
	parameter (br1=0.10711)
	parameter (br2=1.9386)
	parameter (br3=0.49338)

	e=ae
	em=aem
	a=aa
	v=e/em
        if(em.le.pm) go to 1
c        vmin = ppi / em
c        vmax=(2.*prm*em*em-2.*prm*pm*pm)/(2.*prm*em+pm*pm+prm*prm)/em
	vmin = ( ( prm + ppi )**2 - prm * prm ) / 2. / prm / em
	vmax = 1. - pm / em
        if(v.lt.vmin.or.v.gt.vmax) go to 1
        em1=em-e
        if(em1.le.pm) go to 1

	qmin = pm * pm * v * v / ( 1. - v )
	qmax = 2. * prm * em * v - ( ( prm + ppi )**2 - prm * prm )
	qa2 = dlog10(qmin)
	dqa = 0.01
	s = 0.
 3	qa1 = qa2 + dqa / 2.
	q2 = 10.**qa1
	if (q2.ge.qmax) go to 5
	x = q2 / 2. / prm / e
	if (x.ge.1) go to 4

	if (x.le.0.0014) sha = A**(-0.1)
	if (x.gt.0.0014.and.x.le.0.04) 
	1    sha = A**( 0.069 * dlog10(x) + 0.097)
	if (x.ge.0.04) sha = 1.

	px = 1. - 1.85 * x + 2.45 * x * x - 2.35 * x * x * x + 
	1    x * x * x * x

	t = dlog(dlog((q2+q02)/plam)/dlog(q02/plam))
	W2 = prm * prm + 2. * prm * e - q2
	xp = ( q2 + pmp ) / ( q2 + pmp + W2 - prm * prm )
	xr = ( q2 + pmr ) / ( q2 + pmr + W2 - prm * prm )
	cr = cr1 + cr2 * t**cr3
	par = par1 + par2 * t**par3
	br = br1 + br2 * t**br3
	bp = bp1 + bp2 * t**bp3
	cp = cp1 + ( cp1 - cp2 ) * ( 1. / ( 1. + t**cp3 ) - 1. )
	pap = pap1 + ( pap1 - pap2 ) * ( 1. / ( 1. + t**pap3 ) - 1. )
	F21 = cp * xp**pap * ( 1. - x )**bp
	F22 = cr * xr**par * ( 1. - x )**br
	F2p = q2 / ( q2 + pm0 ) * ( F21 + F22 )

	F2 = sha * A / 2. * ( 1. + px ) * F2p
	if(A.le.1.5) F2 = A / 2. * ( 1. + px ) * F2p
	if(A.le.1.1) F2 = A * F2p

	The = 1. + 12. * q2 / ( q2 + 1. ) * 
	1    0.125 * 0.125 / ( 0.125 * 0.125 + x * x )
	R = 0.635 / dlog(q2/0.2/0.2) * The + 0.5747 / q2 -
	1    0.3534 / ( q2 * q2 + 0.09)

	if(R.lt.0.) R=0.

	sigma = 4. * pi * pal * pal / q2 / q2 / x * F2 *
	1    ( 1. - v - prm * x * v / 2. / em + 
	2    ( 1. - 2. * pm * pm / q2 ) *
	3    v * v * ( 1. + 4. * prm * prm * x * x / q2 ) /
	4    2. / ( 1. + R )) * q2 / 2. / prm / v / v / em  
	s = s + sigma * q2 * dqa * 2.303
 4	continue
	qa2 = qa2 + dqa
	if (qa2.le.10.) go to 3
 5	continue
	FN = s * v * 6.023e-4 / A * 0.38938

	IF(FN.LE.0.) FN=0.
	GO TO 2
1	FN=0.
2	CONTINUE
	afn=fn
	RETURN
	END


C	CALCULATION OF CROSS-SECTION OF DELTA-ELECTRON PRODUCTION BY MUON
C	INPUT PARAMETERS: ENERGY OF DELTA-ELECTRON E, MUON ENERGY EM,
C	NUCLEAR CHARGE Z, NUCLEAR MASS A
C	OUTPUT PARAMETER: VALUE OF V*DSIGMA/DV, WHERE V=E/EM
C
	SUBROUTINE SIGE(E,EM,Z,A,FE)
	fmmu=0.105655
	fme=0.000511
	gam=em/fmmu
	beta=sqrt(gam*gam-1.)/gam
	V=E/EM
	emax=2.*fme*(gam*gam-1)/(1.+2.*gam*fme/fmmu+fme*fme/fmmu/fmmu)
	vmax=emax/em
	IF(V.GE.VMAX) GO TO 1
	dfe=1./137./3.1416*(2.*alog(1.+2.*e/fme)*alog((2.*e*(1.-e/emax))/
     *	(fme*(1.-e/em)))-1./2.*(alog((2.*e*(1.-e/emax))/(fme*(1.-e/em))))
     *	**2-3./2.*(alog(1.+2.*e/fme))**2+1./2.*alog(1.+2.*e/fme)*
     *	alog((2.*em*(em-e)*fme)/(fmmu*fmmu*e)))
	FE1=1.535E-4*Z/A/beta/beta/E*
     *  (1.-beta*beta*E/emax+E*E/EM/EM/2.)
	fe=fe1*(1.+dfe)
	IF(FE.LE.0.) FE=0.
	GO TO 2
1	FE=0.
2	CONTINUE
	RETURN
	END


C	CALCULATION OF CROSS-SECTION OF DELTA-ELECTRON PRODUCTION BY MUON
C	INPUT PARAMETERS: ENERGY OF DELTA-ELECTRON E, MUON ENERGY EM,
C	NUCLEAR CHARGE Z, NUCLEAR MASS A
C	OUTPUT PARAMETER: VALUE OF V*DSIGMA/DV, WHERE V=E/EM
C
	SUBROUTINE SIGE1(E,EM,Z,A,FE)
	fmmu=0.105658
	fme=0.000511
	gam=em/fmmu
	beta=sqrt(gam*gam-1.)/gam
	V=E/EM
	emax=2.*fme*(gam*gam-1)/(1.+2.*gam*fme/fmmu+fme*fme/fmmu/fmmu)
	vmax=emax/em
	IF(V.GE.VMAX) GO TO 1
	dfe=1./137./3.1416*(2.*alog(1.+2.*e/fme)*alog((2.*e*(1.-e/emax))/
     *	(fme*(1.-e/em)))-1./2.*(alog((2.*e*(1.-e/emax))/(fme*(1.-e/em))))
     *	**2-3./2.*(alog(1.+2.*e/fme))**2+1./2.*alog(1.+2.*e/fme)*
     *	alog((2.*em*(em-e)*fme)/(fmmu*fmmu*e)))
	FE1=1.535E-4*Z/A/beta/beta/E*
     *  (1.-beta*beta*E/emax+E*E/EM/EM/2.)
	fe=fe1*dfe
***	fe=fe1
	IF(FE.LE.0.) FE=0.
	GO TO 2
1	FE=0.
2	CONTINUE
	RETURN
	END


