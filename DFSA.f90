!============================================================================================
!    DFSA - A Derivative-free Simulated Annealing Method for
!    bound constrained global optimization
!    Copyright (C) 2011  G.Liuzzi, S.Lucidi, V.Piccialli
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    G. Liuzzi, S. Lucidi, V. Piccialli, A. Sotgiu. A magnetic resonance device designed via 
!    global optimization techniques, Mathematical Programming, 101: 339-364 (2004)
!
!============================================================================================
module vincoli
     double precision, allocatable :: eps(:),constr(:),epsiniz(:),constr_old(:)
	 double precision viol,violz,viol_old
end module vincoli

module mod_best
	real*8, parameter				:: tolfeas = 1.d-4
	real*8							:: minfeas, fminfeas, bestamm
	real*8, allocatable			:: xminfeas(:), xbestamm(:)
	real*8							:: tempfeas, tempfob
end module mod_best

module mod_bounds
	integer :: n1,me1,mi1
end module mod_bounds

SUBROUTINE DFSA(rseed)
	use mod_bounds
	use vincoli
	use mod_best

	IMPLICIT NONE

	INTEGER						:: N, ME, MI, NF, numnf, NUMFAL, ISTOP
	INTEGER						:: LS, NPC, NLM, ICONTFAL, IDFAL, NUMITER
	INTEGER						:: NUMCAS, K, I, J, MAXLS, MODLS, nminloc
	INTEGER						:: MAXITER, MAXNF, PRINT_LEVEL
	INTEGER						:: RSEED

	DOUBLE PRECISION,  ALLOCATABLE	:: X(:), XOTT(:)
	DOUBLE PRECISION,  ALLOCATABLE	:: Z(:), VETT(:), XPROP(:)
	DOUBLE PRECISION,  ALLOCATABLE	:: XINF(:), XSUP(:)
	DOUBLE PRECISION,  ALLOCATABLE	:: constrain(:), cl(:), cu(:)
	DOUBLE PRECISION				:: TOL, fglob, violmax2
	DOUBLE PRECISION				:: RFAL, RNUMFAL, RMED, RMED2, FOTT, F
	DOUBLE PRECISION				:: FPROP, VAR2, VAR1, ZETA, TRES, DIFF
	DOUBLE PRECISION				:: TCOEFF

	REAL							:: RRAND
	real							:: time_begin, time_end, timetot

	CHARACTER ( LEN = 40 )			:: nomefun
	CHARACTER ( LEN = 2  )			:: CSTEP


	call setdim(N,Mi,Me)

	n1  = n
	me1 = me
	mi1 = mi

	allocate(xbestamm(n))
	allocate(constr(mi+me),eps(mi+me), epsiniz(mi+me),constr_old(mi+me))
	ALLOCATE(XINF(N),XSUP(N))
	ALLOCATE(X(N),XOTT(N),Z(N),VETT(N),XPROP(N))

	call FUNCTINIT(N,X,XINF,XSUP,nomefun,fglob)


!-----------------------------------------------
!	SET TOLERANCES AND BUDGETS
!-----------------------------------------------
	TOL         = 1.0D-4
	MAXITER     = 100000000
	maxnf       = (ceiling((1.3 * 500*(me+mi)*n)/100)+1)*100
	NUMFAL		= 0
	PRINT_LEVEL = 0
	TCOEFF      = 1.0D0

	fminfeas    = 1.d+30
	minfeas     = 1.d+30
	bestamm     = 1.d+30

!
!	-------
!	 Iniz.
!	-------
!
	
		LS       = 0
		NF       = 0
		NPC      = 0
		NLM      = 0
		RFAL     = 0.0D0
		RNUMFAL  = 0.0D0
		ICONTFAL = 0
		IDFAL    = 0  	
		NUMITER  = 0
		nminloc  = 0 

		RMED     = 0.D0
		RMED2    = 0.D0

!	-------------------------------------
!	  Generazione del campione originale
!	-------------------------------------
		FOTT     = 1.D+32
!--------------------------------------------------------
!	modificato il 23/9/2003
!--------------------------------------------------------
		NUMCAS   = 30
!		NUMCAS   = 1
!--------------------------------------------------------
		IF(PRINT_LEVEL >= 0) THEN
			WRITE(*,2070)
			WRITE(*,2080) 
			WRITE(1,2070)
			WRITE(1,2080) 
		ENDIF

		call cpu_time(time_begin)

		DO K=1,NUMCAS

			DO I=1,N
				!CALL RANDOM_NUMBER(RRAND)
				call randp(rseed)
				VETT(I)=DBLE(rseed)/2147483647.d0
				XPROP(I)=XINF(I)+VETT(I)*(XSUP(I)-XINF(I))
			ENDDO

			call fun_viol(n,mi,me,xprop,tempfeas)
			call funob(n,xprop,tempfob)
			if(fott > 1.d+30) then
				minfeas = tempfeas
				fminfeas = tempfob
				write(*,*) 'update minfeas ',fminfeas,minfeas
				if(minfeas <= tolfeas) then
					bestamm = fminfeas
					xbestamm = xprop
				else
					bestamm = 1.d+30
				endif
			else
				if(tempfeas < minfeas) then
					minfeas  = tempfeas
					fminfeas = tempfob
					write(*,*) 'update minfeas ',fminfeas,minfeas
				elseif((tempfeas <= minfeas).and.(tempfob < fminfeas)) then
					minfeas  = tempfeas
					fminfeas = tempfob
					write(*,*) 'update minfeas ',fminfeas,minfeas
				endif
				if((tempfeas <= tolfeas).and.(tempfob < bestamm)) then
					bestamm  = tempfob
					xbestamm = xprop
					write(*,*) 'update bestamm ',bestamm,tempfeas
				endif
			endif

		   call funct_vinc(n,mi,me,xprop,fprop,constr)
		   
		   do i = 1,mi
			 if(max(0.d0,constr(i)) < 1.d-0) then
				eps(i) = 1.d-3
			 else
				eps(i) = 1.d-1
			 endif
		   enddo
	  	   do i = 1,me
			 if(abs(constr(mi+i)) < 1.d-0) then
				eps(mi+i) = 1.d-3
			 else
				eps(mi+i) = 1.d-1
			 endif
		   enddo

			CALL FUNCT(N,XPROP,FPROP)
			
			NPC=NPC+1
			NF=NF+1
			RMED=((NPC-1.0D0)*RMED+FPROP)/NPC

			IF(NPC.LE.2) THEN
				RMED2=RMED2+FPROP**2
			ELSE
				RMED2=((NPC-2.)*RMED2+FPROP**2)/(NPC-1.0D0)
			ENDIF

			IF(FPROP.LT.FOTT) THEN
				CALL AGGIORNO_OTTIMO(N,XPROP,FPROP,XOTT,FOTT)
			ENDIF

			IF(PRINT_LEVEL >= 0) THEN
				WRITE(*,2130) NUMITER, NF, FPROP,'--'
				WRITE(1,2130) NUMITER, NF, FPROP,'--'
			ENDIF
		ENDDO

		IF(PRINT_LEVEL >= 2) THEN
			write(*,*)'fine generazione random iniziale nf = ', nf
			write(1,*)'fine generazione random iniziale nf = ', nf
		ENDIF

		VAR2=RMED2-(RMED**2)*(NPC/(NPC-1.))

	IF(PRINT_LEVEL >= 0) THEN
		WRITE(*,2130) NUMITER, NF, FOTT,'IN'
		WRITE(1,2130) NUMITER, NF, FOTT,'IN'
	ENDIF			

	DO I=1,N
		X(I)=XOTT(I)
	ENDDO

	MAXLS=100000000      
	MODLS=1

!	--------------------
!	  Ciclo principale
!	--------------------
	IF(PRINT_LEVEL >= 0) THEN
		WRITE(*,2070)
		WRITE(*,2080) 
		WRITE(1,2070)
		WRITE(1,2080) 
	ENDIF
2070 FORMAT( '    ITER      FOTT      SALTO       ALFA       ')
2080 FORMAT( '-----------------------------------------------')
!                 123456  +1.2345E-02    OK     +1.2345E-02  
!2090 FORMAT(1X,  I6, 2X,  I12,  2X,  ES11.4, 4X, A2, 5X,  ES11.4, 1X,   ES11.4, 1X, ES11.4 )
2130 FORMAT(1X,  I6, 2X,  I12,  2X,  ES11.4, 4X, A2 )

100	FORMAT(/,1X,' GENERAZ.  N = ', I9,/)
101	FORMAT(/,1X,'  Attuale valore ottimo = ',D13.6,/)
!102	FORMAT(3(:,1X,'XOTT(',I1,') =',D13.6,1X))
!109	FORMAT(3(:,1X,'XOTT(',I2,') =',D13.6,1X))
6	FORMAT(/,1X,' NF=',I15,4X,' NPC=',I7,4X,' NLM=',I7,4X,' LS=',I7,/)
96	FORMAT(/,1X,'   MED = ',D13.6,'    VAR = ',D13.6,/)

	CSTEP ='  '

8	CONTINUE


	call cpu_time(time_end)
	timetot = time_end - time_begin

	IF((NF >= MAXNF).OR.(((bestamm-fglob)/(1.d0+abs(fglob)))<=tol)) THEN
	!IF(.OR.(NUMITER >= MAXITER).OR.(NF >= MAXNF)) THEN
		IF(PRINT_LEVEL >= 0) THEN	 
			WRITE(1,100) NUMITER
			WRITE(*,100) NUMITER
			WRITE(1,*) FOTT
			WRITE(*,*) FOTT
			WRITE(1,101) FOTT
			WRITE(*,101) FOTT

			WRITE(1,6) NF,NPC,NLM,LS
			WRITE(6,6) NF,NPC,NLM,LS

			VAR1=DSQRT(VAR2)
			WRITE(6,96) RMED,VAR1
			WRITE(1,96) RMED,VAR1
		ENDIF

		if(bestamm < 1.d+30) then
			if(print_level > 0) write(*,*) 'computing violmax2'
			violmax2=0.0d0
			call fun_viol(n,mi,me,xbestamm,violmax2)
		else
			violmax2 = 1.d+30
		endif

		write(3,950) nomefun, n, mi, me,  timetot, bestamm, violmax2, fminfeas, minfeas, nf, nminloc
		write(*,950) nomefun, n, mi, me,  timetot, bestamm, violmax2, fminfeas, minfeas, nf, nminloc 
950		FORMAT(a40,3(' & ',i4),5(' & ', es16.8),2(' & ',i10),'\\')

		GO TO 170
	ENDIF

	!IF (FMIN /= FOTT) PAUSE

	IF(PRINT_LEVEL >= 1) THEN


		WRITE(*,*)'/---------------------------------\'
		WRITE(1,*)'/---------------------------------\'
		WRITE(*,*)' FOTT = ',FOTT
		WRITE(1,*)' FOTT = ',FOTT
!		WRITE(*,*) XOTT
!		WRITE(1,*) XOTT
		WRITE(*,*)'\---------------------------------/'
		WRITE(1,*)'\---------------------------------/'
		WRITE(*,*)
		WRITE(1,*)
	ENDIF


!--------------------------------------

	NUMITER = NUMITER + 1
	IF(PRINT_LEVEL >= 0) THEN
		IF(MOD(NUMITER,30)==0) THEN
			WRITE(*,2070)
			WRITE(*,2080) 
			WRITE(1,2070)
			WRITE(1,2080) 
		ENDIF
	ENDIF

!--------------------------------------
!	genera zeta random tra 0 e 1
!	e un punto xprop su cui calcolare
!	la f.ob.
!--------------------------------------
	!CALL RANDOM_NUMBER(RRAND)
	call randp(rseed)
	zeta=dble(rseed)/2147483647.d0


	DO I=1,N
		!CALL RANDOM_NUMBER(RRAND)
		call randp(rseed)
		vett(i)=dble(rseed)/2147483647.d0
		XPROP(I)=XINF(I)+VETT(I)*(XSUP(I)-XINF(I))
	ENDDO	

!	---------------------------------------------------------
!	30-09-2015 
!	Aggiunta stima dell'eps nel punto corrente xprop
!	---------------------------------------------------------

	   call funct_vinc(n,mi,me,xprop,fprop,constr)
	   
	   do i = 1,mi
		 if(max(0.d0,constr(i)) < 1.d-0) then
			eps(i) = 1.d-3
		 else
			eps(i) = 1.d-1
		 endif
	   enddo
  	   do i = 1,me
		 if(abs(constr(mi+i)) < 1.d-0) then
			eps(mi+i) = 1.d-3
		 else
			eps(mi+i) = 1.d-1
		 endif
	   enddo

	CALL FUNCT(N,XPROP,FPROP)

	NPC=NPC+1
	NF=NF+1

!103	FORMAT(/,1X,'      Proposta di salto')
!105	FORMAT(/,1X,'   Zeta =',D13.6,'   Tres =',D13.6,/)

!	----------------------------
!	  Calcolo media, varianza e
!	  la temperatura del SA
!	----------------------------

	RMED  = ((NPC-1.)*RMED+FPROP)/NPC
	RMED2 = ((NPC-2.)*RMED2+FPROP**2)/(NPC-1.)
	VAR2  = RMED2-(RMED**2)*(NPC/(NPC-1.))

	DIFF  = RMED-FOTT

	TRES  = DEXP(-TCOEFF*(RFAL/DIFF)*(FPROP-FOTT))

	IF(FPROP.LE.FOTT) TRES=1.0D0

!--------------------------------------
!	se zeta e' minore o uguale a tres
!	allora accetta il punto e fa
!	partire una min. locale
!--------------------------------------
!	IF(.FALSE.) THEN
	IF(ZETA.LE.TRES) THEN

		NUMFAL = 0
		IF(PRINT_LEVEL >= 2) THEN
			WRITE(1,405) 
			WRITE(*,405) 
		ENDIF
405		FORMAT(1X,'    Proposta di salto accettata ',/)
		DO I=1,N
			X(I)=XPROP(I)
		ENDDO

		NLM     = NLM+1

!-----------------------------------------------------------------------------
! fai la ricerca locale
!-----------------------------------------------------------------------------

		numnf = 0
		call sd_box(n,x,f,xinf,xsup,1.d-4,min(5000,maxnf-nf),numnf,-1,istop, 1)
		nf = nf + numnf
		nminloc = nminloc + 1

		call funob(n,x,tempfob)
		call fun_viol(n,mi,me,x,tempfeas)

		!write(*,*) ' '
		!write(*,*) ' -----------------------------------------------------------------------------'
		!write(*,*) ' fine minimizzazione locale  --- amm =',tempfeas, '  fo =',tempfob
		!write(*,*) ' --------------------------------------------------------------------------------'
		!write(*,*) ' '

		if(tempfeas < minfeas) then
			minfeas  = tempfeas
			fminfeas = tempfob
			!xminfeas = xtemp
			write(*,*) 'update minfeas ',fminfeas,minfeas
		elseif((tempfeas <= minfeas).and.(tempfob < fminfeas)) then
			minfeas  = tempfeas
			fminfeas = tempfob
			!xminfeas = xtemp
			write(*,*) 'update minfeas ',fminfeas,minfeas
		endif
		if((tempfeas <= tolfeas).and.(tempfob < bestamm)) then
			bestamm  = tempfob
			xbestamm = x
			write(*,*) 'update bestamm ',bestamm,tempfeas
		endif

		!CALL INTERFACCIA(N,X,F,ALFA,D,XINF,XSUP,FOTT,XOTT,ISTOP,MAXNF,fglob)

		IF(F.LT.FOTT) THEN
			CALL AGGIORNO_OTTIMO(N,X,F,XOTT,FOTT)
			ISTOP = 1
		ENDIF

	ELSE
!--------------------------------------
!	altrimenti (zeta > tres) scarta
!	il punto proposto e aggiorna
!	il cont. numfal
!--------------------------------------
		NUMFAL   = NUMFAL   + 1

		RNUMFAL  = RNUMFAL  + 1.0D0
		ICONTFAL = ICONTFAL + 1
		RFAL     = RFAL     + 1.0D0/RNUMFAL
		IDFAL    = 1
	ENDIF

	GO TO 8

170 CONTINUE

	deallocate(xbestamm)
	deallocate(constr,eps,epsiniz,constr_old)
	DEALLOCATE(XINF,XSUP)
	DEALLOCATE(X,XOTT,Z,VETT,XPROP)

END SUBROUTINE DFSA

!==============================================================================================
!==============================================================================================
!==============================================================================================

SUBROUTINE AGGIORNO_OTTIMO(N,XPROP,FPROP,XOTT,FOTT)
	IMPLICIT NONE

	INTEGER,	  INTENT(IN)		:: N
	DOUBLE PRECISION, INTENT(IN)		:: XPROP(N), FPROP
	DOUBLE PRECISION, INTENT(OUT)		:: XOTT(N),  FOTT
	INTEGER					:: I

	XOTT = XPROP
	FOTT = FPROP

	RETURN

!1000 FORMAT(1X,'XOTT(',I3,') = ',ES11.4)
!1010 FORMAT(1X,'FOTT      = ',ES11.4)
!1020 FORMAT(1X,ES11.4)

	RETURN

END SUBROUTINE AGGIORNO_OTTIMO

subroutine funct(n,x,f)
	implicit none
	integer n
	real*8  x(n),f

	call funct_pen(n,x,f)

	return
end subroutine funct

subroutine funct_pen(n,x,f)

    use vincoli
    use mod_bounds
	implicit none
	integer n,i
	real*8  x(n),f,fob,fmax

	!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer :: me, mi
	!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	me=me1
	mi=mi1

    call funct_vinc(n,mi, me, x,fob,constr) 

	fmax = 0.d0
	viol = 0.d0

	do i = 1,mi
		fmax = fmax + (max(0.d0,constr(i)/eps(i)))
		viol=max(viol,constr(i))
	enddo

	do i = 1, me
		fmax = fmax + (abs(constr(mi+i)/eps(mi+i)))
		viol=max(viol,abs(constr(mi+i)))
	enddo

	f = fob + fmax

	return
end subroutine funct_pen


subroutine funct_vinc(n, m, p, x, fob, constr)

	implicit none

	integer			 :: n, m, p
	double precision :: x(n), fob, constr(m+p)
	
	call funob(n, x, fob)
	if (m >0 ) 	call fconstrin(n,m,x,constr(1:m))
	if (p >0 ) call fconstreq(n,p,x,constr(m+1:m+p))


	return
end subroutine funct_vinc


subroutine fun_viol(n, m, p, x, viol)

	implicit none

	integer			 :: n, m, p, i
	double precision :: x(n), fob, constr(m+p), viol
	
	viol = 0.0d0
	if (m >0 ) 	call fconstrin(n,m,x,constr(1:m))
	if (p >0 ) call fconstreq(n,p,x,constr(m+1:m+p))
	do i = 1,m
		viol=max(viol,constr(i))
	enddo

	do i = 1, p
		viol=max(viol,abs(constr(m+i)))
	enddo


	return
end subroutine fun_viol

subroutine randp(ix)
!c     ------------------------------------------------------------------
!c     randp: Portable pseudo-random number generator.
!c            Reference: L. Schrage, "A More Portable Fortran
!c            Random Number Generator", ACM Transactions on
!c            Mathematical Software, June, 1979.
!c     ------------------------------------------------------------------
      implicit none

      integer*4 a,p,ix,b15,b16,xhi,xalo,leftlo,fhi,k
      data a/16807/,b15/32768/,b16/65536/,p/2147483647/

!c     ------------------------------------------------------------------
      xhi=ix/b16
      xalo=(ix-xhi*b16)*a
      leftlo=xalo/b16
      fhi=xhi*a+leftlo
      k=fhi/b15
      ix=(((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
      if (ix.lt.0) ix=ix+p
!c     ------------------------------------------------------------------
      return
end

