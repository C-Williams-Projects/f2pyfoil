C***********************************************************************
C    Module:  xoper.f
C 
C    Copyright (C) 2000 Mark Drela 
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************
C
      SUBROUTINE FCPMIN
C------------------------------------------------
C     Finds minimum Cp on dist for cavitation work
C------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      XCPMNI = X(1)
      XCPMNV = X(1)
      CPMNI = CPI(1)
      CPMNV = CPV(1)
C
      DO I = 2, N + NW
        IF(CPI(I) .LT. CPMNI) THEN
         XCPMNI = X(I)
         CPMNI = CPI(I)
        ENDIF
        IF(CPV(I) .LT. CPMNV) THEN
         XCPMNV = X(I)
         CPMNV = CPV(I)
        ENDIF
      ENDDO
C

      IF(LVISC)THEN
        CPMN  = CPMNV
        XCPMN = XCPMNV
      ELSE
        CPMN  = CPMNI
        XCPMN = XCPMNI
C
        CPMNV  = CPMNI
        XCPMNV = XCPMNI
      ENDIF
C
      RETURN
      END ! FCPMIN

      SUBROUTINE MRSHOW(LM,LR)
      INCLUDE 'XFOIL.INC'
      LOGICAL LM, LR
C
C---- Interactive display removed for library use
C      IF(LM .OR. LR) WRITE(*,*)
C      IF(LM) THEN
C       IF(MATYP.EQ.1) WRITE(*,1100) MINF1
C       IF(MATYP.EQ.2) WRITE(*,1100) MINF1, ' / sqrt(CL)'
C       IF(MATYP.EQ.3) WRITE(*,1100) MINF1, ' / CL'
C      ENDIF
C      IF(LR) THEN
C       IF(RETYP.EQ.1) WRITE(*,1200) INT(REINF1)
C       IF(RETYP.EQ.2) WRITE(*,1200) INT(REINF1), ' / sqrt(CL)'
C       IF(RETYP.EQ.3) WRITE(*,1200) INT(REINF1), ' / CL'
C      ENDIF
C
      RETURN
C
 1100 FORMAT(1X,'M  =' , F10.4, A)
 1200 FORMAT(1X,'Re =' , I10  , A)
      END ! MRSHOW

      SUBROUTINE MHINGE
C----------------------------------------------------
C     Calculates the hinge moment of the flap about
C     (XOF,YOF) by integrating surface pressures.
C----------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      IF(.NOT.LFLAP) THEN
C
        CALL GETXYF(X,XP,Y,YP,S,N, TOPS,BOTS,XOF,YOF)
        LFLAP = .TRUE.
C
      ELSE
C
C------ find top and bottom y at hinge x location
        TOPS = XOF
        BOTS = S(N) - XOF
        CALL SINVRT(TOPS,XOF,X,XP,S,N)      
        CALL SINVRT(BOTS,XOF,X,XP,S,N)      
C
      ENDIF
C
      TOPX = SEVAL(TOPS,X,XP,S,N)
      TOPY = SEVAL(TOPS,Y,YP,S,N)
      BOTX = SEVAL(BOTS,X,XP,S,N)
      BOTY = SEVAL(BOTS,Y,YP,S,N)
C
C
      HMOM = 0.
      HFX  = 0.
      HFY  = 0.
C
C---- integrate pressures on top and bottom sides of flap
      DO 20 I=2, N
        IF(S(I-1).GE.TOPS .AND. S(I).LE.BOTS) GO TO 20
C
         DX = X(I) - X(I-1)
         DY = Y(I) - Y(I-1)
         XMID = 0.5*(X(I)+X(I-1)) - XOF
         YMID = 0.5*(Y(I)+Y(I-1)) - YOF
         IF(LVISC) THEN
          PMID = 0.5*(CPV(I) + CPV(I-1))
         ELSE
          PMID = 0.5*(CPI(I) + CPI(I-1))
         ENDIF
         HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
         HFX  = HFX  - PMID* DY
         HFY  = HFY  + PMID* DX
   20 CONTINUE
C
C---- find S(I)..S(I-1) interval containing s=TOPS
      DO I=2, N
        IF(S(I).GT.TOPS) GO TO 31
      ENDDO
C
   31 CONTINUE
C---- add on top surface chunk TOPS..S(I-1),  missed in the DO 20 loop.
      DX = TOPX - X(I-1)
      DY = TOPY - Y(I-1)
      XMID = 0.5*(TOPX+X(I-1)) - XOF
      YMID = 0.5*(TOPY+Y(I-1)) - YOF
      IF(S(I) .NE. S(I-1)) THEN
       FRAC = (TOPS-S(I-1))/(S(I)-S(I-1))
      ELSE
       FRAC = 0.
      ENDIF
      IF(LVISC) THEN
       TOPP = CPV(I)*FRAC + CPV(I-1)*(1.0-FRAC)
       PMID = 0.5*(TOPP+CPV(I-1))
      ELSE
       TOPP = CPI(I)*FRAC + CPI(I-1)*(1.0-FRAC)
       PMID = 0.5*(TOPP+CPI(I-1))
      ENDIF
      HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
      HFX  = HFX  - PMID* DY
      HFY  = HFY  + PMID* DX
C
C---- add on inside flap surface contribution from hinge to top surface
      DX = XOF - TOPX
      DY = YOF - TOPY
      XMID = 0.5*(TOPX+XOF) - XOF
      YMID = 0.5*(TOPY+YOF) - YOF
      HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
      HFX  = HFX  - PMID* DY
      HFY  = HFY  + PMID* DX
C
C---- find S(I)..S(I-1) interval containing s=BOTS
      DO I=N, 2, -1
        IF(S(I-1).LT.BOTS) GO TO 41
      ENDDO
C
   41 CONTINUE
C---- add on bottom surface chunk BOTS..S(I),  missed in the DO 20 loop.
      DX = X(I) - BOTX
      DY = Y(I) - BOTY
      XMID = 0.5*(BOTX+X(I)) - XOF
      YMID = 0.5*(BOTY+Y(I)) - YOF
      IF(S(I) .NE. S(I-1)) THEN
       FRAC = (BOTS-S(I-1))/(S(I)-S(I-1))
      ELSE
       FRAC = 0.
      ENDIF
      IF(LVISC) THEN
       BOTP = CPV(I)*FRAC + CPV(I-1)*(1.0-FRAC)
       PMID = 0.5*(BOTP+CPV(I))
      ELSE
       BOTP = CPI(I)*FRAC + CPI(I-1)*(1.0-FRAC)
       PMID = 0.5*(BOTP+CPI(I))
      ENDIF
      HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
      HFX  = HFX  - PMID* DY
      HFY  = HFY  + PMID* DX
C
C---- add on inside flap surface contribution from hinge to bottom surface
      DX = BOTX - XOF
      DY = BOTY - YOF
      XMID = 0.5*(BOTX+XOF) - XOF
      YMID = 0.5*(BOTY+YOF) - YOF
      HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
      HFX  = HFX  - PMID* DY
      HFY  = HFY  + PMID* DX
C
C---- add on TE base thickness contribution
      DX = X(1) - X(N)
      DY = Y(1) - Y(N)
      XMID = 0.5*(X(1)+X(N)) - XOF
      YMID = 0.5*(Y(1)+Y(N)) - YOF
      IF(LVISC) THEN
       PMID = 0.5*(CPV(1)+CPV(N))
      ELSE
       PMID = 0.5*(CPI(1)+CPI(N))
      ENDIF
      HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
      HFX  = HFX  - PMID* DY
      HFY  = HFY  + PMID* DX
C
      RETURN
      END ! MHINGE

      SUBROUTINE SPECAL
C-----------------------------------
C     Converges to specified alpha.
C-----------------------------------
      INCLUDE 'XFOIL.INC'
      REAL MINF_CLM, MSQ_CLM
C
C---- calculate surface vorticity distributions for alpha = 0, 90 degrees
      IF(.NOT.LGAMU .OR. .NOT.LQAIJ) CALL GGCALC
C
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
C
C---- superimpose suitably weighted  alpha = 0, 90  distributions
      DO 50 I=1, N
        GAM(I)   =  COSA*GAMU(I,1) + SINA*GAMU(I,2)
        GAM_A(I) = -SINA*GAMU(I,1) + COSA*GAMU(I,2)
   50 CONTINUE
      PSIO = COSA*GAMU(N+1,1) + SINA*GAMU(N+1,2)
C
      CALL TECALC
      CALL QISET
C
C---- set initial guess for the Newton variable CLM
      CLM = 1.0
C
C---- set corresponding  M(CLM), Re(CLM)
      CALL MRCL(CLM,MINF_CLM,REINF_CLM)
      CALL COMSET
C
C---- set corresponding CL(M)
      CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &            CL,CM,CDP, CL_ALF,CL_MSQ)
C
C---- iterate on CLM
      DO 100 ITCL=1, 20
C
        MSQ_CLM = 2.0*MINF*MINF_CLM
        DCLM = (CL - CLM)/(1.0 - CL_MSQ*MSQ_CLM)
C
        CLM1 = CLM
        RLX = 1.0
C
C------ under-relaxation loop to avoid driving M(CL) above 1
        DO 90 IRLX=1, 12
C
          CLM = CLM1 + RLX*DCLM
C
C-------- set new freestream Mach M(CLM)
          CALL MRCL(CLM,MINF_CLM,REINF_CLM)
C
C-------- if Mach is OK, go do next Newton iteration
          IF(MATYP.EQ.1 .OR. MINF.EQ.0.0 .OR. MINF_CLM.NE.0.0) GO TO 91
C
          RLX = 0.5*RLX
   90   CONTINUE
   91   CONTINUE
C
C------ set new CL(M)
        CALL COMSET
        CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &              CL,CM,CDP,CL_ALF,CL_MSQ)
C
        IF(ABS(DCLM).LE.1.0E-6) GO TO 110
C
  100 CONTINUE
      IF(.NOT.LQUIET) WRITE(*,*) 'SPECAL:  Minf convergence failed'
  110 CONTINUE
C
C---- set final Mach, CL, Cp distributions, and hinge moment
      CALL MRCL(CL,MINF_CL,REINF_CL)
      CALL COMSET
      CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &            CL,CM,CDP, CL_ALF,CL_MSQ)
      CALL CPCALC(N,QINV,QINF,MINF,CPI)
      IF(LVISC) THEN
       CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV)
       CALL CPCALC(N+NW,QINV,QINF,MINF,CPI)
      ELSE
       CALL CPCALC(N,QINV,QINF,MINF,CPI)
      ENDIF
      IF(LFLAP) CALL MHINGE
C
      RETURN
      END ! SPECAL
 
      SUBROUTINE SPECCL
C-----------------------------------------
C     Converges to specified inviscid CL.
C-----------------------------------------
      INCLUDE 'XFOIL.INC'
C
C---- calculate surface vorticity distributions for alpha = 0, 90 degrees
      IF(.NOT.LGAMU .OR. .NOT.LQAIJ) CALL GGCALC
C
C---- set freestream Mach from specified CL -- Mach will be held fixed
      CALL MRCL(CLSPEC,MINF_CL,REINF_CL)
      CALL COMSET
C
C---- current alpha is the initial guess for Newton variable ALFA
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
      DO 10 I=1, N
        GAM(I)   =  COSA*GAMU(I,1) + SINA*GAMU(I,2)
        GAM_A(I) = -SINA*GAMU(I,1) + COSA*GAMU(I,2)
   10 CONTINUE
      PSIO = COSA*GAMU(N+1,1) + SINA*GAMU(N+1,2)
C
C---- get corresponding CL, CL_alpha, CL_Mach
      CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &            CL,CM,CDP, CL_ALF,CL_MSQ)
C
C---- Newton loop for alpha to get specified inviscid CL
      DO 100 ITAL=1, 20
C
        DALFA = (CLSPEC - CL) / CL_ALF
        RLX = 1.0
C
        ALFA = ALFA + RLX*DALFA
C
C------ set new surface speed distribution
        COSA = COS(ALFA)
        SINA = SIN(ALFA)
        DO 40 I=1, N
          GAM(I)   =  COSA*GAMU(I,1) + SINA*GAMU(I,2)
          GAM_A(I) = -SINA*GAMU(I,1) + COSA*GAMU(I,2)
   40   CONTINUE
        PSIO = COSA*GAMU(N+1,1) + SINA*GAMU(N+1,2)
C
C------ set new CL(alpha)
        CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &              CL,CM,CDP,CL_ALF,CL_MSQ)
C
        IF(ABS(DALFA).LE.1.0E-6) GO TO 110
  100 CONTINUE
      IF(.NOT.LQUIET) WRITE(*,*) 'SPECCL:  CL convergence failed'
  110 CONTINUE
C
C---- set final surface speed and Cp distributions
      CALL TECALC
      CALL QISET
      IF(LVISC) THEN
       CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV)
       CALL CPCALC(N+NW,QINV,QINF,MINF,CPI)
      ELSE
       CALL CPCALC(N,QINV,QINF,MINF,CPI)
      ENDIF
      IF(LFLAP) CALL MHINGE
C
      RETURN
      END ! SPECCL

      SUBROUTINE VISCAL(NITER1)
C----------------------------------------
C     Converges viscous operating point
C----------------------------------------
      INCLUDE 'XFOIL.INC'
C
C---- convergence tolerance
      DATA EPS1 / 1.0E-4 /
C
      NITER = NITER1
C
C---- calculate wake trajectory from current inviscid solution if necessary
      IF(.NOT.LWAKE) THEN
       CALL XYWAKE
      ENDIF
C
C---- set velocities on wake from airfoil vorticity for alpha=0, 90
      CALL QWCALC
C
C---- set velocities on airfoil and wake for initial alpha
      CALL QISET
C
      IF(.NOT.LIPAN) THEN
C
       IF(LBLINI) CALL GAMQV
C
C----- locate stagnation point arc length position and panel index
       CALL STFIND
C
C----- set  BL position -> panel position  pointers
       CALL IBLPAN
C
C----- calculate surface arc length array for current stagnation point location
       CALL XICALC
C
C----- set  BL position -> system line  pointers
       CALL IBLSYS
C
      ENDIF
C
C---- set inviscid BL edge velocity UINV from QINV
      CALL UICALC
C
      IF(.NOT.LBLINI) THEN
C
C----- set initial Ue from inviscid Ue
       DO IBL=1, NBL(1)
         UEDG(IBL,1) = UINV(IBL,1)
       ENDDO
C
       DO IBL=1, NBL(2)
         UEDG(IBL,2) = UINV(IBL,2)
       ENDDO
C
      ENDIF
C
      IF(LVCONV) THEN
C----- set correct CL if converged point exists
       CALL QVFUE
       IF(LVISC) THEN
        CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV)
        CALL CPCALC(N+NW,QINV,QINF,MINF,CPI)
       ELSE
        CALL CPCALC(N,QINV,QINF,MINF,CPI)
       ENDIF
       CALL GAMQV
       CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &             CL,CM,CDP, CL_ALF,CL_MSQ)
       CALL CDCALC
      ENDIF
C
C---- set up source influence matrix if it doesn't exist
      IF(.NOT.LWDIJ .OR. .NOT.LADIJ) CALL QDCALC
C
C---- Newton iteration for entire BL solution
      IF(NITER.EQ.0) NITER= 50
C      WRITE(*,*)
C      WRITE(*,*) 'Solving BL system ...'
      DO 1000 ITER=1, NITER
C
C------ fill Newton system for BL variables
        CALL SETBL
C
C------ solve Newton system with custom solver
        CALL BLSOLV
C
C------ update BL variables
        CALL UPDATE
C
        IF(LALFA) THEN
C------- set new freestream Mach, Re from new CL
         CALL MRCL(CL,MINF_CL,REINF_CL)
         CALL COMSET
        ELSE
C------- set new inviscid speeds QINV and UINV for new alpha
         CALL QISET
         CALL UICALC
        ENDIF
C
C------ calculate edge velocities QVIS(.) from UEDG(..)
        CALL QVFUE
C
C------ set GAM distribution from QVIS
        CALL GAMQV
C
C------ relocate stagnation point
        CALL STMOVE
C
C------ set updated CL,CD
        CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &              CL,CM,CDP,CL_ALF,CL_MSQ)
        CALL CDCALC
C
C------ display changes and test for convergence
        IF(.NOT.LQUIET) THEN
         IF(RLX.LT.1.0)
     &    WRITE(*,2000) ITER, RMSBL, RMXBL, VMXBL,IMXBL,ISMXBL,RLX
         IF(RLX.EQ.1.0)
     &    WRITE(*,2010) ITER, RMSBL, RMXBL, VMXBL,IMXBL,ISMXBL
         CDP = CD - CDF
         WRITE(*,2020) ALFA/DTOR, CL, CM, CD, CDF, CDP
        ENDIF
C
        IF(RMSBL .LT. EPS1) THEN
         LVCONV = .TRUE.
         AVISC = ALFA
         MVISC = MINF
         GO TO 90
        ENDIF
C
 1000 CONTINUE
      IF(.NOT.LQUIET) WRITE(*,*) 'VISCAL:  Convergence failed'
C
   90 CONTINUE
      CALL CPCALC(N+NW,QINV,QINF,MINF,CPI)
      CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV)
      IF(LFLAP) CALL MHINGE
      RETURN
C....................................................................
 2000   FORMAT
     &   (/1X,I3,'   rms: ',E10.4,'   max: ',E10.4,3X,A1,' at ',I4,I3,
     &     '   RLX:',F6.3)
 2010   FORMAT
     &   (/1X,I3,'   rms: ',E10.4,'   max: ',E10.4,3X,A1,' at ',I4,I3)
 2020   FORMAT
     &   ( 1X,3X,'   a =', F7.3,'      CL =',F8.4  /
     &     1X,3X,'  Cm =', F8.4, '     CD =',F9.5,
     &           '   =>   CDf =',F9.5,'    CDp =',F9.5)
      END ! VISCAL
