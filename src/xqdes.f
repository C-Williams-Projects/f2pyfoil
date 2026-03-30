C***********************************************************************
C    Module:  xqdes.f
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
      SUBROUTINE SPLQSP(KQSP)
C------------------------------------------------------
C     Splines Qspec(s).  The end intervals are treated
C     specially to avoid Gibbs-type problems from 
C     blindly splining to the stagnation point.
C------------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
C---- usual spline with natural end BCs
      CALL SPLIND(QSPEC(2,KQSP),QSPECP(2,KQSP),SSPEC(2),NSP-2,
     &            -999.0,-999.0)
C
ccC---- pseudo-monotonic spline with simple secant slope calculation
cc      CALL SPLINA(QSPEC(2,KQSP),QSPECP(2,KQSP),SSPEC(2),NSP-2)
C
C---- end intervals are splined separately with natural BCs at
C     the trailing edge and matching slopes at the interior points
C
      I = 1
      CALL SPLIND(QSPEC(I,KQSP),QSPECP(I,KQSP),SSPEC(I),2,
     &            -999.0,QSPECP(I+1,KQSP))
C
      I = NSP-1
      CALL SPLIND(QSPEC(I,KQSP),QSPECP(I,KQSP),SSPEC(I),2,
     &            QSPECP(I,KQSP),-999.0)
C
      RETURN
      END


      SUBROUTINE SMOOQ(KQ1,KQ2,KQSP)
C--------------------------------------------
C     Smooths Qspec(s) inside target segment
C--------------------------------------------
      INCLUDE 'XFOIL.INC'
C
cC---- calculate smoothing coordinate
ccc      IF(NSP.EQ.NC1) THEN
cC
cC------ mapping inverse: use circle plane coordinate
c        I = 1
c        W8(I) = 0.0
c        DO 10 I=2, NSP
c          SINW = 2.0*SIN( 0.25*(WC(I)+WC(I-1)) )
c          SINWE = SINW**(1.0-AGTE)
cC
c          DSDW = SINWE * EXP( REAL(0.5*(PIQ(I)+PIQ(I-1)) ))
c          W8(I) = W8(I-1) + (WC(I)-WC(I-1))/DSDW
c   10   CONTINUE
c        DO 11 I=1, NSP
c          W8(I) = W8(I)/W8(NSP)
c 11     CONTINUE
cC
cC------ do not smooth first and last intervals in circle plane
c        KQ1 = MAX(IQ1,2)
c        KQ2 = MIN(IQ2,NSP-1)
cC
ccc      ELSE
C
C------ mixed inverse: use arc length coordinate
        DO 15 I=1, NSP
          W8(I) = SSPEC(I)
   15   CONTINUE
C
ccc      ENDIF
C
C
      IF(KQ2-KQ1 .LT. 2) THEN
       IF(.NOT.LQUIET)
     &  WRITE(*,*) 'Segment is too short.  No smoothing possible.'
       RETURN
      ENDIF
C
C---- set smoothing length ( ~ distance over which data is smeared )
      SMOOL = 0.002*(W8(NSP) - W8(1))
CCC   CALL ASKR('Enter Qspec smoothing length^',SMOOL)
C
C---- set up tri-diagonal system for smoothed Qspec
      SMOOSQ = SMOOL**2
      DO 20 I=KQ1+1, KQ2-1
        DSM = W8(I  ) - W8(I-1)
        DSP = W8(I+1) - W8(I  )
        DSO = 0.5*(W8(I+1) - W8(I-1))
C
        W1(I) =  SMOOSQ * (         - 1.0/DSM) / DSO
        W2(I) =  SMOOSQ * ( 1.0/DSP + 1.0/DSM) / DSO  +  1.0
        W3(I) =  SMOOSQ * (-1.0/DSP          ) / DSO
   20 CONTINUE
C
C---- set fixed-Qspec end conditions
      W2(KQ1) = 1.0
      W3(KQ1) = 0.0
C
      W1(KQ2) = 0.0
      W2(KQ2) = 1.0
C
      IF(LQSLOP) THEN
C----- also enforce slope matching at endpoints
       I = KQ1 + 1
       DSM = W8(I  ) - W8(I-1)
       DSP = W8(I+1) - W8(I  )
       DS  = W8(I+1) - W8(I-1)
       W1(I) = -1.0/DSM - (DSM/DS)/DSM
       W2(I) =  1.0/DSM + (DSM/DS)/DSM + (DSM/DS)/DSP
       W3(I) =                         - (DSM/DS)/DSP
       QSPP1 = W1(I)*QSPEC(I-1,KQSP)
     &       + W2(I)*QSPEC(I  ,KQSP)
     &       + W3(I)*QSPEC(I+1,KQSP)
C
       I = KQ2 - 1
       DSM = W8(I  ) - W8(I-1)
       DSP = W8(I+1) - W8(I  )
       DS  = W8(I+1) - W8(I-1)
       W1(I) =                           (DSP/DS)/DSM
       W2(I) = -1.0/DSP - (DSP/DS)/DSP - (DSP/DS)/DSM
       W3(I) =  1.0/DSP + (DSP/DS)/DSP
       QSPP2 = W1(I)*QSPEC(I-1,KQSP)
     &       + W2(I)*QSPEC(I  ,KQSP)
     &       + W3(I)*QSPEC(I+1,KQSP)
C
       QSPEC(KQ1+1,KQSP) = QSPP1
       QSPEC(KQ2-1,KQSP) = QSPP2
      ENDIF
C
C
C---- solve for smoothed Qspec array
      CALL TRISOL(W2(KQ1),W1(KQ1),W3(KQ1),QSPEC(KQ1,KQSP),(KQ2-KQ1+1))
C
C
cc      IF(LQSYM) THEN
cc        DO 40 I=KQ1+1, KQ2-1
cc          QSPEC(NSP-I+1,KQSP) = -QSPEC(I,KQSP)
cc 40     CONTINUE
cc      ENDIF
C
      RETURN
      END
 

      FUNCTION QINCOM(QC,QINF,TKLAM)
C-------------------------------------
C     Sets incompressible speed from
C     Karman-Tsien compressible speed
C-------------------------------------
C
      IF(TKLAM.LT.1.0E-4 .OR. ABS(QC).LT.1.0E-4) THEN
C----- for nearly incompressible case or very small speed, use asymptotic
C      expansion of singular quadratic formula to avoid numerical problems
       QINCOM = QC/(1.0 - TKLAM)
      ELSE
C----- use quadratic formula for typical case
       TMP = 0.5*(1.0 - TKLAM)*QINF/(QC*TKLAM)
       QINCOM = QINF*TMP*(SQRT(1.0 + 1.0/(TKLAM*TMP**2)) - 1.0)
      ENDIF
      RETURN
      END 

      SUBROUTINE GAMQSP(KQSP)
C------------------------------------------------
C     Sets Qspec(s,k) from current speed Q(s).
C------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      ALQSP(KQSP) = ALGAM
      CLQSP(KQSP) = CLGAM
      CMQSP(KQSP) = CMGAM
C
      DO 10 I=1, NSP
        QSPEC(I,KQSP) = QGAMM(I)
 10   CONTINUE
C
C---- zero out Qspec DOFs
      QDOF0 = 0.0
      QDOF1 = 0.0
      QDOF2 = 0.0
      QDOF3 = 0.0
C
      CALL SPLQSP(KQSP)
C
C---- reset target segment endpoints
      IF(.NOT.LIQSET) THEN
       IQ1 = 1
       IQ2 = NSP
      ENDIF
C
      RETURN
      END

      SUBROUTINE SYMQSP(KQSP)
C-----------------------------------------
C     Forces symmetry of Qspec(KQSP) array
C-----------------------------------------
      INCLUDE 'XFOIL.INC'
C
      ALQSP(KQSP) = 0.
      CLQSP(KQSP) = 0.
      CMQSP(KQSP) = 0.
C
      SSPMID = 0.5*(SSPEC(NSP) - SSPEC(1))
      DO 10 I=1, (NSP+1)/2
        SSPEC(I) = SSPMID + 0.5*(SSPEC(I)      - SSPEC(NSP-I+1)  )
        QSPEC(I,KQSP) =     0.5*(QSPEC(I,KQSP) - QSPEC(NSP-I+1,KQSP))
 10   CONTINUE
C
      DO 15 I=(NSP+1)/2+1, NSP
        SSPEC(I)      = -SSPEC(NSP-I+1)      + 2.0*SSPMID
        QSPEC(I,KQSP) = -QSPEC(NSP-I+1,KQSP)
 15   CONTINUE
C
C---- zero out Qspec DOFs
      QDOF0 = 0.0
      QDOF1 = 0.0
      QDOF2 = 0.0
      QDOF3 = 0.0
C
      CALL SPLQSP(KQSP)
C
C      WRITE(*,1000) KQSP
 1000 FORMAT(/' Qspec',I2,'  made symmetric')
C
      RETURN
      END

      SUBROUTINE MIXED(KQSP,NITERQ)
C-------------------------------------------------
C     Performs a mixed-inverse calculation using 
C     the specified surface speed array QSPEC.
C-------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
C---- distance of internal control point ahead of sharp TE
C-    (fraction of smaller panel length adjacent to TE)
      BWT = 0.1
C
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
      CALL SCALC(X,Y,S,N)
C
C---- zero-out and set DOF shape functions
      DO 1 I=1, N
        QF0(I) = 0.0
        QF1(I) = 0.0
        QF2(I) = 0.0
        QF3(I) = 0.0
    1 CONTINUE
C
C---- set DOF shape functions and specified speed
      DO 2 I=IQ1, IQ2
        FS = (S(I)-S(IQ1)) / (S(IQ2)-S(IQ1))
CCC        QF0(I) = (1.0-FS)**2
CCC        QF1(I) = FS**2
        QF0(I) = 1.0 - FS
        QF1(I) = FS
        IF(LCPXX) THEN
         QF2(I) = EXP(-5.0*     FS )
         QF3(I) = EXP(-5.0*(1.0-FS))
        ELSE
         QF2(I) = 0.0
         QF3(I) = 0.0
        ENDIF
        GAM(I) = QSPEC(I,KQSP) + QDOF0*QF0(I) + QDOF1*QF1(I)
     &                         + QDOF2*QF2(I) + QDOF3*QF3(I)
    2 CONTINUE
C
   99 CONTINUE
C
C---- perform Newton iterations on the new geometry
      DO 1000 ITER=1, NITERQ
C
      DO 3 I=1, N+5
        DO 31 J=1, N+5
          Q(I,J) = 0.
   31   CONTINUE
    3 CONTINUE
C
C---- calculate normal direction vectors along which the nodes move
      CALL NCALC(X,Y,S,N,NX,NY)
C
C---- go over all nodes, setting up  Psi = Psi0  equations
      DO 20 I=1, N
        CALL PSILIN(I,X(I),Y(I),NX(I),NY(I),PSI,PSI_N,.TRUE.,.FALSE.)
C
        DZDN(I) = DZDN(I) + PSI_N
C
C------ fill columns for specified geometry location
        DO 201 J=1, IQ1-1
          Q(I,J) = Q(I,J) + DZDG(J)
  201   CONTINUE
C
C------ fill columns for specified surface speed location
        DO 202 J=IQ1, IQ2
          Q(I,J) = Q(I,J) + DZDN(J)
  202   CONTINUE
C
C------ fill columns for specified geometry location
        DO 203 J=IQ2+1, N
          Q(I,J) = Q(I,J) + DZDG(J)
  203   CONTINUE
C
C------ set residual
        DQ(I) = PSIO - PSI
C
C------ fill global unknown columns
        Q(I,N+1) = Q(I,N+1) - 1.0
        Q(I,N+2) = Q(I,N+2) + Z_QDOF0
        Q(I,N+3) = Q(I,N+3) + Z_QDOF1
        Q(I,N+4) = Q(I,N+4) + Z_QDOF2
        Q(I,N+5) = Q(I,N+5) + Z_QDOF3
   20 CONTINUE
C
C---- set up Kutta condition
      DQ(N+1) = -( GAM(1) + GAM(N) )
      CALL GAMLIN(N+1,1,1.0)
      CALL GAMLIN(N+1,N,1.0)
C
      IF(SHARP) THEN
C----- set zero internal velocity in TE corner 
C
C----- set TE bisector angle
       AG1 = ATAN2(-YP(1),-XP(1)    )
       AG2 = ATANC( YP(N), XP(N),AG1)
       ABIS = 0.5*(AG1+AG2)
       CBIS = COS(ABIS)
       SBIS = SIN(ABIS)
C
C----- minimum panel length adjacent to TE
       DS1 = SQRT( (X(1)-X(2)  )**2 + (Y(1)-Y(2)  )**2 )
       DS2 = SQRT( (X(N)-X(N-1))**2 + (Y(N)-Y(N-1))**2 )
       DSMIN = MIN( DS1 , DS2 )
C
C----- control point on bisector just ahead of TE point
       XBIS = XTE - BWT*DSMIN*CBIS
       YBIS = YTE - BWT*DSMIN*SBIS
ccc       write(*,*) xbis, ybis
C
C----- set velocity component along bisector line
       CALL PSILIN(0,XBIS,YBIS,-SBIS,CBIS,PSI,QBIS,.FALSE.,.TRUE.)
C
CCC--- RES = DQDGj*Gamj + DQDMj*Massj + QINF*(COSA*CBIS + SINA*SBIS)
       RES = QBIS
C
       DO J=1, N+5
         Q(N,J) = 0.
       ENDDO
C
C----- dRes/dgamj
       DO J=1, N
         CALL GAMLIN(N,J, DQDG(J) )
         Q(N,J) = DQDG(J)
       ENDDO
C
C----- dRes/dPsio
       Q(N,N+1) = 0.
C
C----- -dRes/dUinf
       DQ(N) = -RES
      ENDIF
C
C---- pinned IQ1 point condition
      Q(N+2,IQ1) = 1.0
      DQ(N+2) = 0.0
C
C---- pinned IQ2 point condition
      Q(N+3,IQ2) = 1.0
      DQ(N+3) = 0.0
C
      IF(IQ1.GT.1 .AND. LCPXX) THEN
C----- speed regularity IQ1 condition
       RES = GAM(IQ1-1)      - 2.0*  GAM(IQ1)      +   GAM(IQ1+1)
     &  - (QSPEC(IQ1-1,KQSP) - 2.0*QSPEC(IQ1,KQSP) + QSPEC(IQ1+1,KQSP) )
       CALL GAMLIN(N+4,IQ1-1, 1.0)
       CALL GAMLIN(N+4,IQ1  ,-2.0)
       CALL GAMLIN(N+4,IQ1+1, 1.0)
       DQ(N+4) = -RES
      ELSE
C----- zero DOF condition
       Q(N+4,N+4) = 1.0
       DQ(N+4) = -QDOF2
      ENDIF
C
      IF(IQ2.LT.N .AND. LCPXX) THEN
C----- speed regularity IQ2 condition
       RES = GAM(IQ2-1)      - 2.0*  GAM(IQ2)      +   GAM(IQ2+1)
     &  - (QSPEC(IQ2-1,KQSP) - 2.0*QSPEC(IQ2,KQSP) + QSPEC(IQ2+1,KQSP) )
       CALL GAMLIN(N+5,IQ2-1, 1.0)
       CALL GAMLIN(N+5,IQ2  ,-2.0)
       CALL GAMLIN(N+5,IQ2+1, 1.0)
       DQ(N+5) = -RES
      ELSE
C----- zero DOF condition
       Q(N+5,N+5) = 1.0
       DQ(N+5) = -QDOF3
      ENDIF
C
      CALL GAUSS(IQX,N+5,Q,DQ,1)
C
      INMAX = 0
      IGMAX = 0
      DNMAX = 0.0
      DGMAX = 0.0
C
C---- update surface speed GAM before target segment
      DO 100 I=1, IQ1-1
        GAM(I) = GAM(I) + DQ(I)
        IF(ABS(DQ(I)) .GT. ABS(DGMAX)) THEN
         DGMAX = DQ(I)
         IGMAX = I
        ENDIF
  100 CONTINUE
C
C---- update panel nodes inside target segment
      DO 110 I=IQ1, IQ2
        X(I) = X(I) + NX(I)*DQ(I)
        Y(I) = Y(I) + NY(I)*DQ(I)
        IF(ABS(DQ(I)) .GT. ABS(DNMAX)) THEN
         DNMAX = DQ(I)
         INMAX = I
        ENDIF
  110 CONTINUE
C
C---- update surface speed GAM after target segment
      DO 120 I=IQ2+1, N
        GAM(I) = GAM(I) + DQ(I)
        IF(ABS(DQ(I)) .GT. ABS(DGMAX)) THEN
         DGMAX = DQ(I)
         IGMAX = I
        ENDIF
  120 CONTINUE
C
C---- update gloabal variables
      PSIO  = PSIO  + DQ(N+1)
      QDOF0 = QDOF0 + DQ(N+2)
      QDOF1 = QDOF1 + DQ(N+3)
      QDOF2 = QDOF2 + DQ(N+4)
      QDOF3 = QDOF3 + DQ(N+5)
C
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
      CALL SCALC(X,Y,S,N)
C
C---- set correct surface speed over target segment including DOF contributions
      DO 140 I=IQ1, IQ2
        GAM(I) = QSPEC(I,KQSP) + QDOF0*QF0(I) + QDOF1*QF1(I)
     &                         + QDOF2*QF2(I) + QDOF3*QF3(I)
  140 CONTINUE
C
C---- update everything else
      CALL TECALC
      CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &            CL,CM,CDP, CL_ALF,CL_MSQ)
      IF(.NOT.LQUIET) WRITE(*,2000) DNMAX,INMAX,DGMAX,IGMAX,CL
     &             ,DQ(N+2),DQ(N+3)
     &             ,DQ(N+4),DQ(N+5)
 2000 FORMAT(/' dNmax =',E10.3,I4,'   dQmax =',E10.3,I4,'    CL =',F7.4
     &       /' dQf1  =',E10.3,4X,'   dQf2  =',E10.3
     &       /' dQf3  =',E10.3,4X,'   dQf4  =',E10.3)
C
      IF(ABS(DNMAX).LT.5.0E-5 .AND. ABS(DGMAX).LT.5.0E-4) THEN
       RETURN
      ENDIF
C
 1000 CONTINUE
      IF(.NOT.LQUIET)
     & WRITE(*,*) 'Not quite converged.  Can EXEC again if necessary.'
      RETURN
C
      END
      
      SUBROUTINE GAMLIN(I,J,COEF)
C-------------------------------------------------------------------
C     Adds on Jacobian entry for point I due to node speed GAM at J.
C     GAM is either a local unknown if outside target segment,
C     or dependent on global Qspec DOF's if inside target segment.
C-------------------------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      IF(J.GE.IQ1 .AND. J.LE.IQ2) THEN
C----- inside target segment
       Q(I,N+2) = Q(I,N+2) + COEF*QF0(J)
       Q(I,N+3) = Q(I,N+3) + COEF*QF1(J)
       Q(I,N+4) = Q(I,N+4) + COEF*QF2(J)
       Q(I,N+5) = Q(I,N+5) + COEF*QF3(J)
      ELSE
C----- outside target segment
       Q(I,J) = Q(I,J) + COEF
      ENDIF
      RETURN
      END
