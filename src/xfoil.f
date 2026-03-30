C***********************************************************************
C    Module:  xfoil.f
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
      SUBROUTINE INIT
C---------------------------------------------------
C     Variable initialization/default routine.
C     See file XFOIL.INC for variable description.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      PI = 4.0*ATAN(1.0)
      HOPI = 0.50/PI
      QOPI = 0.25/PI
      DTOR = PI/180.0
C
C---- default Cp/Cv (air)
      GAMMA = 1.4
      GAMM1 = GAMMA - 1.0
C
C---- set unity freestream speed
      QINF = 1.0
C
C---- initialize freestream Mach number to zero
      MATYP = 1
      MINF1 = 0.
C
      ALFA = 0.0
      COSA = 1.0
      SINA = 0.0
C
      DO 10 I=1, IQX
        GAMU(I,1) = 0.
        GAMU(I,2) = 0.
        GAM(I) = 0.
        GAM_A(I) = 0.
   10 CONTINUE
      PSIO = 0.
C
      CL = 0.
      CM = 0.
      CD = 0.
C
      SIGTE = 0.0
      GAMTE = 0.0
      SIGTE_A = 0.
      GAMTE_A = 0.
C
      DO 20 I=1, IZX
        SIG(I) = 0.
   20 CONTINUE
C
      NQSP = 0
      DO 30 K=1, IPX
        ALQSP(K) = 0.
        CLQSP(K) = 0.
        CMQSP(K) = 0.
        DO 302 I=1, IBX
          QSPEC(I,K) = 0.
 302    CONTINUE
 30   CONTINUE
C
      AWAKE = 0.0
      AVISC = 0.0
C
      YIMAGE = -10.0
      LIMAGE = .FALSE.
C
      LGAMU  = .FALSE.
      LQINU  = .FALSE.
      LVISC  = .FALSE.
      LWAKE  = .FALSE.
      LBLINI = .FALSE.
      LIPAN  = .FALSE.
      LQAIJ  = .FALSE.
      LADIJ  = .FALSE.
      LWDIJ  = .FALSE.
      LCPXX  = .FALSE.
      LQVDES = .FALSE.
      LQSPEC = .FALSE.
      LQREFL = .FALSE.
      LVCONV = .FALSE.
      LCPREF = .FALSE.
      LFOREF = .FALSE.
      LPFILE = .FALSE.
      LPFILX = .FALSE.
      LBFLAP = .FALSE.
      LFLAP  = .FALSE.
      LEIW   = .FALSE.
      LSCINI = .FALSE.
      LCMINP = .FALSE.
      LHMOMP = .FALSE.
      LFREQP = .TRUE.
      LQUIET = .FALSE.
C
      LGSAME = .FALSE.
C
C
C
C---- input airfoil will not be normalized
      LNORM = .FALSE.
C
C---- airfoil will not be forced symmetric
      LQSYM = .FALSE.
      LGSYM = .FALSE.
C
C---- endpoint slopes will be matched
      LQSLOP = .TRUE.
      LGSLOP = .TRUE.
      LCSLOP = .TRUE.
C
C
C
C---- buffer and current airfoil flap hinge coordinates
      XBF = 0.0
      YBF = 0.0
      XOF = 0.0
      YOF = 0.0
C
C
C                                               n
C---- circle plane array size (257, or largest 2  + 1 that will fit array size)
      ANN = LOG(FLOAT((2*IQX)-1))/LOG(2.0)
      NN = INT( ANN + 0.00001 )
      NC1 = 2**NN + 1
      NC1 = MIN( NC1 , 257 )
C
C---- default paneling parameters
      NPAN = 160
      CVPAR = 1.0
      CTERAT = 0.15
      CTRRAT = 0.2
C
C---- default paneling refinement zone x/c endpoints
      XSREF1 = 1.0
      XSREF2 = 1.0
      XPREF1 = 1.0
      XPREF2 = 1.0
C
C
      LAECEN = .FALSE.
C
C---- default Cm reference location
      XCMREF = 0.25
      YCMREF = 0.
C
C---- default viscous parameters
      RETYP = 1
      REINF1 = 0.
      ACRIT = 9.0
      XSTRIP(1) = 1.0
      XSTRIP(2) = 1.0
      XOCTR(1) = 1.0
      XOCTR(2) = 1.0
      YOCTR(1) = 0.
      YOCTR(2) = 0.
      WAKLEN = 1.0

      ICCOR = 1   ! default: Karman-Tsien (existing behaviour)
C
C---- set BL calibration parameters
      CALL BLPINI
C
C---- Newton iteration limit
      ITMAX = 10
C
C
C---- drop tolerance for BL system solver
      VACCEL = 0.01
C
C---- inverse-mapping auto-filter level
      FFILT = 0.0
C
C
C
      NNAME  = 32
      NAME   = '                                '
CCC             12345678901234567890123456789012
C
C---- MSES domain parameters (not used in XFOIL)
      ISPARS = ' -2.0  3.0  -2.5  3.5'
C
C---- set MINF, REINF, based on current CL-dependence
      CALL MRCL(1.0,MINF_CL,REINF_CL)
C
C---- set various compressibility parameters from MINF
      CALL COMSET
C
      RETURN
      END ! INIT


      SUBROUTINE MRCL(CLS,M_CLS,R_CLS)
C-------------------------------------------
C     Sets actual Mach, Reynolds numbers
C     from unit-CL values and specified CLS
C     depending on MATYP,RETYP flags.
C-------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL M_CLS
C
      CLA = MAX( CLS , 0.000001 )
C
      IF(RETYP.LT.1 .OR. RETYP.GT.3) THEN
        IF(.NOT.LQUIET) THEN
         WRITE(*,*) 'MRCL:  Illegal Re(CL) dependence trigger.'
         WRITE(*,*) '       Setting fixed Re.'
        ENDIF
        RETYP = 1
      ENDIF
      IF(MATYP.LT.1 .OR. MATYP.GT.3) THEN
        IF(.NOT.LQUIET) THEN
         WRITE(*,*) 'MRCL:  Illegal Mach(CL) dependence trigger.'
         WRITE(*,*) '       Setting fixed Mach.'
        ENDIF
        MATYP = 1
      ENDIF
C
C
      IF(MATYP.EQ.1) THEN
C
        MINF  = MINF1
        M_CLS = 0.
C
      ELSE IF(MATYP.EQ.2) THEN
C
        MINF  =  MINF1/SQRT(CLA)
        M_CLS = -0.5*MINF/CLA
C
      ELSE IF(MATYP.EQ.3) THEN
C
        MINF  = MINF1
        M_CLS = 0.
C
      ENDIF
C
C
      IF(RETYP.EQ.1) THEN
C
        REINF = REINF1
        R_CLS = 0.
C
      ELSE IF(RETYP.EQ.2) THEN
C
        REINF =  REINF1/SQRT(CLA)
        R_CLS = -0.5*REINF/CLA
C
      ELSE IF(RETYP.EQ.3) THEN
C
        REINF =  REINF1/CLA
        R_CLS = -REINF /CLA
C
      ENDIF
C
C
      IF(MINF .GE. 0.99) THEN
        IF(.NOT.LQUIET) THEN
         WRITE(*,*)
         WRITE(*,*) 'MRCL: CL too low for chosen Mach(CL) dependence'
         WRITE(*,*) '      Aritificially limiting Mach to  0.99'
        ENDIF
        MINF = 0.99
        M_CLS = 0.
      ENDIF
C
      RRAT = 1.0
      IF(REINF1 .GT. 0.0) RRAT = REINF/REINF1
C
      IF(RRAT .GT. 100.0) THEN
        IF(.NOT.LQUIET) THEN
         WRITE(*,*)
         WRITE(*,*) 'MRCL: CL too low for chosen Re(CL) dependence'
         WRITE(*,*) '      Aritificially limiting Re to ',REINF1*100.0
        ENDIF
        REINF = REINF1*100.0
        R_CLS = 0.
      ENDIF
C
      RETURN
      END ! MRCL

      SUBROUTINE COMSET
      INCLUDE 'XFOIL.INC'
C
C---- Prandtl-Glauert compressibility factor (common to all corrections)
      BETA     = SQRT(1.0 - MINF**2)
      BETA_MSQ = -0.5/BETA
C
C=========================================================
C     Set TKLAM and TKL_MSQ.
C     These are used by the BL solver (xbl.f, GET_BL etc.)
C     for the velocity <-> Cp back-transformation.
C=========================================================
C
      IF(ICCOR .EQ. 0) THEN
C
C----   Prandtl-Glauert: linear; velocity ratio equals inviscid ratio
        TKLAM   = 0.0
        TKL_MSQ = 0.0
C
      ELSE IF(ICCOR .EQ. 1) THEN
C
C----   Karman-Tsien (original formulation)
        TKLAM   = MINF**2 / (1.0 + BETA)**2
        TKL_MSQ =     1.0 / (1.0 + BETA)**2
     &      - 2.0*TKLAM / (1.0 + BETA) * BETA_MSQ
C
      ELSE
C
C----   Laitone: scale KT TKLAM by the adiabatic factor (1 + (g-1)/2 * M^2)
C       This keeps the BL velocity-transformation structurally consistent
C       with the Laitone Cp denominator.
C
C       d(TKLAM)/d(Minf^2):
C         FLAIT = 1 + 0.5*GAMM1*Minf^2
C         d(FLAIT)/d(Minf^2) = 0.5*GAMM1
C         TKLAM = Minf^2*FLAIT / (1+beta)^2
C         => TKL_MSQ = (FLAIT + GAMM1*Minf^2) / (1+beta)^2
C                    - 2*TKLAM/(1+beta)*BETA_MSQ
C
        FLAIT   = 1.0 + 0.5*GAMM1*MINF**2
        TKLAM   = MINF**2 * FLAIT / (1.0 + BETA)**2
        TKL_MSQ = (FLAIT + GAMM1*MINF**2) / (1.0 + BETA)**2
     &      - 2.0*TKLAM / (1.0 + BETA) * BETA_MSQ
C
      ENDIF
C
C=========================================================
C     Set CPCBFAC and CPCBFQ for CPCALC / CLCALC.
C     Cp formula:  DEN = BETA + CPCBFAC*CPINC
C                  CP  = CPINC / DEN
C=========================================================
C
      IF(ICCOR .EQ. 0) THEN
C
C----   Prandtl-Glauert: DEN = BETA  (no CPINC dependence)
        CPCBFAC = 0.0
        CPCBFQ  = 0.0
C
      ELSE IF(ICCOR .EQ. 1) THEN
C
C----   Karman-Tsien
C       BFAC     = 0.5*M^2 / (1+beta)
C       dBFAC/d(M^2) = 0.5/(1+beta) - BFAC*BETA_MSQ/(1+beta)
C
        CPCBFAC = 0.5*MINF**2 / (1.0 + BETA)
        CPCBFQ  = 0.5 / (1.0 + BETA)
     &          - CPCBFAC * BETA_MSQ / (1.0 + BETA)
C
      ELSE
C
C----   Laitone
C       BFAC     = 0.5*M^2 * (1 + 0.5*(g-1)*M^2) / beta
C       dBFAC/d(M^2) = (1 + (g-1)*M^2) / (2*beta)
C                    - BFAC * BETA_MSQ / beta
C
        CPCBFAC = 0.5*MINF**2 * (1.0 + 0.5*GAMM1*MINF**2) / BETA
        CPCBFQ  = (1.0 + GAMM1*MINF**2) / (2.0*BETA)
     &          - CPCBFAC * BETA_MSQ / BETA
C
      ENDIF
C
C---- set sonic pressure coefficient and speed (isentropic)
      IF(MINF.EQ.0.0) THEN
       CPSTAR = -999.0
       QSTAR  = 999.0
      ELSE
       CPSTAR = 2.0 / (GAMMA*MINF**2)
     &        * (( (1.0 + 0.5*GAMM1*MINF**2)
     &            /(1.0 + 0.5*GAMM1        ))**(GAMMA/GAMM1) - 1.0)
       QSTAR  = QINF/MINF
     &        * SQRT( (1.0 + 0.5*GAMM1*MINF**2)
     &               /(1.0 + 0.5*GAMM1        ) )
      ENDIF
C
      RETURN
      END ! COMSET


      SUBROUTINE CPCALC(N,Q,QINF,MINF,CP)
C---------------------------------------------
C     Sets compressible Cp from speed array Q.
C     Correction type is selected via CPCBFAC
C     pre-computed in COMSET (see CPCOR.INC).
C
C     PG:      DEN = BETA
C     KT:      DEN = BETA + [M^2/(2(1+beta))]   * CPINC
C     Laitone: DEN = BETA + [M^2(1+(g-1)/2*M^2)
C                            /(2*beta)]           * CPINC
C---------------------------------------------
      DIMENSION Q(N), CP(N)
      REAL MINF
      INCLUDE 'QUIET.INC'
      INCLUDE 'CPCOR.INC'
C
      LOGICAL DENNEG
C
      BETA   = SQRT(1.0 - MINF**2)
C
      DENNEG = .FALSE.
C
      DO 20 I=1, N
        CPINC  = 1.0 - (Q(I)/QINF)**2
        DEN    = BETA + CPCBFAC*CPINC
        CP(I)  = CPINC / DEN
        IF(DEN .LE. 0.0) DENNEG = .TRUE.
   20 CONTINUE
C
      IF(DENNEG .AND. .NOT.LQUIET) THEN
       WRITE(*,*)
       WRITE(*,*) 'CPCALC: Local speed too large. ',
     &            'Compressibility corrections invalid.'
      ENDIF
C
      RETURN
      END ! CPCALC

 
      SUBROUTINE CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF,
     &                  XREF,YREF,
     &                  CL,CM,CDP, CL_ALF,CL_MSQ)
C-----------------------------------------------------------
C     Integrates surface pressures to get CL, CM, CDP.
C     Calculates dCL/dAlpha and dCL/d(Minf^2) for Newton.
C
C     Compressibility correction coefficients CPCBFAC and
C     CPCBFQ are set by COMSET (see CPCOR.INC).
C-----------------------------------------------------------
      DIMENSION X(N),Y(N), GAM(N), GAM_A(N)
      REAL MINF
      INCLUDE 'CPCOR.INC'
C
      SA = SIN(ALFA)
      CA = COS(ALFA)
C
      BETA     = SQRT(1.0 - MINF**2)
      BETA_MSQ = -0.5/BETA
C
C---- Cp correction coefficient and its Minf^2 sensitivity.
C     Previously these were hardcoded for Karman-Tsien only;
C     now they come from COMSET via CPCOR.INC so that PG,
C     KT and Laitone are all handled without branching here.
      BFAC     = CPCBFAC
      BFAC_MSQ = CPCBFQ
C
      CL     = 0.0
      CM     = 0.0
      CDP    = 0.0
      CL_ALF = 0.0
      CL_MSQ = 0.0
C
      I = 1
      CGINC    = 1.0 - (GAM(I)/QINF)**2
      CPG1     = CGINC/(BETA + BFAC*CGINC)
      CPG1_MSQ = -CPG1/(BETA + BFAC*CGINC)
     &           *(BETA_MSQ + BFAC_MSQ*CGINC)
C
      CPI_GAM  = -2.0*GAM(I)/QINF**2
      CPC_CPI  = (1.0 - BFAC*CPG1) / (BETA + BFAC*CGINC)
      CPG1_ALF = CPC_CPI*CPI_GAM*GAM_A(I)
C
      DO 10 I=1, N
        IP = I+1
        IF(I.EQ.N) IP = 1
C
        CGINC    = 1.0 - (GAM(IP)/QINF)**2
        CPG2     = CGINC/(BETA + BFAC*CGINC)
        CPG2_MSQ = -CPG2/(BETA + BFAC*CGINC)
     &             *(BETA_MSQ + BFAC_MSQ*CGINC)
C
        CPI_GAM  = -2.0*GAM(IP)/QINF**2
        CPC_CPI  = (1.0 - BFAC*CPG2) / (BETA + BFAC*CGINC)
        CPG2_ALF = CPC_CPI*CPI_GAM*GAM_A(IP)
C
        DX = (X(IP) - X(I))*CA + (Y(IP) - Y(I))*SA
        DY = (Y(IP) - Y(I))*CA - (X(IP) - X(I))*SA
        DG = CPG2 - CPG1
C
        AX = (0.5*(X(IP)+X(I))-XREF)*CA + (0.5*(Y(IP)+Y(I))-YREF)*SA
        AY = (0.5*(Y(IP)+Y(I))-YREF)*CA - (0.5*(X(IP)+X(I))-XREF)*SA
        AG = 0.5*(CPG2 + CPG1)
C
        DX_ALF = -(X(IP) - X(I))*SA + (Y(IP) - Y(I))*CA
        AG_ALF =  0.5*(CPG2_ALF + CPG1_ALF)
        AG_MSQ =  0.5*(CPG2_MSQ + CPG1_MSQ)
C
        CL     = CL     + DX* AG
        CDP    = CDP    - DY* AG
        CM     = CM     - DX*(AG*AX + DG*DX/12.0)
     &                  - DY*(AG*AY + DG*DY/12.0)
C
        CL_ALF = CL_ALF + DX*AG_ALF + AG*DX_ALF
        CL_MSQ = CL_MSQ + DX*AG_MSQ
C
        CPG1     = CPG2
        CPG1_ALF = CPG2_ALF
        CPG1_MSQ = CPG2_MSQ
   10 CONTINUE
C
      RETURN
      END ! CLCALC



      SUBROUTINE CDCALC
      INCLUDE 'XFOIL.INC'
C
      SA = SIN(ALFA)
      CA = COS(ALFA)
C
      IF(LVISC .AND. LBLINI) THEN
C
C----- set variables at the end of the wake
       THWAKE = THET(NBL(2),2)
       URAT   = UEDG(NBL(2),2)/QINF
       UEWAKE = UEDG(NBL(2),2) * (1.0-TKLAM) / (1.0 - TKLAM*URAT**2)
       SHWAKE = DSTR(NBL(2),2)/THET(NBL(2),2)
C
C----- extrapolate wake to downstream infinity using Squire-Young relation
C      (reduces errors of the wake not being long enough)
       CD = 2.0*THWAKE * (UEWAKE/QINF)**(0.5*(5.0+SHWAKE))
C
      ELSE
C
       CD = 0.0
C
      ENDIF
C
C---- calculate friction drag coefficient
      CDF = 0.0
      DO 20 IS=1, 2
        DO 205 IBL=3, IBLTE(IS)
          I  = IPAN(IBL  ,IS)
          IM = IPAN(IBL-1,IS)
          DX = (X(I) - X(IM))*CA + (Y(I) - Y(IM))*SA
          CDF = CDF + 0.5*(TAU(IBL,IS)+TAU(IBL-1,IS))*DX * 2.0/QINF**2
 205    CONTINUE
 20   CONTINUE
C
      RETURN
      END ! CDCALC

      SUBROUTINE ADDGEOM(XIN, YIN, NIN)
      INCLUDE 'XFOIL.INC'
C------------------------------------------------------
C     Loads airfoil from passed-in coordinate arrays
C     into the buffer airfoil and does initial
C     processing (orientation check, normalisation,
C     spline, geometry parameters, paneling).
C------------------------------------------------------
      INTEGER NIN
      DIMENSION XIN(NIN), YIN(NIN)
C
      DATA ANGTOL / 40.0 /
C
      IF(NIN.LE.1) THEN
       IF(.NOT.LQUIET) WRITE(*,*) 'ADDGEOM: No airfoil data supplied.'
       RETURN
      ENDIF
C
C---- copy incoming coordinates into global buffer airfoil arrays
      NB = NIN
      DO 10 I=1, NB
        XB(I) = XIN(I)
        YB(I) = YIN(I)
   10 CONTINUE
C
C---- calculate airfoil area assuming counterclockwise ordering
      AREA = 0.0
      DO 50 I=1, NB
        IP = I+1
        IF(I.EQ.NB) IP = 1
        AREA = AREA + 0.5*(YB(I)+YB(IP))*(XB(I)-XB(IP))
   50 CONTINUE
C
      IF(AREA.GE.0.0) THEN
       LCLOCK = .FALSE.
C       WRITE(*,1010) NB
      ELSE
C----- if area is negative (clockwise order), reverse coordinate order
       LCLOCK = .TRUE.
C       WRITE(*,1011) NB
       DO 55 I=1, NB/2
         XTMP = XB(NB-I+1)
         YTMP = YB(NB-I+1)
         XB(NB-I+1) = XB(I)
         YB(NB-I+1) = YB(I)
         XB(I) = XTMP
         YB(I) = YTMP
   55  CONTINUE
      ENDIF
C
      IF(LNORM) THEN
       CALL NORM(XB,XBP,YB,YBP,SB,NB)
C       WRITE(*,1020)
      ENDIF
C
      CALL SCALC(XB,YB,SB,NB)
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
C
      CALL GEOPAR(XB,XBP,YB,YBP,SB,NB, W1,
     &            SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
     &            EI11BA,EI22BA,APX1BA,APX2BA,
     &            EI11BT,EI22BT,APX1BT,APX2BT,
     &            THICKB,CAMBRB )
C
      XBLE = SEVAL(SBLE,XB,XBP,SB,NB)
      YBLE = SEVAL(SBLE,YB,YBP,SB,NB)
      XBTE = 0.5*(XB(1) + XB(NB))
      YBTE = 0.5*(YB(1) + YB(NB))
C
C      WRITE(*,1050) XBLE,YBLE, CHORDB,
C     &              XBTE,YBTE
C
C---- set MSES domain parameters
      XINL = XBLE - 2.0*CHORDB
      XOUT = XBLE + 3.0*CHORDB
      YBOT = YBLE - 2.5*CHORDB
      YTOP = YBLE + 3.5*CHORDB
      XINL = AINT(20.0*ABS(XINL/CHORDB)+0.5)/20.0 * SIGN(CHORDB,XINL)
      XOUT = AINT(20.0*ABS(XOUT/CHORDB)+0.5)/20.0 * SIGN(CHORDB,XOUT)
      YBOT = AINT(20.0*ABS(YBOT/CHORDB)+0.5)/20.0 * SIGN(CHORDB,YBOT)
      YTOP = AINT(20.0*ABS(YTOP/CHORDB)+0.5)/20.0 * SIGN(CHORDB,YTOP)
      WRITE(ISPARS,1005) XINL, XOUT, YBOT, YTOP
 1005 FORMAT(1X, 4F8.2)
C
C---- wipe out old flap hinge location
      XBF = 0.0
      YBF = 0.0
      LBFLAP = .FALSE.
C
C---- panel, copy to current airfoil, and check panel corner angles
      CALL PANGEN(.TRUE.)
      CALL ABCOPY(.TRUE.)
      CALL CANG(X,Y,N,0, IMAX,AMAX)
      IF(ABS(AMAX).GT.ANGTOL .AND. .NOT.LQUIET) THEN
       WRITE(*,1081) AMAX, IMAX
      ENDIF
C
      RETURN
C
 1010 FORMAT(/' Number of input coordinate points:', I4
     &       /' Counterclockwise ordering')
 1011 FORMAT(/' Number of input coordinate points:', I4
     &       /' Clockwise ordering')
 1020 FORMAT(/' Airfoil has been normalized')
 1050 FORMAT(/'  LE  x,y  =', 2F10.5,'  |   Chord =',F10.5
     &       /'  TE  x,y  =', 2F10.5,'  |'                 )
 1081 FORMAT(
     &  /' WARNING: Poor input coordinate distribution'
     &  /'          Excessive panel angle', F7.1,'  at i =', I4
     &  /'          Repaneling with PANE and/or PPAR suggested'
     &  /'           (doing GDES,CADD before repaneling _may_'
     &  /'            improve excessively coarse LE spacing' )
C
      END ! ADDGEOM
 
      SUBROUTINE NACA(IDES1)
      INCLUDE 'XFOIL.INC'
C
C---- number of points per side
      NSIDE = IQX/3
C
      IF(IDES1 .LE. 0) THEN
       IF(.NOT.LQUIET) WRITE(*,*) 'NACA: No designation supplied.'
       RETURN
      ELSE
       IDES = IDES1
      ENDIF
C
      ITYPE = 0
      IF(IDES.LE.25099) ITYPE = 5
      IF(IDES.LE.9999 ) ITYPE = 4
C
      IF(ITYPE.EQ.0) THEN
       IF(.NOT.LQUIET) WRITE(*,*) 'This designation not implemented.'
       RETURN
      ENDIF
C
      IF(ITYPE.EQ.4) CALL NACA4(IDES,W1,W2,W3,NSIDE,XB,YB,NB,NAME)
      IF(ITYPE.EQ.5) CALL NACA5(IDES,W1,W2,W3,NSIDE,XB,YB,NB,NAME)
      CALL STRIP(NAME,NNAME)
C
C---- see if routines didn't recognize designator
      IF(IDES.EQ.0) RETURN
C
      LCLOCK = .FALSE.
C
      XBF = 0.0
      YBF = 0.0
      LBFLAP = .FALSE.
C
      CALL SCALC(XB,YB,SB,NB)
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
C
      CALL GEOPAR(XB,XBP,YB,YBP,SB,NB, W1,
     &            SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
     &            EI11BA,EI22BA,APX1BA,APX2BA,
     &            EI11BT,EI22BT,APX1BT,APX2BT,
     &            THICKB,CAMBRB )
C
C      WRITE(*,1200) NB
 1200 FORMAT(/' Buffer airfoil set using', I4,' points')
C
C---- set paneling
      CALL PANGEN(.TRUE.)
ccc      CALL PANPLT
C
      RETURN
      END ! NACA


      SUBROUTINE PANGEN(SHOPAR)
C---------------------------------------------------
C     Set paneling distribution from buffer airfoil
C     geometry, thus creating current airfoil.
C 
C     If REFINE=True, bunch points at x=XSREF on
C     top side and at x=XPREF on bottom side
C     by setting a fictitious local curvature of
C     CTRRAT*(LE curvature) there.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      LOGICAL SHOPAR
C
      IF(NB.LT.2) THEN
       IF(.NOT.LQUIET) WRITE(*,*)'PANGEN: Buffer airfoil not available.'
       N = 0
       RETURN
      ENDIF
C
C---- Number of temporary nodes for panel distribution calculation
C       exceeds the specified panel number by factor of IPFAC.
      IPFAC = 3
C
C---- number of airfoil panel points
      N = NPAN
C
cC---- number of wake points
c      NW = NPAN/8 + 2
c      IF(NW.GT.IWX) THEN
c       WRITE(*,*)
c     &  'Array size (IWX) too small.  Last wake point index reduced.'
c       NW = IWX
c      ENDIF
C
C---- set arc length spline parameter
      CALL SCALC(XB,YB,SB,NB)
C
C---- spline raw airfoil coordinates
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
C
C---- normalizing length (~ chord)
      SBREF = 0.5*(SB(NB)-SB(1))
C
C---- set up curvature array
      DO I = 1, NB
        W5(I) = ABS( CURV(SB(I),XB,XBP,YB,YBP,SB,NB) ) * SBREF
      ENDDO
C
C---- locate LE point arc length value and the normalized curvature there
      CALL LEFIND(SBLE,XB,XBP,YB,YBP,SB,NB)
      CVLE = ABS( CURV(SBLE,XB,XBP,YB,YBP,SB,NB) ) * SBREF
C
C---- check for doubled point (sharp corner) at LE
      IBLE = 0
      DO I = 1, NB-1
        IF(SBLE.EQ.SB(I) .AND. SBLE.EQ.SB(I+1)) THEN
         IBLE = I
C         WRITE(*,*)
C         WRITE(*,*) 'Sharp leading edge'
         GO TO 21
        ENDIF
      ENDDO
 21   CONTINUE
C
C---- set LE, TE points
      XBLE = SEVAL(SBLE,XB,XBP,SB,NB)
      YBLE = SEVAL(SBLE,YB,YBP,SB,NB)
      XBTE = 0.5*(XB(1)+XB(NB))
      YBTE = 0.5*(YB(1)+YB(NB))
      CHBSQ = (XBTE-XBLE)**2 + (YBTE-YBLE)**2
C
C---- set average curvature over 2*NK+1 points within Rcurv of LE point
      NK = 3
      CVSUM = 0.
      DO K = -NK, NK
        FRAC = FLOAT(K)/FLOAT(NK)
        SBK = SBLE + FRAC*SBREF/MAX(CVLE,20.0)
        CVK = ABS( CURV(SBK,XB,XBP,YB,YBP,SB,NB) ) * SBREF
        CVSUM = CVSUM + CVK
      ENDDO
      CVAVG = CVSUM/FLOAT(2*NK+1)
C
C---- dummy curvature for sharp LE
      IF(IBLE.NE.0) CVAVG = 10.0
C
C---- set curvature attraction coefficient actually used
      CC = 6.0 * CVPAR
C
C---- set artificial curvature at TE to bunch panels there
      CVTE = CVAVG * CTERAT
      W5(1)  = CVTE
      W5(NB) = CVTE
C
C
C**** smooth curvature array for smoother panel size distribution  ****
C
CCC      CALL ASKR('Enter curvature smoothing length/c^',SMOOL)
CCC      SMOOL = 0.010
C
C---- set smoothing length = 1 / averaged LE curvature, but 
C-    no more than 5% of chord and no less than 1/4 average panel spacing
      SMOOL = MAX( 1.0/MAX(CVAVG,20.0) , 0.25 /FLOAT(NPAN/2) )
C
      SMOOSQ = (SMOOL*SBREF) ** 2
C
C---- set up tri-diagonal system for smoothed curvatures
      W2(1) = 1.0
      W3(1) = 0.0
      DO I=2, NB-1
        DSM = SB(I) - SB(I-1)
        DSP = SB(I+1) - SB(I)
        DSO = 0.5*(SB(I+1) - SB(I-1))
C
        IF(DSM.EQ.0.0 .OR. DSP.EQ.0.0) THEN
C------- leave curvature at corner point unchanged
         W1(I) = 0.0
         W2(I) = 1.0
         W3(I) = 0.0
        ELSE
         W1(I) =  SMOOSQ * (         - 1.0/DSM) / DSO
         W2(I) =  SMOOSQ * ( 1.0/DSP + 1.0/DSM) / DSO  +  1.0
         W3(I) =  SMOOSQ * (-1.0/DSP          ) / DSO
        ENDIF
      ENDDO
C
      W1(NB) = 0.0
      W2(NB) = 1.0
C
C---- fix curvature at LE point by modifying equations adjacent to LE
      DO I=2, NB-1
        IF(SB(I).EQ.SBLE .OR. I.EQ.IBLE .OR. I.EQ.IBLE+1) THEN
C------- if node falls right on LE point, fix curvature there
         W1(I) = 0.
         W2(I) = 1.0
         W3(I) = 0.
         W5(I) = CVLE
        ELSE IF(SB(I-1).LT.SBLE .AND. SB(I).GT.SBLE) THEN
C------- modify equation at node just before LE point
         DSM = SB(I-1) - SB(I-2)
         DSP = SBLE    - SB(I-1)
         DSO = 0.5*(SBLE - SB(I-2))
C
         W1(I-1) =  SMOOSQ * (         - 1.0/DSM) / DSO
         W2(I-1) =  SMOOSQ * ( 1.0/DSP + 1.0/DSM) / DSO  +  1.0
         W3(I-1) =  0.
         W5(I-1) = W5(I-1) + SMOOSQ*CVLE/(DSP*DSO)
C
C------- modify equation at node just after LE point
         DSM = SB(I) - SBLE
         DSP = SB(I+1) - SB(I)
         DSO = 0.5*(SB(I+1) - SBLE)
         W1(I) =  0.
         W2(I) =  SMOOSQ * ( 1.0/DSP + 1.0/DSM) / DSO  +  1.0
         W3(I) =  SMOOSQ * (-1.0/DSP          ) / DSO
         W5(I) = W5(I) + SMOOSQ*CVLE/(DSM*DSO)
C
         GO TO 51
        ENDIF
      ENDDO
   51 CONTINUE
C
C---- set artificial curvature at bunching points and fix it there
      DO I=2, NB-1
C------ chord-based x/c coordinate
        XOC = (  (XB(I)-XBLE)*(XBTE-XBLE)
     &         + (YB(I)-YBLE)*(YBTE-YBLE) ) / CHBSQ
C
        IF(SB(I).LT.SBLE) THEN
C------- check if top side point is in refinement area
         IF(XOC.GT.XSREF1 .AND. XOC.LT.XSREF2) THEN
          W1(I) = 0.
          W2(I) = 1.0
          W3(I) = 0.
          W5(I) = CVLE*CTRRAT
         ENDIF
        ELSE
C------- check if bottom side point is in refinement area
         IF(XOC.GT.XPREF1 .AND. XOC.LT.XPREF2) THEN
          W1(I) = 0.
          W2(I) = 1.0
          W3(I) = 0.
          W5(I) = CVLE*CTRRAT
         ENDIF
        ENDIF
      ENDDO
C
C---- solve for smoothed curvature array W5
      IF(IBLE.EQ.0) THEN
       CALL TRISOL(W2,W1,W3,W5,NB)
      ELSE
       I = 1
       CALL TRISOL(W2(I),W1(I),W3(I),W5(I),IBLE)
       I = IBLE+1
       CALL TRISOL(W2(I),W1(I),W3(I),W5(I),NB-IBLE)
      ENDIF
C
C---- find max curvature
      CVMAX = 0.
      DO I=1, NB
        CVMAX = MAX( CVMAX , ABS(W5(I)) )
      ENDDO
C
C---- normalize curvature array
      DO I=1, NB
        W5(I) = W5(I) / CVMAX
      ENDDO
C
C---- spline curvature array
      CALL SEGSPL(W5,W6,SB,NB)
C
C---- Set initial guess for node positions uniform in s.
C     More nodes than specified (by factor of IPFAC) are 
C     temporarily used  for more reliable convergence.
      NN = IPFAC*(N-1)+1
C
C---- ratio of lengths of panel at TE to one away from the TE
      RDSTE = 0.667
      RTF = (RDSTE-1.0)*FLOAT(IPFAC) + 1.0
C
      IF(IBLE.EQ.0) THEN
C
       DSAVG = (SB(NB)-SB(1))/(FLOAT(NN-3) + 2.0*RTF)
       SNEW(1) = SB(1)
       DO I=2, NN-1
         SNEW(I) = SB(1) + DSAVG * (FLOAT(I-2) + RTF)
       ENDDO
       SNEW(NN) = SB(NB)
C
      ELSE
C
       NFRAC1 = (N * IBLE) / NB
C
       NN1 = IPFAC*(NFRAC1-1)+1
       DSAVG1 = (SBLE-SB(1))/(FLOAT(NN1-2) + RTF)
       SNEW(1) = SB(1)
       DO I=2, NN1
         SNEW(I) = SB(1) + DSAVG1 * (FLOAT(I-2) + RTF)
       ENDDO
C
       NN2 = NN - NN1 + 1
       DSAVG2 = (SB(NB)-SBLE)/(FLOAT(NN2-2) + RTF)
       DO I=2, NN2-1
         SNEW(I-1+NN1) = SBLE + DSAVG2 * (FLOAT(I-2) + RTF)
       ENDDO
       SNEW(NN) = SB(NB)
C
      ENDIF
C
C---- Newton iteration loop for new node positions
      DO 10 ITER=1, 20
C
C------ set up tri-diagonal system for node position deltas
        CV1  = SEVAL(SNEW(1),W5,W6,SB,NB)
        CV2  = SEVAL(SNEW(2),W5,W6,SB,NB)
        CVS1 = DEVAL(SNEW(1),W5,W6,SB,NB)
        CVS2 = DEVAL(SNEW(2),W5,W6,SB,NB)
C
        CAVM = SQRT(CV1**2 + CV2**2)
        IF(CAVM .EQ. 0.0) THEN
          CAVM_S1 = 0.
          CAVM_S2 = 0.
        ELSE
          CAVM_S1 = CVS1 * CV1/CAVM
          CAVM_S2 = CVS2 * CV2/CAVM
        ENDIF
C
        DO 110 I=2, NN-1
          DSM = SNEW(I) - SNEW(I-1)
          DSP = SNEW(I) - SNEW(I+1)
          CV3  = SEVAL(SNEW(I+1),W5,W6,SB,NB)
          CVS3 = DEVAL(SNEW(I+1),W5,W6,SB,NB)
C
          CAVP = SQRT(CV3**2 + CV2**2)
          IF(CAVP .EQ. 0.0) THEN
            CAVP_S2 = 0.
            CAVP_S3 = 0.
          ELSE
            CAVP_S2 = CVS2 * CV2/CAVP
            CAVP_S3 = CVS3 * CV3/CAVP
          ENDIF
C
          FM = CC*CAVM + 1.0
          FP = CC*CAVP + 1.0
C
          REZ = DSP*FP + DSM*FM
C
C-------- lower, main, and upper diagonals
          W1(I) =      -FM  +  CC*               DSM*CAVM_S1
          W2(I) =  FP + FM  +  CC*(DSP*CAVP_S2 + DSM*CAVM_S2)
          W3(I) = -FP       +  CC* DSP*CAVP_S3
C
C-------- residual, requiring that
C         (1 + C*curv)*deltaS is equal on both sides of node i
          W4(I) = -REZ
C
          CV1 = CV2
          CV2 = CV3
          CVS1 = CVS2
          CVS2 = CVS3
          CAVM    = CAVP
          CAVM_S1 = CAVP_S2
          CAVM_S2 = CAVP_S3
  110   CONTINUE
C
C------ fix endpoints (at TE)
        W2(1) = 1.0
        W3(1) = 0.0
        W4(1) = 0.0
        W1(NN) = 0.0
        W2(NN) = 1.0
        W4(NN) = 0.0
C
        IF(RTF .NE. 1.0) THEN
C------- fudge equations adjacent to TE to get TE panel length ratio RTF
C
         I = 2
         W4(I) = -((SNEW(I) - SNEW(I-1)) + RTF*(SNEW(I) - SNEW(I+1)))
         W1(I) = -1.0
         W2(I) =  1.0 + RTF
         W3(I) =      - RTF
C
         I = NN-1
         W4(I) = -((SNEW(I) - SNEW(I+1)) + RTF*(SNEW(I) - SNEW(I-1)))
         W3(I) = -1.0
         W2(I) =  1.0 + RTF
         W1(I) =      - RTF
        ENDIF
C
C
C------ fix sharp LE point
        IF(IBLE.NE.0) THEN
         I = NN1
         W1(I) = 0.0
         W2(I) = 1.0
         W3(I) = 0.0
         W4(I) = SBLE - SNEW(I)
        ENDIF
C
C------ solve for changes W4 in node position arc length values
        CALL TRISOL(W2,W1,W3,W4,NN)
C
C------ find under-relaxation factor to keep nodes from changing order
        RLX = 1.0
        DMAX = 0.0
        DO I=1, NN-1
          DS  = SNEW(I+1) - SNEW(I)
          DDS = W4(I+1) - W4(I)
          DSRAT = 1.0 + RLX*DDS/DS
          IF(DSRAT.GT.4.0) RLX = (4.0-1.0)*DS/DDS
          IF(DSRAT.LT.0.2) RLX = (0.2-1.0)*DS/DDS
          DMAX = MAX(ABS(W4(I)),DMAX)
        ENDDO
C
C------ update node position
        DO I=2, NN-1
          SNEW(I) = SNEW(I) + RLX*W4(I)
        ENDDO
C
CCC        IF(RLX.EQ.1.0) WRITE(*,*) DMAX
CCC        IF(RLX.NE.1.0) WRITE(*,*) DMAX,'    RLX =',RLX
        IF(ABS(DMAX).LT.1.E-3) GO TO 11
   10 CONTINUE
      IF(.NOT.LQUIET)
     & WRITE(*,*) 'Paneling convergence failed.  Continuing anyway...'
C
   11 CONTINUE
C
C---- set new panel node coordinates
      DO I=1, N
        IND = IPFAC*(I-1) + 1
        S(I) = SNEW(IND)
        X(I) = SEVAL(SNEW(IND),XB,XBP,SB,NB)
        Y(I) = SEVAL(SNEW(IND),YB,YBP,SB,NB)
      ENDDO
C
C
C---- go over buffer airfoil again, checking for corners (double points)
      NCORN = 0
      DO 25 IB=1, NB-1
        IF(SB(IB) .EQ. SB(IB+1)) THEN
C------- found one !
C
         NCORN = NCORN+1
         XBCORN = XB(IB)
         YBCORN = YB(IB)
         SBCORN = SB(IB)
C
C------- find current-airfoil panel which contains corner
         DO 252 I=1, N
C
C--------- keep stepping until first node past corner
           IF(S(I) .LE. SBCORN) GO TO 252
C
C---------- move remainder of panel nodes to make room for additional node
            DO 2522 J=N, I, -1
              X(J+1) = X(J)
              Y(J+1) = Y(J)
              S(J+1) = S(J)
 2522       CONTINUE
            N = N+1
C
            IF(N .GT. IQX-1)
     &       STOP 'PANEL: Too many panels. Increase IQX in XFOIL.INC'
C
            X(I) = XBCORN
            Y(I) = YBCORN
            S(I) = SBCORN
C
C---------- shift nodes adjacent to corner to keep panel sizes comparable
            IF(I-2 .GE. 1) THEN
             S(I-1) = 0.5*(S(I) + S(I-2))
             X(I-1) = SEVAL(S(I-1),XB,XBP,SB,NB)
             Y(I-1) = SEVAL(S(I-1),YB,YBP,SB,NB)
            ENDIF
C
            IF(I+2 .LE. N) THEN
             S(I+1) = 0.5*(S(I) + S(I+2))
             X(I+1) = SEVAL(S(I+1),XB,XBP,SB,NB)
             Y(I+1) = SEVAL(S(I+1),YB,YBP,SB,NB)
            ENDIF
C
C---------- go on to next input geometry point to check for corner
            GO TO 25
C
  252    CONTINUE
        ENDIF
   25 CONTINUE
C
      CALL SCALC(X,Y,S,N)
      CALL SEGSPL(X,XP,S,N)
      CALL SEGSPL(Y,YP,S,N)
      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
C
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
      CHORD  = SQRT( (XTE-XLE)**2 + (YTE-YLE)**2 )
C
C---- calculate panel size ratios (user info)
      DSMIN =  1000.0
      DSMAX = -1000.0
      DO 40 I=1, N-1
        DS = S(I+1)-S(I)
        IF(DS .EQ. 0.0) GO TO 40
          DSMIN = MIN(DSMIN,DS)
          DSMAX = MAX(DSMAX,DS)
   40 CONTINUE
C
      DSMIN = DSMIN*FLOAT(N-1)/S(N)
      DSMAX = DSMAX*FLOAT(N-1)/S(N)
ccc      WRITE(*,*) 'DSmin/DSavg = ',DSMIN,'     DSmax/DSavg = ',DSMAX
C
C---- set various flags for new airfoil
      LGAMU = .FALSE.
      LQINU = .FALSE.
      LWAKE = .FALSE.
      LQAIJ = .FALSE.
      LADIJ = .FALSE.
      LWDIJ = .FALSE.
      LIPAN = .FALSE.
      LBLINI = .FALSE.
      LVCONV = .FALSE.
      LSCINI = .FALSE.
      LQSPEC = .FALSE.
      LGSAME = .FALSE.
C
      IF(LBFLAP) THEN
       XOF = XBF
       YOF = YBF
       LFLAP = .TRUE.
      ENDIF
C
C---- determine if TE is blunt or sharp, calculate TE geometry parameters
      CALL TECALC
C
C---- calculate normal vectors
      CALL NCALC(X,Y,S,N,NX,NY)
C
C---- calculate panel angles for panel routines
      CALL APCALC
C
      IF(.NOT.SHARP) THEN
       GAP = SQRT((X(1)-X(N))**2 + (Y(1)-Y(N))**2)
      ENDIF
 1090 FORMAT(/1X,A,F9.5)
C
C      IF(SHOPAR) WRITE(*,1100) NPAN, CVPAR, CTERAT, CTRRAT,
C     &                         XSREF1, XSREF2, XPREF1, XPREF2
 1100 FORMAT(/' Paneling parameters used...'
     &       /'   Number of panel nodes      ' , I4
     &       /'   Panel bunching parameter   ' , F6.3
     &       /'   TE/LE panel density ratio  ' , F6.3
     &       /'   Refined-area/LE panel density ratio   ' , F6.3
     &       /'   Top    side refined area x/c limits ' , 2F6.3
     &       /'   Bottom side refined area x/c limits ' , 2F6.3)
C
      RETURN
      END ! PANGEN

      SUBROUTINE SETPAN(NPAN_IN, CVPAR_IN, CTERAT_IN, CTRRAT_IN,
     &                  XSREF1_IN, XSREF2_IN, XPREF1_IN, XPREF2_IN)
      INCLUDE 'XFOIL.INC'
C
      REAL CVPAR_IN, CTERAT_IN, CTRRAT_IN
      REAL XSREF1_IN, XSREF2_IN, XPREF1_IN, XPREF2_IN
      INTEGER NPAN_IN
C
      IF(NB.LE.1) THEN
       IF(.NOT.LQUIET) WRITE(*,*)'SETPAN: Buffer airfoil not available.'
       RETURN
      ENDIF
C
C---- copy input arguments into global common variables
      NPAN   = NPAN_IN
      CVPAR  = CVPAR_IN
      CTERAT = CTERAT_IN
      CTRRAT = CTRRAT_IN
      XSREF1 = XSREF1_IN
      XSREF2 = XSREF2_IN
      XPREF1 = XPREF1_IN
      XPREF2 = XPREF2_IN
C
C---- clamp panel count to array limit
      IF(NPAN .GT. IQX-6) THEN
       NPAN = IQX - 6
       IF(.NOT.LQUIET) WRITE(*,1000) NPAN
 1000  FORMAT(1X,'SETPAN: Panel nodes reduced to array limit:',I4)
      ENDIF
C
C---- generate panel distribution and report max corner angle
      CALL PANGEN(.FALSE.)
      IF(N.GT.0) CALL CANG(X,Y,N,1,IMAX,AMAX)
C
      RETURN
      END ! SETPAN

      SUBROUTINE TECALC
C-------------------------------------------
C     Calculates total and projected TE gap 
C     areas and TE panel strengths.
C-------------------------------------------
      INCLUDE 'XFOIL.INC'
C
C---- set TE base vector and TE bisector components
      DXTE = X(1) - X(N)
      DYTE = Y(1) - Y(N)
      DXS = 0.5*(-XP(1) + XP(N))
      DYS = 0.5*(-YP(1) + YP(N))
C
C---- normal and streamwise projected TE gap areas
      ANTE = DXS*DYTE - DYS*DXTE
      ASTE = DXS*DXTE + DYS*DYTE
C
C---- total TE gap area
      DSTE = SQRT(DXTE**2 + DYTE**2)
C
      SHARP = DSTE .LT. 0.0001*CHORD
C
      IF(SHARP) THEN
       SCS = 1.0
       SDS = 0.0
      ELSE
       SCS = ANTE/DSTE
       SDS = ASTE/DSTE
      ENDIF
C
C---- TE panel source and vorticity strengths
      SIGTE = 0.5*(GAM(1) - GAM(N))*SCS
      GAMTE = -.5*(GAM(1) - GAM(N))*SDS
C
      SIGTE_A = 0.5*(GAM_A(1) - GAM_A(N))*SCS
      GAMTE_A = -.5*(GAM_A(1) - GAM_A(N))*SDS
C
      RETURN
      END ! TECALC
