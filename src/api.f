C***********************************************************************
C    Module:  api.f
C
C    Copyright (C) 2026
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
C---- API subroutines for f2py Python interface
C
C=======================================================================
      SUBROUTINE INIT()
C---------------------------------------
C     Initialises all XFOIL variables.
C     Must be called before any other API
C     subroutine.
C---------------------------------------
      CALL INITIALISE
      END SUBROUTINE INIT

C=======================================================================
      SUBROUTINE Quiet(setting)
C---------------------------------------------------
C     Enables/disables solver iteration output.
C     setting = .TRUE.  suppresses WRITE output
C     setting = .FALSE. allows WRITE output
C---------------------------------------------------
      INCLUDE 'QUIET.INC'
      LOGICAL setting
Cf2py intent(in) :: setting
C
      LQUIET = setting
C
      END SUBROUTINE Quiet

      SUBROUTINE SETGEOM(X_in, Y_in, N_in)
      INCLUDE 'XFOIL.INC'
      INCLUDE 'QUIET.INC'
C------------------------------------------------------
C     Loads airfoil from passed-in coordinate arrays
C     into the buffer airfoil and does initial
C     processing (orientation check, normalisation,
C     spline, geometry parameters, paneling).
C------------------------------------------------------
      INTEGER N_in
      REAL*8 X_in(N_in), Y_in(N_in)
Cf2py intent(in) :: X_in, Y_in, N_in
Cf2py integer intent(hide),depend(X_in) :: N_in = len(X_in)
C
      DATA ANGTOL / 40.0 /
C
      IF(N_in.LE.1) THEN
       IF(.NOT.LQUIET) WRITE(*,*) 'No airfoil data supplied.'
       RETURN
      ELSEIF(N_in.GT.IBX) THEN
       IF(.NOT.LQUIET) WRITE(*,*) 'Too many points supplied.'
       RETURN
      ELSE
       IF(.NOT.LQUIET) WRITE(*,*) 'Airfoil set using', N_in, 'points.'
      ENDIF
C
C---- copy incoming coordinates into global buffer airfoil arrays
      NB = N_in
      DO 10 I=1, NB
        XB(I) = X_in(I)
        YB(I) = Y_in(I)
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
      IF(.NOT.LQUIET) WRITE(ISPARS,1005) XINL, XOUT, YBOT, YTOP
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
      END SUBROUTINE SETGEOM
      
C=======================================================================
      SUBROUTINE SETNACA(desig)
C---------------------------------------------------
C     Generates and panels a NACA 4 or 5 digit
C     airfoil from a character designation string.
C     e.g. desig = "0012" -> NACA 0012
C          desig = "23015" -> NACA 23015
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INCLUDE 'QUIET.INC'
      CHARACTER*(10) desig
Cf2py character*10 desig
Cf2py intent(in) :: desig
C
      INTEGER IDES
C
C---- read integer designation from string
C     "0012" -> 12, "2412" -> 2412, "23015" -> 23015
      READ(desig, *) IDES
C
C
C---- number of points per side
      NSIDE = IQX/3
C
      IF(IDES .LE. 0) THEN
       IF(.NOT.LQUIET) WRITE(*,*) 'NACA: No designation supplied.'
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
      IF(.NOT.LQUIET) WRITE(*,1200) NB
 1200 FORMAT(/' Buffer airfoil set using', I4,' points')
C
C---- set paneling
      CALL PANGEN(.TRUE.)
ccc      CALL PANPLT
C
      RETURN
C
      END SUBROUTINE SETNACA
C
C=======================================================================
      SUBROUTINE SETCST(LEM, TE, Au, NAu, Al, NAl, N1, N2)
C---------------------------------------------------
C     CST airfoil parameterisation entry point.
C     Builds a CST airfoil, loads the buffer arrays
C     XB/YB/NB, calls PANGEN, and resets solution flags.
C
C     Python signature (NAu and NAl hidden by f2py):
C         xf.setcst(LEM, TE, Au, Al, N1, N2)
C---------------------------------------------------
      INTEGER NAu, NAl
      REAL*8 Au(NAu), Al(NAl)
      REAL*8 LEM, TE, N1, N2
Cf2py intent(in) :: LEM, TE
Cf2py intent(in) :: Au
Cf2py integer intent(hide),depend(Au) :: NAu = len(Au)
Cf2py intent(in) :: Al
Cf2py integer intent(hide),depend(Al) :: NAl = len(Al)
Cf2py intent(in) :: N1, N2
C
      CALL CST(LEM, TE, Au, NAu, Al, NAl, N1, N2)
C
      END SUBROUTINE SETCST
C

C=======================================================================
      SUBROUTINE SETPARSEC(rle_up, Xup, Zup, Zxxup,
     &                     rle_lo, Xlo, Zlo, Zxxlo,
     &                     Zte, dZte, alpha_te_deg, beta_te_deg)
C---------------------------------------------------
C     Generates a PARSEC aerofoil from twelve geometric parameters,
C     loads it into the XFOIL buffer arrays, and panels it.
C     Wrapper around PARSEC in parsec.f.
C
C     Input: (all lengths normalised to unit chord)
C       rle_up       upper-surface geometric LE radius
C       Xup          upper crest abscissa  (0 < Xup < 1)
C       Zup          upper crest ordinate
C       Zxxup        upper crest curvature (d^2z/dx^2, typically < 0)
C       rle_lo       lower-surface geometric LE radius
C       Xlo          lower crest abscissa  (0 < Xlo < 1)
C       Zlo          lower crest ordinate
C       Zxxlo        lower crest curvature (typically > 0)
C       Zte          TE mid-ordinate   = 0.5*(z_upper + z_lower) at x=1
C       dZte         TE thickness gap  = z_upper - z_lower at x=1 (>= 0)
C       alpha_te_deg TE direction angle (degrees)
C       beta_te_deg  TE wedge angle     (degrees)
C---------------------------------------------------
      REAL*8 rle_up, Xup, Zup, Zxxup
      REAL*8 rle_lo, Xlo, Zlo, Zxxlo
      REAL*8 Zte, dZte, alpha_te_deg, beta_te_deg
Cf2py intent(in) :: rle_up, Xup, Zup, Zxxup
Cf2py intent(in) :: rle_lo, Xlo, Zlo, Zxxlo
Cf2py intent(in) :: Zte, dZte, alpha_te_deg, beta_te_deg
C
      CALL PARSEC(rle_up, Xup, Zup, Zxxup,
     &            rle_lo, Xlo, Zlo, Zxxlo,
     &            Zte, dZte, alpha_te_deg, beta_te_deg)
C
      END SUBROUTINE SETPARSEC
C

C=======================================================================
      SUBROUTINE FLOWCONS(Re, Ma)
C---------------------------------------------------
C     Sets Reynolds and Mach numbers.
C     Setting Re > 0 automatically enables viscous
C     mode (LVISC = .TRUE.).
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 Re, Ma
Cf2py intent(in) :: Re, Ma
C
      REAL*8 MINF_CLM, REINF_CLM
C
      REINF1 = Re
      MINF1  = Ma
      RETYP  = 1
      MATYP  = 1
      LVISC  = (Re .GT. 0.0)
      LVCONV = .FALSE.
C
      CALL MRCL(1.0, MINF_CLM, REINF_CLM)
      CALL COMSET
C
      END SUBROUTINE FLOWCONS
C
C=======================================================================
      SUBROUTINE SETCCOR(icor)
C---------------------------------------------------
C     Selects compressibility correction method.
C     icor = 0  Prandtl-Glauert
C     icor = 1  Karman-Tsien  (default)
C     icor = 2  Laitone
C---------------------------------------------------
      INCLUDE 'CPCOR.INC'
      INTEGER icor
Cf2py intent(in) :: icor
      ICCOR = icor
      CALL COMSET
      END SUBROUTINE SETCCOR
C
C=======================================================================
      SUBROUTINE MAXITER(Iter)
C---------------------------------------------------
C     Sets the maximum number of viscous iterations.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INTEGER Iter
Cf2py intent(in) :: Iter
C
      ITMAX = Iter
C
      END SUBROUTINE MAXITER
C
C=======================================================================
      SUBROUTINE SETBLRLX(DHI_in, DLO_in, DNUE_in)
C---------------------------------------------------
C     Sets BL solver relaxation parameters.
C     DHI_in:  high limit for normalized change (default 1.5)
C     DLO_in:  low limit for normalized change (default -0.5)
C     DNUE_in: Ue normalization divisor (default 0.25)
C
C     Tighter values (e.g. 0.8, -0.3, 0.10) reduce
C     sensitivity to platform math library differences
C     at the cost of more iterations.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 DHI_in, DLO_in, DNUE_in
Cf2py intent(in) :: DHI_in, DLO_in, DNUE_in
C
      RLXPAR_DHI  = DHI_in
      RLXPAR_DLO  = DLO_in
      RLXPAR_DNUE = DNUE_in
C
      END SUBROUTINE SETBLRLX
C
C=======================================================================
      SUBROUTINE SETTRIP(Xtrip1, Xtrip2)
C---------------------------------------------------
C     Sets transition parameters.
C     Xtrip1: upper surface trip x/c location
C     Xtrip2: lower surface trip x/c location
C     Ncrit:  log of critical amplification ratio
C
C     Xtrip > 1.0 means no forced transition.
C     Negative Xtrip = -s/s_side arc length fraction.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 Xtrip1, Xtrip2
Cf2py intent(in) :: Xtrip1, Xtrip2
C
      XSTRIP(1) = Xtrip1
      XSTRIP(2) = Xtrip2
C
      END SUBROUTINE SETTRIP
C=======================================================================
      SUBROUTINE SETNCRIT(Ncrit)
C---------------------------------------------------
C     Sets transition parameters.
C     Xtrip1: upper surface trip x/c location
C     Xtrip2: lower surface trip x/c location
C     Ncrit:  log of critical amplification ratio
C
C     Xtrip > 1.0 means no forced transition.
C     Negative Xtrip = -s/s_side arc length fraction.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 Ncrit
Cf2py intent(in) :: Ncrit
C
      ACRIT(1)  = Ncrit
      ACRIT(2)  = Ncrit
C
      END SUBROUTINE SETNCRIT
C=======================================================================
      SUBROUTINE SETNCRIT12(Ncrit1, Ncrit2)
C---------------------------------------------------
C     Sets transition parameters.
C     Xtrip1: upper surface trip x/c location
C     Xtrip2: lower surface trip x/c location
C     Ncrit:  log of critical amplification ratio
C
C     Xtrip > 1.0 means no forced transition.
C     Negative Xtrip = -s/s_side arc length fraction.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 Ncrit1, Ncrit2
Cf2py intent(in) :: Ncrit1, Ncrit2
C
      ACRIT(1)  = Ncrit1
      ACRIT(2)  = Ncrit2
C
      END SUBROUTINE SETNCRIT12

C=======================================================================
      SUBROUTINE PANEL()
C---------------------------------------------------
C     Sets paneling parameters and re-panels the
C     airfoil if geometry has been loaded.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
C---- re-panel only if buffer airfoil exists
      IF(NB .GT. 0) THEN
        CALL PANGEN(.TRUE.)
C
C------ reset all dependent flags
        LGAMU  = .FALSE.
        LQAIJ  = .FALSE.
        LWAKE  = .FALSE.
        LBLINI = .FALSE.
        LIPAN  = .FALSE.
        LVCONV = .FALSE.
      ENDIF
C
      END SUBROUTINE PANEL

C=======================================================================
      SUBROUTINE REPANEL(NPan_in, CVPar_in, CTERat_in, CTRRat_in,
     &                   XSREF1_in, XSREF2_in, XPREF1_in, XPREF2_in)
C---------------------------------------------------
C     Sets paneling parameters and re-panels the
C     airfoil if geometry has been loaded.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INTEGER NPan_in
      REAL*8 CVPar_in, CTERat_in, CTRRat_in
      REAL*8 XSREF1_in, XSREF2_in, XPREF1_in, XPREF2_in
Cf2py intent(in) :: NPan_in, CVPar_in, CTERat_in, CTRRat_in
Cf2py intent(in) :: XSREF1_in, XSREF2_in, XPREF1_in, XPREF2_in
C
      NPAN   = NPan_in
      CVPAR  = CVPar_in
      CTERAT = CTERat_in
      CTRRAT = CTRRat_in
      XSREF1 = XSREF1_in
      XSREF2 = XSREF2_in
      XPREF1 = XPREF1_in
      XPREF2 = XPREF2_in
C---- re-panel only if buffer airfoil exists
      IF(NB .GT. 0) THEN
        CALL PANGEN(.TRUE.)
C
C------ reset all dependent flags
        LGAMU  = .FALSE.
        LQAIJ  = .FALSE.
        LWAKE  = .FALSE.
        LBLINI = .FALSE.
        LIPAN  = .FALSE.
        LVCONV = .FALSE.
      ENDIF
C
      END SUBROUTINE REPANEL

C=======================================================================
      SUBROUTINE RESETBL()
C---------------------------------------------------
C     Resets boundary layer initialisation flags.
C     Call between sequential analyses to force
C     fresh BL initialisation.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      LBLINI = .FALSE.
      LIPAN  = .FALSE.
C
      END SUBROUTINE RESETBL

C=======================================================================
      SUBROUTINE ALPHA(alpha_deg,
     &                 cl_out, cd_out, cm_out,
     &                 xtu_out, xtb_out, conv_out)
C---------------------------------------------------
C     Analyses airfoil at specified angle of attack.
C     Input:  alpha_deg  angle of attack (degrees)
C     Output: cl_out     lift coefficient
C             cd_out     drag coefficient
C             cm_out     moment coefficient
C             xtu_out    upper surface transition x/c
C             xtb_out    lower surface transition x/c
C             conv_out   .TRUE. if converged
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 alpha_deg
      REAL*8 cl_out, cd_out, cm_out, xtu_out, xtb_out
      LOGICAL conv_out
Cf2py intent(in)  :: alpha_deg
Cf2py intent(out) :: cl_out, cd_out, cm_out, xtu_out, xtb_out
Cf2py intent(out) :: conv_out
C
      ALFA  = alpha_deg * DTOR
      QINF  = 1.0
      LALFA = .TRUE.
C
      CALL SPECAL
C
      IF(ABS(ALFA - AWAKE) .GT. 1.0E-5) LWAKE  = .FALSE.
      IF(ABS(ALFA - AVISC) .GT. 1.0E-5) LVCONV = .FALSE.
      IF(ABS(MINF - MVISC) .GT. 1.0E-5) LVCONV = .FALSE.
C
      IF(LVISC) THEN
        CALL VISCAL(ITMAX)
        conv_out = LVCONV
      ELSE
        conv_out = .TRUE.
        XOCTR(1) = 0.0
        XOCTR(2) = 0.0
      ENDIF
C
      cl_out  = CL
      cd_out  = CD
      cm_out  = CM
      xtu_out = XOCTR(1)
      xtb_out = XOCTR(2)
C
      END SUBROUTINE ALPHA

C=======================================================================
      SUBROUTINE ASEQ(alpha_i, alpha_f, n_alpha,
     &                a_arr, cl_arr, cd_arr, cm_arr,
     &                xtu_arr, xtb_arr, conv_arr)
C---------------------------------------------------
C     Runs an angle-of-attack sweep.
C     Input:  alpha_i   start angle (degrees)
C             alpha_f   end angle (degrees)
C             n_alpha   number of points
C     Output: a_arr     angle array (degrees)
C             cl_arr    CL array
C             cd_arr    CD array
C             cm_arr    CM array
C             xtu_arr   upper transition x/c array
C             xtb_arr   lower transition x/c array
C             conv_arr  convergence flag array
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 alpha_i, alpha_f
      INTEGER n_alpha
      REAL*8 a_arr(n_alpha), cl_arr(n_alpha), cd_arr(n_alpha)
      REAL*8 cm_arr(n_alpha), xtu_arr(n_alpha), xtb_arr(n_alpha)
      LOGICAL conv_arr(n_alpha)
Cf2py intent(in)  :: alpha_i, alpha_f, n_alpha
Cf2py intent(out) :: a_arr, cl_arr, cd_arr, cm_arr
Cf2py intent(out) :: xtu_arr, xtb_arr, conv_arr
C
      INTEGER i, itmaxs
      REAL*8 a0, da
C
      a0 = alpha_i * DTOR
      IF(n_alpha .GT. 1) THEN
        da = (alpha_f - alpha_i) / FLOAT(n_alpha - 1) * DTOR
      ELSE
        da = 0.0
      ENDIF
C
      LALFA = .TRUE.
C
      DO i = 1, n_alpha
        ALFA = a0 + da * FLOAT(i - 1)
C
        IF(ABS(ALFA - AWAKE) .GT. 1.0E-5) LWAKE  = .FALSE.
        IF(ABS(ALFA - AVISC) .GT. 1.0E-5) LVCONV = .FALSE.
        IF(ABS(MINF - MVISC) .GT. 1.0E-5) LVCONV = .FALSE.
C
        CALL SPECAL
C
        itmaxs = ITMAX + 5
        IF(LVISC) THEN
          CALL VISCAL(itmaxs)
          conv_arr(i) = LVCONV
        ELSE
          conv_arr(i) = .TRUE.
          XOCTR(1) = 0.0
          XOCTR(2) = 0.0
        ENDIF
C
        a_arr(i)   = ALFA / DTOR
        cl_arr(i)  = CL
        cd_arr(i)  = CD
        cm_arr(i)  = CM
        xtu_arr(i) = XOCTR(1)
        xtb_arr(i) = XOCTR(2)
      ENDDO
C
      END SUBROUTINE ASEQ

C=======================================================================
      SUBROUTINE SCL(cl_val,
     &                   a_out, cl_out, cd_out, cm_out, 
     &                   xtu_out, xtb_out, conv_out)
C---------------------------------------------------
C     Analyses airfoil at specified lift coefficient.
C     Input:  cl_val    target CL
C     Output: a_out     solved angle of attack (degrees)
C             cl_out    actual CL achieved
C             cd_out    drag coefficient
C             cm_out    moment coefficient
C             xtu_out   upper transition x/c
C             xtb_out   lower transition x/c
C             conv_out  .TRUE. if converged
C
C     Note: AWAKE/AVISC/MVISC checks are done AFTER SPECCL
C     because SPECCL solves for ALFA -- it is unknown beforehand.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 cl_val
      REAL*8 a_out, cl_out, cd_out, cm_out, xtu_out, xtb_out
      LOGICAL conv_out
Cf2py intent(in)  :: cl_val
Cf2py intent(out) :: a_out, cl_out, cd_out, cm_out, xtu_out, xtb_out
Cf2py intent(out) :: conv_out
C
      CLSPEC = cl_val
      ALFA   = 0.0
      QINF   = 1.0
      LALFA  = .FALSE.
C
      CALL SPECCL
C
C---- invalidation guards applied after SPECCL has solved for ALFA
      IF(ABS(ALFA - AWAKE) .GT. 1.0E-5) LWAKE  = .FALSE.
      IF(ABS(ALFA - AVISC) .GT. 1.0E-5) LVCONV = .FALSE.
      IF(ABS(MINF - MVISC) .GT. 1.0E-5) LVCONV = .FALSE.
C
      IF(LVISC) THEN
        CALL VISCAL(ITMAX)
        conv_out = LVCONV
      ELSE
        conv_out = .TRUE.
        XOCTR(1) = 0.0
        XOCTR(2) = 0.0
      ENDIF
C
      a_out  = ALFA / DTOR
      cl_out = CL
      cd_out = CD
      cm_out = CM
      xtu_out = XOCTR(1)
      xtb_out = XOCTR(2)
C
      END SUBROUTINE SCL

C=======================================================================
      SUBROUTINE CSEQ(cl_start, cl_end, n_step,
     &                a_arr, cl_arr, cd_arr, cm_arr, xtu_arr, xtb_arr,
     &                conv_arr)
C---------------------------------------------------
C     Runs a lift-coefficient sweep.
C     Input:  cl_start   start CL
C             cl_end     end CL (exclusive: step size = range/n_step)
C             n_step     number of points
C     Output: a_arr      solved alpha array (degrees)
C             cl_arr     actual CL array
C             cd_arr     CD array
C             cm_arr     CM array
C             xtu_arr    upper transition x/c array
C             xtb_arr    lower transition x/c array
C             conv_arr   convergence flag array
C
C     Note: CL step is (cl_end-cl_start)/n_step so the endpoint
C     cl_end is not evaluated (analogous to Python range() semantics).
C     AWAKE/AVISC/MVISC checks are applied after each SPECCL call
C     because SPECCL solves for ALFA.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 cl_start, cl_end
      INTEGER n_step
      REAL*8 a_arr(n_step), cl_arr(n_step), cd_arr(n_step)
      REAL*8 cm_arr(n_step), xtu_arr(n_step), xtb_arr(n_step)
      LOGICAL conv_arr(n_step)
Cf2py intent(in)  :: cl_start, cl_end, n_step
Cf2py intent(out) :: a_arr, cl_arr, cd_arr, cm_arr, xtu_arr, xtb_arr, conv_arr
C
      INTEGER i, itmaxs
      REAL*8 cl0, dcl
C
      cl0  = cl_start
      dcl  = (cl_end - cl_start) / FLOAT(n_step)
C
      LALFA = .FALSE.
C
      DO i = 1, n_step
        CLSPEC = cl0 + dcl * FLOAT(i - 1)
        ALFA   = 0.0
        QINF   = 1.0
C
        CALL SPECCL
C
C------ invalidation guards after SPECCL has solved for ALFA
        IF(ABS(ALFA - AWAKE) .GT. 1.0E-5) LWAKE  = .FALSE.
        IF(ABS(ALFA - AVISC) .GT. 1.0E-5) LVCONV = .FALSE.
        IF(ABS(MINF - MVISC) .GT. 1.0E-5) LVCONV = .FALSE.
C
        itmaxs = ITMAX + 5
        IF(LVISC) THEN
          CALL VISCAL(itmaxs)
          conv_arr(i) = LVCONV
        ELSE
          conv_arr(i) = .TRUE.
          XOCTR(1) = 0.0
          XOCTR(2) = 0.0
        ENDIF
C
        a_arr(i)  = ALFA / DTOR
        cl_arr(i) = CL
        cd_arr(i) = CD
        cm_arr(i) = CM
        xtu_arr(i) = XOCTR(1)
        xtb_arr(i) = XOCTR(2)
C
      ENDDO
C
      END SUBROUTINE CSEQ

C=======================================================================
      SUBROUTINE GETNB(n_out)
C---------------------------------------------------
C     Returns the number of buffer airfoil points
C     for GETBUFFER.  Call this first to allocate
C     the buffers, then pass the result directly
C     to GETBUFFER.
C
C     Output: n_out  exact number of buffer points (= NB)
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INTEGER n_out
Cf2py intent(out) :: n_out
C
      n_out = NB
C
      END SUBROUTINE GETNB

C=======================================================================
      SUBROUTINE GETBUFFER(x_arr, y_arr, n_in)
C---------------------------------------------------
C     Returns current buffer airfoil coordinates.
C     Call GETBUFFER_N first to obtain n_in.
C
C     Input:  n_in   number of buffer points (= NB from GETBUFFER_N)
C     Output: x_arr  x coordinates (1:n_in)
C             y_arr  y coordinates (1:n_in)
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INTEGER n_in
      REAL*8 x_arr(n_in), y_arr(n_in)
Cf2py intent(in)  :: n_in
Cf2py intent(out) :: x_arr, y_arr
C
      INTEGER i
C
      DO i = 1, n_in
        x_arr(i) = XB(i)
        y_arr(i) = YB(i)
      ENDDO
C
      END SUBROUTINE GETBUFFER

C=======================================================================
      SUBROUTINE GETNBL(n_out)
C---------------------------------------------------
C     Returns the total number of boundary layer
C     data points for GET_BL.  Call this first to
C     allocate the buffers, then pass the result
C     directly to GET_BL.
C
C     The count includes airfoil panel nodes (N)
C     and, if a wake solution exists (LWAKE), the
C     wake nodes (NW), giving N + NW total.
C
C     Output: n_out  exact number of BL data points (= N + NW)
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INTEGER n_out
Cf2py intent(out) :: n_out
C
      IF(LWAKE) THEN
        n_out = N + NW
      ELSE
        n_out = N
      ENDIF
C
      END SUBROUTINE GETNBL

C=======================================================================
      SUBROUTINE GETN(n_out)
C---------------------------------------------------
C     Returns the number of panel nodes for GET_CP.
C     Call this first to allocate the buffer, then
C     pass the result directly to GET_CP.
C
C     Output: n_out  exact number of panel nodes (= N)
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INTEGER n_out
Cf2py intent(out) :: n_out
C
      n_out = N
C
      END SUBROUTINE GETN

C=======================================================================
      SUBROUTINE GETGEOM(x_arr, y_arr, n_in)
C---------------------------------------------------
C     Returns current airfoil panel coordinates.
C     Call GETN first to obtain n_in.
C
C     Input:  n_in   number of panel nodes (= N from GETN)
C     Output: x_arr  x coordinates (1:n_in)
C             y_arr  y coordinates (1:n_in)
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INTEGER n_in
      REAL*8 x_arr(n_in), y_arr(n_in)
Cf2py intent(in)  :: n_in
Cf2py intent(out) :: x_arr, y_arr
C
      INTEGER i
C
      DO i = 1, n_in
        x_arr(i) = X(i)
        y_arr(i) = Y(i)
      ENDDO
C
      END SUBROUTINE GETGEOM

C=======================================================================
      SUBROUTINE GETCP(cp_arr, n_in)
C---------------------------------------------------
C     Returns pressure coefficient distribution.
C     Returns CPV if viscous, CPI if inviscid.
C     Call GETN first to obtain n_in.
C
C     Input:  n_in   number of panel nodes (= N from GETN)
C     Output: cp_arr Cp values (1:n_in)
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INTEGER n_in
      REAL*8 cp_arr(n_in)
Cf2py intent(in)  :: n_in
Cf2py intent(out) :: cp_arr
C
      INTEGER i
C
      IF(LVISC) THEN
        DO i = 1, n_in
          cp_arr(i) = CPV(i)
        ENDDO
      ELSE
        DO i = 1, n_in
          cp_arr(i) = CPI(i)
        ENDDO
      ENDIF
C
      END SUBROUTINE GETCP

C=======================================================================
      SUBROUTINE GETCF(cf_arr, n_in)
C---------------------------------------------------
C     Returns pressure coefficient distribution.
C     Returns CPV if viscous, CPI if inviscid.
C     Call GETN first to obtain n_in.
C
C     Input:  n_in   number of panel nodes (= N from GETN)
C     Output: cp_arr Cp values (1:n_in)
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INTEGER n_in
      REAL*8 cf_arr(n_in)
Cf2py intent(in)  :: n_in
Cf2py intent(out) :: cf_arr
C
      INTEGER i
C
      DO 10 i = 1, N
        IS = 1
        IF(GAM(i) .LT. 0.0) IS = 2
C
        IF(LIPAN .AND. LVISC) THEN
          IF(IS .EQ. 1) THEN
            IBL = IBLTE(1) - i + 1
          ELSE
            IBL = IBLTE(2) + i - N
          ENDIF
          CF   = TAU(IBL,IS) / (0.5*QINF**2)
        ELSE
          CF   = 0.0
        ENDIF
C
        cf_arr(i)    = CF
  10  CONTINUE
C
      END SUBROUTINE GETCF

C=======================================================================
      SUBROUTINE GETBL(x_arr, y_arr, ue_arr, dstr_arr,
     &                  thet_arr, tstr_arr, hk_arr, 
     &                  n_in)
C---------------------------------------------------
C     Returns boundary layer data for airfoil panel
C     nodes and (if a wake solution exists) wake nodes.
C     Derived from BLDUMP without file I/O.
C     Call GET_BL_N first to obtain n_in.
C
C     Input:  n_in      number of BL data points (= N+NW from GET_BL_N)
C     Output: s_arr     arc length along surface
C             x_arr     x coordinates
C             y_arr     y coordinates
C             ue_arr    edge velocity Ue/Vinf (Karman-Tsien corrected)
C             dstr_arr  displacement thickness  (Dstar)
C             thet_arr  momentum thickness      (Theta)
C             tstr_arr  energy thickness        (Tstar)
C             hk_arr    kinematic shape factor Hk
C             cf_arr    skin friction coefficient
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
      INTEGER n_in
      REAL*8 x_arr(n_in), y_arr(n_in)
      REAL*8 ue_arr(n_in), dstr_arr(n_in), thet_arr(n_in)
      REAL*8 tstr_arr(n_in),  hk_arr(n_in)    
Cf2py intent(in)  :: n_in
Cf2py intent(out) :: x_arr, y_arr, ue_arr, dstr_arr, thet_arr
Cf2py intent(out) :: tstr_arr, hk_arr
C
      INTEGER i, j, k, IS, IBL
      REAL*8 DS, TH, TS, CF, H, HS, UE, UI, AMSQ, HK, DUMMY
C
C---- update compressibility parameters (HSTINV, TKLAM, etc.)
      CALL COMSET
C
C---- airfoil nodes (indices 1:N)
      DO 10 i = 1, N
        IS = 1
        IF(GAM(i) .LT. 0.0) IS = 2
C
        IF(LIPAN .AND. LVISC) THEN
          IF(IS .EQ. 1) THEN
            IBL = IBLTE(1) - i + 1
          ELSE
            IBL = IBLTE(2) + i - N
          ENDIF
          DS   = DSTR(IBL,IS)
          TH   = THET(IBL,IS)
          TS   = TSTR(IBL,IS)
          IF(TH .EQ. 0.0) THEN
            H  = 1.0
            HS = 1.0
          ELSE
            H  = DS/TH
          ENDIF
        ELSE
          DS   = 0.0
          TH   = 0.0
          TS   = 0.0
          H    = 1.0
        ENDIF
C
        UE   = (GAM(i)/QINF) * (1.0-TKLAM) /
     &         (1.0 - TKLAM*(GAM(i)/QINF)**2)
        AMSQ = UE*UE*HSTINV / (GAMM1*(1.0 - 0.5*UE*UE*HSTINV))
        CALL HKIN(H, AMSQ, HK, DUMMY, DUMMY)
C
        x_arr(i)     = X(i)
        y_arr(i)     = Y(i)
        ue_arr(i)    = UE
        dstr_arr(i)  = DS
        thet_arr(i)  = TH
        tstr_arr(i)  = TS
        hk_arr(i)    = HK
  10  CONTINUE
C
C---- wake nodes (indices N+1:N+NW) if wake solution exists
C     H*, P, m, K are not defined in the wake; set to zero
      IF(LWAKE) THEN
        IS = 2
        DO 20 j = 1, NW
          i   = N + j
          k   = N + j
          IBL = IBLTE(2) + j
          DS  = DSTR(IBL,IS)
          TH  = THET(IBL,IS)
          IF(TH .EQ. 0.0) THEN
            H = 1.0
          ELSE
            H = DS/TH
          ENDIF
          UI   = UEDG(IBL,IS)
          UE   = (UI/QINF) * (1.0-TKLAM) /
     &           (1.0 - TKLAM*(UI/QINF)**2)
          AMSQ = UE*UE*HSTINV / (GAMM1*(1.0 - 0.5*UE*UE*HSTINV))
          CALL HKIN(H, AMSQ, HK, DUMMY, DUMMY)
C
          x_arr(k)     = X(i)
          y_arr(k)     = Y(i)
          ue_arr(k)    = UE
          dstr_arr(k)  = DS
          thet_arr(k)  = TH
          tstr_arr(k)  = TS
          hk_arr(k)    = HK
  20    CONTINUE
      ENDIF
C
      END SUBROUTINE GETBL

C=======================================================================
      SUBROUTINE GETMOREBL(s_arr, x_arr, y_arr, ue_arr, dstr_arr,
     &                  thet_arr, tstr_arr, hk_arr, hstr_arr, 
     &                  n_in)
C---------------------------------------------------
C     Returns boundary layer data for airfoil panel
C     nodes and (if a wake solution exists) wake nodes.
C     Derived from BLDUMP without file I/O.
C     Call GET_BL_N first to obtain n_in.
C
C     Input:  n_in      number of BL data points (= N+NW from GET_BL_N)
C     Output: s_arr     arc length along surface
C             x_arr     x coordinates
C             y_arr     y coordinates
C             ue_arr    edge velocity Ue/Vinf (Karman-Tsien corrected)
C             dstr_arr  displacement thickness  (Dstar)
C             thet_arr  momentum thickness      (Theta)
C             tstr_arr  energy thickness        (Tstar)
C             hk_arr    kinematic shape factor Hk
C             cf_arr    skin friction coefficient
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
      INTEGER n_in
      REAL*8 s_arr(n_in), x_arr(n_in), y_arr(n_in)
      REAL*8 ue_arr(n_in), dstr_arr(n_in), thet_arr(n_in)
      REAL*8 tstr_arr(n_in),  hk_arr(n_in), hstr_arr(n_in)    
Cf2py intent(in)  :: n_in
Cf2py intent(out) :: s_arr, x_arr, y_arr, ue_arr, dstr_arr, thet_arr
Cf2py intent(out) :: tstr_arr, hk_arr, hstr_arr
C
      INTEGER i, j, k, IS, IBL
      REAL*8 DS, TH, TS, CF, H, HS, UE, UI, AMSQ, HK, DUMMY
C
C---- update compressibility parameters (HSTINV, TKLAM, etc.)
      CALL COMSET
C
C---- airfoil nodes (indices 1:N)
      DO 10 i = 1, N
        IS = 1
        IF(GAM(i) .LT. 0.0) IS = 2
C
        IF(LIPAN .AND. LVISC) THEN
          IF(IS .EQ. 1) THEN
            IBL = IBLTE(1) - i + 1
          ELSE
            IBL = IBLTE(2) + i - N
          ENDIF
          DS   = DSTR(IBL,IS)
          TH   = THET(IBL,IS)
          TS   = TSTR(IBL,IS)
          IF(TH .EQ. 0.0) THEN
            H  = 1.0
            HS = 1.0
          ELSE
            H  = DS/TH
            HS = TS/TH
          ENDIF
        ELSE
          DS   = 0.0
          TH   = 0.0
          TS   = 0.0
          H    = 1.0
          HS   = 2.0
        ENDIF
C
        UE   = (GAM(i)/QINF) * (1.0-TKLAM) /
     &         (1.0 - TKLAM*(GAM(i)/QINF)**2)
        AMSQ = UE*UE*HSTINV / (GAMM1*(1.0 - 0.5*UE*UE*HSTINV))
        CALL HKIN(H, AMSQ, HK, DUMMY, DUMMY)
C
        s_arr(i)     = S(i)
        x_arr(i)     = X(i)
        y_arr(i)     = Y(i)
        ue_arr(i)    = UE
        dstr_arr(i)  = DS
        thet_arr(i)  = TH
        tstr_arr(i)  = TS
        hk_arr(i)    = HK
        hstr_arr(i)  = HS
  10  CONTINUE
C
C---- wake nodes (indices N+1:N+NW) if wake solution exists
C     H*, P, m, K are not defined in the wake; set to zero
      IF(LWAKE) THEN
        IS = 2
        DO 20 j = 1, NW
          i   = N + j
          k   = N + j
          IBL = IBLTE(2) + j
          DS  = DSTR(IBL,IS)
          TH  = THET(IBL,IS)
          TS   = TSTR(IBL,IS)
          IF(TH .EQ. 0.0) THEN
            H = 1.0
            HS = 1.0
          ELSE
            H = DS/TH
            HS = TS/TH
          ENDIF
          UI   = UEDG(IBL,IS)
          UE   = (UI/QINF) * (1.0-TKLAM) /
     &           (1.0 - TKLAM*(UI/QINF)**2)
          AMSQ = UE*UE*HSTINV / (GAMM1*(1.0 - 0.5*UE*UE*HSTINV))
          CALL HKIN(H, AMSQ, HK, DUMMY, DUMMY)
C
          s_arr(k)     = S(i)
          x_arr(k)     = X(i)
          y_arr(k)     = Y(i)
          ue_arr(k)    = UE
          dstr_arr(k)  = DS
          thet_arr(k)  = TH
          tstr_arr(k)  = TS
          hk_arr(k)    = HK
          hstr_arr(k)  = HS
  20    CONTINUE
      ENDIF
C
      END SUBROUTINE GETMOREBL

C=======================================================================
      SUBROUTINE GETBLALL(s_arr, x_arr, y_arr, ue_arr, dstr_arr,
     &                  thet_arr, hk_arr, cdis_arr, ct_arr,
     &                   P_arr, m_arr, K_arr,
     &                  n_in)
C---------------------------------------------------
C     Returns boundary layer data for airfoil panel
C     nodes and (if a wake solution exists) wake nodes.
C     Derived from BLDUMP without file I/O.
C     Call GET_BL_N first to obtain n_in.
C
C     Input:  n_in      number of BL data points (= N+NW from GET_BL_N)
C     Output: s_arr     arc length along surface
C             x_arr     x coordinates
C             y_arr     y coordinates
C             ue_arr    edge velocity Ue/Vinf (Karman-Tsien corrected)
C             dstr_arr  displacement thickness  (Dstar)
C             thet_arr  momentum thickness      (Theta)
C              tstr_arr  energy thickness        (Tstar)
C             hk_arr    kinematic shape factor Hk
C             cf_arr    skin friction coefficient
C             cdis_arr  dissipation coefficient
C             ct_arr    sqrt(max shear stress coefficient)
C             P_arr     kinematic momentum flux (= Theta * Ue^2)
C                       zero for wake nodes
C             m_arr     kinematic mass flux      (= Dstar * Ue)
C                       zero for wake nodes
C             K_arr     kinematic energy flux    (= Tstar * Ue^3)
C                       zero for wake nodes
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
      INTEGER n_in
      REAL*8 s_arr(n_in),     x_arr(n_in),    y_arr(n_in)
      REAL*8 ue_arr(n_in),    dstr_arr(n_in), thet_arr(n_in)
      REAL*8 tstr_arr(n_in), hk_arr(n_in)
      REAL*8 cdis_arr(n_in),  ct_arr(n_in)
      REAL*8 P_arr(n_in), m_arr(n_in), K_arr(n_in)
Cf2py intent(in)  :: n_in
Cf2py intent(out) :: s_arr, x_arr, y_arr, ue_arr, dstr_arr, thet_arr
Cf2py intent(out) :: tstr_arr, hk_arr, cdis_arr, ct_arr
Cf2py intent(out) :: P_arr, m_arr, K_arr
C
      INTEGER i, j, k, IS, IBL
      REAL*8 DS, TH, TS, CF, H, HS, UE, UI, AMSQ, HK, CDIS, CT, DUMMY
C
C---- update compressibility parameters (HSTINV, TKLAM, etc.)
      CALL COMSET
C
C---- airfoil nodes (indices 1:N)
      DO 10 i = 1, N
        IS = 1
        IF(GAM(i) .LT. 0.0) IS = 2
C
        IF(LIPAN .AND. LVISC) THEN
          IF(IS .EQ. 1) THEN
            IBL = IBLTE(1) - i + 1
          ELSE
            IBL = IBLTE(2) + i - N
          ENDIF
          DS   = DSTR(IBL,IS)
          TH   = THET(IBL,IS)
          TS   = TSTR(IBL,IS)
          IF(TH .EQ. 0.0) THEN
            H  = 1.0
          ELSE
            H  = DS/TH
          ENDIF
          CDIS = DIS(IBL,IS)  / QINF**3
          CT   = CTAU(IBL,IS)
        ELSE
          DS   = 0.0
          TH   = 0.0
          TS   = 0.0
          H    = 1.0
          CDIS = 0.0
          CT   = 0.0
        ENDIF
C
        UE   = (GAM(i)/QINF) * (1.0-TKLAM) /
     &         (1.0 - TKLAM*(GAM(i)/QINF)**2)
        AMSQ = UE*UE*HSTINV / (GAMM1*(1.0 - 0.5*UE*UE*HSTINV))
        CALL HKIN(H, AMSQ, HK, DUMMY, DUMMY)
C
        s_arr(i)     = S(i)
        x_arr(i)     = X(i)
        y_arr(i)     = Y(i)
        ue_arr(i)    = UE
        dstr_arr(i)  = DS
        thet_arr(i)  = TH
        tstr_arr(i)  = TS
        hk_arr(i)    = HK
        cdis_arr(i)  = CDIS
        ct_arr(i)    = CT
        P_arr(i)     = TH * UE**2
        m_arr(i)     = DS * UE
        K_arr(i)     = TS * UE**3
  10  CONTINUE
C
C---- wake nodes (indices N+1:N+NW) if wake solution exists
C     H*, P, m, K are not defined in the wake; set to zero
      IF(LWAKE) THEN
        IS = 2
        DO 20 j = 1, NW
          i   = N + j
          k   = N + j
          IBL = IBLTE(2) + j
          DS  = DSTR(IBL,IS)
          TH  = THET(IBL,IS)
          TS  = TSTR(IBL,IS)
          IF(TH .EQ. 0.0) THEN
            H = 1.0
          ELSE
            H = DS/TH
          ENDIF
          UI   = UEDG(IBL,IS)
          UE   = (UI/QINF) * (1.0-TKLAM) /
     &           (1.0 - TKLAM*(UI/QINF)**2)
          AMSQ = UE*UE*HSTINV / (GAMM1*(1.0 - 0.5*UE*UE*HSTINV))
          CALL HKIN(H, AMSQ, HK, DUMMY, DUMMY)
          CDIS = DIS(IBL,IS)  / QINF**3
          CT   = CTAU(IBL,IS)
C
          s_arr(k)     = S(i)
          x_arr(k)     = X(i)
          y_arr(k)     = Y(i)
          ue_arr(k)    = UE
          dstr_arr(k)  = DS
          thet_arr(k)  = TH
          tstr_arr(k)  = TS
          hk_arr(k)    = HK
          cdis_arr(k)  = CDIS
          ct_arr(k)    = CT
          P_arr(k)     = 0.0
          m_arr(k)     = 0.0
          K_arr(k)     = 0.0
  20    CONTINUE
      ENDIF
C
      END SUBROUTINE GETBLALL

C=======================================================================
C---- GDES menu wrappers
C=======================================================================

C=======================================================================
      SUBROUTINE SETFLAP(Xcloc, Ycloc, DegDef)
C---------------------------------------------------
C     Adds a flap to the current airfoil geometry.
C     Xcloc, Ycloc: hinge point location (normalised to chord)
C     DegDef: flap deflection angle in degrees (positive deflects trailing edge down)
C---------------------------------------------------
      REAL*8 Xcloc, Ycloc, DegDef
Cf2py intent(in) :: Xcloc, Ycloc, DegDef

      CALL FLAP(Xcloc, Ycloc, DegDef)
      END SUBROUTINE SETFLAP
      
C=======================================================================
      SUBROUTINE GETGEOINFO(chord_out, area_out, radle_out,
     &                      angbte_out, thick_out, cambr_out)
C---------------------------------------------------
C     Returns the full set of geometric parameters
C     for the buffer airfoil.
C     Synchronises buffer and current airfoil via
C     ABCOPY before querying.
C
C     Output: chord_out   chord length
C             area_out    enclosed area
C             radle_out   leading-edge radius
C             angbte_out  trailing-edge included angle (rad)
C             thick_out   maximum thickness
C             cambr_out   maximum camber
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 chord_out, area_out, radle_out
      REAL*8 angbte_out, thick_out, cambr_out
Cf2py intent(out) :: chord_out, area_out, radle_out
Cf2py intent(out) :: angbte_out, thick_out, cambr_out
C
C---- synchronise current airfoil from buffer
      CALL ABCOPY(.FALSE.)
C
C---- recompute geometry parameters on the buffer airfoil
      CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
     &            SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
     &            EI11BA,EI22BA,APX1BA,APX2BA,
     &            EI11BT,EI22BT,APX1BT,APX2BT,
     &            THICKB,CAMBRB )
C
      chord_out  = CHORDB
      area_out   = AREAB
      radle_out  = RADBLE
      angbte_out = ANGBTE
      thick_out  = THICKB
      cambr_out  = CAMBRB
C
      END SUBROUTINE GETGEOINFO

C=======================================================================
      SUBROUTINE GETPANANGLE(imax_out, amax_out)
C---------------------------------------------------
C     Returns the index and magnitude of the maximum
C     panel corner angle across the buffer airfoil.
C     Terminal output is suppressed (IPRINT=0).
C
C     Output: imax_out  index of max corner angle point
C             amax_out  max corner angle (degrees)
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INTEGER imax_out
      REAL*8 amax_out
Cf2py intent(out) :: imax_out, amax_out
C
      CALL CANG(XB,YB,NB, 0, imax_out, amax_out)
C
      END SUBROUTINE GETPANANGLE

C=======================================================================
      SUBROUTINE GETSTRUCT(area_out, slen_out,
     &                     xc_out, yc_out,
     &                     aixx_out, aixxt_out,
     &                     aiyy_out, aiyyt_out,
     &                     aj_out, ajt_out)
C---------------------------------------------------
C     Returns structural cross-section integrals
C     for the buffer airfoil (GDES BEND command).
C     Uses power exponent PEX = 16.0 internally.
C
C     Output: area_out   enclosed area
C             slen_out   perimeter length
C             xc_out     centroid x coordinate
C             yc_out     centroid y coordinate
C             aixx_out   solid Ixx (second moment about x)
C             aixxt_out  skin-weighted Ixx
C             aiyy_out   solid Iyy (second moment about y)
C             aiyyt_out  skin-weighted Iyy
C             aj_out     solid torsional constant
C             ajt_out    skin torsional constant
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 area_out, slen_out, xc_out, yc_out
      REAL*8 aixx_out, aixxt_out, aiyy_out, aiyyt_out
      REAL*8 aj_out, ajt_out
Cf2py intent(out) :: area_out, slen_out, xc_out, yc_out
Cf2py intent(out) :: aixx_out, aixxt_out, aiyy_out, aiyyt_out
Cf2py intent(out) :: aj_out, ajt_out
C
      REAL*8 PEX
      REAL*8 XMIN, XMAX, XEXINT
      REAL*8 YMIN, YMAX, YEXINT
C
      PEX = 16.0
C
      CALL IJSECT(NB, XB, YB, PEX,
     &  area_out, slen_out,
     &  xc_out, XMIN, XMAX, XEXINT,
     &  yc_out, YMIN, YMAX, YEXINT,
     &  aixx_out, aixxt_out,
     &  aiyy_out, aiyyt_out,
     &  aj_out, ajt_out)
C
      END SUBROUTINE GETSTRUCT

C=======================================================================
      SUBROUTINE ANGREFINE(atol_in, x1_in, x2_in)
C---------------------------------------------------
C     Inserts additional buffer airfoil points where
C     the panel corner angle exceeds atol_in (degrees),
C     within the chordwise x-range [x1_in, x2_in].
C     The buffer is re-splined and geometry parameters
C     are recomputed internally.
C
C     Input: atol_in   angle tolerance (degrees)
C            x1_in     start of x-range for refinement
C            x2_in     end of x-range for refinement
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 atol_in, x1_in, x2_in
Cf2py intent(in) :: atol_in, x1_in, x2_in
C
      INTEGER NNEW
      REAL*8 XNEW(2*IBX), YNEW(2*IBX)
C
C---- refine buffer airfoil where angles exceed threshold
      CALL AREFINE(XB,YB,SB,XBP,YBP,NB, atol_in,
     &             2*IBX, NNEW, XNEW, YNEW, x1_in, x2_in)
C
      IF(NNEW .GT. 0) THEN
C------ copy refined coordinates back into buffer arrays
        NB = NNEW
        DO I = 1, NB
          XB(I) = XNEW(I)
          YB(I) = YNEW(I)
        ENDDO
C
C------ re-spline
        CALL SCALC(XB,YB,SB,NB)
        CALL SEGSPL(XB,XBP,SB,NB)
        CALL SEGSPL(YB,YBP,SB,NB)
C
C------ recompute geometry parameters
        CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
     &              SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
     &              EI11BA,EI22BA,APX1BA,APX2BA,
     &              EI11BT,EI22BT,APX1BT,APX2BT,
     &              THICKB,CAMBRB )
      ENDIF
C
      LGSAME = .FALSE.
C
      END SUBROUTINE ANGREFINE

C=======================================================================
      SUBROUTINE INTERP(X2_in, Y2_in, N2_in, frac_in)
C---------------------------------------------------
C     Produces an intermediate buffer airfoil by
C     blending the current airfoil (source 0) with a
C     second user-supplied airfoil (source 1) at a
C     given fraction.
C       frac_in = 0.0 -> current airfoil
C       frac_in = 1.0 -> supplied airfoil
C
C     The second airfoil is splined internally.
C     The result replaces the buffer airfoil.
C
C     Input: X2_in(N2_in)  x-coords of second airfoil
C            Y2_in(N2_in)  y-coords of second airfoil
C            N2_in         number of points
C            frac_in       blending fraction (0..1)
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INTEGER N2_in
      REAL*8 X2_in(N2_in), Y2_in(N2_in)
      REAL*8 frac_in
Cf2py intent(in) :: X2_in, Y2_in, N2_in, frac_in
C
      REAL*8 X2(IBX), Y2(IBX), XP2(IBX), YP2(IBX), S2(IBX)
      REAL*8 SLE2
      REAL*8 XOUT(IBX), YOUT(IBX)
      INTEGER N2, NOUT, I
C
C---- copy second airfoil into local arrays
      N2 = MIN(N2_in, IBX)
      DO I = 1, N2
        X2(I) = X2_in(I)
        Y2(I) = Y2_in(I)
      ENDDO
C
C---- spline the second airfoil
      CALL SCALC(X2,Y2,S2,N2)
      CALL SEGSPL(X2,XP2,S2,N2)
      CALL SEGSPL(Y2,YP2,S2,N2)
      CALL LEFIND(SLE2,X2,XP2,Y2,YP2,S2,N2)
C
C---- interpolate using the buffer airfoil as source 0
      CALL INTER(XB,XBP,YB,YBP,SB,NB,SBLE,
     &           X2,XP2,Y2,YP2,S2,N2,SLE2,
     &           XOUT,YOUT,NOUT, frac_in)
C
C---- copy result into buffer arrays
      NB = NOUT
      DO I = 1, NB
        XB(I) = XOUT(I)
        YB(I) = YOUT(I)
      ENDDO
C
C---- re-spline buffer
      CALL SCALC(XB,YB,SB,NB)
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
C
C---- recompute geometry parameters
      CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
     &            SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
     &            EI11BA,EI22BA,APX1BA,APX2BA,
     &            EI11BT,EI22BT,APX1BT,APX2BT,
     &            THICKB,CAMBRB )
C
      LGSAME = .FALSE.
C
      END SUBROUTINE INTERP

C=======================================================================
C---- GDES Phase 2 wrappers
C=======================================================================

C=======================================================================
      SUBROUTINE ROTATEGEOM(angle_in, in_degrees)
C---------------------------------------------------
C     Rotates all buffer airfoil points about the
C     origin by a specified angle.
C     Corresponds to GDES commands ADEG / ARAD.
C
C     Input: angle_in    rotation angle
C            in_degrees  .TRUE.  -> angle is in degrees
C                        .FALSE. -> angle is in radians
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 angle_in
      LOGICAL in_degrees
Cf2py intent(in) :: angle_in, in_degrees
C
      REAL*8 ALFA_RAD
C
      IF(in_degrees) THEN
        ALFA_RAD = angle_in * DTOR
      ELSE
        ALFA_RAD = angle_in
      ENDIF
C
      CALL ROTATE(XB, YB, NB, ALFA_RAD)
C
C---- re-spline buffer
      CALL SCALC(XB,YB,SB,NB)
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
C
C---- recompute geometry parameters
      CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
     &            SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
     &            EI11BA,EI22BA,APX1BA,APX2BA,
     &            EI11BT,EI22BT,APX1BT,APX2BT,
     &            THICKB,CAMBRB )
C
      LGSAME = .FALSE.
C
      END SUBROUTINE ROTATEGEOM

C=======================================================================
      SUBROUTINE SETTGAP(gapnew_in, doc_in)
C---------------------------------------------------
C     Modifies the trailing-edge gap of the buffer
C     airfoil to a specified value while blending
C     smoothly into the rest of the geometry.
C     Corresponds to GDES command TGAP.
C
C     Input: gapnew_in  desired new TE gap
C            doc_in     blending distance/c (0..1)
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 gapnew_in, doc_in
Cf2py intent(in) :: gapnew_in, doc_in
C
      CALL TGAP(gapnew_in, doc_in)
C
      END SUBROUTINE SETTGAP

C=======================================================================
      SUBROUTINE SETLERAD(rfac_in, doc_in)
C---------------------------------------------------
C     Scales the leading-edge radius of the buffer
C     airfoil by a target ratio while preserving the
C     overall shape away from the leading edge.
C     Corresponds to GDES command LERA.
C
C     Input: rfac_in  approx. new/old LE radius
C                     scaling ratio
C            doc_in   blending distance/c from LE
C                     (clamped to >= 0.001)
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 rfac_in, doc_in
Cf2py intent(in) :: rfac_in, doc_in
C
      CALL LERAD(rfac_in, doc_in)
C
      END SUBROUTINE SETLERAD

C=======================================================================
      SUBROUTINE SCALETHCAM(tfac_in, cfac_in)
C---------------------------------------------------
C     Scales the existing thickness and camber
C     distributions of the buffer airfoil by
C     independent factors.
C     Corresponds to GDES command TFAC.
C
C     Input: tfac_in  thickness scale factor
C            cfac_in  camber    scale factor
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 tfac_in, cfac_in
Cf2py intent(in) :: tfac_in, cfac_in
C
      CALL TCSCAL(tfac_in, cfac_in)
C
      END SUBROUTINE SCALETHCAM

C=======================================================================
      SUBROUTINE SETTHCAM(tnew_in, cnew_in)
C---------------------------------------------------
C     Sets the buffer airfoil to a specified maximum
C     thickness and/or camber value, adjusting the
C     distributions while preserving their shapes.
C     Corresponds to GDES command TSET.
C
C     Input: tnew_in  new max thickness
C                     (use 999.0 to skip / leave unchanged)
C            cnew_in  new max camber
C                     (use 999.0 to skip / leave unchanged)
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 tnew_in, cnew_in
Cf2py intent(in) :: tnew_in, cnew_in
C
      CALL TCSET(tnew_in, cnew_in)
C
      END SUBROUTINE SETTHCAM

C=======================================================================
      SUBROUTINE MOVHIPNT(thpnt_in, chpnt_in)
C---------------------------------------------------
C     Shifts the chordwise location of the thickness
C     and/or camber highpoints of the buffer airfoil
C     while preserving the peak values.
C     Corresponds to GDES command HIGH.
C
C     Input: thpnt_in  new thickness highpoint x/c
C                      (use 0.0 to keep current)
C            chpnt_in  new camber highpoint x/c
C                      (use 0.0 to keep current)
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL*8 thpnt_in, chpnt_in
Cf2py intent(in) :: thpnt_in, chpnt_in
C
      CALL HIPNT(thpnt_in, chpnt_in)
C
      END SUBROUTINE MOVHIPNT

C=======================================================================
      SUBROUTINE SETAIRFOIL()
C---------------------------------------------------
C     Copies the buffer airfoil to the current
C     airfoil, re-splines, and resets all dependent
C     flags. Call after any GDES buffer modifications
C     to make them active for analysis.
C---------------------------------------------------
      CALL ABCOPY(.FALSE.)
C
      END SUBROUTINE SETAIRFOIL