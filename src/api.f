C***********************************************************************
C    Module:  api.f
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
C---- API subroutines for f2py Python interface
C
C=======================================================================
      SUBROUTINE Initialise()
C---------------------------------------
C     Initialises all XFOIL variables.
C     Must be called before any other API
C     subroutine.
C---------------------------------------
      CALL INIT
      END SUBROUTINE Initialise

C=======================================================================
      SUBROUTINE AIRFOIL(X_in, Y_in, N_in)
C-----------------------------------------------
C     Loads airfoil coordinates, generates panels.
C-----------------------------------------------
      INTEGER N_in
      REAL X_in(N_in), Y_in(N_in)
Cf2py intent(in) :: X_in, Y_in, N_in
C
      CALL ADDGEOM(X_in, Y_in, N_in)
C
      END SUBROUTINE AIRFOIL

C=======================================================================
      SUBROUTINE NACA_LOAD(desig)
C---------------------------------------------------
C     Generates and panels a NACA 4 or 5 digit
C     airfoil from a character designation string.
C     e.g. desig = "0012" -> NACA 0012
C          desig = "23015" -> NACA 23015
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      CHARACTER*10 desig
Cf2py intent(in) :: desig
C
      INTEGER IDES
C
C---- read integer designation from string
C     "0012" -> 12, "2412" -> 2412, "23015" -> 23015
      READ(desig, *) IDES
C
C---- call existing NACA routine which handles 4 and 5 digit cases,
C     generates buffer airfoil coords and calls PANGEN
      CALL NACA(IDES)
C
      END SUBROUTINE NACA_LOAD
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
      REAL    Au(NAu), Al(NAl)
      REAL    LEM, TE, N1, N2
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
      REAL rle_up, Xup, Zup, Zxxup
      REAL rle_lo, Xlo, Zlo, Zxxlo
      REAL Zte, dZte, alpha_te_deg, beta_te_deg
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
      SUBROUTINE SETCON(Re, Ma)
C---------------------------------------------------
C     Sets Reynolds and Mach numbers.
C     Setting Re > 0 automatically enables viscous
C     mode (LVISC = .TRUE.).
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL Re, Ma
Cf2py intent(in) :: Re, Ma
C
      REAL MINF_CLM, REINF_CLM
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
      END SUBROUTINE SETCON
C
C=======================================================================
      SUBROUTINE SETCCOR(icor)
C---------------------------------------------------
C     Selects compressibility correction method.
C     icor = 0  Prandtl-Glauert
C     icor = 1  Karman-Tsien  (default)
C     icor = 2  Laitone
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INTEGER icor
Cf2py intent(in) :: icor
      ICCOR = icor
      LVCONV = .FALSE.
      CALL COMSET
      END SUBROUTINE SETCCOR
C
C=======================================================================
      SUBROUTINE SETITER(Iter)
C---------------------------------------------------
C     Sets the maximum number of viscous iterations.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INTEGER Iter
Cf2py intent(in) :: Iter
C
      ITMAX = Iter
C
      END SUBROUTINE SETITER
C
C=======================================================================
      SUBROUTINE TRPARS(Xtrip1, Xtrip2, Ncrit)
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
      REAL Xtrip1, Xtrip2, Ncrit
Cf2py intent(in) :: Xtrip1, Xtrip2, Ncrit
C
      XSTRIP(1) = Xtrip1
      XSTRIP(2) = Xtrip2
      ACRIT     = Ncrit
C
      END SUBROUTINE TRPARS

C=======================================================================
      SUBROUTINE REPANEL(NPan_in, CVPar_in, CTERat_in, CTRRat_in,
     &                   XSREF1_in, XSREF2_in, XPREF1_in, XPREF2_in)
C---------------------------------------------------
C     Sets paneling parameters and re-panels the
C     airfoil if geometry has been loaded.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INTEGER NPan_in
      REAL CVPar_in, CTERat_in, CTRRat_in
      REAL XSREF1_in, XSREF2_in, XPREF1_in, XPREF2_in
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
C
C---- re-panel only if buffer airfoil exists
      IF(NB .GT. 0) THEN
        CALL PANGEN(.FALSE.)
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
      SUBROUTINE Quiet(setting)
C---------------------------------------------------
C     Enables/disables solver iteration output.
C     setting = .TRUE.  suppresses WRITE output
C     setting = .FALSE. allows WRITE output
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      LOGICAL setting
Cf2py intent(in) :: setting
C
      LQUIET = setting
C
      END SUBROUTINE Quiet

C=======================================================================
      SUBROUTINE INIT_BL()
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
      END SUBROUTINE INIT_BL

C=======================================================================
      SUBROUTINE SET_VISCOUS(flag)
C---------------------------------------------------
C     Toggles viscous/inviscid mode without
C     changing Re or Ma.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      LOGICAL flag
Cf2py intent(in) :: flag
C
      LVISC  = flag
      LVCONV = .FALSE.
C
      END SUBROUTINE SET_VISCOUS

C=======================================================================
      SUBROUTINE alpha(alpha_deg,
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
      REAL    alpha_deg
      REAL    cl_out, cd_out, cm_out, xtu_out, xtb_out
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
      ENDIF
C
      cl_out  = CL
      cd_out  = CD
      cm_out  = CM
      xtu_out = XOCTR(1)
      xtb_out = XOCTR(2)
C
      END SUBROUTINE alpha

C=======================================================================
      SUBROUTINE ASeq(alpha_i, alpha_f, n_alpha,
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
      REAL    alpha_i, alpha_f
      INTEGER n_alpha
      REAL    a_arr(n_alpha), cl_arr(n_alpha), cd_arr(n_alpha)
      REAL    cm_arr(n_alpha), xtu_arr(n_alpha), xtb_arr(n_alpha)
      LOGICAL conv_arr(n_alpha)
Cf2py intent(in)  :: alpha_i, alpha_f, n_alpha
Cf2py intent(out) :: a_arr, cl_arr, cd_arr, cm_arr
Cf2py intent(out) :: xtu_arr, xtb_arr, conv_arr
C
      INTEGER i, itmaxs
      REAL a0, da
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
      END SUBROUTINE ASeq

C=======================================================================
      SUBROUTINE CL_spec(cl_val,
     &                   a_out, cl_out, cd_out, cm_out, cp_out,
     &                   conv_out)
C---------------------------------------------------
C     Analyses airfoil at specified lift coefficient.
C     Input:  cl_val    target CL
C     Output: a_out     solved angle of attack (degrees)
C             cl_out    actual CL achieved
C             cd_out    drag coefficient
C             cm_out    moment coefficient
C             cp_out    minimum surface Cp
C             conv_out  .TRUE. if converged
C
C     Note: AWAKE/AVISC/MVISC checks are done AFTER SPECCL
C     because SPECCL solves for ALFA -- it is unknown beforehand.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL    cl_val
      REAL    a_out, cl_out, cd_out, cm_out, cp_out
      LOGICAL conv_out
Cf2py intent(in)  :: cl_val
Cf2py intent(out) :: a_out, cl_out, cd_out, cm_out, cp_out
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
      ENDIF
C
      a_out  = ALFA / DTOR
      cl_out = CL
      cd_out = CD
      cm_out = CM
C
      CALL FCPMIN
      cp_out = CPMN
C
      END SUBROUTINE CL_spec

C=======================================================================
      SUBROUTINE CSeq(cl_start, cl_end, n_step,
     &                a_arr, cl_arr, cd_arr, cm_arr, cp_arr,
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
C             cp_arr     minimum surface Cp array
C             conv_arr   convergence flag array
C
C     Note: CL step is (cl_end-cl_start)/n_step so the endpoint
C     cl_end is not evaluated (analogous to Python range() semantics).
C     AWAKE/AVISC/MVISC checks are applied after each SPECCL call
C     because SPECCL solves for ALFA.
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL    cl_start, cl_end
      INTEGER n_step
      REAL    a_arr(n_step), cl_arr(n_step), cd_arr(n_step)
      REAL    cm_arr(n_step), cp_arr(n_step)
      LOGICAL conv_arr(n_step)
Cf2py intent(in)  :: cl_start, cl_end, n_step
Cf2py intent(out) :: a_arr, cl_arr, cd_arr, cm_arr, cp_arr, conv_arr
C
      INTEGER i, itmaxs
      REAL    cl0, dcl
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
        ENDIF
C
        a_arr(i)  = ALFA / DTOR
        cl_arr(i) = CL
        cd_arr(i) = CD
        cm_arr(i) = CM
C
        CALL FCPMIN
        cp_arr(i) = CPMN
      ENDDO
C
      END SUBROUTINE CSeq

C=======================================================================
      SUBROUTINE GET_CP(cp_arr, nbuf, n_out)
C---------------------------------------------------
C     Returns pressure coefficient distribution.
C     Returns CPV if viscous, CPI if inviscid.
C     Input:  nbuf   size of caller-provided buffer
C     Output: cp_arr Cp values (1:n_out)
C             n_out  number of valid entries (= N)
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INTEGER nbuf, n_out
      REAL    cp_arr(nbuf)
Cf2py intent(in)  :: nbuf
Cf2py intent(out) :: cp_arr, n_out
C
      INTEGER i
C
      n_out = MIN(N, nbuf)
C
      IF(LVISC) THEN
        DO i = 1, n_out
          cp_arr(i) = CPV(i)
        ENDDO
      ELSE
        DO i = 1, n_out
          cp_arr(i) = CPI(i)
        ENDDO
      ENDIF
C
      END SUBROUTINE GET_CP

C=======================================================================
      SUBROUTINE GET_XY(x_arr, y_arr, nbuf, n_out)
C---------------------------------------------------
C     Returns current airfoil panel coordinates.
C     Input:  nbuf   size of caller-provided buffer
C     Output: x_arr  x coordinates (1:n_out)
C             y_arr  y coordinates (1:n_out)
C             n_out  number of valid entries (= N)
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INTEGER nbuf, n_out
      REAL    x_arr(nbuf), y_arr(nbuf)
Cf2py intent(in)  :: nbuf
Cf2py intent(out) :: x_arr, y_arr, n_out
C
      INTEGER i
C
      n_out = MIN(N, nbuf)
C
      DO i = 1, n_out
        x_arr(i) = X(i)
        y_arr(i) = Y(i)
      ENDDO
C
      END SUBROUTINE GET_XY

C=======================================================================
      SUBROUTINE GET_BL(s_arr, x_arr, y_arr, ue_arr, dstr_arr,
     &                  thet_arr, cf_arr, hk_arr, cdis_arr, ct_arr,
     &                  nbuf, n_out)
C---------------------------------------------------
C     Returns boundary layer data for airfoil panel
C     nodes and (if a wake solution exists) wake nodes.
C     Derived from BLDUMP without file I/O.
C
C     Input:  nbuf      size of caller-provided buffers
C     Output: s_arr     arc length along surface
C             x_arr     x coordinates
C             y_arr     y coordinates
C             ue_arr    edge velocity Ue/Vinf (Karman-Tsien corrected)
C             dstr_arr  displacement thickness
C             thet_arr  momentum thickness
C             cf_arr    skin friction coefficient
C             hk_arr    kinematic shape factor Hk
C             cdis_arr  dissipation coefficient
C             ct_arr    sqrt(max shear stress coefficient)
C             n_out     number of valid entries written
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
      INTEGER nbuf, n_out
      REAL    s_arr(nbuf),    x_arr(nbuf),    y_arr(nbuf)
      REAL    ue_arr(nbuf),   dstr_arr(nbuf), thet_arr(nbuf)
      REAL    cf_arr(nbuf),   hk_arr(nbuf)
      REAL    cdis_arr(nbuf), ct_arr(nbuf)
Cf2py intent(in)  :: nbuf
Cf2py intent(out) :: s_arr, x_arr, y_arr, ue_arr, dstr_arr, thet_arr
Cf2py intent(out) :: cf_arr, hk_arr, cdis_arr, ct_arr, n_out
C
      INTEGER i, j, k, IS, IBL, nw_out
      REAL    DS, TH, CF, H, UE, UI, AMSQ, HK, CDIS, CT, DUMMY
C
C---- update compressibility parameters (HSTINV, TKLAM, etc.)
      CALL COMSET
C
C---- airfoil nodes (indices 1:N)
      n_out = MIN(N, nbuf)
C
      DO 10 i = 1, n_out
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
          CF   = TAU(IBL,IS) / (0.5*QINF**2)
          IF(TH .EQ. 0.0) THEN
            H = 1.0
          ELSE
            H = DS/TH
          ENDIF
          CDIS = DIS(IBL,IS)  / QINF**3
          CT   = CTAU(IBL,IS)
        ELSE
          DS   = 0.0
          TH   = 0.0
          CF   = 0.0
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
        s_arr(i)    = S(i)
        x_arr(i)    = X(i)
        y_arr(i)    = Y(i)
        ue_arr(i)   = UE
        dstr_arr(i) = DS
        thet_arr(i) = TH
        cf_arr(i)   = CF
        hk_arr(i)   = HK
        cdis_arr(i) = CDIS
        ct_arr(i)   = CT
  10  CONTINUE
C
C---- wake nodes (indices N+1:N+NW) if wake solution exists
      IF(LWAKE) THEN
        IS     = 2
        nw_out = MIN(NW, nbuf - n_out)
        DO 20 j = 1, nw_out
          i   = N + j
          k   = n_out + j
          IBL = IBLTE(2) + j
          DS  = DSTR(IBL,IS)
          TH  = THET(IBL,IS)
          CF  = 0.0
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
          s_arr(k)    = S(i)
          x_arr(k)    = X(i)
          y_arr(k)    = Y(i)
          ue_arr(k)   = UE
          dstr_arr(k) = DS
          thet_arr(k) = TH
          cf_arr(k)   = CF
          hk_arr(k)   = HK
          cdis_arr(k) = CDIS
          ct_arr(k)   = CT
  20    CONTINUE
        n_out = n_out + nw_out
      ENDIF
C
      END SUBROUTINE GET_BL

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
      REAL Xcloc, Ycloc, DegDef
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
      REAL chord_out, area_out, radle_out
      REAL angbte_out, thick_out, cambr_out
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
      REAL    amax_out
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
      REAL area_out, slen_out, xc_out, yc_out
      REAL aixx_out, aixxt_out, aiyy_out, aiyyt_out
      REAL aj_out, ajt_out
Cf2py intent(out) :: area_out, slen_out, xc_out, yc_out
Cf2py intent(out) :: aixx_out, aixxt_out, aiyy_out, aiyyt_out
Cf2py intent(out) :: aj_out, ajt_out
C
      REAL PEX
      REAL XMIN, XMAX, XEXINT
      REAL YMIN, YMAX, YEXINT
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
      REAL atol_in, x1_in, x2_in
Cf2py intent(in) :: atol_in, x1_in, x2_in
C
      INTEGER NNEW
      REAL XNEW(2*IBX), YNEW(2*IBX)
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
      REAL X2_in(N2_in), Y2_in(N2_in)
      REAL frac_in
Cf2py intent(in) :: X2_in, Y2_in, N2_in, frac_in
C
      REAL X2(IBX), Y2(IBX), XP2(IBX), YP2(IBX), S2(IBX)
      REAL SLE2
      REAL XOUT(IBX), YOUT(IBX)
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
