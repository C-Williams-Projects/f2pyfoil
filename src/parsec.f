C***********************************************************************
C    Module:  parsec.f
C
C    Copyright (C) 2024
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
C***********************************************************************
C
C    PARSEC aerofoil parameterisation for the XFOIL f2py API.
C
C    Reference:
C      Sobieczky, H. (1999) "Parametric Airfoils and Wings",
C      Notes on Numerical Fluid Mechanics 68, pp. 71-88.
C
C    Each surface is represented by the six-term series:
C
C      z(x) = sum_{k=1}^{6}  a_k * x^{(2k-1)/2}
C
C    i.e. basis {x^{1/2}, x^{3/2}, x^{5/2}, x^{7/2}, x^{9/2}, x^{11/2}}.
C    The six coefficients per surface are determined by a 6x6 linear
C    system enforcing:
C
C      (1) TE ordinate    z(1)      = Zte +/- 0.5*dZte
C      (2) Crest ordinate z(Xup/lo) = Zup/lo
C      (3) TE slope       z'(1)     = tan((2*alpha -/+ beta)/2)
C      (4) Crest slope    z'(Xup/lo)= 0
C      (5) Crest curvature z''(Xup/lo) = Zxxup/lo
C      (6) LE radius      a_1       = +/-  sqrt(2 * r_LE)
C
C    Note on constraint (6): the local geometry near x=0 gives
C    z(x) ~ a_1 * sqrt(x), from which the curvature formula yields
C    r_LE = a_1^2 / 2, hence a_1 = sqrt(2 * r_LE).  The factor of 2
C    is required and distinct from some implementations that incorrectly
C    use a_1 = sqrt(r_LE).
C
C    Subroutines
C    -----------
C    PRSCSLV(A, B, X, SING)
C        Solve a 6x6 linear system A*x = b using Gaussian elimination
C        with partial pivoting.
C
C    PRSCOEF(rle_up, Xup, Zup, Zxxup,
C            rle_lo, Xlo, Zlo, Zxxlo,
C            Zte, dZte, alpha_te_deg, beta_te_deg,
C            a_up, a_lo, SING)
C        Build and solve the two PARSEC 6x6 systems.
C
C    PRSCSURF(PSI, NPSI, A, Y)
C        Evaluate one PARSEC surface at NPSI chord stations.
C
C    SETPARSEC(rle_up, Xup, Zup, Zxxup,
C              rle_lo, Xlo, Zlo, Zxxlo,
C              Zte, dZte, alpha_te_deg, beta_te_deg)
C        f2py API entry point.
C
C        Python signature:
C            xf.setparsec(rle_up, Xup, Zup, Zxxup,
C                         rle_lo, Xlo, Zlo, Zxxlo,
C                         Zte, dZte, alpha_te_deg, beta_te_deg)
C
C***********************************************************************


C=======================================================================
      SUBROUTINE PRSCSLV(A, B, X, SING)
C-----------------------------------------------
C     Solve the 6x6 system  A * x = b  using
C     Gaussian elimination with partial pivoting.
C
C     Input:
C       A(6,6)   coefficient matrix (destroyed on exit)
C       B(6)     right-hand side    (destroyed on exit)
C
C     Output:
C       X(6)     solution vector
C       SING     .TRUE. if the matrix is numerically singular
C                (pivot magnitude < 1e-12); X is undefined.
C
C     The 6x6 size is fixed by the PARSEC parameterisation and
C     hard-coded here to keep the routine self-contained with no
C     external LAPACK or XFOIL solver dependency.
C-----------------------------------------------
      INTEGER NP
      PARAMETER (NP = 6)
C
      REAL    A(NP,NP), B(NP), X(NP)
      LOGICAL SING
C
      INTEGER I, J, K, IPIV
      REAL    PIVOT, FAC, TMP
C
      SING = .FALSE.
C
C---- forward elimination with partial pivoting
      DO 30 K = 1, NP
C
C------ find row with largest pivot in column K
        IPIV = K
        DO 10 I = K+1, NP
          IF(ABS(A(I,K)) .GT. ABS(A(IPIV,K))) IPIV = I
   10   CONTINUE
C
C------ swap rows K and IPIV in A and B
        IF(IPIV .NE. K) THEN
          DO 15 J = K, NP
            TMP      = A(K,J)
            A(K,J)   = A(IPIV,J)
            A(IPIV,J) = TMP
   15     CONTINUE
          TMP    = B(K)
          B(K)   = B(IPIV)
          B(IPIV) = TMP
        ENDIF
C
        PIVOT = A(K,K)
        IF(ABS(PIVOT) .LT. 1.0E-12) THEN
          SING = .TRUE.
          RETURN
        ENDIF
C
C------ eliminate entries below the pivot
        DO 25 I = K+1, NP
          FAC = A(I,K) / PIVOT
          DO 20 J = K, NP
            A(I,J) = A(I,J) - FAC * A(K,J)
   20     CONTINUE
          B(I) = B(I) - FAC * B(K)
   25   CONTINUE
C
   30 CONTINUE
C
C---- back substitution
      DO 50 I = NP, 1, -1
        X(I) = B(I)
        DO 40 J = I+1, NP
          X(I) = X(I) - A(I,J) * X(J)
   40   CONTINUE
        X(I) = X(I) / A(I,I)
   50 CONTINUE
C
      RETURN
      END


C=======================================================================
      SUBROUTINE PRSCOEF(rle_up, Xup, Zup, Zxxup,
     &                   rle_lo, Xlo, Zlo, Zxxlo,
     &                   Zte, dZte, alpha_te_deg, beta_te_deg,
     &                   a_up, a_lo, SING)
C-----------------------------------------------
C     Build and solve the two PARSEC 6x6 linear systems to obtain
C     surface coefficients a_up(6) and a_lo(6).
C
C     Input: (all lengths normalised to unit chord)
C       rle_up       geometric upper-surface LE radius (>= 0)
C       Xup          upper crest abscissa  (0 < Xup < 1)
C       Zup          upper crest ordinate
C       Zxxup        upper crest curvature (d^2z/dx^2, typically < 0)
C       rle_lo       geometric lower-surface LE radius (>= 0)
C       Xlo          lower crest abscissa  (0 < Xlo < 1)
C       Zlo          lower crest ordinate  (typically < 0)
C       Zxxlo        lower crest curvature (typically > 0)
C       Zte          TE mid-ordinate  = 0.5*(z_upper_TE + z_lower_TE)
C       dZte         TE thickness gap = z_upper_TE - z_lower_TE  (>= 0)
C       alpha_te_deg TE direction angle (degrees)
C       beta_te_deg  TE wedge angle     (degrees)
C
C     Output:
C       a_up(6)      upper surface PARSEC coefficients
C       a_lo(6)      lower surface PARSEC coefficients
C       SING         .TRUE. if either system is singular
C
C     TE slope convention (Sobieczky 1999):
C       dz/dx at x=1 (upper) = tan( alpha - beta/2 )
C       dz/dx at x=1 (lower) = tan( alpha + beta/2 )
C     where  alpha = alpha_te_deg * pi/180
C            beta  = beta_te_deg  * pi/180
C
C     LE radius constraint:
C       z(x) ~ a_1 * sqrt(x)  near x = 0
C       curvature kappa = (1/2) * a_1^(-2)  =>  r_LE = a_1^2 / 2
C       so  a_1 = +sqrt(2 * r_LE)  (upper)
C           a_1 = -sqrt(2 * r_LE)  (lower)
C-----------------------------------------------
      REAL    rle_up, Xup, Zup, Zxxup
      REAL    rle_lo, Xlo, Zlo, Zxxlo
      REAL    Zte, dZte, alpha_te_deg, beta_te_deg
      REAL    a_up(6), a_lo(6)
      LOGICAL SING
C
      INTEGER NP
      PARAMETER (NP = 6)
C
      REAL    AM(NP,NP), BV(NP)
      REAL    ALPHA, BETA, SL_UP, SL_LO, PI
      REAL    XC, ORD, EXP1, EXP2, EXP3
      INTEGER K
C
      PI    = 4.0 * ATAN(1.0)
      ALPHA = alpha_te_deg * PI / 180.0
      BETA  = beta_te_deg  * PI / 180.0
C
C---- TE slope per Sobieczky 1999
C     upper: tan( alpha - beta/2 )
C     lower: tan( alpha + beta/2 )
      SL_UP = TAN(ALPHA - 0.5*BETA)
      SL_LO = TAN(ALPHA + 0.5*BETA)
C
C--------------------------------------------------------------------
C     Build common matrix structure for one surface at crest XC,
C     fill RHS BV, then solve.  Done twice (upper then lower).
C--------------------------------------------------------------------
C
C---- UPPER SURFACE -----------------------------------------------
C
      XC = Xup
C
      DO 10 K = 1, NP
C
C------ orders: (2K-1)/2  =  K-0.5  (i.e. 0.5, 1.5, ..., 5.5)
        EXP1 = REAL(K) - 0.5          ! order of basis function
        EXP2 = REAL(K) - 1.5          ! order of 1st derivative
        EXP3 = REAL(K) - 2.5          ! order of 2nd derivative
C
C       Row 1: z(1) = sum_k a_k * 1^EXP1 = sum_k a_k
        AM(1,K) = 1.0
C
C       Row 2: z(XC) = sum_k a_k * XC^EXP1
        AM(2,K) = XC**EXP1
C
C       Row 3: z'(1) = sum_k a_k * EXP1 * 1^EXP2 = sum_k a_k * EXP1
        AM(3,K) = EXP1
C
C       Row 4: z'(XC) = sum_k a_k * EXP1 * XC^EXP2
C               EXP2 = K-1.5: for K=1 this is XC^(-0.5), safe for XC > 0
        AM(4,K) = EXP1 * (XC**EXP2)
C
C       Row 5: z''(XC) = sum_k a_k * EXP1*EXP2 * XC^EXP3
C               EXP3 = K-2.5: for K=1 this is XC^(-1.5), K=2 XC^(-0.5)
        AM(5,K) = EXP1 * EXP2 * (XC**EXP3)
C
C       Row 6: a_1 = +sqrt(2*r_le)  -->  1*a_1 + 0*a_2 + ... = rhs
        AM(6,K) = 0.0
   10 CONTINUE
      AM(6,1) = 1.0
C
C---- upper RHS
      BV(1) = Zte + 0.5*dZte            ! z(1)  upper TE ordinate
      BV(2) = Zup                        ! z(Xup)
      BV(3) = SL_UP                      ! z'(1)
      BV(4) = 0.0                        ! z'(Xup) = 0  (crest)
      BV(5) = Zxxup                      ! z''(Xup)
      BV(6) = SQRT(MAX(2.0*rle_up, 0.0)) ! a_1 = +sqrt(2*r_le)
C
      CALL PRSCSLV(AM, BV, a_up, SING)
      IF(SING) RETURN
C
C---- LOWER SURFACE -----------------------------------------------
C
      XC = Xlo
C
      DO 20 K = 1, NP
C
        EXP1 = REAL(K) - 0.5
        EXP2 = REAL(K) - 1.5
        EXP3 = REAL(K) - 2.5
C
        AM(1,K) = 1.0
        AM(2,K) = XC**EXP1
        AM(3,K) = EXP1
        AM(4,K) = EXP1 * (XC**EXP2)
        AM(5,K) = EXP1 * EXP2 * (XC**EXP3)
        AM(6,K) = 0.0
   20 CONTINUE
      AM(6,1) = 1.0
C
C---- lower RHS
      BV(1) = Zte - 0.5*dZte             ! z(1)  lower TE ordinate
      BV(2) = Zlo                         ! z(Xlo)
      BV(3) = SL_LO                       ! z'(1)
      BV(4) = 0.0                         ! z'(Xlo) = 0  (crest)
      BV(5) = Zxxlo                       ! z''(Xlo)
      BV(6) = -SQRT(MAX(2.0*rle_lo, 0.0)) ! a_1 = -sqrt(2*r_le)
C
      CALL PRSCSLV(AM, BV, a_lo, SING)
C
      RETURN
      END


C=======================================================================
      SUBROUTINE PRSCSURF(PSI, NPSI, A, Y)
C-----------------------------------------------
C     Evaluate one PARSEC surface at NPSI chord stations.
C
C     The ordinate at each station is:
C
C       z(psi) = sum_{k=1}^{6}  a_k * psi^{(2k-1)/2}
C
C     Note: z(0) = 0 for any valid coefficient set since every
C     basis power (2k-1)/2 > 0.  z(1) = sum_k a_k.
C
C     Input:
C       PSI(NPSI)   chord stations, psi in [0, 1]
C       NPSI        number of stations
C       A(6)        PARSEC surface coefficients
C
C     Output:
C       Y(NPSI)     surface ordinates
C-----------------------------------------------
      INTEGER NPSI
      REAL    PSI(NPSI), A(6), Y(NPSI)
C
      INTEGER I, K
      REAL    P, EXP1, ZVAL
C
      DO 20 I = 1, NPSI
        P    = PSI(I)
        ZVAL = 0.0
C
        IF(P .LE. 0.0) THEN
C-------- all basis functions vanish at the LE (powers > 0)
          ZVAL = 0.0
        ELSE
          DO 10 K = 1, 6
            EXP1 = REAL(K) - 0.5     ! (2k-1)/2
            ZVAL = ZVAL + A(K) * (P**EXP1)
   10     CONTINUE
        ENDIF
C
        Y(I) = ZVAL
   20 CONTINUE
C
      RETURN
      END


C=======================================================================
      SUBROUTINE SETPARSEC(rle_up, Xup, Zup, Zxxup,
     &                     rle_lo, Xlo, Zlo, Zxxlo,
     &                     Zte, dZte, alpha_te_deg, beta_te_deg)
C-----------------------------------------------
C     f2py API entry point.
C
C     Generates a PARSEC aerofoil from the twelve geometric parameters,
C     loads the coordinates into the XFOIL buffer arrays XB/YB/NB in
C     Selig order (TE-upper -> LE -> TE-lower), and calls PANGEN to
C     create the panelled geometry ready for analysis.
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
C
C     Python call:
C       xf.setparsec(rle_up, Xup, Zup, Zxxup,
C                    rle_lo, Xlo, Zlo, Zxxlo,
C                    Zte, dZte, alpha_te_deg, beta_te_deg)
C
C     Notes:
C       - XFOIL must be initialised (xf.Initialise()) before calling.
C       - All dependent solution flags are reset after geometry loading,
C         consistent with the behaviour of REPANEL in api.f.
C       - If the 6x6 system is singular (e.g. Xup == Xlo == 0) a
C         warning is printed and no geometry is loaded.
C-----------------------------------------------
      INCLUDE 'XFOIL.INC'
      INCLUDE 'QUIET.INC'
C
      REAL    rle_up, Xup, Zup, Zxxup
      REAL    rle_lo, Xlo, Zlo, Zxxlo
      REAL    Zte, dZte, alpha_te_deg, beta_te_deg
C
Cf2py intent(in) :: rle_up, Xup, Zup, Zxxup
Cf2py intent(in) :: rle_lo, Xlo, Zlo, Zxxlo
Cf2py intent(in) :: Zte, dZte, alpha_te_deg, beta_te_deg
C
C---- number of psi stations used to seed the buffer airfoil.
C     101 gives smooth curvature for PANGEN to work with.
C     Must be <= IBX/2  (IBX = 2*IQX = 1000 by default).
      INTEGER NPSI
      PARAMETER (NPSI = 101)
C
      REAL    PSI(NPSI), YU(NPSI), YL(NPSI)
      REAL    a_up(6), a_lo(6)
      LOGICAL SING
      INTEGER I, IB
      REAL    PI
C
      PI = 4.0 * ATAN(1.0)
C
C---- solve the two PARSEC 6x6 systems
      CALL PRSCOEF(rle_up, Xup, Zup, Zxxup,
     &             rle_lo, Xlo, Zlo, Zxxlo,
     &             Zte, dZte, alpha_te_deg, beta_te_deg,
     &             a_up, a_lo, SING)
C
      IF(SING) THEN
        IF(.NOT.LQUIET)
     &    WRITE(*,*) 'SETPARSEC: singular 6x6 system -- '//
     &               'check crest abscissae and LE radii.'
        RETURN
      ENDIF
C
C---- cosine-spaced psi stations on [0, 1]
C     Provides finer clustering near LE and TE, where curvature is
C     highest, giving PANGEN a smooth spline to re-panel from.
      DO 10 I = 1, NPSI
        PSI(I) = 0.5 * (1.0 - COS(PI * REAL(I-1) / REAL(NPSI-1)))
   10 CONTINUE
C
C---- evaluate upper and lower surfaces
      CALL PRSCSURF(PSI, NPSI, a_up, YU)
      CALL PRSCSURF(PSI, NPSI, a_lo, YL)
C
C---- assemble Selig-ordered buffer airfoil
C
C       Segment 1: upper surface TE -> LE   (I = NPSI downto 1)
C       Segment 2: lower surface LE -> TE   (I = 2 upto NPSI)
C       Total points: 2*NPSI - 1  (LE shared, not duplicated)
C
      IB = 0
C
      DO 20 I = NPSI, 1, -1
        IB      = IB + 1
        XB(IB)  = PSI(I)
        YB(IB)  = YU(I)
   20 CONTINUE
C
      DO 30 I = 2, NPSI
        IB      = IB + 1
        XB(IB)  = PSI(I)
        YB(IB)  = YL(I)
   30 CONTINUE
C
      NB = IB
C
C---- generate panels from the buffer airfoil.
C     .FALSE. argument suppresses XFOIL screen output.
      CALL PANGEN(.FALSE.)
C
C---- invalidate all dependent solution quantities,
C     consistent with REPANEL in api.f
      LGAMU  = .FALSE.
      LQAIJ  = .FALSE.
      LWAKE  = .FALSE.
      LBLINI = .FALSE.
      LIPAN  = .FALSE.
      LVCONV = .FALSE.
C
      RETURN
      END
