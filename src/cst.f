C***********************************************************************
C    Module:  cst.f
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
C    CST (Class-Shape Transformation) airfoil parameterisation
C    for the XFOIL f2py API.
C
C    Reference:
C      Kulfan, B.M. (2008) "Universal Parametric Geometry
C      Representation Method", J. Aircraft 45(1), pp. 142-158.
C
C    Subroutines
C    -----------
C    CSTBINO(N, K)
C        Binomial coefficient C(N,K) returned as REAL.
C
C    CSTSURF(PSI, NPSI, A, NA, N1, N2, TE_HALF, LEM, Y)
C        Evaluate one CST surface over NPSI psi stations.
C
C    CST(LEM, TE, Au, NAu, Al, NAl, N1, N2)
C        Computation routine. Builds a CST airfoil, loads it into
C        the XFOIL buffer arrays XB/YB/NB, and calls PANGEN.
C
C        Called by the SETCST wrapper in api.f which carries the
C        f2py directives and forms the Python-facing entry point.
C
C***********************************************************************


C=======================================================================
      REAL FUNCTION CSTBINO(N, K)
C-----------------------------------------------
C     Binomial coefficient C(N,K) returned as REAL.
C
C     Uses the multiplicative formula evaluated from
C     the smaller end to minimise rounding for the
C     moderate degrees (typically <= 20) used in CST.
C
C     Returns 0.0 for K < 0 or K > N.
C-----------------------------------------------
      INTEGER N, K, J, KK
      REAL    VAL
C
      IF(K.LT.0 .OR. K.GT.N) THEN
        CSTBINO = 0.0
        RETURN
      ENDIF
      IF(K.EQ.0 .OR. K.EQ.N) THEN
        CSTBINO = 1.0
        RETURN
      ENDIF
C
C---- exploit symmetry C(N,K) = C(N,N-K) to halve the work
      KK  = MIN(K, N - K)
      VAL = 1.0
      DO 10 J = 1, KK
        VAL = VAL * REAL(N - KK + J) / REAL(J)
   10 CONTINUE
      CSTBINO = VAL
C
      RETURN
      END


C=======================================================================
      SUBROUTINE CSTSURF(PSI, NPSI, A, NA, N1, N2,
     &                   TE_HALF, LEM, Y)
C-----------------------------------------------
C     Evaluate one CST surface at NPSI chord stations.
C
C     The ordinate at each station is:
C
C       y(psi) = C(psi)*S(psi) + TE_HALF*psi + LEM*phi(psi)
C
C     Class function (Kulfan Eq. 2):
C       C(psi) = psi^N1 * (1 - psi)^N2
C
C     Shape function (Bernstein series, Kulfan Eq. 3):
C       S(psi) = sum_{i=0}^{n}  A(i+1) * B_{i,n}(psi)
C       B_{i,n}(psi) = C(n,i) * psi^i * (1-psi)^(n-i)
C       n = NA - 1  (polynomial degree)
C
C     Trailing-edge thickness term:
C       TE_HALF = +0.5*TE  for upper surface
C               = -0.5*TE  for lower surface
C
C     Leading-edge modification (LEM) basis (Kulfan Eq. 8):
C       phi(psi) = psi * (1 - psi)^(n + 0.5)
C       vanishes at both ends so it does not affect LE/TE conditions.
C       Pass LEM = 0.0 to disable.
C
C     Input:
C       PSI(NPSI)   normalised chord stations, psi in [0, 1]
C       NPSI        number of stations
C       A(NA)       CST shape coefficients  A_0 ... A_{NA-1}
C       NA          number of coefficients  (degree n = NA - 1)
C       N1, N2      class-function exponents
C       TE_HALF     signed half trailing-edge gap
C       LEM         leading-edge modification coefficient
C
C     Output:
C       Y(NPSI)     surface ordinates
C-----------------------------------------------
      INTEGER NPSI, NA
      REAL    PSI(NPSI), A(NA), Y(NPSI)
      REAL    N1, N2, TE_HALF, LEM
C
      REAL    CSTBINO
      INTEGER I, J, ND
      REAL    P, PM, CFUNC, SVAL, BVAL, PHIVAL
C
      ND = NA - 1
C
      DO 30 I = 1, NPSI
        P  = PSI(I)
        PM = 1.0 - P
C
C------ class function C(psi) = psi^N1 * (1-psi)^N2
C       Guard against 0^(negative or zero exponent) explicitly
        IF(P .EQ. 0.0) THEN
          IF(N1 .GT. 0.0) THEN
            CFUNC = 0.0
          ELSE
            CFUNC = PM**N2
          ENDIF
        ELSEIF(PM .EQ. 0.0) THEN
          IF(N2 .GT. 0.0) THEN
            CFUNC = 0.0
          ELSE
            CFUNC = P**N1
          ENDIF
        ELSE
          CFUNC = (P**N1) * (PM**N2)
        ENDIF
C
C------ shape function S(psi) via Bernstein basis polynomials
        SVAL = 0.0
        DO 20 J = 0, ND
C
C-------- B_{J,ND}(P) = C(ND,J) * P^J * (1-P)^(ND-J)
C         Handle end cases to avoid 0^0 arithmetic
          IF(J .EQ. 0) THEN
            BVAL = PM**(ND)
          ELSEIF(J .EQ. ND) THEN
            BVAL = P**ND
          ELSE
            BVAL = CSTBINO(ND, J) * (P**J) * (PM**(ND - J))
          ENDIF
          SVAL = SVAL + A(J + 1) * BVAL
   20   CONTINUE
C
C------ LEM basis  phi(psi) = psi * (1-psi)^(n + 0.5)
        IF(P .EQ. 0.0 .OR. PM .EQ. 0.0) THEN
          PHIVAL = 0.0
        ELSE
          PHIVAL = P * (PM**(REAL(ND) + 0.5))
        ENDIF
C
        Y(I) = CFUNC * SVAL + TE_HALF * P + LEM * PHIVAL
C
   30 CONTINUE
C
      RETURN
      END


C=======================================================================
      SUBROUTINE CST(LEM, TE, Au, NAu, Al, NAl, N1, N2)
C-----------------------------------------------
C     CST computation routine.
C
C     Generates a CST aerofoil from the supplied upper and lower
C     surface coefficient vectors, loads the coordinates into the
C     XFOIL buffer arrays XB/YB/NB in Selig order
C     (TE-upper -> LE -> TE-lower), and calls PANGEN to create the
C     panelled geometry ready for analysis.
C
C     Because PANGEN re-distributes the panels, the internal cosine
C     spacing used here is only a smooth initial seed; panel count
C     and clustering are controlled by the existing XFOIL panelling
C     parameters (NPAN, CVPAR, etc.) set via REPANEL.
C
C     Input:
C       LEM         leading-edge modification coefficient.
C                   Set to 0.0 to disable.
C       TE          trailing-edge thickness as a fraction of chord.
C                   Set to 0.0 for a sharp trailing edge.
C       Au(NAu)     upper surface CST coefficients (A_0 ... A_{n_u})
C       NAu         length of Au
C       Al(NAl)     lower surface CST coefficients (A_0 ... A_{n_l})
C       NAl         length of Al
C       N1, N2      class-function exponents.
C                   Canonical aerofoil values: N1=0.5, N2=1.0.
C
C     Notes:
C       - XFOIL must be initialised (xf.Initialise()) before calling.
C       - All dependent solution flags are reset after geometry loading,
C         consistent with the behaviour of REPANEL in api.f.
C       - Upper and lower coefficient vectors may have different lengths,
C         allowing asymmetric polynomial degrees on each surface.
C       - Called by SETCST in api.f (the f2py-facing wrapper).
C-----------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      INTEGER NAu, NAl
      REAL    Au(NAu), Al(NAl)
      REAL    LEM, TE, N1, N2
C
C---- number of psi stations used to seed the buffer airfoil.
C     101 gives smooth curvature for PANGEN to work with.
C     Must be <= IBX/2  (IBX = 2*IQX = 1000 by default).
      INTEGER NPSI
      PARAMETER (NPSI = 101)
C
      REAL    PSI(NPSI), YU(NPSI), YL(NPSI)
      INTEGER I, IB
      REAL    PI
C
      PI = 4.0 * ATAN(1.0)
C
C---- cosine-spaced psi stations on [0, 1]
C     Provides finer point density near LE and TE where curvature
C     is highest, giving PANGEN a good quality spline to re-panel.
      DO 10 I = 1, NPSI
        PSI(I) = 0.5 * (1.0 - COS(PI * REAL(I-1) / REAL(NPSI-1)))
   10 CONTINUE
C
C---- upper surface  (positive TE half-gap)
      CALL CSTSURF(PSI, NPSI, Au, NAu, N1, N2,
     &              0.5*TE, LEM, YU)
C
C---- lower surface  (negative TE half-gap)
      CALL CSTSURF(PSI, NPSI, Al, NAl, N1, N2,
     &             -0.5*TE, LEM, YL)
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
        IB = IB + 1
        XB(IB) = PSI(I)
        YB(IB) = YU(I)
   20 CONTINUE
C
      DO 30 I = 2, NPSI
        IB = IB + 1
        XB(IB) = PSI(I)
        YB(IB) = YL(I)
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
