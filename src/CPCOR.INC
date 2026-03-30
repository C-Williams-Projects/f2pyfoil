C
C---- CPCOR.INC
C     Pre-computed compressibility-correction coefficients.
C     Set by COMSET according to ICCOR.
C
C     The Cp formula used throughout is:
C         DEN    = BETA + CPCBFAC * CPINC
C         CP     = CPINC / DEN
C
C     Correction     CPCBFAC
C     ──────────     ────────────────────────────────
C     PG  (0)        0
C     KT  (1)        0.5*Minf^2 / (1+beta)
C     Laitone (2)    0.5*Minf^2 * (1 + 0.5*(gamma-1)*Minf^2) / beta
C
C     CPCBFQ  = d(CPCBFAC)/d(Minf^2)   -- used in CLCALC sensitivities
C
      REAL CPCBFAC, CPCBFQ
      COMMON /CPCORR/ CPCBFAC, CPCBFQ