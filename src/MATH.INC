C
C----- MATH.INC : shared mathematical constants
C
C   PI    = pi
C   HOPI  = 1/(2*pi)   ("half-over-pi")
C   QOPI  = 1/(4*pi)   ("quarter-over-pi")
C   DTOR  = pi/180     (degrees to radians)
C
      REAL*8 PI, HOPI, QOPI, DTOR
      PARAMETER (PI   = 3.141592653589793238D0)
      PARAMETER (HOPI = 0.5D0  / PI)
      PARAMETER (QOPI = 0.25D0 / PI)
      PARAMETER (DTOR = PI / 180.0D0)
