C***********************************************************************
C    Module:  xgdes.f
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
      SUBROUTINE ABCOPY(LCONF)
      INCLUDE 'XFOIL.INC'
      LOGICAL LCONF
C
      IF(NB.LE.1) THEN
       IF(.NOT.LQUIET)
     &  WRITE(*,*) 'ABCOPY: Buffer airfoil not available.'
       RETURN
      ELSEIF(NB.GT.IQX-5) THEN
       IF(.NOT.LQUIET) THEN
        WRITE(*,*) 'Maximum number of panel nodes  : ',IQX-5
        WRITE(*,*) 'Number of buffer airfoil points: ',NB
        WRITE(*,*) 'Current airfoil cannot be set.'
        WRITE(*,*) 'Try executing PANE at Top Level instead.'
       ENDIF
       RETURN
      ENDIF
      IF(N.NE.NB) LBLINI = .FALSE.
C
      N = NB
      DO 101 I=1, N
        X(I) = XB(I)
        Y(I) = YB(I)
  101 CONTINUE
      LGSAME = .TRUE.
C
      IF(LBFLAP) THEN
       XOF = XBF
       YOF = YBF
       LFLAP = .TRUE.
      ENDIF
C
C---- strip out doubled points
      I = 1
 102  CONTINUE
      I = I+1
      IF(X(I-1).EQ.X(I) .AND. Y(I-1).EQ.Y(I)) THEN
        DO 104 J=I, N-1
          X(J) = X(J+1)
          Y(J) = Y(J+1)
 104    CONTINUE
        N = N-1
      ENDIF
      IF(I.LT.N) GO TO 102
C
      CALL SCALC(X,Y,S,N)
      CALL SEGSPL(X,XP,S,N)
      CALL SEGSPL(Y,YP,S,N)
      CALL NCALC(X,Y,S,N,NX,NY)
      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
      CHORD  = SQRT( (XTE-XLE)**2 + (YTE-YLE)**2 )
      CALL TECALC
      CALL APCALC
C
      LGAMU = .FALSE.
      LQINU = .FALSE.
      LWAKE = .FALSE.
      LQAIJ = .FALSE.
      LADIJ = .FALSE.
      LWDIJ = .FALSE.
      LIPAN = .FALSE.
      LVCONV = .FALSE.
      LSCINI = .FALSE.
CCC      LBLINI = .FALSE.
C
C      IF(LCONF) WRITE(*,1200) N
 1200 FORMAT(/' Current airfoil nodes set from buffer airfoil nodes (',
     &        I4,' )')
C
      RETURN
      END ! ABCOPY

      SUBROUTINE SCLXY(SCL, XORG, YORG)
C-----------------------------------------------------------
C     Scale buffer airfoil by factor SCL about point
C     (XORG, YORG).
C-----------------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      DO 10 I=1, NB
        XB(I) = SCL*(XB(I) - XORG) + XORG 
        YB(I) = SCL*(YB(I) - YORG) + YORG 
   10 CONTINUE
C
      CALL SCALC(XB,YB,SB,NB)
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
C
      CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
     &            SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
     &            EI11BA,EI22BA,APX1BA,APX2BA,
     &            EI11BT,EI22BT,APX1BT,APX2BT,
     &            THICKB,CAMBRB )
C
      RETURN
      END ! SCLXY

      SUBROUTINE FLAP(XBF_IN, YBF_IN, DDEF_IN)
C----------------------------------------------------
C     Modifies buffer airfoil for a deflected flap.
C     Points may be added/subtracted in the flap
C     break vicinity to clean things up.
C
C     Inputs:
C       XBF_IN   flap hinge x coordinate
C       YBF_IN   flap hinge y coordinate  
C       DDEF_IN  flap deflection in degrees (+ down)
C----------------------------------------------------
      INCLUDE 'XFOIL.INC'
      LOGICAL LCHANGE
      LOGICAL INSID
      LOGICAL INSIDE
      LOGICAL LT1NEW,LT2NEW,LB1NEW,LB2NEW
C
      REAL XBF_IN, YBF_IN, DDEF_IN
C
C---- set hinge location into global flap variables
      XBF = XBF_IN
      YBF = YBF_IN
C
      CALL GETXYF(XB,XBP,YB,YBP,SB,NB, TOPS,BOTS,XBF,YBF)
      INSID = INSIDE(XB,YB,NB,XBF,YBF)
C
C      WRITE(*,1050) XBF, YBF
 1050 FORMAT(/' Flap hinge: x,y =', 2F9.5 )
C
      DDEF = DDEF_IN
      RDEF = DDEF*PI/180.0
      IF(RDEF .EQ. 0.0) RETURN
C
      IF(INSID) THEN
        ATOP = MAX( 0.0 , -RDEF )
        ABOT = MAX( 0.0 ,  RDEF )
      ELSE
        CHX = DEVAL(BOTS,XB,XBP,SB,NB) - DEVAL(TOPS,XB,XBP,SB,NB)
        CHY = DEVAL(BOTS,YB,YBP,SB,NB) - DEVAL(TOPS,YB,YBP,SB,NB)
        FVX = SEVAL(BOTS,XB,XBP,SB,NB) + SEVAL(TOPS,XB,XBP,SB,NB)
        FVY = SEVAL(BOTS,YB,YBP,SB,NB) + SEVAL(TOPS,YB,YBP,SB,NB)
        CRSP = CHX*(YBF-0.5*FVY) - CHY*(XBF-0.5*FVX)
        IF(CRSP .GT. 0.0) THEN
C-------- flap hinge is above airfoil
          ATOP = MAX( 0.0 ,  RDEF )
          ABOT = MAX( 0.0 ,  RDEF )
        ELSE
C-------- flap hinge is below airfoil
          ATOP = MAX( 0.0 , -RDEF )
          ABOT = MAX( 0.0 , -RDEF )
        ENDIF
      ENDIF
C
C---- find upper and lower surface break arc length values...
      CALL SSS(TOPS,ST1,ST2,ATOP,XBF,YBF,XB,XBP,YB,YBP,SB,NB,1)
      CALL SSS(BOTS,SB1,SB2,ABOT,XBF,YBF,XB,XBP,YB,YBP,SB,NB,2)
C
C---- ... and x,y coordinates
      XT1 = SEVAL(ST1,XB,XBP,SB,NB)
      YT1 = SEVAL(ST1,YB,YBP,SB,NB)
      XT2 = SEVAL(ST2,XB,XBP,SB,NB)
      YT2 = SEVAL(ST2,YB,YBP,SB,NB)
      XB1 = SEVAL(SB1,XB,XBP,SB,NB)
      YB1 = SEVAL(SB1,YB,YBP,SB,NB)
      XB2 = SEVAL(SB2,XB,XBP,SB,NB)
      YB2 = SEVAL(SB2,YB,YBP,SB,NB)
C
C      WRITE(*,1100) XT1, YT1, XT2, YT2,
C     &              XB1, YB1, XB2, YB2
 1100 FORMAT(/' Top breaks: x,y =  ', 2F9.5, 4X, 2F9.5
     &       /' Bot breaks: x,y =  ', 2F9.5, 4X, 2F9.5)
C
C---- find points adjacent to breaks
      DO 5 I=1, NB-1
        IF(SB(I).LE.ST1 .AND. SB(I+1).GT.ST1) IT1 = I+1
        IF(SB(I).LT.ST2 .AND. SB(I+1).GE.ST2) IT2 = I
        IF(SB(I).LE.SB1 .AND. SB(I+1).GT.SB1) IB1 = I
        IF(SB(I).LT.SB2 .AND. SB(I+1).GE.SB2) IB2 = I+1
    5 CONTINUE
C
      DSAVG = (SB(NB)-SB(1))/FLOAT(NB-1)
C
C---- smallest fraction of s increments i+1 and i+2 away from break point
      SFRAC = 0.33333
C
      IF(ATOP .NE. 0.0) THEN
        ST1P = ST1 + SFRAC*(SB(IT1  )-ST1)
        ST1Q = ST1 + SFRAC*(SB(IT1+1)-ST1)
        IF(SB(IT1) .LT. ST1Q) THEN
C-------- simply move adjacent point to ideal SFRAC location
          XT1NEW = SEVAL(ST1Q,XB,XBP,SB,NB)
          YT1NEW = SEVAL(ST1Q,YB,YBP,SB,NB)
          LT1NEW = .FALSE.
        ELSE
C-------- make new point at SFRAC location
          XT1NEW = SEVAL(ST1P,XB,XBP,SB,NB)
          YT1NEW = SEVAL(ST1P,YB,YBP,SB,NB)
          LT1NEW = .TRUE.
        ENDIF
C
        ST2P = ST2 + SFRAC*(SB(IT2 )-ST2)
        IT2Q = MAX(IT2-1,1)
        ST2Q = ST2 + SFRAC*(SB(IT2Q)-ST2)
        IF(SB(IT2) .GT. ST2Q) THEN
C-------- simply move adjacent point
          XT2NEW = SEVAL(ST2Q,XB,XBP,SB,NB)
          YT2NEW = SEVAL(ST2Q,YB,YBP,SB,NB)
          LT2NEW = .FALSE.
        ELSE
C-------- make new point
          XT2NEW = SEVAL(ST2P,XB,XBP,SB,NB)
          YT2NEW = SEVAL(ST2P,YB,YBP,SB,NB)
          LT2NEW = .TRUE.
        ENDIF
      ENDIF
C
      IF(ABOT .NE. 0.0) THEN
        SB1P = SB1 + SFRAC*(SB(IB1  )-SB1)
        SB1Q = SB1 + SFRAC*(SB(IB1-1)-SB1)
        IF(SB(IB1) .GT. SB1Q) THEN
C-------- simply move adjacent point
          XB1NEW = SEVAL(SB1Q,XB,XBP,SB,NB)
          YB1NEW = SEVAL(SB1Q,YB,YBP,SB,NB)
          LB1NEW = .FALSE.
        ELSE
C-------- make new point
          XB1NEW = SEVAL(SB1P,XB,XBP,SB,NB)
          YB1NEW = SEVAL(SB1P,YB,YBP,SB,NB)
          LB1NEW = .TRUE.
        ENDIF
C
        SB2P = SB2 + SFRAC*(SB(IB2 )-SB2)
        IB2Q = MIN(IB2+1,NB)
        SB2Q = SB2 + SFRAC*(SB(IB2Q)-SB2)
        IF(SB(IB2) .LT. SB2Q) THEN
C-------- simply move adjacent point
          XB2NEW = SEVAL(SB2Q,XB,XBP,SB,NB)
          YB2NEW = SEVAL(SB2Q,YB,YBP,SB,NB)
          LB2NEW = .FALSE.
        ELSE
C-------- make new point
          XB2NEW = SEVAL(SB2P,XB,XBP,SB,NB)
          YB2NEW = SEVAL(SB2P,YB,YBP,SB,NB)
          LB2NEW = .TRUE.
        ENDIF
      ENDIF
C
      SIND = SIN(RDEF)
      COSD = COS(RDEF)
C
C---- rotate flap points about the hinge point (XBF,YBF)
      DO 10 I=1, NB
        IF(I.GE.IT1 .AND. I.LE.IB1) GO TO 10
C
        XBAR = XB(I) - XBF
        YBAR = YB(I) - YBF
C
        XB(I) = XBF  +  XBAR*COSD  +  YBAR*SIND
        YB(I) = YBF  -  XBAR*SIND  +  YBAR*COSD
   10 CONTINUE
C
      IDIF = IT1-IT2-1
      IF(IDIF.GT.0) THEN
C----- delete points on upper airfoil surface which "disappeared".
       NB  = NB -IDIF
       IT1 = IT1-IDIF
       IB1 = IB1-IDIF
       IB2 = IB2-IDIF
       DO 21 I=IT2+1, NB
         SB(I) = SB(I+IDIF)
         XB(I) = XB(I+IDIF)
         YB(I) = YB(I+IDIF)
   21  CONTINUE
      ENDIF
C
      IDIF = IB2-IB1-1
      IF(IDIF.GT.0) THEN
C----- delete points on lower airfoil surface which "disappeared".
       NB  = NB -IDIF
       IB2 = IB2-IDIF
       DO 22 I=IB1+1, NB
         SB(I) = SB(I+IDIF)
         XB(I) = XB(I+IDIF)
         YB(I) = YB(I+IDIF)
   22  CONTINUE
      ENDIF
C
      IF(ATOP .EQ. 0.0) THEN
C
C------ arc length of newly created surface on top of airfoil
        DSNEW = ABS(RDEF)*SQRT((XT1-XBF)**2 + (YT1-YBF)**2)
C
C------ number of points to be added to define newly created surface
        NPADD = INT(1.5*DSNEW/DSAVG + 1.0)
C
C------ skip everything if no points are to be added
        IF(NPADD.EQ.0) GO TO 35
C
C------ increase coordinate array length to make room for the new point(s)
        NB  = NB +NPADD
        IT1 = IT1+NPADD
        IB1 = IB1+NPADD
        IB2 = IB2+NPADD
        DO 30 I=NB, IT1, -1
          XB(I) = XB(I-NPADD)
          YB(I) = YB(I-NPADD)
   30   CONTINUE
C
C------ add new points along the new surface circular arc segment
        DANG = RDEF / FLOAT(NPADD)
        XBAR = XT1 - XBF
        YBAR = YT1 - YBF
        DO 31 IP=1, NPADD
          ANG = DANG*(FLOAT(IP) - 0.5)
          CA = COS(ANG)
          SA = SIN(ANG)
C
          XB(IT1-IP) = XBF  +  XBAR*CA + YBAR*SA
          YB(IT1-IP) = YBF  -  XBAR*SA + YBAR*CA
   31   CONTINUE
C
      ELSE
C
C------ set point in the corner and possibly two adjacent points
        NPADD = 1
        IF(LT2NEW) NPADD = NPADD+1
        IF(LT1NEW) NPADD = NPADD+1
C
        NB  = NB +NPADD
        IT1 = IT1+NPADD
        IB1 = IB1+NPADD
        IB2 = IB2+NPADD
        DO 33 I=NB, IT1, -1
          XB(I) = XB(I-NPADD)
          YB(I) = YB(I-NPADD)
   33   CONTINUE
C
        IF(LT1NEW) THEN
         XB(IT1-1) = XT1NEW
         YB(IT1-1) = YT1NEW
         XB(IT1-2) = XT1
         YB(IT1-2) = YT1
        ELSE
         XB(IT1  ) = XT1NEW
         YB(IT1  ) = YT1NEW
         XB(IT1-1) = XT1
         YB(IT1-1) = YT1
        ENDIF
C
        XBAR = XT2NEW - XBF
        YBAR = YT2NEW - YBF
        IF(LT2NEW) THEN
          XB(IT2+1) = XBF  +  XBAR*COSD + YBAR*SIND
          YB(IT2+1) = YBF  -  XBAR*SIND + YBAR*COSD
        ELSE
          XB(IT2  ) = XBF  +  XBAR*COSD + YBAR*SIND
          YB(IT2  ) = YBF  -  XBAR*SIND + YBAR*COSD
        ENDIF
C
      ENDIF
   35 CONTINUE
C
      IF(ABOT .EQ. 0.0) THEN
C
C------ arc length of newly created surface on bottom of airfoil
        DSNEW = ABS(RDEF)*SQRT((XB1-XBF)**2 + (YB1-YBF)**2)
C
C------ number of points to be added to define newly created surface
        NPADD = INT(1.5*DSNEW/DSAVG + 1.0)
C
C------ skip everything if no points are to be added
        IF(NPADD.EQ.0) GO TO 45
C
C------ increase coordinate array length to make room for the new point(s)
        NB  = NB +NPADD
        IB2 = IB2+NPADD
        DO 40 I=NB, IB2, -1
          XB(I) = XB(I-NPADD)
          YB(I) = YB(I-NPADD)
   40   CONTINUE
C
C------ add new points along the new surface circular arc segment
        DANG = RDEF / FLOAT(NPADD)
        XBAR = XB1 - XBF
        YBAR = YB1 - YBF
        DO 41 IP=1, NPADD
          ANG = DANG*(FLOAT(IP) - 0.5)
          CA = COS(ANG)
          SA = SIN(ANG)
C
          XB(IB1+IP) = XBF  +  XBAR*CA + YBAR*SA
          YB(IB1+IP) = YBF  -  XBAR*SA + YBAR*CA
   41   CONTINUE
C
      ELSE
C
C------ set point in the corner and possibly two adjacent points
        NPADD = 1
        IF(LB2NEW) NPADD = NPADD+1
        IF(LB1NEW) NPADD = NPADD+1
C
        NB  = NB +NPADD
        IB2 = IB2+NPADD
        DO 43 I=NB, IB2, -1
          XB(I) = XB(I-NPADD)
          YB(I) = YB(I-NPADD)
   43   CONTINUE
C
        IF(LB1NEW) THEN
         XB(IB1+1) = XB1NEW
         YB(IB1+1) = YB1NEW
         XB(IB1+2) = XB1
         YB(IB1+2) = YB1
        ELSE
         XB(IB1  ) = XB1NEW
         YB(IB1  ) = YB1NEW
         XB(IB1+1) = XB1
         YB(IB1+1) = YB1
        ENDIF
C
        XBAR = XB2NEW - XBF
        YBAR = YB2NEW - YBF
        IF(LB2NEW) THEN
          XB(IB2-1) = XBF  +  XBAR*COSD + YBAR*SIND
          YB(IB2-1) = YBF  -  XBAR*SIND + YBAR*COSD
        ELSE
          XB(IB2  ) = XBF  +  XBAR*COSD + YBAR*SIND
          YB(IB2  ) = YBF  -  XBAR*SIND + YBAR*COSD
        ENDIF
C
      ENDIF
   45 CONTINUE
C
C---- check new geometry for splinter segments
      STOL = 0.2
      CALL SCHECK(XB,YB,NB, STOL, LCHANGE)
C
C---- spline new geometry
      CALL SCALC(XB,YB,SB,NB)
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
C
      CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
     &            SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
     &            EI11BA,EI22BA,APX1BA,APX2BA,
     &            EI11BT,EI22BT,APX1BT,APX2BT,
     &            THICKB,CAMBRB )
C
      LBFLAP = .TRUE.
C
      IF(LGSYM) THEN
       LGSYM = .FALSE.
      ENDIF
C
      LGSAME = .FALSE.

C---- panel and copy to current airfoil
      CALL PANGEN(.FALSE.)
      CALL ABCOPY(.TRUE.)
C
      RETURN
      END ! FLAP

      SUBROUTINE GETXYF(X,XP,Y,YP,S,N, TOPS,BOTS,XF,YF)
      DIMENSION X(N),XP(N),Y(N),YP(N),S(N)
C
C---- find top and bottom arc length values at hinge x location
      TOPS = S(1) + (X(1) - XF)
      BOTS = S(N) - (X(N) - XF)
      CALL SINVRT(TOPS,XF,X,XP,S,N)
      CALL SINVRT(BOTS,XF,X,XP,S,N)
C
      TOPY = SEVAL(TOPS,Y,YP,S,N)
      BOTY = SEVAL(BOTS,Y,YP,S,N)
C
C      WRITE(*,1000) TOPY, BOTY
 1000 FORMAT(/'  Top    surface:  y =', F8.4
     &       /'  Bottom surface:  y =', F8.4)
C
      RETURN
      END ! GETXYF



