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


       SUBROUTINE ROTATE(X,Y,N,ALFA)
      DIMENSION X(N), Y(N)
C
      SA = SIN(ALFA)
      CA = COS(ALFA)
CCC      XOFF = 0.25*(1.0-CA)
CCC      YOFF = 0.25*SA
      XOFF = 0.
      YOFF = 0.
      DO 8 I=1, N
        XT = X(I)
        YT = Y(I)
        X(I) = CA*XT + SA*YT + XOFF
        Y(I) = CA*YT - SA*XT + YOFF
    8 CONTINUE
C
      RETURN
      END

      SUBROUTINE TGAP(GAPNEW, DOC)
C----------------------------------
C     Used to set buffer airfoil
C     trailing edge gap.
C
C     Input:
C       GAPNEW   desired new TE gap
C       DOC      blending distance/c (0..1)
C----------------------------------
      INCLUDE 'XFOIL.INC'
C
      CALL LEFIND(SBLE,XB,XBP,YB,YBP,SB,NB)
      XBLE = SEVAL(SBLE,XB,XBP,SB,NB)
      YBLE = SEVAL(SBLE,YB,YBP,SB,NB)
      XBTE = 0.5*(XB(1)+XB(NB))
      YBTE = 0.5*(YB(1)+YB(NB))
      CHBSQ = (XBTE-XBLE)**2 + (YBTE-YBLE)**2
C
      DXN = XB(1) - XB(NB)
      DYN = YB(1) - YB(NB)
      GAP = SQRT(DXN**2 + DYN**2)
C
      IF(.NOT.LQUIET) WRITE(*,1000) GAP
 1000 FORMAT(/' Current gap =',F9.5)
C
C---- components of unit vector parallel to TE gap
      IF(GAP.GT.0.0) THEN
       DXU = DXN / GAP
       DYU = DYN / GAP
      ELSE
       DXU = -.5*(YBP(NB) - YBP(1))
       DYU = 0.5*(XBP(NB) - XBP(1))
      ENDIF
C
      DOC = MIN( MAX( DOC , 0.0 ) , 1.0 )
C
      DGAP = GAPNEW - GAP
C
C---- go over each point, changing the y-thickness appropriately
      DO 30 I=1, NB
C
C------ chord-based x/c
        XOC = (  (XB(I)-XBLE)*(XBTE-XBLE)
     &         + (YB(I)-YBLE)*(YBTE-YBLE) ) / CHBSQ
C
C------ thickness factor tails off exponentially away from trailing edge
        IF(DOC .EQ. 0.0) THEN
          TFAC = 0.0
          IF(I.EQ.1 .OR. I.EQ.NB) TFAC = 1.0
        ELSE
          ARG = MIN( (1.0-XOC)*(1.0/DOC-1.0) , 15.0 )
          TFAC = EXP(-ARG)
        ENDIF
C
        IF(SB(I).LE.SBLE) THEN
         XB(I) = XB(I) + 0.5*DGAP*XOC*TFAC*DXU
         YB(I) = YB(I) + 0.5*DGAP*XOC*TFAC*DYU
        ELSE
         XB(I) = XB(I) - 0.5*DGAP*XOC*TFAC*DXU
         YB(I) = YB(I) - 0.5*DGAP*XOC*TFAC*DYU
        ENDIF
   30 CONTINUE
      LGSAME = .FALSE.
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
      END ! TGAP

      SUBROUTINE LERAD(RFAC, DOC)
C----------------------------
C     Changes buffer airfoil
C     leading edge radius.
C
C     Input:
C       RFAC   approx. new/old LE radius scaling ratio
C       DOC    blending distance/c from LE (>= 0.001)
C----------------------------
      INCLUDE 'XFOIL.INC'
C
      DOC = MAX( DOC , 0.001 )
C
      CALL LERSCL(XB,XBP,YB,YBP,SB,NB, DOC,RFAC, W1,W2)
C
      DO 40 I=1, NB
        XB(I) = W1(I)
        YB(I) = W2(I)
   40 CONTINUE
      LGSAME = .FALSE.
C
C---- spline new coordinates
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
C---- find max curvature
      CVMAX = 0.
      DO 6 I=NB/4, (3*NB)/4
        CV = CURV(SB(I),XB,XBP,YB,YBP,SB,NB)
        CVMAX = MAX( ABS(CV) , CVMAX )
    6 CONTINUE
C
      RADIUS = 1.0/CVMAX
C
      IF(.NOT.LQUIET) WRITE(*,1000) RADIUS
 1000 FORMAT(/' New LE radius = ',F7.5)
C
      RETURN
      END ! LERAD

      SUBROUTINE TCSET(TNEW, CNEW)
C------------------------------------------------------
C     Sets the buffer airfoil to a specified maximum
C     thickness and/or camber value.
C
C     Input:
C       TNEW   new max thickness (use 999.0 to skip)
C       CNEW   new max camber    (use 999.0 to skip)
C------------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
C--- find the current buffer airfoil camber and thickness
      CALL GETCAM(XCM,YCM,NCM,XTK,YTK,NTK,
     &            XB,XBP,YB,YBP,SB,NB )
      CALL GETMAX(XCM,YCM,YCMP,NCM,CXMAX,CYMAX)
      CALL GETMAX(XTK,YTK,YTKP,NTK,TXMAX,TYMAX)
C
      IF(.NOT.LQUIET) WRITE(*,1000) 2.0*TYMAX,TXMAX, CYMAX,CXMAX
 1000 FORMAT(/' Max thickness = ',F8.4,'  at x = ',F7.3,
     &       /' Max camber    = ',F8.4,'  at x = ',F7.3/)
C
      CFAC = 1.0
      TFAC = 1.0
      IF(CYMAX.NE.0.0 .AND. CNEW.NE.999.0) CFAC = CNEW / (    CYMAX)
      IF(TYMAX.NE.0.0 .AND. TNEW.NE.999.0) TFAC = TNEW / (2.0*TYMAX)
C
C---- sanity check on scaling factors
      IF(.NOT.LQUIET) THEN
        IF(ABS(TFAC) .GT. 100.0 .OR. ABS(CFAC) .GT. 100.0) THEN
          WRITE(*,1100) TFAC, CFAC
 1100     FORMAT(/' Warning: questionable input...'
     &           /' Implied scaling factors are:', F13.2,' x thickness'
     &           /'                             ', F13.2,' x camber   ')
        ENDIF
      ENDIF
C
      CALL THKCAM(TFAC,CFAC)
C
      RETURN
      END ! TCSET

      SUBROUTINE THKCAM(TFAC,CFAC)
C---------------------------------------------------
C     Changes buffer airfoil thickness and camber
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      CALL LEFIND(SBLE,XB,XBP,YB,YBP,SB,NB)
C
C---This fails miserably with sharp LE foils, tsk,tsk,tsk HHY 4/24/01
C---- set baseline vector normal to surface at LE point
c      DXC = -DEVAL(SBLE,YB,YBP,SB,NB)
c      DYC =  DEVAL(SBLE,XB,XBP,SB,NB)
c      DSC = SQRT(DXC**2 + DYC**2)
c      DXC = DXC/DSC
c      DYC = DYC/DSC
C
C---Rational alternative 4/24/01 HHY
      XLE = SEVAL(SBLE,XB,XBP,SB,NB)
      YLE = SEVAL(SBLE,YB,YBP,SB,NB)
      XTE = 0.5*(XB(1)+XB(NB))
      YTE = 0.5*(YB(1)+YB(NB))
      CHORD = SQRT((XTE-XLE)**2 + (YTE-YLE)**2)
C---- set unit chord-line vector
      DXC = (XTE-XLE) / CHORD
      DYC = (YTE-YLE) / CHORD
C
C---- go over each point, changing the y-thickness appropriately
      DO I=1, NB
C------ coordinates of point on the opposite side with the same x value
        CALL SOPPS(SBOPP, SB(I),XB,XBP,YB,YBP,SB,NB,SBLE)
        XBOPP = SEVAL(SBOPP,XB,XBP,SB,NB)
        YBOPP = SEVAL(SBOPP,YB,YBP,SB,NB)
C
C------ set new y coordinate by changing camber & thickness appropriately
        XCAVG =        ( 0.5*(XB(I)+XBOPP)*DXC + 0.5*(YB(I)+YBOPP)*DYC )
        YCAVG = CFAC * ( 0.5*(YB(I)+YBOPP)*DXC - 0.5*(XB(I)+XBOPP)*DYC )

        XCDEL =        ( 0.5*(XB(I)-XBOPP)*DXC + 0.5*(YB(I)-YBOPP)*DYC )
        YCDEL = TFAC * ( 0.5*(YB(I)-YBOPP)*DXC - 0.5*(XB(I)-XBOPP)*DYC )
C
        W1(I) = (XCAVG+XCDEL)*DXC - (YCAVG+YCDEL)*DYC
        W2(I) = (YCAVG+YCDEL)*DXC + (XCAVG+XCDEL)*DYC
      ENDDO
C
      DO I=1, NB
        XB(I) = W1(I)
        YB(I) = W2(I)
      ENDDO
      LGSAME = .FALSE.
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
      END ! THKCAM



      SUBROUTINE HIPNT(THPNT, CHPNT)
C------------------------------------------------------
C     Changes buffer airfoil
C     thickness and/or camber highpoint.
C
C     Input:
C       THPNT  new thickness highpoint x location
C              (use 0.0 to keep current)
C       CHPNT  new camber highpoint x location
C              (use 0.0 to keep current)
C------------------------------------------------------
      INCLUDE 'XFOIL.INC'
      REAL XFN(5), YFN(5), YFNP(5), SFN(5)
C
C--- Check chordline direction (should be unrotated for camber routines)
C    to function correctly
      XLE = SEVAL(SBLE,XB,XBP,SB,NB)
      YLE = SEVAL(SBLE,YB,YBP,SB,NB)
      XTE = 0.5*(XB(1)+XB(NB))
      YTE = 0.5*(YB(1)+YB(NB))
      AROT = ATAN2(YLE-YTE,XTE-XLE) / DTOR
      IF(ABS(AROT).GT.1.0 .AND. .NOT.LQUIET) THEN
        WRITE(*,*) ' '
        WRITE(*,*) 'Warning: HIGH does not work well on rotated foils'
        WRITE(*,*) 'Current chordline angle: ',AROT
        WRITE(*,*) 'Proceeding anyway...'
      ENDIF
C
C---- find leftmost point location
      CALL XLFIND(SBL,XB,XBP,YB,YBP,SB,NB)
      XBL = SEVAL(SBL,XB,XBP,SB,NB)
      YBL = SEVAL(SBL,YB,YBP,SB,NB)
C
C---- find the current buffer airfoil camber and thickness
      CALL GETCAM(XCM,YCM,NCM,XTK,YTK,NTK,
     &            XB,XBP,YB,YBP,SB,NB )
C
C---- find the max thickness and camber
      CALL GETMAX(XCM,YCM,YCMP,NCM,CXMAX,CYMAX)
      CALL GETMAX(XTK,YTK,YTKP,NTK,TXMAX,TYMAX)
C
      IF(.NOT.LQUIET) WRITE(*,1010) 2.0*TYMAX,TXMAX, CYMAX,CXMAX
 1010 FORMAT(/' Max thickness = ',F8.4,'  at x = ',F7.3,
     &       /' Max camber    = ',F8.4,'  at x = ',F7.3/)
C
      IF (THPNT.LE.0.0) THPNT = TXMAX
      IF (CHPNT.LE.0.0) CHPNT = CXMAX
C
C--- a simple cubic mapping function is used to map x/c to move highpoints
C
C    the assumption is that a smooth function (cubic, given by the old and
C    new highpoint locations) maps the range 0-1 for x/c
C    into the range 0-1 for altered x/c distribution for the same y/c
C    thickness or camber (ie. slide the points smoothly along the x axis)
C
C--- shift thickness highpoint
      IF (THPNT .GT. 0.0) THEN
       XFN(1) = XTK(1)
       XFN(2) = TXMAX
       XFN(3) = XTK(NTK)
       YFN(1) = XTK(1)
       YFN(2) = THPNT
       YFN(3) = XTK(NTK)
       CALL SPLINA(YFN,YFNP,XFN,3)
       DO I = 1, NTK
         XTK(I) = SEVAL(XTK(I),YFN,YFNP,XFN,3)
       ENDDO
      ENDIF
C
C--- shift camber highpoint
      IF (CHPNT .GT. 0.0) THEN
       XFN(1) = XCM(1)
       XFN(2) = CXMAX
       XFN(3) = XCM(NCM)
       YFN(1) = XCM(1)
       YFN(2) = CHPNT
       YFN(3) = XCM(NCM)
       CALL SPLINA(YFN,YFNP,XFN,3)
       DO I = 1, NCM
         XCM(I) = SEVAL(XCM(I),YFN,YFNP,XFN,3)
       ENDDO
      ENDIF
C
C---- Make new airfoil from thickness and camber
C     new airfoil points are spaced to match the original
C--- HHY 4/24/01 got rid of splining vs X,Y vs S (buggy), now spline Y(X)
      CALL SEGSPL(YTK,YTKP,XTK,NTK)
      CALL SEGSPL(YCM,YCMP,XCM,NCM)
C
C
C---- for each orig. airfoil point setup new YB from camber and thickness
      DO 40 I=1, NB
C
C------ spline camber and thickness at original xb points
        YCC = SEVAL(XB(I),YCM,YCMP,XCM,NCM)
        YTT = SEVAL(XB(I),YTK,YTKP,XTK,NTK)
C
C------ set new y coordinate from new camber & thickness
        IF (SB(I) .LE. SBL) THEN
          YB(I) = YCC + YTT
         ELSE
          YB(I) = YCC - YTT
        ENDIF
C---- Add Y-offset for original leftmost (LE) point to camber
        YB(I) = YB(I) + YBL
   40 CONTINUE
      LGSAME = .FALSE.
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
      END ! HIPNT

      SUBROUTINE TCSCAL(TFAC, CFAC)
C------------------------------------------------------
C     Scales buffer airfoil thickness and camber
C     by independent factors.
C
C     Input:
C       TFAC   new/old thickness scale factor
C       CFAC   new/old camber    scale factor
C------------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
C--- find the current buffer airfoil camber and thickness
      CALL GETCAM(XCM,YCM,NCM,XTK,YTK,NTK,
     &            XB,XBP,YB,YBP,SB,NB )
      CALL GETMAX(XCM,YCM,YCMP,NCM,CXMAX,CYMAX)
      CALL GETMAX(XTK,YTK,YTKP,NTK,TXMAX,TYMAX)
C
      IF(.NOT.LQUIET) WRITE(*,1000) 2.0*TYMAX,TXMAX, CYMAX,CXMAX
 1000 FORMAT(/' Max thickness = ',F8.4,'  at x = ',F7.3,
     &       /' Max camber    = ',F8.4,'  at x = ',F7.3/)
C
      CALL THKCAM(TFAC,CFAC)
C
      RETURN
      END ! TCSCAL