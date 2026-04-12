C
C====  XFOIL code global INCLUDE file  =====
C
C------ Primary dimensioning limit parameters
C IQX   number of surface panel nodes + 6
C IWX   number of wake panel nodes
C IPX   number of Qspec(s) distributions
C ISX   number of airfoil sides
C
C------ Derived dimensioning limit parameters
C IBX   number of buffer airfoil nodes
C IMX   number of complex mapping coefficients  Cn
C IZX   number of panel nodes (airfoil + wake)
C IVX   number of nodes along BL on one side of airfoil and wake
C
C NTX   number of points in thickness/camber arrays
C
C
      PARAMETER (IQX=740, IPX=5, ISX=2)
cc      PARAMETER (IQX=1000, IWX=72, IPX=5, ISX=2)
      PARAMETER (IWX=IQX/8 + 2)
      PARAMETER (IBX=4*IQX)
      PARAMETER (IZX=IQX+IWX)
      PARAMETER (IVX=IQX/2 + IWX + 50)
C
      PARAMETER (NTX=2*IBX)
      CHARACTER*48 NAME
      CHARACTER*80 ISPARS
      LOGICAL OK,LIMAGE,
     &        LGAMU,LQINU,SHARP,LVISC,LALFA,LWAKE,
     &        LBLINI,LIPAN,LQAIJ,LADIJ,LWDIJ,LCPXX,LQVDES,LQREFL,
     &        LQSPEC,LVCONV,LCPREF,LCLOCK,LPFILE,LPFILX,
     &        LBFLAP,LFLAP,LEIW,LSCINI,LFOREF,LNORM,LGSAME,
     &        LQSYM ,LGSYM,
     &        LQSLOP,LGSLOP, LCSLOP,
     &        LCMINP, LHMOMP, LFREQP
      LOGICAL LIQSET
      LOGICAL LAECEN
      LOGICAL TFORCE
      REAL*8 NX, NY, MASS, MINF1, MINF, MINF_CL, MVISC
      INTEGER RETYP, MATYP, AIJPIV
      CHARACTER*1 VMXBL
C
C---- dimension temporary work and storage arrays (EQUIVALENCED below)
      REAL*8 W1(6*IQX),W2(6*IQX),W3(6*IQX),W4(6*IQX),
     &       W5(6*IQX),W6(6*IQX),W7(6*IQX),W8(6*IQX)
      REAL*8 BIJ(IQX,IZX), CIJ(IWX,IQX)
C
C---- Explicit REAL*8 declarations for all REAL COMMON block variables
C     (needed so f2py generates NPY_DOUBLE instead of NPY_FLOAT)
      REAL*8 VERSION
      REAL*8 Q, DQ, DZDG, DZDN, DZDM, DQDG, DQDM
      REAL*8 QTAN1, QTAN2
      REAL*8 Z_QINF, Z_ALFA
      REAL*8 Z_QDOF0, Z_QDOF1, Z_QDOF2, Z_QDOF3
      REAL*8 AIJ, DIJ
      REAL*8 QINV, QVIS, CPI, CPV, QINVU, QINV_A
      REAL*8 X, Y, XP, YP, S
      REAL*8 SLE, XLE, YLE, XTE, YTE, CHORD, YIMAGE
      REAL*8 WGAP, WAKLEN
      REAL*8 GAM, GAMU, GAM_A, SIG, APANEL
      REAL*8 SST, SST_GO, SST_GP
      REAL*8 GAMTE, GAMTE_A, SIGTE, SIGTE_A
      REAL*8 DSTE, ANTE, ASTE
      REAL*8 SSPLE, SSPEC, XSPOC, YSPOC, QGAMM
      REAL*8 QSPEC, QSPECP
      REAL*8 ALGAM, CLGAM, CMGAM
      REAL*8 ALQSP, CLQSP, CMQSP
      REAL*8 QF0, QF1, QF2, QF3
      REAL*8 QDOF0, QDOF1, QDOF2, QDOF3, CLSPEC, FFILT
      REAL*8 ADEG, ALFA, AWAKE, AVISC
      REAL*8 XCMREF, YCMREF
      REAL*8 CL, CM, CD, CDP, CDF, CL_ALF, CL_MSQ
      REAL*8 PSIO, CIRC, COSA, SINA, QINF
      REAL*8 GAMMA, GAMM1, TKLAM, TKL_MSQ, CPSTAR, QSTAR
      REAL*8 CPMN, CPMNI, CPMNV, XCPMNI, XCPMNV
      REAL*8 PI, HOPI, QOPI, DTOR
      REAL*8 CVPAR, CTERAT, CTRRAT
      REAL*8 XSREF1, XSREF2, XPREF1, XPREF2
      REAL*8 XB, YB, XBP, YBP, SB, SNEW
      REAL*8 XBF, YBF, XOF, YOF, HMOM, HFX, HFY
      REAL*8 SBLE, CHORDB, AREAB, RADBLE, ANGBTE
      REAL*8 EI11BA, EI22BA, APX1BA, APX2BA
      REAL*8 EI11BT, EI22BT, APX1BT, APX2BT
      REAL*8 THICKB, CAMBRB
      REAL*8 XCM, YCM, SCM, XCMP, YCMP
      REAL*8 XTK, YTK, STK, XTKP, YTKP
      REAL*8 XSSI, UEDG, UINV, THET, DSTR
      REAL*8 CTAU, DELT, USLP, GUXQ, GUXD
      REAL*8 TAU, DIS, CTQ, VTI
      REAL*8 REINF1, REINF, REINF_CL, ACRIT
      REAL*8 XSTRIP, XOCTR, YOCTR, XSSITR, UINV_A
      REAL*8 RMSBL, RMXBL, RLX, VACCEL
      REAL*8 XCADD, YCADD, YCADDP
      REAL*8 XPADD, YPADD, YPADDP
      REAL*8 XCAM, YCAM, YCAMP, PCAM, PCAMP
      REAL*8 VA, VB, VDEL, VM, VZ
C
      COMMON/CR01/ VERSION
      COMMON/CC01/ NAME,ISPARS
      COMMON/QMAT/ Q(IQX,IQX),DQ(IQX),
     &             DZDG(IQX),DZDN(IQX),DZDM(IZX),
     &             DQDG(IQX),DQDM(IZX),QTAN1,QTAN2,
     &             Z_QINF,Z_ALFA,Z_QDOF0,Z_QDOF1,Z_QDOF2,Z_QDOF3
      COMMON/CR03/ AIJ(IQX,IQX),DIJ(IZX,IZX)
      COMMON/CR04/ QINV(IZX),QVIS(IZX),CPI(IZX),CPV(IZX),
     &             QINVU(IZX,2), QINV_A(IZX)
      COMMON/CR05/ X(IZX),Y(IZX),XP(IZX),YP(IZX),S(IZX),
     &             SLE,XLE,YLE,XTE,YTE,CHORD,YIMAGE,
     &             WGAP(IWX),WAKLEN
      COMMON/CR06/ GAM(IQX),GAMU(IQX,2),GAM_A(IQX),SIG(IZX),
     &             NX(IZX),NY(IZX),APANEL(IZX),
     &             SST,SST_GO,SST_GP,
     &             GAMTE,GAMTE_A,
     &             SIGTE,SIGTE_A,
     &             DSTE,ANTE,ASTE
      COMMON/CR07/ SSPLE,
     &             SSPEC(IBX),XSPOC(IBX),YSPOC(IBX),
     &             QGAMM(IBX),
     &             QSPEC(IBX,IPX),QSPECP(IBX,IPX),
     &             ALGAM,CLGAM,CMGAM,
     &             ALQSP(IPX),CLQSP(IPX),CMQSP(IPX),
     &             QF0(IQX),QF1(IQX),QF2(IQX),QF3(IQX),
     &             QDOF0,QDOF1,QDOF2,QDOF3,CLSPEC,FFILT
      COMMON/CI01/ IQ1,IQ2,NSP,NQSP,KQTARG,IACQSP,NC1,NNAME,NPREFIX
      COMMON/CR09/ ADEG,ALFA,AWAKE,MVISC,AVISC,
     &             XCMREF,YCMREF,
     &             CL,CM,CD,CDP,CDF,CL_ALF,CL_MSQ,
     &             PSIO,CIRC,COSA,SINA,QINF,
     &             GAMMA,GAMM1,
     &             MINF1,MINF,MINF_CL,TKLAM,TKL_MSQ,CPSTAR,QSTAR,
     &             CPMN,CPMNI,CPMNV,XCPMNI,XCPMNV
C
      COMMON/CR11/ PI,HOPI,QOPI,DTOR
      COMMON/CR12/ CVPAR,CTERAT,CTRRAT,XSREF1,XSREF2,XPREF1,XPREF2
      COMMON/CI04/ N,NB,NW,NPAN,IST,
     &             ITMAX,RETYP,MATYP,AIJPIV(IQX),
     &             NCM,NTK
C---- CR13 (plot sizing) removed -- all variables were plotting-only
      COMMON/CL01/ OK,LIMAGE,SHARP,
     &             LGAMU,LQINU,LVISC,LALFA,LWAKE,
     &             LBLINI,LIPAN,LQAIJ,LADIJ,LWDIJ,LCPXX,LQVDES,LQREFL,
     &             LQSPEC,LVCONV,LCPREF,LCLOCK,LPFILE,LPFILX,
     &             LBFLAP,LFLAP,LEIW,LSCINI,LFOREF,LNORM,LGSAME,
     &             LQSYM ,LGSYM,
     &             LQSLOP,LGSLOP,LCSLOP,
     &             LIQSET,
     &             LAECEN,
     &             LCMINP, LHMOMP, LFREQP
      COMMON/CR14/ XB(IBX),YB(IBX),
     &             XBP(IBX),YBP(IBX),SB(IBX),SNEW(4*IBX),
     &             XBF,YBF,XOF,YOF,HMOM,HFX,HFY,
     &             SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
     &             EI11BA,EI22BA,APX1BA,APX2BA,
     &             EI11BT,EI22BT,APX1BT,APX2BT,
     &             THICKB,CAMBRB,
     &     XCM(2*IBX),YCM(2*IBX),SCM(2*IBX),XCMP(2*IBX),YCMP(2*IBX),
     &     XTK(2*IBX),YTK(2*IBX),STK(2*IBX),XTKP(2*IBX),YTKP(2*IBX)
C
      COMMON/CR15/ XSSI(IVX,ISX),UEDG(IVX,ISX),UINV(IVX,ISX),
     &             MASS(IVX,ISX),THET(IVX,ISX),DSTR(IVX,ISX),
     &             CTAU(IVX,ISX),DELT(IVX,ISX),TSTR(IVX, ISX),
     &             USLP(IVX,ISX),GUXQ(IVX,ISX),GUXD(IVX,ISX),
     &             TAU(IVX,ISX),DIS(IVX,ISX),CTQ(IVX,ISX),
     &             VTI(IVX,ISX),
     &             REINF1,REINF,REINF_CL,
     &             ACRIT(ISX),
     &             XSTRIP(ISX),XOCTR(ISX),YOCTR(ISX),XSSITR(ISX),
     &             UINV_A(IVX,ISX),
     &             TINDEX(ISX)
      COMMON/CI05/ IBLTE(ISX),NBL(ISX),IPAN(IVX,ISX),ISYS(IVX,ISX),NSYS,
     &             ITRAN(ISX), IDAMP
      COMMON/CL02/ TFORCE(ISX)
      REAL*8 RLXPAR_DHI, RLXPAR_DLO, RLXPAR_DNUE
      COMMON/CR17/ RMSBL,RMXBL,RLX,VACCEL,
     &             RLXPAR_DHI, RLXPAR_DLO, RLXPAR_DNUE
      COMMON/CI06/ IMXBL,ISMXBL
      COMMON/CC03/ VMXBL
C
      COMMON/CR19/
     &     XCADD(NTX), YCADD(NTX), YCADDP(NTX),
     &     XPADD(NTX), YPADD(NTX), YPADDP(NTX),
     &     XCAM(NTX),
     &     YCAM(NTX), YCAMP(NTX),
     &     PCAM(NTX), PCAMP(NTX)
      COMMON/VMAT/ VA(3,2,IZX),VB(3,2,IZX),VDEL(3,2,IZX),
     &             VM(3,IZX,IZX),VZ(3,2)
C
C
C---- save storage space
      EQUIVALENCE (Q(1,1 ),W1(1)), (Q(1,7 ),W2(1)),
     &            (Q(1,13),W3(1)), (Q(1,19),W4(1)),
     &            (Q(1,25),W5(1)), (Q(1,31),W6(1)),
     &            (Q(1,37),W7(1)), (Q(1,43),W8(1))
      EQUIVALENCE (VM(1,1,1),BIJ(1,1)), (VM(1,1,IZX/2),CIJ(1,1))
C
C
C   VERSION     version number of this XFOIL implementation
C   NAME        airfoil name
C   ISPARS      ISES domain parameters  (not used in XFOIL)
C
C   Q(..)       generic coefficient matrix
C   DQ(.)       generic matrix righthand side
C
C   DZDG(.)     dPsi/dGam
C   DZDN(.)     dPsi/dn
C   DZDM(.)     dPsi/dSig
C
C   DQDG(.)     dQtan/dGam
C   DQDM(.)     dQtan/dSig
C   QTAN1       Qtan at alpha =  0 deg.
C   QTAN2       Qtan at alpha = 90 deg.
C
C   Z_QINF      dPsi/dQinf
C   Z_ALFA      dPsi/dalfa
C   Z_QDOF0     dPsi/dQdof0
C   Z_QDOF1     dPsi/dQdof1
C   Z_QDOF2     dPsi/dQdof2
C   Z_QDOF3     dPsi/dQdof3
C
C   AIJ(..)     dPsi/dGam  influence coefficient matrix (factored if LQAIJ=t)
C   BIJ(..)     dGam/dSig  influence coefficient matrix
C   CIJ(..)     dQtan/dGam influence coefficient matrix
C   DIJ(..)     dQtan/dSig influence coefficient matrix
C   QINV(.)     tangential velocity due to surface vorticity
C   QVIS(.)     tangential velocity due to surface vorticity & mass sources
C   QINVU(..)   QINV for alpha = 0, 90 deg.
C   QINV_A(.)   dQINV/dalpha
C
C   X(.),Y(.)   airfoil (1<i<N) and wake (N+1<i<N+NW) coordinate arrays
C   XP(.),YP(.) dX/dS, dY/dS arrays for spline evaluation
C   S(.)        arc length along airfoil (spline parameter)
C   SLE         value of S at leading edge
C   XLE,YLE     leading  edge coordinates
C   XTE,YTE     trailing edge coordinates
C   WGAP(.)     thickness of "dead air" region inside wake just behind TE
C   WAKLEN      wake length to chord ratio
C
C   GAM(.)      surface vortex panel strength array
C   GAMU(.2)    surface vortex panel strength arrays for alpha = 0, 90 deg.
C   GAM_A(.)    dGAM/dALFA
C   SIG(.)      surface and wake mass defect array
C
C   NX(.),NY(.) normal unit vector components at airfoil and wake coordinates
C   APANEL(.)   surface and wake panel angle array (+ counterclockwise)
C
C   SST         S value at stagnation point
C   SST_GO      dSST/dGAM(IST)
C   SST_GP      dSST/dGAM(IST+1)
C
C   GAMTE       vortex panel strength across finite-thickness TE
C   SIGTE       source panel strength across finite-thickness TE
C   GAMTE_A     dGAMTE/dALFA
C   SIGTE_A     dSIGTE/dALFA
C   DSTE        TE panel length
C   ANTE,ASTE   projected TE thickness perp.,para. to TE bisector
C   SHARP       .TRUE.  if  DSTE.EQ.0.0 ,  .FALSE. otherwise
C
C   SSPEC(.)    normalized arc length around airfoil (QSPEC coordinate)
C   XSPOC(.)    x/c at SSPEC points
C   YSPOC(.)    y/c at SSPEC points
C   QSPEC(..)   specified surface velocity for inverse calculations
C   QSPECP(..)  dQSPEC/dSSPEC
C   QGAMM(.)    surface velocity for current airfoil geometry
C   SSPLE       SSPEC value at airfoil nose
C
C   IQ1,IQ2     target segment endpoint indices on Qspec(s) plot
C   NSP         number of points in QSPEC array
C   NQSP        number Qspec arrays
C   IACQSP      1:  ALQSP is prescribed for Qspec arrays
C               2:  CLQSP is prescribed for Qspec arrays
C   NC1         number of circle plane points, must be 2**n - 1
C
C   NNAME       number of characters in airfoil name
C   NPREFIX     number of characters in default filename prefix
C
C   ALQSP(.)    alpha,CL,CM corresponding to QSPEC distributions
C   CLQSP(.)    
C   CMQSP(.)    
C   ALGAM       alpha,CL,CM corresponding to QGAMM distribution
C   CLGAM
C   CMGAM
C
C   QF0(.)      shape function for QSPEC modification
C   QF1(.)        "
C   QF2(.)        "
C   QF3(.)        "
C   QDOF0       shape function weighting coefficient (inverse DOF)
C   QDOF1         "
C   QDOF2         "
C   QDOF3         "
C   CLSPEC      specified CL
C   FFILT       circle-plane mapping filter parameter
C
C   ADEG,ALFA   angle of attack in degrees, radians
C   AWAKE       angle of attack corresponding to wake geometry (radians)
C   AVISC       angle of attack corresponding to BL solution   (radians)
C   MVISC       Mach number corresponding to BL solution
C   CL,CM       current CL and CM calculated from GAM(.) distribution
C   CD          current CD from BL solution
C   CDF         current friction CD from BL solution
C   CL_ALF      dCL/dALFA
C   CL_MSQ      dCL/d(MINF^2)
C
C   PSIO        streamfunction inside airfoil
C   CIRC        circulation
C   COSA,SINA   cos(ALFA), sin(ALFA)
C   QINF        freestream speed    (defined as 1)
C   GAMMA,GAMM1 Gas constant Cp/Cv, Cp/Cv - 1
C   MINF1       freestream Mach number at CL=1
C   MINF        freestream Mach number at current CL
C   MINF_CL     dMINF/dCL
C   TKLAM       Karman-Tsien parameter Minf^2 / [1 + sqrt(1-Minf^2)]^2
C   TKL_MSQ     d(TKLAM)/d(MINF^2)
C   CPSTAR      sonic pressure coefficient
C   QSTAR       sonic speed
C
C   PI          3.1415926...
C   HOPI,QOPI   1/(2 Pi) ,  1/(4 Pi)
C   DTOR        Pi / 180    (degrees to radians conversion factor)
C
C   CVPAR       curvature attraction parameter for airfoil paneling
C               0 = uniform panel node spacing around airfoil
C              ~1 = panel nodes strongly bunched in areas of large curvature
C   CTERAT      TE panel density / LE panel density ratio
C   CTRRAT      local refinement panel density / LE panel density ratio
C   XSREF1-2    suction  side local refinement x/c limits
C   XPREF1-2    pressure side local refinement x/c limits
C
C   N           number of points on airfoil
C   NB          number of points in buffer airfoil array
C   NW          number of points in wake
C   NPAN        default/specified number of points on airfoil
C
C   IST         stagnation point lies between S(IST), S(IST+1)
C   ITMAX       max number of Newton iterations
C
C   RETYP       index giving type of Re variation with CL ...
C            ... 1  Re constant 
C            ... 2  Re ~ 1/sqrt(CL)    (fixed lift)
C            ... 3  Re ~ 1/CL          (fixed lift and dynamic pressure)
C
C   MATYP       index giving type of Ma variation with CL ...
C            ... 1  Ma constant 
C            ... 2  Ma ~ 1/sqrt(CL)    (fixed lift)
C
C   AIJPIV(.)   pivot index array for LU factoring routine
C
C
C   OK          user question response
C   LIMAGE      .TRUE. if image airfoil is present
C   LGAMU       .TRUE. if GAMU  arrays exist for current airfoil geometry
C   LQINU       .TRUE. if QINVU arrays exist for current airfoil geometry
C   LVISC       .TRUE. if viscous option is invoked
C   LALFA       .TRUE. if alpha is specifed, .FALSE. if CL is specified
C   LWAKE       .TRUE. if wake geometry has been calculated
C   LBLINI      .TRUE. if BL has been initialized
C   LIPAN       .TRUE. if BL->panel pointers IPAN have been calculated
C   LQAIJ       .TRUE. if dPsi/dGam matrix has been computed and factored
C   LADIJ       .TRUE. if dQ/dSig matrix for the airfoil has been computed
C   LWDIJ       .TRUE. if dQ/dSig matrix for the wake has been computed
C   LQVDES      .TRUE. if viscous Ue is to be plotted in QDES routines
C   LQSPEC      .TRUE. if Qspec has been initialized
C   LQREFL      .TRUE. if reflected Qspec is to be plotted in QDES routines
C   LVCONV      .TRUE. if converged BL solution exists
C   LCPREF      .TRUE. if reference data is to be plotted on Cp vs x/c plots
C   LCLOCK      .TRUE. if source airfoil coordinates are clockwise
C   LPFILE      .TRUE. if polar file is ready to be appended to
C   LPFILX      .TRUE. if polar dump file is ready to be appended to
C   LBFLAP      .TRUE. if buffer  airfoil flap parameters are defined
C   LFLAP       .TRUE. if current airfoil flap parameters are defined
C   LEIW        .TRUE. if unit circle complex number array is initialized
C   LSCINI      .TRUE. if old-airfoil circle-plane arc length s(w) exists
C   LFOREF      .TRUE. if CL,CD... data is to be plotted on Cp vs x/c plots
C   LNORM       .TRUE. if input buffer airfoil is to be normalized
C   LGSAME      .TRUE. if current and buffer airfoils are identical
C
C   LQSYM       .TRUE. if symmetric Qspec will be enforced
C   LGSYM       .TRUE. if symmetric geometry will be enforced
C   LQSLOP      .TRUE. if modified Qspec(s) segment is to match slopes
C   LGSLOP      .TRUE. if modified geometry segment is to match slopes
C   LCSLOP      .TRUE. if modified camber line segment is to match slopes
C   LCMINP      .TRUE. if min Cp is to be written to polar file for cavitation
C   LHMOMP      .TRUE. if hinge moment is to be written to polar file
C   LFREQP      .TRUE. if individual TS-wave frequencies are to be plotted
C
C   LIQSET      .TRUE. if inverse target segment is marked off in QDES
C
C
C   XB(.),YB(.) buffer airfoil coordinate arrays
C   XBP(.)      dXB/dSB
C   YBP(.)      dYB/dSB
C   SB(.)       spline parameter for buffer airfoil
C   SNEW(.)     new panel endpoint arc length array
C
C   XBF,YBF     buffer  airfoil flap hinge coordinates
C   XOF,YOF     current airfoil flap hinge coordinates
C   HMOM        moment of flap about hinge point
C   HFX         x-force of flap on hinge point
C   HFY         y-force of flap on hinge point
C
C~~~~~~~~~~~~~~ properties of current buffer airfoil
C
C   SBLE         LE tangency-point SB location
C   CHORDB       chord
C   AREAB        area
C   RADBLE       LE radius
C   ANGBTE       TE angle  (rad)
C   
C   EI11BA       bending inertia about axis 1    x^2 dx dy
C   EI22BA       bending inertia about axis 2    y^2 dx dy
C   APX1BA       principal axis 1 angle 
C   APX2BA       principal axis 2 angle 
C
C   EI11BT       bending inertia about axis 1    x^2 t ds
C   EI22BT       bending inertia about axis 2    y^2 t ds
C   APX1BT       principal axis 1 angle 
C   APX2BT       principal axis 2 angle 
C
C   THICKB       max thickness
C   CAMBRB       max camber
C
C~~~~~~~~~~~~~~
C
C   XSSI(..)    BL arc length coordinate array on each surface
C   UEDG(..)    BL edge velocity array
C   UINV(..)    BL edge velocity array without mass defect influence
C   MASS(..)    BL mass defect array  ( = UEDG*DSTR )
C   THET(..)    BL momentum thickness array
C   DSTR(..)    BL displacement thickness array
C   TSTR(..)    BL kin. energy  thickness array
C   CTAU(..)    sqrt(max shear coefficient) array
C               (in laminar regions, log of amplification ratio)
C   TSTR(..)    BL kin. energy  thickness array
C
C   TAU(..)     wall shear stress array                 (for plotting only)
C   DIS(..)     dissipation array                       (for plotting only)
C   CTQ(..)     sqrt(equilibrium max shear coefficient) array (  "  )
C   VTI(..)     +/-1 conversion factor between panel and BL variables
C   UINV_A(..)  dUINV/dalfa array
C 
C   REINF1      Reynolds number  Vinf c / ve  for CL=1
C   REINF       Reynolds number for current CL
C   REINF_CL    dREINF/dCL
C
C   ACRIT       log (critical amplification ratio)
C   XSTRIP(.)   transition trip  x/c locations (if XTRIP > 0),
C               transition trip -s/s_side locations (if XTRIP < 0),
C   XOCTR(.)    actual transition x/c locations
C   YOCTR(.)    actual transition y/c locations
C   XSSITR(.)   actual transition xi locations
C
C   IBLTE(.)    BL array index at trailing edge
C   NBL(.)      max BL array index
C   IPAN(..)    panel index corresponding to BL location
C   ISYS(..)    BL Newton system line number corresponding to BL location
C   NSYS        total number of lines in BL Newton system
C   ITRAN(.)    BL array index of transition interval
C   TFORCE(.)   .TRUE. if transition is forced due to transition strip
C   TINDEX(.)   
C
C   IDAMP    = 0   use original enelope e^n f(H,Rtheta) for all profiles
C            = 1   use modified enelope e^n f(H,Rtheta) for separating profile
C
C   VA,VB(...)  diagonal and off-diagonal blocks in BL Newton system
C   VZ(..)      way-off-diagonal block at TE station line
C   VM(...)     mass-influence coefficient vectors in BL Newton system
C   VDEL(..)    residual and solution vectors in BL Newton system
C
C   RMSBL       rms change from BL Newton system solution
C   RMXBL       max change from BL Newton system solution
C   IMXBL       location of max change
C   ISMXBL      index of BL side containing max change
C   VMXBL       character identifying variable with max change
C   RLX         underrelaxation factor for Newton update
C   VACCEL      parameter for accelerating BL Newton system solution
C               (any off-diagonal element < VACCEL is not eliminated,
C                which speeds up each iteration, but MAY increase
C                iteration count)
C                Can be set to zero for unadulterated Newton method
C
