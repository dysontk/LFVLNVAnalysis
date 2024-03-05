      DOUBLE PRECISION FUNCTION DSIG8(PP,WGT,IMODE)
C     ****************************************************
C     
C     Generated by MadGraph5_aMC@NLO v. 3.5.3, 2023-12-23
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Process: g d~ > w+ u u~ u~ WEIGHTED<=5 @4
C     *   Decay: w+ > e+ ve WEIGHTED<=2
C     Process: g s~ > w+ c c~ c~ WEIGHTED<=5 @4
C     *   Decay: w+ > e+ ve WEIGHTED<=2
C     
C     RETURNS DIFFERENTIAL CROSS SECTION
C     Input:
C     pp    4 momentum of external particles
C     wgt   weight from Monte Carlo
C     imode 0 run, 1 init, 2 reweight, 
C     3 finalize, 4 only PDFs,
C     5 squared amplitude only (never
C     generate events)
C     Output:
C     Amplitude squared and summed
C     ****************************************************
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INCLUDE 'genps.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'maxconfigs.inc'
      INCLUDE 'maxamps.inc'
      DOUBLE PRECISION       CONV
      PARAMETER (CONV=389379.66*1000)  !CONV TO PICOBARNS
      REAL*8     PI
      PARAMETER (PI=3.1415926D0)
C     
C     ARGUMENTS 
C     
      DOUBLE PRECISION PP(0:3,NEXTERNAL), WGT
      INTEGER IMODE
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,ITYPE,LP,IPROC
      DOUBLE PRECISION G1
      DOUBLE PRECISION SX2,DX2
      DOUBLE PRECISION XPQ(-7:7),PD(0:MAXPROC)
      DOUBLE PRECISION DSIGUU,R,RCONF
      INTEGER LUN,ICONF,IFACT,NFACT
      DATA NFACT/1/
      SAVE NFACT
C     
C     STUFF FOR DRESSED EE COLLISIONS
C     
      INCLUDE '../../Source/PDF/eepdf.inc'
      DOUBLE PRECISION EE_COMP_PROD

      INTEGER I_EE
C     
C     STUFF FOR UPC
C     
      DOUBLE PRECISION PHOTONPDFSQUARE
C     
C     EXTERNAL FUNCTIONS
C     
      LOGICAL PASSCUTS
      DOUBLE PRECISION ALPHAS2,REWGT,PDG2PDF,CUSTOM_BIAS
      INTEGER NEXTUNOPEN
C     
C     GLOBAL VARIABLES
C     
      INTEGER          IPSEL
      COMMON /SUBPROC/ IPSEL
C     MINCFIG has this config number
      INTEGER           MINCFIG, MAXCFIG
      COMMON/TO_CONFIGS/MINCFIG, MAXCFIG
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      COMMON/TO_MCONFIGS/MAPCONFIG, ICONFIG
C     Keep track of whether cuts already calculated for this event
      LOGICAL CUTSDONE,CUTSPASSED
      COMMON/TO_CUTSDONE/CUTSDONE,CUTSPASSED

      INTEGER SUBDIAG(MAXSPROC),IB(2)
      COMMON/TO_SUB_DIAG/SUBDIAG,IB
      INCLUDE 'coupl.inc'
      INCLUDE 'run.inc'
      INCLUDE '../../Source/PDF/pdf.inc'
C     Common blocks
C     CHARACTER*7         PDLABEL,EPA_LABEL
C     INTEGER       LHAID
C     COMMON/TO_PDF/LHAID,PDLABEL,EPA_LABEL
C     
C     local
C     
      DOUBLE PRECISION P1(0:3, NEXTERNAL)

C     
C     DATA
C     
      DATA G1/1*1D0/
      DATA SX2,DX2/2*1D0/
C     ----------
C     BEGIN CODE
C     ----------
      DSIG8=0D0

      IF(IMODE.EQ.1)THEN
C       Set up process information from file symfact
        LUN=NEXTUNOPEN()
        NFACT=1
        OPEN(UNIT=LUN,FILE='../symfact.dat',STATUS='OLD',ERR=20)
        DO WHILE(.TRUE.)
          READ(LUN,*,ERR=10,END=10) RCONF, IFACT
          ICONF=INT(RCONF)
          IF(ICONF.EQ.MAPCONFIG(MINCFIG))THEN
            NFACT=IFACT
          ENDIF
        ENDDO
 10     CLOSE(LUN)
        RETURN
 20     WRITE(*,*)'Error opening symfact.dat. No symmetry factor used.'
        RETURN
      ENDIF
C     Continue only if IMODE is 0, 4 or 5
      IF(IMODE.NE.0.AND.IMODE.NE.4.AND.IMODE.NE.5) RETURN


      IF (ABS(LPP(IB(1))).GE.1) THEN
          !LP=SIGN(1,LPP(IB(1)))
        G1=PDG2PDF(LPP(IB(1)),0, IB(1),XBK(IB(1)),DSQRT(Q2FACT(IB(1))))
      ENDIF
      IF (ABS(LPP(IB(2))).GE.1) THEN
          !LP=SIGN(1,LPP(IB(2)))
        SX2=PDG2PDF(LPP(IB(2)),-3, IB(2),XBK(IB(2)),DSQRT(Q2FACT(IB(2))
     $   ))
        DX2=PDG2PDF(LPP(IB(2)),-1, IB(2),XBK(IB(2)),DSQRT(Q2FACT(IB(2))
     $   ))
      ENDIF
      PD(0) = 0D0
      IPROC = 0
      IPROC=IPROC+1  ! g d~ > e+ ve u u~ u~
      PD(IPROC)=G1*DX2
      PD(0)=PD(0)+DABS(PD(IPROC))
      IPROC=IPROC+1  ! g s~ > e+ ve c c~ c~
      PD(IPROC)=G1*SX2
      PD(0)=PD(0)+DABS(PD(IPROC))
      IF (IMODE.EQ.4)THEN
        DSIG8 = PD(0)
        RETURN
      ENDIF
      IF(FRAME_ID.NE.6)THEN
        CALL BOOST_TO_FRAME(PP, FRAME_ID, P1)
      ELSE
        P1 = PP
      ENDIF
      CALL SMATRIX8(P1,DSIGUU)
      IF (IMODE.EQ.5) THEN
        IF (DSIGUU.LT.1D199) THEN
          DSIG8 = DSIGUU*CONV
        ELSE
          DSIG8 = 0.0D0
        ENDIF
        RETURN
      ENDIF
C     Select a flavor combination (need to do here for right sign)
      CALL RANMAR(R)
      IPSEL=0
      DO WHILE (R.GE.0D0 .AND. IPSEL.LT.IPROC)
        IPSEL=IPSEL+1
        R=R-DABS(PD(IPSEL))/PD(0)
      ENDDO

      DSIGUU=DSIGUU*REWGT(PP)

C     Apply the bias weight specified in the run card (default is 1.0)
      DSIGUU=DSIGUU*CUSTOM_BIAS(PP,DSIGUU,8)

      DSIGUU=DSIGUU*NFACT

      IF (DSIGUU.LT.1D199) THEN
C       Set sign of dsig based on sign of PDF and matrix element
        DSIG8=DSIGN(PD(0)*CONV*DSIGUU,DSIGUU*PD(IPSEL))
      ELSE
        WRITE(*,*) 'Error in matrix element'
        DSIGUU=0D0
        DSIG8=0D0
      ENDIF
C     Generate events only if IMODE is 0.
      IF(IMODE.EQ.0.AND.DABS(DSIG8).GT.0D0)THEN
C       Call UNWGT to unweight and store events
        CALL UNWGT(PP,DSIG8*WGT,8)
      ENDIF

      END
C     
C     Functionality to handling grid
C     




      SUBROUTINE PRINT_ZERO_AMP8()

      RETURN
      END

      INTEGER FUNCTION GET_NHEL8(HEL, IPART)
C     if hel>0 return the helicity of particule ipart for the selected
C      helicity configuration
C     if hel=0 return the number of helicity state possible for that
C      particle 
      IMPLICIT NONE
      INTEGER HEL,I, IPART
      INCLUDE 'nexternal.inc'
      INTEGER ONE_NHEL(NEXTERNAL)
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=128)
      INTEGER NHEL(NEXTERNAL,0:NCOMB)
      DATA (NHEL(I,0),I=1,6) / 2, 2, 3, 2, 2, 2/
      DATA (NHEL(I,   1),I=1,7) /-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,   2),I=1,7) /-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,   3),I=1,7) /-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,   4),I=1,7) /-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,   5),I=1,7) /-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,   6),I=1,7) /-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,   7),I=1,7) /-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,   8),I=1,7) /-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,   9),I=1,7) /-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  10),I=1,7) /-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  11),I=1,7) /-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  12),I=1,7) /-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  13),I=1,7) /-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  14),I=1,7) /-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  15),I=1,7) /-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  16),I=1,7) /-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  17),I=1,7) /-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  18),I=1,7) /-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  19),I=1,7) /-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  20),I=1,7) /-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  21),I=1,7) /-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  22),I=1,7) /-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  23),I=1,7) /-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  24),I=1,7) /-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  25),I=1,7) /-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  26),I=1,7) /-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  27),I=1,7) /-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  28),I=1,7) /-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  29),I=1,7) /-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  30),I=1,7) /-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  31),I=1,7) /-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  32),I=1,7) /-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  33),I=1,7) /-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  34),I=1,7) /-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  35),I=1,7) /-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  36),I=1,7) /-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  37),I=1,7) /-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  38),I=1,7) /-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  39),I=1,7) /-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  40),I=1,7) /-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  41),I=1,7) /-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  42),I=1,7) /-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  43),I=1,7) /-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  44),I=1,7) /-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  45),I=1,7) /-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  46),I=1,7) /-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  47),I=1,7) /-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  48),I=1,7) /-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  49),I=1,7) /-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  50),I=1,7) /-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  51),I=1,7) /-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  52),I=1,7) /-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  53),I=1,7) /-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  54),I=1,7) /-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  55),I=1,7) /-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  56),I=1,7) /-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  57),I=1,7) /-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  58),I=1,7) /-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  59),I=1,7) /-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  60),I=1,7) /-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  61),I=1,7) /-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  62),I=1,7) /-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  63),I=1,7) /-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  64),I=1,7) /-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  65),I=1,7) / 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  66),I=1,7) / 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  67),I=1,7) / 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  68),I=1,7) / 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  69),I=1,7) / 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  70),I=1,7) / 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  71),I=1,7) / 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  72),I=1,7) / 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  73),I=1,7) / 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  74),I=1,7) / 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  75),I=1,7) / 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  76),I=1,7) / 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  77),I=1,7) / 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  78),I=1,7) / 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  79),I=1,7) / 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  80),I=1,7) / 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  81),I=1,7) / 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  82),I=1,7) / 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  83),I=1,7) / 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  84),I=1,7) / 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  85),I=1,7) / 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  86),I=1,7) / 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  87),I=1,7) / 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  88),I=1,7) / 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  89),I=1,7) / 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  90),I=1,7) / 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  91),I=1,7) / 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  92),I=1,7) / 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  93),I=1,7) / 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  94),I=1,7) / 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  95),I=1,7) / 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  96),I=1,7) / 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  97),I=1,7) / 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  98),I=1,7) / 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  99),I=1,7) / 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I, 100),I=1,7) / 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I, 101),I=1,7) / 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I, 102),I=1,7) / 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I, 103),I=1,7) / 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I, 104),I=1,7) / 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I, 105),I=1,7) / 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I, 106),I=1,7) / 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I, 107),I=1,7) / 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I, 108),I=1,7) / 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I, 109),I=1,7) / 1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(I, 110),I=1,7) / 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I, 111),I=1,7) / 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I, 112),I=1,7) / 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I, 113),I=1,7) / 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I, 114),I=1,7) / 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I, 115),I=1,7) / 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I, 116),I=1,7) / 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I, 117),I=1,7) / 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I, 118),I=1,7) / 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I, 119),I=1,7) / 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I, 120),I=1,7) / 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I, 121),I=1,7) / 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I, 122),I=1,7) / 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I, 123),I=1,7) / 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I, 124),I=1,7) / 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I, 125),I=1,7) / 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I, 126),I=1,7) / 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I, 127),I=1,7) / 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I, 128),I=1,7) / 1, 1,-1, 1, 1,-1,-1/

      GET_NHEL8 = NHEL(IPART, IABS(HEL))
      RETURN
      END
