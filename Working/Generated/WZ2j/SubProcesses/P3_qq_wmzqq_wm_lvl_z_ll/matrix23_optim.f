      SUBROUTINE SMATRIX23(P,ANS)
C     
C     Generated by MadGraph5_aMC@NLO v. 3.5.3, 2023-12-23
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     MadGraph5_aMC@NLO for Madevent Version
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: d u~ > w- z s s~ WEIGHTED<=6 @3
C     *   Decay: w- > e- ve~ WEIGHTED<=2
C     *   Decay: z > e+ e- WEIGHTED<=2
C     Process: s c~ > w- z d d~ WEIGHTED<=6 @3
C     *   Decay: w- > e- ve~ WEIGHTED<=2
C     *   Decay: z > e+ e- WEIGHTED<=2
C     Process: d u~ > w- z s s~ WEIGHTED<=6 @3
C     *   Decay: w- > mu- vm~ WEIGHTED<=2
C     *   Decay: z > mu+ mu- WEIGHTED<=2
C     Process: s c~ > w- z d d~ WEIGHTED<=6 @3
C     *   Decay: w- > mu- vm~ WEIGHTED<=2
C     *   Decay: z > mu+ mu- WEIGHTED<=2
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      Include 'genps.inc'
      Include 'maxconfigs.inc'
      Include 'nexternal.inc'
      Include 'maxamps.inc'
      INTEGER                 NCOMB         
      PARAMETER (             NCOMB=4)
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=12) 
      INTEGER    NDIAGS
      PARAMETER (NDIAGS=12) 
      INTEGER    THEL
      PARAMETER (THEL=2*NCOMB)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS
C     
C     LOCAL VARIABLES 
C     
      INTEGER CONFSUB(MAXSPROC,LMAXCONFIGS)
      INCLUDE 'config_subproc_map.inc'
      INTEGER NHEL(NEXTERNAL,NCOMB)
      REAL*8 T
      REAL*8 R,SUMHEL,TS(NCOMB)
      INTEGER I,IDEN
      INTEGER JC(NEXTERNAL),II
      REAL*8 XTOT
      INTEGER  J, JJ

      double precision get_channel_cut
      external get_channel_cut

C     
C     GLOBAL VARIABLES
C     
      DOUBLE PRECISION AMP2(MAXAMPS), JAMP2(0:MAXFLOW)
      COMMON/TO_AMPS/  AMP2,       JAMP2


C     
C     INFORMATION TO WRITE THE HELICITY IN THE EVENT --not memory
C      efficient--
C     
      CHARACTER*101         HEL_BUFF
      COMMON/TO_HELICITY/  HEL_BUFF

      INTEGER NB_SPIN_STATE_in(2)
      common /nb_hel_state/ nb_spin_state_in

      REAL*8 POL(2)

      COMMON/TO_POLARIZATION/ POL
      double precision tmin_for_channel
      integer sde_strat    !  1 means standard single diagram enhancement strategy,
C     2 means approximation by the	denominator of the propagator
      common/TO_CHANNEL_STRAT/tmin_for_channel,	sde_strat

      INTEGER          ISUM_HEL
      LOGICAL                    MULTI_CHANNEL
      COMMON/TO_MATRIX/ISUM_HEL, MULTI_CHANNEL
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      COMMON/TO_MCONFIGS/MAPCONFIG, ICONFIG
      INTEGER SUBDIAG(MAXSPROC),IB(2)
      COMMON/TO_SUB_DIAG/SUBDIAG,IB
      
      DATA (NHEL(I,1),I=1,8) /-1, 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,2),I=1,8) /-1, 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,3),I=1,8) /-1, 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,4),I=1,8) /-1, 1,-1, 1,-1, 1, 1,-1/
      DATA IDEN/36/

C     ----------
C     BEGIN CODE
C     ----------

      DO I=1,NEXTERNAL
        JC(I) = +1
      ENDDO

      IF (multi_channel) THEN
        DO I=1,NDIAGS
          AMP2(I)=0D0
        ENDDO
        JAMP2(0)=2
        DO I=1,INT(JAMP2(0))
          JAMP2(I)=0D0
        ENDDO
      ENDIF
      ANS = 0D0
      WRITE(HEL_BUFF,'(20I5)') (0,I=1,NEXTERNAL)
C     Kiran please check if you need this:    
      DO I=1,NCOMB
        TS(I)=0d0
      ENDDO

      call MATRIX23(P ,JC(1), TS)
      DO I=1,NCOMB     
        T=TS(I)  
        DO JJ=1,nincoming
          IF(POL(JJ).NE.1d0.AND.NHEL(JJ,I).EQ.INT(SIGN(1d0,POL(JJ))))
     $      THEN
            T=T*ABS(POL(JJ))*NB_SPIN_STATE_IN(JJ)/2d0    !  NB_SPIN_STATE(JJ)/2d0 is added for polarised beam
          ELSE IF(POL(JJ).NE.1d0)THEN
            T=T*(2d0-ABS(POL(JJ)))*NB_SPIN_STATE_IN(JJ)/2d0
          ENDIF
        ENDDO
        ANS=ANS+DABS(T)
        TS(I)=T
      ENDDO

      IF (ANS.ne.0d0) THEN
        CALL RANMAR(R)
        SUMHEL=0d0
        DO I=1,NCOMB
          SUMHEL=SUMHEL+DABS(TS(I))/ANS
          IF(R.LT.SUMHEL)THEN
            WRITE(HEL_BUFF,'(20i5)')(NHEL(II,I),II=1,NEXTERNAL)
C           Set right sign for ANS, based on sign of chosen helicity
            ANS=DSIGN(ANS,TS(I))
            GOTO 10
          ENDIF
        ENDDO
 10     CONTINUE   
      ENDIF
      IF (MULTI_CHANNEL) THEN
        XTOT=0D0
        DO I=1,LMAXCONFIGS
          J = CONFSUB(23, I)
          if (J.ne.0)then
            if (sde_strat.eq.1)then
              AMP2(J) = AMP2(J) * GET_CHANNEL_CUT(P, I)
            else
              AMP2(J) = GET_CHANNEL_CUT(P, I)
            endif
            XTOT=XTOT+AMP2(J)

          endif
        ENDDO
        IF (XTOT.NE.0D0) THEN
          ANS=ANS*AMP2(SUBDIAG(23))/XTOT
        ELSE IF(ANS.ne.0d0) THEN
          write(*,*) 'Problem in the multi-channeling. All amp2 are'
     $     //' zero but not the total matrix-element'
          stop 1
        ENDIF
      ENDIF
      ANS=ANS/DBLE(IDEN)
      END


      Subroutine  MATRIX23(P,IC, TS)
C     
C     Generated by MadGraph5_aMC@NLO v. 3.5.3, 2023-12-23
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: d u~ > w- z s s~ WEIGHTED<=6 @3
C     *   Decay: w- > e- ve~ WEIGHTED<=2
C     *   Decay: z > e+ e- WEIGHTED<=2
C     Process: s c~ > w- z d d~ WEIGHTED<=6 @3
C     *   Decay: w- > e- ve~ WEIGHTED<=2
C     *   Decay: z > e+ e- WEIGHTED<=2
C     Process: d u~ > w- z s s~ WEIGHTED<=6 @3
C     *   Decay: w- > mu- vm~ WEIGHTED<=2
C     *   Decay: z > mu+ mu- WEIGHTED<=2
C     Process: s c~ > w- z d d~ WEIGHTED<=6 @3
C     *   Decay: w- > mu- vm~ WEIGHTED<=2
C     *   Decay: z > mu+ mu- WEIGHTED<=2
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=12) 
      include 'genps.inc'
      include 'nexternal.inc'
      include 'maxamps.inc'
      INTEGER    NWAVEFUNCS,     NCOLOR
      PARAMETER (NWAVEFUNCS=43, NCOLOR=2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0D0,1D0))
      INTEGER NAMPSO, NSQAMPSO
      PARAMETER (NAMPSO=1, NSQAMPSO=1)
      LOGICAL CHOSEN_SO_CONFIGS(NSQAMPSO)
      DATA CHOSEN_SO_CONFIGS/.TRUE./
      SAVE CHOSEN_SO_CONFIGS
      INTEGER                 NCOMB         
      PARAMETER (             NCOMB=4)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      REAL*8 TS(NCOMB)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,M,N,K
      COMPLEX*16 ZTEMP,TMP_JAMP(11)
      COMPLEX*16 TMP(6)
      REAL*8 CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NCOMB,NGRAPHS), JAMP(NCOLOR,NAMPSO)
      COMPLEX*16 W(6,NWAVEFUNCS)
C     Needed for v4 models
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0d0, 0d0), (1d0, 0d0)/

      double precision fk_ZERO 
      double precision fk_mdl_WW 
      double precision fk_mdl_WZ 
      save fk_ZERO 
      save fk_mdl_WW 
      save fk_mdl_WZ 

      logical first
      data first /.true./
      save first
C     
C     FUNCTION
C     
      INTEGER SQSOINDEX23

C     
C     GLOBAL VARIABLES
C     
      Double Precision amp2(maxamps), jamp2(0:maxflow)
      common/to_amps/  amp2,       jamp2
      include 'coupl.inc'

      double precision tmin_for_channel
      integer sde_strat    !  1 means standard single diagram enhancement strategy,
C     2 means approximation by the	denominator of the propagator
      common/TO_CHANNEL_STRAT/tmin_for_channel,	sde_strat

      double precision small_width_treatment
      common/narrow_width/small_width_treatment
C     
C     COLOR DATA
C     
      DATA (CF(i,  1),i=  1,  2) /9.000000000000000d+00
     $ ,3.000000000000000d+00/
C     1 T(2,1) T(7,8)
      DATA (CF(i,  2),i=  1,  2) /3.000000000000000d+00
     $ ,9.000000000000000d+00/
C     1 T(2,8) T(7,1)
C     ----------
C     BEGIN CODE
C     ----------
      if (first) then
        first=.false.
        IF(ZERO.ne.0d0) fk_ZERO = SIGN(MAX(ABS(ZERO), ABS(ZERO
     $   *small_width_treatment)), ZERO)
        IF(mdl_WW.ne.0d0) fk_mdl_WW = SIGN(MAX(ABS(mdl_WW), ABS(mdl_MW
     $   *small_width_treatment)), mdl_WW)
        IF(mdl_WZ.ne.0d0) fk_mdl_WZ = SIGN(MAX(ABS(mdl_WZ), ABS(mdl_MZ
     $   *small_width_treatment)), mdl_WZ)
      endif
      AMP(:,:) = (0d0,0d0)
            CALL IXXXXX(P(0,1),ZERO,-1,+1*IC(1),W(1,1))  !  count 8
      CALL OXXXXX(P(0,2),ZERO,+1,-1*IC(2),W(1,2))  !  count 20
      CALL OXXXXX(P(0,3),ZERO,-1,+1*IC(3),W(1,3))  !  count 1
      CALL IXXXXX(P(0,4),ZERO,+1,-1*IC(4),W(1,4))  !  count 1
      CALL FFV2_3(W(1,4),W(1,3),GC_100,MDL_MW,FK_MDL_WW,W(1,5))  !  count 10
      CALL IXXXXX(P(0,5),ZERO,+1,-1*IC(5),W(1,6))  !  count 1
      CALL IXXXXX(P(0,5),ZERO,-1,-1*IC(5),W(1,7))  !  count 1
      CALL OXXXXX(P(0,6),ZERO,+1,+1*IC(6),W(1,8))  !  count 1
      CALL OXXXXX(P(0,6),ZERO,-1,+1*IC(6),W(1,9))  !  count 1
      CALL FFV2_4_3(W(1,6),W(1,9),GC_50,GC_59,MDL_MZ,FK_MDL_WZ,W(1,10))  !  count 10
      CALL FFV2_4_3(W(1,7),W(1,8),GC_50,GC_59,MDL_MZ,FK_MDL_WZ,W(1,11))  !  count 10
      CALL OXXXXX(P(0,7),ZERO,+1,+1*IC(7),W(1,12))  !  count 7
      CALL OXXXXX(P(0,7),ZERO,-1,+1*IC(7),W(1,13))  !  count 7
      CALL IXXXXX(P(0,8),ZERO,+1,-1*IC(8),W(1,14))  !  count 7
      CALL IXXXXX(P(0,8),ZERO,-1,-1*IC(8),W(1,15))  !  count 7
      CALL FFV2_2(W(1,1),W(1,5),GC_100,ZERO,FK_ZERO,W(1,16))  !  count 7
      CALL FFV2_3_2(W(1,14),W(1,10),GC_50,GC_58,ZERO,FK_ZERO,W(1,17))  !  count 2
      CALL FFV2_3_2(W(1,14),W(1,11),GC_50,GC_58,ZERO,FK_ZERO,W(1,18))  !  count 2
      CALL FFV2_3_2(W(1,15),W(1,10),GC_50,GC_58,ZERO,FK_ZERO,W(1,19))  !  count 2
      CALL FFV2_3_2(W(1,15),W(1,11),GC_50,GC_58,ZERO,FK_ZERO,W(1,20))  !  count 2
      CALL FFV1P0_3(W(1,16),W(1,2),GC_11,ZERO,FK_ZERO,W(1,21))  !  count 8
      CALL FFV1P1N_1(W(1,13), W(1,21), GC_11, TMP(1))
      call CombineAmp(2,
     & (/1,3/), 
     & (/17,18/),
     & TMP, W, AMP(1,1))
      CALL FFV1P1N_1(W(1,12), W(1,21), GC_11, TMP(1))
      call CombineAmp(2,
     & (/2,4/), 
     & (/19,20/),
     & TMP, W, AMP(1,1))  !  count 1
      CALL FFV1P0_3(W(1,14),W(1,13),GC_11,ZERO,FK_ZERO,W(1,22))  !  count 11
      CALL FFV1P0_3(W(1,15),W(1,12),GC_11,ZERO,FK_ZERO,W(1,23))  !  count 11
      CALL FFV2_5_2(W(1,16),W(1,10),-GC_50,GC_58,ZERO,FK_ZERO,W(1,24))  !  count 2
      CALL FFV2_5_2(W(1,16),W(1,11),-GC_50,GC_58,ZERO,FK_ZERO,W(1,25))  !  count 2
      CALL FFV1P1N_1(W(1,2), W(1,22), GC_11, TMP(1))
      call CombineAmp(2,
     & (/1,3/), 
     & (/24,25/),
     & TMP, W, AMP(1,2))
      CALL FFV1P1N_1(W(1,2), W(1,23), GC_11, TMP(1))
      call CombineAmp(2,
     & (/2,4/), 
     & (/24,25/),
     & TMP, W, AMP(1,2))  !  count 1
      CALL FFV2_5_1(W(1,2),W(1,10),-GC_50,GC_58,ZERO,FK_ZERO,W(1,26))  !  count 4
      CALL FFV2_5_1(W(1,2),W(1,11),-GC_50,GC_58,ZERO,FK_ZERO,W(1,27))  !  count 4
      CALL FFV1P1N_2(W(1,16), W(1,22), GC_11, TMP(1))
      call CombineAmp(2,
     & (/1,3/), 
     & (/26,27/),
     & TMP, W, AMP(1,3))
      CALL FFV1P1N_2(W(1,16), W(1,23), GC_11, TMP(1))
      call CombineAmp(2,
     & (/2,4/), 
     & (/26,27/),
     & TMP, W, AMP(1,3))  !  count 1
      CALL FFV2_3_1(W(1,12),W(1,10),GC_50,GC_58,ZERO,FK_ZERO,W(1,28))  !  count 2
      CALL FFV2_3_1(W(1,12),W(1,11),GC_50,GC_58,ZERO,FK_ZERO,W(1,29))  !  count 2
      CALL FFV2_3_1(W(1,13),W(1,10),GC_50,GC_58,ZERO,FK_ZERO,W(1,30))  !  count 2
      CALL FFV2_3_1(W(1,13),W(1,11),GC_50,GC_58,ZERO,FK_ZERO,W(1,31))  !  count 2
      CALL FFV1P1N_2(W(1,14), W(1,21), GC_11, TMP(1))
      call CombineAmp(2,
     & (/1,3/), 
     & (/30,31/),
     & TMP, W, AMP(1,4))
      CALL FFV1P1N_2(W(1,15), W(1,21), GC_11, TMP(1))
      call CombineAmp(2,
     & (/2,4/), 
     & (/28,29/),
     & TMP, W, AMP(1,4))  !  count 1
      CALL FFV2_3_2(W(1,1),W(1,10),GC_50,GC_58,ZERO,FK_ZERO,W(1,32))  !  count 3
      CALL FFV2_3_2(W(1,1),W(1,11),GC_50,GC_58,ZERO,FK_ZERO,W(1,33))  !  count 3
      CALL FFV2_2(W(1,32),W(1,5),GC_100,ZERO,FK_ZERO,W(1,34))  !  count 2
      CALL FFV2_2(W(1,33),W(1,5),GC_100,ZERO,FK_ZERO,W(1,35))  !  count 2
      CALL FFV1P1N_1(W(1,2), W(1,22), GC_11, TMP(1))
      call CombineAmp(2,
     & (/1,3/), 
     & (/34,35/),
     & TMP, W, AMP(1,5))
      CALL FFV1P1N_1(W(1,2), W(1,23), GC_11, TMP(1))
      call CombineAmp(2,
     & (/2,4/), 
     & (/34,35/),
     & TMP, W, AMP(1,5))  !  count 1
      CALL FFV2_1(W(1,2),W(1,5),GC_100,ZERO,FK_ZERO,W(1,36))  !  count 9
      CALL FFV1P1N_1(W(1,36), W(1,22), GC_11, TMP(1))
      call CombineAmp(2,
     & (/1,3/), 
     & (/32,33/),
     & TMP, W, AMP(1,6))
      CALL FFV1P1N_1(W(1,36), W(1,23), GC_11, TMP(1))
      call CombineAmp(2,
     & (/2,4/), 
     & (/32,33/),
     & TMP, W, AMP(1,6))  !  count 1
      CALL FFV1P0_3(W(1,1),W(1,36),GC_11,ZERO,FK_ZERO,W(1,37))  !  count 8
      CALL FFV1P1N_1(W(1,13), W(1,37), GC_11, TMP(1))
      call CombineAmp(2,
     & (/1,3/), 
     & (/17,18/),
     & TMP, W, AMP(1,7))
      CALL FFV1P1N_1(W(1,12), W(1,37), GC_11, TMP(1))
      call CombineAmp(2,
     & (/2,4/), 
     & (/19,20/),
     & TMP, W, AMP(1,7))  !  count 1
      CALL VVV1_2(W(1,5),W(1,10),GC_53,MDL_MW,FK_MDL_WW,W(1,38))  !  count 3
      CALL VVV1_2(W(1,5),W(1,11),GC_53,MDL_MW,FK_MDL_WW,W(1,39))  !  count 3
      CALL FFV1_2(W(1,1),W(1,22),GC_11,ZERO,FK_ZERO,W(1,40))  !  count 6
      CALL FFV1_2(W(1,1),W(1,23),GC_11,ZERO,FK_ZERO,W(1,41))  !  count 6
      CALL FFV2P1N_1(W(1,2), W(1,38), GC_100, TMP(1))
      call CombineAmp(2,
     & (/1,2/), 
     & (/40,41/),
     & TMP, W, AMP(1,8))
      CALL FFV2P1N_1(W(1,2), W(1,39), GC_100, TMP(1))
      call CombineAmp(2,
     & (/3,4/), 
     & (/40,41/),
     & TMP, W, AMP(1,8))  !  count 1
      CALL FFV2_2(W(1,1),W(1,38),GC_100,ZERO,FK_ZERO,W(1,42))  !  count 2
      CALL FFV2_2(W(1,1),W(1,39),GC_100,ZERO,FK_ZERO,W(1,43))  !  count 2
      CALL FFV1P1N_1(W(1,2), W(1,22), GC_11, TMP(1))
      call CombineAmp(2,
     & (/1,3/), 
     & (/42,43/),
     & TMP, W, AMP(1,9))
      CALL FFV1P1N_1(W(1,2), W(1,23), GC_11, TMP(1))
      call CombineAmp(2,
     & (/2,4/), 
     & (/42,43/),
     & TMP, W, AMP(1,9))  !  count 1
      CALL FFV2_3P1N_1(W(1,36), W(1,10), GC_50, GC_58, TMP(1))
      call CombineAmp(2,
     & (/1,2/), 
     & (/40,41/),
     & TMP, W, AMP(1,10))
      CALL FFV2_3P1N_1(W(1,36), W(1,11), GC_50, GC_58, TMP(1))
      call CombineAmp(2,
     & (/3,4/), 
     & (/40,41/),
     & TMP, W, AMP(1,10))  !  count 1
      CALL FFV2P1N_1(W(1,26), W(1,5), GC_100, TMP(1))
      call CombineAmp(2,
     & (/1,2/), 
     & (/40,41/),
     & TMP, W, AMP(1,11))
      CALL FFV2P1N_1(W(1,27), W(1,5), GC_100, TMP(1))
      call CombineAmp(2,
     & (/3,4/), 
     & (/40,41/),
     & TMP, W, AMP(1,11))  !  count 1
      CALL FFV1P1N_2(W(1,14), W(1,37), GC_11, TMP(1))
      call CombineAmp(2,
     & (/1,3/), 
     & (/30,31/),
     & TMP, W, AMP(1,12))
      CALL FFV1P1N_2(W(1,15), W(1,37), GC_11, TMP(1))
      call CombineAmp(2,
     & (/2,4/), 
     & (/28,29/),
     & TMP, W, AMP(1,12))  !  count 1

      JAMP(:,:)  = (0d0,0d0)
      DO K = 1, NCOMB
        
        JAMP(:,:) = (0D0,0D0)
C       JAMPs contributing to orders ALL_ORDERS=1
        TMP_JAMP(6) = AMP( K,11) +  AMP( K,12)   !  used 2 times
        TMP_JAMP(5) = AMP( K,9) +  AMP( K,10)   !  used 2 times
        TMP_JAMP(4) = AMP( K,7) +  AMP( K,8)   !  used 2 times
        TMP_JAMP(3) = AMP( K,5) +  AMP( K,6)   !  used 2 times
        TMP_JAMP(2) = AMP( K,3) +  AMP( K,4)   !  used 2 times
        TMP_JAMP(1) = AMP( K,1) +  AMP( K,2)   !  used 2 times
        TMP_JAMP(9) = TMP_JAMP(6) +  TMP_JAMP(5)   !  used 2 times
        TMP_JAMP(8) = TMP_JAMP(4) +  TMP_JAMP(3)   !  used 2 times
        TMP_JAMP(7) = TMP_JAMP(2) +  TMP_JAMP(1)   !  used 2 times
        TMP_JAMP(10) = TMP_JAMP(8) +  TMP_JAMP(7)   !  used 2 times
        TMP_JAMP(11) = TMP_JAMP(10) +  TMP_JAMP(9)   !  used 2 times
        JAMP(1,1) = (-1.666666666666667D-01)*TMP_JAMP(11)
        JAMP(2,1) = (5.000000000000000D-01)*TMP_JAMP(11)

        TS(K) = 0.D0 
        DO M = 1, NAMPSO
          DO I = 1, NCOLOR
            ZTEMP = (0.D0,0.D0)
            DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J,M)
            ENDDO
            DO N = 1, NAMPSO

              TS(K) = TS(K) + ZTEMP*DCONJG(JAMP(I,N))

            ENDDO
          ENDDO
        ENDDO
        if(sde_strat.eq.1) then
        
        IF(SDE_STRAT.EQ.1)THEN
          AMP2(1)=AMP2(1)+AMP( K,1)*DCONJG(AMP( K,1))
          AMP2(7)=AMP2(7)+AMP( K,7)*DCONJG(AMP( K,7))
          AMP2(2)=AMP2(2)+AMP( K,2)*DCONJG(AMP( K,2))
          AMP2(3)=AMP2(3)+AMP( K,3)*DCONJG(AMP( K,3))
          AMP2(5)=AMP2(5)+AMP( K,5)*DCONJG(AMP( K,5))
          AMP2(6)=AMP2(6)+AMP( K,6)*DCONJG(AMP( K,6))
          AMP2(8)=AMP2(8)+AMP( K,8)*DCONJG(AMP( K,8))
          AMP2(9)=AMP2(9)+AMP( K,9)*DCONJG(AMP( K,9))
          AMP2(10)=AMP2(10)+AMP( K,10)*DCONJG(AMP( K,10))
          AMP2(11)=AMP2(11)+AMP( K,11)*DCONJG(AMP( K,11))
          AMP2(4)=AMP2(4)+AMP( K,4)*DCONJG(AMP( K,4))
          AMP2(12)=AMP2(12)+AMP( K,12)*DCONJG(AMP( K,12))
        ENDIF

        endif
        Do I = 1, NCOLOR
          DO M = 1, NAMPSO
            DO N = 1, NAMPSO

              Jamp2(i)=Jamp2(i)+DABS(DBLE(Jamp(i,m)*dconjg(Jamp(i,n))))

            enddo
          enddo
        Enddo
      ENDDO

      END


      SUBROUTINE PRINT_ZERO_AMP_23()

      integer i
      i =1
      return
      end
C     Set of functions to handle the array indices of the split orders


      INTEGER FUNCTION SQSOINDEX23(ORDERINDEXA, ORDERINDEXB)
C     
C     This functions plays the role of the interference matrix. It can
C      be hardcoded or 
C     made more elegant using hashtables if its execution speed ever
C      becomes a relevant
C     factor. From two split order indices, it return the
C      corresponding index in the squared 
C     order canonical ordering.
C     
C     CONSTANTS
C     

      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=1)
C     
C     ARGUMENTS
C     
      INTEGER ORDERINDEXA, ORDERINDEXB
C     
C     LOCAL VARIABLES
C     
      INTEGER I, SQORDERS(NSO)
      INTEGER AMPSPLITORDERS(NAMPSO,NSO)
      DATA (AMPSPLITORDERS(  1,i),i=  1,  1) /    1/
      COMMON/AMPSPLITORDERS23/AMPSPLITORDERS
C     
C     FUNCTION
C     
      INTEGER SOINDEX_FOR_SQUARED_ORDERS23
C     
C     BEGIN CODE
C     
      DO I=1,NSO
        SQORDERS(I)=AMPSPLITORDERS(ORDERINDEXA,I)
     $   +AMPSPLITORDERS(ORDERINDEXB,I)
      ENDDO
      SQSOINDEX23=SOINDEX_FOR_SQUARED_ORDERS23(SQORDERS)
      END

      INTEGER FUNCTION SOINDEX_FOR_SQUARED_ORDERS23(ORDERS)
C     
C     This functions returns the integer index identifying the squared
C      split orders list passed in argument which corresponds to the
C      values of the following list of couplings (and in this order).
C     []
C     
C     CONSTANTS
C     
      INTEGER    NSO, NSQSO, NAMPSO
      PARAMETER (NSO=1, NSQSO=1, NAMPSO=1)
C     
C     ARGUMENTS
C     
      INTEGER ORDERS(NSO)
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J
      INTEGER SQSPLITORDERS(NSQSO,NSO)
      DATA (SQSPLITORDERS(  1,i),i=  1,  1) /    2/
      COMMON/SQPLITORDERS23/SQPLITORDERS
C     
C     BEGIN CODE
C     
      DO I=1,NSQSO
        DO J=1,NSO
          IF (ORDERS(J).NE.SQSPLITORDERS(I,J)) GOTO 1009
        ENDDO
        SOINDEX_FOR_SQUARED_ORDERS23 = I
        RETURN
 1009   CONTINUE
      ENDDO

      WRITE(*,*) 'ERROR:: Stopping in function' 
      WRITE(*,*) 'SOINDEX_FOR_SQUARED_ORDERS23'
      WRITE(*,*) 'Could not find squared orders ',(ORDERS(I),I=1,NSO)
      STOP

      END

      SUBROUTINE GET_NSQSO_BORN23(NSQSO)
C     
C     Simple subroutine returning the number of squared split order
C     contributions returned when calling smatrix_split_orders 
C     

      INTEGER    NSQUAREDSO
      PARAMETER  (NSQUAREDSO=1)

      INTEGER NSQSO

      NSQSO=NSQUAREDSO

      END

C     This is the inverse subroutine of SOINDEX_FOR_SQUARED_ORDERS.
C      Not directly useful, but provided nonetheless.
      SUBROUTINE GET_SQUARED_ORDERS_FOR_SOINDEX23(SOINDEX,ORDERS)
C     
C     This functions returns the orders identified by the squared
C      split order index in argument. Order values correspond to
C      following list of couplings (and in this order):
C     []
C     
C     CONSTANTS
C     
      INTEGER    NSO, NSQSO
      PARAMETER (NSO=1, NSQSO=1)
C     
C     ARGUMENTS
C     
      INTEGER SOINDEX, ORDERS(NSO)
C     
C     LOCAL VARIABLES
C     
      INTEGER I
      INTEGER SQPLITORDERS(NSQSO,NSO)
      COMMON/SQPLITORDERS23/SQPLITORDERS	  
C     
C     BEGIN CODE
C     
      IF (SOINDEX.gt.0.and.SOINDEX.le.NSQSO) THEN
        DO I=1,NSO
          ORDERS(I) =  SQPLITORDERS(SOINDEX,I)
        ENDDO
        RETURN
      ENDIF

      WRITE(*,*) 'ERROR:: Stopping function'
     $ //' GET_SQUARED_ORDERS_FOR_SOINDEX23'
      WRITE(*,*) 'Could not find squared orders index ',SOINDEX
      STOP

      END SUBROUTINE

C     This is the inverse subroutine of getting amplitude SO orders.
C      Not directly useful, but provided nonetheless.
      SUBROUTINE GET_ORDERS_FOR_AMPSOINDEX23(SOINDEX,ORDERS)
C     
C     This functions returns the orders identified by the split order
C      index in argument. Order values correspond to following list of
C      couplings (and in this order):
C     []
C     
C     CONSTANTS
C     
      INTEGER    NSO, NAMPSO
      PARAMETER (NSO=1, NAMPSO=1)
C     
C     ARGUMENTS
C     
      INTEGER SOINDEX, ORDERS(NSO)
C     
C     LOCAL VARIABLES
C     
      INTEGER I
      INTEGER AMPSPLITORDERS(NAMPSO,NSO)
      COMMON/AMPSPLITORDERS23/AMPSPLITORDERS
C     
C     BEGIN CODE
C     
      IF (SOINDEX.gt.0.and.SOINDEX.le.NAMPSO) THEN
        DO I=1,NSO
          ORDERS(I) =  AMPSPLITORDERS(SOINDEX,I)
        ENDDO
        RETURN
      ENDIF

      WRITE(*,*) 'ERROR:: Stopping function'
     $ //' GET_ORDERS_FOR_AMPSOINDEX23'
      WRITE(*,*) 'Could not find amplitude split orders index ',SOINDEX
      STOP

      END SUBROUTINE

C     This function is not directly useful, but included for
C      completeness
      INTEGER FUNCTION SOINDEX_FOR_AMPORDERS23(ORDERS)
C     
C     This functions returns the integer index identifying the
C      amplitude split orders passed in argument which correspond to
C      the values of the following list of couplings (and in this
C      order):
C     []
C     
C     CONSTANTS
C     
      INTEGER    NSO, NAMPSO
      PARAMETER (NSO=1, NAMPSO=1)
C     
C     ARGUMENTS
C     
      INTEGER ORDERS(NSO)
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J
      INTEGER AMPSPLITORDERS(NAMPSO,NSO)
      COMMON/AMPSPLITORDERS23/AMPSPLITORDERS
C     
C     BEGIN CODE
C     
      DO I=1,NAMPSO
        DO J=1,NSO
          IF (ORDERS(J).NE.AMPSPLITORDERS(I,J)) GOTO 1009
        ENDDO
        SOINDEX_FOR_AMPORDERS23 = I
        RETURN
 1009   CONTINUE
      ENDDO

      WRITE(*,*) 'ERROR:: Stopping function SOINDEX_FOR_AMPORDERS23'
      WRITE(*,*) 'Could not find squared orders ',(ORDERS(I),I=1,NSO)
      STOP

      END

