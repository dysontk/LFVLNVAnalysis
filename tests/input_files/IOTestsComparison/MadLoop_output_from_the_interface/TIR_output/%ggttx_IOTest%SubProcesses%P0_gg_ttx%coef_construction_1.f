      SUBROUTINE ML5_0_COEF_CONSTRUCTION_1(P,NHEL,H,IC)
C     
C     Modules
C     
      USE ML5_0_POLYNOMIAL_CONSTANTS
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=4)
      INTEGER    NCOMB
      PARAMETER (NCOMB=16)
      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=3)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=44, NLOOPGROUPS=28, NCTAMPS=85)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=129)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=10,NLOOPWAVEFUNCS=93)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      REAL*16     MP__ZERO
      PARAMETER (MP__ZERO=0.0E0_16)
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=2)
C     
C     ARGUMENTS
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*16 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)

      LOGICAL DUMMYFALSE
      DATA DUMMYFALSE/.FALSE./
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
      INCLUDE 'mp_coupl.inc'

      INTEGER HELOFFSET
      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/ML5_0_FILTERS/GOODAMP,GOODHEL,HELOFFSET

      LOGICAL CHECKPHASE
      LOGICAL HELDOUBLECHECKED
      COMMON/ML5_0_INIT/CHECKPHASE, HELDOUBLECHECKED

      INTEGER SQSO_TARGET
      COMMON/ML5_0_SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/ML5_0_SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      INTEGER I_SO
      COMMON/ML5_0_I_SO/I_SO
      INTEGER I_LIB
      COMMON/ML5_0_I_LIB/I_LIB

      COMPLEX*16 AMP(NBORNAMPS)
      COMMON/ML5_0_AMPS/AMP
      COMPLEX*16 W(20,NWAVEFUNCS)
      COMMON/ML5_0_W/W

      COMPLEX*16 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*16 PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/ML5_0_WL/WL,PL

      COMPLEX*16 AMPL(3,NCTAMPS)
      COMMON/ML5_0_AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.LOOP_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

C     Coefficient construction for loop diagram with ID 4
      CALL VVV1L2P0_1(PL(0,0),W(1,5),GC_4,ZERO,ZERO,PL(0,1),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,1))
      CALL VVV1L2P0_1(PL(0,1),W(1,8),GC_4,ZERO,ZERO,PL(0,2),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1,2))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,2),2,4,1,2,1,86,H)
C     Coefficient construction for loop diagram with ID 5
      CALL FFV1L3_1(PL(0,0),W(1,3),GC_5,MDL_MT,MDL_WT,PL(0,3),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,3))
      CALL FFV1L2P0_3(PL(0,3),W(1,4),GC_5,ZERO,ZERO,PL(0,4),COEFS)
      CALL ML5_0_UPDATE_WL_1_0(WL(1,0,1,3),4,COEFS,4,4,WL(1,0,1,4))
      CALL VVV1L2P0_1(PL(0,4),W(1,5),GC_4,ZERO,ZERO,PL(0,5),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,4),4,COEFS,4,4,WL(1,0,1,5))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,5),2,4,2,1,1,87,H)
C     Coefficient construction for loop diagram with ID 6
      CALL FFV1L1P0_3(PL(0,0),W(1,3),GC_5,ZERO,ZERO,PL(0,6),COEFS)
      CALL ML5_0_UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,6))
      CALL FFV1L3_2(PL(0,6),W(1,4),GC_5,MDL_MT,MDL_WT,PL(0,7),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,6),4,COEFS,4,4,WL(1,0,1,7))
      CALL FFV1L1_2(PL(0,7),W(1,5),GC_5,MDL_MT,MDL_WT,PL(0,8),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,7),4,COEFS,4,4,WL(1,0,1,8))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,8),2,4,3,1,1,88,H)
C     Coefficient construction for loop diagram with ID 7
      CALL FFV1L1P0_3(PL(0,0),W(1,6),GC_5,ZERO,ZERO,PL(0,9),COEFS)
      CALL ML5_0_UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,9))
      CALL FFV1L3_2(PL(0,9),W(1,9),GC_5,MDL_MT,MDL_WT,PL(0,10),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,9),4,COEFS,4,4,WL(1,0,1,10))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,10),1,4,4,1,1,89,H)
C     Coefficient construction for loop diagram with ID 8
      CALL VVV1L2P0_1(PL(0,0),W(1,2),GC_4,ZERO,ZERO,PL(0,11),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,11))
      CALL FFV1L3_2(PL(0,11),W(1,4),GC_5,MDL_MT,MDL_WT,PL(0,12),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,11),4,COEFS,4,4,WL(1,0,1,12))
      CALL FFV1L1P0_3(PL(0,12),W(1,6),GC_5,ZERO,ZERO,PL(0,13),COEFS)
      CALL ML5_0_UPDATE_WL_2_0(WL(1,0,1,12),4,COEFS,4,4,WL(1,0,1,13))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,13),2,4,5,1,1,90,H)
C     Coefficient construction for loop diagram with ID 9
      CALL FFV1L2_1(PL(0,0),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,14),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,14))
      CALL FFV1L2P0_3(PL(0,14),W(1,4),GC_5,ZERO,ZERO,PL(0,15),COEFS)
      CALL ML5_0_UPDATE_WL_1_0(WL(1,0,1,14),4,COEFS,4,4,WL(1,0,1,15))
      CALL FFV1L3_1(PL(0,15),W(1,6),GC_5,MDL_MT,MDL_WT,PL(0,16),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,15),4,COEFS,4,4,WL(1,0,1,16))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,16),2,4,6,1,1,91,H)
C     Coefficient construction for loop diagram with ID 10
      CALL FFV1L1P0_3(PL(0,0),W(1,10),GC_5,ZERO,ZERO,PL(0,17),COEFS)
      CALL ML5_0_UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,17))
      CALL FFV1L3_2(PL(0,17),W(1,7),GC_5,MDL_MT,MDL_WT,PL(0,18),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,17),4,COEFS,4,4,WL(1,0,1,18))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,18),1,4,7,1,1,92,H)
C     Coefficient construction for loop diagram with ID 11
      CALL FFV1L3_1(PL(0,11),W(1,3),GC_5,MDL_MT,MDL_WT,PL(0,19),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,11),4,COEFS,4,4,WL(1,0,1,19))
      CALL FFV1L2P0_3(PL(0,19),W(1,7),GC_5,ZERO,ZERO,PL(0,20),COEFS)
      CALL ML5_0_UPDATE_WL_2_0(WL(1,0,1,19),4,COEFS,4,4,WL(1,0,1,20))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,20),2,4,8,1,1,93,H)
C     Coefficient construction for loop diagram with ID 12
      CALL FFV1L1_2(PL(0,0),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,21),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,21))
      CALL FFV1L1P0_3(PL(0,21),W(1,3),GC_5,ZERO,ZERO,PL(0,22),COEFS)
      CALL ML5_0_UPDATE_WL_1_0(WL(1,0,1,21),4,COEFS,4,4,WL(1,0,1,22))
      CALL FFV1L3_2(PL(0,22),W(1,7),GC_5,MDL_MT,MDL_WT,PL(0,23),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,22),4,COEFS,4,4,WL(1,0,1,23))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,23),2,4,9,1,1,94,H)
C     Coefficient construction for loop diagram with ID 13
      CALL VVV1L2P0_1(PL(0,0),W(1,1),GC_4,ZERO,ZERO,PL(0,24),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,24))
      CALL FFV1L3_2(PL(0,24),W(1,4),GC_5,MDL_MT,MDL_WT,PL(0,25),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,24),4,COEFS,4,4,WL(1,0,1,25))
      CALL FFV1L1P0_3(PL(0,25),W(1,10),GC_5,ZERO,ZERO,PL(0,26),COEFS)
      CALL ML5_0_UPDATE_WL_2_0(WL(1,0,1,25),4,COEFS,4,4,WL(1,0,1,26))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,26),2,4,10,1,1,95,H)
C     Coefficient construction for loop diagram with ID 14
      CALL FFV1L3_1(PL(0,24),W(1,3),GC_5,MDL_MT,MDL_WT,PL(0,27),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,24),4,COEFS,4,4,WL(1,0,1,27))
      CALL FFV1L2P0_3(PL(0,27),W(1,9),GC_5,ZERO,ZERO,PL(0,28),COEFS)
      CALL ML5_0_UPDATE_WL_2_0(WL(1,0,1,27),4,COEFS,4,4,WL(1,0,1,28))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,28),2,4,11,1,1,96,H)
C     Coefficient construction for loop diagram with ID 15
      CALL VVV1L2P0_1(PL(0,24),W(1,2),GC_4,ZERO,ZERO,PL(0,29),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,24),4,COEFS,4,4,WL(1,0,1,29))
      CALL FFV1L3_2(PL(0,29),W(1,4),GC_5,MDL_MT,MDL_WT,PL(0,30),COEFS)
      CALL ML5_0_UPDATE_WL_2_1(WL(1,0,1,29),4,COEFS,4,4,WL(1,0,1,30))
      CALL FFV1L1P0_3(PL(0,30),W(1,3),GC_5,ZERO,ZERO,PL(0,31),COEFS)
      CALL ML5_0_UPDATE_WL_3_0(WL(1,0,1,30),4,COEFS,4,4,WL(1,0,1,31))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,31),3,4,12,1,1,97,H)
C     Coefficient construction for loop diagram with ID 16
      CALL FFV1L3_1(PL(0,29),W(1,3),GC_5,MDL_MT,MDL_WT,PL(0,32),COEFS)
      CALL ML5_0_UPDATE_WL_2_1(WL(1,0,1,29),4,COEFS,4,4,WL(1,0,1,32))
      CALL FFV1L2P0_3(PL(0,32),W(1,4),GC_5,ZERO,ZERO,PL(0,33),COEFS)
      CALL ML5_0_UPDATE_WL_3_0(WL(1,0,1,32),4,COEFS,4,4,WL(1,0,1,33))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,33),3,4,13,1,1,98,H)
C     Coefficient construction for loop diagram with ID 17
      CALL VVV1L2P0_1(PL(0,29),W(1,8),GC_4,ZERO,ZERO,PL(0,34),COEFS)
      CALL ML5_0_UPDATE_WL_2_1(WL(1,0,1,29),4,COEFS,4,4,WL(1,0,1,34))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,34),3,4,14,1,1,99,H)
C     Coefficient construction for loop diagram with ID 18
      CALL VVVV1L2P0_1(PL(0,24),W(1,8),W(1,2),GC_6,ZERO,ZERO,PL(0,35)
     $ ,COEFS)
      CALL ML5_0_UPDATE_WL_1_0(WL(1,0,1,24),4,COEFS,4,4,WL(1,0,1,35))
      CALL VVVV3L2P0_1(PL(0,24),W(1,8),W(1,2),GC_6,ZERO,ZERO,PL(0,36)
     $ ,COEFS)
      CALL ML5_0_UPDATE_WL_1_0(WL(1,0,1,24),4,COEFS,4,4,WL(1,0,1,36))
      CALL VVVV4L2P0_1(PL(0,24),W(1,8),W(1,2),GC_6,ZERO,ZERO,PL(0,37)
     $ ,COEFS)
      CALL ML5_0_UPDATE_WL_1_0(WL(1,0,1,24),4,COEFS,4,4,WL(1,0,1,37))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,35),1,4,15,2,1,100,H)
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,36),1,4,15,2,1,101,H)
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,37),1,4,15,2,1,102,H)
C     Coefficient construction for loop diagram with ID 19
      CALL FFV1L2_1(PL(0,27),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,38),COEFS)
      CALL ML5_0_UPDATE_WL_2_1(WL(1,0,1,27),4,COEFS,4,4,WL(1,0,1,38))
      CALL FFV1L2P0_3(PL(0,38),W(1,4),GC_5,ZERO,ZERO,PL(0,39),COEFS)
      CALL ML5_0_UPDATE_WL_3_0(WL(1,0,1,38),4,COEFS,4,4,WL(1,0,1,39))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,39),3,4,16,1,1,103,H)
C     Coefficient construction for loop diagram with ID 20
      CALL FFV1L2_1(PL(0,0),W(1,1),GC_5,MDL_MT,MDL_WT,PL(0,40),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,40))
      CALL FFV1L2P0_3(PL(0,40),W(1,4),GC_5,ZERO,ZERO,PL(0,41),COEFS)
      CALL ML5_0_UPDATE_WL_1_0(WL(1,0,1,40),4,COEFS,4,4,WL(1,0,1,41))
      CALL FFV1L3_1(PL(0,41),W(1,10),GC_5,MDL_MT,MDL_WT,PL(0,42),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,41),4,COEFS,4,4,WL(1,0,1,42))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,42),2,4,17,1,1,104,H)
C     Coefficient construction for loop diagram with ID 21
      CALL FFV1L1_2(PL(0,0),W(1,1),GC_5,MDL_MT,MDL_WT,PL(0,43),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,43))
      CALL FFV1L1P0_3(PL(0,43),W(1,3),GC_5,ZERO,ZERO,PL(0,44),COEFS)
      CALL ML5_0_UPDATE_WL_1_0(WL(1,0,1,43),4,COEFS,4,4,WL(1,0,1,44))
      CALL FFV1L3_2(PL(0,44),W(1,9),GC_5,MDL_MT,MDL_WT,PL(0,45),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,44),4,COEFS,4,4,WL(1,0,1,45))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,45),2,4,18,1,1,105,H)
C     Coefficient construction for loop diagram with ID 22
      CALL VVVV1L2P0_1(PL(0,11),W(1,8),W(1,1),GC_6,ZERO,ZERO,PL(0,46)
     $ ,COEFS)
      CALL ML5_0_UPDATE_WL_1_0(WL(1,0,1,11),4,COEFS,4,4,WL(1,0,1,46))
      CALL VVVV3L2P0_1(PL(0,11),W(1,8),W(1,1),GC_6,ZERO,ZERO,PL(0,47)
     $ ,COEFS)
      CALL ML5_0_UPDATE_WL_1_0(WL(1,0,1,11),4,COEFS,4,4,WL(1,0,1,47))
      CALL VVVV4L2P0_1(PL(0,11),W(1,8),W(1,1),GC_6,ZERO,ZERO,PL(0,48)
     $ ,COEFS)
      CALL ML5_0_UPDATE_WL_1_0(WL(1,0,1,11),4,COEFS,4,4,WL(1,0,1,48))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,46),1,4,19,2,1,106,H)
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,47),1,4,19,2,1,107,H)
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,48),1,4,19,2,1,108,H)
C     Coefficient construction for loop diagram with ID 23
      CALL VVV1L2P0_1(PL(0,44),W(1,2),GC_4,ZERO,ZERO,PL(0,49),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,44),4,COEFS,4,4,WL(1,0,1,49))
      CALL FFV1L3_2(PL(0,49),W(1,4),GC_5,MDL_MT,MDL_WT,PL(0,50),COEFS)
      CALL ML5_0_UPDATE_WL_2_1(WL(1,0,1,49),4,COEFS,4,4,WL(1,0,1,50))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,50),3,4,20,1,1,109,H)
C     Coefficient construction for loop diagram with ID 24
      CALL FFV1L2_1(PL(0,40),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,51),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,40),4,COEFS,4,4,WL(1,0,1,51))
      CALL FFV1L2P0_3(PL(0,51),W(1,4),GC_5,ZERO,ZERO,PL(0,52),COEFS)
      CALL ML5_0_UPDATE_WL_2_0(WL(1,0,1,51),4,COEFS,4,4,WL(1,0,1,52))
      CALL FFV1L3_1(PL(0,52),W(1,3),GC_5,MDL_MT,MDL_WT,PL(0,53),COEFS)
      CALL ML5_0_UPDATE_WL_2_1(WL(1,0,1,52),4,COEFS,4,4,WL(1,0,1,53))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,53),3,4,21,1,1,110,H)
C     Coefficient construction for loop diagram with ID 25
      CALL FFV1L1_2(PL(0,43),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,54),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,43),4,COEFS,4,4,WL(1,0,1,54))
      CALL FFV1L1P0_3(PL(0,54),W(1,3),GC_5,ZERO,ZERO,PL(0,55),COEFS)
      CALL ML5_0_UPDATE_WL_2_0(WL(1,0,1,54),4,COEFS,4,4,WL(1,0,1,55))
      CALL FFV1L3_2(PL(0,55),W(1,4),GC_5,MDL_MT,MDL_WT,PL(0,56),COEFS)
      CALL ML5_0_UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,56))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,56),3,4,22,1,1,111,H)
C     Coefficient construction for loop diagram with ID 26
      CALL VVVV1L2P0_1(PL(0,0),W(1,2),W(1,1),GC_6,ZERO,ZERO,PL(0,57)
     $ ,COEFS)
      CALL ML5_0_UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,57))
      CALL VVVV3L2P0_1(PL(0,0),W(1,2),W(1,1),GC_6,ZERO,ZERO,PL(0,58)
     $ ,COEFS)
      CALL ML5_0_UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,58))
      CALL VVVV4L2P0_1(PL(0,0),W(1,2),W(1,1),GC_6,ZERO,ZERO,PL(0,59)
     $ ,COEFS)
      CALL ML5_0_UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,59))
      CALL VVV1L2P0_1(PL(0,57),W(1,8),GC_4,ZERO,ZERO,PL(0,60),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,57),4,COEFS,4,4,WL(1,0,1,60))
      CALL VVV1L2P0_1(PL(0,58),W(1,8),GC_4,ZERO,ZERO,PL(0,61),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,58),4,COEFS,4,4,WL(1,0,1,61))
      CALL VVV1L2P0_1(PL(0,59),W(1,8),GC_4,ZERO,ZERO,PL(0,62),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,59),4,COEFS,4,4,WL(1,0,1,62))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,60),1,4,23,2,1,112,H)
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,61),1,4,23,2,1,113,H)
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,62),1,4,23,2,1,114,H)
C     Coefficient construction for loop diagram with ID 27
      CALL VVVV1L2P0_1(PL(0,4),W(1,2),W(1,1),GC_6,ZERO,ZERO,PL(0,63)
     $ ,COEFS)
      CALL ML5_0_UPDATE_WL_1_0(WL(1,0,1,4),4,COEFS,4,4,WL(1,0,1,63))
      CALL VVVV3L2P0_1(PL(0,4),W(1,2),W(1,1),GC_6,ZERO,ZERO,PL(0,64)
     $ ,COEFS)
      CALL ML5_0_UPDATE_WL_1_0(WL(1,0,1,4),4,COEFS,4,4,WL(1,0,1,64))
      CALL VVVV4L2P0_1(PL(0,4),W(1,2),W(1,1),GC_6,ZERO,ZERO,PL(0,65)
     $ ,COEFS)
      CALL ML5_0_UPDATE_WL_1_0(WL(1,0,1,4),4,COEFS,4,4,WL(1,0,1,65))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,63),1,4,24,1,1,115,H)
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,64),1,4,24,1,1,116,H)
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,65),1,4,24,1,1,117,H)
C     Coefficient construction for loop diagram with ID 28
      CALL GHGHGL2_1(PL(0,0),W(1,5),GC_4,ZERO,ZERO,PL(0,66),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),1,COEFS,1,1,WL(1,0,1,66))
      CALL GHGHGL2_1(PL(0,66),W(1,8),GC_4,ZERO,ZERO,PL(0,67),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,66),1,COEFS,1,1,WL(1,0,1,67))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,67),2,1,1,1,1,118,H)
C     Coefficient construction for loop diagram with ID 29
      CALL GHGHGL1_2(PL(0,0),W(1,1),GC_4,ZERO,ZERO,PL(0,68),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),1,COEFS,1,1,WL(1,0,1,68))
      CALL GHGHGL1_2(PL(0,68),W(1,2),GC_4,ZERO,ZERO,PL(0,69),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,68),1,COEFS,1,1,WL(1,0,1,69))
      CALL GHGHGL1_2(PL(0,69),W(1,8),GC_4,ZERO,ZERO,PL(0,70),COEFS)
      CALL ML5_0_UPDATE_WL_2_1(WL(1,0,1,69),1,COEFS,1,1,WL(1,0,1,70))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,70),3,1,14,1,1,119,H)
C     Coefficient construction for loop diagram with ID 30
      CALL GHGHGL2_1(PL(0,0),W(1,1),GC_4,ZERO,ZERO,PL(0,71),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),1,COEFS,1,1,WL(1,0,1,71))
      CALL GHGHGL2_1(PL(0,71),W(1,2),GC_4,ZERO,ZERO,PL(0,72),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,71),1,COEFS,1,1,WL(1,0,1,72))
      CALL GHGHGL2_1(PL(0,72),W(1,8),GC_4,ZERO,ZERO,PL(0,73),COEFS)
      CALL ML5_0_UPDATE_WL_2_1(WL(1,0,1,72),1,COEFS,1,1,WL(1,0,1,73))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,73),3,1,14,1,1,120,H)
C     Coefficient construction for loop diagram with ID 31
      CALL FFV1L2_1(PL(0,0),W(1,5),GC_5,ZERO,ZERO,PL(0,74),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,74))
      CALL FFV1L2_1(PL(0,74),W(1,8),GC_5,ZERO,ZERO,PL(0,75),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,74),4,COEFS,4,4,WL(1,0,1,75))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,75),2,4,1,1,4,121,H)
C     Coefficient construction for loop diagram with ID 32
      CALL FFV1L1_2(PL(0,0),W(1,1),GC_5,ZERO,ZERO,PL(0,76),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,76))
      CALL FFV1L1_2(PL(0,76),W(1,2),GC_5,ZERO,ZERO,PL(0,77),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1,77))
      CALL FFV1L1_2(PL(0,77),W(1,8),GC_5,ZERO,ZERO,PL(0,78),COEFS)
      CALL ML5_0_UPDATE_WL_2_1(WL(1,0,1,77),4,COEFS,4,4,WL(1,0,1,78))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,78),3,4,14,1,4,122,H)
C     Coefficient construction for loop diagram with ID 33
      CALL FFV1L2_1(PL(0,0),W(1,1),GC_5,ZERO,ZERO,PL(0,79),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,79))
      CALL FFV1L2_1(PL(0,79),W(1,2),GC_5,ZERO,ZERO,PL(0,80),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,79),4,COEFS,4,4,WL(1,0,1,80))
      CALL FFV1L2_1(PL(0,80),W(1,8),GC_5,ZERO,ZERO,PL(0,81),COEFS)
      CALL ML5_0_UPDATE_WL_2_1(WL(1,0,1,80),4,COEFS,4,4,WL(1,0,1,81))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,81),3,4,14,1,4,123,H)
C     Coefficient construction for loop diagram with ID 34
      CALL FFV1L2_1(PL(0,0),W(1,5),GC_5,MDL_MB,ZERO,PL(0,82),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,82))
      CALL FFV1L2_1(PL(0,82),W(1,8),GC_5,MDL_MB,ZERO,PL(0,83),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,82),4,COEFS,4,4,WL(1,0,1,83))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,83),2,4,25,1,1,124,H)
C     Coefficient construction for loop diagram with ID 35
      CALL FFV1L1_2(PL(0,0),W(1,1),GC_5,MDL_MB,ZERO,PL(0,84),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,84))
      CALL FFV1L1_2(PL(0,84),W(1,2),GC_5,MDL_MB,ZERO,PL(0,85),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,84),4,COEFS,4,4,WL(1,0,1,85))
      CALL FFV1L1_2(PL(0,85),W(1,8),GC_5,MDL_MB,ZERO,PL(0,86),COEFS)
      CALL ML5_0_UPDATE_WL_2_1(WL(1,0,1,85),4,COEFS,4,4,WL(1,0,1,86))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,86),3,4,26,1,1,125,H)
C     Coefficient construction for loop diagram with ID 36
      CALL FFV1L2_1(PL(0,0),W(1,1),GC_5,MDL_MB,ZERO,PL(0,87),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,87))
      CALL FFV1L2_1(PL(0,87),W(1,2),GC_5,MDL_MB,ZERO,PL(0,88),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,87),4,COEFS,4,4,WL(1,0,1,88))
      CALL FFV1L2_1(PL(0,88),W(1,8),GC_5,MDL_MB,ZERO,PL(0,89),COEFS)
      CALL ML5_0_UPDATE_WL_2_1(WL(1,0,1,88),4,COEFS,4,4,WL(1,0,1,89))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,89),3,4,26,1,1,126,H)
C     Coefficient construction for loop diagram with ID 37
      CALL FFV1L2_1(PL(0,0),W(1,5),GC_5,MDL_MT,MDL_WT,PL(0,90),COEFS)
      CALL ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,90))
      CALL FFV1L2_1(PL(0,90),W(1,8),GC_5,MDL_MT,MDL_WT,PL(0,91),COEFS)
      CALL ML5_0_UPDATE_WL_1_1(WL(1,0,1,90),4,COEFS,4,4,WL(1,0,1,91))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,91),2,4,27,1,1,127,H)
C     Coefficient construction for loop diagram with ID 38
      CALL FFV1L1_2(PL(0,54),W(1,8),GC_5,MDL_MT,MDL_WT,PL(0,92),COEFS)
      CALL ML5_0_UPDATE_WL_2_1(WL(1,0,1,54),4,COEFS,4,4,WL(1,0,1,92))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,92),3,4,28,1,1,128,H)
C     Coefficient construction for loop diagram with ID 39
      CALL FFV1L2_1(PL(0,51),W(1,8),GC_5,MDL_MT,MDL_WT,PL(0,93),COEFS)
      CALL ML5_0_UPDATE_WL_2_1(WL(1,0,1,51),4,COEFS,4,4,WL(1,0,1,93))
      CALL ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,93),3,4,28,1,1,129,H)
C     At this point, all loop coefficients needed for (QCD=6), i.e. of
C      split order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 4000

      GOTO 1001
 4000 CONTINUE
      LOOP_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

