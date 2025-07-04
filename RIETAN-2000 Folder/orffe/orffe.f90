      PROGRAM ORFFE
C**********************************************************************C
C     OR FFE. FORTRAN CRYSTALLOGRAPHIC FUNCTION AND ERROR PROGRAM      C
C     ORIGINAL WRITTEN BY BUSING, MARTIN AND LEVY (1964)               C
C     EINE LOKALE VERSION IN KOMBINATION MIT DEM FMLS-PROGRAMM         C
C     MV13 (13.01.83, K. KATO)                                         C
C     UEBERSETZUNG MIT (AUTODBL(DBLPAD),ALC)-OPTION                    C
C     FURTHER MODIFIED BY F. IZUMI, 1988.7.8                           C
C**********************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NCS=512,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      COMMON /CN/ NEWNO(NCS),C1(NCS),C2(10,NCS),INP(10,NCS),NCNPAR(NCS),
     &NCNSTR

      COMMON /LG/ LAUEG
      COMMON /NEWC/ ATOMN
      COMMON /FLAG1/ LWR
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
*     CHARACTER TITLE*72,ATOMN((NNP-14)/10)*9,IIS(2,3,NW)*2
      CHARACTER TITLE*25,ATOMN((NNP-14)/10)*9,IIS(2,3,NW)*2
      CHARACTER IXYZ(7)*2
      CHARACTER TRBUF(3)*20
      CHARACTER*50 FILE9
      CHARACTER EOC*6, EOR*5, EOS*6
*     FOR g77
      CHARACTER*150 LINPAR
      DATA IXYZ/'-Z','-Y','-X','  ','+X','+Y','+Z'/

*     FLAG TO INDICATE THAT DISTANCE +OR- ERROR IS TO BE WRITTEN
      LWR=1
* Changed on 2000.11.28 Izumi {
      CALL INITIAL
*     FILE #6 WAS CHANGED INTO #3 
      OPEN(UNIT=3,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     &  STATUS='SCRATCH')
* }
      REWIND(UNIT=9)
      READ(9,'(A)') TITLE
      WRITE(3,'(/A)') ' INTERATOMIC DISTANCES AND ANGLES IN '//
     &TITLE
*     READ(5,202) IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB
 0202 FORMAT (24I3)
*     MAXIMUM NUMBER OF SYMMETRY OPERATIONS
      NS=NW
*     FIXED VALUES ARE ASSIGNED
      IPM=1
      NP=0
      IAM=1
      IF(NP)00297,00294,00297
* 294 WRITE(3,295)
* 295 FORMAT('0INPUT DATA TO BE READ FROM A SCRATCH FILE')
  294 CONTINUE
      GO TO 00301
  297 WRITE(3,298) NP
00298    FORMAT ('0INPUT DATA TO BE READ FROM CARDS'/
     1   '0NUMBER OF STRUCTURE PARAMETERS IS',I4)
* 301 WRITE(3,302)
*0302 FORMAT (68H0VARIANCE-COVARIANCE MATRIX AND PARAMETER SELECTION INF
*    1ORMATION WILL)
  301 CONTINUE
      IF(IPM)00308,00305,00308
* 305 WRITE(3,306)
*0306    FORMAT (1H+,T69,12H NOT BE USED)
  305 CONTINUE
         GO TO 00310
* 308 WRITE(3,309)
*0309    FORMAT (1H+,T69,8H BE USED)
  308 CONTINUE
      READ(9,*) LAUEG
      I=0
 3081 CONTINUE
         READ(9,'(A)') EOC
         IF (EOC .EQ. 'ENDCON') GO TO 3082
         I=I+1
         BACKSPACE 9
         READ(9,*) NEWNO(I),NCNPAR(I),C1(I),
     &   (C2(J,I),INP(J,I),J=1,NCNPAR(I))
      GO TO 3081
 3082 NCNSTR=I
*     READ THE NUMBER OF PHASES
      BACKSPACE 9
      READ(9,'(6X,I1)') NPHASE
*     WRITE(3,'(A,I3)') '0LAUE GROUP NO. IS',LAUEG
* 310 WRITE(3,311) NS
*0311 FORMAT (28H0NUMBER OF SYMMETRY CARDS ISI3)
  310 CONTINUE
*     WRITE(3,314)
*0314 FORMAT (26H0CELL PARAMETER ERRORS ARE)
*     IF(IAM-1)00316,00319,00322
* 316 WRITE(3,317)
*0317    FORMAT (1H+,T27,15H NOT TO BE USED)
*        GO TO 00401
* 319 WRITE(3,320)
*0320    FORMAT (1H+,T27,42H TO BE READ IN THE FORM OF STANDARD ERRORS)
*        GO TO 00401
* 322 WRITE(3,323)
*0323    FORMAT (1H+,T27,37H TO BE READ IN THE FORM OF A VARIANCE,
*    1      18H-COVARIANCE MATRIX)
00401 IF(NP)01001,01001,00501
C        READ PARAMETERS AND VARIANCE-COVARIANCE MATRIX FROM CARDS
* 501 WRITE(3,502) NV
*0502    FORMAT (44H0ORDER OF THE VARIANCE-COVARIANCE MATRIX FOR
*    1      28H THE STRUCTURE PARAMETERS IS,I4)
  501 CONTINUE
         IF(JXP)00509,00509,00505
* 505 WRITE(3,506) JXP,JX
*0506       FORMAT (41H0PERIOD OF THE POSITION PARAMETERS IN THE
*    1         18H PARAMETER LIST ISI3/22H0POSITION OF THE FIRST
*    2         38H X COORDINATE IN THE PARAMETER LIST ISI3)
  505 CONTINUE
00509    ITF=0
         IF(JBP)00601,00601,00510
* 510 WRITE(3,511) JBP,JB
*0511       FORMAT (46H0PERIOD OF THE TEMPERATURE FACTOR COEFFICIENTS
*    1         25H IN THE PARAMETER LIST ISI3/16H0POSITION OF THE
*    2         44H FIRST TEMPERATURE FACTOR COEFFICIENT IN THE
*    3         18H PARAMETER LIST ISI3)
  510 CONTINUE
  601 READ(9,603) (P(I),I=1,NP)
00603    FORMAT (8F9.4)
         IF(IPM)00801,00907,00801
  801 READ(9,802) (KI1(I),I=1,NP)
00802       FORMAT (72I1)
            NM=(NV*(NV+1))/2
      READ(9,905) SCALE
      READ(9,905) (PM(I),I=1,NM)
  905 FORMAT(6E12.6)
      DO 906 I=1,NM
  906 PM(I)=PM(I)*SCALE
00907 WRITE(3,'(''0PARAMETERS AND KEY-INTEGERS''/1H )')
      WRITE(3,'((1H ,10(F11.6,I2)))') (P(I),KI1(I),I=1,NP)
      GO TO 2115
C        READ PARAMETERS AND VARIANCE-COVARIANCE MATRIX
C        FROM OR FLS TAPE
*1001 READ(9,'(3I6)') ITF,NQ,NP
 1001 READ(9,*) NP
*     WRITE(3,1004) NP
*1004    FORMAT (34H0NUMBER OF STRUCTURE PARAMETERS ISI4)
      DO 1301 J=1,NP/10
         READ(9,'(A9,1X,10F14.0)') ATOMN(J),
     &   (P(I),I=10*(J-1)+1,10*(J-1)+10)
 1301 CONTINUE
*        JX=NQ+4
*        JB=NQ+7
*        IF(ITF-1)01301,01301,01401
*1301 JXP=5
*     JBP=5
*           GO TO 01501
*1401 JXP=10
*     JBP=10
*     READ THE PARAMETER FILE
      ITF=2
      NQ=-2
      JXP=10
      JX=2
      JBP=10
      JB=5
*1501 WRITE(3,506) JXP,JX
*     WRITE(3,511) JBP,JB
 1501 CONTINUE
*     READ(9,'(14A)') (ATOMN(I),I=1,(NP-NQ-2)/JXP)
         IF(IPM)01701,02101,01701
 1701 READ(9,'(150I1)') (KI1(I),I=1,NP)
      READ(9,*) NV
*     WRITE(3,502) NV
      NM=(NV*(NV+1))/2
      READ(9,*) (PM(I),I=1,NM)
C     PUT OUT INPUT PARAMETERS
*2101 WRITE(3,2102)
*2102 FORMAT(11H0INPUT DATA/19H0GENERAL PARAMETERS/1H0,17X,7HP    KI/1H
*    1)
 2101 CONTINUE
      IEND=NQ+2
      DO 2104 I=1,IEND
 2104 WRITE(3,2105) P(I),KI1(I)
 2105 FORMAT(13X,F10.6,I2)
      WRITE(3,2111) NPHASE,
     &  ' A  ATOM      OCPNCY  KI     X    KI     Y    KI     Z    '//
     &  'KI   B/B11  KI    B22   KI    B33   KI    B12   KI    B13   KI'
     &  //'    B23   KI'
 2111 FORMAT(/'0ATOMIC PARAMETERS FOR PHASE #',I1/1H0,A/1H )
      DO 2114 J=1,(NP-NQ-2)/JXP
      IANF=IEND+1
      IEND=IEND+JXP
      WRITE(3,2113) J,ATOMN(J),(P(I),KI1(I),I=IANF,IEND)
*2113 FORMAT(1H ,I3,2X,A,F9.6,I2,9(F10.6,I2))
 2113 FORMAT(1H ,I2,2X,A,F9.6,I2,9(F10.6,I2))
 2114 CONTINUE
C     READ AND PUT OUT CELL PARAMETERS
*2115 READ(5,2202) (A(I),I=1,6)
 2115 READ(9,*) (A(I),I=1,6)
*2202 FORMAT (6F9.4)
      WRITE(3,2205) (A(I),I=1,6)
*2205 FORMAT(16H0CELL PARAMETERS/1H06F11.4)
 2205 FORMAT(/'0CELL PARAMETERS: A, B, C, COS(ALPHA), COS(BETA), COS(GAM
     &MA)'/1H0,6F12.5)
      IF(IAM-1)03401,02401,02701
C        READ STANDARD ERRORS OF CELL PARAMETERS
02401    DO 02402 I=1,21
02402       AM(I)=0.0
*     READ(5,2202) AM(1),AM(7),AM(12),AM(16),AM(19),AM(21)
      READ(9,*) AM(1),AM(7),AM(12),AM(16),AM(19),AM(21)
      WRITE(3,2503) AM(1),AM(7),AM(12),AM(16),AM(19),AM(21)
*2503    FORMAT (49H0STANDARD ERRORS, RESPECTIVELY, OF THE ABOVE CELL
*    1      11H PARAMETERS/1H06F11.4)
 2503 FORMAT ('0STANDARD ERRORS, RESPECTIVELY, OF THE ABOVE CELL PARAMET
     1ERS'/1H0,6F12.5)
         DO 02602 I=1,21
02602       AM(I)=AM(I)*AM(I)
         GO TO 03201
C        READ VARIANCE-COVARIANCE MATRIX FOR CELL PARAMETERS
C2701 READ(5,603) (AM(I),I=1,21)
 2701 READ(5,2702) (AM(I),I=1,21)
 2702 FORMAT(6E9.4/9X,5E9.4/18X,4E9.4/27X,3E9.4/36X,2E9.4/45X,E9.4)
      WRITE(3,2703)
02703    FORMAT (47H0VARIANCE-COVARIANCE MATRIX FOR CELL PARAMETERS)
         IJ=1
         DO 03101 I=1,6
            DO 02902 J=1,6
02902          ROW(J)=0.0
            DO 03003 J=I,6
               ROW(J)=AM(IJ)
03003          IJ=IJ+1
 3101    WRITE(3,3102) (ROW(J),J=1,6)
03102       FORMAT (1H0,6E11.4)
C        COMPUTE CELL PARAMETER INCREMENTS USED TO OBTAIN DERIVATIVES
03201    K=1
         L=6
         DO 03303 I=1,6
      DA(I)=(0.01)*SQRT(AM(K))
            K=K+L
03303       L=L-1
03401 IF(NS)03701,03701,03601
C        READ AND PUT OUT SYMMETRY TRANSFORMATIONS
*3601 READ(5,3602) ((TS(I,J),(IS(K,I,J),K=1,2),I=1,3),J=1,NS)
*3601 READ(9,3602,IOSTAT=IOS)
*    &((TS(I,J),(IS(K,I,J),K=1,2),I=1,3),J=1,NS)
 3601 DO 3603 J=1,NS
         READ(9,'(A)') EOS
         IF (EOS .EQ. 'ENDSYM') GO TO 999
         BACKSPACE 9
         READ(9,3602) (TS(I,J),(IS(K,I,J),K=1,2),I=1,3)
 3602    FORMAT (F11.6,2I2,F11.6,2I2,F11.6,2I2)
 3603 CONTINUE
  999 NS=J-1
      WRITE(3,3604)
 3604 FORMAT(/'0SYMMETRY INFORMATION'/'0  S                X''
     &                      Y''                      Z'''/1H )
      DO 3605 J=1,NS
      DO 3605 I=1,3
      DO 3605 K=1,2
 3605 IIS(K,I,J)=IXYZ(IS(K,I,J)+4)
 
*     WRITE(3,3608) (J,(TS(I,J),(IIS(K,I,J),K=1,2),I=1,3),J=1,NS)
*3608 FORMAT(1H ,I3,F20.6,2A2,F20.6,2A2,F20.6,2A2)
      DO J = 1, NS
         DO I = 1, 3
	    IF (TS(I,J) .EQ. 0.0) THEN
*              TO DELETE '0.000000'
	       WRITE(TRBUF(I),'(A)') ' '
	    ELSE
	       WRITE(TRBUF(I),'(F20.6)') TS(I,J)
	    END IF
	 END DO
	 WRITE(3,'(1H ,I3,3(A20,2A2))') 
     &   J,(TRBUF(I),(IIS(K,I,J),K=1,2),I=1,3)
      END DO
      
      WRITE(3,3690)
 3690 FORMAT(/'0DESIGNATION OF ATOMS:'/'0AN ATOM IS DESIGNATED BY A PAIR
     * OF INTEGER NUMBERS IN PARENTHESES (A,1000*C+S).'/' A IS THE NUMBE
     *R OF THE ATOM IN THE LIST OF ATOMIC PARAMETERS.'/' S IS THE NUMBER
     * OF THE SYMMETRY INFORMATION TO BE APPLIED, AND C DEFINES THE UNIT
     *-CELL TRANSLATIONS AS DESCRIBED BELOW.'/' THE PROGRAMM OBTAINS THE
     * COORDINATES OF AN ATOM IN THE FOLLOWING WAY:'/'0   1. THE INTEGER
     * A IS USED TO COMPUTE THE LOCATION OF THE COORDINATES IN THE PARAM
     *ETER LIST AND X, Y, Z ARE PICKED UP.'/'       IF A=0 THE PROGRAM S
     *ETS X=Y=Z=0.'/'    2. THESE COORDINATES ARE THEN TRANSFORMED TO X'
     *', Y'', Z'' ACCORDING TO THE S-TH SYMMETRY INFORMATION.'/'       I
     *F S=0 NO TRANSFORMATION IS MADE.'/'    3. THE CELL TRANSLATIONS AR
     *E THEN MADE ACCORDING TO THE FOLLOWING TABLE:'/)

      WRITE(3,3695)
      WRITE(3,3696)
      WRITE(3,3697)
      WRITE(3,3698)
      WRITE(3,3699)
 3695 FORMAT('0',T10,'C',T20,'X''''',T30,'Y''''',T40,'Z'''''/'0',T10,
     *'0',T20,'X''',T30,'Y''',T40,'Z'''/' ',T10,'1',T20,'X''-1',T30,
     *'Y''',T40,'Z'''/' ',T10,'2',T20,'X''-1',T30,'Y''+1',T40,'Z'''/
     *' ',T10,'3',T20,'X''-1',T30,'Y''',T40,'Z''+1'/' ',T10,'4',T20,
     *'X''-1',T30,'Y''+1',T40,'Z''+1'/' ',T10,'5',T20,'X''',T30,
     *'Y''-1',T40,'Z'''/' ',T10,'6',T20,'X''+1',T30,'Y''-1',T40,
     *'Z'''/' ',T10,'7',T20,'X''',T30,'Y''-1',T40,'Z''+1'/' ',T10,
     *'8',T20,'X''+1',T30,'Y''-1',T40,'Z''+1'/' ',T10,'9',T20,'
     *X''-1',T30,'Y''-1',T40,'Z'''/' ',T9,'10',T20,'X''-1',T30,'Y''-1',
     *T40,'Z''+1'/' ',T9,'11',T20,'X''',T30,'Y''',T40,'Z''-1'/' ',T9,
     *'12',T20,'X''+1',T30,'Y''',T40,'Z''-1'/' ',T9,'13',T20,'X''',T30,
     *'Y''+1',T40,'Z''-1'/' ',T9,'14',T20,'X''+1',T30,'Y''+1',T40,
     *'Z''-1'/' ',T9,'15',T20,'X''-1',T30,'Y''',T40,'Z''-1'/' ',T9,'16',
     *T20,'X''-1',T30,'Y''+1',T40,'Z''-1'/' ',T9,'17',T20,'X''',T30,
     *'Y''-1',T40,'Z''-1'/' ',T9,'18',T20,'X''+1',T30,'Y''-1',T40,
     *'Z''-1'/' ',T9,'19',T20,'X''-1',T30,'Y''-1',T40,'Z''-1'/' ',T9,
     *'20',T20,'X''+1',T30,'Y''',T40,'Z'''/' ',T9,'21',T20,'X''',T30,
     *'Y''+1',T40,'Z'''/' ',T9,'22',T20,'X''+1',T30,'Y''+1',T40,'Z'''
     */' ',T9,'23',T20,'X''',T30,'Y''',T40,'Z''+1'/' ',T9,'24',T20,
     *'X''+1',T30,'Y''',T40,'Z''+1'/' ',T9,'25',T20,'X''',T30,'Y''+1',
     *T40,'Z''+1'/' ',T9,'26',T20,'X''+1',T30,'Y''+1',T40,'Z''+1')
 3696 FORMAT(' ',T9,'27',T20,'X''-2',T30,'Y''',T40,'Z'''/' ',T9,'28',
     *T20,'X''+2',T30,'Y''',T40,'Z'''/' ',T9,'29',T20,'X''',T30,'Y''-2',
     *T40,'Z'''/' ',T9,'30',T20,'X''',T30,'Y''+2',T40,'Z'''/' ',T9,'31',
     *T20,'X''',T30,'Y''',T40,'Z''-2'/' ',T9,'32',T20,'X''',T30,'Y''',
     *T40,'Z''+2'/' ',T9,'33',T20,'X''-2',T30,'Y''-2',T40,'Z'''/' ',T9,
     *'34',T20,'X''-2',T30,'Y''-1',T40,'Z'''/' ',T9,'35',T20,'X''-2',
     *T30,'Y''+1',T40,'Z'''/' ',T9,'36',T20,'X''-2',T30,'Y''+2',T40,
     *'Z'''/' ',T9,'37',T20,'X''-2',T30,'Y''',T40,'Z''-2'/' ',T9,'38',
     *T20,'X''-2',T30,'Y''',T40,'Z''-1'/' ',T9,'39',T20,'X''-2',T30,
     *'Y''',T40,'Z''+1'/' ',T9,'40',T20,'X''-2',T30,'Y''',T40,'Z''+2'/
     *' ',T9,'41',T20,'X''-1',T30,'Y''-2',T40,'Z'''/' ',T9,'42',T20,
     *'X''-1',T30,'Y''+2',T40,'Z'''/' ',T9,'43',T20,'X''-1',T30,'Y''',
     *T40,'Z''-2'/' ',T9,'44',T20,'X''-1',T30,'Y''',T40,'Z''+2'/' ',T9,
     *'45',T20,'X''',T30,'Y''-2',T40,'Z''-2'/' ',T9,'46',T20,'X''',T30,
     *'Y''-2',T40,'Z''-1'/' ',T9,'47',T20,'X''',T30,'Y''-2',T40,'Z''+1'
     */' ',T9,'48',T20,'X''',T30,'Y''-2',T40,'Z''+2'/' ',T9,'49',T20,
     *'X''',T30,'Y''-1',T40,'Z''-2'/' ',T9,'50',T20,'X''',T30,'Y''-1',
     *T40,'Z''+2'/' ',T9,'51',T20,'X''',T30,'Y''+1',T40,'Z''-2'/' ',T9,
     *'52',T20,'X''',T30,'Y''+1',T40,'Z''+2'/' ',T9,'53',T20,'X''',T30,
     *'Y''+2',T40,'Z''-2'/' ',T9,'54',T20,'X''',T30,'Y''+2',T40,'Z''-1')
 3697 FORMAT(' ',T9,'55',T20,'X''',T30,'Y''+2',T40,'Z''+1'/' ',T9,'56',
     *T20,'X''',T30,'Y''+2',T40,'Z''+2'/' ',T9,'57',T20,'X''+1',T30,
     *'Y''-2',T40,'Z'''/' ',T9,'58',T20,'X''+1',T30,'Y''+2',T40,'Z'''/
     *' ',T9,'59',T20,'X''+1',T30,'Y''',T40,'Z''-2'/' ',T9,'60',T20,
     *'X''+1',T30,'Y''',T40,'Z''+2'/' ',T9,'61',T20,'X''+2',T30,'Y''-2',
     *T40,'Z'''/' ',T9,'62',T20,'X''+2',T30,'Y''-1',T40,'Z'''/' ',T9,
     *'63',T20,'X''+2',T30,'Y''+1',T40,'Z'''/' ',T9,'64',T20,'X''+2',
     *T30,'Y''+2',T40,'Z'''/' ',T9,'65',T20,'X''+2',T30,'Y''',T40,
     *'Z''-2'/' ',T9,'66',T20,'X''+2',T30,'Y''',T40,'Z''-1'/' ',T9,'67',
     *T20,'X''+2',T30,'Y''',T40,'Z''+1'/' ',T9,'68',T20,'X''+2',T30,
     *'Y''',T40,'Z''+2'/' ',T9,'69',T20,'X''-2',T30,'Y''-2',T40,
     *'Z''-2'/' ',T9,'70',T20,'X''-2',T30,'Y''-2',T40,'Z''-1'/' ',T9,
     *'71',T20,'X''-2',T30,'Y''-2',T40,'Z''+1'/' ',T9,'72',T20,'X''-2',
     *T30,'Y''-2',T40,'Z''+2'/' ',T9,'73',T20,'X''-2',T30,'Y''-1',T40,
     *'Z''-2'/' ',T9,'74',T20,'X''-2',T30,'Y''-1',T40,'Z''-1'/' ',T9,
     *'75',T20,'X''-2',T30,'Y''-1',T40,'Z''+1'/' ',T9,'76',T20,'X''-2',
     *T30,'Y''-1',T40,'Z''+2'/' ',T9,'77',T20,'X''-2',T30,'Y''+1',T40,
     *'Z''-2'/' ',T9,'78',T20,'X''-2',T30,'Y''+1',T40,'Z''-1'/' ',T9,
     *'79',T20,'X''-2',T30,'Y''+1',T40,'Z''+1'/' ',T9,'80',T20,'X''-2',
     *T30,'Y''+1',T40,'Z''+2')
 3698 FORMAT(' ',T9,'81',T20,'X''-2',T30,'Y''+2',T40,'Z''-2'/' ',T9,
     *'82',T20,'X''-2',T30,'Y''+2',T40,'Z''-1'/' ',T9,'83',T20,'X''-2',
     *T30,'Y''+2',T40,'Z''+1'/' ',T9,'84',T20,'X''-2',T30,'Y''+2',T40,
     *'Z''+2'/' ',T9,'85',T20,'X''-1',T30,'Y''-2',T40,'Z''-2'/' ',T9,
     *'86',T20,'X''-1',T30,'Y''-2',T40,'Z''-1'/' ',T9,'87',T20,'X''-1',
     *T30,'Y''-2',T40,'Z''+1'/' ',T9,'88',T20,'X''-1',T30,'Y''-2',T40,
     *'Z''+2'/' ',T9,'89',T20,'X''-1',T30,'Y''-1',T40,'Z''-2'/' ',T9,
     *'90',T20,'X''-1',T30,'Y''-1',T40,'Z''+2'/' ',T9,'91',T20,'X''-1',
     *T30,'Y''+1',T40,'Z''-2'/' ',T9,'92',T20,'X''-1',T30,'Y''+1',T40,
     *'Z''+2'/' ',T9,'93',T20,'X''-1',T30,'Y''+2',T40,'Z''-2'/' ',T9,
     *'94',T20,'X''-1',T30,'Y''+2',T40,'Z''-1'/' ',T9,'95',T20,'X''-1',
     *T30,'Y''+2',T40,'Z''+1'/' ',T9,'96',T20,'X''-1',T30,'Y''+2',
     *T40,'Z''+2'/' ',T9,'97',T20,'X''+1',T30,'Y''-2',T40,'Z''-2'/' ',
     *T9,'98',T20,'X''+1',T30,'Y''-2',T40,'Z''-1'/' ',T9,'99',T20,'X''+
     *1',T30,'Y''-2',T40,'Z''+1'/' ',T8,'100',T20,'X''+1',T30,'Y''-2',
     *T40,'Z''+2'/' ',T8,'101',T20,'X''+1',T30,'Y''-1',T40,'Z''-2'/' ',
     *T8,'102',T20,'X''+1',T30,'Y''-1',T40,'Z''+2'/' ',T8,'103',T20,
     *'X''+1',T30,'Y''+1',T40,'Z''-2')
 3699 FORMAT(' ',T8,'104',T20,'X''+1',T30,'Y''+1',T40,'Z''+2'/' ',T8,
     *'105',T20,'X''+1',T30,'Y''+2',T40,'Z''-2'/' ',T8,'106',T20,'X''+
     *1',T30,'Y''+2',T40,'Z''-1'/' ',T8,'107',T20,'X''+1',T30,'Y''+2',
     *T40,'Z''+1'/' ',T8,'108',T20,'X''+1',T30,'Y''+2',T40,'Z''+2'/' ',
     *T8,'109',T20,'X''+2',T30,'Y''-2',T40,'Z''-2'/' ',T8,'110',T20,
     *'X''+2',T30,'Y''-2',T40,'Z''-1'/' ',T8,'111',T20,'X''+2',T30,
     *'Y''-2',T40,'Z''+1'/' ',T8,'112',T20,'X''+2',T30,'Y''-2',T40,
     *'Z''+2'/' ',T8,'113',T20,'X''+2',T30,'Y''-1',T40,'Z''-2'/' ',T8,
     *'114',T20,'X''+2',T30,'Y''-1',T40,'Z''-1'/' ',T8,'115',T20,
     *'X''+2',T30,'Y''-1',T40,'Z''+1'/' ',T8,'116',T20,'X''+2',T30,
     *'Y''-1',T40,'Z''+2'/' ',T8,'117',T20,'X''+2',T30,'Y''+1',T40,
     *'Z''-2'/' ',T8,'118',T20,'X''+2',T30,'Y''+1',T40,'Z''-1'/' ',T8,
     *'119',T20,'X''+2',T30,'Y''+1',T40,'Z''+1'/' ',T8,'120',T20,'X''+2'
     *,T30,'Y''+1',T40,'Z''+2'/' ',T8,'121',T20,'X''+2',T30,'Y''+2',T40,
     *'Z''-2'/' ',T8,'122',T20,'X''+2',T30,'Y''+2',T40,'Z''-1'/' ',T8,
     *'123',T20,'X''+2',T30,'Y''+2',T40,'Z''+1'/' ',T8,'124',T20,'X''+2'
     *,T30,'Y''+2',T40,'Z''+2')
03701 IF(IPM)03801,04001,03801
C        COMPUTE PARAMETER INCREMENTS USED TO OBTAIN DERIVATIVES
03801    K=1
         L=NV
         DO 03903 I=1,NV
      DP(I)=(0.01)*SQRT(PM(K))
            K=K+L
03903       L=L-1
 4001 WRITE(3,'(1H1/1H ,A/1H )') TITLE
C        READ ONE SET OF INSTRUCTIONS
04101    K=1
         DO 04501 I=1,10
            L=K+23
      READ(9,'(A)',END=4701) EOR
      IF (EOR .NE. 'ENDOR') THEN
         BACKSPACE 9
      ELSE
         GO TO 4701
      END IF
**    READ(9,202) (IN(J),J=K,L)
      READ(9,'(24I5)') (IN(J),J=K,L)
            IF(IN(L))04501,04601,04501
04501       K=L
04601    IF(IN(1))04801,04701,04801
C           END OF JOB
*4701 CLOSE(UNIT=9,STATUS='DELETE')
*     FOR RIETAN-ORFFE
*4701 CLOSE(UNIT=9)
 4701 CONTINUE
*     RETURN
      GO TO 9999
04801    IF(IN(1)-INSAVE)04901,05001,04901
C           PUT OUT HEADING FOR NEW TYPE OF FUNCTION
04901       CALL HEDI(IN(1))
05001    INSAVE=IN(1)
         IF(IN(1)-100)05101,05101,05201
C           COMPUTE SINGLE-VALUED FUNCTION
05101       CALL SUB19
            GO TO 04101
C           COMPUTE MULTIPLE-VALUED FUNCTION
05201       CALL SUB21
            GO TO 04101
 9999 CONTINUE
      CALL DSORT
      STOP 
      END
      
C**********************************************************************C
C                            SUB8                                      C
C          COMPUTE ALL DISTANCES BETWEEN TWO ASYMMETRIC UNITS          C
C**********************************************************************C
      SUBROUTINE SUB8
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      COMMON /FLAG1/ LWR
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
C     SET MAXIMUM DISTANCE ACCEPTED
      IF(IN(6))107,107,111
  107    DMAX=4.0
         GO TO 113
  111 DMAX=FLOAT(IN(6))/10.0
  113 NA=IN(2)
      IF(NA)217,217,117
C     START LOOP FOR FIRST ATOM
  117 DO 215 I=1,NA
         IN(2)=I
C        START LOOP FOR SECOND ATOM
        DO 213 J=1,NA
C           TEST TO AVOID DUPLICATION
            IF(IN(3)-IN(5))201,125,201
  125          IF(J-I)201,215,215
  201       IN(4)=J
            NG=0
C           COMPUTE AND TEST DISTANCE
            CALL FUNI(IN(1))
            IF(FX-DMAX)209,209,213
C              COMPUTE ERROR AND PUT OUT RESULTS
  209          F=FX
               CALL SUB13
  213       CONTINUE
C           END LOOP FOR SECOND ATOM
  215    CONTINUE
C        END LOOP FOR FIRST ATOM
* 217 RETURN
  217 LWR=1
      END
C**********************************************************************C
C                              SUB10                                   C
C                   MULTIVALUED FUNCTIONS 7, 8, AND 9                  C
C**********************************************************************C
      SUBROUTINE SUB10
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
C     START LOOP THRU PRINCIPAL AXES
      DO 109 I=1,3
         IN(4)=I
C        COMPUTE AND PUT OUT FUNCTION AND ERROR
  109    CALL SUB19
      RETURN
      END
C**********************************************************************C
C                                SUB11                                 C
C                   MULTIVALUED FUNCTIONS 10 AND 11                    C
C**********************************************************************C
      SUBROUTINE SUB11
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
C     START LOOP THRU PRINCIPAL AXES
      DO 113 I=1,3
         IN(4)=I
C        START LOOP THRU REFERENCE AXES
         DO 113 J=1,3
            IN(5)=J
C           COMPUTE AND PUT OUT FUNCTION AND ERROR
  113       CALL SUB19
      RETURN
      END
C**********************************************************************C
C                              SUB13                                   C
C                    ERROR CALCULATION AND OUTPUT                      C
C**********************************************************************C
      SUBROUTINE SUB13
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      COMMON /FLAG1/ LWR
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
C     PUT OUT FUNCTION DESCRIPTION
      CALL OUTI(IN(1))
      IF(NG)107,113,107
C        PUT OUT ERROR INDICATOR IF NOT ZERO
  107 WRITE(3,109) NG
  109     FORMAT (1H ,51X,3H***,I3)
         GO TO 723
  113 VARA=0.0
      VARP=0.0
      IF(IAM)119,313,119
C        COMPUTE DERIVATIVES WITH RESPECT TO CELL PARAMETERS
  119    DO 211 I=1,6
            IF(DA(I))201,123,201
  123          DFDA(I)=0.0
               GO TO 211
  201          SAVEA=A(I)
               A(I)=A(I)+DA(I)
               CALL FUNI(IN(1))
               A(I)=SAVEA
               DFDA(I)=(FX-F)/DA(I)
  211       CONTINUE
C        COMPUTE VARIANCE BASED ON CELL PARAMETERS
         K=1
         L=6
         DO 311 I=1,6
            IF(DFDA(I))225,221,225
  221          K=K+L
               GO TO 311
  225          C=1.0
               DO 309 J=I,6
                  IF(DFDA(J))305,307,305
  305                VARA=VARA+C*DFDA(I)*DFDA(J)*AM(K)
  307             K=K+1
  309             C=2.0
  311       L=L-1
  313 IF(IPM)315,615,315
C        SELECT DERIVATIVES TO BE COMPUTED
  315    DO 319 I=1,NP
  319       KI2(I)=0
         CALL PREI(IN(1))
C        COMPUTE DERIVATIVES WITH RESPECT TO STRUCTURE PARAMETERS
         J=0
         LNZ=0
         DO 513 I=1,NV
  403       J=J+1
            IF(KI1(J))407,403,407
  407       IF(KI2(J))413,409,413
  409          DFDP(I)=0.0
               GO TO 513
  413          IF(DP(I))501,409,501
  501          SAVEP=P(J)
               P(J)=P(J)+DP(I)
               CALL FUNI(IN(1))
               P(J)=SAVEP
               DFDP(I)=(FX-F)/DP(I)
               LNZ=I
  513       CONTINUE
C        COMPUTE VARIANCE BASED ON STRUCTURE PARAMETERS
         KK=1
         KKD=NV
         DO 613 I=1,LNZ
            IF(DFDP(I))523,612,523
  523          K=KK
               C=1.0
               DO 611 J=I,LNZ
                  IF(DFDP(J))607,609,607
  607             VARP=VARP+C*DFDP(I)*DFDP(J)*PM(K)
  609             K=K+1
  611             C=2.0
  612       KK=KK+KKD
  613       KKD=KKD-1
  615 IF(NG)617,623,617
C        PUT OUT ERROR INDICATOR IF NOT ZERO
 617  WRITE(3,619) F,NG
**619    FORMAT (1H 48X,F9.4,6X,3H***I3)
  619    FORMAT (' ',52X,F9.4,6X,'***',I3)
         GO TO 723
C     COMPUTE STANDARD ERRORS AND PUT OUT RESULTS
  623 E1=SQRT(VARP)
      E=SQRT(VARP+VARA)
      IF(IPM)703,705,703
  703    IF(IAM)707,713,707
  705    IF(IAM)713,719,713
* 707 WRITE(3,709) F,E,E1
  707 IF (LWR .EQ. 1) WRITE(3,709) F,E,E1
**709       FORMAT (1H 48X,F9.4,5H +OR-F7.4,2H (F7.4,1H))
  709       FORMAT (' ',52X,F9.4,' +OR-',F7.4,' (',F7.4,')')
            GO TO 723
  713 WRITE(3,715) F,E
**715       FORMAT (1H 48X,F9.4,5H +OR-F7.4)
  715       FORMAT (' ',52X,F9.4,' +OR-',F7.4)
            GO TO 723
  719 WRITE(3,721) F
**721       FORMAT (1H 48X,F9.4)
  721       FORMAT (' ',52X,F9.4)
  723 RETURN
      END
C**********************************************************************C
C                                SUB19                                 C
C                   FUNCTION AND ERROR CALCULATION                     C
C**********************************************************************C
      SUBROUTINE SUB19
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      NG=0
      CALL FUNI(IN(1))
      F=FX
      CALL SUB13
      RETURN
      END
C**********************************************************************C
C                                SUB21                                 C
C                    COMPUTE MULTIVALUED FUNCTIONS                     C
C**********************************************************************C
      SUBROUTINE SUB21
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
C     UNPACK INSTRUCTION NUMBER
      IH=IN(1)/100
      IN(1)=IN(1)-100*IH
C     TRANSFER TO APPROPRIATE SECTION
      IF(IN(1)-11) 111,111,400
  111 IF(IN(1)-9)113,113,305
  113 IF(IN(1)-6)115,115,215
  115 IF(IN(1)-1)117,117,319
  117    IF(IH-1)119,119,123
C           COMPUTE ALL DISTANCES BETWEEN TWO UNITS
  119       CALL SUB8
            GO TO 319
C           COMPUTE DISTANCES INVOLVING BASIC UNIT AND ALL OTHERS
  123       NA=IN(2)
            NSP=NS+1
            IN(3)=0
**          DO 211 I = 1, 8
            DO 211 I = 1, 29
               DO 211 J=1,NSP
                  IN(2)=NA
**                IN(5)=100*(I-1)+J-1
                  IN(5) = 1000*(I-1) + J - 1
  211             CALL SUB8
            GO TO 319
  215    IF(IH-1)217,217,221
C           COMPUTE FUNCTIONS INVOLVING THREE PRINCIPAL AXES
  217       CALL SUB10
            GO TO 319
C           COMPUTE FUNCTIONS INVOLVING THREE PRINCIPAL AXES
C                                            FOR ALL ATOMS
  221       NA=IN(2)
            DO 301 I=1,NA
               IN(2)=I
  301          CALL SUB10
            GO TO 319
  305    IF(IH-1)307,307,311
C           COMPUTE FOR ALL PRINCIPAL AXES AND ALL REFERENCE AXES
  307       CALL SUB11
            GO TO 319
C           COMPUTE FOR ALL PRINCIPAL AXES, ALL REFERENCE AXES,
C                                              AND ALL ATOMS
  311       NA=IN(2)
            DO 317 I=1,NA
               IN(2)=I
  317          CALL SUB11
  319 RETURN
  400 IF(IN(1)-16) 319,410,319
  410 NA=IN(2)
      DO 411 I=1,NA
      IN(2)=I
  411 CALL SUB19
      GO TO 319
      END
C**********************************************************************C
C                              SETKX                                   C
C                SET KEY WORDS FOR ATOM COORDINATES                    C
C**********************************************************************C
      SUBROUTINE SETKX(I)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     I=IN(K), THE INSTRUCTION INTEGER SPECIFYING THE ATOM NUMBER
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      IF(I)119,119,111
  111 J=JX+JXP*(I-1)
      KI2(J)=1
      KI2(J+1)=1
      KI2(J+2)=1
  119 RETURN
      END
C**********************************************************************C
C                                 SETKB                                C
C                     SET KEY WORDS FOR ATOM BETAS                     C
C**********************************************************************C
      SUBROUTINE SETKB(I)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     I=IN(K), THE INSTRUCTION INTEGER SPECIFYING THE ATOM NUMBER
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      IF(I)119,119,111
  111 J=JB+JBP*(I-1)
      DO 117 K=1,6
      KI2(J)=1
  117 J=J+1
  119 RETURN
      END
C**********************************************************************C
C                               STAOAA                                 C
C                         STORE METRIC TENSOR                          C
C**********************************************************************C
      SUBROUTINE STOAA
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      AA(1,1)=A(1)*A(1)
      AA(2,2)=A(2)*A(2)
      AA(3,3)=A(3)*A(3)
      AA(1,2)=A(1)*A(2)*A(6)
      AA(2,1)=AA(1,2)
      AA(1,3)=A(1)*A(3)*A(5)
      AA(3,1)=AA(1,3)
      AA(2,3)=A(2)*A(3)*A(4)
      AA(3,2)=AA(2,3)
      RETURN
      END
C**********************************************************************C
C                                STOBB                                 C
C                    STORE RECIPROCAL METRIC TENSOR                    C
C**********************************************************************C
      SUBROUTINE STOBB
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      DIMENSION AI(6),CI(6),BBII(3),BBJK(3)
      DO 121 I=1,3
      IF(A(I))111,111,115
  111 NG=6
      GO TO 225
  115 AI(I)=A(I)
      AI(I+3)=A(I)
      CI(I)=A(I+3)
  121 CI(I+3)=A(I+3)
      X=1.0/(AI(1)*AI(2)*AI(3)*(1.0-CI(1)*CI(1)-CI(2)*CI(2)
     1-CI(3)*CI(3)+2.0*CI(1)*CI(2)*CI(3)))
      DO 205 I=1,3
      BBII(I)=X*(1.0-CI(I)*CI(I))*AI(I+1)*AI(I+2)/AI(I)
  205 BBJK(I)=X*AI(I)*(CI(I+1)*CI(I+2)-CI(I))
      BB(1,1)=BBII(1)
      BB(1,2)=BBJK(3)
      BB(1,3)=BBJK(2)
      BB(2,1)=BBJK(3)
      BB(2,2)=BBII(2)
      BB(2,3)=BBJK(1)
      BB(3,1)=BBJK(2)
      BB(3,2)=BBJK(1)
      BB(3,3)=BBII(3)
  225 RETURN
      END
C**********************************************************************C
C                                ATOM                                  C
C                     ATOM COORDINATE SUBROUTINE                       C
C**********************************************************************C
      SUBROUTINE ATOM(I1,Z)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      DIMENSION I(2),X(3),Y(3),Z(3)
      I(1)=IN(I1)
      I(2)=IN(I1+1)
      IF(I(1))109,109,117
  109 X(1)=0.0
      X(2)=0.0
      X(3)=0.0
      GO TO 125
  117 K=JXP*(I(1)-1)+JX
      IF(K+2-NP)119,119,503
  503 NG=5
      GO TO 325
  119 DO 123 J=1,3
      X(J)=P(K)
  123 K=K+1
**125 KC=I(2)/100
**    KS=I(2)-100*KC
  125 KC = I(2)/1000
      KS = I(2) - 1000*KC
      IF(KS-NS)203,203,403
  403 NG=1
      GO TO 325
  203 IF(KS)403,205,213
  205 Y(1)=X(1)
      Y(2)=X(2)
      Y(3)=X(3)
      GO TO 311
  213 DO 215 J=1,3
  215 Y(J)=TS(J,KS)
      DO 309 K=1,3
      DO 307 J=1,2
      L=IS(J,K,KS)
      IF(L)225,307,305
  225 L=-L
      Y(K)=Y(K)-X(L)
      GO TO 307
  305 Y(K)=Y(K)+X(L)
  307 CONTINUE
  309 CONTINUE
**311 KC4=KC/4
**    KC3=KC-4*KC4
**    KC2=KC3/2
**    KC1=KC3-2*KC2
**    Z(1)=Y(1)-FLOAT (KC1)
**    Z(2)=Y(2)-FLOAT (KC2)
**    Z(3)=Y(3)-FLOAT (KC4)

  311 Z(1) = Y(1)
      Z(2) = Y(2)
      Z(3) = Y(3)
  
      SELECT CASE (KC)
         CASE (0)
            RETURN
         CASE (1)
            Z(1) = Y(1) - 1.0
         CASE (2)
            Z(1) = Y(1) - 1.0
	    Z(2) = Y(2) + 1.0
         CASE (3)
            Z(1) = Y(1) - 1.0
            Z(3) = Y(3) + 1.0
         CASE (4)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) + 1.0
         CASE (5)
            Z(2) = Y(2) - 1.0
         CASE (6)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) - 1.0
         CASE (7)
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) + 1.0
         CASE (8)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) + 1.0
         CASE (9)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) - 1.0
         CASE (10)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) + 1.0
         CASE (11)
            Z(3) = Y(3) - 1.0
         CASE (12)
            Z(1) = Y(1) + 1.0
            Z(3) = Y(3) - 1.0
         CASE (13)
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) - 1.0
         CASE (14)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) - 1.0
         CASE (15)
            Z(1) = Y(1) - 1.0
            Z(3) = Y(3) - 1.0
         CASE (16)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) - 1.0
         CASE (17)
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) - 1.0
         CASE (18)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) - 1.0
         CASE (19)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) - 1.0
         CASE (20)
            Z(1) = Y(1) + 1.0
         CASE (21)
            Z(2) = Y(2) + 1.0
         CASE (22)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) + 1.0
         CASE (23)
            Z(3) = Y(3) + 1.0
         CASE (24)
            Z(1) = Y(1) + 1.0
            Z(3) = Y(3) + 1.0
         CASE (25)
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) + 1.0
         CASE (26)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) + 1.0
         CASE (27)
            Z(1) = Y(1) - 2.0
         CASE (28)
            Z(1) = Y(1) + 2.0
         CASE (29)
            Z(2) = Y(2) - 2.0
         CASE (30)
            Z(2) = Y(2) + 2.0
         CASE (31)
            Z(3) = Y(3) - 2.0
         CASE (32)
            Z(3) = Y(3) + 2.0
         CASE (33)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) - 2.0
         CASE (34)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) - 1.0
         CASE (35)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) + 1.0
         CASE (36)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) + 2.0
         CASE (37)
            Z(1) = Y(1) - 2.0
            Z(3) = Y(3) - 2.0
         CASE (38)
            Z(1) = Y(1) - 2.0
            Z(3) = Y(3) - 2.0
         CASE (39)
            Z(1) = Y(1) - 2.0
            Z(3) = Y(3) + 1.0
         CASE (40)
            Z(1) = Y(1) - 2.0
            Z(3) = Y(3) + 2.0
         CASE (41)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) - 2.0
         CASE (42)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) + 2.0
         CASE (43)
            Z(1) = Y(1) - 1.0
            Z(3) = Y(3) - 2.0
         CASE (44)
            Z(1) = Y(1) - 1.0
            Z(3) = Y(3) + 2.0
         CASE (45)
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) - 2.0
         CASE (46)
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) - 1.0
         CASE (47)
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) + 1.0
         CASE (48)
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) + 2.0
         CASE (49)
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) - 2.0
         CASE (50)
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) + 2.0
         CASE (51)
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) - 2.0
         CASE (52)
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) + 2.0
         CASE (53)
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) - 2.0
         CASE (54)
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) - 1.0
         CASE (55)
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) + 1.0
         CASE (56)
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) + 2.0
         CASE (57)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) - 2.0
         CASE (58)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) + 2.0
         CASE (59)
            Z(1) = Y(1) + 1.0
            Z(3) = Y(3) - 2.0
         CASE (60)
            Z(1) = Y(1) + 1.0
            Z(3) = Y(3) + 2.0
         CASE (61)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) - 2.0
         CASE (62)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) - 1.0
         CASE (63)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) + 1.0
         CASE (64)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) + 2.0
         CASE (65)
            Z(1) = Y(1) + 2.0
            Z(3) = Y(3) - 2.0
         CASE (66)
            Z(1) = Y(1) + 2.0
            Z(3) = Y(3) - 1.0
         CASE (67)
            Z(1) = Y(1) + 2.0
            Z(3) = Y(3) + 1.0
         CASE (68)
            Z(1) = Y(1) + 2.0
            Z(3) = Y(3) + 2.0
         CASE (69)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) - 2.0
         CASE (70)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) - 1.0
         CASE (71)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) + 1.0
         CASE (72)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) + 2.0
         CASE (73)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) - 2.0
         CASE (74)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) - 1.0
         CASE (75)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) + 1.0
         CASE (76)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) + 2.0
         CASE (77)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) - 2.0
         CASE (78)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) - 1.0
         CASE (79)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) + 1.0
         CASE (80)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) + 2.0
         CASE (81)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) - 2.0
         CASE (82)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) - 1.0
         CASE (83)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) + 1.0
         CASE (84)
            Z(1) = Y(1) - 2.0
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) + 2.0
         CASE (85)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) - 2.0
         CASE (86)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) - 1.0
         CASE (87)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) + 1.0
         CASE (88)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) + 2.0
         CASE (89)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) - 2.0
         CASE (90)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) + 2.0
         CASE (91)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) - 2.0
         CASE (92)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) + 2.0
         CASE (93)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) - 2.0
         CASE (94)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) - 1.0
         CASE (95)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) + 1.0
         CASE (96)
            Z(1) = Y(1) - 1.0
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) + 2.0
         CASE (97)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) - 2.0
         CASE (98)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) - 1.0
         CASE (99)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) + 1.0
         CASE (100)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) + 2.0
         CASE (101)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) - 2.0
         CASE (102)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) + 2.0
         CASE (103)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) - 2.0
         CASE (104)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) + 2.0
         CASE (105)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) - 2.0
         CASE (106)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) - 1.0
         CASE (107)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) + 1.0
         CASE (108)
            Z(1) = Y(1) + 1.0
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) + 2.0
         CASE (109)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) - 2.0
         CASE (110)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) - 1.0
         CASE (111)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) + 1.0
         CASE (112)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) - 2.0
            Z(3) = Y(3) + 2.0
         CASE (113)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) - 2.0
         CASE (114)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) - 1.0
         CASE (115)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) + 1.0
         CASE (116)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) - 1.0
            Z(3) = Y(3) + 2.0
         CASE (117)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) - 2.0
         CASE (118)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) - 1.0
         CASE (119)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) + 1.0
         CASE (120)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) + 1.0
            Z(3) = Y(3) + 2.0
         CASE (121)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) - 2.0
         CASE (122)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) - 1.0
         CASE (123)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) + 1.0
         CASE (124)
            Z(1) = Y(1) + 2.0
            Z(2) = Y(2) + 2.0
            Z(3) = Y(3) + 2.0
      END SELECT

  325 RETURN
      END
      
C**********************************************************************C
C                                BETA                                  C
C         STORE TRANSFORMED ANISOTROPIC TEMP FACTOR MATRIX             C
C         INS IS ATOM DESCRIPTION, Z IS TRANSFORMED MATRIX             C
C**********************************************************************C
      SUBROUTINE BETA(INS1,Z)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      DIMENSION INS(2),Z(3,3),B1(6),B2(9)
      INS(1)=IN(INS1)
      INS(2)=IN(INS1+1)
      IF(ITF-1)100,111,115
  100 IF(JBP-6)111,115,115
  111 NG=4
      GO TO 423
  115 KS= MOD (INS(2),100)
      IF(KS)121,119,119
  119 IF(KS-NS)125,125,121
  121 NG=1
      GO TO 423
  125 IF(INS(1))211,201,207
  201 DO 203 I=1,6
  203 B1(I)=0.0
      GO TO 221
  207 J=JB+JBP*(INS(1)-1)
      IF(J+5-NP)215,215,211
  211 NG=5
      GO TO 423
  215 DO 219 I=1,6
      B1(I)=P(J)
  219 J=J+1
  221 B2(1)=B1(1)
      B2(2)=B1(4)
      B2(3)=B1(5)
      B2(4)=B1(2)
      B2(6)=B1(6)
      B2(9)=B1(3)
      DO 421 I=1,3
      DO 419 J=I,3
      IF(KS)121,313,319
  313 M=I*J
      B3=B2(M)
      GO TO 415
  319 B3=0.0
      DO 413 K=1,2
      DO 411 L=1,2
      M=IS(K,I,KS)*IS(L,J,KS)
      IF(M)407,411,403
  403 B3=B3+B2(M)
      GO TO 411
  407 M=-M
      B3=B3-B2(M)
  411 CONTINUE
  413 CONTINUE
  415 Z(I,J)=B3
      Z(J,I)=B3
  419 CONTINUE
  421 CONTINUE
  423 RETURN
      END
C**********************************************************************C
C                              HEDI                                    C
C              SELECT THE HED SUBROUTINE TO BE ENTERED                 C
C**********************************************************************C
      SUBROUTINE HEDI(I)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     HEADING 1
  101 FORMAT(34H0INTERATOMIC DISTANCE IN ANGSTROMS)
C     HEADING 2
  102 FORMAT(//47H0BOND ANGLE IN DEGREES.  CENTRAL ATOM IS VERTEX)
C     HEADING 3
  103 FORMAT(58H0DIHEDRAL ANGLE BETWEEN PLANES EACH DEFINED BY THREE ATO
     1MS)
C     HEADING 4
  104 FORMAT(45H0DIFFERENCE BETWEEN TWO INTERATOMIC DISTANCES)
C     HEADING 5
  105 FORMAT(35H0DIFFERENCE BETWEEN TWO BOND ANGLES)
C     HEADING 6
  106 FORMAT(27H0SUM OF SEVERAL BOND ANGLES)
C     HEADING 7
  107 FORMAT(72H0RMS COMPONENT OF THERMAL DISPLACEMENT ALONG PRINCIPAL A
     1XIS R. ANGSTROMS/25H0            ATOM       R)
C     HEADING 8
  108 FORMAT(63H0ANGLE BETWEEN PRINCIPAL AXIS R AND VECTOR DEFINED BY TW
     1O ATOMS/40H0            ATOM       R         VECTOR)
C     HEADING 9
  109 FORMAT('0RMS COMPONENT OF THERMAL DISPLACEMENT ALONG PRINCIPAL'//
     1'AXIS R PROJECTED ON VECTOR DEFINED BY TWO ATOMS. ANGSTROMS'/
     2'0H0         ATOM       R         VECTOR')
C     HEADING 10
 1100 FORMAT(85H0ANGLE BETWEEN PRINCIPAL AXIS R AND AXIS I OF CARTESIAN
     1SYSTEM DEFINED BY TWO VECTORS/45H0            ATOM     R  I   DEFI
     2NING VECTORS)
C     HEADING 11
  111 FORMAT(101H0RMS COMPONENT OF THERMAL DISPLACEMENT ALONG PRINCIPAL
     1AXIS R PROJECTED ON AXIS I OF CARTESIAN SYSTEM/34H DEFINED BY TWO
     2VECTORS. ANGSTROMS/45H0            ATOM     R  I   DEFINING VECTOR
     3S)
C     HEADING 12
  112 FORMAT(83H0RMS COMPONENT OF THERMAL DISPLACEMENT IN DIRECTION DEFI
     1NED BY TWO ATOMS. ANGSTROMS/40H0            ATOM                 V
     2ECTOR)
C     HEADING 13
  113 FORMAT(51H0RMS RADIAL THERMAL DISPLACEMENT OF ATOM. ANGSTROMS)
C     HEADING 14
  114 FORMAT(88H0INTERATOMIC DISTANCE AVERAGED OVER THERMAL MOTION. SECO
     1ND ATOM ASSUMED TO RIDE ON FIRST)
C     HEADING 15
  115 FORMAT(87H0INTERATOMIC DISTANCE AVERAGED OVER THERMAL MOTION. ATOM
     1S ASSUMED TO MOVE INDEPENDENTLY)
C     HEADING (16)
  116 FORMAT('0MEAN ISOTROPIC THERMAL PARAMETER')
      K= MOD (I,100)
      IF(K)160,160,4
    4 IF(K-16) 8,8,160
    8 GO TO (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,260),K
   10 WRITE (3,101)
         GO TO 160
   20 WRITE (3,102)
         GO TO 160
   30 WRITE (3,103)
         GO TO 160
   40 WRITE (3,104)
         GO TO 160
   50 WRITE (3,105)
         GO TO 160
   60 WRITE (3,106)
         GO TO 160
   70 WRITE (3,107)
         GO TO 160
   80 WRITE (3,108)
         GO TO 160
   90 WRITE (3,109)
         GO TO 160
  100 WRITE (3,1100)
         GO TO 160
  110 WRITE (3,111)
         GO TO 160
  120 WRITE (3,112)
         GO TO 160
  130 WRITE (3,113)
         GO TO 160
  140 WRITE (3,114)
         GO TO 160
  150 WRITE (3,115)
  160 RETURN
  260 WRITE(3,116)
      GO TO 160
      END
C**********************************************************************C
C                               PREI                                   C
C              SELECT THE PRE SUBROUTINE TO BE ENTERED                 C
C**********************************************************************C
      SUBROUTINE PREI(IKY)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
    4 IF(IKY)160,160,6
    6 IF(IKY-16) 8,8,160
    8 GO TO (10,20,30,40,30,60,70,80,80,100,100,120,70,140,140,70),IKY
   10 CONTINUE
C     PRELIMINARY SUBROUTINE 1              (PRE1)
C     SET KEY WORDS FOR INTERATOMIC DISTANCE
      CALL SETKX(IN(2))
      CALL SETKX(IN(4))
         GO TO 160
   20 CONTINUE
C     PRELIMINARY SUBROUTINE 2              (PRE2)
      DO 127 I=2,6,2
  127 CALL SETKX(IN(I))
         GO TO 160
   30 CONTINUE
C     PRELIMINARY SUBROUTINE 3              (PRE3)
      DO 137 I=2,12,2
  137 CALL SETKX(IN(I))
         GO TO 160
   40 CONTINUE
C     PRELIMINARY SUBROUTINE 4              (PRE4)
      DO 147 I=2,8,2
  147 CALL SETKX(IN(I))
         GO TO 160
   60 CONTINUE
C     PRELIMINARY SUBROUTINE 6              (PRE6)
      J=IN(2)*6+1
      DO 169 I=3,J,2
  169 CALL SETKX(IN(I))
         GO TO 160
   70 CONTINUE
C     PRELIMINARY SUBROUTINE 7              (PRE7)
      CALL SETKB(IN(2))
         GO TO 160
   80 CONTINUE
C     PRELIMINARY SUBROUTINE 8              (PRE8)
      CALL SETKB(IN(2))
      DO 189 I=5,7,2
  189 CALL SETKX(IN(I))
         GO TO 160
  100 CONTINUE
C     PRELIMINARY SUBROUTINE 10             (PRE10)
      CALL SETKB(IN(2))
      DO 109 I=6,12,2
  109 CALL SETKX(IN(I))
         GO TO 160
  120 CONTINUE
C     PRELIMINARY SUBROUTINE 12             (PRE12)
      CALL SETKB(IN(2))
      DO 129 I=4,6,2
  129 CALL SETKX(IN(I))
         GO TO 160
  140 CONTINUE
C     PRELIMINARY SUBROUTINE 14             (PRE14)
      DO 149 I=2,4,2
      CALL SETKX(IN(I))
  149 CALL SETKB(IN(I))
  160 RETURN
      END
C**********************************************************************C
C                                FUNI                                  C
C              SELECT THE FUN SUBROUTINE TO BE ENTERED                 C
C**********************************************************************C
      SUBROUTINE FUNI(IKY)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      DIMENSION V1(3),V2(3),V3(3),V4(3),V5(3),V6(3)
      DIMENSION X1(3),X2(3),X3(3),X4(3),X5(3),X6(3)
      DIMENSION W(3,3)
      DIMENSION CC(2)
      CALL SETA(A)
      CALL SETP(P)
      IF(IKY)6,6,5
    5 IF(IKY-16) 8,8,6
    6    NG=11
         GO TO 160
    8 GO TO (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,260),IKY
   10 CONTINUE
C     COMPUTE INTERATOMIC DISTANCE                    (FUN1)
      FX=FUND(2)
         GO TO 160
   20 CONTINUE
C     BOND ANGLE SUBROUTINE                           (FUN2)
      FX=FUNA(2)
         GO TO 160
   30 CONTINUE
C     DIHEDRAL ANGLE SUBROUTINE                       (FUN3)
      CALL STOAA
      CALL STOBB
      IF(NG)207,113,207
  113 CALL ATOM ( 2,X1)
      CALL ATOM ( 4,X2)
      CALL ATOM ( 6,X3)
      CALL ATOM ( 8,X4)
      CALL ATOM (10,X5)
      CALL ATOM (12,X6)
      IF(NG)207,119,207
  119 CALL DIFV (X2,X1,V1)
      CALL DIFV (X3,X1,V2)
      CALL DIFV (X5,X4,V3)
      CALL DIFV (X6,X4,V4)
      CALL NORM(V1,V2,V5)
      CALL NORM(V3,V4,V6)
      FX=ACOS(COSVV(V5,V6))*57.29577951
  207    GO TO 160
   40 CONTINUE
C     DIFFERENCE BETWEEN BOND DISTANCES               (FUN4)
      FX=FUND(2)-FUND(6)
         GO TO 160
   50 CONTINUE
C     DIFFERENCE BETWEEN BOND ANGLES                  (FUN5)
      FX=FUNA(8)-FUNA(2)
         GO TO 160
   60 CONTINUE
C     SUM OF BOND ANGLES                              (FUN6)
      N=IN(2)
      FX=0.0
      DO 111 J=1,N
  111 FX=FX+FUNA(6*J-3)
         GO TO 160
   70 CONTINUE
C     RMS PRINCIPAL DISPLACEMENT                      (FUN7)
      CALL FUNB(W,Z,FX)
         GO TO 160
   80 CONTINUE
C     ANGLE BETWEEN PRINCIPAL AXIS AND VECTOR         (FUN8)
      CALL FUNC(C,Z)
      FX=ACOS(C)*57.29577951
         GO TO 160
   90 CONTINUE
C     PRINCIPAL AXIS PROJECTED ON VECTOR              (FUN9)
      CALL FUNC(C,Z)
      FX=C*Z
         GO TO 160
  100 CONTINUE
C     ANGLE BETWEEN PRINCIPAL AND CARTESIAN AXES      (FUN10)
      CALL FUNX(C,Z)
      FX=ACOS(C)*57.29577951
         GO TO 160
  110 CONTINUE
C     PRINCIPAL AXIS PROJECTED ON CARTESIAN AXIS      (FUN11)
      CALL FUNX(C,Z)
      FX=C*Z
         GO TO 160
  120 CONTINUE
C     RMS DISPLACEMENT IN GIVEN DIRECTION             (FUN12)
      CALL FUNXI(2,4,XISQ)
      FX=SQRT (XISQ)
         GO TO 160
  130 CONTINUE
C     RMS RADIAL DISPLACEMENT                         (FUN13)
      CALL FUNR(2,RSQ)
      FX=SQRT (RSQ)
         GO TO 160
  140 CONTINUE
C     MEAN BOND DISTANCE ASSUMING RIDING              (FUN14)
      CALL FUNCR(CC,R)
      FX=R+(CC(2)-CC(1))/(2.0*R)
         GO TO 160
  150 CONTINUE
C     MEAN INTERATOMIC DISTANCE ASSUMING INDEPENDENT MOTION  (FUN15)
      CALL FUNCR(CC,R)
      FX=R+(CC(2)+CC(1))/(2.0*R)
  160 RETURN
C     ISOTROPIC TEMPERATURE FACTOR FROM ANISOTROPIC ONE
  260 CONTINUE
      CALL FUNBE
      GO TO 160
      END
C**********************************************************************C
C                              OUTI                                    C
C            SELECT THE OUT SUBROUTINE TO BE ENTERED                   C
C**********************************************************************C
      SUBROUTINE OUTI(IKY)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
    4 IF(IKY)160,160,6
    6 IF(IKY-16) 8,8,160
    8 GO TO (10,20,30,40,30,60,70,80,80,100,100,120,130,10,10,260),IKY
   10 CONTINUE
C     PUT OUT DESCRIPTION OF INTERATOMIC DISTANCE
C     OUTPUT DESCRIPTION 1                  (OUT1)
*       WRITE            (3,117)(IN(I),I=2,5)
* 117 FORMAT(12H0          (I2,1H,I3,3H) (I2,1H,I3,1H))
      CALL FRCO(IN,TS,IS,P)
         GO TO 160
   20 CONTINUE
C     OUTPUT DESCRIPTION 2                  (OUT2)
          WRITE            (3, 127) (IN(I), I=2,7)
**127 FORMAT(12H0          (I2,1H,I3,3H) (I2,1H,I3,3H) (I2,1H,I3,1H))
C 127 FORMAT(12H0          (I2,1H,I5,3H) (I2,1H,I5,3H) (I2,1H,I5,1H))
  127 FORMAT('0',10X,'(',I2,',',I5,') (',I2,',',I5,') (',I2,',',I5,')')
         GO TO 160
   30 CONTINUE
C     OUTPUT DESCRIPTION 3                  (OUT3)
          WRITE            (3, 137)(IN(I),I=2,13)
**137 FORMAT(12H0          (I2,1H,I3,3H) (I2,1H,I3,3H) (I2,1H,I3,1H)/12H
**   1           (I2,1H,I3,3H) (I2,1H,I3,3H) (I2,1H,I3,1H))
C 137 FORMAT(12H0          (I2,1H,I5,3H) (I2,1H,I5,3H) (I2,1H,I5,1H)/12H
C    1           (I2,1H,I5,3H) (I2,1H,I5,3H) (I2,1H,I5,1H))
  137 FORMAT('0',10X,'(',I2,',',I5,') (',I2,',',I5,') (',I2,',',I5,')'/
     1      11X,'(',I2,',',I5,') (',I2,',',I5,') (',I2,',',I5,')')
         GO TO 160
   40 CONTINUE
C     OUTPUT DESCRIPTION 4                  (OUT4)
        WRITE            (3,147)(IN(I),I=2,9)
**147 FORMAT(12H0          (I2,1H,I3,3H) (I2,1H,I3,1H)/12H           (I2
**   1,1H,I3,3H) (I2,1H,I5,1H))
C 147 FORMAT(12H0          (I2,1H,I5,3H) (I2,1H,I5,1H)/12H           (I2
C    1,1H,I5,3H) (I2,1H,I5,1H))
  147 FORMAT('0',10X,'(',I2,',',I5,') (',I2,',',I5,')'/
     1       11X,'(',I2,',',I5,') (',I2,',',I5,')')
         GO TO 160
   60 CONTINUE
C     OUTPUT DESCRIPTION 6                  (OUT6)
      J=IN(2)*6+2
          WRITE            (3, 169) (IN(I), I=3,J)
**169 FORMAT(12H0          (I2,1H,I3,3H) (I2,1H,I3,3H) (I2,1H,I3,1H)/(12
**   1H           (I2,1H,I3,3H) (I2,1H,I3,3H) (I2,1H,I3,1H)))
C 169 FORMAT(12H0          (I2,1H,I5,3H) (I2,1H,I5,3H) (I2,1H,I5,1H)/(12
C    1H           (I2,1H,I5,3H) (I2,1H,I5,3H) (I2,1H,I5,1H)))
  169 FORMAT('0',10X,'(',I2,',',I5,') (',I2,',',I5,') (',I2,',',I5,')'/
     1       (11X,'(',I2,',',I5,') (',I2,',',I5,') (',I2,',',I5,')'))
         GO TO 160
   70 CONTINUE
C     OUTPUT DESCRIPTION 7                  (OUT7)
          WRITE            (3, 177) (IN(I), I=2,4)
**177 FORMAT(12H0          (I2,1H,I3,6H)     I1)
C 177 FORMAT(12H0          (I2,1H,I5,6H)     I1)
  177 FORMAT('0',10X,'(',I2,',',I5,')',5X,I1)
         GO TO 160
   80 CONTINUE
C     OUTPUT DESCRIPTION 8                  (OUT8)
          WRITE            (3, 187) (IN(I), I=2,8)
**187 FORMAT(12H0          (I2,1H,I3,6H)     I1,5H    (I2,1H,I3,3H) (I2,
**   11H,I3,1H))
C 187 FORMAT(12H0          (I2,1H,I5,6H)     I1,5H    (I2,1H,I5,3H) (I2,
C    11H,I5,1H))
  187 FORMAT('0',10X,'(',I2,',',I5,')',5X,I1,4X,'(',I2,',',I5,') (',
     1    I2,',',I5,')')
         GO TO 160
  100 CONTINUE
C     OUTPUT DESCRIPTION 10                 (OUT10)
          WRITE            (3,1007) (IN(I), I=2,13)
**1007 FORMAT(12H0          (I2,1H,I3,4H)   I1,I3,4H   (I2,1H,I3,3H) (I2,
**   11H,I3,1H)/30H                             (I2,1H,I3,3H) (I2,1H,I3,
**   21H))
C1007 FORMAT(12H0          (I2,1H,I5,4H)   I1,I5,4H   (I2,1H,I5,3H) (I2,
C    11H,I5,1H)/30H                             (I2,1H,I5,3H) (I2,1H,I5,
C    21H))
 1007 FORMAT('0',10X,'(',I2,',',I5,')',3X,I1,I5,3X,'(',I2,',',I5,') (',
     1       I2,',',I5,')'/29X,'(',I2,',',I5,') (',I2,',',I5,')')
         GO TO 160
  120 CONTINUE
C     OUTPUT DESCRIPTION 12                 (OUT12)
          WRITE            (3,1207) (IN(I), I=2,7)
**1207 FORMAT(12H0          (I2,1H,I3,12H)          (I2,1H,I3,3H) (I2,1H,
**   1I3,1H))
C1207 FORMAT(12H0          (I2,1H,I5,12H)          (I2,1H,I5,3H) (I2,1H,
C    1I5,1H))
 1207 FORMAT('0',10X,'(',I2,',',I5,')',10X,'(',I2,',',I5,') (',I2,',',
     1        I5,')')
         GO TO 160
  130 CONTINUE
C     OUTPUT DESCRIPTION 13                 (OUT13)
          WRITE            (3,1307)(IN(I), I=2,3)
**1307 FORMAT(12H0          (I2,1H,I3,1H))
C1307 FORMAT(12H0          (I2,1H,I5,1H))
 1307 FORMAT('0',10X,'(',I2,',',I5,')')
  160 RETURN
C     OUTPUT DESCRIPTION 16
  260 CONTINUE
      WRITE(3,2260) (IN(I),I=2,3)
**2260 FORMAT(12H0          (I2,1H,I3,1H))
C2260 FORMAT(12H0          (I2,1H,I5,1H))
 2260 FORMAT('0',10X,'(',I2,',',I5,')')
      GO TO 160
      END

************************************************************************

      SUBROUTINE FRCO(IN,TS,IS,P)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*     PRINT OUT THE FRACTIONAL COORDINATES OF TWO ATOMS
      PARAMETER(NNP=2048,NATM=10000,NW=384)
      INTEGER IN(*),IS(2,3,NW),NOATM(2,NATM)
      DOUBLE PRECISION TS(3,NW),P(*),XYZ(6),COORD(6,NATM)
      CHARACTER ATOMN((NNP-14)/10)*9
      COMMON /NEWC/ ATOMN
      COMMON /FLAG1/ LWR
      DATA IATM,TOL/1,0.0001D0/
      SAVE NOATM,COORD,IATM

**    CALL CALCD(P((IN(2)-1)*10+2), IN(3)/100, MOD(IN(3),100),
**   &  TS, IS, XYZ(1))
      CALL CALCD(P((IN(2)-1)*10+2), IN(3)/1000, MOD(IN(3),1000),
     &  TS, IS, XYZ(1))
     
**    CALL CALCD(P((IN(4)-1)*10+2), IN(5)/100, MOD(IN(5),100),
**   &  TS, IS, XYZ(4))
      CALL CALCD(P((IN(4)-1)*10+2), IN(5)/1000, MOD(IN(5),1000),
     &  TS, IS, XYZ(4))

      LWR=0
      IF (IN(2) .EQ. IN(4)) THEN
*        SKIP IF THE TWO ATOMS OVERLAP WITH EACH OTHER
         DO 10 J=1,3
            IF (ABS(XYZ(J)-XYZ(J+3)) .GT. TOL) GO TO 1
   10    CONTINUE
         RETURN
      END IF

*     SAVE ATOM # AND FRACTIONAL COORDINATES OF THE TWO ATOMS
    1 NOATM(1,IATM)=IN(2)
      NOATM(2,IATM)=IN(4)
      DO 15 J=1,6
         COORD(J,IATM)=XYZ(J)
   15 CONTINUE

*     CHECK WHETER OR NOT THE PRESENT BOND PAIR HAS ALREADY BEEN WRITTEN
      DO 60 I=1,IATM-1
         IF (NOATM(1,I) .EQ. IN(2) .AND. NOATM(2,I) .EQ. IN(4)) THEN
            DO 20 J=1,6
               IF (ABS(XYZ(J)-COORD(J,I)) .GT. TOL) GO TO 3
   20       CONTINUE
            RETURN
         END IF
    3    IF (NOATM(1,I) .EQ. IN(4) .AND. NOATM(2,I) .EQ. IN(2)) THEN
            DO 30 J=1,3
               IF (ABS(XYZ(J)-COORD(J+3,I)) .GT. TOL) GO TO 60
   30       CONTINUE
            DO 40 J=4,6
               IF (ABS(XYZ(J)-COORD(J-3,I)) .GT. TOL) GO TO 60
   40       CONTINUE
            RETURN
         END IF
   60 CONTINUE
      IATM=IATM+1
      IF (IATM .GT. NATM) WRITE(3,'(/A)') 'NATM IS TOO SMALL'
      
      LWR=1
      WRITE (3,117) (IN(I),I=2,5),ATOMN(IN(2)),(XYZ(I),I=1,3),
     &  ATOMN(IN(4)),(XYZ(I),I=4,6)
**117 FORMAT(12H0          (I2,1H,I3,3H) (I2,1H,I3,1H),
**   &  24X,A9,' (',F7.4,',',F8.4,',',F8.4,');   ',A9,' (',
**   &  F7.4,',',F8.4,',',F8.4,')')
C  117 FORMAT(12H0          (I2,1H,I5,3H) (I2,1H,I5,1H),
C     &  24X,A9,' (',F7.4,',',F8.4,',',F8.4,');   ',A9,' (',
C     &  F7.4,',',F8.4,',',F8.4,')')
  117 FORMAT('0',10X,'(',I2,',',I5,') (',I2,',',I5,')',
     &  24X,A9,' (',F7.4,',',F8.4,',',F8.4,');   ',A9,' (',
     &  F7.4,',',F8.4,',',F8.4,')')
      END

************************************************************************

      SUBROUTINE CALCD(P,C,S,TS,IS,XYZ)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*     CALCULATE FRACTIONAL COORDINATES OF AN ATOM
      PARAMETER (NW=384)
      DOUBLE PRECISION P(*),TS(3,NW),XYZ(3)
      INTEGER C,S,IS(2,3,NW)

*     ROTATION AND TRANSLATION
      IF (S .EQ. 0) THEN
         XYZ(1)=P(1)
         XYZ(2)=P(2)
         XYZ(3)=P(3)
      ELSE
         DO 20 I=1,3
            XYZ(I) = TS(I,S) + SIGN(1,IS(1,I,S)) * P(ABS(IS(1,I,S)))
            IF (IS(2,I,S) .NE. 0) THEN
               XYZ(I) = XYZ(I) + SIGN(1,IS(2,I,S)) * P(ABS(IS(2,I,S)))
            END IF
   20    CONTINUE
      END IF
      
*     CELL TRANSLATION
** MODIFICATION BEGIN
      SELECT CASE (C)
         CASE (0)
            RETURN
         CASE (1)
            XYZ(1) = XYZ(1) - 1.0
         CASE (2)
            XYZ(1) = XYZ(1) - 1.0
	    XYZ(2) = XYZ(2) + 1.0
         CASE (3)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (4)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (5)
            XYZ(2) = XYZ(2) - 1.0
         CASE (6)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) - 1.0
         CASE (7)
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (8)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (9)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) - 1.0
         CASE (10)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (11)
            XYZ(3) = XYZ(3) - 1.0
         CASE (12)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (13)
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (14)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (15)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (16)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (17)
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (18)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (19)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (20)
            XYZ(1) = XYZ(1) + 1.0
         CASE (21)
            XYZ(2) = XYZ(2) + 1.0
         CASE (22)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) + 1.0
         CASE (23)
            XYZ(3) = XYZ(3) + 1.0
         CASE (24)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (25)
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (26)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (27)
            XYZ(1) = XYZ(1) - 2.0
         CASE (28)
            XYZ(1) = XYZ(1) + 2.0
         CASE (29)
            XYZ(2) = XYZ(2) - 2.0
         CASE (30)
            XYZ(2) = XYZ(2) + 2.0
         CASE (31)
            XYZ(3) = XYZ(3) - 2.0
         CASE (32)
            XYZ(3) = XYZ(3) + 2.0
         CASE (33)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) - 2.0
         CASE (34)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) - 1.0
         CASE (35)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) + 1.0
         CASE (36)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) + 2.0
         CASE (37)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (38)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (39)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (40)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (41)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) - 2.0
         CASE (42)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) + 2.0
         CASE (43)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (44)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (45)
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (46)
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (47)
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (48)
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (49)
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (50)
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (51)
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (52)
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (53)
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (54)
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (55)
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (56)
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (57)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) - 2.0
         CASE (58)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) + 2.0
         CASE (59)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (60)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (61)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) - 2.0
         CASE (62)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) - 1.0
         CASE (63)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) + 1.0
         CASE (64)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) + 2.0
         CASE (65)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (66)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (67)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (68)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (69)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (70)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (71)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (72)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (73)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (74)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (75)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (76)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (77)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (78)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (79)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (80)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (81)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (82)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (83)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (84)
            XYZ(1) = XYZ(1) - 2.0
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (85)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (86)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (87)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (88)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (89)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (90)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (91)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (92)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (93)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (94)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (95)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (96)
            XYZ(1) = XYZ(1) - 1.0
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (97)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (98)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (99)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (100)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (101)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (102)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (103)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (104)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (105)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (106)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (107)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (108)
            XYZ(1) = XYZ(1) + 1.0
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (109)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (110)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (111)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (112)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) - 2.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (113)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (114)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (115)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (116)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) - 1.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (117)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (118)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (119)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (120)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) + 1.0
            XYZ(3) = XYZ(3) + 2.0
         CASE (121)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) - 2.0
         CASE (122)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) - 1.0
         CASE (123)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) + 1.0
         CASE (124)
            XYZ(1) = XYZ(1) + 2.0
            XYZ(2) = XYZ(2) + 2.0
            XYZ(3) = XYZ(3) + 2.0
      END SELECT
** MODIFICATION END
      END

C**********************************************************************C
C                             FUNA                                     C
C            ANGLE SUBROUTINE USED BY FUN2, FUN5, FUN6                 C
C**********************************************************************C
      FUNCTION FUNA(I)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      DIMENSION X1(3),X2(3),X3(3),V1(3),V2(3)
      CALL STOAA
      CALL ATOM(I,  X1)
      CALL ATOM(I+2,X2)
      CALL ATOM(I+4,X3)
      IF(NG)123,117,123
  117 CALL DIFV(X1,X2,V1)
      CALL DIFV(X3,X2,V2)
      FUNA=ACOS(COSVV(V1,V2))*57.29577951
  123 RETURN
      END
C**********************************************************************C
C                             FUNB                                     C
C               SET UP MATRIX AND GET EIGENVALUE                       C
C**********************************************************************C
      SUBROUTINE FUNB(W,Z,Z1)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      DIMENSION B(3,3),W(3,3),Y(3)
      CALL STOAA
      CALL BETA(2,B)
      IF(NG)123,113,123
  113 CALL MM(B,AA,W)
      CALL EIGVAL(W,Y)
      I=IN(4)
      Z=Y(I)
      Z1=SQRT (Z*0.0506605918)
  123 RETURN
      END
C**********************************************************************C
C                             FUNC                                     C
C             COS ANGLE OF PRINCIPAL AXIS AND VECTOR                   C
C**********************************************************************C
      SUBROUTINE FUNC(C,Z)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      DIMENSION W(3,3),X1(3),X2(3),V1(3),V2(3)
      CALL FUNB(W,Y,Z)
      IF(NG)125,111,125
  111 CALL EIGVEC(W,Y,V1)
      IF(NG)125,115,125
  115 CALL ATOM(5,X1)
      CALL ATOM(7,X2)
      IF(NG)125,121,125
  121 CALL DIFV(X2,X1,V2)
      C=COSVV(V1,V2)
  125 RETURN
      END
C**********************************************************************C
C                               FUND                                   C
C            DISTANCE SUBROUTINE USED BY FUN1 AND FUN4                 C
C**********************************************************************C
      FUNCTION FUND(I)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      DIMENSION X1(3),X2(3),V(3)
      CALL STOAA
      CALL ATOM(I,X1)
      CALL ATOM(I+2,X2)
      CALL DIFV(X2,X1,V)
      FUND=SQRT (VMV(V,AA,V))
      RETURN
      END
C**********************************************************************C
C                               FUNX                                   C
C           COS ANGLE OF PRINCIPAL AND CARTESIAN AXES                  C
C**********************************************************************C
      SUBROUTINE FUNX(C,Z)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      DIMENSION W(3,3),V(3),V1(3),V2(3),AX(3,3),X1(3),X2(3),X3(3),X4(3)
      DIMENSION AX3(3)
      CALL FUNB(W,Y,Z)
      CALL STOBB
      IF(NG)207,113,207
  113 CALL EIGVEC(W,Y,V)
      IF(NG)207,117,207
  117 CALL ATOM( 6,X1)
      CALL ATOM( 8,X2)
      CALL ATOM(10,X3)
      CALL ATOM(12,X4)
      IF(NG)207,123,207
  123 CALL DIFV(X2,X1,V1)
      CALL DIFV(X4,X3,V2)
      CALL AXES(V1,V2,AX)
      IN5=IN(5)
      DO 125 I=1,3
  125 AX3(I)=AX(I,IN5)
      C=COSVV(V,AX3)
  207 RETURN
      END
C**********************************************************************C
C                               FUNR                                   C
C                 MEAN SQUARE RADIAL DISPLACEMENT                      C
C**********************************************************************C
      SUBROUTINE FUNR(I,RSQ)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      DIMENSION B(3,3),BAA(3,3)
      CALL STOAA
      CALL BETA(I,B)
      CALL MM(B,AA,BAA)
      RSQ=TRASE(BAA)*0.0506605918
      RETURN
      END
C**********************************************************************C
C                               FUNCR                                  C
C             COMPUTE QUANTITIES FOR MEAN BOND DISTANCE                C
C**********************************************************************C
      SUBROUTINE FUNCR(C,R)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      DIMENSION C(2),X1(3),X2(3),V(3)
      DO 115 I=1,2
      CALL FUNR(2*I,RSQ)
      CALL FUNXI(2*I,2,XISQ)
  115 C(I)=RSQ-XISQ
      CALL ATOM(2,X1)
      CALL ATOM(4,X2)
      CALL DIFV(X2,X1,V)
      R=SQRT (VMV(V,AA,V))
      RETURN
      END
C**********************************************************************C
C                               FUNXI                                  C
C            MEAN SQUARE DISPLACEMENT IN GIVEN DIRECTION               C
C**********************************************************************C
      SUBROUTINE FUNXI(I,J,XISQ)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      DIMENSION B(3,3),X1(3),X2(3),V(3),BAA(3,3),AABAA(3,3)
      CALL STOAA
      CALL BETA(I,B)
      CALL ATOM(J,X1)
      CALL ATOM(J+2,X2)
      IF(NG)207,117,207
  117 CALL DIFV(X2,X1,V)
      D=VMV(V,AA,V)
      IF(D)123,123,201
  123 NG=10
      GO TO 207
  201 CALL MM(B,AA,BAA)
      CALL MM(AA,BAA,AABAA)
      XISQ=VMV(V,AABAA,V)*0.0506605918/D
  207 RETURN
      END
C**********************************************************************C
C                              FUNBE                                   C
C                          MV13 13.01.1983                             C
C**********************************************************************C
      SUBROUTINE FUNBE
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      DIMENSION X(3,3)
      CALL STOAA
      CALL BETA(2,X)
      FX=0.0
      DO 20 I=1,3
      DO 10 J=1,3
      FX=FX+X(I,J)*AA(I,J)
   10 CONTINUE
   20 CONTINUE
      FX=4.0*FX/3.0
      RETURN
      END
C**********************************************************************C
C                              NORM                                    C
C           STORE A VECTOR Z NORMAL TO VECTORS X AND Y                 C
C**********************************************************************C
      SUBROUTINE NORM(X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      DIMENSION X(3),Y(3),Z(3),X1(6),Y1(6),Z1(3)
      DO 115 I=1,3
      X1(I)=X(I)
      X1(I+3)=X(I)
      Y1(I)=Y(I)
  115 Y1(I+3)=Y(I)
      DO 119 I=1,3
  119 Z1(I)=X1(I+1)*Y1(I+2)-X1(I+2)*Y1(I+1)
      CALL MV(BB,Z1,Z)
      RETURN
      END
C**********************************************************************C
C                               AXES                                   C
C               STORE THREE MUTUALLY PERPENDICULAR                     C
C            VECTORS X(I,1), X(I,2), AND X(I,3) GIVEN                  C
C                      VECTORS U AND V.                                C
C**********************************************************************C
      SUBROUTINE AXES(U,V,X)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION U(3),V(3),X(3,3),X2(3),X3(3)
      DO 113 I=1,3
  113 X(I,1)=U(I)
      CALL NORM(U,V,X2)
      CALL NORM(U,X2,X3)
      DO 114 I=1,3
      X(I,2)=X2(I)
  114 X(I,3)=X3(I)
      RETURN
      END
C**********************************************************************C
C                             EIGVAL                                   C
C                 FIND EIGENVALUES Y OF MATRIX W                       C
C**********************************************************************C
      SUBROUTINE EIGVAL(W,Y)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON NG
      DIMENSION W(3,3),Y(3),X(3),Z(6,6)
      DOUBLE PRECISION A, B, A3, B2, A27, B4
      DO 119 J=1,3
      DO 119 I=1,3
      Z1=W(I,J)
      Z(I,J)=Z1
      Z(I+3,J)=Z1
      Z(I,J+3)=Z1
  119 Z(I+3,J+3)=Z1
      P=0.0
      Q=0.0
      R=0.0
      DO 207 I=1,3
      P=P-Z(I,I)
      Q=Q+Z(I,I)*Z(I+1,I+1)-Z(I,I+1)*Z(I+1,I)
  207 R=R+Z(3,I)*Z(2,I+1)*Z(1,I+2)-Z(1,I)*Z(2,I+1)*Z(3,I+2)
      P3=P/3.0
      A=Q-P*P3
      B=2.0*P3*P3*P3-Q*P3+R
      B2=B/2.0
      A3=A/3.0
      B4=B2*B2
      A27=A3*A3*A3
      IF(B4+A27)303,303,215
  215 IF(B4+1.0020*A27)220,220,225
  220 A27=-B4
      GO TO 303
  225 NG=7
      GO TO 317
  303 PHI3=DASIN(DSQRT(1.0D0+(B4/A27)))/3.0D0
      IF(B)307,305,307
  305 B=0.0
  307 C=-DSIGN(2.0D0*DSQRT(-A3),B)
      X(1)=C*COS (PHI3)
      X(2)=C*COS (PHI3+4.188790205)
      X(3)=C*COS (PHI3+2.094395103)
      IF(B)311,313,313
  311 HOLD=X(1)
      X(1)=X(3)
      X(3)=HOLD
  313 DO 315 I=1,3
  315 Y(I)=X(I)-P3
  317 RETURN
      END
C**********************************************************************C
C                              EIGVEC                                  C
C                 COMPUTE EIGENVECTOR Z OF MATRIX                      C
C                      W GIVEN EIGENVALUE Y                            C
C**********************************************************************C
      SUBROUTINE EIGVEC(W,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON NG
      DIMENSION W(3,3),X(6,6),Z(3),P(3)
      DO 123 J=1,3
      DO 123 I=1,3
      X1=W(I,J)
      X(I,J)=X1
      X(I+3,J)=X1
      X(I,J+3)=X1
  123 X(I+3,J+3)=X1
      Y1=Y
      DO 209 I=1,3
      X(I,I)=X(I,I)-Y1
      X(I+3,I)=X(I+3,I)-Y1
      X(I,I+3)=X(I,I+3)-Y1
  209 X(I+3,I+3)=X(I+3,I+3)-Y1
      S1=0.0
      DO 307 I=1,3
      S=0.0
      DO 223 J=1,3
      PJ=X(I,J+1)*X(I+1,J+2)-X(I,J+2)*X(I+1,J+1)
      P(J)=PJ
  223 S=S+PJ*PJ
      IF(S-S1)307,307,301
  301 S1=S
      DO 305 J=1,3
  305 Z(J)=P(J)
  307 CONTINUE
      IF(S1)311,311,313
  311 NG=8
  313 RETURN
      END
C**********************************************************************C
C                              MM                                      C
C           MULTIPLY TWO MATRICES Z(3,3)=X(3,3)*Y(3,3)                 C
C**********************************************************************C
      SUBROUTINE MM(X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(3,3),Y(3,3),Z(3,3)
      DO 117 I=1,3
      DO 117 K=1,3
      Z(I,K)=0.0
      DO 117 J=1,3
  117 Z(I,K)=Z(I,K)+X(I,J)*Y(J,K)
      RETURN
      END
C**********************************************************************C
C                               MV                                     C
C                 MATRIX * VECTOR Z(3)=X(3,3)*Y(3)                     C
C**********************************************************************C
      SUBROUTINE MV(X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(3,3),Y(3),Z(3)
      DO 113 I=1,3
      Z(I)=0.0
      DO 113 J=1,3
  113 Z(I)=Z(I)+X(I,J)*Y(J)
      RETURN
      END
C**********************************************************************C
C                                 VM                                   C
C           TRANSPOSED VECTOR TIMES MATRIX Z(3)=X(3)*Y(3,3)            C
C**********************************************************************C
      SUBROUTINE VM(X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(3),Y(3,3),Z(3)
      DO 115 J=1,3
      Z(J)=0.0
      DO 115 I=1,3
  115 Z(J)=Z(J)+X(I)*Y(I,J)
      RETURN
      END
C**********************************************************************C
C                                VV                                    C
C              TRANSPOSED VECTOR * VECTOR VV=X(3)*Y(3)                 C
C**********************************************************************C
      FUNCTION VV(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(3),Y(3)
      VV=0.0
      DO 111 I=1,3
  111 VV=VV+X(I)*Y(I)
      RETURN
      END
C**********************************************************************C
C                                 VMV                                  C
C     TRANSPOSED VECTOR * MATRIX * VECTOR   VMV=W(3)*X(3,3)*Y(3)       C
C**********************************************************************C
      FUNCTION VMV(W,X,Y)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION W(3),X(3,3),Y(3),Z(3)
      CALL MV(X,Y,Z)
      VMV=VV(W,Z)
      RETURN
      END
C**********************************************************************C
C                              DIFV                                    C
C                 VECTOR - VECTOR Z(3)=X(3)-Y(3)                       C
C**********************************************************************C
      SUBROUTINE DIFV(X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(3),Y(3),Z(3)
      DO 111 I=1,3
  111 Z(I)=X(I)-Y(I)
      RETURN
      END
C**********************************************************************C
C                               SUMV                                   C
C          COMPUTE THE SUM OF TWO VECTORS Z(3)=X(3)+Y(3)               C
C**********************************************************************C
      SUBROUTINE SUMV(X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(3),Y(3),Z(3)
      DO 111 I=1,3
  111 Z(I)=X(I)+Y(I)
      RETURN
      END
C**********************************************************************C
C                              COSVV                                   C
C              COSINE OF ANGLE BETWEEN VECTORS X AND Y                 C
C**********************************************************************C
      FUNCTION COSVV(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER(NNV=2048,NNP=2048,NW=384)
      COMMON  NG,IPM,NP,IAM,NS,NV,JXP,JX,JBP,JB,NM,SCALE,ITF,NQ,
     &INSAVE,DMAX,NA,F,FX,LNZ,KK,KKD,SAVEA,SAVEP,VARA,VARP,E1,E,IH,NSP,
     &PM,DP,DFDP,AM,DA,DFDA,IN,TS,IS,A,AA,BB,P,ROW,KI1,KI2
      INTEGER*2 KI1,KI2
      DIMENSION DP(NNV),DFDP(NNV),AM(21),DA(6)
      DIMENSION DFDA(6),IN(231),TS(3,NW),IS(2,3,NW),A(6),AA(3,3),BB(3,3)
      DIMENSION P(NNP),ROW(6),KI1(NNP),KI2(NNP)
      DIMENSION PM(NNV*(NNV+1)/2)
      DIMENSION X(3),Y(3)
      D=SQRT (VMV(X,AA,X)*VMV(Y,AA,Y))
      IF(D)111,111,115
  111 NG=9
      GO TO 117
  115 COSVV=VMV(X,AA,Y)/D
  117 RETURN
      END
C**********************************************************************C
C                              TRASE                                   C
C                     COMPUTE TRACE OF MATRIX X                        C
C**********************************************************************C
      FUNCTION TRASE(X)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(3,3)
      TRASE=0.0
      DO 111 I=1,3
  111 TRASE=TRASE+X(I,I)
      RETURN
      END
C**********************************************************************C
C                              SETP                                    C
C              SET THE VALUES OF STRUCTURE PARAMETERS                  C
C**********************************************************************C
      SUBROUTINE SETP(P)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NCS=512)
      COMMON /CN/ NEWNO(NCS),C1(NCS),C2(10,NCS),INP(10,NCS),NCNPAR(NCS),
     &NCNSTR
      DOUBLE PRECISION P(*)

      DO 50 I=1,NCNSTR
         T=C1(I)
         DO 20 J=1,NCNPAR(I)
            T=T+C2(J,I)*P(INP(J,I))
   20    CONTINUE
         P(NEWNO(I))=T
   50 CONTINUE
      END
C**********************************************************************C
C                              SETA                                    C
C               SET THE VALUES OF LATTICE CONSTANTS                    C
C**********************************************************************C
      SUBROUTINE SETA(A)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /LG/ LAUEG
      DOUBLE PRECISION A(*)

      IF (LAUEG .EQ. 6 .OR. LAUEG .EQ. 7) THEN
         A(2)=A(1)
      ELSE IF (LAUEG .EQ. 8 .OR. LAUEG .EQ. 10) THEN
         A(2)=A(1)
         A(3)=A(1)
         A(5)=A(4)
         A(6)=A(4)
      ELSE IF (LAUEG .EQ. 9 .OR. LAUEG .EQ. 11 .OR. LAUEG .EQ. 12
     &.OR. LAUEG .EQ. 13) THEN
         A(2)=A(1)
      ELSE IF (LAUEG .GE. 14) THEN
         A(2)=A(1)
         A(3)=A(1)
      ENDIF
      END
