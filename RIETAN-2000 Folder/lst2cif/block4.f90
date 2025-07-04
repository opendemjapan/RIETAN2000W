C
C--------------------------------------------------------------
C     CRYSTAL DATA
C--------------------------------------------------------------
C
CAdditional information:
C  For ATOM SPAACIES 
C          Lines: from 145 to 168 for the new version of RIETAN2001
C                 from 170 to 190 for the RIETAN98
C =============================================================

      SUBROUTINE WRITE_BLOCK4(INP,OUT,TMP,ID_TRANS)

C =============================================================
C       NPH - MAXIMUM NUMBER OF PHASES 
C       MAT - MAXIMUM NUMBER OF ATOMS 
C -------------------------------------------------------------
C    NPHASE - NUMBER OF PHASES
C     PHAME - NAME OF PHASE
C        SG - SPACE GROUP NAME
C       SGV - SPACE GROUP DATABASE FILE (I - SPGRI / A - SPGRA)
C       SGN - SPACE GROUP NUMBER
C       SGS - SPACE GROUP SETTING
C      SYNG - SYMMETRY
C        PG - POINT GROUP NAME
C -------------------------------------------------------------
C         P - a,b,calpha,beta,gamma, Volume   
C        SP - 
C -------------------------------------------------------------
C     NATOM - NUMBER OF ATOMS
C     LATOM - LABELS OF ATOMS    
C    NATOMS - NUMBER OF ATOM SPECIES
C     ATOMT - TYPE OF ATOMS (ATOMS SPECIES)   
C    CATOMT - TYPE OF ATOMS    
C     ATOMN - NUMBER OF ATOM SPECIES IN CELL
C         M - MULTIPLICITY
C      PARA - g, x, y, z, Biso   
C     EPARA - 
C      CTMP - 
C -------------------------------------------------------------
C        BB - ANISOTROPIC THERMAL PARAMETERS Uij and Biso
C       EBB - 
C     CTMP2 - 
C  THERM_ID - TYPE OF THERMAL PARAMETERS (Uani / Biso)
C	IF THERM_ID = 'Biso' THEN ID_THERMAL_TYPE=0
C	IF THERM_ID = 'Uani' THEN ID_THERMAL_TYPE=1
C -------------------------------------------------------------
C   NCONSTR - NUMBER OF LINEAR CONSTRAINTS
C     SCON1 - LEFT CONSTRAIN
C     SCON2 - RIGHT CONSTRAIN   
C     CONL1  - ATOM LABEL (LEFT CONSTRAIN)
C     CONL2  - ATOM LABEL (RIGHT CONSTRAIN)
C =============================================================
      INTEGER INP,OUT,TMP,DB,ID_TRANS
	CHARACTER CMD*300,PAR*80, UPD*20
      CHARACTER DATE_TIME*16
	INTEGER ID(5)
      INTEGER NPH,MAT
      PARAMETER (NPH=8,MAT=100)
      INTEGER NPHASE
      CHARACTER PHNAME(NPH)*80,SG(NPH)*15,SGV(NPH)*1,SYNG(NPH)*15,
     &          PG(NPH)*10
      INTEGER SGN(NPH), SGS(NPH)
      REAL P(7,NPH), SP(7,NPH) 
      CHARACTER LATOM(MAT,NPH)*9,ATOMT(MAT,NPH)*4,CATOMT(MAT)*4
	CHARACTER CTMP*12,CTMP2*9
	INTEGER NATOMS
	CHARACTER ATOMN(MAT,NPH)*15 
      INTEGER NATOM(NPH), M(MAT,NPH)
      REAL PARA(7,MAT,NPH),EPARA(7,MAT,NPH)
      REAL BB(13,MAT,NPH),EBB(13,MAT,NPH)
      CHARACTER THERM_ID(MAT,NPH)*4 
	INTEGER ID_TERM(MAT)
	CHARACTER*16 SCON1,SCON2
	CHARACTER*9 CONL1,CONL2


      CALL GET_DATE_TIME(DATE_TIME)
C -------------------------------------------------------------
C     GET NUMBER OF PHASES
C -------------------------------------------------------------
      CMD= 'PHASE_BEGIN'
      CALL FINDSTRING(0,TMP,CMD,11,RESULT)
      CALL SKIPLINE(TMP,1)
      READ(TMP,*,ERR=991) NPHASE

C -------------------------------------------------------------
C     GET PHNAME,SG,SGV,SGN,SGS,SYNG,PG
C -------------------------------------------------------------
      DO I=1,NPHASE
         CMD='Phase #'
         CALL FINDSTRING(1,TMP,CMD,7,RESULT)
         CALL SKIPLINE(TMP,1)
         READ(TMP,'(A)',ERR=991) PHNAME(I)
         READ(TMP,'(A)',ERR=991) SG(I)
         READ(TMP,'(A)',ERR=991) SGV(I)
         READ(TMP,*,ERR=991) SGN(I),SGS(I)
         READ(TMP,'(A)',ERR=991) SYNG(I)
         READ(TMP,'(A)',ERR=991) PG(I)
      END DO
      GOTO 1
991   CALL ENDPROG('ERROR IN <PHASE> BLOCK')

C -------------------------------------------------------------
C     GET UNIT CELL PARAMETERS
C -------------------------------------------------------------
  1   CONTINUE 
      CMD='CELL_BEGIN'
      CALL FINDSTRING(0,TMP,CMD,10,RESULT)
      CALL SKIPLINE(TMP,1)
      CMD='Phase #'
      DO I=1,NPHASE
         CALL FINDSTRING(1,TMP,CMD,7,RESULT)
         CALL SKIPLINE(TMP,1)
         READ(TMP,'(3F10.5,3F10.4,F12.4)',ERR=992) (P(J,I), J=1,7)
         READ(TMP,'(3F10.5,3F10.4,F12.4)',ERR=992) (SP(J,I), J=1,7)
      END DO
      GOTO 2
992   CALL ENDPROG('ERROR IN <CELL> BLOCK')

C -------------------------------------------------------------
C     GET atom, m, g, x, y, z, Biso
C -------------------------------------------------------------
  2   CONTINUE 
      CMD= 'ATOM_BEGIN'
      CALL FINDSTRING(0,TMP,CMD,10,RESULT)
      CALL SKIPLINE(TMP,1)
      DO I=1,NPHASE
         CALL SKIPLINE(TMP,1)
         READ(TMP,*,ERR=993) NATOM(I)
         DO J=1,NATOM(I)
            READ(TMP,'(A,I3,F9.4,3F10.5,F8.3)',ERR=993) 
     &          LATOM(J,I),M(J,I),(PARA(K,J,I), K=1,5)
            READ(TMP,'(A,F9.4,3F10.5,F8.3)',ERR=993) 
     &          CTMP, (EPARA(K,J,I),K=1,5)
         ENDDO
      ENDDO
      GOTO 3
993   CALL ENDPROG('ERROR IN <ATOM> BLOCK')

C -------------------------------------------------------------
C     GET ATOMS SPECIES
C -------------------------------------------------------------
  3   CONTINUE 
C =================================================================
C       FOR NEW VERSION OF RIETAN20001
C =================================================================
      REWIND(INP)
      DO I=1,NPHASE
	  J = 0
31	  PAR = ''
        READ(INP,'(A)') PAR
        IF ((PAR(37:37).EQ.':').AND.(INDEX(PAR(24:35),'/').GT.0)) THEN
	    DO K=1,NATOM(I)
	      I1 = INDEX(PAR(23:34),LATOM(K,I)(1:LEN_TRIM(LATOM(K,I))))
	      IF(I1.GT.0) THEN
	        J = J + 1
	        I1 = INDEX(PAR,'/') + 1
	        I2 = INDEX(PAR,':') - 1
* Corrected by R. Dilanian on 2001.06.07 { 
C	  	    ATOMT(K,I) = PAR(I1:I2)
	  	    ATOMT(J,I) = PAR(I1:I2)
* }
	        GOTO 32
	      ENDIF
	    ENDDO
32        CONTINUE
	  ENDIF
	  IF(J.LT.NATOM(I)) GOTO 31
	ENDDO
      CLOSE(UNIT=INP)

C =================================================================
C       FOR RIETAN98
C =================================================================
CC      CMD= 'ATYPE_BEGIN'
CC      CALL FINDSTRING(0,TMP,CMD,11,RESULT)
CC      CALL SKIPLINE(TMP,1)
CC      DO I=1,NPHASE
CC         CALL SKIPLINE(TMP,1)
CC         READ(TMP,*,ERR=994) NATOMS
CC         DO J=1,NATOMS
CC           READ(TMP,'(A4,A12)',ERR=994) CATOMT(J),ATOMN(J,I)
CC         ENDDO
CC	   DO J=1,NATOM(I)
CC	     DO K=1,NATOMS
CC	        I1=INDEX(LATOM(J,I),CATOMT(K)(1:LEN_TRIM(CATOMT(K))))
CC	        IF(I1.GT.0) ATOMT(J,I)=CATOMT(K) 
CC		 ENDDO 
CC	   ENDDO
CC      ENDDO
CC      GOTO 4
CC994   CALL ENDPROG('ERROR IN <ATOME> BLOCK')
    
C -------------------------------------------------------------
C     GET ANISOTROPIC TERMAL PARAMETERS
C -------------------------------------------------------------
  4   CONTINUE 
      CMD= 'THERMAL_BEGIN'
      CALL FINDSTRING(0,TMP,CMD,13,RESULT)
      CALL SKIPLINE(TMP,1)
      DO I=1,NPHASE
         CALL SKIPLINE(TMP,1)
         READ(TMP,*,ERR=995) NATOM(I)
         DO J=1,NATOM(I)
            READ(TMP,'(A,12F10.5,F8.3)',ERR=995) 
     &          LATOM(J,I),(BB(K,J,I),K=1,13)
            READ(TMP,'(A9,12F10.5,F8.3)',ERR=995) 
     &          CTMP2,(EBB(K,J,I),K=1,13)
         ENDDO
      ENDDO
      GOTO 5
995   CALL ENDPROG('ERROR IN <THERMAL> BLOCK')
  5   CONTINUE 
C -------------------------------------------------------------
C     GET LINEAR CONSTRAINTS FOR g,x,y,z,Biso,Uij
C -------------------------------------------------------------
      CMD= 'CONSTR_BEGIN'
      CALL FINDSTRING(0,TMP,CMD,12,RESULT)
      CALL SKIPLINE(TMP,1)
	READ(TMP,*) NCONSTR
	IF (NCONSTR.GT.0) THEN
	 DO K=1,NCONSTR
	    READ(TMP,'(A)') PAR
	    I1=INDEX(PAR,'=')
	    SCON1 = PAR(1:I1-1)
	    SCON2 = PAR(I1+1:)
	    I11 = INDEX(SCON1,'(')
	    I12 = INDEX(SCON1,',')
	    I21 = INDEX(SCON2,'(')
	    I22 = INDEX(SCON2,',')
	    CONL1=SCON1(I11+1:I12-1)
	    CONL2=SCON2(I21+1:I22-1)
	    K1=0
	    IF(INDEX(SCON1,'g').GT.0) K1 = 1
	    IF(INDEX(SCON1,'x').GT.0) K1 = 2
	    IF(INDEX(SCON1,'y').GT.0) K1 = 3
	    IF(INDEX(SCON1,'z').GT.0) K1 = 4
	    IF(INDEX(SCON1,'B').GT.0) K1 = 5
	    IF(INDEX(SCON1,'U11').GT.0) K1 = 7
	    IF(INDEX(SCON1,'U22').GT.0) K1 = 8
	    IF(INDEX(SCON1,'U33').GT.0) K1 = 9
	    IF(INDEX(SCON1,'U12').GT.0) K1 = 10
	    IF(INDEX(SCON1,'U13').GT.0) K1 = 11
	    IF(INDEX(SCON1,'U23').GT.0) K1 = 12
	    IF(INDEX(SCON2,'g').GT.0) K2 = 1
	    IF(INDEX(SCON2,'x').GT.0) K2 = 2
	    IF(INDEX(SCON2,'y').GT.0) K2 = 3
	    IF(INDEX(SCON2,'z').GT.0) K2 = 4
	    IF(INDEX(SCON2,'B').GT.0) K2 = 5
	    IF(INDEX(SCON2,'U11').GT.0) K2 = 7
	    IF(INDEX(SCON2,'U22').GT.0) K2 = 8
	    IF(INDEX(SCON2,'U33').GT.0) K2 = 9
	    IF(INDEX(SCON2,'U12').GT.0) K2 = 10
	    IF(INDEX(SCON2,'U13').GT.0) K2 = 11
	    IF(INDEX(SCON2,'U23').GT.0) K2 = 12
	    DO I1=1,NPHASE
	      DO J1=1,NATOM(I1)
             IF(LATOM(J1,I1)(1:LEN_TRIM(LATOM(J1,I1))) .EQ.
     &       CONL1(1:LEN_TRIM(CONL1))) GOTO 51
	      ENDDO
	    ENDDO
51       CONTINUE
	    DO I2=1,NPHASE
	      DO J2=1,NATOM(I2)
 	        IF(LATOM(J2,I2)(1:LEN_TRIM(LATOM(J2,I2))).EQ.
     &       CONL2(1:LEN_TRIM(CONL2))) GOTO 52
	      ENDDO
	    ENDDO
52        CONTINUE
          SELECT CASE(K1)
* Corrected by R. Dilanian on 2001.04.12 { 
C	      CASE(0)
C	         CALL ENDPROG('ERROR IN <CONSTR> BLOCK')
*}
	      CASE(1,2,3,4,5)
	         EPARA(K1,J1,I1) = EPARA(K2,J2,I2)
	      CASE(7,8,9,10,11,12)
	         EBB(K1,J1,I1) = EBB(K2,J2,I2)
	    END SELECT
	 ENDDO
	ENDIF

C -------------------------------------------------------------
C     WRITE PARAMETERS TO CIF FILE
C -------------------------------------------------------------
      WRITE(OUT,900)
      WRITE(OUT,902)
      WRITE(OUT,900)
      WRITE(OUT,400)
      WRITE(OUT,900)
      DO I=1,NPHASE
        WRITE(OUT,904)
        WRITE(OUT,401) I
        WRITE(OUT,900)
        WRITE(OUT,402)
        WRITE(OUT,403) DATE_TIME(1:10),I
        WRITE(OUT,900)
        WRITE(OUT,404) PHNAME(I)(1:LEN_TRIM(PHNAME(I)))
C ------------------------------------------------------- A
	  CALL UPDE(P(1,I),SP(1,I),1,UPD,IS)
        WRITE(OUT,405) UPD(IS:)
C ------------------------------------------------------- B
	  CALL UPDE(P(2,I),SP(2,I),1,UPD,IS)
        WRITE(OUT,406) UPD(IS:)
C ------------------------------------------------------- C
	  CALL UPDE(P(3,I),SP(3,I),1,UPD,IS)
        WRITE(OUT,407) UPD(IS:)
C ------------------------------------------------------- ALPHA
	  IS=1
	  IF((P(4,I).EQ.90.0).OR.(P(4,I).EQ.120.0)) THEN
	    IF(P(4,I).EQ.90.0) UPD = '90.0'
	    IF(P(4,I).EQ.120.0) UPD = '120.0'
	  ELSE
	    CALL UPDE(P(4,I),SP(4,I),2,UPD,IS)
	  ENDIF
        WRITE(OUT,408) UPD(IS:)
C ------------------------------------------------------- BETA
	  IS=1
	  IF((P(5,I).EQ.90.0).OR.(P(5,I).EQ.120.0)) THEN
	    IF(P(5,I).EQ.90.0) UPD = '90.0'
	    IF(P(5,I).EQ.120.0) UPD = '120.0'
	  ELSE
	    CALL UPDE(P(5,I),SP(5,I),2,UPD,IS)
	  ENDIF
        WRITE(OUT,409) UPD(IS:)
C ------------------------------------------------------- GAMMA
	  IS=1
	  IF((P(6,I).EQ.90.0).OR.(P(6,I).EQ.120.0)) THEN
	    IF(P(6,I).EQ.90.0) UPD = '90.0'
	    IF(P(6,I).EQ.120.0) UPD = '120.0'
	  ELSE
	    CALL UPDE(P(6,I),SP(6,I),2,UPD,IS)
	  ENDIF
        WRITE(OUT,410) UPD(IS:)
C ------------------------------------------------------- VOLUME
        CALL UPDE(P(7,I),SP(7,I),3,UPD,IS)
        WRITE(OUT,411) UPD(IS:)

        WRITE(OUT,412)
        WRITE(OUT,413) SYNG(I)
        WRITE(OUT,414) SG(I)(1:LEN_TRIM(SG(I)))
        WRITE(OUT,415) SGN(I)

C -------------------------------------------- _symmetry_equiv_pos_as_xyz  
        WRITE(OUT,900)
        WRITE(OUT,903)
        WRITE(OUT,416)
        WRITE(OUT,417)
        CALL OPENDATABASE(SGV(I),DB)
 6      READ(DB,'(5I3,A)',END=5) (ID(J),J=1,5),PAR
        IF (ID(1).EQ.SGN(I) .AND. ID(2).EQ.SGS(I)) THEN 
          GOTO 7
        ELSE
          CALL SKIPLINE(DB,ID(5)+1)
          GOTO 6
        END IF
 7      CALL SKIPLINE(DB,1)
        CALL SIMPOS(OUT,DB,ID(3),ID(4),ID(5),SG(I)(1:LEN_TRIM(SG(I)))
     &	          ,ID_TRANS)         
        CLOSE(UNIT=DB)

C --------------------------------------------------------------------------
C    GET _thermal_displace_type AND Biso(ESD) [IF _thermal_displace_type=Uani]
C --------------------------------------------------------------------------
        ID_THERMAL_TYPE=0 
        DO J=1,NATOM(I) 
          IF (PARA(5,J,I).EQ.0.0) THEN
             THERM_ID(J,I)='Uani'                
             PARA(5,J,I)=BB(13,J,I)
             EPARA(5,J,I)=EBB(13,J,I)
	       ID_TERM(J) = 1
	       ID_THERMAL_TYPE = 1
          END IF
        END DO
        DO J=1,NATOM(I)
          IF (BB(13,J,I).EQ.0.0) THEN 
             THERM_ID(J,I)='Biso'
             BB(13,J,I)=PARA(5,J,I)
             EBB(13,J,I)=EPARA(5,J,I)             
	       ID_TERM(J) = 0
          END IF 
        END DO
        DO J=1,NATOM(I)
          IF (THERM_ID(J,I).EQ.'Uani') THEN
            CALL BESD(P(1,I),P(2,I),P(3,I),P(4,I),P(5,I),P(6,I),
     &       SP(1,I),SP(2,I),SP(3,I),SP(4,I),SP(5,I),SP(6,I),
     &       BB(1,J,I),BB(2,J,I),BB(3,J,I),BB(4,J,I),BB(5,J,I),
     &       BB(6,J,I),EBB(1,J,I),EBB(2,J,I),EBB(3,J,I),EBB(4,J,I),
     &       EBB(5,J,I),EBB(6,J,I),BB(13,J,I),EBB(13,J,I))
            EPARA(5,J,I)=EBB(13,J,I) 
          END IF
        END DO                   
C ------------------------ atom_label,M,g,x,y,z,thermal_type,Biso,atom_type  

C ================================================================
C                            FORMAT					   MAX CHAR
C ================================================================
C ------------------------------------------- SPACE    --  3 CHAR           
C   atom_label : xxxxxxxxx                             --  6 CHAR
C ------------------------------------------- SPACE    --  1 CHAR           
C            M : xxx                                   --  3 CHAR
C ------------------------------------------- SPACE    --  1 CHAR           
C            g : x.xxxx(x)                             --  9 CHAR
C ------------------------------------------- SPACE    --  1 CHAR           
C            x : x.xxxxx(x)                            -- 11 CHAR
C ------------------------------------------- SPACE    --  1 CHAR           
C            y : x.xxxxx(x)                            -- 11 CHAR
C ------------------------------------------- SPACE    --  1 CHAR           
C            z : x.xxxxx(x)                            -- 11 CHAR
C ------------------------------------------- SPACE    --  1 CHAR           
C thermal_type : xxxx                                  --  4 CHAR
C ------------------------------------------- SPACE    --  1 CHAR           
C         Biso : xxx.xxx(x)                            -- 10 CHAR
C ------------------------------------------- SPACE    --  1 CHAR           
C    atom_type : xxxx                                  --  4 CHAR
C ================================================================
C  TOTAL                                               -- 80 CHAR   
C ================================================================
        WRITE(OUT,900)
        WRITE(OUT,903)
        WRITE(OUT,418)
        WRITE(OUT,419)
        WRITE(OUT,420)
        WRITE(OUT,421)
        WRITE(OUT,422)
        WRITE(OUT,423)
        WRITE(OUT,424)
        WRITE(OUT,425)
        WRITE(OUT,426)
        DO J=1,NATOM(I)
	    PAR = ''
C ------------------------------------------------------- atom_label
          PAR(4:9) = LATOM(J,I)
C ------------------------------------------------------- M
          WRITE(PAR(11:13),'(I3)') M(J,I)
C ------------------------------------------------------- g
          CALL UPDE(PARA(1,J,I),EPARA(1,J,I),4,UPD,IS)
* Corrected by R. Dilanian on 2001.04.13 { 
C..........................................................
C         CHECK g<1.0 OR g=1.0
C..........................................................
          L1 = 16
* Corrected by F. Izumi on 2001.06.12 { 
C         IF (INDEX(UPD,'1.').GT.0) L1 = 18 
          IF (INDEX(UPD,'1.') .GT. 0 .OR. INDEX(UPD,'0.') .GT. 0) 
     &    L1 = 15 
* }		  
C..........................................................
*}
	    IF(EPARA(1,J,I).EQ.0.0) THEN
	       K=INDEX(UPD,'.')
	       DO L=LEN_TRIM(UPD),K+2,-1
	         IF(UPD(L:L).NE.'0') GOTO 8
	       ENDDO
8            PAR(L1:23) = UPD(IS:L)
          ELSE
            PAR(L1:23) = UPD(IS:)
	    ENDIF
C ------------------------------------------------------- x
          CALL UPDE(PARA(2,J,I),EPARA(2,J,I),1,UPD,IS)
* Corrected by R. Dilanian on 2001.04.12 { 
C..........................................................
C         CHECK NEGATIVE OR NOT
C..........................................................
          L1 = 26
          IF ((INDEX(UPD,'-').GT.0).OR.(INDEX(UPD,'1.').GT.0)) L1 = 25 
C..........................................................
*}
	    IF(EPARA(2,J,I).EQ.0.0) THEN
	       K=INDEX(UPD,'.')
	       DO L=LEN_TRIM(UPD),K+2,-1
	         IF(UPD(L:L).NE.'0') GOTO 9
	       ENDDO
9            PAR(L1:35) = UPD(IS:L)
          ELSE
            PAR(L1:35) = UPD(IS:)
	    ENDIF
C ------------------------------------------------------- y
          CALL UPDE(PARA(3,J,I),EPARA(3,J,I),1,UPD,IS)
* Corrected by R. Dilanian on 2001.04.12 { 
C..........................................................
C         CHECK NEGATIVE OR NOT
C..........................................................
          L1 = 38
          IF ((INDEX(UPD,'-').GT.0).OR.(INDEX(UPD,'1.').GT.0)) L1 = 37 
C..........................................................
*}
	    IF(EPARA(3,J,I).EQ.0.0) THEN
	       K=INDEX(UPD,'.')
	       DO L=LEN_TRIM(UPD),K+2,-1
	         IF(UPD(L:L).NE.'0') GOTO 10
	       ENDDO
10           PAR(L1:47) = UPD(IS:L)
          ELSE
            PAR(L1:47) = UPD(IS:)
	    ENDIF
C ------------------------------------------------------- z
          CALL UPDE(PARA(4,J,I),EPARA(4,J,I),1,UPD,IS)
* Corrected by R. Dilanian on 2001.04.12 { 
C..........................................................
C         CHECK NEGATIVE OR NOT
C..........................................................
          L1 = 50
          IF ((INDEX(UPD,'-').GT.0).OR.(INDEX(UPD,'1.').GT.0)) L1 = 49 
C..........................................................
*}
	    IF(EPARA(4,J,I).EQ.0.0) THEN
	       K=INDEX(UPD,'.')
	       DO L=LEN_TRIM(UPD),K+2,-1
	         IF(UPD(L:L).NE.'0') GOTO 11
	       ENDDO
11            PAR(L1:59) = UPD(IS:L)
          ELSE
            PAR(L1:59) = UPD(IS:)
	    ENDIF
C ------------------------------------------------------- thermal_type
          PAR(61:64) = THERM_ID(J,I)
C ------------------------------------------------------- Biso
          L1 = 67
          CALL UPDE(PARA(5,J,I),EPARA(5,J,I),5,UPD,IS)
          IF (INDEX(UPD,'-').GT.0) L1 = 66 
 	    IF(EPARA(5,J,I).EQ.0.0) THEN
	       K=INDEX(UPD,'.')
	       DO L=LEN_TRIM(UPD),K+2,-1
	         IF(UPD(L:L).NE.'0') GOTO 12
	       ENDDO
12            PAR(L1:) = UPD(IS:L)
          ELSE
            PAR(L1:) = UPD(IS:)
	    ENDIF
C ------------------------------------------------------- atom_type
          PAR(77:) = ATOMT(J,I)
          WRITE(OUT,'(A)') PAR
C ------------------------------------------------------- 
 	  ENDDO

	  IF(ID_THERMAL_TYPE.EQ.1) THEN
C -------------------------------------------- atom_label,Uij  

C ================================================================
C                            FORMAT					   MAX CHAR
C ================================================================
C ------------------------------------------- SPACE    --  3 CHAR           
C   atom_label : xxxxxxxxx                             --  9 CHAR
C ------------------------------------------- SPACE    --  1 CHAR           
C          U11 : xx.xxxx(x)                            -- 10 CHAR
C ------------------------------------------- SPACE    --  1 CHAR           
C          U22 : xx.xxxx(x)                            -- 10 CHAR
C ------------------------------------------- SPACE    --  1 CHAR           
C          U33 : xx.xxxx(x)                            -- 10 CHAR
C ------------------------------------------- SPACE    --  1 CHAR           
C          U12 : xx.xxxx(x)                            -- 10 CHAR
C ------------------------------------------- SPACE    --  1 CHAR           
C          U13 : xx.xxxx(x)                            -- 10 CHAR
C ------------------------------------------- SPACE    --  1 CHAR           
C          U23 : xx.xxxx(x)                            -- 10 CHAR
C ================================================================
C  TOTAL                                               -- 78 CHAR   
C ================================================================
          WRITE(OUT,900)
          WRITE(OUT,903)
          WRITE(OUT,427)
          WRITE(OUT,428)
          WRITE(OUT,429)
          WRITE(OUT,430)
          WRITE(OUT,431)
          WRITE(OUT,432)
          WRITE(OUT,433)
          DO J=1,NATOM(I)
	      IF(ID_TERM(J).EQ.1) THEN
C ===============================================================
C   FOR Uij AND ESD(Uij): change format from F10.5 to F10.4   
              DO I1=7,12
	          IF((EBB(I1,J,K)*100000).LT.10) EBB(I1,J,K)=0.00010
	        ENDDO
C ===============================================================
	        PAR = ''
C ------------------------------------------------------- _atom_site_aniso_label
              PAR(4:12) = LATOM(J,I)
C ------------------------------------------------------- U11,U22,U33,U12,U13,U23
	        I1 = 14
			I2 = I1 + 11 
	        DO IM=1,6
                CALL UPDE(BB(6+IM,J,I),EBB(6+IM,J,I),1,UPD,IS)
	          IF(EBB(6+IM,J,I).EQ.0.0) THEN
	            K=INDEX(UPD,'.')
	            DO L=LEN_TRIM(UPD),K+2,-1
	              IF(UPD(L:L).NE.'0') GOTO 13
	            ENDDO
13                PAR(I1:I2) = UPD(IS:L)
                ELSE
                  PAR(I1:I2) = UPD(IS:)
	          ENDIF
	          I1 = I2
	          I2 = I1 + 11
	        ENDDO
              WRITE(OUT,'(A)') PAR
	      ENDIF
	    ENDDO
        ENDIF
	ENDDO


C **************************************************************
C--------------------------------------------------------------------
C     FORMATS FOR CRYSTAL DATA
C--------------------------------------------------------------------

 9999 CONTINUE
 400  FORMAT('# CRYSTAL DATA') 
 401  FORMAT('data_RIETAN_phase_',I1) 
 402  FORMAT('_pd_block_id') 
 403  FORMAT(6X,"'",A,'|','PHASE_0',I1,'|','..creator_name..','|',
     &       '..instr_name..',"'")

 404  FORMAT('_pd_phase_name',25X,A) 
 405  FORMAT('_cell_length_a',25X,A) 
 406  FORMAT('_cell_length_b',25X,A) 
 407  FORMAT('_cell_length_c',25X,A) 
 408  FORMAT('_cell_angle_alpha',22X,A) 
 409  FORMAT('_cell_angle_beta',23X,A) 
 410  FORMAT('_cell_angle_gamma',22X,A) 
 411  FORMAT('_cell_volume',27X,A) 
 412  FORMAT('_cell_formula_units',20X,'?') 
 413  FORMAT('_symmetry_cell_setting',17X,A) 
 414  FORMAT('_symmetry_space_group_name_H-M',9X,"'",A,"'") 
 415  FORMAT('_symmetry_Int_Tables_number',12X,I3)
  
 416  FORMAT(3X,'_symmetry_equiv_pos_site_id') 
 417  FORMAT(3X,'_symmetry_equiv_pos_as_xyz') 

 418  FORMAT(3X,'_atom_site_label') 
 419  FORMAT(3X,'_atom_site_symmetry_multiplicity') 
 420  FORMAT(3X,'_atom_site_occupancy') 
 421  FORMAT(3X,'_atom_site_fract_x') 
 422  FORMAT(3X,'_atom_site_fract_y') 
 423  FORMAT(3X,'_atom_site_fract_z') 
 424  FORMAT(3X,'_atom_site_thermal_displace_type') 
* Corrected by F. Izumi on 2001.04.10 { 
C425  FORMAT(3X,'_atom_site_U_iso_or_equiv') 
 425  FORMAT(3X,'_atom_site_B_iso_or_equiv')
* }
* Corrected by R. Dilanian on 2001.06.07 { 
C426  FORMAT(3X,'_atom_type_symbol') 
 426  FORMAT(3X,'_atom_site_type_symbol') 
* }
 427  FORMAT(3X,'_atom_site_aniso_label') 
 428  FORMAT(3X,'_atom_site_aniso_U11') 
 429  FORMAT(3X,'_atom_site_aniso_U22') 
 430  FORMAT(3X,'_atom_site_aniso_U33') 
 431  FORMAT(3X,'_atom_site_aniso_U12') 
 432  FORMAT(3X,'_atom_site_aniso_U13') 
 433  FORMAT(3X,'_atom_site_aniso_U23') 
 

C 433  FORMAT(3X,'_atom_type_symbol') 
c 434  FORMAT(3X,'_atom_type_number_in_cell') 

 435  FORMAT('#',1X,A,1X,A,1X,A) 
C--------------------------------------------------------------------
C
 900  FORMAT(' ')
 901  FORMAT(';')
 902  FORMAT('#',70('='))
 903  FORMAT('loop_')
 904  FORMAT('#',70('-'))
C
C **************************************************************
      END



C ---------------------------------------------------------------
C ---------------------------------------------------------------
      SUBROUTINE UPDE(X,Y,F,PAR,S)
C ---------------------------------------------------------------
C ---------------------------------------------------------------
	CHARACTER PAR1*20,PAR*20,FORMS*7
	INTEGER F,S
     
C--------------------------------------------------------------------
C       1  F10.5 - a,b,c ; x,y,z ; Uij
C       2  F10.4 - alpha,beta,gamma
C       3  F12.4 - volume
C       4  F9.4  - g
C       5  F8.3  - Biso
C--------------------------------------------------------------------
	SELECT CASE(F)
	  CASE(1)
	     CALL ROUNDF(F,5,FORMS)
	  CASE(2,3,4)
	     CALL ROUNDF(F,4,FORMS)
	  CASE(5)
	     CALL ROUNDF(F,3,FORMS)
	END SELECT

	IF(Y.EQ.0) GOTO 4
	WRITE(PAR1,FORMS) Y
	I1 = INDEX(PAR1,'.')+1
	DO I=I1,LEN(PAR1)
	  IF(PAR1(I:I).NE.'0') GOTO 1
	ENDDO
1	I2 = I-I1+1

	CALL ROUNDF(F,I2,FORMS)
	WRITE(PAR1,FORMS) Y
	I1 = INDEX(PAR1,'.')+1
	DO I=I1,LEN(PAR1)
	  IF(PAR1(I:I).NE.'0') GOTO 2
	ENDDO
2	I2 = I-I1+1

	CALL ROUNDF(F,I2,FORMS)
	WRITE(PAR1,FORMS) Y
	I1 = INDEX(PAR1,'.')+1
	DO I=I1,LEN(PAR1)
	  IF(PAR1(I:I).NE.'0') GOTO 3
	ENDDO
3     CONTINUE

4	WRITE(PAR,FORMS) X
	IF(Y.EQ.0) GOTO 5
	I2=LEN_TRIM(PAR)
	PAR(I2+1:I2+1)='('
	PAR(I2+3:I2+3)=')'
	PAR(I2+2:I2+2)=PAR1(I:I)
5     CONTINUE   
	DO I=1,LEN(PAR)
	  IF(PAR(I:I).NE.' ') GOTO 6
	ENDDO
6     S = I

	END

C ---------------------------------------------------------------
C ---------------------------------------------------------------
	SUBROUTINE ROUNDF(F,IROUND,FORM)
C ---------------------------------------------------------------
C ---------------------------------------------------------------
	CHARACTER*7 FORM1(5),FORM2(4),FORM3(4),FORM4(4),FORM5(3),FORM
	INTEGER F
	DATA FORM1 /'(F10.1)','(F10.2)','(F10.3)','(F10.4)','(F10.5)'/
	DATA FORM2 /'(F10.1)','(F10.2)','(F10.3)','(F10.4)'/
	DATA FORM3 /'(F12.1)','(F12.2)','(F12.3)','(F12.4)'/
	DATA FORM4 /'(F9.1)','(F9.2)','(F9.3)','(F9.4)'/
	DATA FORM5 /'(F8.1)','(F8.2)','(F8.3)'/

	SELECT CASE(F)
	  CASE(1)
	    FORM=FORM1(IROUND)
	  CASE(2)
	    FORM=FORM2(IROUND)
	  CASE(3)
	    FORM=FORM3(IROUND)
	  CASE(4)
	    FORM=FORM4(IROUND)
	  CASE(5)
	    FORM=FORM5(IROUND)
	END SELECT

	END


C--------------------------------------------------------------------
C--------------------------------------------------------------------
      SUBROUTINE BESD(X1,X2,X3,X4,X5,X6,S1,S2,S3,S4,
     &                S5,S6,B11,B22,B33,B12,B13,B23,
     &                S7,S8,S9,S10,S11,S12,BB,EBB)
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C     
      REAL DB(12),SX(12)
C
      PI=3.1415927
      SX(1)=S1
      SX(2)=S2
      SX(3)=S3
      SX(4)=S4
      SX(5)=S5
      SX(6)=S6
      SX(7)=S7
      SX(8)=S8
      SX(9)=S9
      SX(10)=S10
      SX(11)=S11
      SX(12)=S12
      DB(1)=2.0*B11*X1+B12*X2*COS(X6*PI/180.0)+B13*X3*COS(X5*PI/180.0)
      DB(2)=2.0*B22*X2+B12*X1*COS(X6*PI/180.0)+B23*X3*COS(X4*PI/180.0)
      DB(3)=2.0*B33*X3+B13*X1*COS(X5*PI/180.0)+B23*X2*COS(X4*PI/180.0)
      DB(4)=-B23*X2*X3*SIN(X4*PI/180.0)
      DB(5)=-B13*X1*X3*SIN(X5*PI/180.0)
      DB(6)=-B12*X1*X2*SIN(X6*PI/180.0)
      DB(7)=X1*X1
      DB(8)=X2*X2
      DB(9)=X3*X3
      DB(10)=X1*X2*COS(X6*PI/180.0)
      DB(11)=X1*X3*COS(X5*PI/180.0)
      DB(12)=X2*X3*COS(X4*PI/180.0)	  
      SESDB=0.0
      DO J=1,12
         SESDB=SESDB+(DB(J)*SX(J))**2
      END DO
      SESDB=4.0*SESDB/3.0
      EBB=SQRT(SESDB)
C
      RETURN 
      END


C ---------------------------------------------------------------
C ---------------------------------------------------------------

      SUBROUTINE  SIMPOS(F_ID,DB_ID,LAUE,INV,NPOSM,SGN,ID_TRANS)        
C ---------------------------------------------------------------
C ---------------------------------------------------------------
C      F_ID - OUTPUT FILE  
C      DB_ID - DATABASE FILE  
C      LAUE - LAUE GROUP NUMBER  
C      INV - INVERSION
C            INV=0 - NO INVERSION             
C            INV=1 - ADD INVERSION  
C      NPOSM - NUMBER OF XYZ POSITIONS
C      SGN - SPACE GROUP NUMBER
C      ID_TRANS - LATTICE TRANSLATIONS
C            ID_TRANS=0 - ASK ABOUT             
C            ID_TRANS=1 - ADD LATTICE TRANSLATIONS             
C            ID_TRANS=2 - NO LATTICE TRANSLATIONS             
C ---------------------------------------------------------------
      INTEGER ID_TRANS  
      INTEGER MAXSIM,MAXSIM2
      PARAMETER (MAXSIM=48,MAXSIM2=96)
      INTEGER F_ID, DB_ID, NPOSM, INV, LAUE
      CHARACTER*15 SGN
      CHARACTER*60 PAR,CMD1,CMD2,CMD3
	CHARACTER L(3)*10,PAR1*10,T(3)*5,C(2,3)*2,S(2,3)*1
	CHARACTER*5 TRANS(6,3,MAXSIM2),TR_B(6,3,MAXSIM2)
	INTEGER I,J,K,N,M
      CHARACTER*5 TRANS1(8),TRANS2(8),TRANS3(8),TRANS4(8),TRANS5(8)
C
* Modified by F. Izumi on 2001.02.21 {
* The following lines seems to be correct but causes errors in
* Absoft Pro Fortran
C     DATA TRANS1 /'     ','1/4  ','3/4  ','1/2  ','1/3  ','2/3  ',
C	&'1/6  ','5/6  '/
      DATA (TRANS1(J),J=1,8) /'     ','1/4  ','3/4  ','1/2  ','1/3  ',
     &  '2/3  ','1/6  ','5/6  '/
C [ -1 ]
C     DATA TRANS2 /'     ','3/4  ','1/4  ','1/2  ','2/3  ','1/3  ',
C	&'5/6  ','1/6  '/
      DATA (TRANS2(J),J=1,8) /'     ','3/4  ','1/4  ','1/2  ','2/3  ',
     &  '1/3  ','5/6  ','1/6  '/
C [ +1/2 ]				 
C     DATA TRANS3 /'1/2  ','3/4  ','1/4  ','     ','5/6  ','1/6  ',
C	&'2/3  ','1/3  '/
      DATA (TRANS3(J),J=1,8) /'1/2  ','3/4  ','1/4  ','     ','5/6  ',
     &  '1/6  ','2/3  ','1/3  '/
C [ +1/3 ]				 
C     DATA TRANS4 /'1/3  ','7/12 ','1/12 ','5/6  ','2/3  ','     ',
C	&'1/2  ','1/6  '/
      DATA (TRANS4(J),J=1,8) /'1/3  ','7/12 ','1/12 ','5/6  ','2/3  ',
     &  '     ','1/2  ','1/6  '/
C [ +2/3 ]
C      DATA TRANS5 /'2/3  ','11/12','5/12 ','1/6  ','     ','1/3  ',
C	&'5/6  ','1/2  '/
      DATA (TRANS5(J),J=1,8) /'2/3  ','11/12','5/12 ','1/6  ','     ',
     &  '1/3  ','5/6  ','1/2  '/
* }
C 

	DO I=1,NPOSM
	  DO J=1,3
	    L(J) = ''
		T(J) = ''
	    DO K=1,3
	      C(K,J) = '' 
            S(K,J) = ''
	    ENDDO
	  ENDDO
        PAR = '' 
        READ(DB_ID,'(A)') PAR

	  I1 = INDEX(PAR,',')
	  I2 = INDEX(PAR,',',BACK=.TRUE.)
	  L(1) = PAR(1:I1-1)
	  L(2) = PAR(I1+1:I2-1)
	  L(3) = PAR(I2+1:)
	  PAR1 = ''
	  DO J=1,3
	     PAR1 = L(J)
	     L(J) = ''
	     K = 1
	     DO M=1,LEN_TRIM(PAR1)
	        IF(PAR1(M:M).NE.' ') THEN 
	          L(J)(K:K) = PAR1(M:M)
	          K = K + 1
	        ENDIF
	     ENDDO
	  ENDDO
	  DO J=1,3
	     PAR1 = ''
	     I1 = INDEX(L(J),'/')
	     IF(I1.GT.0) THEN
		     T(J)=L(J)(I1-1:I1+1)
	         PAR1 = L(J)
	         IF(I1.EQ.2) THEN 
	            L(J) = PAR1(I1+2:)
	         ELSE 
	            L(J) = PAR1(1:I1-3)
	         ENDIF
	     ENDIF
	     N = 0
	     M = 0
	     DO K=1,LEN_TRIM(L(J))
	        IF((L(J)(K:K).EQ.'-') .OR. (L(J)(K:K).EQ.'+')) THEN
	          N = N + 1 
	          M = 0
	          S(N,J) = L(J)(K:K)
	        ELSE
	          IF (N .EQ. 0) N = 1 
	          M = M+1
                C(N,J)(M:M)=L(J)(K:K)
	        ENDIF
	     ENDDO
	     IF (S(1,J).EQ.'+') S(1,J) = ''
	  ENDDO
	  DO J=1,3
	     TRANS(1,J,I) = S(1,J)
	     TRANS(2,J,I) = C(1,J)
	     TRANS(3,J,I) = S(2,J)
	     TRANS(4,J,I) = C(2,J)
	     IF (T(J).NE.'') THEN
		      TRANS(5,J,I) = '+'
	     ELSE
		      TRANS(5,J,I) = ''
	     ENDIF
	     TRANS(6,J,I) = T(J)
	  ENDDO
C --------------------------------------------------------------------
C              CHECK INVERSION
C --------------------------------------------------------------------
        IF (INV.EQ.1) THEN
	    DO J=1,3
	      IF(S(1,J).EQ.'-') THEN
	 	    TRANS(1,J,I+NPOSM) = ''
	      ELSE
		    TRANS(1,J,I+NPOSM) = '-'
	      ENDIF
	      IF(S(2,J).EQ.'-') THEN
		  	TRANS(3,J,I+NPOSM) = '+'
	      ELSE
		 	IF(S(2,J).EQ.'+') THEN
			   TRANS(3,J,I+NPOSM) = '-'
	        ELSE
			   TRANS(3,J,I+NPOSM) = S(2,J)
	        ENDIF
	      ENDIF
            TRANS(6,J,I+NPOSM) = TRANS(6,J,I)
	      DO K=1,8
	        IF(INDEX(T(J),TRANS1(K)).GT.0) THEN
	          TRANS(6,J,I+NPOSM) = TRANS2(K)
	          GOTO 1
	        ENDIF
	      ENDDO
1	      TRANS(2,J,I+NPOSM) = TRANS(2,J,I)
	      TRANS(4,J,I+NPOSM) = TRANS(4,J,I)
	      TRANS(5,J,I+NPOSM) = TRANS(5,J,I)
	    ENDDO
	  ENDIF
	ENDDO
C --------------------------------------------------------------------
C			   WRITE DATA TO FILE
C --------------------------------------------------------------------
      IF (INV.EQ.1) THEN
        IDIM = 2*NPOSM
	ELSE
        IDIM = NPOSM
	ENDIF     
	DO I=1,IDIM
	  WRITE(CMD1,'(6A)')(TRANS(K,1,I)(1:LEN_TRIM(TRANS(K,1,I))),K=1,6) 
	  WRITE(CMD2,'(6A)')(TRANS(K,2,I)(1:LEN_TRIM(TRANS(K,2,I))),K=1,6) 
	  WRITE(CMD3,'(6A)')(TRANS(K,3,I)(1:LEN_TRIM(TRANS(K,3,I))),K=1,6) 
	  WRITE(F_ID,100) I,CMD1(1:LEN_TRIM(CMD1)),',',
     &                    CMD2(1:LEN_TRIM(CMD2)),',',
     &					CMD3(1:LEN_TRIM(CMD3))
	ENDDO
C
C ************************************************************************
C
C -------------------------------------------------------------------------
C                      ADD BRAVE TRANSLATIONS (F,I,A,B,C,R)
C -------------------------------------------------------------------------
      ID_SG=0
      IF (SGN(1:1).EQ.'F') ID_SG=1
      IF (SGN(1:1).EQ.'I') ID_SG=2
      IF (SGN(1:1).EQ.'A') ID_SG=3
      IF (SGN(1:1).EQ.'B') ID_SG=4
      IF (SGN(1:1).EQ.'C') ID_SG=5
      IF (SGN(1:1).EQ.'R') THEN
          IF (LAUE.EQ.9 .OR. LAUE.EQ.11) ID_SG=6
      END IF
      IF (ID_SG.EQ.0) GOTO 999 

	IF(ID_TRANS.EQ.0) THEN
	  WRITE(6,*) 'DO YOU WONT TO ADD LATTICE TRANSLATIONS? 
     & (1-YES/2-NO)'
	  READ(*,*) ID_TRANS
	ENDIF
      IF(ID_TRANS.NE.1) GOTO 999   


C -------------------------------------------------------------------------
C                      FOR (I,A,B,C)
C -------------------------------------------------------------------------
	DO I=1,IDIM
        DO J=1,3  
           DO K=1,6
             TR_B(K,J,I) = TRANS(K,J,I)
	     ENDDO
        ENDDO
        SELECT CASE (ID_SG)
C------------------------------------------- { I : 1/2,1/2,1/2 }
          CASE (2)
	      DO J=1,3
              DO K=1,8
	          IF(INDEX(TR_B(6,J,I),TRANS1(K)).GT.0) THEN
	            TR_B(6,J,I) = TRANS3(K)
	            GOTO 2
	          ENDIF
	        ENDDO
2	        IF(INDEX(TR_B(6,J,I),'/').GT.0) THEN
                     TR_B(5,J,I)='+' 
	        ELSE
                     TR_B(5,J,I)='' 
	        ENDIF
            ENDDO
C------------------------------------------- { A : 0,1/2,1/2 }
	    CASE (3)
	      DO J=2,3
              DO K=1,8
	          IF(INDEX(TR_B(6,J,I),TRANS1(K)).GT.0) THEN
	            TR_B(6,J,I) = TRANS3(K)
	            GOTO 3
	          ENDIF
	        ENDDO
3	        IF(INDEX(TR_B(6,J,I),'/').GT.0) THEN
                     TR_B(5,J,I)='+' 
	        ELSE
                     TR_B(5,J,I)='' 
	        ENDIF
            ENDDO
C------------------------------------------- { B : 1/2,0,1/2 }
	    CASE (4)
	      DO J=1,3,2
              DO K=1,8
	          IF(INDEX(TR_B(6,J,I),TRANS1(K)).GT.0) THEN
	            TR_B(6,J,I) = TRANS3(K)
	            GOTO 4
	          ENDIF
	        ENDDO
4	        IF(INDEX(TR_B(6,J,I),'/').GT.0) THEN
                     TR_B(5,J,I)='+' 
	        ELSE
                     TR_B(5,J,I)='' 
	        ENDIF
            ENDDO
C------------------------------------------- { C : 1/2,1/2,0 }
	    CASE (5)
	      DO J=1,2
              DO K=1,8
	          IF(INDEX(TR_B(6,J,I),TRANS1(K)).GT.0) THEN
	            TR_B(6,J,I) = TRANS3(K)
	            GOTO 5
	          ENDIF
	        ENDDO
5	        IF(INDEX(TR_B(6,J,I),'/').GT.0) THEN
                     TR_B(5,J,I)='+' 
	        ELSE
                     TR_B(5,J,I)='' 
	        ENDIF
            ENDDO

	  END SELECT
C --------------------------------------------------------------------
C			   WRITE DATA TO FILE
C --------------------------------------------------------------------
	  IF((ID_SG.GT.1).AND.(ID_SG.LT.6)) THEN
	    WRITE(CMD1,'(6A)')(TR_B(K,1,I)(1:LEN_TRIM(TR_B(K,1,I))),K=1,6)
	    WRITE(CMD2,'(6A)')(TR_B(K,2,I)(1:LEN_TRIM(TR_B(K,2,I))),K=1,6) 
	    WRITE(CMD3,'(6A)')(TR_B(K,3,I)(1:LEN_TRIM(TR_B(K,3,I))),K=1,6) 

	    WRITE(F_ID,100) (I+IDIM),CMD1(1:LEN_TRIM(CMD1)),',',
     &                       CMD2(1:LEN_TRIM(CMD2)),',',
     &   	                   CMD3(1:LEN_TRIM(CMD3))
	  ENDIF
	ENDDO
C -------------------------------------------------------------------------
C                      FOR (F)
C -------------------------------------------------------------------------
      IF(ID_SG.EQ.1) THEN
	DO M=1,3
	  DO I=1,IDIM
          DO J=1,3  
             DO K=1,6
               TR_B(K,J,I) = TRANS(K,J,I)
	       ENDDO
          ENDDO
          SELECT CASE (M)
C------------------------------------------- { F : 1/2,1/2,0 }
             CASE (1)
	         DO J=1,2
                 DO K=1,8
	             IF(INDEX(TR_B(6,J,I),TRANS1(K)).GT.0) THEN
	               TR_B(6,J,I) = TRANS3(K)
	               GOTO 6
	             ENDIF
	           ENDDO
6	           IF(INDEX(TR_B(6,J,I),'/').GT.0) THEN
                    TR_B(5,J,I)='+' 
	           ELSE
                    TR_B(5,J,I)='' 
	           ENDIF
               ENDDO
C------------------------------------------- { F : 1/2,0,1/2 }
             CASE (2)
	         DO J=1,3,2
                 DO K=1,8
	             IF(INDEX(TR_B(6,J,I),TRANS1(K)).GT.0) THEN
	               TR_B(6,J,I) = TRANS3(K)
	               GOTO 7
	             ENDIF
	           ENDDO
7	           IF(INDEX(TR_B(6,J,I),'/').GT.0) THEN
                    TR_B(5,J,I)='+' 
	           ELSE
                    TR_B(5,J,I)='' 
	           ENDIF
               ENDDO
C------------------------------------------- { F : 0,1/2,1/2 }
             CASE (3)
	         DO J=2,3
                 DO K=1,8
	             IF(INDEX(TR_B(6,J,I),TRANS1(K)).GT.0) THEN
	               TR_B(6,J,I) = TRANS3(K)
	               GOTO 8
	             ENDIF
	           ENDDO
8	           IF(INDEX(TR_B(6,J,I),'/').GT.0) THEN
                    TR_B(5,J,I)='+' 
	           ELSE
                    TR_B(5,J,I)='' 
	           ENDIF
               ENDDO
	    END SELECT
C --------------------------------------------------------------------
C			   WRITE DATA TO FILE
C --------------------------------------------------------------------
	    WRITE(CMD1,'(6A)')(TR_B(K,1,I)(1:LEN_TRIM(TR_B(K,1,I))),K=1,6)
	    WRITE(CMD2,'(6A)')(TR_B(K,2,I)(1:LEN_TRIM(TR_B(K,2,I))),K=1,6) 
	    WRITE(CMD3,'(6A)')(TR_B(K,3,I)(1:LEN_TRIM(TR_B(K,3,I))),K=1,6) 
	    WRITE(F_ID,100) (I+M*IDIM), CMD1(1:LEN_TRIM(CMD1)),',',
     &                       CMD2(1:LEN_TRIM(CMD2)),',',
     &	    			   CMD3(1:LEN_TRIM(CMD3))
 	  ENDDO
	ENDDO
	ENDIF
C -------------------------------------------------------------------------
C                      FOR (R: HEXAGONAL SETTING)
C -------------------------------------------------------------------------
      IF(ID_SG.EQ.6) THEN
	DO M=1,2
	  DO I=1,IDIM
          DO J=1,3  
             DO K=1,6
               TR_B(K,J,I) = TRANS(K,J,I)
	       ENDDO
          ENDDO
          SELECT CASE (M)
C------------------------------------------- { R : 2/3,1/3,1/3 }
             CASE (1)
               DO K=1,8
                 IF(INDEX(TR_B(6,1,I),TRANS1(K)).GT.0) THEN
	             TR_B(6,1,I) = TRANS5(K)
	             GOTO 9
	           ENDIF
	         ENDDO
9	         IF(INDEX(TR_B(6,1,I),'/').GT.0) THEN
                   TR_B(5,1,I)='+' 
	         ELSE
                   TR_B(5,1,I)='' 
	         ENDIF
			 DO J=2,3
                 DO K=1,8
	             IF(INDEX(TR_B(6,J,I),TRANS1(K)).GT.0) THEN
	               TR_B(6,J,I) = TRANS4(K)
	               GOTO 10
	             ENDIF
	           ENDDO
10	           IF(INDEX(TR_B(6,J,I),'/').GT.0) THEN
                    TR_B(5,J,I)='+' 
	           ELSE
                    TR_B(5,J,I)='' 
	           ENDIF
               ENDDO
C------------------------------------------- { R : 1/3,2/3,2/3 }
             CASE (2)
               DO K=1,8
                 IF(INDEX(TR_B(6,1,I),TRANS1(K)).GT.0) THEN
	             TR_B(6,1,I) = TRANS4(K)
	             GOTO 11
	           ENDIF
	         ENDDO
11	         IF(INDEX(TR_B(6,1,I),'/').GT.0) THEN
                   TR_B(5,1,I)='+' 
	         ELSE
                   TR_B(5,1,I)='' 
	         ENDIF
			 DO J=2,3
                 DO K=1,8
	             IF(INDEX(TR_B(6,J,I),TRANS1(K)).GT.0) THEN
	               TR_B(6,J,I) = TRANS5(K)
	               GOTO 12
	             ENDIF
	           ENDDO
12	           IF(INDEX(TR_B(6,J,I),'/').GT.0) THEN
                    TR_B(5,J,I)='+' 
	           ELSE
                    TR_B(5,J,I)='' 
	           ENDIF
               ENDDO
	    END SELECT
C --------------------------------------------------------------------
C			   WRITE DATA TO FILE
C --------------------------------------------------------------------
	    WRITE(CMD1,'(6A)')(TR_B(K,1,I)(1:LEN_TRIM(TR_B(K,1,I))),K=1,6)
	    WRITE(CMD2,'(6A)')(TR_B(K,2,I)(1:LEN_TRIM(TR_B(K,2,I))),K=1,6) 
	    WRITE(CMD3,'(6A)')(TR_B(K,3,I)(1:LEN_TRIM(TR_B(K,3,I))),K=1,6) 
	    WRITE(F_ID,100)(I+M*IDIM),CMD1(1:LEN_TRIM(CMD1)),',',
     &                      CMD2(1:LEN_TRIM(CMD2)),',',
     &	    		      CMD3(1:LEN_TRIM(CMD3))
 	  ENDDO
	ENDDO
	ENDIF

C
C ************************************************************************
C
999	IDIM=IDIM

100   FORMAT(3X,I3,3X,5A)


	END
