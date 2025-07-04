C *******************************************************************
	SUBROUTINE READLST(INP,OUT,IERR_ID)
C Read data from *.lst file -----------------------------------------
C
C *******************************************************************	
      INTEGER INP,OUT
	LOGICAL IERR_ID
* Corrected by F. Izumi on 2001.02.20 
C     CHARACTER CMD*300,ATOML*6
      CHARACTER CMD*300,ATOML*9
* Correction end
      CHARACTER PAR1*300, PAR2*300
      INTEGER(2) METHOD,ID
      REAL PARAM(11)
      INTEGER SCANL
	LOGICAL RESULT

      INTEGER NPH
      PARAMETER (NPH=8)
      INTEGER IATOM(NPH), ITATOM(NPH)
      CHARACTER*3 CONSTR(11)
      INTEGER CL(11) 
      DATA CONSTR /'g','x','y','z','B','B11','B22','B33','B12','B13',
     &'B23'/
      DATA CL /1,1,1,1,1,3,3,3,3,3,3/

C ------------------------------------------------------------------
C ---- CHECK OUTPUT MODE (NPRINT)
C  Current version of lst2cif work with NPRINT = 0 (Minimal output)
C ------------------------------------------------------------------
      REWIND(INP)
	CMD='NPRINT = '
	CALL FINDSTRING(0,INP,CMD,9,RESULT)
	IF(RESULT) THEN
	    GOTO 8888
	ELSE
	   WRITE(6,*) 'WRONG INPUT FILE'
         WRITE(6,*) ' NPRINT(OUTPUT MODE) > 0'
	   GOTO 9999
	ENDIF
8888	CONTINUE
C
C ------------------------------------------------------------------
C ---- CHECK RIETVELD  ANALYSIS METHOD
C ---- METHOD=0     TOF neutron diffraction
C ---- METHOD=1     X-ray and neutron diffraction
C ------------------------------------------------------------------
      WRITE(OUT,'(A)') 'GENERAL_BEGIN'
      J=-1
      I=-1 
      REWIND(INP)
	CMD='TOF neutron diffraction'
	CALL FINDSTRING(0,INP,CMD,23,RESULT)
	IF (RESULT) THEN
	   METHOD=0
	ELSE   
	  CMD='X-ray and neutron diffraction'
	  CALL FINDSTRING(0,INP,CMD,29,RESULT)
	  IF (RESULT) THEN
	     METHOD=1
	  ELSE
	     WRITE(6,*) 'WRONG INPUT FILE'
           WRITE(6,*) 'UNKNOWN RIETVELD ANALYSIS METHOD'
	     GOTO 9999
	  ENDIF
	ENDIF
      WRITE(OUT,'(I2)') METHOD                           ! METHOD
C ------------------------------------------------------------------
C ---- CHECK NBEAM      
C ---- NBEAM = 0 Neutron powder diffraction.
C ---- NBEAM = 1 Conventional X-ray powder diffraction with characteristic X rays.
C ---- NBEAM = 2 Synchrotron X-ray powder diffraction.
C ------------------------------------------------------------------
      IF (METHOD.EQ.1) THEN
	  CMD='NBEAM = '
	  CALL FINDSTRING(1,INP,CMD,8,RESULT)
	  IF (RESULT) THEN
 	     READ(INP,'(A)') CMD
	     I = INDEX(CMD,'NBEAM = ')
           WRITE(OUT,*) CMD(I+8:I+8)                    ! NBEAM
	  ELSE
	     WRITE(6,*) 'WRONG INPUT FILE'
		 WRITE(6,*) 'UNKNOWN BEAM TYPE'
	     GOTO 9999
	  ENDIF
C ----Wavelength
	  CMD='Wavelength'
	  CALL FINDSTRING(0,INP,CMD,10,RESULT)
	  IF (RESULT) THEN
	     READ(INP,'(A)') CMD
           WRITE(OUT,*) CMD(INDEX(CMD,'=')+2:INDEX(CMD,'=')+8)
	  ELSE
	     CMD='WAVELENGTH'
	     CALL FINDSTRING(0,INP,CMD,10,RESULT)
	     IF (RESULT) THEN
	        READ(INP,'(A)') CMD
              WRITE(OUT,*) CMD(INDEX(CMD,'=')+2:INDEX(CMD,'=')+8)
           ELSE
	       WRITE(6,*) 'WRONG INPUT FILE'
		 WRITE(6,*) 'UNKNOWN WAVELENGTH'
	       GOTO 9999
	     ENDIF
	  ENDIF
      END IF
C ------------------------------------------------------------------
C ---- CHECK RADIUS AND DENSITY
C ------------------------------------------------------------------
      CMD='RADIUS'
      CALL FINDSTRING(0,INP,CMD,6,RESULT)
      IF (RESULT) THEN
	  READ(INP,'(A)') CMD
	  WRITE(OUT,*) CMD(INDEX(CMD,'RADIUS')+8:INDEX(CMD,'RADIUS')+13)
	  WRITE(OUT,*) CMD(INDEX(CMD,'DENSTY')+8:INDEX(CMD,'DENSTY')+13)
	ELSE
	  WRITE(OUT,*) '0'
	  WRITE(OUT,*) '0'
      END IF
C ------------------------------------------------------------------
C ---- CHECK ATOMS
C ----No.  Atom    b/fm      SIGMAA/m**2   SIGMAI/m**2   CATOM/mol  At.wt.       ATOMNO/m**(-3)
C ------------------------------------------------------------------
C      CMD='Coefficients for analytic approximations'
      CMD='Atom'
      CALL FINDSTRING(0,INP,CMD,4,RESULT)
      IF (RESULT) THEN
	  READ(INP,'(A)') CMD
        CALL SKIPBLANK(INP,RESULT)
	  IF (RESULT) THEN
	     WRITE(6,*) 'WRONG INPUT FILE'
           WRITE(6,*) 'UNKNOWN TYPE OF ATOMS'
	     GOTO 9999
	  ENDIF
C	  READ(INP,'(A)') CMD
        ILINE=0
1	  READ(INP,'(A)') CMD
	  IF (LEN_TRIM(CMD).GT.0) THEN 
          J=SCANL(CMD,'0123456789',100,10)
          IF (J.GT.0) THEN
	       ILINE=ILINE+1 		
             BACKSPACE(INP)
C             READ(INP,*) ID,ATOML,(PARAM(I),I=1,11)
C             WRITE(OUT,100) ID,ATOML,(PARAM(I),I=1,11)
             READ(INP,*) ID,ATOML
             WRITE(OUT,100) ID,ATOML
             GOTO 1
	    ENDIF
	  ENDIF
      ENDIF
      IF (ILINE.EQ.0) THEN
	   WRITE(6,*) 'WRONG INPUT FILE'
	   WRITE(6,*) 'UNKNOWN NUMBER OF ATOMS'
	   GOTO 9999
	ENDIF
      WRITE(OUT,'(I3)') ILINE                 ! NUMBER OF ATOMS
      WRITE(OUT,'(A)') '###'
C ------------------------------------------------------------------
C IF METHOD=1 THEN: 2-THETA (MIN/MAX); STEP; NSTEP
C ------------------------------------------------------------------
      IF (METHOD.EQ.1) THEN 
         CMD='2-theta(min) ='
         CALL FINDSTRING(0,INP,CMD,14,RESULT)
	   IF (RESULT) THEN
           READ(INP,'(A)') CMD
           I1=INDEX(CMD,'2-theta(min) =')
           I2=INDEX(CMD,'2-theta(max) =')
           WRITE(OUT,'(A)') CMD(I1+15:I2-1)		! 2-THETA(MIN)
           WRITE(OUT,'(A)') CMD(I2+15:I2+15+10)	! 2-THETA(MAX)
         ELSE
           WRITE(OUT,'(A)') '0.0'
           WRITE(OUT,'(A)') '0.0'
         END IF
         CMD='STEP ='
         CALL FINDSTRING(1,INP,CMD,6,RESULT)
         IF (RESULT) THEN
            READ(INP,'(A)') CMD
            I1=INDEX(CMD,'STEP =')
            I2=INDEX(CMD,'NSTEP =')
            WRITE(OUT,'(A)') CMD(I1+7:I2-1)        ! STEP
            WRITE(OUT,'(A)') CMD(I2+8:I2+8+10)	 ! NSTEP
         ELSE
           WRITE(OUT,'(A)') '0.0'
           WRITE(OUT,'(A)') '0.0'
         END IF
      END IF
      WRITE(OUT,'(A)') 'GENERAL_END'
C ------------------------------------------------------------------
C ---- CHECK NUMBER OF PHASES
C ------------------------------------------------------------------
      I=0
      NPHASE=0
      CMD='Phase #'
 2    CALL FINDSTRING(I,INP,CMD,7,RESULT)
      IF (RESULT) THEN
         I=1
         NPHASE=NPHASE+1
	   CALL SKIPLINE(INP,1)
         GOTO 2
      ENDIF 
      IF (NPHASE.EQ.0) THEN
	   WRITE(6,*) 'WRONG INPUT FILE'
	   WRITE(6,*) 'UNKNOWN NUMBER OF PHASES'
	   GOTO 9999
	ELSE
         WRITE(OUT,'(A)') 'PHASE_BEGIN'
         WRITE(OUT,'(I3)') NPHASE                ! NUMBER OF PHASES      
	ENDIF
      I=0
      CMD='Phase #'
	DO II=1,NPHASE
        CALL FINDSTRING(I,INP,CMD,7,RESULT)
        IF (RESULT) THEN
           READ(INP,'(A)') PAR1
           I1=INDEX(PAR1,'Phase #')
           WRITE(OUT,'(A)') PAR1(I1:I1+9)    
           WRITE(OUT,'(A)') PAR1(I1+10:I1+90)       ! PAHSE NAME
           READ(INP,'(A)') PAR1         
           I1=INDEX(PAR1,'Space group:')
           I2=INDEX(PAR1,'(')
           I3=INDEX(PAR1,',')
           WRITE(OUT,'(A)') PAR1(I1+13:I2-1)        ! SPACE GROUP 
           WRITE(OUT,'(A)') PAR1(I2+6:I3-1)         ! VOLUME
           I4=INDEX(PAR1,'Setting #')
           IF (I4.GT.0) THEN
              PAR2(1:1)=PAR1(I4+9:I4+9)
              WRITE(OUT,'(3A)') PAR1(I3+1:I4-3),' ',PAR2(1:1)    ! SPACE GROUP NUMBER, SETTING
           ELSE
              I4=INDEX(PAR1,')')
              WRITE(OUT,'(2A)') PAR1(I2+8:I4-1),' 1'             ! SPACE GROUP NUMBER, SETTING   
           END IF
           READ(INP,'(A)') PAR1         
           I1=INDEX(PAR1,':')
           I2=INDEX(PAR1,',')
           WRITE(OUT,'(A)') PAR1(I1+2:I2-1)         ! SINGONY 
           WRITE(OUT,'(A)') PAR1(I2+2:I2+7)         ! POINT GROUP 
 	     CALL SKIPLINE(INP,1)
           I=1
        ELSE
	    WRITE(6,*) 'WRONG INPUT FILE'
	    WRITE(6,*) 'UNKNOWN SPACE GROUP FOR PHASE',II
	    GOTO 9999
        END IF
	ENDDO
      WRITE(OUT,'(A)') 'PHASE_END'
C ------------------------------------------------------------------
C ---- UNIT CELL PARAMETERS
C ------------------------------------------------------------------      
      CMD='Lattice parameters'
      WRITE(OUT,'(A)') 'CELL_BEGIN'
      DO I=1,NPHASE
        IF (I.EQ.1) THEN
           J=0
        ELSE
           J=1
        END IF 
        CALL FINDSTRING(J,INP,CMD,18,RESULT)
        IF (RESULT) THEN
	     CALL SKIPLINE(INP,1)
           CALL SKIPBLANK(INP,RESULT)
	     IF (RESULT) THEN
	        WRITE(6,*) 'WRONG INPUT FILE'
	        WRITE(6,*) 'UNKNOWN UNIT CELL PARAMETERS FOR PHASE',I
	        GOTO 9999
	     ELSE
		   WRITE(OUT,'(A,I2)') 'Phase #',I                    !Phase No
	       READ(INP,'(A)') PAR1 
             II=SCANL(PAR1,'abc',100,3)
		   IF (II.GT.0) THEN 
               READ(INP,*)  P_A,P_B,P_C,P_AL,P_BE,P_GA,P_VO
               WRITE(OUT,200) P_A,P_B,P_C,P_AL,P_BE,P_GA,P_VO     !a,b,c,alpha,beta,gamma,Volume
               READ(INP,'(A)') PAR1
 3	         K=INDEX(PAR1,'-')
               IF (K.GT.0) THEN
                  PAR1(K:K)='0'
                  GOTO 3
               END IF           
               READ(PAR1,*) SP_A,SP_B,SP_C,SP_AL,SP_BE,SP_GA,SP_VO       !ESD for a,b,c,alpha,beta,gamma,Volume
               WRITE(OUT,200) SP_A,SP_B,SP_C,SP_AL,SP_BE,SP_GA,SP_VO
	       ELSE
	         WRITE(6,*) 'WRONG INPUT FILE'
	         WRITE(6,*) 'UNKNOWN UNIT CELL PARAMETERS FOR PHASE',I
	         GOTO 9999
	       ENDIF
	     ENDIF
	  ELSE
	    WRITE(6,*) 'WRONG INPUT FILE'
	    WRITE(6,*) 'UNKNOWN UNIT CELL PARAMETERS FOR PHASE',I
	    GOTO 9999
        END IF 
      END DO
      WRITE(OUT,'(A)') 'CELL_END'
C ------------------------------------------------------------------
C ---- ATOMIC POSITIONS
C ------------------------------------------------------------------      
      CMD='Structure parameters'
      J=0
      DO I=1,NPHASE
        CALL FINDSTRING(J,INP,CMD,20,RESULT)
        IF (RESULT) THEN
           J=1
	     CALL SKIPLINE(INP,2)
           CALL SKIPBLANK(INP,RESULT)
	     IF (RESULT) THEN
	        WRITE(6,*) 'WRONG INPUT FILE'
	        WRITE(6,*) 'UNKNOWN STRUCTURE PARAMETERS FOR PHASE',I
	        GOTO 9999
	     ELSE
	       READ(INP,'(A)') PAR1
             II=SCANL(PAR1,'xyz',100,3)
	       IF (II.GT.0) THEN
                IATOM(I)=0
 4	          READ(INP,'(A)') PAR1
                IF (LEN_TRIM(PAR1).EQ.0) GOTO 44
 	          READ(INP,'(A)') PAR1
                IATOM(I)=IATOM(I)+1
                GOTO 4 
44              CONTINUE
	       ELSE
	          WRITE(6,*) 'WRONG INPUT FILE'
	          WRITE(6,*) 'UNKNOWN STRUCTURE PARAMETERS FOR PHASE',I
	          GOTO 9999
	       ENDIF
	     ENDIF
        ELSE
	    WRITE(6,*) 'WRONG INPUT FILE'
	    WRITE(6,*) 'UNKNOWN STRUCTURE PARAMETERS FOR PHASE',I
	    GOTO 9999
        END IF
      END DO
      J=0 
      WRITE(OUT,'(A)') 'ATOM_BEGIN'
      DO I=1,NPHASE 
         WRITE(OUT,'(A,I2)') 'Phase #',I              !Phase No
         WRITE(OUT,'(I3)') IATOM(I)                   !Number of atoms
         CALL FINDSTRING(J,INP,CMD,20,RESULT)
         READ(INP,'(A)') PAR1
* Corrected by F. Izumi on 2001.02.18 
C        I1=INDEX(PAR1,'Structure parameters')
         I1=INDEX(PAR1,'Structure parameters') + 2
* Correction end
         I2=I1+LEN(ATOML)
         J=1
	   CALL SKIPLINE(INP,1)
         CALL SKIPBLANK(INP,RESULT)
	   CALL SKIPLINE(INP,1)
         DO K=1,IATOM(I)
   	      READ(INP,'(A)') PAR1
   	      READ(INP,'(A)') PAR2
            DO L=I2+1,LEN(PAR1)
               I3=INDEX(PAR1(L:L),'-')
               IF ((I3.GT.0).AND.(PAR1(L+1:L+1).NE.'0')) PAR1(L:L)='0'
			 I3=INDEX(PAR2(L:L),'-')
               IF (I3.GT.0) PAR2(L:L)='0'
            END DO 
            READ(PAR1(I1:I2),'(A)') ATOML
            READ(PAR1(I2+1:),*) INEQ,G,SP,X,Y,Z,BEQ
            WRITE(OUT,300) ATOML,INEQ,G,X,Y,Z,BEQ           !Atom,g,x,y,z,Beq
            READ(PAR2(I2+1:),*) INEQ,SG,SP,SX,SY,SZ,SBEQ  
            WRITE(OUT,400) SG,SX,SY,SZ,SBEQ                 !ESD for g,x,y,z,Beq 
         END DO
      END DO 
      WRITE(OUT,'(A)') 'ATOM_END'
C ------------------------------------------------------------------
C ---- ATOMS TYPES
C ------------------------------------------------------------------      
      WRITE(OUT,'(A)') 'ATYPE_BEGIN'
      CMD='Number and weight of each species in the unit cell,'
      J=0
      DO I=1,NPHASE
        I1=0
        CALL FINDSTRING(J,INP,CMD,51,RESULT)
        IF (RESULT) THEN
           J=1
           CALL SKIPLINE(INP,1)
           CALL SKIPBLANK(INP,RESULT)
	     IF (RESULT) THEN
	        WRITE(6,*) 'WRONG INPUT FILE'
	        WRITE(6,*) 'UNKNOWN NUMBER OF SPECIES IN THE UNIT CELL'
	        GOTO 9999
	     ELSE
	        READ(INP,'(A)') PAR1
              II=SCANL(PAR1,'Atom',100,4)
	        IF (II.GT.0) THEN
  5	           READ(INP,'(A)') PAR1
                 IF (INDEX(PAR1,'---').LE.0) THEN
                    I1=I1+1
                    GOTO 5
	           ENDIF
	        ELSE
	           WRITE(6,*) 'WRONG INPUT FILE'
	           WRITE(6,*) 'UNKNOWN NUMBER OF SPECIES IN THE UNIT CELL'
	           GOTO 9999
              END IF
              ITATOM(I)=I1
	     ENDIF
	  ELSE
	    WRITE(6,*) 'WRONG INPUT FILE'
	    WRITE(6,*) 'UNKNOWN NUMBER OF SPECIES IN THE UNIT CELL'
	    GOTO 9999
        ENDIF
      END DO  
      J=0
      DO I=1,NPHASE
        CALL FINDSTRING(J,INP,CMD,51,RESULT)
        IF (RESULT) THEN
           WRITE(OUT,'(A,I2)') 'Phase #',I                     !Phase No
           WRITE(OUT,*) ITATOM(I)                              ! Number of Atoms by Type
           J=1
           CALL SKIPLINE(INP,1)
           CALL SKIPBLANK(INP,RESULT)
           READ(INP,'(A)') PAR1
           I1=INDEX(PAR1,'Atom')
  6	     READ(INP,'(A)') PAR1
           IF (INDEX(PAR1,'---').LE.0) THEN
              WRITE(OUT,'(2A)') PAR1(I1:I1+2),PAR1(I1+3:I1+3+15)    !Atom Type ; Number of Atoms in Unit Cell
              GOTO 6
           END IF
        END IF
      END DO  
      WRITE(OUT,'(A)') 'ATYPE_END'
C ------------------------------------------------------------------
C ---- THERMAL PARAMETERS
C ------------------------------------------------------------------      
	CMD='beta11  beta22  beta33'
      WRITE(OUT,'(A)') 'THERMAL_BEGIN'
      J=0
      DO I=1,NPHASE
        WRITE(OUT,'(A,I2)') 'Phase #',I                     !Phase No
        WRITE(OUT,'(I3)') IATOM(I)                          !Number of atoms 
        CALL FINDSTRING(J,INP,CMD,22,RESULT)
        IF (RESULT) THEN
	      J=1
		  CALL SKIPLINE(INP,1)
            DO K=1,IATOM(I)
	        CALL SKIPBLANK(INP,RESULT)
	        IF (RESULT) THEN
	          WRITE(6,*) 'WRONG INPUT FILE'
	          WRITE(6,*) 'UNKNOWN THERMAL PARAMETERS FOR ATOM ',K,
     &		     	     ' FOR PHASE ', I
	          GOTO 9999
	        ENDIF
   	        READ(INP,'(A)') PAR1
	        CALL SKIPBLANK(INP,RESULT)
	        IF (RESULT) THEN
	          WRITE(6,*) 'WRONG INPUT FILE'
	          WRITE(6,*) 'UNKNOWN ESD VALUE FOR ATOM ',K,
     &		     	     ' FOR PHASE ', I
	          GOTO 9999
	        ENDIF
		    READ(INP,'(A)') PAR2
	        CALL SKIPBLANKCH(PAR1,200,I1)
	        I2=I1+LEN(ATOML)
              DO L=I2+1,LEN(PAR1)
                 I3=INDEX(PAR1(L:L),'-')
                 IF (I3.GT.0) PAR1(L:L)='0'
                 I3=INDEX(PAR2(L:L),'-')
                 IF (I3.GT.0) PAR2(L:L)='0'
              END DO 
              READ(PAR1(I1:I2),'(A)') ATOML
              READ(PAR1(I2+1:),*) B11,B22,B33,B12,B13,B23,
     &                            U11,U22,U33,U12,U13,U23,BEQ
              READ(PAR2(I2+1:),*) SB11,SB22,SB33,SB12,SB13,SB23,
     &                             SU11,SU22,SU33,SU12,SU13,SU23,SBEQ
              WRITE(OUT,500) ATOML,                                    !B11,B22,B33,B12,B13,B23 and
     &                       B11/100000.0,B22/100000.0,B33/100000.0,   !U11,U22,U33,U12,U13,U23 and
     &                       B12/100000.0,B13/100000.0,B23/100000.0,   !Beq
     &                       U11/100000.0,U22/100000.0,U33/100000.0,
     &                       U12/100000.0,U13/100000.0,U23/100000.0,
     &                       BEQ
              WRITE(OUT,600)                                           !ESD for
     &                      SB11/100000.0,SB22/100000.0,SB33/100000.0,  !B11,B22,B33,B12,B13,B23 and
     &                      SB12/100000.0,SB13/100000.0,SB23/100000.0,  !U11,U22,U33,U12,U13,U23 and
     &                      SU11/100000.0,SU22/100000.0,SU33/100000.0,  !Beq 
     &                      SU12/100000.0,SU13/100000.0,SU23/100000.0,
     &                      SBEQ
            END DO
	  ELSE	    
	    WRITE(6,*) 'WRONG INPUT FILE'
	    WRITE(6,*) 'UNKNOWN THERMAL PARAMETERS FOR PHASE',I
	    GOTO 9999
        END IF
      END DO
      WRITE(OUT,'(A)') 'THERMAL_END'
C ------------------------------------------------------------------
C ---- R - FACTORS
C ------------------------------------------------------------------      
      WRITE(OUT,'(A)') 'RFACTOR_BEGIN'
      CMD='Rwp ='
      CALL FINDSTRING(0,INP,CMD,5,RESULT)
	 IF (RESULT) THEN
	   READ(INP,'(A)') PAR1
	   IF ((INDEX(PAR1,'Rwp').GT.0).AND.(INDEX(PAR1,'Rp').GT.0).AND.
     &       (INDEX(PAR1,'RR').GT.0).AND.(INDEX(PAR1,'Re').GT.0)) THEN
            I2=INDEX(PAR1,'=')
            READ(PAR1(I2+1:),*) RWP             
            PAR1=PAR1(I2+1:)
            I2=INDEX(PAR1,'=')
            READ(PAR1(I2+1:),*) RP
            PAR1=PAR1(I2+1:)
            I2=INDEX(PAR1,'=')
            READ(PAR1(I2+1:),*) RR
            PAR1=PAR1(I2+1:)
            I2=INDEX(PAR1,'=')
            READ(PAR1(I2+1:),*) RE
            WRITE(OUT,'(4F6.2)') RWP,RP,RR,RE          !Rwp,Rp,RR,Re
            DO I=1,NPHASE
              WRITE(OUT,'(A,I2)') 'Phase #',I         !Phase No
	        CALL SKIPBLANK(INP,RESULT)
	        IF (RESULT) THEN
	           WRITE(6,*) 'WRONG INPUT FILE'
	           WRITE(6,*) 'UNKNOWN R-FACTORS'
	           GOTO 9999
	        ENDIF
              CALL SKIPLINE(INP,1)
              READ(INP,'(A)') PAR1
              I2=INDEX(PAR1,'=')
              READ(PAR1(I2+1:),*) RI             
              PAR1=PAR1(I2+1:)
              I2=INDEX(PAR1,'=')
              READ(PAR1(I2+1:),*) RF
              WRITE(OUT,'(2F6.2)') RI,RF              !Ri,Rf
	      ENDDO
	   ELSE
	      WRITE(6,*) 'WRONG INPUT FILE'
	      WRITE(6,*) 'UNKNOWN R-FACTORS'
	      GOTO 9999
	   ENDIF
	ENDIF
	WRITE(OUT,'(A)') 'RFACTOR_END'


C ------------------------------------------------------------------
C ---- LINEAR CONSTRAINTS
C ------------------------------------------------------------------      
      CMD='Linear constraints'
      NCONSTR=0
      PAR2=''
      CALL FINDSTRING(0,INP,CMD,18,RESULT)
      IF (RESULT) THEN
         CALL SKIPLINE(INP,1)
 10	   READ(INP,'(A)') PAR1  
         I1=INDEX(PAR1,'A(') 
         IF (I1.GT.0) THEN 
            J=1     
            DO I=1,LEN(PAR1)
               IF (PAR1(I:I).NE.' ' .AND. PAR1(I:I).NE.';') THEN
                     PAR2(J:J)=PAR1(I:I)
                     J=J+1
               ELSE
                  IF (PAR1(I:I).EQ.';' .OR. I.EQ.LEN(PAR1)) THEN
                     NCONSTR=NCONSTR+1 
                     J=1
                     PAR2=''    
                  END IF
               END IF
            END DO
            GOTO 10
         END IF
      END IF 
      WRITE(OUT,'(A)') 'CONSTR_BEGIN'
      WRITE(OUT,*) NCONSTR                                                ! Number of Constraints
      IF (NCONSTR.GT.0) THEN   
         CMD='Linear constraints'
         CALL FINDSTRING(0,INP,CMD,18,RESULT)
         CALL SKIPLINE(INP,1)
 11	   PAR1=''
         READ(INP,'(A)') PAR1  
* Corrected by R. Dilanian on 2001.04.13 { 
         I1=INDEX(PAR1,'A(') 
         IF (I1.GT.0) THEN 
            J=1     
            DO I=1,LEN(PAR1)
               IF (PAR1(I:I).NE.' ' .AND. PAR1(I:I).NE.';') THEN
                     PAR2(J:J)=PAR1(I:I)
                     J=J+1
               ELSE
                  IF (PAR1(I:I).EQ.';' .OR. I.EQ.LEN(PAR1)) THEN
                     WRITE(OUT,'(A)') PAR2(1:J)         !Constraint
                     J=1
                     PAR2=''    
                  END IF
               END IF
            END DO
            GOTO 11
         END IF
*}
      END IF 
      WRITE(OUT,'(A)') 'CONSTR_END'  
	IERR_ID=.TRUE.
	GOTO 13
9999	IERR_ID=.FALSE.
  13  CONTINUE
C ----------------------FINISH---------------------------------------
C 100  FORMAT(I2,3X,A5,9F11.6,2F8.3)
 100  FORMAT(I2,3X,A5)
 200  FORMAT(3F10.5,3F10.4,F12.4)
 300  FORMAT(A,I3,F9.4,3F10.5,F8.3)	
* Corrected by F. Izumi on 2001.02.20 
C400  FORMAT(9X,F9.4,3F10.5,F8.3)	
 400  FORMAT(12X,F9.4,3F10.5,F8.3)	
* Correction end
 500  FORMAT(A,12F10.5,F8.3)	
* Corrected by F. Izumi on 2001.02.20 
C600  FORMAT(6X,12F10.5,F8.3)	
 600  FORMAT(9X,12F10.5,F8.3)	
* Correction end

C 100  FORMAT(I2,3X,A5,0P,F8.3,2(5X,1PE9.3),4X,0P,F7.3,F13.4,5X,1PE10.4)

* Corrected by F. Izumi on 2001.02.25
*     *.lst will be read in block4.f90 to input label/species
C     CLOSE(UNIT=INP)
      CLOSE(UNIT=OUT)

	END
