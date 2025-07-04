C **************************************************************************
	SUBROUTINE SKIPBLANK(INP,END_ID)
C --- SKIP BLANK LINES IN THE FILE
C --- INP - FILE UNIT NUMBER
C **************************************************************************
	INTEGER INP
	CHARACTER CMD*300
	LOGICAL END_ID 

	END_ID=.FALSE.
1	CMD=''
      READ(INP,'(A)',END=9) CMD
	I=LEN_TRIM(CMD)
	IF (I.EQ.0) GOTO 1
      BACKSPACE(INP)
	GOTO 2
9     END_ID=.TRUE.
2     CONTINUE

      END

C **************************************************************************
      SUBROUTINE SKIPLINE(INP,LCOUNT)
C --- SKIP LINES IN THE FILE
C --- INP - FILE UNIT NUMBER
C --- LCOUNT - NUMBER OF LINES
C **************************************************************************
      INTEGER INP,LCOUNT

      DO I=1,LCOUNT
        READ(INP,*,END=99)
      ENDDO
      
	
 99   RETURN
      END

C **************************************************************************
      SUBROUTINE FINDSTRING(REW_ID,INP,STRING,STRINGL,RESULT)
C --- FIND STRING IN THE FILE
C --- STRINGL - LENGTH OF STRING
C --- INP - FILE UNIT NUMBER
C --- IF REW_ID=0 THEN REWIND FILE
C ---   			ELSE CONTINUE
C **************************************************************************
      CHARACTER STRING*300
      INTEGER STRINGL, REW_ID
	LOGICAL RESULT
	CHARACTER PAR*300	

C
      RESULT=.FALSE.
      IF (REW_ID.EQ.0) REWIND(INP)
 1    READ(INP,'(A)',END=999) PAR
      I=INDEX(PAR,STRING(1:STRINGL))
      IF (I.GT.0) THEN
	 BACKSPACE(INP)
	 RESULT=.TRUE.
	 GOTO 3
      END IF
      GOTO 1
C
999   CONTINUE
3     RETURN
      END


C **************************************************************************
      INTEGER FUNCTION SCANL(STRING,SUBSTRING,LSTRING,LSUBSTRING)
C --- SCAN STRING AND FINDE CHARACTER FROM SUBSTRING
C --- LSTRING - LENGTH OF STRING
C --- LSUBSTRING -  LENGTH OF SUBSTRING
C **************************************************************************
      CHARACTER*(*) STRING, SUBSTRING
      INTEGER LSTRING,LSUBSTRING      

      SCANL=0
      DO I=1,LSTRING
         DO J=1,LSUBSTRING
            IF(STRING(I:I).EQ.SUBSTRING(J:J)) THEN
		    SCANL=I
			GOTO 9999
	      ENDIF
         END DO
      END DO

9999  RETURN
      END  
        

C **************************************************************************
      SUBROUTINE SKIPBLANKCH(STRING, LSTRING,N)
C --- SKIP ' 'CHARACTER IN THE STRING
C --- LSTRING - LENGTH OF STRING 
C **************************************************************************
      CHARACTER*(*) STRING
      INTEGER LSTRING,N

      N=0
	K=0
	DO I=1,LSTRING
	 K=K+1 
	 J=LEN_TRIM(STRING(I:I))
	IF (J.GT.0) GOTO 4
	ENDDO
	GOTO 5
4     N=K
5	RETURN
	END	


C **********************************************************************
	SUBROUTINE ENDPROG(MSG)
C **************************************************************************
	CHARACTER*(*) MSG

	WRITE(6,'(5X,A)') MSG

	STOP
	END