      SUBROUTINE INITIAL
*     OPEN UNIT 5 AND PRINT DATE
*     Cf. Language Reference Manual of Digital Fortran, E-16.
	  USE DFLIB
      CHARACTER FNAME*50
      COMMON /MATER/ FNAME
      
      CALL GETARG(1,FNAME)

* Modified by F. Izumi on 2001.04.14 {
*     Check the extension
      DO J = 50, 1, -1
         IF (FNAME(J:J) .NE. ' ') EXIT
      END DO
      IF (J .LT. 5) STOP 'Too short file name'
      IF (FNAME(J-3:J) .NE. '.xyz') 
     &STOP 'The extension in the file name should be ''xyz'''
* }
      OPEN(UNIT=9,FILE=FNAME,STATUS='OLD')
      END

************************************************************************

      SUBROUTINE OPEN4(LEXIST)
      CHARACTER FNAME*50
      LOGICAL LEXIST
      COMMON /MATER/ FNAME

      IDOT = INDEX(FNAME,'.')
      INQUIRE(FILE=FNAME(1:IDOT)//'ffe',EXIST=LEXIST)
      IF (.NOT. LEXIST) THEN
         OPEN(4,FILE=FNAME(1:IDOT)//'ffe',STATUS='NEW')
         REWIND 4
      END IF
      END
