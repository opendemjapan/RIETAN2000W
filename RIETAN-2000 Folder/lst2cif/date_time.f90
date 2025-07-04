      SUBROUTINE GET_DATE_TIME(DATE_TIME)
*     Refer to 9-46 in Compaq Fortran - Launguage Reference Manual
      CHARACTER CDATE*8, CTIME*10
      CHARACTER DATE_TIME*16

      CALL DATE_AND_TIME(CDATE,CTIME)
* Corrected by F. Izumi on 2001.02.21 
C     WRITE(DATE_TIME,100) CDATE(1:4),CDATE(5:6),CDATE(7:8),
C    &                     CTIME(1:2),CTIME(3:4)
C
C100  FORMAT(A,'-',A,'-',A,'T',A,':',A)
      WRITE(DATE_TIME,100) CDATE(1:4),CDATE(5:6),CDATE(7:8),CTIME(1:2),
     &  CTIME(3:4)
  100 FORMAT(A,'-',A,'-',A,' ',A,':',A)
* Correction end
      END 
