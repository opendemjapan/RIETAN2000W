      SUBROUTINE JOBEND(MESSAG)
C     PRINT OUT THE CPU TIME AND END OF JOB MESSAGE
      REAL TARRAY(2)
      CHARACTER*(*) MESSAG

      IF(MESSAG(1:1).NE.' ') WRITE(6,'(//11X,A)') MESSAG
      WRITE(6,'(///7X,A//)') '*** End of job ***'
      SECONDS = ELAPTIME(TARRAY)
      WRITE(6,200) NINT(SECONDS)/60, MOD(SECONDS,60.0)
  200 FORMAT(' ',10X,'Elapsed time =',I4,' min',F6.1,' s'//
     &' ',10X,'--- RIETAN-2000.',
*Rev 1.1 2003.05.16 Izumi
     &'  Copyleft 2000-2003 by F. Izumi ---'//)
      STOP
      END
