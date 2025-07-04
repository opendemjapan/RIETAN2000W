C *******************************************************************
C     WRITE DATA TO CIF FILE
C *******************************************************************
* Corrected by F. Izumi on 2001.02.25 
C     SUBROUTINE WRITECIF(OUT,TMP) 
C     SUBROUTINE WRITECIF(INP,OUT,TMP) 
* Corrected by R. Dilanian on 2001.03.01 
      SUBROUTINE WRITECIF(INP,OUT,TMP,ID_TRANS) 
C
      INTEGER INP,TMP,OUT,ID_TRANS
* Correction end
C
C--------------------------------------------------------------------
      CALL WRITE_BLOCK0(OUT)        ! MAIN BLOCK
	WRITE(6,*) 'BLOCK0 WRITE --------------------- OK'
      CALL WRITE_BLOCK1(OUT)        ! SUBMISSION DETAILS   
	WRITE(6,*) 'BLOCK1 WRITE --------------------- OK'
      CALL WRITE_BLOCK2(OUT)        ! PROCESSING SUMMARY
	WRITE(6,*) 'BLOCK2 WRITE --------------------- OK'
      CALL WRITE_BLOCK3(OUT)        ! TITLE, AUTHOR LIST AND TEXT
	WRITE(6,*) 'BLOCK3 WRITE --------------------- OK'

* Corrected by F. Izumi on 2001.02.25 
C     CALL WRITE_BLOCK4(OUT,TMP)    ! CRYSTAL DATA
C      CALL WRITE_BLOCK4(INP,OUT,TMP)    ! CRYSTAL DATA
* Corrected by R. Dilanian on 2001.03.01 
      CALL WRITE_BLOCK4(INP,OUT,TMP,ID_TRANS)    ! CRYSTAL DATA
	WRITE(6,*) 'BLOCK4 WRITE --------------------- OK'
* Correction end
      CALL WRITE_BLOCK5(OUT,TMP)    ! POWDER SPECIMEN AND EXPERIMENTAL DATA
	WRITE(6,*) 'BLOCK5 WRITE --------------------- OK'
      CALL WRITE_BLOCK6(OUT,TMP)    ! REFINEMENT DATA
	WRITE(6,*) 'BLOCK6 WRITE --------------------- OK'

      WRITE(OUT,900)
C     WRITE(OUT,900)
      WRITE(OUT,100)
C--------------------------------------------------------------------
C
 100  FORMAT('#',15('--eof'),'--#')
 900  FORMAT(' ')
 901  FORMAT(';')
 902  FORMAT('#',70('='))
 903  FORMAT('loop_')
C
      END
