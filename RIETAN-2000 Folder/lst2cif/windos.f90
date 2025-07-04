*     Files used in this program

* Corrected by F. Izumi on 2001.04.14 {
C     Unit No.    File		                  var
C     1           *.lst (not "rietan output")	  INP
C     2           *.tmp (scratch file)		  TMP
C     3           *.cif                           OUT
C     4           spgri or spgra
* }

*     An environment variable DATABASE (absolute path for database
*     files) has to be set in the batch file.

C *******************************************************************
	PROGRAM LST2CIF

C CREATING CIF FILE FROM RIETAN'S OUTPUT FILES:
C     *.lst  - OUTPUT FROM RIETAN-2000

*     For Fortran Power Station
C	USE MSFLIB
*     For Visual Fortran
      USE DFLIB

C
C *******************************************************************	
      CHARACTER*80 INPFNAME,OUTFNAME,TEMPFNAME
      INTEGER INP,TMP,OUT
C	INTEGER*2 N1,N2,STATUS1,STATUS2
	LOGICAL RESULT


	WRITE(6,*) '**************************************************'
	WRITE(6,*) '**                   lst2cif                    **'
      WRITE(6,*) '**   A FORTRAN PROGRAM FOR CONVERTION OF THE    **'
	WRITE(6,*) '**     RIETVELD ANALYSIS OUTPUT (*.lst) FILE    **'
	WRITE(6,*) '**            TO THE CIF FORMAT FILE            **'
	WRITE(6,*) '**                                              **'
      WRITE(6,*) '**     COPYLEFT (C) 2001, BY R.A. DILANIAN      **'
      WRITE(6,*) '**              ALL RIGHTS RELEASED             **'
      WRITE(6,*) '**                                              **'
	WRITE(6,*) '**************************************************'

* Corrected by F. Izumi on 2001.04.14 {
*      ID_TRANS - LATTICE TRANSLATIONS
*            ID_TRANS=0 - ASK ABOUT             
*            ID_TRANS=1 - ADD LATTICE TRANSLATIONS             
*            ID_TRANS=2 - NO LATTICE TRANSLATIONS             
* Command line to run LST2CIF:
*      LST2CIF [INP FILE] [OUT FILE] [ID_TRANS]
* or
*      LST2CIF 
*     ID_TRANS is fixed at 1 (CrystalMaker is not used on Windows)
      ID_TRANS = 1                             
* }

C	N1=1
C	N2=2
C	CALL GETARG(N1,INPFNAME,STATUS1)
C	CALL GETARG(N2,OUTFNAME,STATUS2)
C	IF (STATUS1.GT.0) THEN
C        IF (STATUS2.LT.0) THEN
C	      OUTFNAME=INPFNAME(1:INDEX(INPFNAME,'.',BACK=.TRUE.))//'cif' 
C	   ENDIF
C	ELSE
C        WRITE(6,*) 'INPUT FILE NAME'
C	   READ(*,*) INPFNAME
C	   WRITE(6,*) 'OUTPUT FILE NAME'
C	   READ(*,*) OUTFNAME
C	ENDIF

      CALL GETARG(1,INPFNAME)
* Modified by F. Izumi on 2001.02.19 
      TEMPFNAME=INPFNAME(1:INDEX(INPFNAME,'.'))//'tmp'
	OUTFNAME=INPFNAME(1:INDEX(INPFNAME,'.'))//'cif'
* Modification end

      CALL OPENLST(INP,TMP,INPFNAME,TEMPFNAME)
      CALL READLST(INP,TMP,RESULT)
	IF (RESULT) THEN
*     Reopen the temporary file
        OPEN(UNIT=TMP,FILE=TEMPFNAME,STATUS='OLD')
        CALL OPENCIF(OUT,OUTFNAME)
* Modified by F. Izumi on 2001.04.14 
        CALL WRITECIF(INP,OUT,TMP,ID_TRANS) 
        CLOSE(UNIT=OUT)
* Modified by F. Izumi on 2001.02.19 
C       CLOSE(UNIT=TMP)
        CLOSE(UNIT=TMP,STATUS='DELETE')
* Modification end
	ENDIF

      END

C *******************************************************************
C     OPEN *.LST FILE AND SCRATCH FILE
C *******************************************************************
      SUBROUTINE OPENLST(F_ID1,F_ID2,FNAME1,FNAME2)
* Modified by F. Izumi on 2001.02.18 
C*	USE MSFLIB
C     USE DFLIB
C
C     CHARACTER(80) FNAME1,FNAME2
      CHARACTER*80 FNAME1,FNAME2
* Modification end

      INTEGER F_ID1,F_ID2
C     INTEGER(2) N1

      F_ID1=1
      F_ID2=2
C 	N1=1
C     CALL GETARG(N1,FNAME1)
C     FNAME2='tmp.out'
      OPEN(UNIT=F_ID1,FILE=FNAME1,STATUS='OLD')
      OPEN(UNIT=F_ID2,FILE=FNAME2,STATUS='UNKNOWN')

      END

C *******************************************************************
C     OPEN *.CIF FILE
C *******************************************************************
      SUBROUTINE OPENCIF(F_ID,FNAME)
* Modified by F. Izumi on 2001.02.18 
C*	USE MSFLIB
C     USE DFLIB
C
C     CHARACTER(80) FNAME
      CHARACTER*80 FNAME
* Modification end

      INTEGER F_ID
C     INTEGER(2) N2

      F_ID=3
C 	N2=2
C     CALL GETARG(N2,FNAME)
      OPEN(UNIT=F_ID,FILE=FNAME,STATUS='UNKNOWN')
C
      END

C *******************************************************************
C     OPEN DATABASE FILE (SPGRA or SPGRI)
C *******************************************************************
      SUBROUTINE OPENDATABASE(D_ID,F_ID)
* Modified by F. Izumi on 2001.02.18 
C*    USE MSFLIB
      USE DFPORT
* Modification end
C
      INTEGER F_ID
* Modified by F. Izumi on 2001.02.18 
* A too short environment variable
C     CHARACTER D_ID*1,EVALUE*50 
      CHARACTER D_ID*1,EVALUE*250
* Modification end

      F_ID=4
      CALL GETENV('DATABASE',EVALUE)
* Corrected by F. Izumi on 2001.02.18 
* Environment variable DATABASE may contain a space
      DO J = 250, 1, -1
	   IF (EVALUE(J:J) .NE. ' ') EXIT
      END DO
      IF (D_ID .EQ. 'I') THEN
C        OPEN(UNIT=F_ID,FILE=EVALUE(1:INDEX(EVALUE,' ')-1)//'/spgri',
C    &   STATUS='OLD')
         OPEN(UNIT=F_ID,FILE=EVALUE(:J)//'\spgri',STATUS='OLD')
      ELSE
C        OPEN(UNIT=F_ID,FILE=EVALUE(1:INDEX(EVALUE,' ')-1)//'/spgra',
C    &   STATUS='OLD')
        OPEN(UNIT=F_ID,FILE=EVALUE(:J)//'\spgra',STATUS='OLD')
* Correction end
      END IF
      END
