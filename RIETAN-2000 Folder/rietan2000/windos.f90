      SUBROUTINE INITIAL
*     OPEN UNIT 5 AND PRINT DATE
*     Cf. Language Reference Manual of Compaq Fortran, E-16.
      USE DFLIB
      CHARACTER FILE5*50,CYMD*8,HMS*10

      CALL GETARG(1,FILE5)
      OPEN(UNIT=5,FILE=FILE5,STATUS='OLD')

*     Example: 08:31:41, 16-FEB-1999
*     Get the date and time
*     Cf. Language Reference Manual of Compaq Fortran, 9-46.
*     CYMD: CCYYMMDD
*     HMS: hhmmss.sss
      CALL DATE_AND_TIME(CYMD,HMS)
*     Example: 08:31:41, 1999.02.18
      WRITE(6,'(//1H ,85X,11A)') HMS(1:2),':',HMS(3:4),':',HMS(5:6),
     &  ', ',CYMD(1:4),'.',CYMD(5:6),'.',CYMD(7:8)
      END

************************************************************************

      REAL FUNCTION ELAPTIME(TARRAY)
*     Cf. Language Reference Manual of Compaq Fortran, E-1.
      USE DFPORT
      REAL TARRAY(2)

      ELAPTIME = ETIME(TARRAY)
      END

************************************************************************

      SUBROUTINE OPENSPGR(VOLUME,NFILE)
*     OPEN A SPACE GROUP FILE
*     Cf. Language Reference Manual of Compaq Fortran, E-1.
      USE DFPORT
      CHARACTER VOLUME*1, EVALUE*50

      CALL GETENV('RIETAN',EVALUE)
*Rev 1.0n1 2001.01.15 Izumi {
      DO J = 50, 1, -1
	   IF (EVALUE(J:J) .NE. ' ') EXIT
      END DO
      IF (VOLUME .EQ. 'I') THEN
         NFILE=1
*        OPEN(UNIT=NFILE,FILE=EVALUE(1:INDEX(EVALUE,' ')-1)//'\spgri',
*    &   STATUS='OLD')
         OPEN(UNIT=NFILE,FILE=EVALUE(1:J)//'\spgri',STATUS='OLD')
      ELSE
         NFILE=11
*        OPEN(UNIT=NFILE,FILE=EVALUE(1:INDEX(EVALUE,' ')-1)//'\spgra',
*    &   STATUS='OLD')
         OPEN(UNIT=NFILE,FILE=EVALUE(1:J)//'\spgra',STATUS='OLD')
      END IF
* }
      END

************************************************************************

      SUBROUTINE OPENASFDC
*     OPEN FILE asfdc
*     Cf. Language Reference Manual of Compaq Fortran, E-1.
      USE DFPORT
      CHARACTER EVALUE*50

      CALL GETENV('RIETAN',EVALUE)
*Rev 1.0n1 2001.01.15 Izumi {
      DO J = 50, 1, -1
	   IF (EVALUE(J:J) .NE. ' ') EXIT
      END DO
*     OPEN(UNIT=2,FILE=EVALUE(1:INDEX(EVALUE,' ')-1)//'\asfdc',
*    &  STATUS='OLD')
      OPEN(UNIT=2,FILE=EVALUE(1:J)//'\asfdc',STATUS='OLD')
* }
      END

************************************************************************

      SUBROUTINE OPENFILE(NARG,FILENAME)
*     OPEN A FILE WHOSE ORDER IN THE ARGUMENT IS NARG
*     Cf. Language Reference Manual of Compaq Fortran, E-16.
      USE DFLIB
      CHARACTER*50 FILENAME

*     NARG: ARGUMENT ORDER NUMBER
      CALL GETARG(NARG,FILENAME)

      SELECT CASE (NARG)
         CASE (1)
*           *.ins
*           FILE #5 NEED NOT BE REOPENED IN THE CASE OF Visual Fortran
*           DO NOTHING
            RETURN
         CASE (2)
*           *.int
            NUNIT = 3
         CASE (3)
*           *.bkg
            NUNIT = 8
         CASE (4)
*           *.pat
            NUNIT = 20
         CASE (5)
*           *.hkl
            NUNIT = 21
         CASE (6)
*           *.xyz
            NUNIT = 9
         CASE (7)
*           *.mem
            NUNIT = 30
         CASE (8)
*           *.ffe
            NUNIT = 10
         CASE (9)
*           *.fba
            NUNIT = 32
         CASE (10)
*           *.ffi
            NUNIT = 22
         CASE (11)
*           *.ffo
            NUNIT = 23
*Rev 1.07 2002.08.22 Izumi {
         CASE (12)
*           *.vcs
	      NUNIT = 9998
* }
      END SELECT

      OPEN(UNIT=NUNIT,FILE=FILENAME,STATUS='UNKNOWN')
      END

************************************************************************

*Rev 1.07 2002.08.22 Izumi {
      SUBROUTINE UPDATE_LPP(A,SIGMAA)
*     Update lattice and positional parameters in a VICS file, *.vcs
*     Cf. Language Reference Manual of Compaq Fortran, E-16.
      USE DFLIB
      PARAMETER (NB=7000,NT=999,NS=48,NAP=150,NPH=8)
      REAL A(*),SIGMAA(*)
      INTEGER R,RX,RY,RZ,IG(NAP)
      CHARACTER FILE12*50,LINE*80,LINE2*55
      LOGICAL VICS
      CHARACTER PARNAM*60,PHNAME*25
      COMMON /CC/ APR(NT)
      COMMON /PAR/ PARNAM(NT),PHNAME(NPH)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /RT/ T(3,NS,NPH),HT(NS,NB),NSYM(NPH),NSITE(NPH),
     &  IDSYM(NS,NAP,NPH),R(3,3,NS,NPH),RX(NS,NB),RY(NS,NB),RZ(NS,NB)

      CALL GETARG(12,FILE12)
*     Check whether or not *.vcs (VICS file) exists in the current directory
      INQUIRE(FILE=FILE12,EXIST=VICS)
      IF (.NOT. VICS) RETURN

      NVCS = 9998
      CALL OPENFILE(12,FILE12)

*     Open a scratch file
      NSCR = 9999
      OPEN(UNIT=NSCR,STATUS='SCRATCH')

*     Copy *.vcs into the scratch file
*     Both EXIT and CYCLE can be used within DO ... END DO loops
*     Refer to 7-14 in Compaq Fortran Language Reference Manual
      DO
         READ(NVCS,'(A)',END=1) LINE
         CALL PRLINE(LINE,NSCR)
      END DO
    1 REWIND NVCS
      REWIND NSCR
     
      DO
         READ(NSCR,'(A)',END=9) LINE
         CALL PRLINE(LINE,NVCS)
         IF (LINE(1:5) .EQ. 'CELLP') EXIT
      END DO

      LA = KPHB(1) + NAX
*     Lattice parameters
      WRITE(NVCS,'(F10.5,5F11.5)') (APR(J), J = LA, LA+5)
*     Erros in the lattice parameters
      WRITE(NVCS,'(F10.5,5F11.5)') (SIGMAA(J), J = LA, LA+5)

      DO I = 1, 3
*        Skip three lines
         READ(NSCR,'(A)') LINE
      END DO
	WRITE(NVCS,'(A)') 'STRUC'
 
*     Update g, x, y, and z with their esd's

*     Only the first phase is dealt with
      DO J = 1, NSITE(1)
         READ(NSCR,'(A)') LINE
*        DETERMINE THE PARAMETER NUMBERS OF OCCUPATION FACTORS
         SELECT CASE (J)
            CASE (1)
               IG(J) = KPHB(1) + NG1X
            CASE DEFAULT
               IG(J) = IG(J-1) + NPSITE(J-1,1)
         END SELECT
*        Write site number, element, site name, g, x, y, and z
         WRITE(LINE2,'(A,F8.4,3F11.5)') LINE(1:14),
     &   (A(JJ), JJ = IG(J), IG(J) + 3)
	   WRITE(NVCS,'(A)') LINE2
*        Read the line of esd's
	   READ(NSCR,'(A)') LINE
*        Write esd's of x, y, and g
         WRITE(NVCS,'(F33.5,2F11.5)')
     &   (SIGMAA(JJ), JJ = IG(J) + 1, IG(J) + 3)
      END DO
      
      DO
         READ(NSCR,'(A)',END=2) LINE
         CALL PRLINE(LINE,NVCS)
      END DO

    2 CLOSE(UNIT=NVCS)
      CLOSE(UNIT=NSCR)
	RETURN
    9 CALL JOBEND('*.vcs does not include the line of CELLP')
      END
* }

************************************************************************

      SUBROUTINE XATT
*     A dummy subroutine to do nothing (called in SUBROUTINE OUT20)
      END

