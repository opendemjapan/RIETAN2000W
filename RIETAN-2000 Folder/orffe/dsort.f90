      SUBROUTINE DSORT
*     Sort, print out interatomic distances for atoms in the
*     asymmetric unit, and add instructions 2 to the tail of *.xyz
*     File #5 was changed into #3
*      Input: FILE #3 (Standard input = "orffe output")
*      Input: FILE #9 (*.xyz)
*     Output: FILE #6 (Standard output)
*     Output: FILE #4 (Output file with no carriage control, charcters at 1st columns)
*     Output: FILE #9999 (Scra
*     MDIST: Maximum interatomic distances
*     MCN: Maximum coordination number
*     Modified by F. izumi on 2001.04.23 {
*     PARAMETER (MDIST=10000,MCN=500,NAP=150)
      PARAMETER (MDIST=100000,MCN=5000,NAP=150)
* }
      INTEGER IN(4,MDIST),NUM(MCN),IOUT(MCN),ISITE(3,NAP),
     &  IBOND(2,MCN,NAP)
      REAL XYZ(6,MDIST),DIST(3,MDIST)
      CHARACTER ATOMN(2,MDIST)*9,LINE*137,FNAME*50,LINEXYZ*150
      LOGICAL LEXIST
       
      REWIND 3
      CALL OPEN4(LEXIST)

*     Read the number of atoms in the asymmetric unit
      DO
         READ(3,'(A)') LINE
         CALL WR46(LINE, LEXIST)
         IF (LINE(1:3) .EQ. '  1') THEN
            DO J = 2, 500
               READ(3, '(A)') LINE
               CALL WR46(LINE, LEXIST)
               READ(LINE, '(I3)') NATOM
               IF (NATOM .EQ. 0) GO TO 1
            END DO
         END IF
      END DO
    1 NATOM = J - 1
  
      DO
         READ(3,'(A)') LINE
         CALL WR46(LINE, LEXIST)
         IF (LINE(1:35) .EQ.  '0INTERATOMIC DISTANCE IN ANGSTROMS ')
     &   EXIT
      END DO

*     UNIT = 9: *.xyz
      REWIND 9
*     A scratech file to add instructions 2 to *.xyz
      NSCR = 9999
      OPEN(UNIT = NSCR, STATUS = 'SCRATCH')
      DO 
         READ(9, '(A)', END = 2) LINEXYZ
         CALL PRLINE(LINEXYZ, NSCR)
      END DO
      
    2 J=1
      DO
         READ(3,100,END=6) (IN(I,J), I = 1, 4), ATOMN(1,J),
     &   (XYZ(I,J),I=1,3), ATOMN(2,J), (XYZ(I,J), I = 4, 6)
  100    FORMAT(12X,I2,1X,I5,3X,I2,1X,I5,25X,A9,2X,F7.4,1X,F8.4,1X,F8.4,
     &   5X,A9,2X,F7.4,1X,F8.4,1X,F8.4)
         IF (IN(1,J) .EQ. 0) EXIT
         READ(3,'(53X,F9.4,5X,F7.4,2X,F7.4)') (DIST(I,J), I = 1, 3)
*        SKIP IF THE INTERATOMIC DISTANCE IS ZERO
         IF (DIST(1,J) .EQ. 0.0) CYCLE

*        Swap a pair of atoms if the 2nd atom is within the asymmetric unit
         IF (IN(4,J) .EQ. 0) THEN
            IN(1,J+1) = IN(3,J)
            IN(2,J+1) = IN(4,J)
            IN(3,J+1) = IN(1,J)
            IN(4,J+1) = IN(2,J)
            ATOMN(1,J+1) = ATOMN(2,J)
            ATOMN(2,J+1) = ATOMN(1,J)
            XYZ(1,J+1) = XYZ(4,J)
            XYZ(2,J+1) = XYZ(5,J)
            XYZ(3,J+1) = XYZ(6,J)
            XYZ(4,J+1) = XYZ(1,J)
            XYZ(5,J+1) = XYZ(2,J)
            XYZ(6,J+1) = XYZ(3,J)
            DIST(1,J+1) = DIST(1,J)
            DIST(2,J+1) = DIST(2,J)
            DIST(3,J+1) = DIST(3,J)
            J=J+1
         END IF
         J = J + 1
         IF (J. GT. MDIST - 1) THEN
            WRITE(6, '(A)') ' INCREASE MDIST'
            STOP
         END IF
      END DO
    6 NDIST = J - 1

*      ISITE(K,I)  K = 1: A for site I
*                  K = 2: 1000C + S for site I
*                  K = 3: Total number of bonds for site I
*                      I: Site number
*     BOND(K,J,I)  K = 1: A for the Jth atoms bonded to the atoms at site I
*                  K = 2: 1000C + S for the Jth atoms bonded to the atom at site I
*                      J: Number for an atom bonded to the atom at site I
*                      I: Site number
*     IORDER: Sseirs numbers for interatomic distances and bond angles
      IORDER = 1
      DO I = 1, NATOM
         INUM = 1
         DO JJ = 1, NDIST
            IF (IN(1,JJ) .EQ. I) THEN
               NUM(INUM) = JJ
*              IOUT = 0: This distance has not been printed out
*              IOUT = 1: This distance has already been printed out
               IOUT(INUM) = 0
               INUM = INUM + 1
               IF (INUM .GT. MCN) THEN
                  WRITE(6,'(A)') ' INCREASE MCN'
                  STOP
               END IF
            END IF
         END DO
         INUM = INUM - 1
         IF (INUM .EQ. 0) CYCLE

         DO J = 1, MCN
*           Extremely large distance (dummy)
            DMIN = 10000.0
            DO JJ = 1, INUM 
               IF (IOUT(JJ) .EQ. 0 .AND. DIST(1,NUM(JJ)) .LT. DMIN) THEN
                  DMIN = DIST(1,NUM(JJ))
                  JOUT = JJ
               END IF
            END DO
         
            IF (ABS(DMIN - 10000.0) .LT. 0.1) THEN
*              End of distances for this site
               WRITE(LINE, '(1H0, 129A1)') ('-', JJ = 1, 129)
               CALL WR46(LINE, LEXIST)
               ISITE(3,I) = J - 1
               EXIT
            END IF
            
            WRITE(LINE,150) IORDER, (IN(JJ, NUM(JOUT)), JJ = 1, 4),
     &      ATOMN(1,NUM(JOUT)), (XYZ(JJ, NUM(JOUT)), JJ = 1, 3),
     &      ATOMN(2,NUM(JOUT)), (XYZ(JJ, NUM(JOUT)), JJ = 4, 6)
            CALL WR46(LINE, LEXIST)
  150       FORMAT('0',I5,5X,'(',I2,1H,I5,3H) (I2,1H,I5,1H),24X, A9,
     &      ' (',F7.4,',',F8.4,',',F8.4,');   ',A9,' (',F7.4,',',
     &      F8.4,',',F8.4,')')
     
            WRITE(LINE,'(1H 52X,F9.4,5H +OR-F7.4,2H (F7.4,1H))') 
     &      (DIST(JJ,NUM(JOUT)), JJ = 1, 3) 
            CALL WR46(LINE, LEXIST)
            IOUT(JOUT) = 1
               
            IF (J .EQ. 1) THEN
               ISITE(1,I) = IN(1,NUM(JOUT))
               ISITE(2,I) = IN(2,NUM(JOUT))
            END IF
            IBOND(1,J,I) = IN(3,NUM(JOUT))
            IBOND(2,J,I) = IN(4,NUM(JOUT))
            
            IORDER = IORDER + 1
         END DO
      END DO

      IF (.NOT. LEXIST) WRITE(4, '()')
      WRITE(6, '()')
      LANGLE = 0
      DO
         READ(3,'(A)', END = 8) LINE
         IF (LINE(12:12) .EQ. '(') THEN
            WRITE(LINE, '(A, I5, 5X, A)') '0', IORDER, LINE(12:)
            CALL WR46(LINE, LEXIST)
            IORDER = IORDER + 1
            LANGLE = 1
         ELSE
            CALL WR46(LINE, LEXIST)         
         END IF
      END DO
    8 CLOSE(UNIT = 3)
      IF (.NOT. LEXIST) CLOSE(UNIT = 4)
      CLOSE(UNIT = 6)
*     IF Instruction 2 has already been included in *.xyz, automatic 
*     generation of instructions 2 is skipped
      IF (LANGLE .EQ. 1) RETURN
      
*     A series of instructions 2 to calculate bond angle J1-I-J3
      DO I = 1, NATOM
         DO J1 = 1, ISITE(3,I)
            DO J3 = J1 + 1, ISITE(3,I) 
               WRITE(NSCR,'(7I5)') 2, IBOND(1,J1,I), IBOND(2,J1,I),
     &         ISITE(1,I), ISITE(2,I), IBOND(1,J3,I), IBOND(2,J3,I)
            END DO
         END DO
      END DO
      
*     Copy the whole content of the scratch file to *.xyz
      REWIND 9
      REWIND NSCR
      DO
         READ(NSCR, '(A)', END = 9) LINEXYZ
         CALL PRLINE(LINEXYZ, 9)
      END DO
    9 CLOSE(UNIT = 9)
      CLOSE(UNIT = NSCR)
      END
      
************************************************************************

      SUBROUTINE WR46(LINE, LEXIST)
      CHARACTER*(*) LINE
      LOGICAL LEXIST
      
      IF (.NOT. LEXIST) WRITE(4, '(2A)') ' ',LINE(2:)
      WRITE(6, '(A)') LINE
      END
      
************************************************************************

      SUBROUTINE PRLINE(LINE, IOUT)
*     DETERMINE THE LINE LENGTH AND PRINT OUT LINE
      CHARACTER*150 LINE
      
      DO J = 150, 1, -1
         IF (LINE(J:J) .NE. ' ') EXIT
      END DO
      IF (J .EQ. 0) THEN
         WRITE(IOUT, '()')
      ELSE
         WRITE(IOUT, '(A)') LINE(1:J)
      END IF
      END
