      SUBROUTINE SLESUM(X,N,A,SUM)
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1985-08-02 SLESUM Lawson  Initial code.
C
C     THIS SUBROUTINE EVALUATES THE SUM OF
C          A(J) * P(J) FOR J = 0,...,N,
C     WHERE P(J)'S  ARE LEGENDRE POLYNOMIALS OF DEGREE J.
C
C     THE RECURSION FORMULA IS :
C     B(K) = X+B(K+1)*(2*K+1)/(K+1)-B(K+2)*(K+1)/(K+2)+A(K)
C
C     C.L.LAWSON & S.CHAN, JPL, 1983 JUNE 9
C
C     -------------------------------------------------------------
C     SUBROUTINE ARGUMENTS
C     --------------------
C     X     ARGUMENT OF LEG POLYS, X SHOULD BE NON-NEGATIVE.
C     N     SUM IS TO INCLUDE LEGENDRE POLYS OF DEGREE ZERO
C           THRU N.
C     A()   A,...,A(N) CONTAIN COEFFS TO BE USED IN
C           FORMING THE SUM.
C     SUM   SUM OF COMBINATION
C
C     -------------------------------------------------------------
C
      REAL A(0:N), X
      DATA ZERO,ONE,TWO / 0.E0, 1.E0, 2.E0 /
C
      IF (N .GE. 0) GO TO 5
        SUM = ZERO
      RETURN
C
   5  C1 = REAL(N+1)
      C3 = C1 + C1 - ONE
      B1 = ZERO
      B  = A(N)
C
      DO 10 K = N-1,0,-1
C
C     C1 = K + 1, C3 = 2K + 1
C
        C2 = C1
        C1 = C1 - ONE
        C3 = C3 - TWO
        B2 = B1
        B1 = B
        B  = X * B1 * C3 / C1  -  B2 * C1 / C2  +  A(K)
  10  CONTINUE
      SUM = B
      RETURN
      END

************************************************************************

      subroutine AMACH(MODE, I, I1, R1, D1)
c>> 1992-04-07 AMACH  Oken    Removed ^Z at EOF (error found by VAX compile)
c>> 1992-02-20 AMACH  Snyder  Added Cray-YMP stuff, q.v.
c>> 1990-06-11 AMACH  Snyder  Added Apollo DN-10000 stuff, q.v.
c>> 1990-12-14 AMACH  Lawson  Changed to eliminate ENTRY statements.
c>> 1990-08-21 AMACH  Krogh   No test was getting done for bad machine.
c>> 1990-02-28 AMACH  Krogh   Correct missing DOUBLE PRECISION AMSUB1
c>> 1989-08-14 AMACH  Krogh   Parameterized everything -- Massive change
c>> 1989-03-30 AMACH  Snyder  Correct missing "/" line 921
c>> 1989-01-30 AMACH  Snyder  Incorporate more constants from NETLIB.
C>> 1988-05-19 AMACH  Lawson  Initial code.
c File AMACH.FOR contains user-callable functions I1MACH, D1MACH, and
c R1MACH, plus second-level subroutines AMACH, AMTEST, and AMSUB1.
c Appropriate lines must be switched between comment and non-comment
c status when this code is moved to a different computer system.
c     These changes can be done with any text editor, however the "c++"
c lines permit automation of the change using the MARVEL processor.
c Note that when the MARVEL processor activates a line it shifts
c Columns 2-72 to 1-71 and puts a blank in Column 72.  When it inactiv-
c ates a line it shifts Columns 1-71 to 2-72 and puts a C in Column 1.
c     The possible choices using MARVEL are:
c              c++ SET SYS = IEEE
c              c++ SET SYS = AMDAHL 
c              c++ SET SYS = APOLLO-10000
c              c++ SET SYS = BUR1700 
c              c++ SET SYS = BUR5700 
c              c++ SET SYS = BUR67-7700 
c              c++ SET SYS = CDC60-7000 
c              c++ SET SYS = CONVEXC-1 
c              c++ SET SYS = CRAY1 
c              c++ SET SYS = CRAY1-SD (Sngl prec.arith. used for dble.)
c              c++ SET SYS = CRAY1-64 (64 bit integers)
c              c++ SET SYS = CRAY1-SD-64 (64 bit int, SP used for DP)
c              c++ SET SYS = CRAY-YMP
c              c++ SET SYS = DG-S2000 
c              c++ SET SYS = HARRIS220 
c              c++ SET SYS = HON600-6000 
c              c++ SET SYS = HON-DPS-8-70 
c              c++ SET SYS = IBM360-370 
c              c++ SET SYS = INTERDATA-8-32 
c              c++ SET SYS = PDP10-KA 
c              c++ SET SYS = PDP10-KB 
c              c++ SET SYS = PDP11 
c              c++ SET SYS = PRIME50 
c              c++ SET SYS = SEQ-BAL-8000 
c              c++ SET SYS = UNIVAC 
c              c++ SET SYS = VAX 
c     The current choice is:
c++ SET SYS = IEEE
c
C  I/O UNIT NUMBERS:
C
C    IM1 = I1MACH( 1) = THE STANDARD INPUT UNIT.
C    IM2 = I1MACH( 2) = THE STANDARD OUTPUT UNIT.
C    IM3 = I1MACH( 3) = THE STANDARD PUNCH UNIT.
C    IM4 = I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
C
C  WORDS:
C
C    IM5 = I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
C    IM6 = I1MACH( 6) = THE NUMBER OF CHARACTERS/INTEGER STORAGE UNIT.
C
C  INTEGERS:
C
C    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
C
C               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
C
C    IM7 = I1MACH( 7) = A, THE BASE.
C    IM8 = I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
C    IM9 = I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
C
C  FLOATING-POINT NUMBERS:
C
C    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
C    BASE-B FORM
C
C               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
C               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
C
C    IM10 = I1MACH(10) = B, THE BASE.
C
C  SINGLE-PRECISION:
C
C    IM11 = I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
C    IM12 = I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
C    IM13 = I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
C
C  DOUBLE-PRECISION:
C
C    IM14 = I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
C    IM15 = I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
C    IM16 = I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
C
C  CONVERSION FROM FUNCTIONAL TO STRUCTURAL FLOATING POINT CONSTANTS
C
C    IM17 = CONSTANT SUCH THAT IM14 + IM17 = ACTUAL NUMBER OF BASE-B
C           DIGITS IN DOUBLE PRECISION, USED FOR CHECKING THAT CORRECT
C           VERSION OF THIS PROGRAM IS INSTALLED.  (SEE DEFINITION OF
C           DM6, AND THE USE OF DM6 IN CALLING AMTEST.)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF PARAMETER STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF
C  IM1 - IM4 SHOULD BE CHECKED FOR CONSISTENCY
C  WITH THE LOCAL OPERATING SYSTEM.
c     -----------------------------------------------------------------
c     Original design and code due to P. A. Fox, A. D. Hall, and
c     N. L. Schryer, Bell Laboratories.  See ACM TOMS, 4,(1978),177-188.
c     Adapted to Univac 1100 by Kris Stewart, JPL, 7/30/81.
c     Adapted for the JPL MATH77 library by C. L. Lawson and F. T. Krogh
c     Sept, 1987.
c     1989-08-14 AMACH  Krogh   Parameterized everything. Major changes.
C     1990 Dec. CLL reorganized code to avoid using ENTRY statements
c     for functions of different types.  Also added save statements.
c     -----------------------------------------------------------------
c     On the first call to this function, tests are done to verify that
c     IM10 and IM14 are not grossly wrong for the host environment.
c     This gives some protection against using the wrong version of this
c     subprogram.
c     -----------------------------------------------------------------
      integer MODE, I, I1
      real R1
      double precision D1, TEST
c
      integer IMACH(16)
      integer IM1, IM2, IM3, IM4, IM5, IM6, IM7, IM8, IM9, IM10, IM11,
     1   IM12, IM13, IM14, IM15, IM16, IM17
      real             RMACH(5), RM1, RM2, RM3, RM4, RM5,
     1                 RMA, RMB, RBASE
      double precision DMACH(5), DM1, DM2, DM3, DM4, DM5, DM6,
     1                 DMA, DMB, DBASE
      save TEST, IMACH, RMACH, DMACH
C     -----------------------------------------------------------------
C     Machine constants for IEEE standard binary floating-point
c     processors.  This includes PC's and work-stations using the
c     Intel 8087, 80287, 80387, ... processors or the
c     Motorola 68881, 68882, ... processors.
c     Note:  We are setting the "most negative exponent" (IMACH(12) and
c     IMACH(15)) to be the exponent of the smallest normalized number.
c     An IEEE processor actually handles smaller numbers before
c     underflowing, however these "unnormalized" numbers have
c     diminished precision.
c
c++ Code for SYS = IEEE is ACTIVE
       PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
       PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
       PARAMETER (IM10 =2, IM11 =24, IM12 =-125, IM13 =128)
       PARAMETER (IM14 =53, IM15 =-1021, IM16 =1024, IM17=0)
C     -----------------------------------------------------------------
c++ Code for SYS = AMDAHL is INACTIVE
CC     MACHINE CONSTANTS FOR AMDAHL MACHINES.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
C      PARAMETER (IM10 =16, IM11 =6, IM12 =-64, IM13 =63)
C      PARAMETER (IM14 =14, IM15 =-64, IM16 =63, IM17=0)
C      -----------------------------------------------------------------
c++ Code for SYS = APOLLO-10000 is INACTIVE
cc     MACHINE CONSTANTS FOR APOLLO DN-10000 MACHINES.
cc     The only difference from IEEE is IM13.  This difference has
cc     nothing to do with the arithmetic or representation used by the
cc     machine.  It is caused by a bug in the compiler:  The right-hand
cc     side of RM2 (below) is apparently evaluated in double precision.
cc     When the compiler is ready to store the resulting value into its
cc     internal data structures, it compares it to an incorrect value
cc     of the overflow limit.  It appears the incorrect value has the
cc     correct exponent, but the fraction is 1.5 instead of 2-2**(-p),
cc     where p is the precision in bits.  You can get the correct result
cc     by changing IM13 to 128, changing RM2 from a parameter to a
cc     variable, and changing the parameter statement that assigns a
cc     value to RM2 into an ordinary assignment statement.
CC
c      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
c      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
c      PARAMETER (IM10 =2, IM11 =24, IM12 =-125, IM13 =127)
c      PARAMETER (IM14 =53, IM15 =-1021, IM16 =1024, IM17 =0)
CC     -----------------------------------------------------------------
c++ Code for SYS = BUR1700 is INACTIVE
CC     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
CC
C      PARAMETER (IM1 =7, IM2 =2, IM3 =2, IM4 =2)
C      PARAMETER (IM5 =36, IM6 =4, IM7 =2, IM8 =33)
C      PARAMETER (IM10 =2, IM11 =24, IM12 =-256, IM13 =255)
C      PARAMETER (IM14 =60, IM15 =-256, IM16 =255, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = BUR5700 is INACTIVE
CC     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =48, IM6 =6, IM7 =2, IM8 =39)
C      PARAMETER (IM10 =8, IM11 =13, IM12 =-50, IM13 =76)
C      PARAMETER (IM14 =26, IM15 =-50, IM16 =76, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = BUR67-7700 is INACTIVE
CC     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =48, IM6 =6, IM7 =2, IM8 =39)
C      PARAMETER (IM10 =8, IM11 =13, IM12 =-50, IM13 =76)
C      PARAMETER (IM14 =26, IM15 =-32754, IM16 =32780, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = CDC60-7000 is INACTIVE
CC     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =60, IM6 =10, IM7 =2, IM8 =48)
C      PARAMETER (IM10 =2, IM11 =47, IM12 =-929, IM13 =1070)
C      PARAMETER (IM14 =94, IM15 =-929, IM16 =1069, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = CONVEXC-1 is INACTIVE
CC     MACHINE CONSTANTS FOR CONVEX C-1.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
C      PARAMETER (IM10 =2, IM11 =24, IM12 =-128, IM13 =127)
C      PARAMETER (IM14 =53, IM15 =-1024, IM16 =1023, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = CRAY1 is INACTIVE
CC     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
C      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =46)
C      PARAMETER (IM10 =2, IM11 =47, IM12 =-8189, IM13 =8190)
C      PARAMETER (IM14 =94, IM15 =-8099, IM16 =8190, IM17=2)
CC     -----------------------------------------------------------------
c++ Code for SYS = CRAY-YMP is INACTIVE
CC     MACHINE CONSTANTS FOR THE CRAY YMP
CC     Cray claims the overflow exponent (IM13 and IM16) is 8189, and
CC     the underflow exponent (IM12 and IM15) is -8189, but these values
CC     don't seem to work in cf77:  the underflow limit underflows, and
CC     the overflow limit overflows when using Cray's values.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
C      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =46)
C      PARAMETER (IM10 =2, IM11 =47, IM12 =-8188, IM13 =8189)
C      PARAMETER (IM14 =94, IM15 =-8188, IM16 =8189, IM17=2)
CC     -----------------------------------------------------------------
c++ Code for SYS = CRAY1-SD is INACTIVE
CC     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3, WHEN DOUBLE
CC     PRECISION IS TO USE SINGLE PRECISION ARITHMETIC.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
C      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =46)
C      PARAMETER (IM10 =2, IM11 =47, IM12 =-8189, IM13 =8190)
C      PARAMETER (IM14 =47, IM15 =-8189, IM16 =8190, IM17=1)
CC     -----------------------------------------------------------------
c++ Code for SYS = CRAY1-64 is INACTIVE
CC     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
C      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =63)
C      PARAMETER (IM10 =2, IM11 =47, IM12 =-8189, IM13 =8190)
C      PARAMETER (IM14 =94, IM15 =-8099, IM16 =8190, IM17=2)
CC     -----------------------------------------------------------------
c++ Code for SYS = CRAY1-SD-64 is INACTIVE
CC     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3, WHEN DOUBLE
CC     PRECISION IS TO USE SINGLE PRECISION ARITHMETIC.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
C      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =63)
C      PARAMETER (IM10 =2, IM11 =47, IM12 =-8189, IM13 =8190)
C      PARAMETER (IM14 =47, IM15 =-8189, IM16 =8190, IM17=1)
CC     -----------------------------------------------------------------
c++ Code for SYS = DG-S2000 is INACTIVE
CC     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
CC
C      PARAMETER (IM1 =11, IM2 =12, IM3 =8, IM4 =10)
C      PARAMETER (IM5 =16, IM6 =2, IM7 =2, IM8 =15)
C      PARAMETER (IM10 =16, IM11 =6, IM12 =-64, IM13 =63)
C      PARAMETER (IM14 =14, IM15 =-64, IM16 =63, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = HARRIS220 is INACTIVE
CC     MACHINE CONSTANTS FOR THE HARRIS 220, SLASH 6, SLASH 7.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =0, IM4 =6)
C      PARAMETER (IM5 =24, IM6 =3, IM7 =2, IM8 =23)
C      PARAMETER (IM10 =2, IM11 =23, IM12 =-127, IM13 =127)
C      PARAMETER (IM14 =38, IM15 =-127, IM16 =127, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = HON600-6000 is INACTIVE
CC     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =43, IM4 =6)
C      PARAMETER (IM5 =36, IM6 =6, IM7 =2, IM8 =35)
C      PARAMETER (IM10 =2, IM11 =27, IM12 =-127, IM13 =127)
C      PARAMETER (IM14 =63, IM15 =-127, IM16 =127, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = HON-DPS-8-70 is INACTIVE
CC     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =43, IM4 =6)
C      PARAMETER (IM5 =36, IM6 =4, IM7 =2, IM8 =35)
C      PARAMETER (IM10 =2, IM11 =27, IM12 =-127, IM13 =127)
C      PARAMETER (IM14 =63, IM15 =-127, IM16 =127, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = IBM360-370 is INACTIVE
CC     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
CC     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
C      PARAMETER (IM10 =16, IM11 =6, IM12 =-64, IM13 =63)
C      PARAMETER (IM14 =14, IM15 =-64, IM16 =63, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = INTERDATA-8-32 is INACTIVE
CC     MACHINE CONSTANTS FOR THE INTERDATA 8/32
CC     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
CC 
C      PARAMETER (IM1 =5, IM2 =6, IM3 =6, IM4 =6)
C      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
C      PARAMETER (IM10 =16, IM11 =6, IM12 =-64, IM13 =62)
C      PARAMETER (IM14 =14, IM15 =-64, IM16 =62, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = PDP10-KA is INACTIVE
CC     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =36, IM6 =5, IM7 =2, IM8 =35)
C      PARAMETER (IM10 =2, IM11 =27, IM12 =-128, IM13 =127)
C      PARAMETER (IM14 =54, IM15 =-101, IM16 =127, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = PDP10-KB is INACTIVE
CC     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
CC
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =36, IM6 =5, IM7 =2, IM8 =35)
C      PARAMETER (IM10 =2, IM11 =27, IM12 =-128, IM13 =127)
C      PARAMETER (IM14 =62, IM15 =-128, IM16 =127, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = PDP11 is INACTIVE
CC     MACHINE CONSTANTS FOR PDP-11 FORTRAN'S SUPPORTING
CC     16-BIT INTEGER ARITHMETIC.
C
C      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C      PARAMETER (IM5 =16, IM6 =2, IM7 =2, IM8 =15)
C      PARAMETER (IM10 =2, IM11 =24, IM12 =-127, IM13 =127)
C      PARAMETER (IM14 =56, IM15 =-127, IM16 =127, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = PRIME50 is INACTIVE
CC     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
CC     WITH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
CC     SUPPLIED BY IGOR BRAY.
C
C      PARAMETER (IM1 =1, IM2 =1, IM3 =2, IM4 =1)
C      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
C      PARAMETER (IM10 =2, IM11 =23, IM12 =-127, IM13 =127)
C      PARAMETER (IM14 =47, IM15 =-32895, IM16 =32637, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = SEQ-BAL-8000 is INACTIVE
CC     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
CC
C      PARAMETER (IM1 =0, IM2 =0, IM3 =7, IM4 =0)
C      PARAMETER (IM5 =32, IM6 =1, IM7 =2, IM8 =31)
C      PARAMETER (IM10 =2, IM11 =24, IM12 =-125, IM13 =128)
C      PARAMETER (IM14 =53, IM15 =-1021, IM16 =1024, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = UNIVAC is INACTIVE
CC     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
CC
CC     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 1
CC     WHICH IS APPROPRIATE FOR THE UNIVAC-FTN SYSTEM.
CC     IF YOU HAVE THE UNIVAC-FOR SYSTEM, SET IT TO 7.
CC     IM6 = 4 for FTN (4 chars per word), 6 for FOR (6 chars per word).
Cc
C      PARAMETER (IM1 =5, IM2 =6, IM3 =1, IM4 =6)
C      PARAMETER (IM5 =36, IM6 =4, IM7 =2, IM8 =35)
C      PARAMETER (IM10 =2, IM11 =27, IM12 =-128, IM13 =127)
C      PARAMETER (IM14 =60, IM15 =-1024, IM16 =1023, IM17=0)
CC     -----------------------------------------------------------------
c++ Code for SYS = VAX is INACTIVE
C     MACHINE CONSTANTS for the VAX/VMS
C     and for PDP-11 FORTRAN SUPPORTING 32-BIT INTEGER ARITHMETIC.
C
C     PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
C     PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
C     PARAMETER (IM10 =2, IM11 =24, IM12 =-127, IM13 =127)
C     PARAMETER (IM14 =56, IM15 =-127, IM16 =127, IM17=0)
c++ end
C     -----------------------------------------------------------------
C
      PARAMETER (IM9 = 2 * (2**(IM8-1) - 1) + 1)
C
C Real parameters
C
C  RM1 = R1MACH(1) = B**(EMIN-1), The smallest positive number, i.e.,
c                    the underflow limit.
C  RM2 = R1MACH(2) = B**EMAX*(1 - B**(-T)), The largest number, i.e.,
c                    the overflow limit.
C  RM3 = R1MACH(3) = B**(-T), The smallest relative spacing, i.e., the
c                    difference between 1.0 and the next smaller number.
C  RM4 = R1MACH(4) = B**(1-T), The largest relative spacing, i.e., the
c                     difference between 1.0 and the next larger number.
C  RM5 = R1MACH(5) = LOG10(B).  When B = 2 this value is
c              Log10(2) = 0.30102_99956_63981_19521_37388_94724
C
C Parameter RMA and RMB are selected so that for values of the base =
C 2, 8, 16, 10, RMA has the values 1, 3, 4, 0, and RMB has the values 0,
C 0, 0, 1.  These values are used in computing RM5.
C $$$$ Note that if other bases are to be supported, the calculation of
C $$$$ RMA and RMB will have to be generalized.
C
      PARAMETER (RMA = ((IM10 - 10) * (-3 + ((IM10 - 2) * (-77 +
     1    12 * (IM10 - 8))) / 14)) / 24)
      PARAMETER (RMB = ((IM10 - 2) * (IM10 - 8) * (16 - IM10)) / 96)
      PARAMETER (RBASE = IM10)
C
C     Weird subterfuges below are NECESSARY to compute DM1 and DM2 on
C     some systems.  DON'T SIMPLIFY THEM.  We compute RM1 and RM2 using
C     these subterfuges so it will be clear we're computing the REAL and
C     DOUBLE PRECISION characteristics in the same way.
      PARAMETER (RM1 = (RBASE**(IM12/2)) * (RBASE**(IM12-IM12/2-1)))
      PARAMETER (RM2 = RBASE**(IM13-IM11) * ((RBASE**IM11 - RBASE)
     1               + (RBASE - 1.0E0)))
      PARAMETER (RM3 = RBASE**(-IM11))
      PARAMETER (RM4 = RBASE**(1-IM11))
c     PARAMETER (RM5 = RMA*0.30102 99956 63981 19521 37388 94724E0+RMB)
      PARAMETER (RM5 = RMA*0.301029995663981195213738894724E0+RMB)
C
C Double precision paramters -- (Defined like the real ones.)
C
      PARAMETER (DMA = ((IM10 - 10) * (-3 + ((IM10 - 2) * (-77 +
     1    12 * (IM10 - 8))) / 14)) / 24)
      PARAMETER (DMB = ((IM10 - 2) * (IM10 - 8) * (16 - IM10)) / 96)
      PARAMETER (DBASE = IM10)
C
C     Weird subterfuges below are NECESSARY to compute DM1 and DM2 on
C     some systems.  DON'T SIMPLIFY THEM.
      PARAMETER (DM1 = (DBASE**(IM15/2)) * (DBASE**(IM15-IM15/2-1)))
      PARAMETER (DM2 = DBASE**(IM16-IM14) * ((DBASE**IM14 - DBASE)
     1               + (DBASE - 1.0D0)))
      PARAMETER (DM3 = DBASE**(-IM14))
      PARAMETER (DM4 = DBASE**(1-IM14))
c     PARAMETER (DM5 = DMA*0.30102 99956 63981 19521 37388 94724D0+DMB)
      PARAMETER (DM5 = DMA*0.301029995663981195213738894724D0+DMB)
C DM6 and TEST are used in checking that the correct constants have been
C selected.
      PARAMETER (DM6 = DBASE**(-IM14-IM17))
      data TEST / 0.D0 /
C
c     DATA IMACH / IM1, IM2, IM3, IM4, IM5, IM6, IM7, IM8, IM9, IM10,
c    1   IM11, IM12, IM13, IM14, IM15, IM16 /
c     DATA RMACH / RM1, RM2, RM3, RM4, RM5 /
c     DATA DMACH / DM1, DM2, DM3, DM4, DM5 /
C     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if (TEST .eq. 0.0D0) then
         IMACH(1) = IM1
         IMACH(2) = IM2
         IMACH(3) = IM3
         IMACH(4) = IM4
         IMACH(5) = IM5
         IMACH(6) = IM6
         IMACH(7) = IM7
         IMACH(8) = IM8
         IMACH(9) = IM9
         IMACH(10) = IM10
         IMACH(11) = IM11
         IMACH(12) = IM12
         IMACH(13) = IM13
         IMACH(14) = IM14
         IMACH(15) = IM15
         IMACH(16) = IM16
         RMACH(1) = RM1
         RMACH(2) = RM2
         RMACH(3) = RM3
         RMACH(4) = RM4
         RMACH(5) = RM5
         DMACH(1) = DM1
         DMACH(2) = DM2
         DMACH(3) = DM3
         DMACH(4) = DM4
         DMACH(5) = DM5
         CALL AMTEST (TEST, DM6)
      ENDIF

      if (MODE .eq. 0) then
         I1=IMACH(I)
      else if (MODE .eq. 1) then
         R1=RMACH(I)
c                                  Here we assume MODE = 2.
      else
         D1=DMACH(I)
      endif
      return
      end
c     ==================================================================
      integer function I1MACH(I)
      integer I, I1
      real R1
      double precision D1
      IF (I .LT. 1  .OR.  I .GT. 16) THEN
         PRINT*,'I1MACH.. Bad argument: I =',I
         STOP 'I1MACH error'
      END IF
      call AMACH (0, I, I1, R1, D1)
      I1MACH = I1
      return
      end
c     ==================================================================
c
      real function R1MACH(I)
      integer I, I1
      real R1
      double precision D1
      IF (I .lt. 1  .or.  I .gt. 5) THEN
         print*,'R1MACH.. Bad argument: I = ',I
         stop 'R1MACH error'
      END IF
      call AMACH (1, I, I1, R1, D1)
      R1MACH = R1
      RETURN
      end
c     ==================================================================
c
      double precision function D1MACH(I)
      integer I, I1
      real R1
      double precision D1
      IF (I .lt. 1  .or.  I .gt. 5) THEN
         print*,'D1MACH.. Bad argument: I = ',I
         stop 'D1MACH error'
      END IF
      call AMACH (2, I, I1, R1, D1)
      D1MACH = D1
      RETURN
      END
c     ==================================================================
c
      SUBROUTINE AMTEST (TEST, D6)
c Verifies that D6 is an appropriate values for DM6.
      DOUBLE PRECISION AMSUB1, D6, TEST
      TEST = AMSUB1(1.D0 + D6)
C
C The comparison with 1.875E0*D6 in the line below is to guard
C against the possibility that TEST is > 0 as a result of rounding
C up in the addition of D6 to 1.
C
      IF ((TEST .eq. 0.D0) .or. (TEST .gt. 1.875D0*D6)) THEN
         TEST = (D6 + D6) + 1.D0
         IF (AMSUB1(TEST) .ne. 0.D0) RETURN
      END IF
* Commented out by F. Izumi on 2000.12.30 because this error occurs
* whenever function erf is calculated
C     print*,'AMACH has bad parameters for current environment.'
C     stop
      WRITE(6,'(//11X,A)') 
     &'AMACH has bad parameters for current environment (probably, no pr
     &oblem).'
      END
c     ==================================================================
c
      DOUBLE PRECISION FUNCTION AMSUB1 (TEST1)
      DOUBLE PRECISION TEST1
C     Returns the value of TEST1 - 1.
      AMSUB1 = TEST1 - 1.0D0
      RETURN
      END

c++   Set VERSION = F
c++ CODE for VERSION = C is inactive
c%% static FILE *c_handle[2], *scratch_file;
c%% static char *c_fname[2]={"MESSF-xx", "MESSF-xx"};
c%% char *ctmp;
c++ END
      subroutine MESS(MACT, TEXT, IDAT)
c     .  Copyright (C) 1991, California Institute of Technology.
c     .  All rights reserved.  U. S. Government sponsorship under
c     .  NASA contract NAS7-918 is acknowledged.
c>> 1993-05-19 MESS  Krogh  Changed TEXT to array of character strings.
c>> 1993-04-14 MESS  Krogh  Fixes for conversion to C. (C%% comments.)
c>> 1993-03-10 MESS  Krogh  Broke into smaller pieces.
c>> 1992-12-02 MESS  Krogh  Added save statment in the block data subpr.
c>> 1992-07-13 MESS  Krogh  Add checks in heading set up.
c>> 1992-07-12 MESS  Krogh  Fixed so $$ prints a single $ in TEXT.
c>> 1992-07-12 MESS  Krogh  Set out of bound inputs to limit values.
c>> 1992-07-12 MESS  Krogh  Fixed so output works to alternate files.
c>> 1992-07-12 MESS  Krogh  Added integer declarations for parameters.
c>> 1992-06-24 MESS  Krogh  More blanks allowed on break of long lines.
c>> 1992-06-10 MESS  Krogh  Minor fix to vector output.
c>> 1992-05-27 MESS  Krogh  Fixed bug on line width setting.
c>> 1992-05-14 MESS  Krogh  Put common blocks in save statement.
c>> 1992-05-11 MESS  Krogh  Added label to assigned go to & a comment.
c>> 1992-04-08 MESS  Krogh  Unused labels 60, 220 and 320 removed.
c>> 1992-03-20 MESS  Krogh  Changed status on open to SCRATCH.
c>> 1992-03-10 MESS  Krogh  1 line below label 690 changed max to min.
c>> 1992-02-05 MESS  Krogh  Fixed bugs in printing matrix labels.
c>> 1992-01-29 MESS  Krogh  Added UMESS and multiple print option.
c>> 1991-12-09 MESS  Krogh  Fine tuning of vector output.
c>> 1991-10-10 MESS  Krogh  Insure no stop if stop level = 9.
c>> 1991-06-26 MESS  Krogh  Initial Code.
c Processes Messages -- Actions are controlled by MACT().
c This routine is intended for use primarily by other library routines.
c Users of library routines may want to use values of MACT from MERET-
c MESUNI, and may have an interest in using it to print messages
c from their own software.
c This routine has companion routines that are called with the same
c three arguments, plus one additional argument.  This argument is
c referred to here as FDAT since actions specified here can result
c in returns to print data from FDAT.  The name FDAT is used because
c this other routine will ordinarily print floating point data, but
c it could also print other kinds of data, e.g. logical.  At present
c only SMESS and DMESS are defined which are for single and double
c precision floating point data.
c MACT is a vector specifying sequentially the actions desired.
c Some of these actions require more than one location, in which
c case the next action follows the last datum required by the
c previous action.  Internal variables together with default
c values in parentheses which are used to keep track of locations
c are as follows:
c  NTEXT  (1)   The next text output starts at TEXT(NTEXT).
c  NIDAT  (1)   The next output from IDAT starts at IDAT(NIDAT).
c  NFDAT  (1)   The next output from FDAT starts at FDAT(NFDAT).
c  NMDAT  (1)   The next output from MDAT starts at MDAT(NMDAT), where
c               MDAT is defined by actions MEMDA1-MEMDA5 below, and
c               NMDAT is set to one at end of every text output.
c An action which uses data pointed to by one of the above will cause
c the pointer to be incrmented to one past the last location used.  An
c exception is NMDAT which when it reaches 5 is not incremented and the
c value pointed to is incremented instead.
c Actions are encoded by values starting in MACT(1) as follows.
c (Arguments required are given in parentheses at the start of
c description.  These arguments follow the action index.  The next
c action follows the last argument for the preceding action.  Action
c indices have been selected so that it is easy to add new functionality
c without affecting codes using an earlier version.  Where bounds are
c indicated for an argument, if the argument is outside the bounds it is
c treated as if it had the value for the bound violated.)
c MESUNI=10  (0 .le. K10 .le. 99) Set the unit to use for a scratch
c            file.  The default unit for a scratch file is 30.  If a
c            scratch file is needed, (only needed here if a table
c            exceeds the line length), and unit 30 can not be opened as
c            a new scratch file, then units 29, 28, ..., will be tried
c            until an acceptable unit is found.  Library routines may
c            use this file but must be sure that the use does not
c            conflict with the printing of tables here, or the use by
c            any other library routines.  If K10 is 0, a scratch unit is
c            assumed not to be available, and tables with long lines
c            will be printed with each line on multiple lines.
c MEHEAD=11  (0 .le. K11 .le. 1) Defines the print that surrounds an
c            error message.  K11=0 gives nothing, and 1 gives the first
c            4 characters in TEXT repeated 18 times.  If this is not
c            used, one gets 72 $'s.  (To get a blank line use 1 with
c            TEXT = '    '.)
c MEDDIG=12  (-50 .le. K12 .le. 50) Set default digits to print for
c            floating point.  If K12 > 0 then K12 significant digits
c            will be printed, if K12 < 0, then -K12 digits will be
c            printed after the decimal point, and if K12 = 0, the
c            default will be used, which is the full machine precision.
c            Setting or getting this value will only work properly if
c            the action is taken by calling SMESS or DMESS as
c            appropriate.
c MEMLIN=13  (39 .le. K13 .le. 500) Set message line length to K13.
c            (Default is 128.)
c MEELIN=14  (39 .le. K14 .le. 500) Set error message line length to
c            K14. (Default is 79)
c MEMUNI=15  (-99 .le. K15 .le. 99) Messages go to unit K15.  If K15 = 0
c            (default), 'print' is used.  If K15 < 0, messages go to
c            both 'print' and to unit abs(K15).
c MEEUNI=16  (-99 .le. K16 .le. 99) As for MEMUNI, except for Error
c            Messages.
c MESCRN=17  (0 .le. K17 .le. 100000000) Set number of lines to print to
c            standard output before pausing for "go" from user.  Default
c            is 0, which never stops.
c MEDIAG=18  (0 .le. K18 .le. 1000000000) Set the diagnostic level
c            desired.  This routine makes no use of K18.  It merely
c            serves as a place to set it and to answer inquiries on its
c            value.  It is intended to be set by users of library
c            software.  Library packages that make use of this number
c            are expected to use it as described below.  If K18 = 0 (the
c            default), no diagnostic print is being requested.  Else m =
c            mod(K18,256) determines whether a package will do
c            diagnostic printing.  Associated with a library package is
c            a number L which must be a power of 2 < 129, and which
c            should be mentioned in the documentation for the package.
c            If the bit logical or(m,L) = m then diagnosics output for
c            the routine with the associated value of L is activated.
c            The value of L should have been selected by the following
c            somewhat vague rules.  Let base 2 log(L) = 2*i + j, where j
c            is 0 or 1.  Select i = level of the library package, where
c            the level is 0 if no other library routine that is likely
c            to be used with the package could reasonably be expected to
c            want any embedded diagnostics, and otherwise is
c            min(4, I+1), where I is the maximum level for any library
c            routine which is likely to be used with the package.
c            Select j = 0 if the user is relatively unlikely to want
c            diagnostics, and j = 1, if this is a routine for which
c            considering its level the user is relatively likely to want
c            diagnostic output.  The next 8 bits, mod(K18/256, 256), may
c            be used by the library routine to select the actual output
c            that is to be given.  These bits may be ignored,  but if
c            they are used, the lowest order bits should correspond to
c            less voluminous output that is more likely to be requested.
c            Finally, K18 / (2**16) may be used to give a count on how
c            many times to print the diagnostics that are given.  This
c            count may be interpreted by library routines in slightly
c            different ways, but when used it should serve to turn off
c            all output after a certain limit is reached.  By
c            convention, if this is 0 there is no upper bound on the
c            count.
c MEMAXE=19  (0 .le. K19 .le. 1000000000) Set the maximum error value.
c            When retrieving this value, it is the maximum value seen
c            for 10000*s + 1000*p + i, where s, p, and i are the stop
c            and print levels, and the index on the last error message
c            processed, respectively.  See MEEMES below.
c MESTOP=20  (0 .le. K20 .le. 8) Set the stop level for error messages.
c            If an error message has a stop index > min(K20, 8), the
c            program is stopped after processing the message.  The
c            default value is K20=3.
c MEPRNT=21  (0 .le. K21 .le. 8) Set the print level for error messages.
c            If an error message has a print index > K21, or the message
c            is going to stop when finished, information in an error
c            message is processed, else all the actions including
c            printing are skipped.  (MESTOP controls stopping.)  The
c            default value is MEPRNT = 3.
c An action index of -i, for i < METDIG, will return in the location
c ordinarily used for Ki the current default value for the internal
c variable set by Ki.  In the case of MESUNI, if the scratch unit has
c not been opened, it will be opened before returning the unit number.
c
c METDIG=22  (-50 .le. K22 .le. 50) As for MEDDIG, except the value here
c            is temporary, lasting until the return, or next use of this
c            action.  If 0, the internal value for K12 is used instead.
c MENTXT=23  (1 .le. K23 .le. 10000000) Set value of NTEXT to K23.
c MEIDAT=24  (1 .le. K24 .le. 1000000000) Set value of NIDAT to K24.
c MEFDAT=25  (1 .le. K25 .le. 1000000000) Set value of NFDAT to K25.
c MEMDAT=26  (1 .le. K26 .le. 5) Set value of NMDAT to K26.
c MEMDA1=27  (K27) set MDAT(1) to K27.  See description of NMDAT above.
c MEMDA2=28  (K28) set MDAT(2) to K28.
c MEMDA3=29  (K29) set MDAT(3) to K29.
c MEMDA4=30  (K30) set MDAT(4) to K30.
c MEMDA5=31  (K31) set MDAT(5) to K31.
c METABS=32  (1 .le. K32 .le. 100) set spacing for tabs to K32.
c MECONT=50  Exit, but no print of current print buffer.  The error or
c            diagnostic message is to be continued immediately.
c MERET=51   All done with diagnostic or error message, complete
c            processing and return, or for some error messages stop.
c MEEMES=52  (K52, L52, M52) Start an error message with severity level
c            K52,index for the error of L52, and message text starting
c            at TEXT(M52).  If M52 is 0, message text starts at
c            TEXT(NTEXT), and if M52 < 0, no message text is
c            printed as part of this action.  Library routines should
c            set K52 = 10*s + p, where s is the stop level desired, and
c            p the print level, and should have 10 > p .ge. s .ge. 0.
c            We offer the following guidelines as a yardstick for
c            setting the value of s.
c   = 9  User has ignored warning that program was going to be stopped.
c   = 8  Program has no way to continue.
c   = 7  User has given no indication of knowing that functionality of
c        results is reduced.  (E.g. not enough space for some result.)
c   = 6  Program could continue but with reduced functionality.
c   = 5  Results far worse than user expected to want.
c   = 4  User has given no indication of knowing that results do not
c        meet requested or expected accuracy.
c   = 3  Warning is given that program will be stopped without some
c        kind of response from the calling program.
c   = 2  Program is not delivering requested or expected accuracy.
c   = 1  Some kind of problem that user could correct with more care in
c        coding or in problem formulation.
c   = 0  Message is for information of uncritical nature.
c            Print levels might be counted down so that warnings given
c            several times are no longer given, or be counted up so
c            that a warning is only given after a certain threshold is
c            reached.  Levels should be selected with the understanding
c            that the default is to print only levels above 3.
c METEXT=53  Print TEXT, starting at TEXT(NTEXT).  Print ends
c            with the last character preceding the first '$'.  Special
c            actions are determined by the character following the '$'.
c            Except as noted, the '$' and the single character which
c            follows are not printed.  In the text below, "to continue",
c            means to continue print of TEXT with the next character
c            until the next "$".  Except for the one case noted, NTEXT
c            is set to point to the second character after the "$".
c            Note, letters must be in upper case.  Possibilities are:
c      B  Break text, but don't start a new line.
c      E  End of text and line.
c      R  Break text, don't change the value of NTEXT.  Thus next
c         text Repeats the current.
c      N  Start a New line, and continue.
c      I  Print IDAT(NIDAT), set NIDAT=NIDAT+1, and continue.
c      J  As for I above, except use the last integer format
c         defined by a "$(", see below.
c      F  Print FDAT(NFDAT), set NFDAT=NFDAT+1, and continue.
c      G  As for F above, except use the last floating format
c         defined by a "$(", see below.
c      M  Print MDAT(NMDAT), set NMDAT=NMDAT+1, and continue.
c      H  Marks terminator for column amd row Headings, see table,
c         vector, and matrix output below.  This causes enough blanks to
c         be generated to keep column headings centered over their
c         columns.  After the blanks are generated, text is continued
c         until the next '$'.  This is not to be used except inside
c         column or row headings.  The last row or column should be
c         terminated with a '$E' or if appropriate, a '$#' for a row or
c         column label.
c      (  Starts the definition of a format for integer or floating
c         point output.  The format may not contain a "P" field, and
c         must require no more than 12 characters for floating point
c         (e.g. "(nnEww.ddEe)", where each of the lower case letters
c         represents a single digit), and no more than 7 characters for
c         integer output.  Following the ")" that ends the format, if
c         the next character is not a "$" then "$J" or "$G" type output
c         is done, see above.  In either case processing of TEXT then
c         continues.
c      T  Tab.
c      #  Used in matrix row or column labels this prints the current
c         row or column index, respectively, ends the text for the
c         current row or column, and resets the text pointer to where
c         it started.
c      $  a single '$' is printed, continue till the next '$'.
c      C  Only used by PMESS which deletes it and the preceding '$'.
c         Used at the end of a line to indicate continued text.
c   other Don't use this -- the '$' is ignored, but new features may
c         change the action.  (E.g. $P might be added to get a prompt.)
c ME????=54  Not used.
c METABL=55  (K55, L55, M55, N55)  Note this action automatically
c            returns when done, further locations in MACT are not
c            examined.  This action prints a heading and/or data that
c            follows a heading.  If K55 is 1, then the heading text
c            starting in TEXT(NTEXT) is printed.  This text
c            should contain embedded "$H"'s to terminate columns of the
c            heading.  If there is no heading on a column, use " $H".
c            Note the leading blank.  If the heading is to continue
c            over k columns, begin the text with "$H" repeated k-1
c            times with no other embedded characters.  The very last
c            column must be terminated with "$E" rather than "$H".
c            After tabular data is printed, K55 is incremented by 1,
c            and compared with L55.  If K55 > L55, K55 is reset to 1,
c            and if the data that was to be printed had lines that were
c            too long, data saved in the scratch file is printed using
c            the headings for the columns that would not fit on the
c            first pass.  Note that only one line of tabular data can
c            be printed on one call to this subroutine.
c            M55 gives the number of columns of data associated with the
c            heading.
c            N55 is a vector containing M55 entries.  The k-th integer
c            in N55 defines the printing action for the k-th column
c            of the table.  Let such an integer have a value defined by
c            rr + 100 * (t + 10 * (dd + 100 * ww)), i.e. wwddtrr, where
c            0 .le. rr,dd,ww < 100, and 0 .le. t < 10.
c      rr    The number of items to print.
c      t     The type of output.
c            1  Print text starting at TEXT(NTEXT), rr = 01.
c            2  Print the value of K55, rr = 01.
c            3  Print integers starting at IDAT(NIDAT).
c            4  Print starting at FDAT(NFDAT), using an F format.
c            5  Print starting at FDAT(NFDAT), using an E format.
c            6  Print starting at FDAT(NFDAT), using an G format.
c      dd    Number of digits after the decimal point.
c      ww    The total number of column positions used by the column,
c            including the space used to separate this column from the
c            preceding one.  This must be big enough so that the column
c            headings will fit without overlap.
c MEIVEC=57  (K57) Print IDAT as a vector with K57 entries.  The vector
c            output starts on the current line even if the current line
c            contains text.  This is useful for labeling the vector.
c            The vector starts at IDAT(NIDAT).
c            If K57 < 0,  indices printed in the labels for the vector
c            start at at NIDAT, and entries from NIDAT to -K57 are
c            printed.
c MEIMAT=58  (K58, L58, M58, I58, J58) Print IDAT as a matrix with K58
c            declared rows, L58 actual rows, and M58 columns.  If L58<0,
c            the number of actual rows is mod(-L58, 100000), and the
c            starting row index is -L58 / 100000.  Similarly for M58<0.
c            TEXT(I58) starts the text for printing row labels.  If
c            I58 < 0, no row labels are printed.  If I58 = 0, it is as
c            if it pointed to text containing "Row $E".  Any "$" in a
c            row or column label must be followed by "H" or
c            "E" which terminates the text for the label.  In the case
c            of $H, text for the next label follows immediately,
c            in the case of $E the current row index is printed in place
c            of the $E and the next label uses the same text.  J58 is
c            treated similarly to I58, except for column labels, and
c            with "Row $E" replaced with "Col $E".  The matrix starts at
c            IDAT(NIDAT), and NIDAT points one past the end of the
c            matrix when finished.
c MEJVEC=59  (K59) As for MEIVEC, except use format set with $(.
c MEJMAT=60  (K60, L60, M60, I60, J60) As for MEIMAT, except use the
c            format set with $(.
c MEFVEC=61  (K61) As for MEIVEC, except print FDAT as a vector with
c            K61 entries.  The vector starts at FDAT(NFDAT).
c MEFMAT=62  (K62, L62, M62, I62, J62) As for action MEIMAT, but
c            instead print FDAT, and use NFDAT in place of NIDAT.
c MEGVEC=63  (K63) As for MEFVEC, except use format set with $(.
c MEGMAT=64  (K64, L64, M64, I64, J64) As for MEIMAT, except use the
c            format set with $(.
c
c++ CODE for VERSION = C is inactive
c      parameter (LNMSG = 128)
c      parameter (LNERR = 79)
c      integer  kciwid, kccwid, kcrwid, lfprec, lgprec
c      common /MESSCC/ kciwid, kccwid, kcrwid, lfprec, lgprec
c      integer kc
c++ END
c
c ************** Parameters Defining Actions (See Above) ***************
c
      integer   MESUNI, MEHEAD, MEDDIG, MEMLIN, MEELIN, MEMUNI, MEEUNI,
     1  MESCRN, MEDIAG, MEMAXE, MESTOP, MEPRNT, METDIG, MENTXT, MEIDAT,
     2  MEFDAT, MEMDAT, MEMDA1, MEMDA2, MEMDA3, MEMDA4, MEMDA5, METABS,
     3  MECONT, MERET , MEEMES, METEXT, METABL, MERES3, MEIVEC,
     4  MEIMAT, MEJVEC, MEJMAT, MEFVEC, MEFMAT, MEGVEC, MEGMAT, MEMAXI,
     5  MEGBAS, MEVBAS, MEVLAS
c Parameters for changing the environment.
      parameter (MESUNI=10,MEHEAD=11,MEDDIG=12,MEMLIN=13,MEELIN=14,
     1 MEMUNI=15,MEEUNI=16,MESCRN=17,MEDIAG=18,MEMAXE=19,MESTOP=20,
     2 MEPRNT=21,METDIG=22,MENTXT=23,MEIDAT=24,MEFDAT=25,MEMDAT=26,
     3 MEMDA1=27,MEMDA2=28,MEMDA3=29,MEMDA4=30,MEMDA5=31,METABS=32)
c Parameters for actions.
      parameter (MECONT=50, MERET=51,MEEMES=52,METEXT=53,
     1 METABL=55,MERES3=56,MEIVEC=57,MEIMAT=58,MEJVEC=59,MEJMAT=60,
     2 MEFVEC=61,MEFMAT=62,MEGVEC=63,MEGMAT=64)
c Parameter derived from those above.
      parameter (MEMAXI=64,MEGBAS=49,MEVBAS=10,MEVLAS=32)
c
c ************************** Internal Variables ************************
c
c BUF    Character string holding characters to be output.
c C      Used for temp. storage of a character.
c DOLS   A heading/trailing string of characters, default = $'s.
c ERMSG  Default part of error message.
c EUNIT  Unit number for output of error messages.
c FDAT   Formal array, containing floating point data to output.  Only
c   appears external to this subroutine.
c FMTC   Format for integer output of matrix column headings.
c FMTF   Format for floating point or other output.
c FMTG   Format set by user for floating point or other output.
c FMTI   Character string holding format for integer output.
c FMTIM  Equivalenced to FMTR, FMTC.
c FMTJ   Format set by user for integer output.
c FMTR   Value of FMTI for format of row indices.
c FMTT   Format to be stored in FMTJ or FMTG.
c GETW   Set true if still need to get width for a format.
c GOTFMT Set .true. if format has been set by user.  When printing
c   tables, set true when heading has been output.
c I      Index of current action from MACT.
c ICHAR0 Value of ICHAR('0')
c ICOL   Current column index in matrix output.
c IDAT   Formal array, containing integer data to output.
c IMAG   Magnitude of integer to output, with negative sign if integer
c   is < 0.
c INCM   Array giving amount of space used by the options.
c INERR  0 if not processing an error message, 1 if printing an error
c   message, -1 if in an error message that is not being printed, and >1
c   if printing an error message that stops.  Set to -2 when the error
c   message is supposed to stop.
c IOUT   Integer to be output.
c IRC    = 1 for rows, = 2 for columns when determining labels for
c   matrix output.
c IROW   Row index for matrix output.  Also used in table output to
c    count lines for printing line index on read from scratch unit.
c IROW1  Starting row index for matrix output.
c ITEXT  Index of the element of TEXT use for the next text output.
c ITXTSV Saved value of NTEXT when doing matrix output.
c IVAR   Integer array, that is equivalenced to a number of integer
c   variables that can be set by the user.
c IWF    Width to be used in a floating pt. (or other) format.
c IWG    Value like IWF for user set format.
c J      Used as a temporary index.
c JJ      Used as a temporary index.
c K      Used as a temporary index.
c K      Used as a temporary index.
c K1     Used as a temporary index.
c K2     Used as a temporary index.
c KDF    Current number of digits to print for floating point.  See
c   description of MACT(*) = METDIG above.
c KDFDEF Current default for KDF, see description of MACT(*) = MEDDIG.
c KDI    Number of digits used to print last integer.
c KDIAG  Not directly referenced.  Holds place in IVAR for reference
c   from outside.  See comments above.
c KDILAB Length for printing index in vector output.
c KDJ    As for KDI, except for format set by user.
c KK     Temporary index.
c KLINE  Count of number of things to print on a line.  (In table print
c   is the number to print for one spec field.)
c KNT    In vector output gives the current index for output.
c KOLWID Length to use for printing a column heading.  Equivalenced to
c   MAXWID(2).
c KP     Index from error action input for the print action.
c KRES1  Holds place in common block for future use.
c KS     Index from error action input for the stop action.
c KSCRN  Number of lines to "print" before pausing.
c KSHIFT Amount to shift column heading before printing.
c KSPEC  Defines action after looking for character in TEXT.  (Also
c   used as a temporary index.)
c   1.  $B   Break the text here continue on same line.
c   2.  $E   Break text, print what is in BUF.
c   3.  $R   Break text, continue on same line, NTEXT set to repeat the
c            current text.
c   4.  $N   Print BUF, continue with following text.
c   5.  $I   Print IDAT(NIDAT), continue TEXT.
c   6.  $F   Print FDAT(NFDAT), continue TEXT.
c   7.  $M   Print MDAT(NMDAT), continue TEXT.
c   8.  $J   As for $I, but with user format.
c   9.  $G   As for $F, but with user format.
c  10.  $(   Set a user format.
c  11.  $T   Tab.
c  12.       Set when done with an action.
c  13.       Set when done with boiler plate text for an error message.
c   0. Other Ignore the "$", continue with TEXT.
c KT     Used for logic in output of headings.
c        = 1 Output table headings.
c        = 2 Get row/column widths for matrix output.
c        = 3 Output column headings for matrix output.
c LASKNT In vector output value index for last element to print.
c LASTI  Last index for matrix output, or for finding values that
c   determine format.
c LBUF   Position of characters in BUF, usually the last to print.
c LBUF1  Start of text to shift when shifting text in BUF to the left.
c LBUF2  End of text to shift when shifting text in BUF to the left.
c LENBUF Parameter giving the number of character in BUF.
c LENLIN Gives number of character in output lines.
c LENOUT Length of output for table or vector/matrix output.
c LENTXT Length of character array elements in TEXT.
c LENTRY Tells what to do on entry (and sometimes other places.)
c   = 1  The value on first entry.
c   = 2  A previous entry is to be continued.
c   = 3  A non printing error message is to be continued
c   = 4  Just done output from inside a METEXT action.
c   = 5  Got "maximum" value for entries in a vector.
c   = 6  Got "maximum" value for entries in a matrix.
c   = 7  Vector (either print or get format for label indices.)
c   = 8  Matrix (either print or get format for label indices.)
c   = 9  Output of data in a table.
c LHEAD  If = 0 no print of DOLS, else DOLS printed in error messages.
c LINERR Gives LENLIN for error messages.
c LINMSG Gives LENLIN for diagnostic messages.
c LINSTR Space at start of line for label in vector and matrix output.
c   Equivalenced to MAXWID(1).
c LNERR  Parameter giving the default value for LENERR, only in XMESS.
c LNMSG  Parameter giving the default value for LINMSG, only in XMESS.
c LOCBEG Index of first item in vector and matrix output.
c LPRINT For error messages with a print level .le. LPRINT nothing is
c   printed (unless the message would result in a stop).
c LSTOP  As for LPRINT, except with the stop level, and stopping.
c LSTRT  Starting character position in BUF for storing next characters.
c LTEXT  Length of heading text in TEXT.
c M      Index for the current action.
c MACT   Formal integer array specifying the actions, see above.
c MAXERR Value save in IVAR for user to get indication of last (really
c   the maximum error seen so far = 1000 * (10*stop + print) + index.
c MAXWID Equivalenced to  LINSTR and KOLWID.
c MBNDHI Upper bounds for inputs to IVAR.
c MBNDLO Lower bounds for inputs to IVAR.
c MDAT   Array where user can store integers in IVAR for later output.
c   Also used to store column indices when tables are saved on scratch
c   unit.
c
c The following parameter names starting with ME define actions
c   which the user can request.  They have been documented in the
c   description above, except for the ones defined just below.
c MEGBAS is 1 less than the smallest action that involves something
c   other than just storing or retrieving a value.
c MEMAXI is the largest action which a user can request.
c MEVBAS is the smallest action index, used to set the starting index in
c   IVAR.
c MEVLAS is the largest index for a variable in IVAR.
c MECONT,  MEDDI,  MEELI,  MEEME, MEEUNI, MEFDAT, MEFMAT, MEFVEC,
c MEGBAS, MEGMAT, MEGVEC, MEHEAD, MEIDAT, MEIMAT, MEIVEC, MEJMAT,
c MEJVEC, MEMAXE, MEMAXI, MEMDA1, MEMDA2, MEMDA3, MEMDA4, MEMDA5,
c MEMDAT, MEMLIN, MEMUNI, MENTXT, MEPRNT, MESCRN, MERES1, MERES2,
c MERES3, MERET, MESTOP, MESUNI, METABL, METDIG, METEXT
c MPT    Current pointer to data for matrix or vector output.
c MTEXT  Equivalenced to MTEXTR and MTEXTC.
c MTEXTC TEXT(MTEXTC) starts text for printing column labels.
c MTEXTR TEXT(MTEXTR) starts text for printing row labels.
c MUNIT  Output unit used for messages that aren't in an error message.
c NCOL   Number of columns for matrix output, 0 for vector output,
c   count of column left for table output.
c NDIM   Distance between columns for matrix output.
c NFDAT  Index of next item in FDAT to print.
c NIDAT  Index of next item in IDAT to print.
c NLINE  Maximum number of data items to print on a line for matrix and
c   vector output.  If the scratch file is used for table output,
c   NLINE gives the original end of the buffer.
c NMDAT  Pointer to next thing to print from MDAT.
c NROCO  Equivalenced to (NROW, NCOL).  Used in matrix output.
c NROW   Number of rows for matrix output.    When printing tables,
c   MDAT(NROW) gives place where line was split.  (=0 if not split)
c NTEXT  Index inside an element of TEXT for the next text output.
c NTEXTR Value of NTEXT to use if get a $R.
c NTXTSV Saved value of NTEXT when doing matrix output.
c OUNIT  Index of the current output unit.
c SC     Parameter for special character used to introduce actions.
c   Default value is '$'.  If this is changed the "$"'s in comments
c   should be changed to the new value of the character.  (Note that
c   SC = '\' is not portable.)
c SCRNAM Name of file constructed for error output or message output.
c SUNIT  Index for the scratch unit, -1 if not yet assigned.
c TEXT   Formal argument giving the character string from which all text
c   is taken.
c UMESS  Name of subroutine called that does nothing, but which may be
c   modified by the user to cause different actions to be taken.
c XARG   If .true., output data is not integer, and a return is made to
c   print data from FDAT.
c XARGOK Set .true. if call is from program that will print data from
c   FDAT.
c XMESS  Name of the block data subprogram that initializes common data.
c
c ************************** Variable Declarations *********************
c
      integer    MACT(*), IDAT(*)
      character  TEXT(*)*(*)
c
      integer    I, ICOL, INCM(MECONT:MEIMAT), INERR, IOUT, IROW, IROW1,
     1    ITEXTR, ITXTSV, J, JJ, K, K1, K2, KDILAB, KK, KNT, KOLWID, KP,
     2    KS, LASKNT, LBUF1, LBUF2, LENBUF, LINSTR, M,
     3    MBNDHI(MEVBAS:MEVLAS-5), MBNDLO(MEVBAS:MEVLAS-5), MTEXT(2),
     4    MTEXTC, MTEXTR, NLINE, NROCO(2), NTEXTR, NTXTSV
      integer MESSGS
      logical   GETW
      character ERMSG*63, ERMSG1*27
      character SC, C, SCRNAM*12
      parameter (SC = '$')
      save  I, ICOL, INERR, IROW, IROW1, ITXTSV, KDILAB, KNT, LASKNT, M,
     1   MTEXT, NLINE, NTXTSV
C     COMMENTED OUT BY KUMAZAWA
C     save /CMESSI/, /CMESSC/
      equivalence (MTEXT(1), MTEXTR), (MTEXT(2), MTEXTC)
c
c ************************** Data from common block ********************
c
c++ CODE for VERSION = F is active
      external XMESS
c++ END
      parameter (LENBUF = 250)
      logical          XARG, GOTFMT, XARGOK
      integer          EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS), IMAG,
     1   ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ, KLINE, KSCRN,
     2   KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN, LENOUT,
     3   LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT, LSTOP,
     4   LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL, NDIM,
     5   NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
c
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*15, FMTG*15,
     1  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,
     1   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT,
     2   NFDAT, NMDAT, MDAT, TABSPA, ICHAR0, IMAG, IRC, ITEXT, IWF,
     3   IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI, LBUF, LENLIN,
     4   LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT, MAXWID, MPT,
     5   NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
c
      equivalence (NROCO, NROW)
      equivalence (MAXWID(1), LINSTR), (MAXWID(2), KOLWID)
c ************************** End of stuff from common block ************
c
      data INERR / 0 /
      data ERMSG /
     1' reports error: Stop level = x, Print level = y, Error index = '/
      data ERMSG1 / ': Print level = y, Index = ' /
c                 50  51, 52  53  54 55 56 57 58
      data INCM /  1,  1,  4,  1,  2, 0, 0, 2, 6 /
      data MBNDLO /  0, 0, -50,  39,  39, -99, -99,         0,
     1           0,          0, 0, 0, -50,        1,          1,
     2           1, 1,   1 /
      data MBNDHI / 99, 1,  50, 500, 500,  99,  99, 100000000,
     1  1000000000, 1000000000, 8, 8,  50, 10000000, 1000000000,
     2  1000000000, 5, 100 /
c
c ************************* Start of Executable Code *******************
c
c
      go to (5, 10, 20, 850, 1160, 1620, 1130, 1530, 960), LENTRY
      ICHAR0 = ICHAR('0')
c++ CODE for VERSION = C is inactive
c      LSTOP = 3
c      LHEAD = 1
c      LPRINT = 3
c      LINMSG = LNMSG
c      LINERR = LNERR
c      LENLIN = LNMSG
c      SUNIT = -1
c      TABSPA = 6
c      KDJ = 7
c%%    for (kc=0; kc < 71; kc++)  cmessc.dols[kc] = '$';
c      FMTI = '%*d'
c      FMTJ = '%*d\0'
c      FMTG = '%*.*E\0'
c++ END
c                             First entry for a message
    5 LBUF = 0
c                             Usual continuation entry
   10 I = 1
      NTEXT = 1
      ITEXT = 1
      LENTXT = len(TEXT(1))
      NIDAT = 1
      NFDAT = 1
      NMDAT = 1
      go to 120
c                     Continuation entry when have a non printing error
c Skip all actions -- Inside non-printing error message.
   20 I = 1
   30 K = MACT(I)
      if (K .le. MERET) then
         if (K .eq. MERET) go to 120
         if (K .eq. MECONT) return
         if (K .le. -MEGBAS) go to 180
         I = I + 2
      else
         if (K .gt. MEIMAT) then
            if (K .gt. MEMAXI) go to 180
            K = MEIVEC + mod(K - MEIVEC, 2)
         end if
         I = I + INCM(K)
      end if
      go to 30
c
c Print BUF
   40 call MESSPR
c                             Usual place to end an action request.
  100 I = I + INCM(M)
c                             Pick up the next action request
  120 M = MACT(I)
      if (M .gt. MEGBAS) go to 140
      I = I + 2
      if (abs(M) .gt. MEVLAS) go to 180
      if (M .gt. 0) then
         IVAR(M) = MACT(I-1)
         if (M .ge. MEMDA1) then
            if (M .le. MEMDA5) go to 120
            M = M - 5
         end if
         if (IVAR(M) .lt. MBNDLO(M)) then
            IVAR(M) = MBNDLO(M)
         else if (IVAR(M) .gt. MBNDHI(M)) then
            IVAR(M) = MBNDHI(M)
         end if
c            MEHEAD, MEDDIG, MEMLIN, MEELIN, MEMUNI, MEEUNI
         go to (122,    124,    126,    126,    128,    128), M - MESUNI
         if (M .ne. MENTXT) go to 120
         ITEXT = (NTEXT-1) / LENTXT
         NTEXT = NTEXT - LENTXT*ITEXT
         ITEXT = ITEXT + 1
         go to 120
  122    if (LHEAD .ne. 0) then
         end if
         go to 120
  124    KDF = KDFDEF
         go to 120
  126    LENLIN = LINMSG
         go to 120
  128    if (IVAR(M) .ne. 0) then
c++ CODE for VERSION = F is active
            K = abs(IVAR(M))
            inquire (K, OPENED=XARG)
            if (XARG) then
               open (UNIT=K, STATUS='UNKNOWN')
            else
               write(SCRNAM, '(A, I2.2, A)') 'MESSF_', K, '.tmp'
               open (UNIT=K, STATUS='NEW', FILE=SCRNAM)
            end if
c++ CODE for VERSION = C is inactive
c%%          c_fname[m-15][6] = k / 10 + '0';
c%%          c_fname[m-15][7] = k % 10 + '0';
c%%          if (strcmp(&c_fname[16-m][6], &c_fname[m-15][6]))
c%%             c_handle[m-15] = fopen(c_fname[m-15],"w");
c%%          else
c%%             c_handle[m-15] = c_handle[16-m];
c++ END
         end if
         OUNIT = MUNIT
         go to 120
      end if
      if (M .eq. -MESUNI) then
c++ CODE for VERSION = F is active
         if (SUNIT .le. 0) SUNIT = MESSGS()
c++ CODE for VERSION = C is inactive
c%%      if (cmessi.sunit == -1L) {
c%%          scratch_file = tmpfile();
c%%          cmessi.sunit = 1L;}
c++ END
      end if
C
      MACT(I-1) = IVAR(-M)
      go to 120
c  ME ..    CONT  RET EMES ETXT ???? TABL
  140 go to (170, 200, 310, 400, 180, 910, 180), M-MEGBAS
      if (M .le. MEGMAT) go to 1000
      go to 180
c
c Action MECONT -- Continue message on next entry
  170 LENTRY = 2
      return
c
c Some kind of error in message specification.
  180 continue
c++ CODE for VERSION = F is active
      BUF(1:57) =
     1   'Actions in MESS terminated due to error in usage of MESS.'
c++ CODE for VERSION = C is inactive
c%%   memcpy(cmessc.buf,
c%%   "Actions in MESS terminated due to error in usage of MESS.",57);
c++ END
      LBUF = 57
c
c Action MERET -- Finish a message.
  200 LENTRY = 1
      J = INERR
      INERR = 0
      if (J .ge. 2) INERR = -2
      if (J .gt. 0) go to 330
c                       Finish print before exit.
      call MESSPR
      return
c
c Action MEEMES -- Start an error message
  310 LENTRY = 3
      call UMESS(TEXT, MACT(I+1), IVAR)
      IMAG = max( 0, min(999, MACT(I+2)))
      K = MACT(I+1)
      MAXERR = max(MAXERR, 1000*K + IMAG)
      KS = K / 10
      KP = K - 10 * KS
      if (KS .le. min(LSTOP, 8)) then
         if (KP .le. LPRINT) then
            INERR = -1
            go to 20
         end if
         INERR = 1
      else
         INERR = 2
      end if
      OUNIT = EUNIT
      LENLIN = LINERR
c                        Output a blank line.
      BUF(1:1) = ' '
      LBUF = 1
  330 call MESSPR
c                        Put out line of $'s
      if (LHEAD .ne. 0) then
         LBUF = min(len(DOLS), LENLIN)
         BUF = DOLS
c++ CODE for VERSION = F is active
         if (INERR.lt.0) BUF(5:37)=' Fatal error -- Program stopped. '
c++ CODE for VERSION = C is inactive
c%%      if (inerr < 0L)
c%%       memcpy(&cmessc.buf[4]," Fatal error -- Program stopped. ",34);
c++ END
         call MESSPR
      end if
      if (INERR .le. 0) then
c                                 Just finished an error message
         if (INERR .ne. 0) stop 'Stopped in MESS.'
         OUNIT = MUNIT
         LENLIN = LINMSG
         return
      end if
c                     Just starting an error message get program name
      NTEXTR = 0
      go to 410
c                     Got the program name in BUF.
  370 LBUF = min(LBUF, 40)
      if (KS .eq. 0) then
         ERMSG1(17:17) = char(KP + ICHAR0)
c++ CODE for VERSION = F is active
         BUF(LBUF+1:LBUF+len(ERMSG1)) = ERMSG1
c++ CODE for VERSION = C is inactive
c%%       memcpy(&cmessc.buf[cmessi.lbuf], ermsg1, strlen(ermsg1));
c++ END
         LBUF = LBUF + len(ERMSG1)
      else
         ERMSG(30:30) = char(KS + ICHAR0)
         ERMSG(47:47) = char(KP + ICHAR0)
c++ CODE for VERSION = F is active
         BUF(LBUF+1:LBUF+len(ERMSG)) = ERMSG
c++ CODE for VERSION = C is inactive
c%%       memcpy(&cmessc.buf[cmessi.lbuf], ermsg, strlen(ermsg));
c++ END
         LBUF = LBUF + len(ERMSG)
      end if
      LSTRT = LBUF + 1
      call MESSFI
      LBUF = LBUF + KDI
c++ CODE for VERSION = F is active
      write (BUF(LSTRT:LBUF), FMTI) IMAG
c++ CODE for VERSION = C is inactive
c%%   sprintf(&cmessc.buf[cmessi.lstrt-1L], "%*d",
c%%           messcc.kciwid, cmessi.imag);
c++ END
c          Finish up the start error message action.
      if (MACT(I+3) .lt. 0) go to 40
      if (MACT(I+3) .ne. 0) then
         ITEXT = (MACT(I+3)-1) / LENTXT
         NTEXT = MACT(I+3) - LENTXT*ITEXT
         ITEXT = ITEXT + 1
      end if
      KSPEC = 13
      go to 480
c                  Take care of any left over print from error header
  390 if (LBUF .ne. 0) call MESSPR
c
c Action METEXT -- Print string from TEXT
  400 LENTRY = 4
      NTEXTR = NTEXT
      ITEXTR = ITEXT
c                  Continue with print from TEXT
  410 K1 = LENBUF - LBUF + NTEXT
      K2 = min(K1, LENTXT)
      LSTRT = LBUF + 1
c++ CODE for VERSION = F is active
      K = index(TEXT(ITEXT)(NTEXT:K2), SC)
c++ CODE for VERSION = C is inactive
c%%       if ((ctmp=memchr(TEXT(cmessi.itext-1L,cmessi.ntext-1), SC,
c%%          k2 - cmessi.ntext + 1)) == NULL)
c%%             k = 0;
c%%       else
c%%             k = ctmp - TEXT(cmessi.itext-1L,cmessi.ntext-1) + 1;
c++ END
      if (K .eq. 0) then
         LBUF = LSTRT + K2 - NTEXT
c++ CODE for VERSION = F is active
         BUF(LSTRT:LBUF) = TEXT(ITEXT)(NTEXT:K2)
c++ CODE for VERSION = C is inactive
c%%       memcpy(&cmessc.buf[cmessi.lstrt-1L], TEXT(cmessi.itext-1L,
c%%         cmessi.ntext-1), k2 - cmessi.ntext + 1L);
c++ END
         if (K1 .ge. LENTXT) then
            ITEXT = ITEXT + 1
            NTEXT = 1
            if (LBUF .le. LENLIN) go to 410
            K1 = 1
         end if
         NTEXT = NTEXT + K1 - 1
         KSPEC = 12
         if (ITEXT - ITEXTR .lt. 4000) go to 480
         KSPEC = 2
         go to 430
      end if
      LBUF = LBUF + K - 1
c++ CODE for VERSION = F is active
      if (K .ge. 2) BUF(LSTRT:LBUF) = TEXT(ITEXT)(NTEXT:NTEXT+K-2)
c++ CODE for VERSION = C is inactive
c%%   if (k >= 2) memcpy(&cmessc.buf[cmessi.lstrt-1],
c%%     TEXT(cmessi.itext-1L, cmessi.ntext-1), k - 1L);
c++ END
      NTEXT = NTEXT + K + 1
      if (NTEXT .gt. LENTXT) then
         ITEXT = ITEXT + 1
         if (NTEXT .eq. LENTXT + 1) then
            C = TEXT(ITEXT-1)(LENTXT:LENTXT)
            NTEXT = 1
         else
            C = TEXT(ITEXT)(1:1)
            NTEXT = 2
         end if
      else
         C = TEXT(ITEXT)(NTEXT-1:NTEXT-1)
      end if
      if (NTEXTR .eq. 0) then
         if (LENTRY .eq. 3) go to 370
         go to 1510
      end if
      KSPEC = index('BERNIMFJG(T', C)
  430 if (LBUF .gt. LENLIN) go to 480
c              1   2   3   4   5   6   7   8   9  10  11  12, 13
c              B   E   R   N   I   M   F   J   G   (   T done end err
      go to (455,480,450,460,700,680,900,700,900,600,690,410,390), KSPEC
c                        No match, continue with the text.
  440 LBUF = LBUF + 1
      BUF(LBUF:LBUF) = C
      go to 410
c                        Reset NTEXT for $R
  450 NTEXT = NTEXTR
      ITEXT = ITEXTR
c                             Done with METEXT action.
  455 NMDAT = 1
      go to 100
c           At this point want to output all in BUF
  460 do 470 LBUF = LBUF, 1, -1
         if (BUF(LBUF:LBUF) .ne. ' ') go to 480
  470 continue
  480 LBUF2 = LBUF
      if (LBUF2 .eq. 0) then
         LBUF = 1
         BUF(1:1) = ' '
      else if (LBUF .gt. LENLIN) then
         do 485 K = LENLIN+1, LENLIN/3, -1
            if (BUF(K:K) .eq. ' ') then
               LBUF = K - 1
               go to 490
            end if
  485    continue
         LBUF = LENLIN
      end if
  490 LBUF1 = LBUF
      call MESSPR
      if (LBUF1 .ge. LBUF2) then
c                       The entire buffer has been printed.
         if (KSPEC .le. 2) go to 455
         if (KSPEC .ne. 4) go to 430
         go to 410
      end if
c                       Remove trailing blanks
      do 510 LBUF1 = LBUF1+1, LBUF2
         if (BUF(LBUF1:LBUF1) .ne. ' ') go to 520
  510 continue
c                       Shift the contents of the buffer.
  520 LBUF = LBUF2-LBUF1+1
      LSTRT = 1
  530 if (LBUF .ge. LBUF1) then
c                              Take care of overlap.
         K = 2*LBUF1 - LSTRT
c++ CODE for VERSION = F is active
         BUF(LSTRT:LBUF1-1) = BUF(LBUF1:K-1)
c++ CODE for VERSION = C is inactive
c%%  memcpy(&cmessc.buf[cmessi.lstrt-1],&cmessc.buf[lbuf1-1L],k-lbuf1);
c++ END
         LSTRT = LBUF1
         LBUF1 = K
         go to 530
      end if
c++ CODE for VERSION = F is active
      if (LBUF .ge. LSTRT) BUF(LSTRT:LBUF) = BUF(LBUF1:LBUF2)
c++ CODE for VERSION = C is inactive
c%%  if (cmessi.lbuf>=cmessi.lstrt) memcpy(&cmessc.buf[cmessi.lstrt-1],
c%%       &cmessc.buf[lbuf1-1L], lbuf2-lbuf1+1);
c++ END
      go to 430
c
c Get information on user format
  600 KSPEC = 8
c              I,   i,   F,   f,   E,   e,   G,   g
      go to (604, 604, 601, 601, 602, 602, 602, 602),
     1   index('IiFfEeGf',TEXT(ITEXT)(NTEXT:NTEXT))
      go to 180
  601 continue
c++ CODE for VERSION = F is active
      FMTG='(   99F'
      go to 603
c++ END
  602 continue
c++ CODE for VERSION = F is active
      FMTG='(1P,99'//TEXT(ITEXT)(NTEXT:NTEXT)
c++ END
  603 KSPEC = 9
  604 IMAG = 0
      GETW = .true.
c++ CODE for VERSION = C is inactive
c%%   strcpy(cmessc.fmtg, "%*.*E\0");
c      FMTG(5:5) = TEXT(ITEXT)(NTEXT:NTEXT)
c%%   messcc.lgprec = 0;
c++ END
      K = NTEXT
  606 continue
         NTEXT = NTEXT + 1
         if (NTEXT .gt. LENTXT) then
            ITEXT = ITEXT + 1
            NTEXT = 1
         end if
c++ CODE for VERSION = F is active
         FMTT(NTEXT-K:NTEXT-K) = TEXT(ITEXT)(NTEXT:NTEXT)
c++ END
         JJ = ichar(TEXT(ITEXT)(NTEXT:NTEXT)) - ICHAR0
         if (GETW) then
            if ((JJ .ge. 0) .and. (JJ .le. 9)) then
               IMAG = 10*IMAG + JJ
            else
               if (TEXT(ITEXT)(NTEXT:NTEXT) .eq. ')')  go to 610
               if (TEXT(ITEXT)(NTEXT:NTEXT) .ne. '.')  go to 180
               GETW = .false.
            end if
         else
            if (TEXT(ITEXT)(NTEXT:NTEXT) .eq. ')') go to 610
            if ((JJ .lt. 0) .or. (JJ .gt. 9)) go to 180
c++ CODE for VERSION = C is inactive
c%%         messcc.lgprec = 10*messcc.lgprec + jj;
c++ END
         end if
      go to 606
c
  610 NTEXT = NTEXT + 1
      if (NTEXT .gt. LENTXT) then
         ITEXT = ITEXT + 1
         NTEXT = 1
      end if
c++ CODE for VERSION = F is active
      if (KSPEC .eq. 8) then
         KDJ = IMAG
         FMTJ(5:7) = FMTT
      else
         IWG = IMAG
         FMTG(8:15) = FMTT
      end if
c++ CODE for VERSION = C is inactive
c%%   if (cmessi.kspec == 8)
c%%       cmessi.kdj = cmessi.imag;
c%%   else
c%%       cmessi.iwg = cmessi.imag;
c++ END
      if (TEXT(ITEXT)(NTEXT:NTEXT) .eq. SC) go to 410
      if (KSPEC .eq. 8) go to 700
      if (XARGOK) return
      go to 440
c
c                         Print from MDAT
  680 IOUT = MDAT(NMDAT)
      if (NMDAT .ge. 6) then
         MDAT(NMDAT) = MDAT(NMDAT) + 1
      else
         NMDAT = NMDAT + 1
      end if
      go to 720
c
c                         Process a tab
  690 LSTRT = LBUF + 1
      LBUF = min(LBUF + TABSPA - mod(LBUF, TABSPA), LENLIN+1)
c++ CODE for VERSION = F is active
      BUF(LSTRT:LBUF) = ' '
c++ CODE for VERSION = C is inactive
c%%   for (kc=cmessi.lstrt-1; kc<cmessi.lbuf; kc++) cmessc.buf[kc]=' ';
c++ END
      go to 850
c                         Print from IDAT
  700 IOUT = IDAT(NIDAT)
      NIDAT = NIDAT + 1
  720 LSTRT = LBUF + 1
      IMAG = IOUT
      if (KSPEC .ge. 8) then
         LBUF = LBUF + KDJ
c++ CODE for VERSION = F is active
      write (BUF(LSTRT:LBUF), FMTJ) IOUT
c++ CODE for VERSION = C is inactive
c%%   sprintf(&cmessc.buf[cmessi.lstrt-1], "%*d", cmessi.kdj, iout);
c++ END
         go to 850
      end if
c
c                Get format for integer output.
  800 call MESSFI
      LBUF = LBUF + KDI
c++ CODE for VERSION = F is active
      write (BUF(LSTRT:LBUF), FMTI) IOUT
c++ CODE for VERSION = C is inactive
c%%  sprintf(&cmessc.buf[cmessi.lstrt-1], "%*d", messcc.kciwid, iout);
c++ END
c                         Entry here to check line after numeric output.
  850 if (LBUF .le. LENLIN) go to 410
      KSPEC = 12
      go to 480
c
c                          Take care of output for extra argument.
  900 if (XARGOK) return
      go to 180
c
c Action METABL -- Start a table
  910 GOTFMT = MACT(I+1) .ne. 1
      if (.not. GOTFMT) then
         IROW = 0
         KOLWID = 0
      end if
      LENTRY = 9
      if (LBUF .ne. 0) call MESSPR
  920 BUF = ' '
      NROW = 1
      NCOL = MACT(I+3)
      ICOL = I + 3
  940 ICOL = ICOL + 1
      JJ = MACT(ICOL)
      KLINE = mod(JJ, 100)
      LENOUT = JJ / 100000
      NCOL = NCOL - max(KLINE, 1)
      if (GOTFMT) then
c                                 Print the data
         LSTRT = LBUF + 1
         LBUF = min(LBUF + KLINE * LENOUT, LENBUF)
         JJ = JJ / 100
         KK = mod(JJ, 10)
c              Text,   I   I',   F    E    G
         go to (948, 941, 941, 943, 945, 944), KK
         go to 180
c                             Integer output
  941    continue
c++ CODE for VERSION = F is active
         FMTI(5:5) = char(LENOUT / 10 + ichar0)
         FMTI(6:6) = char(mod(LENOUT, 10) + ichar0)
c++ END
         if (KK .eq. 3) then
c++ CODE for VERSION = F is active
            write (BUF(LSTRT:LBUF), FMTI) MACT(I+1)
c++ CODE for VERSION = C is inactive
c%%         sprintf(&cmessc.buf[cmessi.lstrt-1], "%*d",
c%%            cmessi.lenout, mact[i]);
c++ END
            go to 960
         end if
c                            Regular integer output
c++ CODE for VERSION = F is active
         write (BUF(LSTRT:LBUF), FMTI) (IDAT(K), K = NIDAT,
     1      NIDAT+KLINE-1)
         NIDAT = NIDAT + KLINE
c++ CODE for VERSION = C is inactive
c%%  kk = cmessi.nidat;
c%%  for (cmessi.nidat=kk; cmessi.nidat<kk+cmessi.kline; cmessi.nidat++)
c%%     sprintf(&cmessc.buf[cmessi.lstrt+cmessi.lenout*(cmessi.nidat
c%%       - kk) - 1], "%*d", cmessi.lenout, idat[cmessi.nidat-1]);
c++ END
         go to 960
c                           Various floating point output
  943    continue
c++ CODE for VERSION = F is active
         FMTF = '(   99F  .  )'
         go to 946
c++ END
  944    continue
c++ CODE for VERSION = F is active
         FMTF = '(1P,99G  .  )'
         go to 946
c++ END
  945    continue
c++ CODE for VERSION = F is active
         FMTF = '(1P,99E  .  )'
c++ END
  946    JJ = mod(JJ/10, 100)
c++ CODE for VERSION = F is active
         FMTF(8:8) = char(ICHAR0 + LENOUT / 10)
         FMTF(9:9) = char(ICHAR0 + mod(LENOUT, 10))
         FMTF(11:11) = char(ICHAR0 + JJ / 10)
         FMTF(12:12) = char(ICHAR0 + mod(JJ, 10))
c++ CODE for VERSION = C is inactive
c%%      strcpy(cmessc.fmtf, "%*.*E\0");
c        IWF = LENOUT
c        lfprec = JJ
c++ END
         if (.not. XARGOK) go to 180
         MPT = NFDAT
         NFDAT = NFDAT + KLINE
         return
c                           Text output
  948    K1 = NTEXT + LBUF - LSTRT
c++ CODE for VERSION = F is active
         BUF(LSTRT:LBUF) = TEXT(ITEXT)(NTEXT:K1-1)
c++ CODE for VERSION = C is inactive
c%%    memcpy(&cmessc.buf[cmessi.lstrt-1], TEXT(cmessi.itext-1,
c%%       cmessi.ntext -1), k1 - cmessi.ntext);
c++ END
         NTEXT = K1
      else
c                                 Print the heading
         KT = 1
         call MESSMH(TEXT)
         if (KT .lt. 0) go to 180
      end if
  960 if ((LBUF .le. MDAT(NROW)) .and. (NCOL .gt. 0)) go to 940
      if (NROW .eq. 1) then
         JJ = LBUF
         LBUF = MDAT(1)
         call MESSPR
         LBUF = JJ
      else
         if (IROW .eq. 0) then
            if (NROW .eq. 2) then
c++ CODE for VERSION = F is active
               if (SUNIT .le. 0) SUNIT = MESSGS()
               rewind(SUNIT)
c++ CODE for VERSION = C is inactive
c%%        if (cmessi.sunit == -1) {
c%%           scratch_file = tmpfile();
c%%           cmessi.sunit = 1;}
c%%        rewind(scratch_file);
c++ END
            end if
         end if
c++ CODE for VERSION = F is active
         write(SUNIT) BUF(5:MDAT(NROW))
c++ CODE for VERSION = C is inactive
c%%       fwrite(&cmessc.buf[4], cmessi.mdat[cmessi.nrow]-4, 1,
c%%          scratch_file);
c++ END
      end if
      if (LBUF .gt. MDAT(NROW)) then
c++ CODE for VERSION = F is active
         BUF(5:LBUF - MDAT(NROW) + 4) = BUF(MDAT(NROW)+1:LBUF)
c++ CODE for VERSION = C is inactive
c%%   memcpy(&cmessc.buf[4], &cmessc.buf[cmessi.mdat[cmessi.nrow-1]],
c%%     cmessi.lbuf - cmessi.mdat[cmessi.nrow-1]);
c++ END
         LBUF = LBUF - MDAT(NROW) + 4
         NROW = NROW + 1
         if (.not. GOTFMT) then
            if (NROW .gt. 5) go to 180
            MDAT(NROW) = LBUF
         end if
         if (NCOL .eq. 0) go to 960
         go to 940
      end if
      LBUF = 0
      if (.not. GOTFMT) then
         GOTFMT = .true.
         IROW = IROW - 1
         go to 920
      end if
      MACT(I+1) = MACT(I+1) + 1
      if (MACT(I+1) .le. MACT(I+2)) go to 999
      MACT(I+1) = 1
      if (NROW .eq. 1) go to 999
c++ CODE for VERSION = F is active
      endfile SUNIT
c++ CODE for VERSION = C is inactive
c%%    fputc(EOF, scratch_file);
c++ END
      KK = 1
  994 KK = KK + 1
      if (KK .gt. NROW) go to 999
c++ CODE for VERSION = F is active
      rewind(SUNIT)
c++ CODE for VERSION = C is inactive
c%%   rewind(scratch_file);
c++ END
      IROW = -1
      K = KK
  995 LBUF = 5
      do 996 J = 2, K
         if (J .eq. K) LBUF = MDAT(KK)
c++ CODE for VERSION = F is active
         read(SUNIT, END = 994) BUF(5:LBUF)
c++ CODE for VERSION = C is inactive
c%%       if (fread(&cmessc.buf[4], cmessi.lbuf-4, 1,
c%%         scratch_file) == 0) goto L_994;
c++ END
  996 continue
      K = NROW
      IROW = IROW + 1
      if (IROW .ne. 0) then
c++ CODE for VERSION = F is active
         write(BUF(1:4), '(I4)') mod(IROW, 10000)
c++ CODE for VERSION = C is inactive
c%%      sprintf(&cmessc.buf[cmessi.lstrt-1], "%4d",  irow%10000);
c++ END
      else
c++ CODE for VERSION = F is active
         BUF(1:4) = ' '
c++ CODE for VERSION = C is inactive
c%%    for (kc=0; kc < 4; kc++) cmessc.buf[kc] = ' ';
c++ END
      end if
      call MESSPR
      go to 995
  999 LENTRY = 1
      return
c
c                          Get started with vector or matrix output
 1000 LOCBEG = NIDAT
      XARG = M .gt. MEJMAT
      if (XARG) then
         M = M - 4
         LOCBEG = NFDAT
         if (.not. XARGOK) go to 40
      end if
      GOTFMT = M .gt. MEIMAT
      if (GOTFMT) M = M - 2
      MPT = LOCBEG
      if (M .eq. MEIMAT) go to 1300
c                           Take care of setup for vector output
      KNT = 0
      LASKNT = MACT(I+1)
      if (LASKNT .le. 0) then
         LASKNT = -LASKNT
         KNT = LOCBEG - 1
         if (LASKNT .le. KNT) go to 40
      end if
      IMAG = LASKNT
      LASTI = LOCBEG + LASKNT - 1 - KNT
      NCOL = 0
c                          Get format for label output.
      call MESSFI
c++ CODE for VERSION = F is active
      FMTR = FMTI
c++ CODE for VERSION = C is inactive
c%%   messcc.kcrwid = messcc.kciwid;
c++ END
      KDILAB = KDI+1
      LINSTR = 2*KDILAB+1
      if (XARG) then
         if (.not. GOTFMT) go to 1150
         IWF = IWG
         FMTF = FMTG
c++ CODE for VERSION = C is inactive
c%%      cmessi.iwf = cmessi.iwg;
c%%      messcc.lfprec = messcc.lgprec;
c++ END
         go to 1160
      end if
      call MESSFD(IDAT)
c                                          After integer format
      LENOUT = KDI
      NIDAT = LASTI + 1
c                                          Common code continues here
 1080 NLINE = (LENLIN - LINSTR + 1) / LENOUT
      if (LBUF .eq. 0) go to 1090
      K = max(LINSTR, LBUF+1)
      if (((K-LINSTR)/LENOUT + (LENLIN-K+1)/LENOUT) .lt. NLINE) K = K +
     1   LENOUT - mod(K-LINSTR, LENOUT)
      KLINE = (LENLIN - K + 1) / LENOUT
      if (KLINE .lt. min(LASKNT-KNT, NLINE/2)) go to 1085
      LINSTR = K - LENOUT * ((K - LINSTR) / LENOUT)
      if (KLINE .ge. LASKNT-KNT)  then
         KLINE = LASKNT - KNT
         K = LBUF + 1
      end if
      KNT = KNT + KLINE
c++ CODE for VERSION = F is active
      BUF(LBUF+1:K) = ' '
c++ CODE for VERSION = C is inactive
c%%    for (kc=cmessi.lbuf; kc < k; kc++) cmessc.buf[kc] = ' ';
c++ END
      LBUF = K
      go to 1110
 1085 call MESSPR
 1090 BUF = ' '
c++ CODE for VERSION = F is active
      write (BUF(1:KDILAB), FMTR) KNT+1
c++ CODE for VERSION = C is inactive
c%%  sprintf(cmessc.buf, "%*d", messcc.kcrwid, knt+1);
c++ END
      BUF(KDILAB:KDILAB) = '-'
      KLINE = min(NLINE, LASKNT - KNT)
      KNT = KNT + KLINE
c++ CODE for VERSION = F is active
      write (BUF(KDILAB+1:2*KDILAB), FMTR) KNT
c++ CODE for VERSION = C is inactive
c%%  sprintf(&cmessc.buf[kdilab], "%*d", messcc.kcrwid, knt);
c++ END
c++ CODE for VERSION = F is active
      BUF(2*KDILAB:LINSTR-1) = ':'
c++ CODE for VERSION = C is inactive
c%%    cmessc.buf[kdilab*2L-1] = ':';
c%%    for (kc=kdilab*2L; kc < *linstr-1; kc++) cmessc.buf[kc] = ' ';
c++ END
      LBUF = LINSTR
 1110 LSTRT = LBUF
      LBUF = LBUF + LENOUT * KLINE - 1
      if (XARG) return
c                                    Integer output
c++ CODE for VERSION = F is active
      write (BUF(LSTRT:LBUF), FMTI) (IDAT(K), K = MPT, MPT+KLINE-1)
c++ CODE for VERSION = C is inactive
c%%   for (k=cmessi.mpt; k<=cmessi.mpt+cmessi.kline-1; k++)
c%%  sprintf(&cmessc.buf[cmessi.lstrt+messcc.kciwid*(k-cmessi.mpt)-1],
c%%      "%*d", messcc.kciwid, idat[k-1]);
c++ END
      MPT = MPT + KLINE
c
c                                     Entry here after vector output.
 1130 if (MPT .le. LASTI) go to 1085
      go to 40
c                                          Get other format
 1150 LENTRY = 5
      return
c                                          After other format
 1160 LENOUT = IWF
      LENTRY = 7
      NFDAT = LASTI + 1
      go to 1080
c
c                           Take care of setup for matrix output
 1300 NDIM = MACT(I+1)
      if (NDIM .le. 0) go to 40
      ICOL = 1
      IROW1 = 1
      NROW = MACT(I+2)
      if (NROW .le. 0) then
         if (NROW .eq. 0) go to 40
         IROW1 = -NROW / 100000
         NROW = -NROW - 99999 * IROW1 - 1
      end if
      NCOL = MACT(I+3)
      if (NCOL .le. 0) then
         if (NCOL .eq. 0) go to 40
         ICOL = -NCOL / 100000
         NCOL = -NCOL - 99999 * IROW1 - 1
      end if
      NTXTSV = NTEXT
      ITXTSV = ITEXT
      IRC = 1
c                        Compute widths for row and column labels
 1320 MAXWID(IRC) = 0
      MTEXT(IRC) = MACT(I+IRC+3)
      IMAG = NROCO(IRC)
      KLINE = IMAG
 1330 NTEXT = MTEXT(IRC)
      if (NTEXT .ge. 0) then
         if (NTEXT .eq. 0) then
            LTEXT = 5
         else
c                        Go get row/column widths
            KT = 2
            call MESSMH(TEXT)
            if (KT .lt. 0) then
               MTEXT(IRC) = 0
               go to 1330
            end if
         end if
         call MESSFI
         MAXWID(IRC) = max(MAXWID(IRC), LTEXT + KDI)
c++ CODE for VERSION = F is active
         FMTIM(IRC) = FMTI
c++ CODE for VERSION = C is inactive
c%%      if (cmessi.irc == 1)
c%%         messcc.kcrwid = cmessi.kdi;
c%%      else
c%%         messcc.kccwid = cmessi.kdi;
c++ END
      end if
      IRC = IRC + 1
      if (IRC .eq. 2) go to 1320
c                 Widths for Row and column titles have been computed.
      KSHIFT = 1
      LASTI = LOCBEG + NROW - IROW1
      if (XARG) then
         if (.not. GOTFMT) go to 1610
c++ CODE for VERSION = F is active
         IWF = IWG
         FMTF = FMTG
c++ CODE for VERSION = C is inactive
c%%      cmessi.iwf = cmessi.iwg;
c%%      messcc.lfprec = messcc.lgprec;
c++ END
         go to 1620
      end if
      call MESSFD(IDAT)
c
      If (KDI .ge. KOLWID) then
         LENOUT = KDI
      else
         KSHIFT = (KOLWID - KDI + 2) /2
         LENOUT = KOLWID
c++ CODE for VERSION = F is active
         FMTI(5:5) = char(ICHAR0 + KOLWID / 10)
         FMTI(6:6) = char(ICHAR0 + mod(KOLWID, 10))
c++ CODE for VERSION = C is inactive
c%%  messcc.kciwid = *kolwid;
c++ END
      end if
      NIDAT = NIDAT + NDIM*NCOL
c                              Continue with commmon code
 1390 NLINE = (LENLIN - LINSTR) / LENOUT
      if (LBUF .le. LINSTR) go to 1420
 1400 call MESSPR
 1420 IROW = IROW1
      KLINE = min(NLINE, NCOL-ICOL+1)
c                       Output column labels (if any)
      if (MTEXTC .lt. 0) go to 1480
      NTEXT = MTEXTC
      IMAG = ICOL
      KT = 3
      call MESSMH(TEXT)
      if (KT .lt. 0) go to 180
c                       Return from output of column labels.
 1470 MTEXTC = NTEXT
 1480 ICOL = ICOL + KLINE
 1490 call MESSPR
c
c                      Output row labels (if any)
      if (MTEXTR .lt. 0) go to 1520
      if (MTEXTR .eq. 0) then
c++ CODE for VERSION = F is active
         BUF(LBUF+1:LBUF+4) = 'Row '
c++ CODE for VERSION = C is inactive
c%%       memcpy(&cmessc.buf[cmessi.lbuf],"Row ", 4);
c++ END
         LBUF = LBUF + 4
         go to 1515
      end if
      NTEXT = MTEXTR
      ITEXT = (NTEXT-1) / LENTXT
      NTEXT = NTEXT - ITEXT * LENTXT
      ITEXT = ITEXT + 1

c                     Go get text for row label
      NTEXTR = 0
      go to 410
c                     Return from getting text for row label
 1510 if (C .ne. '#') then
         MTEXTR = NTEXT + LENTXT * (ITEXT-1)
c++ CODE for VERSION = F is active
         BUF(LBUF+1:LINSTR) = ' '
c++ CODE for VERSION = C is inactive
c%%    for (kc=cmessi.lbuf; kc < *linstr; kc++) cmessc.buf[kc] = ' ';
c++ END
         go to 1520
      end if
 1515 continue
c++ CODE for VERSION = F is active
      write (BUF(LBUF+1:LINSTR), FMTR) IROW
c++ CODE for VERSION = C is inactive
c%%   sprintf(&cmessc.buf[cmessi.lbuf], "%*d", messcc.kcrwid, irow);
c%%    for (kc=cmessi.lbuf+messcc.kcrwid;
c%%       kc < *linstr; kc++) cmessc.buf[kc] = ' ';
c++ END
 1520 LSTRT = LINSTR + 1
      LBUF = LINSTR + LENOUT*KLINE
      LASTI = MPT + NDIM * KLINE - 1
      if (XARG) return
c                                    Integer output
c++ CODE for VERSION = F is active
         write (BUF(LSTRT:LBUF), FMTI) (IDAT(K), K=MPT,LASTI,NDIM)
c++ CODE for VERSION = C is inactive
c%% for (k=cmessi.mpt; k<=cmessi.lasti; k+=cmessi.ndim)
c%%  sprintf(&cmessc.buf[cmessi.lstrt + messcc.kciwid*(k-cmessi.mpt)/
c%%     cmessi.ndim - 1], "%*d", messcc.kciwid, idat[k-1]);
c++ END
c
c                                     Entry here after matrix output.
 1530 MPT = MPT + 1
      IROW = IROW + 1
c
      if (IROW .le. NROW) go to 1490
      if (ICOL .gt. NCOL) then
         NTEXT = NTXTSV
         ITEXT = ITXTSV
         go to 40
      end if
      MPT = NDIM*(ICOL-1) + 1
      MTEXTR = MACT(I+4)
      call MESSPR
      LBUF = 1
      BUF(1:1) = ' '
      go to 1400
c                                Need to get format for matrix print.
 1610 LENTRY = 6
      return
c                                Entry after got format for matrix print
 1620 If (IWF .ge. KOLWID) then
         LENOUT = IWF
      else
         KSHIFT = (KOLWID - IWF + 2) /2
         LENOUT = KOLWID
c++ CODE for VERSION = F is active
         write (FMTF(7:8), '(I2)') KOLWID
c++ CODE for VERSION = C is inactive
c%%      cmessi.iwf = *kolwid;
c%%      strcpy(cmessc.fmtf, "%*.*E\0");
c++ END
      end if
      NFDAT = NFDAT + NDIM*NCOL
      LENTRY = 8
      go to 1390
      end

      subroutine MESSFD(IDAT)
c Get the format for data to be printed in vectors and arrays.
c
c ************** Variable only used here *******************************
c
c K      Temporary index.
c J      Temporary index.
c IDAT   Input array to MESS
c IMAX   Used when computing largest integer in array.
c IMIN   Used when computing smallest integer in array.
c
      integer J, K, IDAT(*), IMAX, IMIN
c
c For comments on other variables, see the listing for MESS.
      integer   LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF = 250)
      parameter (MEVBAS = 10)
      parameter (MEVLAS = 32)
      logical          XARG, GOTFMT, XARGOK
      integer          EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS), IMAG,
     1   ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ, KLINE, KSCRN,
     2   KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN, LENOUT,
     3   LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT, LSTOP,
     4   LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL, NDIM,
     5   NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
c
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*15, FMTG*15,
     1  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,
     1   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT,
     2   NFDAT, NMDAT, MDAT, TABSPA, ICHAR0, IMAG, IRC, ITEXT, IWF,
     3   IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI, LBUF, LENLIN,
     4   LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT, MAXWID, MPT,
     5   NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
c
C     COMMENTED OUT BY KUMAZAWA
C     save /CMESSI/, /CMESSC/
c
      if (GOTFMT) then
         KDI = KDJ
c++ CODE for VERSION = F is active
         FMTI = FMTJ
c++ CODE for VERSION = C is inactive
c%%      messcc.kciwid = cmessi.kdj;
c++ END
         return
      end if
      K = 1
      IMAX = IDAT(LOCBEG)
      IMIN = IMAX
   10 do 20 J = LOCBEG, LASTI
         IMAX = max(IMAX, IDAT(J))
         IMIN = MIN(IMIN, IDAT(J))
   20 continue
      if (NCOL .ne. 0) then
         K = K + 1
         LOCBEG = LOCBEG + NDIM
         LASTI = LASTI + NDIM
         if (K .le. NCOL) go to 10
      end if
      IMAG = IMAX
      if (IMAG + 10*IMIN .le. 0) IMAG = IMIN
      IMAG = 10 * IMAG
      call MESSFI
      return
      end


      subroutine MESSFI
c Get the format for the integer IMAG.
c
c ************** Variable only used here *******************************
c
c K, KD are used in determining number of characters needed to represent
c       IMAG.
c
      integer K, KD
c
c For comments on other variables, see the listing for MESS.
      integer   LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF = 250)
      parameter (MEVBAS = 10)
      parameter (MEVLAS = 32)
      logical          XARG, GOTFMT, XARGOK
      integer          EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS), IMAG,
     1   ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ, KLINE, KSCRN,
     2   KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN, LENOUT,
     3   LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT, LSTOP,
     4   LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL, NDIM,
     5   NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
c
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*15, FMTG*15,
     1  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,
     1   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT,
     2   NFDAT, NMDAT, MDAT, TABSPA, ICHAR0, IMAG, IRC, ITEXT, IWF,
     3   IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI, LBUF, LENLIN,
     4   LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT, MAXWID, MPT,
     5   NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
c

C     COMMENTED OUT BY KUMAZAWA
C     save /CMESSI/, /CMESSC/
c
      KD = 0
      K = 1
      if (IMAG .lt. 0) then
         IMAG = -IMAG
         KD = KD + 1
      end if
   10 K = 10 * K
      KD = KD + 1
      if (IMAG .ge. K) go to 10
      if (KD .ne. KDI) then
         KDI = KD
c++ CODE for VERSION = F is active
         FMTI(5:5) = char(ICHAR0 + KDI / 10)
         FMTI(6:6) = char(ICHAR0 + mod(KDI, 10))
c++ CODE for VERSION = C is inactive
c%%      messcc.kciwid = cmessi.kdi;
c++ END
      end if
      return
      end

c++ CODE for VERSION = F is active
      integer function MESSGS()
c                                 Get a scratch unit assigned.
      integer J
c
      MESSGS = 31
   10 MESSGS = MESSGS - 1
      if (MESSGS .eq. 0) stop 'Could not assign scratch unit in MESS.'
      open (MESSGS, STATUS='SCRATCH', ACCESS='SEQUENTIAL',
     1    FORM='UNFORMATTED', IOSTAT=J)
      if (J .ne. 0) go to 10
      return
      end
c++ END

      subroutine MESSMH(TEXT)
c Processing of multiple headings:
c
c J     Used as a temporary index.
c K     Used as a temporary index.
c KK    Used as a temporary index.
c KT    Used for logic in output of headings.  Set <0 on exit if there
c       is an input error.
c     KT = 1 Output table headings.  (Set to -1 on fatal error.)
c     KT = 2 Get row/column widths for matrix output.  (Set to -2 if
c            error results in no headings.)
c     KT = 3 Output column headings.  (Set to -1 on fatal error.)
c L     Used as a temporary index.
c TEXT  Original input character vector.
c
      integer J, K, KK, L
      character*(*)  TEXT(*)
      character SC, C
      parameter (SC = '$')
c For comments on other variables, see the listing for MESS.
      integer   KOLWID, LINSTR, LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF = 250)
      parameter (MEVBAS = 10)
      parameter (MEVLAS = 32)
      logical          XARG, GOTFMT, XARGOK
      integer          EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS), IMAG,
     1   ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ, KLINE, KSCRN,
     2   KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN, LENOUT,
     3   LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT, LSTOP,
     4   LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL, NDIM,
     5   NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
c
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*15, FMTG*15,
     1  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,
     1   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT,
     2   NFDAT, NMDAT, MDAT, TABSPA, ICHAR0, IMAG, IRC, ITEXT, IWF,
     3   IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI, LBUF, LENLIN,
     4   LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT, MAXWID, MPT,
     5   NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
c
      equivalence (MAXWID(1), LINSTR), (MAXWID(2), KOLWID)

C     COMMENTED OUT BY KUMAZAWA
C     save /CMESSI/, /CMESSC/
c++ CODE for VERSION = C is inactive
c      integer kc
c++ END
c
      if (NTEXT .ne. 0) then
         ITEXT = (NTEXT-1) / LENTXT
         NTEXT = NTEXT - ITEXT * LENTXT
         ITEXT = ITEXT + 1
      end if
      do 300 J = 1, max(1,KLINE)
         if (NTEXT .eq. 0) then
            K = KOLWID
            go to 210
         end if
         LTEXT = 0
  110    continue
c++ CODE for VERSION = F is active
         L = index(TEXT(ITEXT)(NTEXT:LENTXT), SC)
c++ CODE for VERSION = C is inactive
c%%       ctmp=memchr(TEXT(cmessi.itext-1L,cmessi.ntext-1), SC,
c%%          cmessi.lentxt - cmessi.ntext + 1);
c%%       if (ctmp == NULL)
c%%             l = 0;
c%%       else
c%%             l = ctmp - TEXT(cmessi.itext-1L,cmessi.ntext-1) + 1;
c++ END
         if (L .eq. 0) then
            LTEXT = LTEXT + LENTXT - NTEXT + 1
            if (LTEXT .lt. 80) then
               ITEXT = ITEXT + 1
               NTEXT = 1
               go to 110
            end if
            LTEXT = 0
            if (KT .eq. 3) go to 310
            go to 160
         end if
         NTEXT = NTEXT + L + 1
         LTEXT = L + LTEXT - 1
         if (NTEXT .gt. LENTXT) then
            ITEXT = ITEXT + 1
            if (NTEXT .eq. LENTXT + 1) then
               C = TEXT(ITEXT-1)(LENTXT:LENTXT)
               NTEXT = 1
            else
               C = TEXT(ITEXT)(1:1)
               NTEXT = 2
            end if
         else
            C = TEXT(ITEXT)(NTEXT-1:NTEXT-1)
         end if
         if (C .eq. 'H') go to (180, 190, 200), KT
         if (C .eq. 'E') go to (180, 310, 200), KT
         if (C .eq. '#') go to (140, 310, 200), KT
         if (KT .ne. 1) go to 160
  140    LTEXT = LTEXT + 2
         go to 110
  160    KT = -KT
         go to 310
c
  180    KOLWID = KOLWID + LENOUT
         if (LTEXT .eq. 0) go to 300
         KK = KOLWID-LTEXT
         if (KK .lt. 0) stop
     1   'Stopped in MESS -- Column width too small in a heading.'
         if (XARG)  KK = 1 + KK/2
         LSTRT = LBUF + KK + 1
         LBUF = LBUF + KOLWID
         if (LBUF .le. LENLIN) MDAT(NROW) = LBUF
         KOLWID = 0
         go to 220
c
c                                  Set up column widths
  190    MAXWID(IRC) = max(MAXWID(IRC), LTEXT)
         go to 300
c
c                                  Output matrix column
  200    K = KOLWID
         if (C .ne. '#') K = LTEXT
  210    KK = LENOUT - KOLWID
         if (J .eq. 1) then
c                        Special setup for the first column.
            if (XARG) KK = (KK + 1) / 2
            KK = KK + KSHIFT + LINSTR - LBUF
         end if
         KK = KK + KOLWID - K
         LSTRT = LBUF + KK + 1
         LBUF = LSTRT + K - 1
c                                  Set initial blanks
  220    continue
c++ CODE for VERSION = F is active
         if (KK .gt. 0) BUF(LSTRT-KK:LSTRT-1) = ' '
c++ CODE for VERSION = C is inactive
c%%      if (kk > 0) for (kc=cmessi.lstrt-kk-1; kc<cmessi.lstrt-1; kc++)
c%%         cmessc.buf[kc] = ' ';
c++ END
c                                  Move characters
         if (NTEXT .eq. 0) then
c++ CODE for VERSION = F is active
            BUF(LSTRT:LSTRT+3) = 'Col '
c++ CODE for VERSION = C is inactive
c%%       memcpy(&cmessc.buf[cmessi.lstrt-1],"Col ", 4);
c++ END
            C = '#'
            LSTRT = LSTRT+4
         else
            K = NTEXT - LTEXT - 2
            if (K .le. 0) then
               KK = max(0, 3-NTEXT)
c++ CODE for VERSION = F is active
               BUF(LSTRT:LSTRT-K-KK)=TEXT(ITEXT-1)(LENTXT+K:LENTXT-KK)
c++ CODE for VERSION = C is inactive
c%%       memcpy(&cmessc.buf[cmessi.lstrt-1L], TEXT(cmessi.itext-2L,
c%%         cmessi.lentxt+k-1), -k-kk+1L);
c++ END
               LSTRT = LSTRT-K-KK+1
               K = 1
            end if
            if (NTEXT .gt. 3) then
c++ CODE for VERSION = F is active
               BUF(LSTRT:LSTRT+NTEXT-K-3) = TEXT(ITEXT)(K:NTEXT-3)
c++ CODE for VERSION = C is inactive
c%%       memcpy(&cmessc.buf[cmessi.lstrt-1L], TEXT(cmessi.itext-1L,
c%%         k-1), cmessi.ntext-k-2L);
c++ END
               LSTRT = LSTRT + NTEXT - K - 2
            end if
         end if
         if (C .eq. '#') then
c                                  Output column index
c++ CODE for VERSION = F is active
            write (BUF(LSTRT:LBUF), FMTC) IMAG + J - 1
c++ CODE for VERSION = C is inactive
c%%         sprintf(&cmessc.buf[cmessi.lstrt-1], "%*d ",
c%%           cmessi.lbuf-cmessi.lstrt, cmessi.imag+j-1);
c++ END
            if (NTEXT .ne. 0) NTEXT = K
            go to 300
         end if
c                                  Set trailing blanks
c++ CODE for VERSION = F is active
         if (LSTRT .le. LBUF) BUF(LSTRT:LBUF) = ' '
c++ CODE for VERSION = C is inactive
c%%      if (cmessi.lstrt <= cmessi.lbuf)
c%%           for (kc=cmessi.lstrt-1; kc < cmessi.lbuf; kc++)
c%%              cmessc.buf[kc] = ' ';
c++ END
  300 continue
  310 return
      end

      subroutine MESSPR
c Prints the buffer for MESS
c
c ************** Variable only used here *******************************
c
c NSCRN  Number of lines currently on CRT from messages.
c
      integer   NSCRN
      save      NSCRN
c
c For comments on other variables, see the listing for MESS.
      integer   LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF = 250)
      parameter (MEVBAS = 10)
      parameter (MEVLAS = 32)
      logical          XARG, GOTFMT, XARGOK
      integer          EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS), IMAG,
     1   ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ, KLINE, KSCRN,
     2   KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN, LENOUT,
     3   LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT, LSTOP,
     4   LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL, NDIM,
     5   NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
c
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*15, FMTG*15,
     1  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,
     1   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT,
     2   NFDAT, NMDAT, MDAT, TABSPA, ICHAR0, IMAG, IRC, ITEXT, IWF,
     3   IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI, LBUF, LENLIN,
     4   LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT, MAXWID, MPT,
     5   NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
c
C     COMMENTED OUT BY KUMAZAWA
C     save /CMESSI/, /CMESSC/
      data NSCRN / 0 /
c
      if (LBUF .ne. 0) then
         if (OUNIT .le. 0) then
            if (KSCRN .gt. 0) then
               if (NSCRN .ge. KSCRN) then
c++ CODE for VERSION = F is active
                  print '('' Type "Enter" to continue'')'
c++ CODE for VERSION = C is inactive
c%%               fprintf( stdout, " Type 'Enter' to continue\n" );
c++ END
                  read (*, *)
                  NSCRN = 0
               end if
               NSCRN = NSCRN + 1
            end if
c++ CODE for VERSION = F is active
            print '(1X, A)', BUF(1:LBUF)
c++ CODE for VERSION = C is inactive
c%%      fprintf(stdout, "%.*s\n", cmessi.lbuf, cmessc.buf);
c++ END
            if (OUNIT .eq. 0) go to 10
         end if
c++ CODE for VERSION = F is active
         write (abs(OUNIT), '(A)') BUF(1:LBUF)
c++ CODE for VERSION = C is inactive
c%%      fprintf(c_handle[labs(cmessi.ounit)-1], "%.*s\n", cmessi.lbuf,
c%%         cmessc.buf);
c++ END
   10    LBUF = 0
      end if
      return
      end

      subroutine MESSFT(MACT, FTEXT)
c  Prints FTEXT, which contains a Fortran character string, and then
c  call MESSS to do the actions in MACT.  Actions in MACT can not do
c  anything other than actions that reference MACT.
c  This routine intended for use by library subroutines getting text in
c  the form of a Fortran character string.
c
      integer MACT(*)
      character FTEXT*(*)
c
      integer K, IDAT(1), MECONT
c++ CODE for VERSION = F is active
      character TEXT(1)*1
c++ CODE for VERSION = C is inactive
c      character TEXT(1)*2
c++ END
      parameter (MECONT=50)
c
      parameter (LENBUF = 250)
      parameter (MEVBAS = 10)
      parameter (MEVLAS = 32)
      logical          XARG, GOTFMT, XARGOK
      integer          EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS), IMAG,
     1   ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ, KLINE, KSCRN,
     2   KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN, LENOUT,
     3   LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT, LSTOP,
     4   LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL, NDIM,
     5   NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
c
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*15, FMTG*15,
     1  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,
     1   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT,
     2   NFDAT, NMDAT, MDAT, TABSPA, ICHAR0, IMAG, IRC, ITEXT, IWF,
     3   IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI, LBUF, LENLIN,
     4   LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT, MAXWID, MPT,
     5   NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
c
      K = len(FTEXT)
      NTEXT = 1
      if (K .ne. 0) then
         if (FTEXT(1:1) .eq. '0') then
            NTEXT = 2
            K = K - 1
            if (LBUF .eq. 0) then
               BUF(1:1) = ' '
               LBUF = 1
            end if
         end if
         call MESSPR
         LBUF = K
c++ CODE for VERSION = F is active
         BUF(1:K) = FTEXT(NTEXT:NTEXT+K-1)
c++ CODE for VERSION = C is inactive
c%%      memcpy(cmessc.buf, &ftext[cmessi.ntext-1], k);
c++ END
      end if
      ICHAR0 = ICHAR('0')
      LENTRY = 2
      if (MACT(1) .ne. MECONT) call MESS(MACT, TEXT, IDAT)
      return
      end

c++ CODE for VERSION = F is active
      block data XMESS
c For comments on these variables, see the listing for MESS.
c
      integer   LNERR, LNMSG, LENBUF, MEVBAS, MEVLAS
      parameter (LNERR = 79)
      parameter (LNMSG = 128)
      parameter (LENBUF = 250)
      parameter (MEVBAS = 10)
      parameter (MEVLAS = 32)
      logical          XARG, GOTFMT, XARGOK
      integer          EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS), IMAG,
     1   ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ, KLINE, KSCRN,
     2   KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN, LENOUT,
     3   LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT, LSTOP,
     4   LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL, NDIM,
     5   NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
c
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*15, FMTG*15,
     1  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,
     1   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT,
     2   NFDAT, NMDAT, MDAT, TABSPA, ICHAR0, IMAG, IRC, ITEXT, IWF,
     3   IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI, LBUF, LENLIN,
     4   LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT, MAXWID, MPT,
     5   NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
c
      character   DOLE(72)
      equivalence (DOLS,DOLE)
C     COMMENTED OUT BY KUMAZAWA
C     save /CMESSI/, /CMESSC/
c
      data EUNIT,KDFDEF,MAXERR,LSTOP,LHEAD,LPRINT,LINMSG,LINERR,MUNIT
     1   /     0,     0,     0,    3,    1,     3, LNMSG, LNERR,    0 /
      data KDF, LENLIN, LENTRY, OUNIT, SUNIT, TABSPA, KDI, KDJ
     1  /    0,  LNMSG,      0,     0,    -1,      6,   0,   6 /
      data (DOLE(I), I = 1, 72) / 72*'$' /
      data FMTI, FMTJ / '(99Ixx)', '(99I06)' /
      data FMTG / '(1P,99Exx.xx)  ' /
      end
c++ END

      subroutine SILUP (X, Y, NTAB, XT, YT, NDEG, LUP, IOPT, EOPT)
c     .  Copyright (C) 1989, California Institute of Technology.
c     .  All rights reserved.  U. S. Government sponsorship under
c     .  NASA contract NAS7-918 is acknowledged.
c>> 1993-04-28 SILUP  Krogh  Additions for Conversion to C.
c>> 1992-05-27 SILUP  Krogh  Fixed bug in error estimate.
c>> 1992-04-08 SILUP  Krogh  Removed unused labels 510 and 2080.
c>> 1991-10-17 SILUP  Krogh  Initial Code.
c
c Polynomial Interpolation with look up.
c Design/Code by  Fred T. Krogh, Jet Propulsion Laboratory, Pasadena, CA
c Last revision: May 7, 1991.
c
c In addition to doing standard 1 dimensional interpolation, this
c subroutine supports efficient interpolation of several functions
c defined at the same values of the independent variable and supports
c SILUPM, a subroutine for doing multidimensional interpolation.  Error
c estimates and derivatives of the interpolant can be obtained via
c options.
c Algorithms used are described in "Efficient Algorithms for Polynomial
c Interpolation and Numerical Differentiation", by Fred T. Krogh, Math.
c of Comp. Vol. 24, #109 (Jan. 1970), pp. 185-190.
c
c     *************     Formal Arguments     ***************************
c
c X      Independent variable where value of interpolant is desired.
c Y      Value of interpolant, computed by this subroutine.  The
c        interpolant is always a piecewise polynomial.
c NTAB   Number of points in the table.
c XT     Array of independent variable values.  Must be monotone
c        increasing or monotone decreasing.  If XT(I) = XT(I+1) for some
c        I, then the XT(J)'s that are used in the interpolation will
c        either have all J's .le. I, or all J's .ge. I+1.  If the XT's
c        are equally spaced an option allows one to provide only XT(1),
c        and the increment between the XT's.
c YT     Array of dependent variable values.  Y(XT(I)) = YT(I).
c NDEG   If NDEG < 2 or odd, then it gives the degree of the polynomial
c        used.  Else the polynomial used is a linear combination of two
c        polynomials of degree NDEG, such that the resulting polynomial
c        is of degree NDEG+1 and has a continuous first derivative.
c        Typically accuracy will improve as NDEG increases up to some
c        point where it will start to get worse because of either
c        rounding errors or the inherant instability of high degree
c        polynoimial interpolation. If more than MAXDEG-th degree is
c        desired, parameter MAXDEG must be changed below.  It is
c        currently 15.  If X is so close to the end of the table that
c        the same number of points can not be selected on both sides of
c        x the degree of the interpolant is NDEG exactly.  When
c        extrapolating, the degree used is max(2, 2*(NDEG/2)), where the
c        divide is truncates to the nearest integer.
c LUP    Defines the type of look up method.  (Changed if LUP < 1.)
c    < 0   Use a sequential search starting with an index of -LUP, and
c          set LUP to -k on exit where k minimizes abs(X-XT(k)).
c    = 0   As for < 0, except start with a binary search.
c    = 1   Use a binary search.
c    = 2   Start a sequential search with an index = [1.5 + 
c          (NTAB-1) * (X-XT(1)) / (XT(NTAB)-XT(1))].
c    = 3   YT(k) corresponds to XT(1) + (k-1)*XT(2), no search is needed
c          to do the look up and XT can have dimension 2.
c    = 4   Internal information connected with X-XT(k) values used in
c          in the last interpolation is reused.  (Only use if there are
c          no intervening calls to SILUP.)  No options should be
c          specified and only YT should be different.  Intended for
c          interpolating components after the first of a vector valued
c          function.
c IOPT   IOPT(1) is used to return a status as follows.
c    -9    An option index is out of range.
c    -8    NTAB is outside of allowed limits.
c    -7    NDEG is outside of allowed limits.
c    -6    Option 6 (use old parameters), used when not ready for it.
c    -5    Option 3 (compute derivatives), has requested more than
c          MAXDEG derivatives.
c    -4    LUP = 3 and XT(2) = 0.
c    -3    XT(1) = XT(NTAB), and NTAB is not 1.
c    -2    Only one table entry available, req. err. est. not computed.
c    -1    The accuracy requested was not obtained.
c     0    Nothing special to flag.
c     1    X was outside the domain of the table, extrapolation used.
c     2    NTAB is so small, it restricted the degree of the polynomial.
c
c        Starting with IOPT(2) options are specified by integers in the
c        range 0-6, followed in some cases by integers providing
c        argument(s).  Each option, together with its options its
c        arguments if any is followed in IOPT by the next option (or 0).
c     0    No more options; this must be last in the option list.
c     1    An error estimate is to be returned in EOPT(1).
c     2    (Argument: K2)  K2 gives the polynomial degree to use when
c          extrapolating.
c     3    (Argument: K3, L3) Save (k-th derivative of interpolating
c          polynomial) / (k!) in EOPT(K3+K-1) for k = 1, 2, ..., L3.
c          These values are the coefficients of the polynomial in the
c          monomial basis expanded about X.  One must have 0<L3<NDEG+1.
c     4    (Argument K4) The absolute and relative errors expected in
c          YT entries are specified in EOPT(K4) and EOPT(K4+1)
c          respectively.  The values provided here are used in
c          estimating the error in the interpolation which is stored
c          in EOPT(1).
c     5    (Argument K5, L5) Do the interpolation to the accuracy
c          requested by the absolute error tolerance specified in
c          EOPT(K5) and the relative error tolerance in EOPT(K5+1)
c          respectively.  An attempt is made to keep the final error
c          < EOPT(K5) + EOPT(K5+1) * (abs(YT(i1) +abs(YT(i2)), where
c          i1 and i2 are indices for table values close to X.
c          Interpolation is done as for the default case, except that
c          in this case NDEG gives the maximal degree polynomial to use
c          in the interpolation, and polynomial interpolation of even
c          degree can happen.  (The form that gives the C1 continuity
c          is not available in this case, and in fact continuity of the
c          interpolant itself is not to be expected.  The actual degree
c          used in doing the interpolation is stored in the space for
c          the argument L5.  An error estimate is returned in EOPT(1).
c     6    (Argument K6) Do not use point k in the interpolation if
c          YT(k) = EOPT(K6).
c  >253    Used for calls from SILUPM, a number of variables in the
c          common block have been set.  If this is used, it must be the
c          first option, and common variables must be set.  The rest
c          of IOPT is not examined for options.
c EOPT   Array used to return an error estimate and also used for
c        options.
c     EOPT(1)  if an error estimate is returned, this contains a crude
c              estimate of the error in the interpolation.
c
c      ************     External Procedures      ***********************
c
c R1MACH Returns system parameters.
c SMESS  Prints error messages.
c
c      ************     Variables referenced     ***********************
c
c BADPT  In common CSILUP.  YT(K) is not to be used if YT(K) = BADPT,
c        and LBADPT = .true.
c CY     Array giving the YT values used to construct the interpolant.
c  Also used to store divided differences.
c DX     Array giving the differences X - XT(I), where I runs through
c  the points that are used in the interpolation.
c DXS    Saved value of DX when computing derivatives
c E1     Error estimate from one iteration ago (for variable order) / 8
c E2     Error estimate from this iteration (for variable order)
c E2L    Used in computing error estimate for variable order.
c EBND   If estimated error is < EBND then order is suff. high.
c EBNDI  Internal value for bounding error.
c EBNDR  Used in getting part of EBND due to relative accuracy request.
c EF     Factor used in estimating errors for variable order.
c EN     Used in computing error estimate for variable order.
c EOPT   Formal argument, see above.
c EPSR   Relative error level (R1MACH(4)).
c ERRDAT Holds floating point data for error messages.
c ERREST Estimate of the error in the interpolation.
c GETERR Logical variable that is .true. if we are getting an error est.
c H      In indexed lookup contains the difference between XT values.
c I      Temporary index.
c IDIR   When LBADPT = .true. this is 1 if XT values are increasing and
c        is -1 when they are decreasing.
c IDX    In common CSILUP.  Gives indices of points selected for the
c        interpolation in order selected.
c IIFLG  Internal value for IOPT(1).
c ILI    Lower (sometimes upper) bound on the XT indices that will be
c  used in the interpolation.
c INDXED Logical variable that is .true. if we have an indexed XT.
c IOPT   Formal argument, see above.
c IUI    As for ILI, except an upper (sometimes lower) bound.
c K      Index usually into YT and XT, into CY for derivatives.
c KAOS   In common CSILUP.  Keeps track of state on variable order.
c    = 0   No variable order -- this value only set in SILUPM.
c    = 1   First time
c    = 2   Second time
c    = 3   Just had an increase or a very weak decrease in the error.
c    = 4   Just had a strong decrease in the error.
c KEXTRP In common CSILUP.  Gives degree to use when extrapolating.
c KGO    Indicates type of interpolation as follows:
c    1     Standard polynomial interpolation.
c    2     Get error estimate after 1.
c    3     Compute an interpolant with some desired accuracy.
c    4     Compute an interpolant with a continuous derivative.
c  >10     Same as if KGO were 10 smaller, except XT values have already
c          been selected.
c KGOC   In common CSILUP.  Value saved for -KGO.  If YT contains values
c        of Y computed elsewhere in the order specified in IDX, then
c        KGOC should be set to abs(KGOC) before calling SILUP.  If -2,
c        an extra value was computed for an error estimate.
c KK     Temporary index used in selecting data points to use.
c L      Temporary index used in searching the XT array.
c LBADPT Logical variable set = .true. if checking for bad points.
c LDERIV In common CSILUP. Loc. where first deriv. is stored in EOPT().
c LEXERR In common CSILUP. Loc. where absolute and relative error
c        information on Y is stored in EOPT.
c LEXIT  In common CSILUP.  Defines action after finding location in XT.
c    0     Don't compute anything for Y, just return. (Used for SILUPM.)
c    1     Usual case, just interpolate.
c   >1     Compute LEXIT-1 derivatives of the interpolant.
c LINC   Logical variable used in the sequential search.  Set .true. if
c  the values in XT are increasing with I, and set .false. otherwise.
c LNDEG  Location in IOPT to save degree when variable order is used.
c LOPT   Index of the last option.
c LUP    Formal argument, see above.
c MACT   Array used for error message actions, see SMESS for details.
c MAXDEG Parameter giving the maximum degree polynomial interpolation
c  supported.  MAXDEG must be odd and > 2.
c MESS   Message program called from SMESS to print error messages.
c MEEMES Parameter giving value for printing an error message in MESS.
c MEIVEC Parameter giving value for printing an integer vector in MESS.
c MEMDA1 Parameter giving value for indicating value in MACT for MESS.
c MEMDA2 Parameter giving value for indicating value in MACT for MESS.
c MEMDAT Parameter giving value for specifying MACT index for MESS.
c MENTXT Parameter giving value for specifying location in MTXTAA to
c  start print in MESS.
c MERET  Parameter giving value to indicate no more actions for MESS.
c METEXT Parameter giving value to tell MESS to print from MTXTAA.
c LTXTxx Parameter names of this form were generated by PMESS in making
c        up text for error messages and are used to locate various parts
c        of the message.
c MLOC   Array giving starting loctions for error message text.
c MTXTAA Character array holding error message text for MESS.
c N      Used in logic for deciding bounds.  If N = 0, ILI can not be
c  reduced further.  If N = 3*NTABI, IUI can not be increased further.
c NDEG   Formal argument, see above.
c NDEGE  Degree of polynomial to be used.  Starts out = NDEG.
c NDEGEC In common CSILUP.  Degree actually used, saved in common.
c NDEGI  Degree of polynomial up to which y & differences are computed.
c NDEGQ  Used when KGO > 10.  In this case If L .ge. NDEGQ special
c        action is needed.  (Usually means quitting.)
c NTAB   Formal argument, see above.
c NTABI  Internal value of NTAB, = number of points in XT and YT.
c PI     Array use to store interpolation coefficients
c PID    Temporary storage used when computing derivatives.
c R1MACH Function to get parameters of floating point arithmetic.
c SMESS  Program calling MESS to print error messages.
c TP1    Used for temporary accumulation of values and temp. storage.
c TP2    Used for temporary storage.
c X      Formal argument, see above.
c XI     Internal value for X, the place where interpolation is desired.
c XL     Value of "left hand" X when selecting indexed points.
c XT     Formal argument, see above.
c XU     Value of "right hand" X when selecting indexed points.
c Y      Formal argument, see above.
c YL     Last value of Y when selecting the order.
c YT     Formal argument, see above.
c YTORD  Logical variable set .true. when YT values are ordered on entry
c
c     *************     Formal Variable Declarations     ***************
c
      integer          IOPT(*), LUP, NDEG, NTAB
      real             X, XT(*), Y, YT(*), EOPT(*)
c
c     *************     Common Block and Parameter     *****************
c
c                               MAXDEG must be odd and > 2.
      parameter (MAXDEG = 15)
c
      logical          GETERR
      integer          KAOS, KEXTRP, KGOC,LDERIV,LEXERR,LEXIT,NDEGEC,
     1                 IDX(0:MAXDEG+1)
      real             BADPT, DX(0:MAXDEG+1), EBND, EBNDR
      common / CSILUP / BADPT, DX, EBND, EBNDR, KAOS, KEXTRP, KGOC,
     1                  LDERIV, LEXERR,LEXIT,NDEGEC,IDX,GETERR
c
c     *************     Local Variables     ****************************
c
      logical          INDXED, LBADPT, LINC, YTORD
      integer          I, IDIR, IIFLG, ILI, IUI, K, KGO,
     1                 KK, L, LNDEG, LOPT, N, NDEGE, NDEGI, NDEGQ, NTABI
      real             CY(0:MAXDEG+1), DXS, E1, E2, E2L, EBNDI,
     1                 EF, EM, EPSR, ERREST, H, PI(0:MAXDEG+1),
     2                 PID(0:MAXDEG), TP1, TP2, XI, XL, XU, YL
      real             R1MACH
      save             EPSR
c
c ************************ Error Message Stuff and Data ****************
c
c Parameter defined below are all defined in the error message program
c SMESS.
c
      parameter (MENTXT =23)
      parameter (MEMDAT =26)
      parameter (MEMDA1 =27)
      parameter (MEMDA2 =28)
      parameter (MERET =51)
      parameter (MEEMES =52)
      parameter (METEXT =53)
      parameter (MEIVEC =57)
c
      real             ERRDAT(2)
      integer MLOC(9), MACT(14)
c
c ********* Error message text ***************
c[Last 2 letters of Param. name]  [Text generating message.]
cAA SILUP$B
cAB Estimated error = $F, Requested error = $F.$E
cAC NTAB = 1, no error estimate computed.$E
cAD XT(1) = XT(NTAB=$M) = $F ??$E
cAE LUP = 3, and XT(2) = 0.$E
cAF Order of requested derivative, $M, is > bound of $M.$E
cAG LUP > 3, with unprepared common block.$E
cAH NDEG = $M is not in the allowed interval of [0, $M].$E
cAI NTAB = $M is not in the allowed interval of [1, $M].$E
cAJ IOPT($M) = $M is not a valid option.$E
cAK IOPT(1:$M):$B
      integer LTXTAA,LTXTAB,LTXTAC,LTXTAD,LTXTAE,LTXTAF,LTXTAG,LTXTAH,
     * LTXTAI,LTXTAJ,LTXTAK
      parameter (LTXTAA=  1,LTXTAB=  8,LTXTAC= 53,LTXTAD= 92,LTXTAE=121,
     * LTXTAF=146,LTXTAG=200,LTXTAH=240,LTXTAI=294,LTXTAJ=348,
     * LTXTAK=386)
      character MTXTAA(2) * (199)
      data MTXTAA/'SILUP$BEstimated error = $F, Requested error = $F.$EN
     *TAB = 1, no error estimate computed.$EXT(1) = XT(NTAB=$M) = $F ??$
     *ELUP = 3, and XT(2) = 0.$EOrder of requested derivative, $M, is > 
     *bound of $M.$E','LUP > 3, with unprepared common block.$ENDEG = $M
     * is not in the allowed interval of [0, $M].$ENTAB = $M is not in t
     *he allowed interval of [1, $M].$EIOPT($M) = $M is not a valid opti
     *on.$EIOPT(1:$M):$B'/
c
c                      1 2       3 4       5 6 7 8      9 10     11
      data MACT / MEMDA1,0, MEMDA2,0, MEEMES,0,0,0, MERET,1, METEXT,
     1   MEIVEC,0, MERET /
c            12 13    14
c
      data MLOC / LTXTAB, LTXTAC, LTXTAD, LTXTAE, LTXTAF, LTXTAG,
     1            LTXTAH, LTXTAI, LTXTAJ /
c
      data EPSR / 0.E0 /
c
c
c      ************     Start of executable code     *******************
c
      IIFLG = 0
      LBADPT = .false.
      ILI = -LUP
      if (ILI .le. -4) go to 1400
      LEXERR = 0
      NTABI = NTAB
      if ((NTABI .le. 0) .or. (NTABI .gt. 9999999)) then
         IIFLG = -8
         MACT(2) = NTABI
         MACT(4) = 9999999
         go to 2120
      end if
      NDEGI = NDEG
      NDEGE = NDEGI
      if ((NDEGI .lt. 0) .or. (NDEGI .gt. MAXDEG)) then
         IIFLG = -7
         MACT(2) = NDEGI
         MACT(4) = MAXDEG
         go to 2120
      end if
      XI = X
      KGO = 1
      I = 1
      if (IOPT(2) .ge. 254) then
         LBADPT = IOPT(2) .eq. 255
         IDIR = 0
         if (KAOS .le. 0) go to 50
         LNDEG = 0
         KGO = 3
         KAOS = 1
         NDEGE = IOPT(3)
         NDEGI = min(MAXDEG, NDEGE+1)
         NDEGE = -NDEGI
         go to 60
      end if
      GETERR = .false.
      KEXTRP = -1
      LEXIT = 1
c                                      Loop to take care of options
   10 I = I + 1
         LOPT = IOPT(I)
         if (LOPT .ne. 0) then
            go to (1930, 1500, 1600, 1900, 1950, 1970), LOPT
            IIFLG = -9
            MACT(2) = I
            MACT(4) = LOPT
            MACT(13) = min(I, 100)
            MACT(9) = MEMDAT
            go to 2120
         end if
   50 if (KGO .eq. 1) then
         if (NDEGI .ge. 2) then
            if (mod(NDEGI, 2) .eq. 0) then
               NDEGE = NDEGE + 1
               KGO = 4
            end if
         end if
      end if
   60 if (KEXTRP .lt. 0) then
         KEXTRP = max(NDEGE-1, min(NDEGE, 2))
         if (KEXTRP .lt. 0) KEXTRP = NDEGI
      end if
      if (ILI .gt. 0) go to 200
      if (ILI + 2) 1300, 1200, 100
c
c                                      Binary search, then sequential
  100    continue
c                         In binary search XT(ILI) .le. XI .le.  XT(IUI)
c                         or we are extrapolating
            ILI = 1
            IUI = NTABI
            if (XT(NTABI) - XT(1)) 110, 1220, 130
  110       ILI = NTABI
            L = 1
  120       IUI = L
  130       L = (IUI - ILI) / 2
            if (L .eq. 0) go to 210
            L = ILI + L
            if (XT(L) .gt. XI) go to 120
            ILI = L
            go to 130
c
c                                      Sequential search, then exit
  200    continue
            if ((ILI .le. 0) .or. (ILI .gt. NTABI)) go to 100
            if (XT(NTABI) .eq. XT(1)) go to 1220
  210       INDXED = .false.
            LINC = XT(NTABI) .gt. XT(1)
            if ((XT(ILI) .gt. XI) .eqv. LINC) go to 230
  220       if (ILI .eq. NTABI) go to 1240
            ILI = ILI + 1
            if ((XT(ILI) .lt. XI) .eqv. LINC) go to 220
            N = 2*ILI
            if (abs(XT(ILI-1)-XI) .lt. abs(XT(ILI) - XI)) then
               ILI = ILI - 1
               N = N - 1
            end if
            go to 240
  230       if (ILI .eq. 1) go to 1240
            ILI = ILI - 1
            if ((XT(ILI) .gt. XI) .eqv. LINC) go to 230
            N = 2*ILI + 1
            if (abs(XT(ILI+1)-XI) .lt. abs(XT(ILI) - XI)) then
               ILI = ILI + 1
               N = N + 1
            end if
  240       if (LUP .le. 0) LUP = -ILI
  250       DX(0) = XI - XT(ILI)
c                                  Get bounding indices and interpolate
  260       IUI = ILI
            K = ILI
            L = -1
            if (LEXIT .eq. 0) then
               NDEGI = 0
               if (GETERR) then
                  if (KGO .ne. 3) then
                     if (KGO .eq. 1) NDEGE = NDEGE + 1
                     KEXTRP = min(KEXTRP+1, NDEGE)
                  else
                     NDEGE = -NDEGE
                  end if
               end if
               go to 280
            end if
c                                  Just got index for next XT
  270       TP2 = YT(K)
            if (LBADPT) then
c                                  Check if point should be discarded
               if (TP2 .eq. BADPT) then
                  if (IDIR .eq. 0) then
                     IDIR = 1
                     if (INDXED) then
                        if (XT(2) .lt. 0.E0) IDIR = -IDIR
                     else if (XT(1) .gt. XT(2)) then
                        IDIR = -IDIR
                     end if
                  end if
                  N = N - IDIR * int(sign(1.E0, DX(L+1)))
                  go to 1000
               end if
            end if
  280       L = L + 1
            IDX(L) = K
  290       if (L .eq. 0) then
                  PI(1) = DX(0)
                  TP1 = TP2
            else
               PI(L+1) = PI(L) * DX(L)
               if (L .le. NDEGI) then
c                                 Get divided differences & interpolate.
                  do 300 I = 1, L
                     TP2 = (TP2 - CY(I-1)) / (DX(I-1) - DX(L))
  300             continue
                  TP1 = TP1 + PI(L) * TP2
               end if
            end if
            CY(L) = TP2
            if (L .lt. NDEGE) go to 1000
  340       if (LEXIT .eq. 0) go to 2110
            go to (500, 600, 900, 800), KGO
            if (L .ge. NDEGQ) go to (500,600,900,800), KGO-10
c                        Already got the points selected.
  400       L = L + 1
            if (YTORD) then
               TP2 = YT(L+1)
            else
               TP2 = YT(IDX(L))
            end if
            go to 290
c Got simple interpolated value.
  500       Y = TP1
            if (.not. GETERR) go to  2000
            NDEGI = NDEGI + 1
            KGO = KGO + 1
            if (NDEGE .ge. 0) go to 1000
            go to 400
c Got info. for error estimate.
  600       PI(0) = 1.E0
            if (L .eq. 0) then
               IIFLG = -2
            else
               ERREST=1.5E0*(abs(TP1-Y)+.03125E0*abs(PI(L-1)*CY(L-1)))
               L = L - 1
            end if
            go to  2000
c C1 interpolant.
  800       do 810 I = 1, L-1
               TP2 = (TP2 - CY(I-1)) / (DX(I-1) - DX(L))
  810       continue
            TP2 = (TP2 - CY(L-1)) / (DX(0) - DX(1))
            CY(L) = TP2
            PI(L) = PI(L-1) * DX(0)
            Y = TP1 + PI(L) * TP2
            if (GETERR) ERREST = 1.5E0 * abs(PI(L-1)) * abs(((DX(L-1) *
     1         (DX(0) - DX(1)) / (DX(L-1) - DX(L)) - DX(0)) * TP2) +
     2         abs(.03125E0*CY(L-1)))
            go to 2000
c    Variable order, check estimated error and for convergence.
  900       continue
            if (KAOS .ge. 3) go to 930
            if (KAOS .eq. 2) go to 920
c                                        First time
            KAOS = 2
            E2L = abs(TP1)
            E2 = EBND + 1.E30
            go to 950
c                                        Second time
  920       KAOS = 4
            EBNDI = .66666E0*(EBND + EBNDR*(abs(CY(0))+abs(YT(IDX(1)))))
            EF = DX(0)
            if (LEXIT .gt. 1) then
               EF = DX(1) - DX(0)
            end if
            E2 = abs(EF * TP2)
            EM = max(.75E0, E2 / (E1 + E2 + 1.E-20))
            go to 950
c                                        Usual case            
  930       E2L = E2
            E1 = E2L * (5.E0 * EM / real(L))
            EF = EF * DX(L-1)
            E2 = abs(EF*TP2)
            EM = 0.5E0 * EM + E2 / (E2L + E2 + 1.E-20)
            if (E2 .ge. E1) then
               if (KAOS .eq. 3) then
c                          Apparently diverging, so quit.
                  L = L - 1
                  go to 960
               else
                  KAOS = 3
               end if
            else
               KAOS = 4
            end if
  950       YL = TP1
            if (E2L + E2 .gt. EBNDI) then
               if ((L+NDEGE .lt. 0) .and. (IIFLG .ne. 2)) go to 970
               if (KGO .eq. 13) then
c                                 May need early exit to get next point.
                  IOPT(2) = 0
                  if (L .lt. NDEG) go to 2110
               end if
            end if
  960       TP2 = 1.5E0
            if (LEXIT .gt. 1) TP2 = 1.5E0 * abs(DX(0) / (DX(1) - DX(0)))
            ERREST = TP2 * (E2  + .0625E0 * E2L)
            if (LNDEG .ne. 0) IOPT(LNDEG) = L
            Y = YL
            go to 2000
  970       if (KGO .gt. 10) go to 400
c 
 1000       KK = min(N - IUI - ILI, 2) - 1
c In this section of code, KK=: 1, decrease ILI only; 0, increase IUI;
c -1, decrease ILI; and <-1, increase IUI only
            if (abs(KK) .eq. 1) then
               ILI = ILI - 1
               K = ILI
               if (ILI .ne. 0) then
                  if (INDXED) then
                     XL = XL + H
                     DX(L+1) = XL
                     go to 270
                  end if
                  DX(L+1) = XI - XT(K)
                  if (XT(ILI+1) .ne. XT(ILI)) go to 270
               end if
               if (KK .eq. 1) go to 1100
               N = 0
            else
               IUI = IUI + 1
               K = IUI
               if (IUI .le. NTABI) then
                  if (INDXED) then
                     XU = XU - H
                     DX(L+1) = XU
                     go to 270
                  end if
                  DX(L+1) = XI - XT(K)
                  if (XT(IUI-1) .ne. XT(IUI)) go to 270
               end if
               if (KK .ne. 0) go to 1100
               N = 3*NTABI
            end if
            if (KGO .lt. 4) go to 1000
            KGO = 1
            NDEGE = NDEGE - 1
            if (L .lt. NDEGE) go to 1000
            go to 340
c                                     No more data accessible.
 1100       NDEGI = min(NDEGI, L)
            NDEGE = 0
            IIFLG = 2
            go to 340
c
c                                     Secant start, then use sequential
 1200    continue
            if (XT(1) .eq. XT(NTABI)) go to 1220
            ILI = max(1, min(NTABI, int(1.5E0+real(NTABI-1)*(XI-XT(1))/
     1         (XT(NTABI) - XT(1)))))
         go to 210
c
c                         Special cases
c                                  1 entry in XT
 1220    ILI = 1
         IUI = 1
         KGO = 1
         K = 1
         if (NDEGE .ne. 0) IIFLG = 2
         NDEGE = 0
         if (NTABI .eq. 1) go to 250
c                         Error -- XT(1) .eq. XT(NTAB), and NTAB .ne. 1
         IIFLG = -3
         MACT(2) = NTABI
         ERRDAT(1) = XT(1)
         KGO = 0
         go to 2120
c                                  Extrapolating
 1240    IIFLG = 1
         N = 6 * (ILI - 1)
         if (KGO .ge. 4) KGO = 1
         NDEGI = KEXTRP
         NDEGE = sign(NDEGI, NDEGE)
         if (INDXED) go to 260
         go to 240
c
c                                      Index search, then exit
 1300    continue
         INDXED = .true.
         H = XT(2)
         if (H .eq. 0) then
            IIFLG = -4
            go to 2120
         end if 
         TP1 = 1.E0 + (XI - XT(1)) / H
         ILI = min(max(1, int(TP1+0.5E0)), NTABI)
         XU = (TP1 - real(ILI)) * H
         XL = XU
         DX(0) = XU
         N = ILI + ILI
         if (TP1 .gt. real(ILI))  N = N + 1
         if (TP1 .lt. 1.E0) go to 1240
         if (TP1 .le. real(NTABI)) go to 260
         go to 1240
c
c                                 Already set up
 1400    KGO = KGOC
         if (KGO .eq. 0) then
            IIFLG = -6
            go to 2120
         end if
         YTORD = KGO .gt. 0
         KGO = abs(KGO)
         if (KGO .lt. 10) KGO = KGO + 10
         if (KGO .eq. 12) KGO = 11
         NDEGI = NDEGEC
         if (GETERR) then
            if (KGO .eq. 11) NDEGI = NDEGI - 1
         end if
         NDEGQ = NDEGI
         if (KGO .eq. 13) NDEGQ = 0
         NDEGE = -NDEGI
         L = -1
         if (KGO .eq. 14) NDEGI = NDEGI - 1
         go to 400
c
c                                 Set number of points for extrapolation
 1500    continue
            I = I + 1
            KEXTRP = min(IOPT(I), NDEGE)
            go to 10
c
c                                 Get derivatives of interpolant
 1600    continue
            I = I + 2
            LDERIV = IOPT(I-1)
            LEXIT = IOPT(I) + 1
            go to 10
c
c                                     Get expected errors in YT.
1900     continue
            I = I + 1
            LEXERR = IOPT(I)
c
c                                     Set to get error estimate.
 1930    continue
            GETERR = .true.
            go to 10
c
c                                     Set for automatic order selection.
1950     continue
            I = I + 2
            LNDEG = I
            EBND = max(0.E0, EOPT(IOPT(I-1)))
            EBNDR = max(0.E0, EOPT(IOPT(I-1)+1))
            KGO = 3
            KAOS = 1
            NDEGI = min(MAXDEG, NDEGE+1)
            NDEGE = -NDEGI
         go to 1930
c
c                                    Set to ignore special points
1970     continue
            LBADPT = .true.
            I = I + 1
            BADPT = EOPT(IOPT(I))
            IDIR = 0
            go to 10
c
c                                    Put any new options just above here
c                    End of the loop
c
 2000 if (LEXIT .lt. 2) go to 2100
c                       Compute derivatives of Y
      if (KGO .gt. 3) then
         DXS = DX(L-1)
         DX(L-1) = DX(0)
      end if
      N = LEXIT - 1
      if (N .gt. L) then
        if (N .gt. MAXDEG) then
           MACT(2) = N
           MACT(4) = MAXDEG
           N = MAXDEG
           IIFLG = -5
        end if
        do 2020 I = L+1, N
           EOPT(I+LDERIV-1) = 0.E0
 2020   continue
        N = L
      end if
      TP1 = CY(1)
      PID(0) = 1.E0
      do 2030 I = 1, L-1
         PID(I) = PI(I) + DX(I) * PID(I-1)
         TP1 = TP1 + PID(I) * CY(I+1)
 2030 continue
      EOPT(LDERIV) = TP1
      do 2070 K = 2, N
         TP1 = CY(K)
         do 2060 I = 1, L-K
            PID(I) = PID(I) + DX(I+K-1) * PID(I-1)
            TP1 = TP1 + PID(I) * CY(I+K)
 2060    continue
         EOPT(K+LDERIV-1) = TP1
 2070 continue
      if (KGO .gt. 3) DX(L-1) = DXS
c                                    Save info. and return
 2100 continue
      if (GETERR) then
         if (EPSR .eq. 0.E0) then
            EPSR = R1MACH(4)
         end if
         TP1 = EPSR
         TP2 = 0.E0
         if (LEXERR .ne. 0) then
            TP2 = max(0.E0, EOPT(LEXERR))
            TP1 = max(TP1, EOPT(LEXERR+1))
         end if
         if (IIFLG .eq. 2) ERREST = 32.E0 * ERREST
         EOPT(1) = ERREST + TP2 + TP1*(abs(CY(0)) + abs(DX(0)*CY(1)))
         if ((KGO .eq. 3) .or. (KGO .eq. 13)) then
            if (EBNDI .ne. 0.E0) then
               if (EOPT(1) .gt. EBND) then
                  ERRDAT(1) = EOPT(1)
                  ERRDAT(2) = EBND
                  IIFLG = -1
               end if
            end if
         end if
      end if
 2110 KGOC = -KGO
      NDEGEC = L
 2120 IOPT(1) = IIFLG
      if (IIFLG .ge. 0) return
c Take care of error messages, IIFLG < -2 should stop.
      MACT(6) = 88
      if (IIFLG .ge. -2) MACT(6) = 24 - IIFLG
      if (IOPT(2) .ge. 254)  then
         if (IIFLG .eq. -1) return
         MACT(6) = 28
      end if
      MACT(7) = -IIFLG
      MACT(8) = MLOC(-IIFLG)
      call SMESS(MACT, MTXTAA, IOPT, ERRDAT)
      MACT(9) = MERET
      return
c                       End of Interpolation routine -- SILUP
      end

      subroutine SMESS (MACT, TEXT, IDAT, FDAT)
c     .  Copyright (C) 1991, California Institute of Technology.
c     .  All rights reserved.  U. S. Government sponsorship under
c     .  NASA contract NAS7-918 is acknowledged.
c++   Set VERSION = F
c>> 1993-05-14 SMESS Krogh  Changed TEXT to array of character strings.
c>> 1993-04-14 SMESS Krogh  Fixes for conversion to C. (C%% comments.)
c>> 1992-07-12 SMESS Krogh  Fixed so negative KDFDEF works.
c>> 1992-05-27 SMESS Krogh  Initialized LDFDEF in a data statement.
c>> 1992-05-14 SMESS Krogh  Put common blocks in save statement.
c>> 1992-04-28 SMESS Krogh  Corrected minor error in floating pt. format
c>> 1992-02-07 SMESS Krogh  Initial Code.
c
c Processes Messages -- Actions are controlled by MACT().  See
c comment is subroutine MESS.  This program is for the extra
c argument of type real.
c
c BUF    In common CMESSC, see MESS.
c DOLS   In common for intitialization, not used here.  See MESS.
c EUNIT  In common for intitialization, not used here.  See MESS.
c FDAT   Formal argument -- gives floating point data to print.
c FMAG   Magnitude of floating point number to output, with negative sign
c   if floating point number is < 0.
c FMIN   Minimum floating point number.
c FMTF   In common CMESSC, format for printing floating point number.
c FMTG   In common CMESSC, user format to use in place of FMTF.
c FNMAX  Maximum negative floating point number.
c FOUT   Floating point number to be output.
c FSMA   Smallest postitive floating point number.
c I1GETD Parameter giving the value to pass to I1MACH to get the number
c   of base b digits in a floating point number.  (=11 for single
c   precision code and =14 for double precision code.)
c I1MACH External function giving integer info. about the environment.
c ICOL   In common CMESSI, see MESS.
c ID     Number of decimal digits for floating point format statement.
c IDAT   Integer data -- passed to MESS.
c IVAR   In common CMESSI, see MESS.
c IWF    In common CMESSI, see MESS.
c IWG    In common CMESSI, see MESS.
c J      Temporary index.
c K      Temporary index.
c KAP    Number of extra 0's after the decimal point in "F" format.  If
c        < 0, -KAP gives the number of extra digits to the left of the
c        decimal point.  KAP depends on abs(smallest number printed).
c KBP    Number of extra digits before the decimal point required by the
c        largest number to be printed.
c KDF    In common CMESSI, see MESS.
c KDFDEF In common CMESSI, see MESS.
c KDIAG  In common CMESSI, not used here, see MESS.
c KEXE   Extra space required for E format.
c KF     In common CMESSI, see MESS.
c KLINE  In common CMESSI, see MESS.
c KSCRN  In common CMESSI, see MESS.
c KRES1  In common CMESSI, see MESS.
c KSPEC  In common CMESSI, see MESS.
c LASTER In common CMESSI, not used here, see MESS.
c LASTI  In common CMESSI, see MESS.
c LBUF   In common CMESSI, see MESS.
c LDFDEF Value of KDFDEF for this routine.  (Saved)
c LENBUF In common CMESSI, see MESS.
c LENLIN In common CMESSI, not used here, see MESS.
c LENTRY In common CMESSI, see MESS.
c LHEAD  In common CMESSI, not used here, see MESS.
c LINERR In common CMESSI, not used here, see MESS.
c LINMSG In common CMESSI, not used here, see MESS.
c LOCBEG In common CMESSI, see MESS.
c LPRINT In common CMESSI, not used here, see MESS.
c LSTOP  In common CMESSI, not used here, see MESS.
c LSTRT  In common CMESSI, see MESS.
c MACT   Formal argument, see MESS.
c MDAT   In common CMESSI, not used here, see MESS.
c MEMDA5 In common CMESSI, see MESS.
c MESS   Program called for most of the message processing.
c MPT    In common CMESSI, see MESS.
c MUNIT  In common CMESSI, not used here, see MESS.
c NCOL   In common CMESSI, see MESS.
c NDIM   In common CMESSI, see MESS.
c NFDAT  In common CMESSI, see MESS.
c NIDAT  In common CMESSI, not used here, see MESS.
c NMDAT  In common CMESSI, not used here, see MESS.
c NTEXT  In common CMESSI, not used here, see MESS.
c OUNIT  In common CMESSI, not used here, see MESS.
c R1MACH External func. giving floating pt. info. about the environment.
c SUNIT  In common CMESSI, not used here, see MESS.
c TEXT   Formal argument, passed to MESS, see there.
c XARGOK In common CMESSI, see MESS.
c
c I1GETD is 11 for real, and 14 for double precision.
      integer I1GETD
      parameter (I1GETD = 11)
      integer          MACT(*), IDAT(*)
      real             FDAT(*)
      character        TEXT(*)*(*)
      integer          I1MACH, ICOL, ID, J, K, KAP, KBP, KEXE, LDFDEF
      real             FMAG, FMIN, FOUT, FNMAX, FSMA, R1MACH
      save LDFDEF
C     COMMENTED OUT BY KUMAZAWA
C     save /CMESSI/, /CMESSC/
c++ CODE for VERSION = C is inactive
c      integer  kciwid, kccwid, kcrwid, lfprec, lgprec
c      common /MESSCC/ kciwid, kccwid, kcrwid, lfprec, lgprec
c++ END
c
c ************************** Data from common block ********************
c
c For comments on these variables, see the listing for MESS.
c
      integer   LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF = 250)
      parameter (MEVBAS = 10)
      parameter (MEVLAS = 32)
      logical          XARG, GOTFMT, XARGOK
      integer          EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS), IMAG,
     1   ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ, KLINE, KSCRN,
     2   KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN, LENOUT,
     3   LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT, LSTOP,
     4   LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL, NDIM,
     5   NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
c
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*15, FMTG*15,
     1  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,
     1   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT,
     2   NFDAT, NMDAT, MDAT, TABSPA, ICHAR0, IMAG, IRC, ITEXT, IWF,
     3   IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI, LBUF, LENLIN,
     4   LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT, MAXWID, MPT,
     5   NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
c
C     COMMENTED OUT BY KUMAZAWA
C     save /CMESSI/, /CMESSC/

c
      data LDFDEF / 0 /
c
c ************************* Start of Executable Code *******************
c
      XARGOK = .true.
      if (LDFDEF .eq. 0) then
         LDFDEF = 1 + int(real(I1MACH(I1GETD)) * R1MACH(5))
      end if
      KDFDEF = LDFDEF
      KDF = KDFDEF
   10 call MESS (MACT, TEXT, IDAT)
      go to (20, 100, 100, 200, 300, 400), LENTRY-3
      XARGOK = .false.
      LDFDEF = KDFDEF
      return
c                                      Print from FDAT
   20 J = LBUF + 1
      FOUT = FDAT(NFDAT)
      NFDAT = NFDAT + 1
      if (KSPEC .ge. 8) then
         LBUF = LBUF + IWG
c++ CODE for VERSION = F is active
         write (BUF(J:LBUF), FMTG) FOUT
c++ CODE for VERSION = C is inactive
c%%      sprintf(&cmessc.buf[j-1], cmessc.fmtg, cmessi.iwg,
c%%         messcc.lgprec, fout);
c++ END
         go to 10
      end if
      FMAG = FOUT
      FSMA = abs(FOUT)
      IWF = 1
c                                      Get the format.
   40 KEXE = 3
      KBP = 0
      KAP = 0
      if (FMAG .ne. 0.E0) then
         if (FMAG .lt. 0.E0) then
            FMAG = -FMAG
            IWF = IWF + 1
         end if
         if (FSMA .ge. 1.E0) then
            KAP = -log10(FSMA) - 1.000001E0
         else
            KAP = -log10(FSMA) + 1.E-6
         end if
         if (FMAG .ne. FSMA) then
            if (FMAG .ge. 1.E0) KBP = 1 + int(log10(FMAG) + 1.E-6)
         else
            KBP = max(-KAP, 0)
         end if
         K = max(KAP+1, KBP)
         if (K .ge. 10) KEXE = 3 + int(log10(real(K)) + 1.e-5)
      end if
      if (KDF .eq. 0) KDF = KDFDEF
      if (KDF .gt. 0) then
         if (KBP + max(KAP, 0) .ge. KEXE) go to 50
         ID = max(KDF + KAP, 0)
         IWF = IWF + ID + KBP
c If LENTRY is 5 or 6 (Vector or Matrix output), and KBP is 0, need 1
c extra since with the extra space for a blank an extra 0 is printed.
         if (LENTRY - 2*KBP .gt. 4) IWF = IWF + 1
      else
         if (KBP .gt. KEXE) go to 50
         ID = -KDF
         IWF = IWF + KBP + ID
      end if
c++ CODE for VERSION = F is active
      write (FMTF, '(6H(0P99F,I2,1H.,I2,1H))') IWF,ID
c++ CODE for VERSION = C is inactive
c%%    strcpy(cmessc.fmtf, "%*.*f\0");
c      lfprec = id
c++ END
      go to 60
   50 ID = KDF - 1
      if (ID .lt. 0) ID = -ID - 1
c++ CODE for VERSION = F is active
      IWF = IWF + ID + KEXE + 1
      write (FMTF, '(6H(1P99E,I2,1H.,I2,1HE,I1,1H))') IWF,ID,KEXE-2
c++ CODE for VERSION = C is inactive
c      iwf = iwf + id + 5
c%%    strcpy(cmessc.fmtf, "%*.*E\0");
c      lfprec = id
c++ END
   60 if (LENTRY .ne. 4) go to 10
c
      LBUF = LBUF + IWF
c++ CODE for VERSION = F is active
      write (BUF(J:LBUF),FMTF) FOUT
c++ CODE for VERSION = C is inactive
c%%      sprintf(&cmessc.buf[j-1], cmessc.fmtf, cmessi.iwf,
c%%        messcc.lfprec, fout);
c++ END
      go to 10
c                                     Get format for a vector or matrix
  100 if (KDF .eq. 0) KDF = KDFDEF
      ICOL = 1
      FMAG = FDAT(LOCBEG)
      FMIN = FMAG
      FSMA = 1.E0
      FNMAX = -1.E0
  110 do 120 J = LOCBEG, LASTI
         if (FDAT(J) .le. 0.E0) then
            FMIN = min(FMIN, FDAT(J))
            if (FDAT(J) .ne. 0.E0) FNMAX = max(FNMAX, FDAT(J))
         else
            FMAG = max(FMAG, FDAT(J))
            FSMA = min(FSMA, FDAT(J))
         end if
  120 continue
      if (NCOL .ne. 0) then
         ICOL = ICOL + 1
         LOCBEG = LOCBEG + NDIM
         LASTI = LASTI + NDIM
         if (ICOL .le. NCOL) go to 110
      end if
      if (FMIN .lt. 0.E0) then
         FMAG = min(FMIN, -0.1E0*FMAG)
         FSMA = min(-FNMAX, FSMA)
      end if
      if (FSMA .ge. 1.E0) then
         FSMA = abs(FMAG)
      else if (abs(FMAG) .lt. 1.E0) then
         FMAG = sign(FSMA, FMAG)
      end if
      IWF = 2
      go to 40
c                                    Floating point vector output
  200 continue
c++ CODE for VERSION = F is active
      write (BUF(LSTRT:LBUF), FMTF) (FDAT(K), K = MPT, MPT+KLINE-1)
c++ CODE for VERSION = C is inactive
c%%      for (k=cmessi.mpt; k<cmessi.mpt+cmessi.kline; k++)
c%%         sprintf(&cmessc.buf[cmessi.lstrt+cmessi.iwf*(k-cmessi.mpt)
c%%         -1], cmessc.fmtf, cmessi.iwf, messcc.lfprec, fdat[k-1]);
c++ END
      MPT = MPT + KLINE
      go to 10
c                                    Floating point matrix output
  300 continue
c++ CODE for VERSION = F is active
      write (BUF(LSTRT:LBUF), FMTF) (FDAT(K), K = MPT, LASTI, NDIM)
c++ CODE for VERSION = C is inactive
c%%      j = 0;
c%%      for (k=cmessi.mpt; k<=cmessi.lasti; k+=cmessi.ndim)
c%%         sprintf(&cmessc.buf[cmessi.lstrt+cmessi.iwf*(j++)-1],
c%%            cmessc.fmtf, cmessi.iwf, messcc.lfprec, fdat[k-1]);
c++ END
      go to 10
c                                    Table output
  400 continue
c++ CODE for VERSION = F is active
      write (BUF(LSTRT:LBUF), FMTF) (FDAT(K), K = MPT, MPT+KLINE-1)
c++ CODE for VERSION = C is inactive
c%%      for (k=cmessi.mpt; k<cmessi.mpt+cmessi.kline; k++)
c%%         sprintf(&cmessc.buf[cmessi.lstrt+cmessi.iwf*(k-cmessi.mpt)
c%%         -1], cmessc.fmtf, cmessi.iwf, messcc.lfprec, fdat[k-1]);
c++ END
      go to 10
      end

      subroutine UMESS(TEXT, MACT, IVAR)
c
c>> 1992-01-29 UMESS  Krogh  Initial code.
c
c Dummy routine called by the error message routine MESS.  This gives
c the default action when an error message is requested, which is to do
c nothing here and take the default actions in MESS.  One can change the
c actions taken by changing the values in MACT or IVAR.  Comments here
c refer to things defined in the listing for MESS.
c Variables in the calling sequence are defined as follows:
c
c TEXT  Same as TEXT passed into MESS.  The characters up to the first
c       "$" sign are the name of the subroutine that made the call to
c       MESS.
c MACT  Contents of MACT connected with the error message.  The first
c       location is the value of K52 = 10*s + p, where s is the stop
c       level specified for the error message, and p the print level.
c       The next is the value of L52, which is an error index associated
c       with an error message, and the next is the value of M52.  These
c       values or a subset of them together with TEXT should serve to
c       uniquely identify the error message.  The values of MACT(1)
c       (=K52), or  MACT(2) (=L52) can be changed and the new values
c       will be used if a return is made to MESS.  (For example MACT(1)
c       = MACT(1) - mod(MACT(1), 10), would turn off printing of a
c       message.)  Note that if MACT(1) > 90, and it is replaced with a
c       value which results in MESS returning rather than stopping,
c       unpredictable results are to be expected.
c IVAR  An array containing the parameters defining default values for
c       MESS.  The parameter values below can be used to access data in
c       this array.  For example, IVAR(MESTOP) gives the current stop
c       value, so one can tell (using s from MACT(1) above) if an error
c       is going to stop, and perhaps do some clean up before the stop
c       occurs.  Values in IVAR can be changed, but this should be done
c       with great care, since if MESS doesn't stop, the next error or
c       diagnostic message will use the new value, and diagnostic
c       messages do not call this routine.
c
      integer   MESUNI, MEHEAD, MEDDIG, MEMLIN, MEELIN, MEMUNI, MEEUNI,
     1  MESCRN, MEDIAG, MEMAXE, MESTOP, MEPRNT, METDIG, MENTXT, MEIDAT,
     2  MEFDAT, MEMDAT, MEMDA1, MEMDA2, MEMDA3, MEMDA4, MEMDA5, METABS,
     3  MECONT, MERET , MEEMES, METEXT, METXTF, METABL, MERES3, MEIVEC,
     4  MEIMAT, MEJVEC, MEJMAT, MEFVEC, MEFMAT, MEGVEC, MEGMAT, MEMAXI,
     5  MEGBAS, MEVBAS, MEVLAS, MESLAS
      parameter (MESUNI=10,MEHEAD=11,MEDDIG=12,MEMLIN=13,MEELIN=14,
     1 MEMUNI=15,MEEUNI=16,MESCRN=17,MEDIAG=18,MEMAXE=19,MESTOP=20,
     2 MEPRNT=21,METDIG=22,MENTXT=23,MEIDAT=24,MEFDAT=25,MEMDAT=26,
     3 MEMDA1=27,MEMDA2=28,MEMDA3=29,MEMDA4=30,MEMDA5=31,METABS=32)
c
      character*(*) TEXT(*)
      parameter (MEVBAS = 10, MEVLAS = 32)
      integer   MACT(*), IVAR(MEVBAS:MEVLAS)
c
c Code for your actions goes here.
c
      return
      end

************************************************************************

      real function SGAMMA(X)
C>> 1991-10-21 SGAMMA CLL Eliminated DGAM1 as a separate subroutine.
C>> 1991-01-16 SGAMMA Lawson  Replaced D2MACH/R2MACH with DGAM1.
C>> 1985-08-02 SGAMMA Lawson  Initial code.
C 
C-----------------------------------------------------------------------
C 
C  THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A double precision
C      ARGUMENT X. PERMITS NEGATIVE AS WELL AS POSITIVE X. NOTE
C      THAT THE GAMMA FUNCTION HAS POLES AT ZERO AND AT NEGATIVE
C      ARGUMENTS. COMPUTATION IS BASED ON AN ALGORITHN OUTLINED IN
C      W.J.CODY, 'AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
C      FUNCTIONS', LECTURE NOTES IN MATHEMATICS, 506, NUMERICAL ANALYSIS
C      DUNDEE, 1975, G. A. WATSON (ED.),SPRINGER VERLAG, BERLIN, 1976.
C      THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
C      FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS. COEFFICIENTS
C      FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
C      THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM HART, ET. AL.,
C      COMPUTER APPROXIMATIONS, WILEY AND SONS, NEW YORK, 1968.
C      LOWER ORDER APPROXIMATIONS CAN BE SUBSTITUTED FOR THESE ON
C      MACHINES WITH LESS PRECISE ARITHMETIC.
C 
C  Designed & programmed by W.J.CODY, Argonne National Lab.,1982.
C  Minor changes for the JPL library by C.L.LAWSON & S.CHAN,JPL,1983.
C 
C***********************************************************************
C 
C  EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
C 
C  EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
C           1.0 + EPS .GT. 1.0  (EPS = [D/R]1MACH(4).)
C  XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER.
c           XINF = [D/R]1MACH(2).
C  XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
c           both XMININ and 1/XMININ are representable.  This is
c           the larger of [D/R]1MACH(1) and 1.0 / [D/R]1MACH(2).
C  XGBIG  - A value such that    Gamma(XGBIG) = 0.875 * XINF.
c           (Computed and used in [D/S]GAMMA.)
C  XLBIG  - A value such that LogGamma(XLBIG) = 0.875 * XINF.
c           (Computed and used in [D/S]LGAMA.)
C 
C      Values of XINF, XGBIG, and XLBIG for some machines:
C 
c        XINF              XGBIG     XLBIG       Machines
c 
c  2**127  = 0.170e39      34.81  0.180e37     Vax SP & DP; Unisys SP
c  2**128  = 0.340e39      35.00  0.358e37     IEEE SP
c  2**252  = 0.723e76      57.54  0.376e74     IBM30xx DP
c  2**1023 = 0.899e308    171.46  0.112e306    Unisys DP
c  2**1024 = 0.180e309    171.60  0.2216e306   IEEE DP
c  2**1070 = 0.126e323    177.78  0.1501e320   CDC/7600 SP
c  2**8191 = 0.550e2466   966.94  0.8464e2462  Cray SP & DP
c 
C***********************************************************************
C 
C  ERROR RETURNS
C 
C  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
C     WHEN OVERFLOW WOULD OCCUR. THE COMPUTATION IS BELIEVED
C     TO BE FREE OF UNDERFLOW AND OVERFLOW.
C 
C  OTHER SUBPROGRAM REQUIRED
C 
C      real,exp,log,sin
C 
C 
C  AUTHOR:W. J. CODY
C         APPLIED MATHMATICS DIVISION
C         ARGONNE NATIONAL LABORATORY
C         ARGONNE, IL 60439
C 
C  LATEST MODIFICATION by Cody: MAY 18, 1982
C 
c     ------------------------------------------------------------------
c--   D version uses DGAMMA, D1MACH, DERM1, DERV1, dble
c--   S version uses SGAMMA, R1MACH, SERM1, SERV1, real
C     ------------------------------------------------------------------
      real R1MACH
      real C(7), CONST, DEL, EPS, F, FACT, FP, HALF
      real ONE,P(8), PI,Q(8), RES,C1
      real SUM,TEMP, TWELVE,TWO
      real X,X1, X2, XGBIG,XDEN,XINF,XMININ,XNUM
      real Y,Y1,YSQ,Z,ZERO
      integer I,J,N
      logical PARITY
C 
      save EPS, XGBIG, XMININ, XINF
C 
      parameter( ONE = 1.0E0, HALF = 0.5E0, TWO = 2.0e0)
      parameter( ZERO = 0.0E0, TWELVE = 12.0E0)
C 
C                      C1 = LOG base e of SQRT(2*PI)
C 
      parameter( C1 = 0.9189385332046727417803297E0)
C 
      parameter( PI = 3.1415926535897932384626434E0)
C 
      data XINF/0.0E0/
C 
C-----------------------------------------------------------------------
C  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
C     APPROXIMATION OVER (1,2).
C-----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,
     *       -3.79804256470945635097577E+2,6.29331155312818442661052E+2,
     *       8.66966202790413211295064E+2,-3.14512729688483675254357E+4,
     *       -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
      DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,
     *      -1.01515636749021914166146E+3,-3.10777167157231109440444E+3,
     *        2.25381184209801510330112E+4,4.75584627752788110767815E+3,
     *      -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
C-----------------------------------------------------------------------
C  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
C-----------------------------------------------------------------------
      DATA C/-1.910444077728E-03,8.4171387781295E-04,
     *     -5.952379913043012E-04,7.93650793500350248E-04,
     *     -2.777777777777681622553E-03,8.333333333333333331554247E-02,
     *      5.7083835261E-03/
C-----------------------------------------------------------------------
C 
      IF (XINF .EQ. ZERO) THEN
        EPS = R1MACH(4)
        XINF = R1MACH(2)
         if(R1MACH(1) * R1MACH(2) .ge. ONE) then
            XMININ = R1MACH(1)
         else
            XMININ = ONE / R1MACH(2)
         endif
c 
c                         Compute XGBIG
c 
c        XGBIG will satisfy Gamma(XGBIG) = 0.875 * XINF.
c        Use a Newton iteration and the following approximation for
c        the gamma function:
c        log(gamma(x)) ~ (x - .5)*log(x) - x + 0.5 * log(2 * PI)
c 
         TEMP = log(0.875e0 * XINF)
         CONST = HALF * log(TWO * PI) - TEMP
         X1 = TEMP * 0.34e0
         do 40 J=1,7
            F = (X1-HALF)*log(X1) - X1 + CONST
            FP = ((X1-HALF)/X1)  + log(X1) - ONE
            DEL = -F/FP
            X2 = X1+DEL
            if(abs(DEL) .lt. 0.5e-5 * X2) go to 45
            X1 = X2
   40    continue
   45    continue
         XGBIG = X2
      END IF
      PARITY = .FALSE.
      FACT = ONE
      N = 0
      Y = X
      IF (Y .GT. ZERO) GO TO 200
C-----------------------------------------------------------------------
C  ARGUMENT IS NEGATIVE OR ZERO
C-----------------------------------------------------------------------
      Y = -X
      J = INT(Y)
      RES = Y - real(J)
      IF (RES .EQ. ZERO) GO TO 700
      IF (J .NE. (J/2)*2) PARITY = .TRUE.
      FACT = -PI / sin(PI*RES)
      Y = Y + ONE
C-----------------------------------------------------------------------
C  ARGUMENT IS POSITIVE
C-----------------------------------------------------------------------
  200 IF (Y .LT. EPS) GO TO 650
      IF (Y .GE. TWELVE) GO TO 300
      Y1 = Y
      IF (Y .GE. ONE) GO TO 210
C-----------------------------------------------------------------------
C  0.0 .LT. ARGUMENT .LT. 1.0
C-----------------------------------------------------------------------
      Z = Y
      Y = Y + ONE
      GO TO 250
C-----------------------------------------------------------------------
C  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
C-----------------------------------------------------------------------
  210 N= int(Y) - 1
      Y = Y - real(N)
      Z = Y - ONE
C-----------------------------------------------------------------------
C  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
C-----------------------------------------------------------------------
  250 XNUM = ZERO
      XDEN = ONE
      DO 260 I = 1, 8
         XNUM = (XNUM + P(I)) * Z
         XDEN = XDEN * Z + Q(I)
  260 CONTINUE
      RES = (XNUM / XDEN + HALF) + HALF
      IF (Y .EQ. Y1) GO TO 900
      IF (Y1 .GT. Y) GO TO 280
C-----------------------------------------------------------------------
C  ADJUST RESULT FOR CASE 0.0 .LT. ARGUMENT .LT. 1.0
C-----------------------------------------------------------------------
      RES = RES / Y1
      GO TO 900
C-----------------------------------------------------------------------
C  ADJUST RESULT FOR CASE 2.0 .LT. 12.0
C-----------------------------------------------------------------------
  280 DO 290 I = 1, N
         RES = RES * Y
         Y = Y + ONE
  290 CONTINUE
      GO TO 900
C-----------------------------------------------------------------------
C  EVALUATE FOR ARGUMENT .GE. 12.0
C-----------------------------------------------------------------------
  300 IF (Y .GT. XGBIG) GO TO 720
      YSQ = Y * Y
      SUM = C(7)
      DO 350 I = 1, 6
         SUM = SUM / YSQ + C(I)
  350 CONTINUE
      SUM = ((SUM/Y + C1) - Y) + (Y - HALF) * log(Y)
      RES = exp(SUM)
      GO TO 900
C-----------------------------------------------------------------------
C  ARGUMENT .LT. EPS
C-----------------------------------------------------------------------
  650 IF (Y .LT. XMININ) GO TO 740
      RES = ONE / Y
      GO TO 900
C-----------------------------------------------------------------------
C  RETURN FOR SINGULARITIES,EXTREME ARGUMENTS, ETC.
C-----------------------------------------------------------------------
  700 CALL SERM1('SGAMMA',1,0,'POLE AT 0 AND NEG INTEGERS',
     *           'X',X,'.')
      GO TO 780
C 
  720 CALL SERM1('SGAMMA',2,0,'X SO LARGE VALUE WOULD OVERFLOW',
     *           'X',X,',')
      CALL SERV1('LIMIT',XGBIG,'.')
      GO TO 780
C 
  740 CALL SERM1('SGAMMA',3,0,'X TOO NEAR TO A SINGULARITY.'//
     *           'VALUE WOULD OVERFLOW.','X',X,'.')
C 
  780 SGAMMA = XINF
      GO TO 950
C-----------------------------------------------------------------------
C  FINAL ADJUSTMENTS AND RETURN
C-----------------------------------------------------------------------
  900 IF (PARITY) RES = -RES
      IF (FACT .NE. ONE) RES = FACT / RES
      SGAMMA = RES
  950 RETURN
      END

************************************************************************

      SUBROUTINE SERM1(SUBNAM,INDIC,LEVEL,MSG,LABEL,VALUE,FLAG)
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1985-09-25 SERM1  Lawson  Initial code.
C
      CHARACTER*(*) SUBNAM,MSG,LABEL
      CHARACTER*1 FLAG
      CALL ERMSG(SUBNAM,INDIC,LEVEL,MSG,',')
      CALL SERV1(LABEL,VALUE,FLAG)
C
      RETURN
      END

************************************************************************

      SUBROUTINE SERV1(LABEL,VALUE,FLAG)
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1987-09-23 SERV1  Lawson  Initial code.
C
C     ------------------------------------------------------------
C     SUBROUTINE ARGUMENTS
C     --------------------
C     LABEL     An identifing name to be printed with VALUE.
C
C     VALUE     A real number to be printed.
C
C     FLAG      See write up for FLAG in ERMSG.
C
C     ------------------------------------------------------------
C
      COMMON/M77ERR/IDELTA,IALPHA
      INTEGER IDELTA,IALPHA
      CHARACTER*(*) LABEL
      CHARACTER*1 FLAG
      SAVE /M77ERR/
C
      IF (IALPHA.GE.-1) THEN
        WRITE (*,1002) LABEL,VALUE
        IF (FLAG.EQ.'.') CALL ERFIN
      ENDIF
      RETURN
C
 1002 FORMAT(3X,A,' = ',G16.7)
      END

************************************************************************

      SUBROUTINE ERFIN
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1985-09-23 ERFIN  Lawson  Initial code.
C
      COMMON/M77ERR/IDELTA,IALPHA
      SAVE /M77ERR/
C
      PRINT 1003
      IF (IALPHA.GE.2) STOP
      RETURN
 1003 FORMAT(1X,72('$')/' ')
      END

************************************************************************

      SUBROUTINE ERMSG(SUBNAM,INDIC,LEVEL,MSG,FLAG)
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1992-10-20 ERMSG  WV Snyder  added ERLSET, ERLGET
C>> 1985-09-25 ERMSG  Lawson  Initial code.
C
C     --------------------------------------------------------------
C
C     Four entries: ERMSG, ERMSET, ERLGET, ERLSET
C     ERMSG initiates an error message. This subr also manages the
C     saved value IDELOC and the saved COMMON block M77ERR to
C     control the level of action. This is intended to be the
C     only subr that assigns a value to IALPHA in COMMON.
C     ERMSET resets IDELOC & IDELTA.  ERLGET returns the last value
C     of LEVEL passed to ERMSG.  ERLSET sets the last value of LEVEL.
C     ERLSET and ERLGET may be used together to determine the level
C     of error that occurs during execution of a routine that uses
C     ERMSG.
C
C     --------------------------------------------------------------
C     SUBROUTINE ARGUMENTS
C     --------------------
C     SUBNAM   A name that identifies the subprogram in which
C              the error occurs.
C
C     INDIC    An integer printed as part of the mininal error
C              message. It together with SUBNAM can be used to
C              uniquely identify an error.
C
C     LEVEL    The user sets LEVEL=2,0,or -2 to specify the
C              nominal action to be taken by ERMSG. The
C              subroutine ERMSG contains an internal variable
C              IDELTA, whose nominal value is zero. The
C              subroutine will compute IALPHA = LEVEL + IDELTA
C              and proceed as follows:
C              If (IALPHA.GE.2)        Print message and STOP.
C              If (IALPHA=-1,0,1)      Print message and return.
C              If (IALPHA.LE.-2)       Just RETURN.
C
C     MSG      Message to be printed as part of the diagnostic.
C
C     FLAG     A single character,which when set to '.' will
C              call the subroutine ERFIN and will just RETURN
C              when set to any other character.
C
C     --------------------------------------------------------------
C
C     C.Lawson & S.Chan, JPL, 1983 Nov
C
C     ------------------------------------------------------------------
      INTEGER OLDLEV
      COMMON/M77ERR/IDELTA,IALPHA
      CHARACTER*(*) SUBNAM,MSG
      CHARACTER*1 FLAG
      SAVE/M77ERR/,IDELOC,OLDLEV
      DATA IDELOC/0/, OLDLEV /0/
      OLDLEV = LEVEL
      IDELTA = IDELOC
      IALPHA = LEVEL + IDELTA
      IF (IALPHA.GE.-1) THEN
c
c            Setting FILE = 'CON' works for MS/DOS systems.
c
c
        WRITE (*,1001) SUBNAM,INDIC
        WRITE (*,*) MSG
        IF (FLAG.EQ.'.') CALL ERFIN
      ENDIF
      RETURN
C
 1001 FORMAT('0',72('$')/' SUBPROGRAM ',A,' REPORTS ERROR NO. ',I4)
C
C
      ENTRY ERMSET(IDEL)
      IDELTA=IDEL
      IDELOC=IDEL
      RETURN
C
C
      ENTRY ERLSET (LEVEL)
      OLDLEV = LEVEL
      RETURN
C
C     
      ENTRY ERLGET (LEVEL)
      LEVEL = OLDLEV
      RETURN
      END
	  
      subroutine IERM1(SUBNAM,INDIC,LEVEL,MSG,LABEL,VALUE,FLAG)
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1990-01-18 CLL Added Integer stmt for VALUE.  Typed all variables.
C>> 1985-08-02 IERM1  Lawson  Initial code.
C
      integer INDIC, LEVEL, VALUE
      character*(*) SUBNAM,MSG,LABEL
      character*1 FLAG
      call ERMSG(SUBNAM,INDIC,LEVEL,MSG,',')
      call IERV1(LABEL,VALUE,FLAG)
C
      return
      end
	  
      SUBROUTINE IERV1(LABEL,VALUE,FLAG)
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1985-09-20 IERV1  Lawson  Initial code.
C
C     ------------------------------------------------------------
C     SUBROUTINE ARGUMENTS
C     --------------------
C     LABEL     An identifing name to be printed with VALUE.
C
C     VALUE     A integer to be printed.
C
C     FLAG      See write up for FLAG in ERMSG.
C
C     ------------------------------------------------------------
C
      COMMON/M77ERR/IDELTA,IALPHA
      INTEGER IDELTA,IALPHA,VALUE
      CHARACTER*(*) LABEL
      CHARACTER*1 FLAG
      SAVE /M77ERR/
C
      IF (IALPHA.GE.-1) THEN
        WRITE (*,1002) LABEL,VALUE
        IF (FLAG .EQ. '.') CALL ERFIN
      ENDIF
      RETURN
C
 1002 FORMAT(3X,A,' = ',I5)
      END

      subroutine SFMIN(X,XORF,MODE,TOL)
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1987-12-09 SFMIN  Lawson  Initial code.
c
C     Finds a local minimum of a function f(x) in the closed interval
c     between A and B.
c     The function f(x) is defined by code in the calling program.
c     The user must initially provide A, B, and a tolerance, TOL.
c     The function will not be evaluated at A or B unless the
c     minimization search leads to one of these endpoints.
c
c     This subroutine uses reverse communication, i.e., it returns to
c     the calling program each time it needs to have f() evaluated at a
c     new value of x.
c     ------------------------------------------------------------------
c                   Subroutine Arguments
c
c     MODE, X, XORF  [all are inout]
c        On the initial call to this subroutine to solve a
c        new problem, the user must set MODE = 0, and must set
c           X = A
c           XORF = B
c           TOL = desired tolerance
c        A and B denote endpoints defining a closed
c        interval in which a local minimum is to be found.
c        Permit A < B or A > B or A = B.  If A = B the solution, x = A,
c        will be found immediately.
c     TOL [in]  is an absolute tolerance on the uncertainty in the final
c        estimate of the minimizing abcissa, X1.  Recommend TOL .ge. 0.
c        This subr will set TOLI = TOL if TOL > 0., and otherwise sets
c        TOLI = 0..  The operational tolerance at any trial abcissa, X,
c        will be
c            TOL2 = 2 * abs(X) * sqrt(D1MACH(4)) + (2/3) * C3
c        where C3 = D1MACH(4) ** 4.
c
c        On each return this subr will set MODE to a value in the range
c        [1:3] to indicate the action needed from the calling program or
c        the status on termination.
c
c        = 1 means the calling program must evaluate
c        f(X), store the value in XORF, and then call this
c        subr again.
c
c        = 2 means normal termination.  X contains the point of the
c        minimum and XORF contains f(X).
c
c        = 3 as for 2 but the requested accuracy could not be obtained.
c
c        = 4 means termination on an error condition:
c            MODE on entry not in [0:1].
c     ------------------------------------------------------------------
c     This algorithm is due to Richard. P. Brent.  It is presented as
c     an ALGOL procedure, LOCALMIN, in his 1973 book,
c     "Algorithms for minimization without derivatives", Prentice-Hall.
c     Published as subroutine FMIN in Forsythe, Malcolm, and Moler,
c     "Computer Methods for Mathematical Computations", Prentice-Hall,
c     1977.
c     The current subroutine adapted from F.,M.& M. by C. L. Lawson, and
c     F. T. Krogh, JPL, Oct 1987, for use in Fortran 77 in the JPL
c     MATH77 library.  The changes improve performance when a minimum is
c     at an endpoint, and change the user interface to reverse
c     communication.  Unlike the original, the changed algorithm may
c     evaluate the function at one of the end points.
c     ------------------------------------------------------------------
c     Subprograms referenced: R1MACH, IERM1
c     ------------------------------------------------------------------
C     THE METHOD USED IS A COMBINATION OF  GOLDEN  SECTION  SEARCH  AND
C  SUCCESSIVE PARABOLIC INTERPOLATION.  CONVERGENCE IS NEVER MUCH SLOWER
C  THAN  THAT  FOR  A  FIBONACCI SEARCH.  IF  F  HAS A CONTINUOUS SECOND
C  DERIVATIVE WHICH IS POSITIVE AT THE MINIMUM (WHICH IS NOT  AT  AX  OR
C  BX),  THEN  CONVERGENCE  IS  SUPERLINEAR, AND USUALLY OF THE ORDER OF
C  ABOUT  1.324....
C      THE FUNCTION  F  IS NEVER EVALUATED AT TWO POINTS CLOSER TOGETHER
C  THAN  EPS*ABS(FMIN) + (TOLI/3), WHERE EPS IS APPROXIMATELY THE SQUARE
C  ROOT  OF  THE  RELATIVE  MACHINE  PRECISION.  IF F IS A UNIMODAL
C  FUNCTION AND THE COMPUTED VALUES OF F ARE ALWAYS UNIMODAL WHEN
C  SEPARATED BY AT LEAST  EPS*ABS(X) + (TOLI/3), THEN  FMIN APPROXIMATES
C  THE ABCISSA OF THE GLOBAL MINIMUM OF F ON THE INTERVAL [AX,BX] WITH
C  AN ERROR LESS THAN  3*EPS*ABS(FMIN) + TOLI.  IF F IS NOT UNIMODAL,
C  THEN FMIN MAY APPROXIMATE A LOCAL, BUT PERHAPS NON-GLOBAL, MINIMUM TO
C  THE SAME ACCURACY. (Comments from F., M. & M.)
c     ------------------------------------------------------------------
      real X,XORF,TOL, R1MACH
      real A,B,C,D,E,EPS,XM,P,Q,R,TOL1,TOL2,U,V,W
      real FU,FV,FW,FX,XI, TOLI, C3
      real ASAVE, BSAVE
      real HALF, ONE, TWO, THREE, FIVE, ZERO, PT95
      integer MODE, NEXT, IC, NP
      parameter(HALF = 0.5E0, ONE = 1.0E0, TWO = 2.0E0, THREE = 3.0E0)
      parameter(FIVE = 5.0E0, ZERO = 0.0E0, PT95 = 0.95E0)
      save
      data EPS / ZERO /, NP / 0 /
c     ------------------------------------------------------------------
      IF (MODE .LT. 0) THEN
         NP = -1 - MODE
         RETURN
      END IF
      IF (NP .GT. 0) THEN
         NP = NP - 1
         PRINT 10, X, XORF, IC
   10    FORMAT(' In SFMIN -- X=',1PE15.8,' XORF=',E15.8, ' IC=', I3)
      END IF
      if(MODE .eq. 1) then
         go to (301,302), NEXT
      endif
      if( MODE .ne. 0) then
c                                           Error:  MODE not 0 or 1
         call IERM1('SFMIN',MODE,0,
     *   'MODE not in the range [0:1]','MODE',MODE,'.')
         MODE = 4
         return
      endif
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C               C is the squared inverse of the "Golden Ratio", 1.618...
C               C = 0.381966
c
      C = HALF*(THREE - sqrt(FIVE))
C
C           EPS IS THE SQUARE ROOT OF THE RELATIVE MACHINE  PRECISION.
C
      if (EPS .eq. ZERO) EPS = sqrt(R1MACH(4))
C
C  INITIALIZATION
C
      A = X
      B = XORF
      if( A .gt. B) then
         XI = A
         A = B
         B = XI
      endif
c                           Now we have A .le. B
      ASAVE = A
      BSAVE = B
      TOLI = TOL
      if(TOLI .le. ZERO) TOLI = ZERO
      C3 = TOLI/THREE + EPS ** 4
      if (C3 .eq. ZERO) C3 = EPS ** 2
      V = A + C*(B - A)
      W = V
      XI = V
      IC = 0
      E = ZERO
c                   Return to calling prog for computation of FX = f(XI)
      X = XI
      MODE = 1
      NEXT = 1
      return
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  301 continue
      NEXT = 2
      FX = XORF
      FV = FX
      FW = FX
C                                         MAIN LOOP STARTS HERE
   20 XM = HALF*(A + B)
      TOL1 = EPS*abs(XI) + C3
      TOL2 = TWO*TOL1
C                            CHECK STOPPING CRITERION
C
      if (abs(XI - XM) .le. (TOL2 - HALF*(B - A))) go to 90
C
C                            IS GOLDEN-SECTION NECESSARY ?
C
      if (abs(E) .le. TOL1) go to 40
C
C                            FIT PARABOLA
C
      R = (XI - W)*(FX - FV)
      Q = (XI - V)*(FX - FW)
      P = (XI - V)*Q - (XI - W)*R
      Q = TWO*(Q - R)
      if (Q .gt. ZERO) P = -P
      Q =  abs(Q)
      R = E
      E = D
C
C  IS PARABOLA ACCEPTABLE
C
      if (abs(P) .ge. abs(HALF*Q*R)) go to 40
      if (P .le. Q*(A - XI)) go to 40
      if (P .ge. Q*(B - XI)) go to 40
C
C  A PARABOLIC INTERPOLATION STEP
C
      D = P/Q
      U = XI + D
C
C  F MUST NOT BE EVALUATED TOO CLOSE TO A OR B
C
      if ((U - A) .lt. TOL2 .or. (B - U) .lt. TOL2)
     *    D = sign(TOL1, XM - XI)
      go to 50
C
C  A GOLDEN-SECTION STEP
C
   40 IC = IC + 1
      if (IC .gt. 3) then
         if (A .eq. ASAVE) then
            if (IC .eq. 4) then
               E = A + (EPS * abs(A) + C3) - XI
            else
               if (B .ne. W) go to 45
               E = A - XI
            end if
         else if (B .eq. BSAVE) then
            if (IC .eq. 4) then
               E = B - (EPS * abs(B) + C3) - XI
            else
               if (A .ne. W) go to 45
               E = B - XI
            end if
         else
            IC = -99
            go to 45
         end if
         D = E
         go to 50
      end if
   45 if (XI .ge. XM) then
         E = A - XI
      else
         E = B - XI
      endif
      D = C*E
C
C  F MUST NOT BE EVALUATED TOO CLOSE TO XI
C
   50 if (abs(D) .ge. TOL1) then
         U = XI + D
      else
         if (B - A .gt. THREE * TOL1) TOL1 = TOL2
         U = min(B, max(A, XI + sign(PT95*TOL1, D)))
      endif
c
c                    Return to calling prog for computation of FU = f(U)
c                    Returning with MODE = 1 and NEXT = 2
      X = U
      return
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  302 continue
      FU = XORF
C
C  UPDATE  A, B, V, W, AND XI
C
      if (FU .le. FX) then
         if (FU .eq. FX) then
            if (B - min(XI, U) .gt. max(XI, U) - A) then
               B = max(XI, U)
            else
               A = min(XI, U)
            end if
         else
            if (U .ge. XI) then
               A = XI
            else
               B = XI
            endif
         end if
         V = W
         FV = FW
         W = XI
         FW = FX
         XI = U
         FX = FU
         go to 20
      endif
c
      if (U .lt. XI) then
         A = U
      else
         B = U
      endif
      if (FU .le. FW  .or.  W .eq. XI) then
         V = W
         FV = FW
         W = U
         FW = FU
      elseif (FU .le. FV  .or.  V .eq. XI  .or.  V .eq. W) then
         V = U
         FV = FU
      endif
      go to 20
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   90 continue
         X = XI
         XORF = FX
         MODE = 2
         if ((B - A .gt. THREE*TOLI) .and. (TOLI .ne. ZERO)) MODE = 3
         return
         end

