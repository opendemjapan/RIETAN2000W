      SUBROUTINE MINI(YFIT,B,XINT,DEG,NTERMS,NSTEP,X)
C     POWELL'S ALGORITHM FOR UNCONSTRAINED OPTIMIZATION WITHOUT THE USE
C     OF DERIVATIVES (CONJUGATE DIRECTION METHOD)
C     M. J. D. POWELL, COMPUTER J., 7 (1964) 155.
      PARAMETER (NB=7000,NT=999,NR=400,NPH=8,NAP=150)
      DIMENSION YFIT(*),W(NR),SECND(NR),B(*),XINT(*),DEG(*),X(*)
      COMMON /NUML/ NUPDT
      COMMON /ONE/ Y(NR),S(NR),FX,FY
      COMMON /FL0/ NC,TK,CONV,NCONV,FINC,SUMPEN,SUMP1,FACTOR
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /OLD/ AOLD(NT),AOUT(NT)
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /POW/ MITER,STEP,ACC
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /TWO/ DIRECT(NR,NR),BEFORE(NR),FIRST(NR)
      COMMON /THREE/ N,NFUNCT,NDRV,ITER,INDIC,NPRINT
      EQUIVALENCE (W,SECND)
      
      IF (NMODE .EQ. 4 .OR. NMODE .EQ. 5) CALL SCYPEAK(0)
      ICONVG=1
      ITER=0
      NTRY=1
      N1=N-1
      DO 2 I=1,N
      DO 1 J=1,N
      DIRECT(J,I)=0.0
    1 CONTINUE
      IF(FIRST(I).EQ.0.0) THEN
         IF(X(I).EQ.0.0) THEN
            DIRECT(I,I)=STEP
         ELSE
            DIRECT(I,I)=ABS(X(I))*STEP
         ENDIF
      ELSE
         IF(X(I).EQ.0.0) THEN
            DIRECT(I,I)=FIRST(I)
         ELSE
            DIRECT(I,I)=ABS(X(I))*FIRST(I)
         ENDIF
      ENDIF
    2 CONTINUE
      STEP=1.0
      STEPA=STEP
  100 FX=AOFFN(DEG,XINT,YFIT,NPTS,NSTEP,B,0,IDPOS)
      NFUNCT=NFUNCT+1
      WRITE(6,'(/11X,A)') 'Calculation using initial parameters'
      IF(NC.GT.0) THEN
         WRITE(6,2010) FX,FX-SUMPEN,SUMPEN
 2010    FORMAT(' ',10X,'AOF =',1P,G13.6,4X,'OF =',1P,G13.6,4X,'PF =',
     &   1P,G13.6)
      ELSE
         WRITE(6,'(11X,A,1P,G13.6)') 'OF =',FX
      ENDIF
      DO 85 J=1,NTERMS
         AOLD(J)=B(J)
   85 CONTINUE
      GO TO 301
    3 ITER=ITER+1
      CALL PRPARA(YFIT,B,XINT,DEG,NTERMS,NPTS,NSTEP,X)
      IF(ITER.EQ.MITER) RETURN
  301 F1=FX
      DO 4 I=1,N
    4 BEFORE(I)=X(I)
      SUM=0.0
      DO 9 I=1,N
      DO 5 J=1,N
    5 S(J)=DIRECT(J,I)*STEP
      CALL SEARCH(YFIT,B,XINT,DEG,NTERMS,NSTEP,X)
      A=FX-FY
      IF(A-SUM) 7,7,6
    6 ISAVE=I
      SUM=A
    7 DO 8 J=1,N
    8 X(J)=Y(J)
    9 FX=FY
      F2=FX
      DO 10 I=1,N
   10 W(I)=2.0*X(I)-BEFORE(I)
      CALL FUN(YFIT,W,F3,B,XINT,DEG,NTERMS,NSTEP,*999)
  999 A=F3-F1
      IF(A) 11,19,19
   11 A=2.0*(F1-2.0*F2+F3)*((F1-F2-SUM)/A)**2
      IF(A-SUM) 12,19,19
   12 IF(ISAVE-N) 13,15,15
   13 DO 14 I=ISAVE,N1
      II=I+1
      DO 14 J=1,N
   14 DIRECT(J,I)=DIRECT(J,II)
   15 A=0.0
      DO 16 J=1,N
      DIRECT(J,N)=X(J)-BEFORE(J)
   16 A=DIRECT(J,N)**2+A
      A=1.0/SQRT(A)
      DO 17 J=1,N
      DIRECT(J,N)=DIRECT(J,N)*A
   17 S(J)=DIRECT(J,N)*STEP
      CALL SEARCH(YFIT,B,XINT,DEG,NTERMS,NSTEP,X)
      FX=FY
      DO 18 I=1,N
   18 X(I)=Y(I)
   19 CALL TEST(F1,FX,BEFORE,X,FLAG,N,ACC)
      IF(FLAG) 22,22,20
   20 IF(F1-FX) 121,120,120
  121 STEP=-0.4*SQRT(ABS(F1-FX))
      GO TO 123
  120 STEP=0.4*SQRT(F1-FX)
  123 IF(STEPA-STEP) 21,3,3
   21 STEP=STEPA
      GO TO 3
   22 GO TO (23,24), ICONVG
   23 RETURN
   24 GO TO (25,27), NTRY
   25 NTRY=2
      DO 26 I=1,N
      FIRST(I)=X(I)
   26 X(I)=X(I)+ACC*10.0
      FFIRST=FX
      GO TO 100
   27 FSECND=FX
      A=0.0
      DO 28 I=1,N
      SECND(I)=X(I)
      S(I)=FIRST(I)-SECND(I)
   28 A=A+S(I)**2
      IF(A) 23,23,29
   29 A=STEP/SQRT(A)
      DO 30 I=1,N
   30 S(I)=S(I)*A
      CALL SEARCH(YFIT,B,XINT,DEG,NTERMS,NSTEP,X)
      CALL TEST(FFIRST,FY,FIRST,Y,FLAG,N,ACC)
      IF(FLAG) 32,32,31
   31 CALL TEST(FSECND,FY,SECND,Y,FLAG,N,ACC)
      IF(FLAG) 32,32,34
   32 DO 33 I=1,N
   33 X(I)=Y(I)
      FX=FY
      RETURN
   34 A=A/STEP
      DO 35 I=1,N
      DIRECT(I,1)=(FIRST(I)-SECND(I))*A
   35 FIRST(I)=SECND(I)
      GO TO 3
      END

************************************************************************

      SUBROUTINE SEARCH(YFIT,B,XINT,DEG,NTERMS,NSTEP,X)
C     DSC-POWELL UNIDIMENSIONAL SEARCH
      PARAMETER (NR=400)
      DIMENSION YFIT(*),B(*),XINT(*),DEG(*),X(*)
      COMMON /FL0/ NC,TK,CONV,NCONV,FINC,SUMPEN,SUMP1,FACTOR
      COMMON /ONE/ Y(NR),S(NR),FX,FY
      COMMON /THREE/ N,NFUNCT,NDRV,ITER,INDIC,NPRINT
      COMMON /TWO/ H(NR,NR),DELG(NR),GX(NR)
      IEXIT=0
      NTOL=0
      FTOL=.001
      FTOL2=FTOL/100.0
      FA=FX
      FB=FX
      FC=FX
      DA=0.0
      DB=0.0
      DC=0.0
      K=-2
      M=0
      STEP=1.0
      D=STEP
    1 DO 2 I=1,N
    2 Y(I)=X(I)+D*S(I)
      CALL FUN(YFIT,Y,F,B,XINT,DEG,NTERMS,NSTEP,*999)
      K=K+1
      IF(F-FA) 5,3,6
    3 DO 4 I=1,N
    4 Y(I)=X(I)+DA*S(I)
      FY=FA
      IF(NPRINT.EQ.2) WRITE(6,2100)
 2100 FORMAT(/11X,'Search failed.  Function value independent of search
     &direction.')
      GO TO 326
    5 FC=FB
      FB=FA
      FA=F
      DC=DB
      DB=DA
      DA=D
      D=2.0*D+STEP
      GO TO 1
    6 IF(K) 7,8,9
    7 FB=F
      DB=D
      D=-D
      STEP=-STEP
      GO TO 1
    8 FC=FB
      FB=FA
      FA=F
      DC=DB
      DB=DA
      DA=D
      GO TO 21
    9 DC=DB
      DB=DA
      DA=D
      FC=FB
      FB=FA
      FA=F
   10 D=0.5*(DA+DB)
      DO 11 I=1,N
   11 Y(I)=X(I)+D*S(I)
      CALL FUN(YFIT,Y,F,B,XINT,DEG,NTERMS,NSTEP,*999)
   12 IF((DC-D)*(D-DB)) 15,13,18
   13 DO 14 I=1,N
   14 Y(I)=X(I)+DB*S(I)
      FY=FB
      IF(IEXIT.EQ.1) GO TO 32
      IF(NPRINT.EQ.2) WRITE(6,2200)
 2200 FORMAT(/11X,'Search failed.  Location of minimum limited by roundi
     &ng.')
      GO TO 325
   15 IF(F-FB) 16,13,17
   16 FC=FB
      FB=F
      DC=DB
      DB=D
      GO TO 21
   17 FA=F
      DA=D
      GO TO 21
   18 IF(F-FB) 19,13,20
   19 FA=FB
      FB=F
      DA=DB
      DB=D
      GO TO 21
   20 FC=F
      DC=D
   21 A=FA*(DB-DC)+FB*(DC-DA)+FC*(DA-DB)
      IF(A) 22,30,22
   22 D=0.5*((DB*DB-DC*DC)*FA+(DC*DC-DA*DA)*FB+(DA*DA-DB*DB)*FC)/A
      IF((DA-D)*(D-DC)) 13,13,23
   23 DO 24 I=1,N
   24 Y(I)=X(I)+D*S(I)
      CALL FUN(YFIT,Y,F,B,XINT,DEG,NTERMS,NSTEP,*999)
      IF(ABS(FB)-FTOL2) 25,25,26
   25 A=1.0
      GO TO 27
   26 A=1.0/FB
   27 IF((ABS(FB-F)*A)-FTOL) 28,28,12
   28 IEXIT=1
      IF(F-FB) 29,13,13
   29 FY=F
      GO TO 32
   30 IF(M) 31,31,13
   31 M=M+1
      GO TO 10
   32 DO 99 I=1,N
      IF(Y(I).NE.X(I)) GO TO 325
   99 CONTINUE
      GO TO 33
  325 IF(NTOL.NE.0 .AND. NPRINT.EQ.2) WRITE(6,3000) NTOL
 3000 FORMAT(/11X,'Tolerance reduced ',I1,' time(s)')
  326 IF(FY.LT.FX) RETURN
      IF(S(1).NE.-GX(1) .OR. FY.LT.FX) RETURN
      IF(NPRINT.EQ.2) THEN
         WRITE(6,5100) ITER,NFUNCT,NDRV
         IF(NC.GT.0) THEN
            WRITE(6,5200) FY,FY-SUMPEN,SUMPEN
         ELSE
            WRITE(6,'(11X,A,1P,G13.6/)') 'OF =',FY
         ENDIF
      ENDIF
 5100 FORMAT(/11X,'Search failed on a gradient step.  Search terminated.
     &'//11X,'ITER   =',I3,5X,'NFUNCT =',I4,5X,'NDRV   =',I3)
 5200 FORMAT(11X,'AOF =',1P,G13.6,4X,'OF =',G13.6,4X,'PF =',G13.6/)
      RETURN
   33 IF(NTOL.EQ.5) GO TO 34
      IEXIT=0
      NTOL=NTOL+1
      FTOL=FTOL/10.0
      GO TO 12
   34 IF(NPRINT.EQ.2) WRITE(6,2000)
 2000 FORMAT(/11X,'A point better than the entering point cannot be foun
     &d')
      RETURN
  999 DO 91 I=1,N
         Y(I)=X(I)
   91 CONTINUE
      FY=FX
      END

************************************************************************

      SUBROUTINE TEST(FI,FF,RI,RF,FLAG,N,ACC)
      DIMENSION RI(*),RF(*)
      FLAG=2.0
      IF(ABS(FI)-ACC) 2,2,1
    1 IF(ABS((FI-FF)/FI)-ACC) 3,3,7
    2 IF(ABS(FI-FF)-ACC) 3,3,7
    3 DO 6 I=1,N
      IF(ABS(RI(I))-ACC) 5,5,4
    4 IF(ABS((RI(I)-RF(I))/RI(I))-ACC) 6,6,7
    5 IF(ABS(RI(I)-RF(I))-ACC) 6,6,7
    6 CONTINUE
      FLAG=-2.0
    7 END

************************************************************************

      SUBROUTINE FUN(YFIT,Z,FX,A,XINT,DEG,NTERMS,NSTEP,*)
C     FUNCTION TO BE MINIMIZED BY POWELL'S METHOD
      PARAMETER (NT=999)
      DIMENSION Z(*),A(*),XINT(*),DEG(*),YFIT(*)
      COMMON /OLD/ AOLD(NT),AOUT(NT)
      COMMON /THREE/ N,NFUNCT,NDRV,ITER,INDIC,NPRINT
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      
C     COPY CURRENT VALUES OF VARIABLE PARAMETERS INTO ARRAY A
      IDSAME=1
      DO J = 1, NTERMS
         IF(ID(J).EQ.1) THEN
            A(J)=Z(IR(J))
            IF(A(J).NE.AOLD(J)) IDSAME=0
         ENDIF
      END DO
      IF(IDSAME.EQ.1) RETURN 1
      CALL REFPOW(A,NTERMS)
C
      FX = AOFFN(DEG,XINT,YFIT,NPTS,NSTEP,A,1,IDPOS)
      NFUNCT=NFUNCT+1
      END

************************************************************************

      SUBROUTINE PRPARA(YFIT,A,XINT,DEG,NTERMS,NPTS,NSTEP,X)
C     PRINT THE CURRENT VALUES OF ARRAYS A AND ID
      CHARACTER PARNAM*60,PHNAME*25
      PARAMETER (NT=999,NR=400,NPH=8,NAP=150)
C     DIMENSION YFIT(*),A(*),XINT(*),DEG(*),X(*),IDLAT(10)
      DIMENSION YFIT(*),A(*),XINT(*),DEG(*),X(*)
      COMMON /CC/ APR(NT)
      COMMON /FL0/ NC,TK,CONV,NCONV,FINC,SUMPEN,SUMP1,FACTOR
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /OLD/ AOLD(NT),AOUT(NT)
      COMMON /ONE/ Y(NR),S(NR),FX,FY
      COMMON /PAR/ PARNAM(NT),PHNAME(NPH)
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /THREE/ N,NFUNCT,NDRV,ITER,INDIC,NPRINT
      COMMON /UCC/ IDPEAK(NPH),IDSF(NPH),IDPRF(NPH),IDPO(NPH)
      PARAMETER (NB=7000)
      COMMON /SAVEYP/ SYPEAK(NB),SAVEA(NT),NTOTPAR
      
      DO J = 1, NTERMS
         IF(ID(J).EQ.1) A(J)=X(IR(J))
      END DO
      DO L = 1, NPHASE
         SELECT CASE (NMODE)
            CASE(0,5)
               IDPEAK(L)=1
               IDSF(L)=1
               IDPRF(L)=1
            CASE(2:4)
               IDPEAK(L)=1
               IDPRF(L)=1
            CASE(6)
               IDPRF(L)=1
         END SELECT
      END DO
      IF (NMODE .EQ. 4 .OR. NMODE .EQ. 5) THEN
*        RESOTRE THE PARAMETER AND YPEAK GIVING THE MINIMUM 
*        SUM OF SQUARES
*        RESTORE A
         DO J = 1, NTOTPAR
            A(J) = SAVEA(J)
         END DO
      END IF
      FX = AOFFN(DEG,XINT,YFIT,NPTS,NSTEP,A,1,IDPOS)
      IF (NMODE .EQ. 4 .OR. NMODE .EQ. 5) CALL SCYPEAK(0)
 
      ENTRY PRFIN(A,NTERMS)
      WRITE(6,100) ITER,NFUNCT
  100 FORMAT(//' ',10X,'ITER   =',I3/' ',10X,'NFUNCT =',I4)
      IF(NC.GT.0) THEN
         WRITE(6,110) FX,FX-SUMPEN,SUMPEN
  110    FORMAT(' ',10X,'AOF    =',1P,G13.6,4X,'OF    =',1P,G13.6,4X,
     &   'PF    =',1P,G13.6)
      ELSE
         WRITE(6,'(11X,A,1P,G13.6)') 'OF     =',FX
      ENDIF
      CALL SETPAR(A)
      CALL RFACTR(A,DEG,YFIT,XINT,NPTS,0)
      IF(NPRINT.NE.0) THEN
         WRITE(6,'(/12X,A)') 'No.    A(old)    +   Delta.A   =  '//
     &   'A(refined)    ID'
         CALL APRNTD(A,NTERMS)
         DO JJ = 1, NTERMS
            IF(ID(JJ) .NE. 0 .OR. LATVAR(JJ) .EQ. 1) THEN
               DELTA=APR(JJ)-AOUT(JJ)
            ELSE
               DELTA=0.0
            ENDIF
            WRITE(6,'(11X,I3,2X,1P,3(G13.6,1X),I4,3X,A)')
     &      JJ,AOUT(JJ),DELTA,APR(JJ),ID(JJ),PARNAM(JJ)
            AOUT(JJ)=APR(JJ)
         END DO
      ENDIF
      END

************************************************************************

      SUBROUTINE SETPAR(A)
C     UPDATE THE METRIC TENSOR AND CONSTRAINED PARAMETERS
      DIMENSION A(*)
      PARAMETER (NT=999,NCS=400,NPH=8,NAP=150)
      COMMON /CNSTR/ C1(NCS),C2(10,NCS),IN(10,NCS),NCNPAR(NCS),
     &  NACNS(NCS),NCNSTR
      COMMON /H/ NQ(NPH),NCENTR(NPH),LAUEG(NPH)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      LP(J)=KPHB(L)+J
      P(J)=A(LP(J))

      K=1
*     CONSTRAINTS FOR PPP'S
      K = 1
      DO J = NPPPX, NPPPX + NRELAX*NPRIM - 1
         IF (ID(J) .EQ. 2) THEN
            A(J) = A(IN(1,K))
            K = K + 1
         END IF
      END DO
      
      DO L = 1, NPHASE
         IF(LAUEG(L) .EQ. 6 .OR. LAUEG(L) .EQ. 7) THEN
            A(LP(NBX))=P(NAX)
         ELSE IF(LAUEG(L) .EQ. 8 .OR. LAUEG(L) .EQ. 10) THEN
            A(LP(NBX))=P(NAX)
            A(LP(NCX))=P(NAX)
            A(LP(NBETX))=P(NALPX)
            A(LP(NGAMX))=P(NALPX)
         ELSE IF(LAUEG(L) .EQ. 9 .OR. LAUEG(L) .EQ. 11 .OR. 
     &   LAUEG(L) .EQ. 12 .OR. LAUEG(L) .EQ. 13) THEN
            A(LP(NBX))=P(NAX)
            A(LP(NGAMX))=0.5*P(NAX)
         ELSE IF(LAUEG(L).GE.14) THEN
            A(LP(NBX))=P(NAX)
            A(LP(NCX))=P(NAX)
         END IF
         DO J = LP(1), KPHE(L)
            IF (ID(J) .EQ. 2) THEN
               A(J) = C1(K)
               DO MM = 1, NCNPAR(K)
                  A(J) = A(J) + C2(MM,K)*A(IN(MM,K))
               END DO
               K = K+1
            END IF
         END DO
      END DO
      END

************************************************************************

      SUBROUTINE OVERWR(A,NTERMS,NSPLBL,LABEL,LPAR,NLINE)
*     UPDATE ARRAY A IN THE INPUT FILE
      PARAMETER (NT=999)
      INTEGER LPAR(*)
      INTEGER NSPLBL(*)
      REAL A(*)
      CHARACTER LABEL(*)*25,FILE5*50
      CHARACTER*80 LINE
      COMMON /CC/ APR(NT)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN

*     READ THE INPUT FILE
      INP = 5
      REWIND INP
*     SCRATCH FILE STORING *.ins FOR USE DURING UPDATING PARAMETERS
      ICOPY = 55
      OPEN (ICOPY,STATUS='SCRATCH',ACCESS='SEQUENTIAL')
      DO WHILE (.TRUE.)
         READ(INP,'(A)',END=1) LINE
         CALL PRLINE(LINE,ICOPY)
      END DO

*     REOPEN FILE INP IN THE CASE OF ONLY LS FORTRAN
    1 CALL OPENFILE(1,FILE5)
      REWIND INP
      REWIND ICOPY
      LFIRST = 0
      IA = 0
      DO WHILE (.TRUE.)
         READ(ICOPY,'(A)') LINE
         CALL CHKLBL(LINE,LFIRST,ID2,IA)
         IF (LFIRST .EQ. 0) THEN
            CALL PRLINE(LINE,INP)
         ELSE
            BACKSPACE ICOPY
            EXIT
         END IF
      END DO

      CALL APRNTD(A,NTERMS)
*     WRITE LABEL, A, AND ID
      CALL NEWPAR(NLINE,LPAR,NSPLBL,LABEL,APR,IDSAVE,INP,ICOPY)
 
      DO WHILE (.TRUE.)
         READ(ICOPY,'(A)',END=9) LINE
         CALL PRLINE(LINE,INP)
      END DO
    9 CLOSE(UNIT=ICOPY)
      END

************************************************************************

      SUBROUTINE PRLINE(LINE,INP)
*     DETERMINE THE LINE LENGTH AND PRINT OUT LINE
      CHARACTER*80 LINE
      DO J = 80, 1, -1
         IF (LINE(J:J) .NE. ' ') EXIT
      END DO
      IF (J .EQ. 0) THEN
         WRITE(INP,'()')
      ELSE
         WRITE(INP,'(A)') LINE(1:J)
      END IF
      END

************************************************************************

      SUBROUTINE RFACTR(A,DEG,YFIT,XINT,NPTS,LRB)
C     CALCULATE AND PRINT OUT RELIABILITY INDICES
      PARAMETER (NT=999,NB=7000,NP=80000,NPH=8,NAP=150)
      REAL A(*),DEG(*),YFIT(*),XINT(*)
      DOUBLE PRECISION SUM(9),SUMPH(4,NPH),YOPROF,YCPROF,YCIND
      DIMENSION ESCIO(NPH),NUMESC(NPH)
      CHARACTER PARNAM*60,PHNAME*25,TOP1*8,FILE20*50
      EXTERNAL PRFINT
      COMMON /A/ IH(NB),IK(NB),IL(NB),L12(NB),NOPH(NB)
      COMMON /ABSVR/ ABSORP(NP),VARLEN(NP),ROUGH(NP),NSURFR
      COMMON /B/ U(NB),NREFL
      COMMON /C/ BG(NP),BGINC(NP),DEGNOR(NP)
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /OBSCAL/ OBSI(NB),CALCI(NB),LREFINT(NB)
      COMMON /DETECT/ VARINV(NP)
      COMMON /EXCREG/ DEGEXC(2,50),NSKIP
      COMMON /G/ I
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /ORDER/ IPOINT(NB)
      COMMON /I/ RDEG(2,NB)
      INTEGER WIDTH,HEIGHT,YMIN,YMAX,LBGPLOT,LDEL,OFFSETD,MSIZE,FSIZE,
     &  LSIZE,OFFSET,INDREF
      COMMON /IGOR/ WIDTH,HEIGHT,YMIN,YMAX,LBGPLOT,LDEL,OFFSETD,MSIZE,
     &  FSIZE,LSIZE,OFFSET(NPH),INDREF
      COMMON /INTFMT/ TOP1
      COMMON /MC/ CTHM,NTRAN,PCOR
      COMMON /PAR/ PARNAM(NT),PHNAME(NPH)
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /RLS/ NRANGE
      COMMON /T/ FF(NB),COEF(NPH)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
      COMMON /UNCONS/ NUCREF,IPHASE(90),IHKL(3,90),NOPPP(90),LPPP(NB),
     &  RCUT(2,90),NCUT
      COMMON /YOKO/ ABSCISSA(NP),NPOINTS,AWRY
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      DATA RD/0.01745329/
      
      DO J = 1, 9
         SUM(J) = 0.0D0
      END DO
*     XINT: RAW INTENSITY
*     YFIT: INTENSITY CALCULATED WITH THE MODEL FUNCTION
*     VARINV: 1/XINT (SINGLE COUNTER) AND N/XINT (N COUNTERS)
*     BG: BACKGROUND
      DO KK = 1, NPTS
         DELTA = XINT(KK) - YFIT(KK)
         IF (VARINV(1) .EQ. 0.0) THEN
            SUM(1) = SUM(1) + DELTA**2/XINT(KK)
            SUM(2) = SUM(2) + XINT(KK)
         ELSE
            SUM(1) = SUM(1) + DELTA**2*VARINV(KK)
            SUM(2) = SUM(2) + VARINV(KK)*XINT(KK)**2
            SUM(9) = SUM(9) + XINT(KK)
         END IF
         SUM(3) = SUM(3) + ABS(DELTA)
         SUM(4) = SUM(4) + ABS(XINT(KK) - BG(KK))
         IF (KK .EQ. 1) CYCLE
         IF (VARINV(1) .EQ. 0.0) THEN
            SUM(8) = SUM(8) + DELTA**2/XINT(KK)
            SUM(5) = SUM(5) + (DELTA/SQRT(XINT(KK)) -
     &      (XINT(KK-1) - YFIT(KK-1))/SQRT(XINT(KK-1)))**2
         ELSE
            SUM(8) = SUM(8) + DELTA**2*VARINV(KK)
            SUM(5) = SUM(5) + (DELTA*SQRT(VARINV(KK)) -
     &      (XINT(KK-1) - YFIT(KK-1))*SQRT(VARINV(KK-1)))**2
         END IF
         SUM(6) = SUM(6) + (DELTA - (XINT(KK-1) - YFIT(KK-1)))**2
         SUM(7) = SUM(7) + DELTA**2
      END DO
      RWP = SQRT(SUM(1)/SUM(2))
      IF (VARINV(1) .EQ. 0.0) THEN 
         RP = SUM(3)/SUM(2)
      ELSE
         RP = SUM(3)/SUM(9)
      END IF
      RR = SUM(3)/SUM(4)
      RE = SQRT((NPTS - NRFN)/SUM(2))
*     Durbin-Watson statistic, d, according to Hill & Flack (1987)
      DWDS1 = SUM(5)/SUM(8)
*     Durbin-Watson statistic, d, according to 
*     Schwarzenbach et al., Acta Crystallogr., Sect. A, 45 (1989) 63.
      DWDS2 = SUM(6)/SUM(7)

      IF (LRB .NE. 0) WRITE(6,100)
  100 FORMAT(//' ',10X,'Reliability factors, goodness-of-fit indicator, 
     &and Durbin-Watson statistic'//)
      WRITE(6,110) RWP, RP, RR, RE, RWP/RE, DWDS1, DWDS2
  110 FORMAT(' ',10X,'Rwp =',2P,F6.2,5X,'Rp =',F6.2,5X,'RR =',F6.2,
     &5X,'Re =',F6.2,5X,'S =',0P,F7.4,5X,'d1 =',F7.4,5X,'d2 =',F7.4)
      IF (LRB .EQ. 0) RETURN

      IF (NPAT .EQ. 5) THEN
         CALL OPENFILE(4,FILE20)
         REWIND 20
         WRITE(20,'(A)') 'IGOR'
         IF (INDREF .EQ. 1) 
     &   WRITE(20,'(A/A)') 'WAVES/O xrefl, yrefl','BEGIN'
      END IF

*     COPY NPTS INTO NPOINTS AND DEG INTO VARINV, WHICH WILL BE USED 
*     IN FUNCTION PRFINT
      NPOINTS = NPTS
      DO J = 1, NPTS
         ABSCISSA(J) = DEG(J)
      END DO

*     LREFINT(I) = 0: FULL OBSERVED (POSITIVE) AND CALCULATED PROFILES 
*                     ARE OBTAINED.
*     LREFINT(I) = 1: PART OF THE OBSERVED PROFILE IS LACKING AS A 
*                     RESULT OF EXCLUDING A 2-THETA REGION.
*     LREFINT(I) = 2: FULL OBSERVED (NEGATIVE) AND CALCULATED PROFILES
*                     ARE OBTAINED.
*     LREFINT(I) = 3: PART OF THE PROFILE EXCEEDS 2-THETA(MAX).
      DO I = 1, NREFL
         LREFINT(I) = 0
         OBSI(I) = 0.0
*        NUMERICAL INTEGRATION OF THE BRAGG INTENSITY FROM RDEG(1,I) 
*        TO RDEG(2,I)
*        RDEG: 2-THETA REGION FOR CALCULATING THE PROFILE (RADIAN)
         CALL QROMB(PRFINT, RDEG(1,I), RDEG(2,I), CALCI(I), A)
      END DO
      
*     DETERMINE THE EXCLUDED REGION WITH THE HIGHEST 2-THETA
*     NSKIP: NUMBER OF EXCLUDED REGIONS
*     DEGEXC: EXCLUDED 2-THETA REGION (DEGREE)
      LAST = NSKIP
      DO KK = 1, NSKIP - 1
         IF (DEGEXC(1,KK) .GT. DEGEXC(1,NSKIP)) LAST = KK
      END DO

*     LEXCEND = 1: THE HIGH 2-THETA OF RAW INTENSITY DATA IS CUT.
*     THE PROGRAM ASSUMES THAT PART OF INTENSITY DATA HAVE BEEN EXCLUDED
*     WITH THE LAST DEGEXC.
      IF (DEG(NPTS) .LT. DEGEXC(1,LAST)*RD) THEN
         LEXCEND = 1
      ELSE
         LEXCEND = 0
      END IF
      
      ISTART = 1
      DO IPH = 1, NPHASE
         DO J = 1, 4
            SUMPH(J,IPH) = 0.0D0
         END DO

         ESCIO(IPH)= 0.0
         NUMESC(IPH)= 0
      END DO
      DO 20 J = 1, NREFL
         I = IPOINT(J)
         IF (LEXCEND .EQ. 1 .AND. RDEG(2,I) .GT. DEG(NPTS)) THEN
            DO JJ = J, NREFL
               LREFINT(IPOINT(JJ)) = 3
            END DO
            EXIT
         END IF
	 
         DO KK = 1, NSKIP
            IF ((RDEG(1,I) .GT. DEGEXC(1,KK)*RD .AND. 
     &      RDEG(1,I) .LT. DEGEXC(2,KK)*RD) .OR.
     &      (RDEG(2,I) .GT. DEGEXC(1,KK)*RD .AND. 
     &      RDEG(2,I) .LT. DEGEXC(2,KK)*RD) .OR.
     &      (RDEG(1,I) .LT. DEGEXC(1,KK)*RD .AND.
     &      RDEG(2,I) .GT. DEGEXC(2,KK)*RD)) THEN
               LREFINT(I) = 1
*              In the present version, the observed intensity is set 
*              equal to the calculated one in *.mem for a reflection 
*              whose profile is included in an excluded region.
               OBSI(I) = CALCI(I)
               GO TO 20
            END IF
         END DO

*        YOPROF: INTEGRATION OF THE OBSERVED OVERLAPPED PROFILE
*        YCPROF: INTEGRATION OF THE CALCULATED OVERLAPPED PROFILE
         YOPROF = 0.0D0
         YCPROF = 0.0D0
         YCIND = 0.0D0

*        SUBROUTINE HRPDINT IS CALLED WHEN THE 8TH CHARACTER IS '$'
*        FOR EVERY FORMAT
         IF (TOP1(1:4) .EQ. 'OVLP' .OR. TOP1(8:8) .EQ. '$') THEN
            CALL HRPDINT(DEG,XINT,YFIT,NPTS,ISTART,YOPROF,YCPROF,YCIND,
     &      A,1)
         ELSE
            CALL AREAVS(DEG,XINT,YFIT,NPTS,ISTART,YOPROF,YCPROF,YCIND,
     &      LFLAG,A,1)
*           LFLAG = 1: REACHED REFLECTIONS WHOSE PROFILES ARE PARTIALLY 
*                      LACKING
            IF (LFLAG .EQ. 1) EXIT
         END IF
         IF (YOPROF .GT. 0.0D0 .AND. YCPROF .GT. 0.0D0 .AND. 
     &   CALCI(I) .GT. 0.0D0) THEN
            OBSI(I) = YOPROF*CALCI(I)/YCPROF
            ESCIO(NOPH(I)) = ESCIO(NOPH(I)) + YCIND/OBSI(I)
            NUMESC(NOPH(I)) = NUMESC(NOPH(I)) + 1
            SUMPH(1,NOPH(I)) = SUMPH(1,NOPH(I)) + 
     &                         ABS(OBSI(I) - CALCI(I))
            SUMPH(2,NOPH(I)) = SUMPH(2,NOPH(I)) + OBSI(I)
*Rev 1.0z 2001.05.02 Izumi {
*           No RF can be obtained in Le Bail refinement or individual
*           profile fitting for the 1st phase because FF(I) is not calculated.
            IF (NMODE .LT. 4 .OR. NOPH(I) .GE. 2) THEN
*              Call PRFINT to obtain AWRY at the peak position
*              AWRY = AWR*YPEAK(I)
*              AWR = (absorption factor)*(irradiation width)*
*                    (surface roughness)
               PROFBA = PRFINT(PEAK(I),A)
*              CALCI(I) = AWRY = AWR*s*U(I)*FLP(I)*PROR1(I)*FF(I),
*              where s is the scale factor.
               FFOBS = OBSI(I)/(AWRY/FF(I))
               FFCALC = CALCI(I)/(AWRY/FF(I))
               IF (FFOBS .GT. 0.0 .AND. FFCALC .GT. 0.0) THEN
                  SUMPH(3,NOPH(I)) = SUMPH(3,NOPH(I)) +
     &                               ABS(SQRT(FFOBS) - SQRT(FFCALC))
                  SUMPH(4,NOPH(I)) = SUMPH(4,NOPH(I)) + SQRT(FFOBS)
               END IF
            END IF
* }
         ELSE IF (YOPROF .LE. 0.0D0) THEN
            LREFINT(I) = 2
         END IF
   20 CONTINUE

      IF (NPAT .EQ. 5 .AND. INDREF .EQ. 1) WRITE(20,'(A)') 'END'

*Rev 1.0z 2001.05.04 Izumi {
      DO IPH = 1, NPHASE
         IF (NMODE .LT. 4 .OR. IPH .GE. 2) THEN      
            WRITE(6,300) PHNAME(IPH), SUMPH(1,IPH)/SUMPH(2,IPH), 
     &      SUMPH(3,IPH)/SUMPH(4,IPH), ESCIO(IPH)/FLOAT(NUMESC(IPH))
  300       FORMAT(/11X,A/11X,'RI =',2P,F6.2,5X,'RF =',F6.2,5X,
     &      'E(SCIO) = ',1P,G13.6)
         ELSE
            WRITE(6,310) PHNAME(IPH), SUMPH(1,IPH)/SUMPH(2,IPH), 
     &      ESCIO(IPH)/FLOAT(NUMESC(IPH))
  310       FORMAT(/11X,A/11X,'RI =',2P,F6.2,5X,'E(SCIO) = ',1P,G13.6)
         END IF
      END DO
* }
      END
 
************************************************************************
 
      SUBROUTINE INTEGINT(A,DEG,YFIT,XINT,NPTS)
*     UPDATE INTEGRATED INTENSITIES, YPEAK, ON THE BASIS OF NEW PARAMETERS
*     USING RIETVELD'S PROCEDURE
      PARAMETER (NB=7000,NP=80000,NPH=8,NAP=150,NS=48)
      REAL A(*),DEG(*),YFIT(*),XINT(*)
      DOUBLE PRECISION YOPROF,YCPROF,YCIND
      CHARACTER TOP1*8
      EXTERNAL PRFINT
      COMMON /A/ IH(NB),IK(NB),IL(NB),L12(NB),NOPH(NB)
      COMMON /B/ U(NB),NREFL
      COMMON /OBSCAL/ OBSI(NB),CALCI(NB),LREFINT(NB)
      COMMON /DETECT/ VARINV(NP)
      COMMON /EXCREG/ DEGEXC(2,50),NSKIP
      COMMON /G/ I
      COMMON /ORDER/ IPOINT(NB)
      COMMON /I/ RDEG(2,NB)
      COMMON /INTFMT/ TOP1
      COMMON /MC/ CTHM,NTRAN,PCOR
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      INTEGER R,RX,RY,RZ
      COMMON /RT/ TV(3,NS,NPH),HT(NS,NB),NSYM(NPH),NSITE(NPH),
     &  IDSYM(NS,NAP,NPH),R(3,3,NS,NPH),RX(NS,NB),RY(NS,NB),RZ(NS,NB)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
      COMMON /YOKO/ ABSCISSA(NP),NPOINTS,AWRY
      DATA RD/0.01745329/
      
*     COPY NPTS INTO NPOINTS AND DEG INTO VARINV, WHICH WILL BE USED 
*     IN SUBROUTINE PRFINT
      NPOINTS = NPTS
      DO J = 1, NPTS
         ABSCISSA(J) = DEG(J)
      END DO

*     LREFINT(I) = 0: FULL OBSERVED (POSITIVE) AND CALCULATED PROFILES 
*                     ARE OBTAINED.
*     LREFINT(I) = 1: PART OF THE OBSERVED PROFILE IS LACKING AS A 
*                     RESULT OF EXCLUDING A 2-THETA REGION.
*     LREFINT(I) = 2: FULL OBSERVED (NEGATIVE) AND CALCULATED PROFILES
*                     ARE OBTAINED.
*     LREFINT(I) = 3: PART OF THE PROFILE EXCEEDS 2-THETA(MAX).
      DO I = 1, NREFL
         IF (NOPH(I) .GT. 1) CYCLE
         LREFINT(I) = 0
         OBSI(I) = 0.0
*        NUMERICAL INTEGRATION OF THE BRAGG INTENSITY FROM RDEG(1,I) 
*        TO RDEG(2,I)
*        RDEG: 2-THETA REGION FOR CALCULATING THE PROFILE (RADIAN)
         CALL QROMB(PRFINT, RDEG(1,I), RDEG(2,I), CALCI(I), A)
      END DO
      
*     DETERMINE THE EXCLUDED REGION WITH THE HIGHEST 2-THETA
*     NSKIP: NUMBER OF EXCLUDED REGIONS
*     DEGEXC: EXCLUDED 2-THETA REGION (DEGREE)
      LAST = NSKIP
      DO KK = 1, NSKIP - 1
         IF (DEGEXC(1,KK) .GT. DEGEXC(1,NSKIP)) LAST = KK
      END DO

*     LEXCEND = 1: THE HIGH 2-THETA OF RAW INTENSITY DATA IS CUT.
*     THE PROGRAM ASSUMES THAT PART OF INTENSITY DATA HAVE BEEN EXCLUDED
*     WITH THE LAST DEGEXC.
      IF (DEG(NPTS) .LT. DEGEXC(1,LAST)*RD) THEN
         LEXCEND = 1
      ELSE
         LEXCEND = 0
      END IF
      
      IF (R12 .EQ. 0.0) THEN
         MONO = 1
      ELSE
         MONO = 0
      END IF
      ISTART = 1
      DO 20 J = 1, NREFL
         I = IPOINT(J)
         IF (NOPH(I) .GT. 1) CYCLE
         IF (LEXCEND .EQ. 1 .AND. RDEG(2,I) .GT. DEG(NPTS)) THEN
            DO JJ = J, NREFL
               LREFINT(IPOINT(JJ)) = 3
            END DO
            EXIT
         END IF
	 
         DO KK = 1, NSKIP
            IF ((RDEG(1,I) .GT. DEGEXC(1,KK)*RD .AND. 
     &      RDEG(1,I) .LT. DEGEXC(2,KK)*RD) .OR.
     &      (RDEG(2,I) .GT. DEGEXC(1,KK)*RD .AND. 
     &      RDEG(2,I) .LT. DEGEXC(2,KK)*RD) .OR.
     &      (RDEG(1,I) .LT. DEGEXC(1,KK)*RD .AND.
     &      RDEG(2,I) .GT. DEGEXC(2,KK)*RD)) THEN
               LREFINT(I) = 1
*              In the present version, the observed intensity is set 
*              equal to the calculated one in *.mem for a reflection 
*              whose profile is included in an excluded region.
               OBSI(I) = CALCI(I)
               GO TO 20
            END IF
         END DO

*        YOPROF: INTEGRATION OF THE OBSERVED OVERLAPPED PROFILE
*        YCPROF: INTEGRATION OF THE CALCULATED OVERLAPPED PROFILE
         YOPROF = 0.0D0
         YCPROF = 0.0D0
         YCIND = 0.0D0

*        SUBROUTINE HRPDINT IS CALLED WHEN THE 8TH CHARACTER IS '$'
*        FOR EVERY FORMAT
         IF (TOP1(1:4) .EQ. 'OVLP' .OR. TOP1(8:8) .EQ. '$') THEN
            CALL HRPDINT(DEG,XINT,YFIT,NPTS,ISTART,YOPROF,YCPROF,YCIND,
     &      A,0)
         ELSE
            CALL AREAVS(DEG,XINT,YFIT,NPTS,ISTART,YOPROF,YCPROF,YCIND,
     &      LFLAG,A,0)
*           LFLAG = 1: REACHED REFLECTIONS WHOSE PROFILES ARE PARTIALLY 
*                      LACKING
            IF (LFLAG .EQ. 1) EXIT
         END IF

         SELECT CASE (MONO)
            CASE (0)
*              CHARACTERISTIC X RAYS
               IF (YOPROF .LE. 0.0) THEN
                  OBSI(I) = 0.0
               ELSE IF (YCPROF .GT. 0.0D0 .AND. CALCI(I) .GT. 0.0D0) 
     &         THEN
                  OBSI(I) = YOPROF*CALCI(I)/YCPROF
               ELSE 
                  OBSI(I) = 0.0
                  LREFINT(I) = 2
               END IF
            CASE (1)
*              MONOCHROMATIC BEAM
               IF (YOPROF .LE. 0.0) THEN
                  YPEAK(I) = 1.0E-10
               ELSE IF (YCPROF .GT. 0.0D0 .AND. CALCI(I) .GT. 0.0D0) 
     &         THEN
                  OBSI(I) = YOPROF*CALCI(I)/YCPROF
*                 UPDATE ARRAY YPEAK
                  YPEAK(I) = YPEAK(I)*OBSI(I)/CALCI(I)
                  IF (YPEAK(I) .LE. 0.0) YPEAK(I) = 1.0E-10
               ELSE 
                  YPEAK(I) = 1.0E-10
                  LREFINT(I) = 2
               END IF
         END SELECT
   20 CONTINUE

      IF (MONO .EQ. 1) RETURN
*     CHARACTERISTIC X RAYS
      DO J = 1, NREFL
         I = IPOINT(J)
         IF (NOPH(I) .GT. 1 .OR. L12(I) .NE. 0) CYCLE
*        SUM OF INTEGRATED INTENSITIES FOR K-ALPHA1 AND K-ALPHA2
         DO JJ = I+1, NREFL
            IF (L12(JJ) .EQ. I) EXIT
         END DO
*        UPDATE ARRAY YPEAK
         IF (OBSI(I) .EQ. 0.0) THEN
            YPEAK(I) = 1.0E-10
         ELSE
            YPEAK(I) = 
     &      YPEAK(I)*(OBSI(I)+OBSI(JJ))/(CALCI(I)+CALCI(JJ))
         END IF
         YPEAK(JJ) = YPEAK(I)*R12*FLP(JJ)/FLP(I)
      END DO
      END
 
************************************************************************
 
      SUBROUTINE QROMB(FUNC,A,B,SS,AA)
*     ROMBERG INTEGRATION
      INTEGER JMAX,JMAXP,K,KM
      REAL A,B,FUNC,SS,EPS,AA(*)
      EXTERNAL FUNC
*     EPS WAS INCREASED FROM 1.0E-6
*Rev 1.0x 2001.04.20 Izumi {
*     PARAMETER (EPS=9.0E-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
      PARAMETER (EPS=2.0E-5, JMAX=25, JMAXP=JMAX+1, K=5, KM=K-1)
* }
      INTEGER J
      REAL DSS,H(JMAXP),S(JMAXP)
      H(1)=1.
      DO J=1,JMAX
        CALL TRAPZD(FUNC,A,B,S(J),J,AA)
        IF (J .GE. K) THEN
          CALL POLINT(H(J-KM),S(J-KM),K,0.,SS,DSS)
          IF (ABS(DSS).LE.EPS*ABS(SS)) RETURN
        ENDIF
        S(J+1)=S(J)
        H(J+1)=0.25*H(J)
      END DO
      CALL JOBEND('Too many steps in SUBROUTINE QROMB')
      END
 
************************************************************************
 
      SUBROUTINE TRAPZD(FUNC,A,B,S,N,AA)
      INTEGER N
      REAL A,B,S,FUNC,AA(*)
      EXTERNAL FUNC
      INTEGER IT,J
      REAL DEL,SUM,TNM,X
      IF (N.EQ.1) THEN
        S=0.5*(B-A)*(FUNC(A,AA)+FUNC(B,AA))
      ELSE
        IT=2**(N-2)
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5*DEL
        SUM=0.
        DO J=1,IT
          SUM=SUM+FUNC(X,AA)
          X=X+DEL
        END DO
        S=0.5*(S+(B-A)*SUM/TNM)
      ENDIF
      END

************************************************************************

      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      INTEGER N,NMAX
      REAL DY,X,Y,XA(N),YA(N)
      PARAMETER (NMAX=10)
      INTEGER I,M,NS
      REAL DEN,DIF,DIFT,HO,HP,W,C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO I = 1, N
         DIFT=ABS(X-XA(I))
         IF (DIFT.LT.DIF) THEN
            NS=I
            DIF=DIFT
         END IF
         C(I)=YA(I)
         D(I)=YA(I)
      END DO
      Y=YA(NS)
      NS=NS-1
      DO M = 1, N-1
         DO I = 1, N-M
            HO=XA(I)-X
            HP=XA(I+M)-X
            W=C(I+1)-D(I)
            DEN=HO-HP
            IF(DEN.EQ.0.) CALL JOBEND('Failure in SUBROUTINE POLINT')
            DEN=W/DEN
            D(I)=HP*DEN
            C(I)=HO*DEN
         END DO
         IF (2*NS.LT.N-M)THEN
            DY=C(NS+1)
         ELSE
            DY=D(NS)
            NS=NS-1
         END IF
         Y=Y+DY
      END DO
      END

************************************************************************

      REAL FUNCTION PRFINT(X,A)
*     CALCULATE THE PROFILE INTENSITY AT X
      INTEGER IOPT(2)
      REAL A(*),EOPT(2)
      PARAMETER (NB=7000,NP=80000,NPH=8)
      COMMON /ABSVR/ ABSORP(NP),VARLEN(NP),ROUGH(NP),NSURFR
      COMMON /G/ I
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /YOKO/ ABSCISSA(NP),NPOINTS,AWRY
      COMMON /MC/ CTHM,NTRAN,PCOR
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD

*     INTERPOLATE THE ABSORPTION FACTOR
      LUP = 0      ! Start with a binary search
      IOPT(2) = 0  ! End of the option list
      CALL SILUP(X,VA,NPOINTS,ABSCISSA,ABSORP,3,LUP,IOPT,EOPT)
*     INTERPOLATE THE IRRADIATION WIDTH
      LUP = 0      ! Start with a binary search
      IOPT(2) = 0  ! End of the option list
      CALL SILUP(X,VW,NPOINTS,ABSCISSA,VARLEN,3,LUP,IOPT,EOPT)
*     INTERPOLATE THE SURFACE ROUGHNESS
      LUP = 0      ! Start with a binary search
      IOPT(2) = 0  ! End of the option list
      CALL SILUP(X,VR,NPOINTS,ABSCISSA,ROUGH,3,LUP,IOPT,EOPT)
*     PRFINT = VA*VW*VR*YPEAK(I)*PRFL(X - PEAK(I))
      AWRY = VA*VW*VR*YPEAK(I)
      PRFINT = AWRY*PRFL(X,A)
      END

************************************************************************

      SUBROUTINE HRPDINT(DEG,XINT,YFIT,NPTS,ISTART,YOPROF,YCPROF,YCIND,
     &  A,L20)            
*     DETERMINE AN INTEGRATED INTENSITY
*     ORIGINALLY WRITTEN FOR INTENSITY DATA TAKEN ON HRPD IN A 
*     OVERLAPPED MODE
      PARAMETER (NB=7000,NP=80000,NPH=8)
      REAL DEG(*),XINT(*),YFIT(*),A(*)
      DOUBLE PRECISION YOPROF,YCPROF,YCIND
      COMMON /ABSVR/ ABSORP(NP),VARLEN(NP),ROUGH(NP),NSURFR
      COMMON /C/ BG(NP),BGINC(NP),DEGNOR(NP)
      COMMON /G/ I
      COMMON /I/ RDEG(2,NB)
      INTEGER WIDTH,HEIGHT,YMIN,YMAX,LBGPLOT,LDEL,OFFSETD,MSIZE,FSIZE,
     &  LSIZE,OFFSET,INDREF
      COMMON /IGOR/ WIDTH,HEIGHT,YMIN,YMAX,LBGPLOT,LDEL,OFFSETD,MSIZE,
     &  FSIZE,LSIZE,OFFSET(NPH),INDREF
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /MC/ CTHM,NTRAN,PCOR
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD

      NBEGIN = 0
*Rev 1.1 2003.05.18 {
*     Suggestion of Dr. Kajimoto on 2003.05.16
      NEND = NPTS
* }
      DO K = ISTART, NPTS
         IF (DEG(K) .GT. RDEG(2,I)) THEN
            NEND = K - 1
            ISTART = NBEGIN
            EXIT
         END IF
         IF (NBEGIN .EQ. 0 .AND. DEG(K) .GE. RDEG(1,I)) NBEGIN = K
         YC = ABSORP(K)*VARLEN(K)*ROUGH(K)*YPEAK(I)*
     &   PRFL(DEG(K),A)
         YCIND = YCIND + YC
         IF (NPAT .EQ. 5 .AND. INDREF .EQ. 1 .AND. L20 .EQ. 1) 
     &   WRITE(20,'(F7.3,1P,E12.5E1)') DEG(K)*57.29578, YC + BG(K)
      END DO
      IF (NPAT .EQ. 5 .AND. INDREF .EQ. 1 .AND. L20 .EQ. 1) 
     &WRITE(20,'(A)') 'NaN NaN'
	    
*     PRIMITIVE ALGORITHM FOR THE TRAPEZOIDAL INTEGRATION
*     FOR VARIABLE STEP WIDTHS
      DO JJ = NBEGIN, NEND - 1
         VSTEP = DEG(JJ+1) - DEG(JJ)
         DELTA1 = XINT(JJ) - BG(JJ)
         DELTA2 = XINT(JJ+1) - BG(JJ+1)
         YOPROF = YOPROF + 0.5*VSTEP*(DELTA1 + DELTA2)
         DELTA1 = YFIT(JJ) - BG(JJ)
         DELTA2 = YFIT(JJ+1) - BG(JJ+1)
         YCPROF = YCPROF + 0.5*VSTEP*(DELTA1 + DELTA2)
      END DO
      END
 
************************************************************************
 
      SUBROUTINE AREAVS(DEG,XINT,YFIT,NPTS,ISTART,YOPROF,YCPROF,YCIND,
     &  LFLAG,A,L20)
*     DETERMINE AN INTEGRATED INTENSITY, TAKING INTO ACCOUNT VARIABLE
*     STEP WIDTHS
      PARAMETER (NB=7000,NP=80000,NPH=8)
      INTEGER IREG(2,2)
      REAL DEG(*),XINT(*),YFIT(*),STEP(2),A(*)
      DOUBLE PRECISION YOPROF,YCPROF,YCIND
      COMMON /ABSVR/ ABSORP(NP),VARLEN(NP),ROUGH(NP),NSURFR
      COMMON /B/ U(NB),NREFL
      COMMON /C/ BG(NP),BGINC(NP),DEGNOR(NP)
      COMMON /OBSCAL/ OBSI(NB),CALCI(NB),LREFINT(NB)
      COMMON /G/ I
      COMMON /I/ RDEG(2,NB)
      INTEGER WIDTH,HEIGHT,YMIN,YMAX,LBGPLOT,LDEL,OFFSETD,MSIZE,FSIZE,
     &  LSIZE,OFFSET,INDREF
      COMMON /IGOR/ WIDTH,HEIGHT,YMIN,YMAX,LBGPLOT,LDEL,OFFSETD,MSIZE,
     &  FSIZE,LSIZE,OFFSET(NPH),INDREF
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /MC/ CTHM,NTRAN,PCOR
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /ORDER/ IPOINT(NB)
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
      DATA RD/0.01745329/

      LFLAG = 0
*     NUMBER OF REGIONS TO BE INTEGRATED
      NREGION = 1
*     IREG(2,2): STARTING AND ENDING POINTS FOR REGIONS WITH EQUAL STEP
*                WIDTHS AND REGION NUMBERS.
*                NOTE THAT IREG(2,1) = IREG(1,2)
      IREG(1,1) = 0
      DO K = ISTART, NPTS
         IF (DEG(K) .GT. RDEG(2,I)) THEN
            IREG(2,NREGION) = K - 1
*           RESET ISTART AT THE STARTING POINT FOR THE PRESENT REFLECTION
            ISTART = IREG(1,1)
            GO TO 2
         ELSE IF (DEG(K) .LT. RDEG(1,I)) THEN
            CYCLE
         END IF

*        STEP(2): STEP WIDTHS FOR REGION NUMBERS 1 AND 2
         IF (IREG(1,1) .EQ. 0) THEN
*           FIRST POINT FOR REFLECTION I
            IREG(1,1) = K
            IF (K + 1 .LE. NPTS .OR. K .EQ. 1) THEN
               STEP(1) = DEG(K+1) - DEG(K)
            ELSE
               STEP(1) = DEG(K) - DEG(K-1)
            END IF
         ELSE IF (ABS(DEG(K) - DEG(K-1) - STEP(NREGION)) .GT.
     &   0.001*RD) THEN
            IREG(2,NREGION) = K - 1
            NREGION = NREGION + 1
            IF (NREGION .GT. 2) CALL JOBEND('Too large NREGION')
            IREG(1,2) = K - 1
            STEP(2) = DEG(K) - DEG(K-1)
         END IF
      END DO

      DO JJ = I, NREFL
         LREFINT(IPOINT(JJ)) = 3
      END DO
      LFLAG = 1
      RETURN
*     THE TOTAL NUMBER OF POINTS FOR INTEGRATION BY THE EXTENDED 
*     SIMPSON'S RULE SHOULD BE 2N + 1 (N: INTEGER)
    2 IF (MOD(IREG(2,1) - IREG(1,1),2) .EQ. 1) THEN
         IREG(1,1) = IREG(1,1) - 1
         IF (IREG(1,1) .LE. 0) IREG(1,1) = IREG(1,1) + 2
      END IF 
      IF (MOD(IREG(2,2) - IREG(1,2),2) .EQ. 1) THEN
         IREG(2,2) = IREG(2,2) + 1
         IF (IREG(2,2) .GT. NPTS) THEN
            IREG(2,2) = IREG(2,2) - 2
*           NEGLECT INTEGRATION FOR ONLY ONE POINT OR LESS
            IF (IREG(2,2) - IREG(1,2) + 1 .LT. 3) 
     &      IREG(2,2) = IREG(1,2) - 1
         END IF
      END IF

      DO JJ = 1, NREGION
         IP = 0
         DO K = IREG(1,JJ), IREG(2,JJ)
*           COEFF: COEFFICIENTS IN THE EXTENDED SIMPSON'S RULE
*           "Numerical Recipes in FORTRAN," 2nd ed., p. 128.
*           STEP(JJ) IS MULTIPLIED TO TAKE INTO ACCOUNT THE STEP WIDTH
            IF (K .EQ. IREG(1,JJ) .OR. K .EQ. IREG(2,JJ)) THEN
*              FIRST OR LAST POINT
               COEFF = 0.3333333*STEP(JJ)
            ELSE IF (MOD(IP,2) .EQ. 1) THEN
*              IP = odd
               COEFF = 1.333333*STEP(JJ)
            ELSE
*              IP = even
               COEFF = 0.6666667*STEP(JJ)
            END IF
            YOPROF = YOPROF + COEFF*(XINT(K) - BG(K))
            YCPROF = YCPROF + COEFF*(YFIT(K) - BG(K))

            YC = ABSORP(K)*VARLEN(K)*ROUGH(K)*YPEAK(I)*
     &      PRFL(DEG(K),A)
            YCIND = YCIND + YC
            IF (NPAT .EQ. 5 .AND. INDREF .EQ. 1 .AND. L20 .EQ. 1) 
     &      WRITE(20,'(F7.3,1P,E12.5E1)') DEG(K)/RD, YC + BG(K)
            IP = IP + 1
         END DO
      END DO
      IF (NPAT .EQ. 5 .AND. INDREF .EQ. 1 .AND. L20 .EQ. 1) 
     &WRITE(20,'(A)') 'NaN NaN'
      END
 
************************************************************************
 
      SUBROUTINE RDINT(THINIT,STEP,NEXC,NSTEP,DEG,XINT)
*     READ STEP-SCANNED INTENSITY DATA
*     (A) INPUT
*         NP: MAXIMUM DATA POINTS FOR XINT
*     (B) OUTPUT
*         THINIT: MINIMUM TWO-THETA
*           STEP: STEP WIDTH
*          NSTEP: NUMBER OF INTENSITY DATA
*            DEG: TWO-THETA
*           XINT: STEP-SCANNED INTENSITY DATA
      PARAMETER (NP=80000)
      INTEGER NTOTAL,NEXC,NSKIP,NSTEP
      REAL THINIT,STEP,DEG(*),XINT(*)
      CHARACTER LINE*80,TOP1*8,TOP2*6
      COMMON /DETECT/ VARINV(NP)
      COMMON /EXCREG/ DEGEXC(2,50),NSKIP
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
      COMMON /INTFMT/ TOP1
      CHARACTER*50 FILE3*50
      DATA RD/0.01745329/

      CALL OPENFILE(2,FILE3)
      REWIND 3
*     CHECK WHETER HRPD DATA OR NOT
      READ(3,'(A)') TOP1
      READ(3,'(A)') TOP2
      REWIND 3
*     VARINV(1) <> 0.0 IF NUMBERS OF COUNTERS ARE GIVEN
      VARINV(1)=0.0

      IF (NBEAM .EQ. 0 .AND. TOP1 .EQ. 'EXPNO = ') THEN
*        #1. HRPD (JAERI, JRR-3M) FORMATS (TWO TYPES EXIST)
         DO J=1,10000
            READ(3,'(A)') LINE
            IF (LINE(1:6) .EQ. 'DETN =') GO TO 1 
         END DO
    1    READ(3,'(A)') LINE
         IF (LINE(1:4) .EQ. '*** ') THEN
*           #1.1 HRPD FORMAT, TYPE 1
            DO J=1,10000
               READ(3,'(A)') LINE
               IF (LINE(1:3) .EQ. '***') GO TO 3
               READ(LINE,'(10F8.3)') (DEG(10*(J-1)+K),K=1,10)
            END DO
    3       NTOTAL=10*(J-2)+10
            READ(3,*) (XINT(J),J=1,NTOTAL)
         ELSE IF (LINE(1:4) .EQ. 'STPN') THEN
*           #1.2 HRPD FORMAT, TYPE 2
            DO J=1,NP
               READ(3,'(5X,F9.3,7X,F9.0)',END=5) DEG(J),XINT(J)
            END DO
            CALL JOBEND('Too small NP value')
    5       NTOTAL=J-1
         END IF
*        STEP = 0.05 DEGREES IN HRPD
         STEP = 0.05
      ELSE IF (NBEAM .EQ. 1 .AND. TOP2 .EQ. '*CLASS') THEN
*        #2. RIGAKU FORMAT (RINT)
*        TAB'S ARE INCLUDED IN FILES IN THE RIGAKU FORMAT
         DO J=1,10000
            READ(3,'(A)') LINE
            IF (LINE(1:6) .EQ. '*START') THEN
               ITOP = INDEX(LINE,'=') + 3
               READ(LINE(ITOP:),*) THINIT
               READ(3,'(A)')
               READ(3,'(A)') LINE
               ITOP = INDEX(LINE,'=') + 3
               READ(LINE(ITOP:),*) STEP
               GO TO 7
            END IF
         END DO
    7    DO J=1,10000
            READ(3,'(A)') LINE
            IF (LINE(1:6) .EQ. '*COUNT') GO TO 8
         END DO
    8    ITOP = INDEX(LINE,'=') + 3
         READ(LINE(ITOP:),*) NTOTAL
         IF (NTOTAL .GT. NP-2) THEN
            WRITE(6,150) NP
            CALL JOBEND(' ')
         END IF
         READ(3,*) (XINT(J),J=1,NTOTAL)
         DO J=1,NTOTAL
            DEG(J)=THINIT+FLOAT(J-1)*STEP
         END DO
      ELSE IF (TOP1(1:7) .EQ. 'GENERAL') THEN
*        #3. GENERAL FORMAT
         READ(3,'(A)') TOP1
         READ(3,*) NTOTAL
         READ(3,*) (DEG(J),XINT(J),J=1,NTOTAL)
         STEP = DEG(2)-DEG(1)
      ELSE IF (TOP1(1:4) .EQ. 'IGOR') THEN
*        #3'. IGOR TEXT FILE
         READ(3,'()') 
         READ(3,'()')
         READ(3,'()')
         DO J = 1, NP
            READ(3,'(A)') LINE
            IF (LINE(1:3) .EQ. 'END') EXIT
            READ(LINE,*) DEG(J),XINT(J)
         END DO
         NTOTAL = J -1
         STEP = DEG(2)-DEG(1)
      ELSE IF (TOP1(1:4) .EQ. 'HB-4') THEN
*        #4, HFIR FORMAT
         DO WHILE (LINE(1:2) .NE. ' 1' .OR. LINE(9:10) .NE. ' 1' .OR.
     &   LINE(17:18) .NE. ' 1' .OR. LINE(25:26) .NE. ' 1' .OR.
     &   LINE(33:34) .NE. ' 1' .OR. LINE(41:42) .NE. ' 1' .OR.
     &   LINE(49:50) .NE. ' 1' .OR. LINE(57:58) .NE. ' 1' .OR.
     &   LINE(65:66) .NE. ' 1' .OR. LINE(73:74) .NE. ' 1')
            READ(3,'(A)') LINE
         END DO
         BACKSPACE 3
         STEP=0.05
*        STARTING 2-THETA = 11.0
         CALL MULTDT(DEG,XINT,VARINV,NTOTAL,STEP,11.0,TOP1)
      ELSE IF (TOP1(1:3) .EQ. 'D2B' .OR. TOP1(1:3) .EQ. 'D1A' .OR.
     &TOP1(1:3) .EQ. 'D20' .OR. TOP1(1:4) .EQ. 'VSCT') THEN
*        DIFFRACTOMETERS AT ILL OR VARIABLE STEP-COUNTING TIMES
         READ(3,'()')
         READ(3,'()')
         READ(3,'(16X,F8.3)') STEP
         READ(3,'(F8.3)') THINIT
         READ(3,'()')
         CALL MULTDT(DEG,XINT,VARINV,NTOTAL,STEP,THINIT,TOP1)
      ELSE IF (TOP1(1:4) .EQ. 'FVFM' .OR. TOP1(1:4) .EQ. 'OVLP') THEN
*        FVFM: FULLY VARIABLE FORMAT WHERE BOTH THE STEP WIDTH AND COUNTING
*              TIME ARE VARIABLE
*        OVLP: THE SAME AS THE FVFM FORMAT.  MEASUED WITH HRPD IN A
*              OVERLAPPING MODE. 
         CALL MULTDT(DEG,XINT,VARINV,NTOTAL,STEP,THINIT,TOP1)
*        THE STEP WIDTH AT THE HIGHEST ANGLE IS USED BECAUSE SKIPPED
*        REGIONS FOR IGOR PRO ARE WIDER THAN 2.0*STEP (CF. SUBROUTINE
*        OUT20).
	 STEP = DEG(NTOTAL)-DEG(NTOTAL-1)
      ELSE IF (TOP1(1:4) .EQ. 'BT-1') THEN
         READ(3,'(A)') TOP2         
         DO WHILE (TOP2(1:5) .NE. 'BANK ')
	    READ(3,'(A)') TOP2  
	 END DO
	 BACKSPACE 3
         READ(3,'(12X,I4,15X,2F10.2)') NTOTAL, START_BT1, STEP
	 START_BT1 = 0.01*START_BT1
	 STEP = 0.01*STEP
         CALL MULTDT(DEG,XINT,VARINV,NTOTAL,STEP,START_BT1,TOP1)
      ELSE
*        ROUTINE CODED BY DR. Xiaohua Yan OF MAC SCIENCE
*        #4, RIETAN FORMAT
*        #5, MAC SCIENCE FORMAT
         NTOTAL = IntRd(3, NP, STEP, DEG, XINT)
         if (NTOTAL .eq. 0) then
            CALL JOBEND('Error during reading intensity data.')
         else if (NTOTAL .lt. 0) then
            WRITE(6,150) -NTOTAL
            CALL JOBEND(' ') 
         endif
      END IF
  150 FORMAT(//11X,'The number of step-scan data has extended the maximu
     &m value (',I5,').')
      CLOSE(UNIT=3,STATUS='KEEP')

*     DELETE INTENSITY DATA TO BE SKIPPED
      NSTEP=0
      DO 60 LL=1,NTOTAL
         X1=DEG(LL)
         IF (XINT(LL) .EQ. 0) CALL JOBEND
     &   ('A zero intensity has been found.  Delete it!')
         IF (NEXC .EQ. 1) THEN
            DO J=1,NSKIP
               IF (X1 .GT. DEGEXC(1,J)-0.0005 .AND. X1 .LT.
     &         DEGEXC(2,J)+0.0005) GO TO 60
            END DO
         END IF
         NSTEP=NSTEP+1
         DEG(NSTEP)=X1*RD
         XINT(NSTEP)=XINT(LL)
         IF (VARINV(1) .NE. 0.0) VARINV(NSTEP)=VARINV(LL)
   60 CONTINUE
      STEP=STEP*RD
      END

************************************************************************

      SUBROUTINE MULTDT(DEG,XINT,VARINV,NTOTAL,STEP,THINIT,TOP1)
*     READ NUMBER OF DETETORS AND INTENSITIES
      PARAMETER (NP=80000)

      REAL DEG(*),XINT(*),VARINV(*)
      CHARACTER*8 TOP1

      DO J=1,NP
         XINT(J)=0.0
      END DO
      IF (TOP1(1:4) .EQ. 'BT-1') THEN
*        DATA MEASURED ON BT-1 AT NIST
*        VARINV: STANDARD DEVIATIONS
         READ(3,'(5(2F8.0))') (XINT(J),VARINV(J),J=1,NTOTAL)
         DO J = 1, NTOTAL
            IF (XINT(J) .EQ. 0.0) THEN
               NTOTAL = J - 1
               EXIT
            END IF
            VARINV(J) = 1.0/VARINV(J)**2
            DEG(J) = THINIT + FLOAT(J-1)*STEP
         END DO
         RETURN
      ELSE IF (TOP1(1:4) .EQ. 'VSCT') THEN
*        VARIABLE STEP-COUNTING TIMES (IKEDA)
         READ(3,'(10(F3.0,F10.3))',END=1) (VARINV(J),XINT(J),J=1,NP)
*Rev 1.09 2002.12.29 Izumi {
      ELSE IF (TOP1(1:4) .EQ. 'FVFM' .OR. TOP1(1:4) .EQ. 'OVLP') THEN
*        FVFM: FULLY VARIABLE FORMAT
*        OVLP: INTENSITY DATA TAKEN ON HRPD WITH AN OVERLAPPING MODE
         READ(3,'()')
         READ(3,*,END=1) (DEG(J),XINT(J),VARINV(J),J=1,NP)
* }
      ELSE
*        D2B OR HB-4 DATA
         READ(3,'(10(F2.0,F6.0))',END=1) (VARINV(J),XINT(J),J=1,NP)
      END IF
      
    1 DO J=NP,1,-1
         IF (XINT(J) .NE. 0.0) EXIT
      END DO
      NTOTAL=J
      DO J=1,NTOTAL
*        VARINV: (NUMBER OF DETECTORS)/Y(OBS)
         IF (XINT(J) .EQ. 0.0) THEN
            VARINV(J) = 1.0
         ELSE
            VARINV(J) = VARINV(J)/XINT(J)
         END IF
         IF (TOP1(1:4) .NE. 'FVFM' .AND. TOP1(1:4) .NE. 'OVLP')
     &   DEG(J)=THINIT+FLOAT(J-1)*STEP
      END DO
      END

************************************************************************

c     integer function IntRd(io, np, step, deg, xint)
*     CODED BY DR. Xiaohua Yan OF MAC SCIENCE AND MODIFIED BY F. IZUMI
c
c     A function to read intensity data for the following formats:
c             1) RIETAN format
c                     first non-comment line: 
c                             NTOTAL   start_two_theta   two_theta_step
c                     other lines: 
c                             intensities
c             2) MAC Science format 
c                     each non-comment line contains a pair of 
c                             two-theta and intensity
c
c             In the data file, if the first column of a line is '*',
c             this line is taken as a comment line

c     input: 
c             io      logical I/O number for this data file
c             np      maximum number of intensity data
c     output:
c             deg     real array, two-theta
c             xint    real array, intensities
c             step    step width, step = deg(2) - deg(1)
c       
c     returned value:
c           ntotal    data have been read in sucessufully
c                0    error during reading the file
c              < 0    the number of data exceeds np - 2 
c
c     note:
c             1. before calling this function, the data file should be opened as unit io
c             2. the returned value can/should be checked to examine whether file reading is successful

      integer function IntRd(io, np, step, deg, xint)

      integer io, np
      real step, deg(*), xint(*)
      integer ntotal
      character buf*80
	 
c     find the first non-comments line
      rewind(io)
   10 continue   
         read (io, '(a80)', err = 8, end = 8) buf
      if (buf(1:1) .eq. '*') goto 10

      read(buf, *, err = 1, end = 1) ntotal, deg(1), step
      if (ntotal .gt. np - 2) goto 8
*     RIETAN format
      read(io, *) (xint(i), i = 1, ntotal)
      do i = 2, ntotal
         deg(i) = deg(1) + float(i-1)*step
      end do
      goto 9

*     MAC Science format
    1 backspace(io)
      ntotal = 0
20    continue
         read (io, '(a80)', err = 8, end = 2) buf 
         if (buf(1:1).ne.'*') then
            ntotal = ntotal + 1
            read (buf, *, err = 8, end = 2) 
     &      deg(ntotal), xint(ntotal)
            if (ntotal .gt. np -2) goto 8
         endif
      goto 20
    2 step = deg(2) - deg(1)
      goto 9

    8 continue        
*     An error occurred.  Reading error or too many intensity data
      if (ntotal .gt. np-2) then
         ntotal = -ntotal
      else
         ntotal = 0
      end if

    9 IntRd = ntotal
      end

************************************************************************

      SUBROUTINE RDMEED
*     READ STRUCTURE FACTORS THAT HAVE RESULTED FROM MEM ANALYSIS
      CHARACTER TOP*1,HEAD*70
      CHARACTER*50 FILE32
      LOGICAL FOREVER

      PARAMETER (NB=7000)
      COMMON /FMEM/ MH(NB),MK(NB),ML(NB),FMEMR(NB),FMEMI(NB),NHKL
      
      CALL OPENFILE(9,FILE32)
      REWIND 32
      
      READ(32,'(A)') HEAD
      IF (HEAD.NE.'$FB-MEM-DATA') CALL JOBEND
     &(FILE32(1:INDEX(FILE32,'.'))//'fba contains an invalid header')
      FOREVER = .TRUE.
      DO WHILE (FOREVER)
         READ(32,'(A)') TOP
	 IF (TOP .NE. '*') THEN
	    BACKSPACE 32
	    FOREVER = .FALSE.
	 END IF
      END DO
      READ(32,*) NHKL
*Rev 1.2 2003.5.28 {
*     FMEMR: REAL PART
*     FMEMI: IMAGINARY PART
* }
      DO J = 1, NHKL
         READ(32,*) MH(J), MK(J), ML(J), FMEMR(J), FMEMI(J)
      END DO
      DO J = NHKL + 1, NB
         READ(32,*,END=9) MH(J), MK(J), ML(J), FMEMR(J), FMEMI(J)
      END DO
    9 NHKL = J - 1
      CLOSE(UNIT=32)
      END

************************************************************************

      SUBROUTINE RDSFF
*     READ INITIAL INTEGRATED INTENSITIES FOR LE BAIL REFINEMENT
      CHARACTER FILE22*50
      PARAMETER (NB=7000,NPH=8)
      COMMON /ABPI/ AP(NB),BP(NB)
      COMMON /FMEM/ MH(NB),MK(NB),ML(NB),FMEMR(NB),FMEMI(NB),NHKL
      character temp*80
          
      CALL OPENFILE(10,FILE22)
      REWIND 22
      
*     SKIP THE COMMENT LINE
      READ(22,'(A)')
      DO I = 1, NB
         READ(22,'(3I4,2F15.4)',END=9) MH(I),MK(I),ML(I),FWHM,AP(I)
      END DO
    9 NHKL = I - 1
      CLOSE(UNIT=22)
      END

************************************************************************

      SUBROUTINE MKFLOR(A,G,TITLE,SIGMAA,NTERMS,NESD)
*     MAKE A FILE FOR ORFFE AND ORTEP-II
      PARAMETER (NT=999,NR=400,NAP=150,NPH=8,NS=48,NB=7000,NCS=400)
      PARAMETER (NOD=NR*(NR+1)/2,NSP=NAP*10)
      REAL A(*),SIGMAA(*),STPAR(NSP),PM(NOD),TV(3)
      DOUBLE PRECISION G(*)
      INTEGER IG(NAP),KI1(NSP),R,RX,RY,RZ,ROR(2,3),NEWNO(NT),NC(NT)
      CHARACTER TITLE*25,PARNAM*60,PHNAME*25,INSTR*72
      CHARACTER*50 FILE9
      COMMON /CC/ APR(NT)
      COMMON /CNSTR/ C1(NCS),C2(10,NCS),IN(10,NCS),NCNPAR(NCS),
     &  NACNS(NCS),NCNSTR
      COMMON /ERMAT/ SI2,SP2
      COMMON /EXC/ LSPSYM(NPH),NSSYM(14,NPH)
      COMMON /FRDA/ NFR,NMEM,MEED(15),EPSD,SCIO,TSCAT,TSCAT1,TSCAT2,
     &  UCLAG,NDA,NFPX(3)
      COMMON /H/ NQ(NPH),NCENTR(NPH),LAUEG(NPH)
      COMMON /LTM/ LISO(NAP,NPH),LMAG(NPH)
      COMMON /MULT/ SM(NAP)
      COMMON /PAR/ PARNAM(NT),PHNAME(NPH)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /RT/ T(3,NS,NPH),HT(NS,NB),NSYM(NPH),NSITE(NPH),
     &  IDSYM(NS,NAP,NPH),R(3,3,NS,NPH),RX(NS,NB),RY(NS,NB),RZ(NS,NB)
      COMMON /SITE/ NOAT(NAP,NPH)
      COMMON /T/ FF(NB),COEF(NPH)
      
      LP(J)=KPHB(NDA)+J
      PS(J)=SIGMAA(LP(J))

      CALL OPENFILE(6,FILE9)
      REWIND 9
      WRITE(9,'(A)') TITLE
*     LAUE GROUP NUMBER, CENTROSYMMETRIC/NONCENTROSYMMETRIC,
*     STEP OF SYMMETRY OPERATIONS
      WRITE(9,'(4I5)') LAUEG(NDA),NCENTR(NDA),2*NINT(COEF(NDA)),
     &  NSSYM(7,NDA)
*     NC(RIETAN PARAM. #) = (RIETAN SERIES # OF CONSTRAINTS)
      IPOINT=0
      DO 10 I=1,NTERMS
         IF (ID(I) .EQ. 2) THEN
            IPOINT=IPOINT+1
            NC(I)=IPOINT
         END IF
   10 CONTINUE

      DO 20 J=1,NSITE(NDA)
*        DETERMINE THE PARAMETER NUMBERS OF OCCUPATION FACTORS
         IF (J .EQ. 1) THEN
            IG(J)=LP(NG1X)
         ELSE
            IG(J)=IG(J-1)+NPSITE(J-1,NDA)
         END IF
*        NEWNO(RIETAN PARAM. #) = PARAM. # OF ORFFE
         DO I = IG(J),IG(J)+LISO(J,NDA)+3
            NEWNO(I)=(J-1)*10+I-IG(J)+1
         END DO
   20 CONTINUE
      NP=10*NSITE(NDA)
      DO 25 I=1,NP
         STPAR(I)=0.0
         KI1(I)=0
   25 CONTINUE

      DO 35 J=1,NSITE(NDA)
         DO 30 I=0,LISO(J,NDA)+3
*           STORE SUTRUCTURE PARAMETERS INCLUDING DUMMY ONES
            STPAR((J-1)*10+1+I)=A(IG(J)+I)
            IF (ID(IG(J)+I) .EQ. 1) THEN
               KI1((J-1)*10+1+I)=1
            ELSE IF (ID(IG(J)+I) .EQ. 2) THEN
               IC=NC(IG(J)+I)
*              WRITE CONSTRAINTS IMPOSED ON STRUCTURE PARAMETERS
               WRITE(9,'(2I5,1P,G14.6,10(1P,G14.6,I5:))')
     &         NEWNO(IG(J)+I),NCNPAR(IC),C1(IC),(C2(K,IC),
     &         NEWNO(IN(K,IC)),K=1,NCNPAR(IC))
            END IF
   30    CONTINUE
         IF (LISO(J,NDA) .EQ. 0) STPAR((J-1)*10+5)=A(LP(NOTX))
*        G ==> OCCUPANCY * MULTIPLICITY / (NUMBER OF GENERAL EQUIVALENT
*              POSITIONS)
         STPAR((J-1)*10+1)=SM(J)
   35 CONTINUE
      WRITE(9,'(A,I1)') 'ENDCON',NDA

*     WRITE STRUCTURE PARAMETERS AND REFINEMENT IDENTIFIERS
      WRITE(9,'(I5)') NP
      DO 40 J=1,NSITE(NDA)
         WRITE(9,'(A,1X,10(1P,G14.6:))') PARNAM(IG(J))(1:9),
     &   (STPAR(K),K=(J-1)*10+1,(J-1)*10+1+MAX(4,LISO(J,NDA)+3))
   40 CONTINUE
      WRITE(9,'(150I1)') (KI1(J),J=1,NP)

*     WRITE THE VARIANCE-COVARIANCE MATRIX
      IPOINT=1
*Rev 1.0l 2001.01.09 Izumi 
      DO J = IG(1), IG(NSITE(NDA))+LISO(NSITE(NDA),NDA)+3
         IF (ID(J) .NE. 1 .OR. INDEX(PARNAM(J),'Magnetic moment') .GT. 
     &   0) CYCLE
         DO K = J, IG(NSITE(NDA))+LISO(NSITE(NDA),NDA)+3
*           ARRAY ELEMENTS FOR REFINABLE STRUCTURE PARAMETERS ARE COPIED
*           CHI-SQUARE/(N-M) IS MULTIPLIED
            IF (ID(K) .NE. 1 .OR. INDEX(PARNAM(K),'Magnetic moment') 
     &      .GT. 0) CYCLE
            IF (NESD .EQ. 1) THEN
*              SCOTT'S METHOD OF CALCULATING E.S.D.'S
               PM(IPOINT)=SI2*G(IR(K)*(IR(K)-1)/2+IR(J))
            ELSE
*              CONVENTIONAL METHOD OF CALCULATING E.S.D.'S
               PM(IPOINT)=SP2*G(IR(K)*(IR(K)-1)/2+IR(J))
            END IF
            IPOINT=IPOINT+1
         END DO
      END DO
* }
      NV=0
      DO 55 J=1,NP
         IF (KI1(J) .EQ. 1) NV=NV+1
   55 CONTINUE
      WRITE(9,'(I5)') NV
      WRITE(9,'(1P,10G14.6)') (PM(I),I=1,NV*(NV+1)/2)

*     A, B, C, COS(ALPHA), COS(BETA), COS(GAMMA) AND THEIR STAND. DEV.
      WRITE(9,'(1P,6G14.6/1P,6G14.6)')
     &(APR(LP(J)),J=NAX,NCX),(COS(APR(LP(J))*0.01745329),J=NALPX,NGAMX),
     &(PS(J),J=NAX,NCX),(PS(J)*SIN(PS(J)*0.01745329),J=NALPX,NGAMX)

*     WRITE SYMMETRY OPERATIONS
      DO 70 ISYM=1,NSYM(NDA)
         DO 65 I=1,3
            ROR(2,I)=0
            IPOINT=1
            DO 60 J=1,3
               IF (R(I,J,ISYM,NDA) .EQ. -1) THEN
                  ROR(IPOINT,I)=-J
                  IPOINT=IPOINT+1
               ELSE IF (R(I,J,ISYM,NDA) .EQ. 1) THEN
                  ROR(IPOINT,I)=J
                  IPOINT=IPOINT+1
               END IF
   60       CONTINUE
            TV(I)=T(I,ISYM,NDA)
   65    CONTINUE
         CALL WRSYOP(TV,ROR,NSSYM(7,NDA),ISYM-1)
         IF (NCENTR(NDA) .EQ. 1) THEN
            DO 68 I=1,3
*              INVERTED POSITION: -X, -Y, -Z
               TV(I)=-TV(I)
               ROR(1,I)=-ROR(1,I)
               ROR(2,I)=-ROR(2,I)
   68       CONTINUE
            CALL WRSYOP(TV,ROR,NSSYM(7,NDA),-1)
         END IF
   70 CONTINUE
      WRITE(9,'(A)') 'ENDSYM'

*     WRITE ORFFE INSTRUCTIONS
   80 CONTINUE
         READ(4,'(A)',END=9) INSTR
         IF (INSTR(1:1) .EQ. '}') GO TO 9
*        CORRECT THE WRONG NUMBER OF ATOMS IN THE PARAMETER LIST
*        FORMATS OF ORFFE INSTRUCTIONS WERE CHANGED (3 ==> 5 COLUMNS)
         READ(INSTR,'(2I5)') NFUNC,NSITES
         IF ((NFUNC .EQ. 101 .OR. NFUNC .EQ. 201 .OR. NFUNC .EQ. 207
     &   .OR. NFUNC .EQ. 208 .OR. NFUNC .EQ. 209 .OR. NFUNC .EQ. 210
     &   .OR. NFUNC .EQ. 211) .AND. NSITES .NE. NSITE(NDA)) THEN
            WRITE(9,'(2I5,A)') NFUNC,NSITE(NDA),INSTR(11:)
         ELSE
            WRITE(9,'(A)') INSTR
         END IF
      GO TO 80
    9 END

************************************************************************

      SUBROUTINE WRSYOP(TV,ROR,NCON,LSYM)
*     WRITE SYMMETRY OPERATIONS COMPATIBLE WITH ORFLS AND ORFFE
      REAL TV(3)
      INTEGER ROR(2,3)

      IF (LSYM .NE. 0) CALL TVROR(TV(1),TV(2),TV(3),ROR,LSYM)
      IF (NCON .EQ. 13) THEN
*        BODY CENTERED
         CALL TVROR(TV(1)+0.5,TV(2)+0.5,TV(3)+0.5,ROR,LSYM)
      ELSE IF (NCON .EQ. 4) THEN
*        A-FACE CENTERED
         CALL TVROR(TV(1),TV(2)+0.5,TV(3)+0.5,ROR,LSYM)
      ELSE IF (NCON .EQ. 5) THEN
*        B-FACE CENTERED
         CALL TVROR(TV(1)+0.5,TV(2),TV(3)+0.5,ROR,LSYM)
      ELSE IF (NCON .EQ. 6) THEN
*        C-FACE CENTERED
         CALL TVROR(TV(1)+0.5,TV(2)+0.5,TV(3),ROR,LSYM)
      ELSE IF (NCON .EQ. 7) THEN
*        ALL-FACE CENTERED
         CALL TVROR(TV(1),TV(2)+0.5,TV(3)+0.5,ROR,LSYM)
         CALL TVROR(TV(1)+0.5,TV(2),TV(3)+0.5,ROR,LSYM)
         CALL TVROR(TV(1)+0.5,TV(2)+0.5,TV(3),ROR,LSYM)
      ELSE IF (NCON .EQ. 14 .OR. NCON .EQ. 15) THEN
*        RHOMBOHEDRALLY CENTERED (HEXAGONAL AXES)
         CALL TVROR(TV(1)+0.6666667,TV(2)+0.3333333,TV(3)+0.3333333,ROR,
     &   LSYM)
         CALL TVROR(TV(1)+0.3333333,TV(2)+0.6666667,TV(3)+0.6666667,ROR,
     &   LSYM)
      END IF
      END

************************************************************************

      SUBROUTINE TVROR(TV1,TV2,TV3,ROR,LSYM)
*     MAKE TV SMALL ENOUGH AND WRITE SYMMETRY OPERATION ON A FILE
      INTEGER ROR(2,3)

*     LSYM=-1: INVERTED POSITION (USED IN MADEL.F)

      WRITE(9,'(3(F11.6,2I2))') POSTR(TV1),ROR(1,1),ROR(2,1),
     &   POSTR(TV2),ROR(1,2),ROR(2,2),POSTR(TV3),ROR(1,3),ROR(2,3)
      END

************************************************************************

      REAL FUNCTION POSTR(X)
*     CONVERT A NEGATIVE COMPONENT OF THE TRANSLATION VECTOR INTO A 
*     POSITIVE ONE

      POSTR = X
      IF (X .LT. 0.0) POSTR = POSTR + 1.0
      END

************************************************************************

      INTEGER FUNCTION NUMPAR(NTERMS,LINE1,LABEL,LPAR,NLINE)
*     RETURN A PARAMETER NUMBER FROM LABEL,NUMBER
      INTEGER NTERMS,LPAR(*),NLINE
      CHARACTER LINE1*(*),LABEL(*)*25,FORM(4)*4,LINE*30,TEMP*30
      COMMON /FLAG/ NPRE
      DATA FORM /'(I1)', '(I2)', '(I3)', '(I4)'/

      CALL ADDCOM(LINE1,LENGTH,LINE,ICOM)
      TEMP=LINE
      IF (ICOM .EQ. 0) THEN
*        INTEGER NUMBER
         READ(LINE,FORM(LENGTH)) N
         IF (N .GT. NTERMS .OR. N .LE. 0) THEN
            CALL JOBEND('Parameter number out of range: '//TEMP)
         END IF
         NUMPAR=N
         RETURN
      END IF
      IF (ICOM .EQ. LENGTH .OR. INDEX(LINE(ICOM+1:),',') .NE. 0 .OR.
     &   INDEX(LINE,'(') .NE. 0 .OR.
     &   INDEX(LINE,')') .NE. 0) THEN
         CALL JOBEND('Bad parameter number: '//TEMP)
      END IF

*     A VARIABLE NAME IS SPECIFIED INSTEAD OF A NUMBER
      NUM=0
      IF (LINE(ICOM+1:) .EQ. 'g') THEN
*        ",g" ==> ",1"
         NUM=1
      ELSE IF (LINE(ICOM+1:) .EQ. 'x') THEN
*        ",x" ==> ",2"
         NUM=2
      ELSE IF (LINE(ICOM+1:) .EQ. 'y') THEN
*        ",y" ==> ",3"
         NUM=3
      ELSE IF (LINE(ICOM+1:) .EQ. 'z') THEN
*        ",z" ==> ",4"
         NUM=4
      ELSE IF (LINE(ICOM+1:) .EQ. 'B' .OR. LINE(ICOM+1:) .EQ. 'beta11'
     &.OR. LINE(ICOM+1:) .EQ. 'B11') THEN
*        ",B", ",beta11", OR ",B11' ==> ",5"
         NUM=5
      ELSE IF (LINE(ICOM+1:) .EQ. 'beta22' .OR. LINE(ICOM+1:) .EQ.
     &'B22') THEN
*        ",beta22" OR ",B22" ==> ",6"
         NUM=6
      ELSE IF (LINE(ICOM+1:) .EQ. 'beta33' .OR. LINE(ICOM+1:) .EQ.
     &'B33') THEN
*        ",beta33" OR ",B33" ==> ",7"
         NUM=7
      ELSE IF (LINE(ICOM+1:) .EQ. 'beta12' .OR. LINE(ICOM+1:) .EQ.
     &'B12') THEN
*        ",beta12" OR ",B12" ==> ",8"
         NUM=8
      ELSE IF (LINE(ICOM+1:) .EQ. 'beta13' .OR. LINE(ICOM+1:) .EQ.
     &'B13') THEN
*        ",beta13" OR ",B13" ==> ",9"
         NUM=9
      ELSE IF (LINE(ICOM+1:) .EQ. 'BETA23' .OR. LINE(ICOM+1:) .EQ.
     &'B23') THEN
*        ",beta23" OR ",B23" ==> ",10"
         NUM=10
      END IF

      IF (NUM .EQ. 0) THEN
*        LABEL,NUMBER
         DO 20 I=ICOM+1,LENGTH
            IF (INDEX('0123456789',LINE(I:I)) .EQ. 0) THEN
               CALL JOBEND('A label is followed by a bad number: '//
     &         TEMP)
            END IF
            READ(LINE(ICOM+1:LENGTH),FORM(LENGTH-ICOM)) NUM
            IF (NUM .LE. 0) THEN
               CALL JOBEND('Negative number in (LABEL,NUMBER): '//TEMP)
            END IF
   20    CONTINUE
      END IF

*     THE PARAMETER NUMBER CONTAINS A LABEL
      IPAR=1
      DO 50 I=1,NLINE
*        A LABEL FOR AN ATOMIC SITE ALWAYS CONTAINS '/' AS A SEPARATOR
         LENLBL=INDEX(LABEL(I),'/')-1
         IF (LENLBL .EQ. -1) LENLBL=LENSTR(LABEL(I))
         IF (LENLBL .NE. 0 .AND. LINE(1:ICOM-1) .EQ.
     &   LABEL(I)(1:LENLBL) .AND. NUM .LE. LPAR(I)) THEN
            IF (LPAR(I) .LE. 6 .AND. LINE(ICOM+1:ICOM+4) .EQ. 'BETA')
     &      THEN
               CALL JOBEND('A label is followed by a bad number: '//
     &         TEMP)
            ELSE IF (LPAR(I) .GE. 7 .AND. LINE(LENGTH:LENGTH) .EQ. 'B')
     &      THEN
               CALL JOBEND('A label is followed by a bad number: '//
     &         TEMP)
            END IF
            NUMPAR=IPAR+NUM-1
            RETURN
         ELSE IF (LENLBL .NE. 0 .AND. LINE(1:ICOM-1) .EQ.
     &   LABEL(I)(1:LENLBL) .AND. NUM .GT. LPAR(I)) THEN
            CALL JOBEND('Too large number in (LABEL,NUMBER): '//TEMP)
         END IF
         IPAR=IPAR+LPAR(I)
   50 CONTINUE
*     Dummy statement to avoid a warning messange in Digital Fortran
      NUMPAR = 0
      CALL JOBEND('Such a label has not been input: '//TEMP(1:ICOM-1))
      END

************************************************************************

      SUBROUTINE ADDCOM(LINE1,LENGTH,LINE,ICOM)
*     ADD ",1" IF ONLY A LABEL HAS BEEN INPUT
      CHARACTER LINE1*(*),LINE*30

      LENGTH=LEN(LINE1)
      LINE=LINE1
      ICOM=INDEX(LINE,',')
      IF (ICOM .GT. 0) RETURN
      DO 30 I=1,LENGTH
         IF (INDEX('0123456789',LINE(I:I)) .EQ. 0) GO TO 1
   30 CONTINUE
*     A PARAMETER NUMBER HAS BEEN INPUT
      RETURN
*     ONLY A LABEL HAS BEEN INPUT
    1 LINE(LENGTH+1:LENGTH+2)=',1'
      ICOM=LENGTH+1
      LENGTH=LENGTH+2
      END

************************************************************************

      SUBROUTINE NEWPAR(NLINE,LPAR,NSPLBL,LABEL,A,ID,INP,ICOPY)
*     WRITE LABELS, FINAL PARAMETERS, AND REFINEMENT IDENTIFIERS
*     INPUT
*        NLINE: NUMBER OF LOGICAL PARAMETER LINES
*         LPAR: NUMBERS OF PARAMETERS IN LOGICAL PARAMETER LINES
*       NSPLBL: STARTING POSITION OF LABEL
*        LABEL: LABELS OF LOGICAL PARAMETER LINES (CHARACTER*25)
*            A: FINAL PARAMETERS
*           ID: REFINEMENT IDENTIFIERS
*          INP: FILE NUMBER FOR THE ORIGINAL INPUT FILE
*        ICOPY: FILE NUMBER FOR THE COPY OF THE INPUT FILE

      INTEGER NLINE,LPAR(*),ID(*)
      INTEGER NSPLBL(*)
      REAL A(*)
      CHARACTER LABEL(*)*25, LINE*80

      J=1
      DO I = 1, NLINE
         DO WHILE (.TRUE.)
            READ(ICOPY,'(A)') LINE
            IF (LINE .EQ. ' ') THEN
               WRITE(INP,'(A)') ' '
               CYCLE
            END IF
            DO IS = 1, 80
               IF (LINE(IS:IS) .NE. ' ') EXIT
            END DO
            LENLBL=LENSTR(LABEL(I))
            IF (IS+LENLBL-1 .GT. 80) CALL JOBEND
     &      ('IS + LENLBL - 1 > 80: '//LINE)
            IF (LINE(IS:IS+LENLBL-1) .EQ. LABEL(I)(1:LENLBL)) EXIT
*           WRITE A COMMENT LINE, A SKIPPED LINE, 'If ...', 
*           'else ...', OR 'end if'
            CALL PRLINE(LINE,INP)
         END DO

         CALL WLAID(NSPLBL(I),LABEL(I),A(J),ID(J),LPAR(I),INP)
         J = J + LPAR(I)

         DO WHILE (.TRUE.)
*           SKIP LINES BELONGING TO THE PRESENT LABLE
            READ(ICOPY,'(A)') LINE
            IF (LINE .EQ. ' ') THEN
               WRITE(INP,'(A)') ' '
               EXIT
            END IF
            DO IS = 1, 80
               IF (LINE(IS:IS) .NE. ' ') EXIT
            END DO
            IF (I .EQ. NLINE .AND. LINE(IS:IS) .EQ. '}') THEN
               CALL PRLINE(LINE,INP)
               RETURN
            ELSE IF (INDEX('0123456789.+-',LINE(IS:IS)) .EQ. 0) THEN
*              NEITHER A(I) NOR ID(I) IS PRESENT AT THE TOP OF THIS LINE
               BACKSPACE ICOPY
               EXIT
            END IF
         END DO
      END DO
      END

************************************************************************

      SUBROUTINE WLAID(NSPLBL,LABEL,A,ID,LPAR,INP)
*     WRITE ONE LINE, PACKING AS MANY PARAMETERS AS POSSIBLE IN IT
      INTEGER ID(*),LPAR
      INTEGER NSPLBL
      REAL A(*)
      CHARACTER LABEL*25,BUF*14,LINE*80
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST

*     WRITE LABEL
      LINE=' '
      LENLBL=LENSTR(LABEL)
      LINE(NSPLBL:NSPLBL+LENLBL-1)=LABEL(1:LENLBL)
      IP = NSPLBL+LENLBL+2

      DO I=1,LPAR
         CALL A2NUM(A,I,ISTART,IEND,BUF)
         IF (IP .GE. 71 .OR. IEND-ISTART+1 .GT. 72-IP+1) THEN
            CALL CONTIN(LINE,IP,LENLBL,INP)
         END IF
         LINE(IP:IP+IEND-ISTART)=BUF(ISTART:IEND)
*        TWO SPACES ARE INSERTED BETWEEN TWO PARAMETERS
         IP=IP+IEND-ISTART+3
      END DO
*     LPAR = LENGTH OF AN ID GROUP
      IF (NMODE .NE. 1 .AND. (IP .GT. 72 .OR. LPAR .GT. 72-IP+1)) THEN
         CALL CONTIN(LINE,IP,LENLBL,INP)
      END IF

*     WRITE REFINEMENT IDENTIFIERS 
      WRITE(LINE(IP:IP+LPAR-1),'(80(I1:))') (ID(I),I=1,LPAR)
      CALL PRLINE(LINE,INP)
      END

************************************************************************

      SUBROUTINE A2NUM(A,I,ISTART,IEND,BUF)
*     CONVERT A(I) INTO A PACKED STRING
      REAL A(*)
      CHARACTER BUF*14

      IF (A(I) .LT. 0.1 .OR. A(I) .GE. 10.0E6) THEN
         WRITE(BUF,'(1P,E12.5)') A(I)
      ELSE
         WRITE(BUF,'(1P,G14.6)') A(I)
      END IF
*     DETERMINE THE START AND END POINTS OF A PARAMETER
      DO 20 ISTART=1,14
         IF (BUF(ISTART:ISTART) .NE. ' ') GO TO 1
   20 CONTINUE
    1 IEND=LENSTR(BUF)

      IE=INDEX(BUF,'E')
*     EXPONENTIAL EXPRESSION
      IF (IE .NE. 0) THEN
         CALL EXPO(BUF,IE,ISTART,IEND)
      ELSE
*        DELETE ZERO'S OF THE TAIL
         DO 30 J=IEND,ISTART+1,-1
            IF (BUF(J:J) .NE. '0') GO TO 6
   30    CONTINUE
    6    IEND=J
      END IF
      IF (BUF(IEND:IEND) .EQ. '.' .AND. IEND .LE. 13) THEN
*        ". " ==> ".0"
         IEND=IEND+1
         BUF(IEND:IEND)='0'
      END IF
      END

************************************************************************

      INTEGER FUNCTION NEXT(LINPAR,IP)
*     RETURN THE VALUE OF THE NEXT POSITION OF A PARAMETER OR ID'S
      CHARACTER*80 LINPAR
*     NEXT: = 0, ONLY SPACES FOLLOW IN THIS LINE.
*           > 0, THE POSITION OF THE FOLLOWING PARAMETER.

      DO 50 J=IP,80
         IF (LINPAR(J:J) .NE. ' ') THEN
            NEXT=J
            RETURN
         END IF
   50 CONTINUE
      NEXT=0
      END

************************************************************************

      SUBROUTINE CONTIN(LINE,IP,LENLBL,INP)
*     WRITE THE PRESENT LINE AND PROCEED TO THE NEXT ONE
      INTEGER IP,LENLBL
      CHARACTER LINE*80

      CALL PRLINE(LINE(1:IP-2),INP)
      LINE=' '
      IF (LENLBL .EQ. 0) THEN
         IP=1
      ELSE
         IP=LENLBL+3
      END IF
      END

************************************************************************

      SUBROUTINE EXPO(BUF,IE,ISTART,IEND)
*     COMPRESS AN EXPONENTIAL EXPRESSION
      INTEGER IE,ISTART,IEND
      CHARACTER BUF*14,BUF2*14,CH*1

      DO 10 IS=IE-1,1,-1
         IF ((BUF(IS:IS) .NE. '0') .OR.
     &   (IS .NE. 1 .AND. BUF(IS-1:IS) .EQ. '.0')) THEN
*           DELETE ZERO'S JUST BEFORE 'E'
            BUF2=BUF(IE:14)
            BUF(IS+1:)=BUF2
            IEND=LENSTR(BUF)
            GO TO 2
         END IF
   10 CONTINUE
    2 IF (IEND-ISTART+1 .GT. 4 .AND. BUF(IEND-3:IEND) .EQ. 'E+00') THEN
*        DELETE 'E+00' AT THE TAIL
         IEND=IEND-4
      ELSE IF (IEND-ISTART+1 .GT. 4 .AND. BUF(IEND-3:IEND-1) .EQ. 'E+0')
     &THEN
*        "E+0N" ==> "EN", WHERE N IS AN INTEGER
         CH=BUF(IEND:IEND)
         IEND=IEND-2
         BUF(IEND:IEND)=CH
      ELSE IF (IEND-ISTART+1 .GT. 4 .AND. BUF(IEND-3:IEND-1) .EQ. 'E-0')
     &THEN
*        "E-0N" ==> "E-N", WHERE N IS AN INTEGER
         CH=BUF(IEND:IEND)
         IEND=IEND-1
         BUF(IEND:IEND)=CH
      END IF
*     '0.0' may be output as '0.0E-NN')
      IF (BUF(ISTART:ISTART+3) .EQ. '0.0E') IEND=ISTART+2
      END

************************************************************************

      INTEGER FUNCTION LENSTR(LABEL)
*     RETURN THE ESSENTIAL LENGTH OF A STRING
      CHARACTER*(*) LABEL

      DO 20 I=LEN(LABEL),1,-1
         IF (LABEL(I:I) .NE. ' ') GO TO 1
   20 CONTINUE
    1 LENSTR=I
      END

************************************************************************

      SUBROUTINE DECPAR(NTERMS,A,ID,NSPLBL,LABEL,LPAR,NLINE,NPHL)
*     READ ((LABEL +) PARAMETER + IDENTIFIER) LINES AND DECODE THEM

*     OUTPUT
*       NTERMS: NUMBER OF PARAMETERS (USED TO CHECK THE CORRESPONDING
*               NUMBER IN THE INPUT)
*            A: PARAMETER ARRAY
*           ID: REFINEMENT IDENTIFIER ARRAY
*       NSPLBL: STARTING POSITION OF LABEL
*        LABEL: LABELS OF LOGICAL PARAMETER LINES (CHARACTER*25)
*         LPAR: NUMBERS OF PARAMETERS IN LOGICAL PARAMETER LINES
*        NLINE: NUMBER OF LOGICAL PARAMETER LINES
*         NPHL: NUMBER OF PHYSICAL PARAMETER LINES

      PARAMETER (MAXLEN=320,MAXLAB=800)
      INTEGER NTERMS,ID(*),LPAR(*),NLINE
      INTEGER NSPLBL(*)
      REAL A(*)
      SAVE LINE
      CHARACTER LABEL(*)*25,LINE*320,LINE80*80,NEXTL*80

      NLINE=1
      NPHL=0
*     POINTER FOR PARAMETER NUMBER
      IPAR=1
      DO ILINE=1,1000
         LINE=' '
         LENGTH=0
         DO I=1,MAXLEN/80
            READ(4,'(A)',END=9) LINE80
            NPHL=NPHL+1
            IF (LINE80(1:7) .EQ. 'ENDPARA') THEN
               NLINE=NLINE-1
               NPHL=NPHL-1
               NTERMS=IPAR-1
               RETURN
            END IF

            READ(4,'(A)',END=9) NEXTL
            BACKSPACE 4
            IP=NEXT(NEXTL,1)
            LAB=INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ',NEXTL(IP:IP))
            IF (LAB .EQ. 0 .AND. I .EQ. MAXLEN/80) THEN
               CALL JOBEND('Continuation of lines is limited to 4'//
     &         ' times.  Divide these lines using two or more labels')
            ELSE
               J=LENSTR(LINE80)
               IF (LENGTH+J+1 .GT. MAXLEN) THEN
                  CALL JOBEND('Too long parameter line: '//LINE(1:80)//
     &            '.....')
               END IF
               LINE(LENGTH+1:LENGTH+J+1)=LINE80(1:J)//' '
               LENGTH=LENGTH+J+1
               IF (LAB .NE. 0) THEN
*                 END OF A PARAMETER GROUP
                  LENGTH=LENGTH-1
                  CALL DECOMP(LINE,LENGTH,NLINE,A,ID,IPAR,NSPLBL,LABEL,
     &            LPAR)
                  NLINE = NLINE + 1
                  IF (NLINE .GT. MAXLAB) THEN
                     CALL JOBEND('The total number of labels have '//
     &               'exceeded MAXLAB')
                  END IF
                  EXIT
               END IF
            END IF
         END DO
      END DO

    9 CALL JOBEND('Missing ''}'' to indicate the end of parameters')
      END

************************************************************************

      SUBROUTINE DECOMP(LINE,LENGTH,NLINE,A,ID,IPAR,NSPLBL,LABEL,LPAR)
*     DECOMPOSE THE NLINE'TH LINE INTO A'S AND ID'S
      PARAMETER (MAXLEN=320,MAXPAR=26,NT=999)
      INTEGER LENGTH,NLINE,ID(*),IPAR,LPAR(*),ISTART(26),IEND(26)
*     MAXPAR: MAXIMUM NUMBER OF PARAMETERS IN A LINE
      REAL A(*)
      INTEGER NSPLBL(*)
      SAVE LINE2
      CHARACTER LINE*(*),LINE2*320,LABEL(*)*25,BUFF*120
      COMMON /FLAG/ NPRE
      COMMON /NUML/ NUPDT
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST

      DO I = 1, MAXLEN
         IF (LINE(I:I) .NE. ' ') EXIT
      END DO
      NSPLBL(NLINE) = I
      DO J = I+1, MAXLEN
         IF (LINE(J:J) .EQ. ' ') EXIT
      END DO
      LABEL(NLINE)=LINE(NSPLBL(NLINE):J-1)
      IPOINT=J+1

    4 CONTINUE
*     DETERMINE THE NUMBER OF PARAMETERS CONTAINED IN THIS LINE
      K=1
*     LSTART: FLAG TO INDICATE THAT A DATA IS BEING INPUT
      LSTART=0
      DO J=IPOINT,LENGTH
         IF (J .EQ. 1) THEN
            ISTART(K)=1
            LSTART=1
         ELSE IF (J .EQ. LENGTH) THEN
            IEND(K)=J
            IF (LSTART .EQ. 0) ISTART(K)=J
            K=K+1
         ELSE IF (LINE(J-1:J-1) .EQ. ' ' .AND. LINE(J:J) .NE. ' ' .AND.
     &   K .GT. MAXPAR) THEN
            CALL JOBEND('Too many parameters are contained under '//
     &      'the same label')
         ELSE IF (LINE(J-1:J-1) .EQ. ' ' .AND. LINE(J:J) .NE. ' ') THEN
            ISTART(K)=J
            LSTART=1
         ELSE IF (LINE(J-1:J-1) .NE. ' ' .AND. LINE(J:J) .EQ. ' ') THEN
            IEND(K)=J-1
            LSTART=0
            K=K+1
         END IF
      END DO

      IF (NMODE .NE. 1) THEN
*        A SERIES OF ID'S IS NOT COUNTED
         LPAR(NLINE)=K-2
      ELSE 
*        ID(I)'S ARE NEGLECTED (ADDED BY KUMAZAWA)
         LPAR(NLINE) = K - 2
         N1 = ISTART(LPAR(NLINE)+1)
         N2 = IEND(LPAR(NLINE)+1)
         IF (INDEX(LINE(N1:N2),'.') .NE. 0) 
     &   LPAR(NLINE) = LPAR(NLINE) + 1
      END IF

*     CHECK THE COINCIDENCE OF THE NUMBER OF PARAMETERS AND THE NUMBER
*     OF ID'S
      IF (NMODE .NE. 1) THEN
         N1=ISTART(LPAR(NLINE)+1)
         N2=IEND(LPAR(NLINE)+1)
         IF (LPAR(NLINE) .GE. 5 .AND. LPAR(NLINE) .NE. N2-N1-4 .AND.
     &   LINE(ISTART(5):ISTART(5)) .EQ. '+') THEN
            WRITE(6,'(//11X,A,I2,A)') 'The numbers of parameters and '//
     &      'ID''s in label #',NLINE,' do not coincide (+B type)'
            CALL JOBEND(' ')
         ELSE IF (LPAR(NLINE) .GE. 5 .AND. LPAR(NLINE) .EQ. N2-N1-4
     &   .AND. LINE(ISTART(5):ISTART(5)) .EQ. '+') THEN
*           THE 5'TH PARAMETER WHOSE FIRST CHARACTER IS '+' IS REGARDED
*           AS AN EQUIVALENT ISOTROPIC DISPLACEMENT PARAMETER, 
*           WHICH WILL BE CONVERTED INTO SIX ANISOTROPIC DISPLACEMENT 
*           PARAMETERS LATER. '+' AT THE TOP OF THE 5'TH PARAMETER IS
*           DELETED, TEMPORARILY.  BETA22, BETA33, BETA12, BETA13, AND
*           BETA23 (ALL ZERO) ARE INSERTED, AND THEN THE NEW LINE IS 
*           DECODED AGAIN.
            IF (NUPDT .EQ. 2) THEN
               CALL JOBEND('+B should not be entered if NUPDT = 2')
            END IF
            LINE2 = LINE(1:ISTART(5)-1) // LINE(ISTART(5)+1:IEND(5)) //
     &              ' 0 0 0 0 0' // LINE(IEND(5)+1:)
            LENGTH=LENGTH+9
            LINE=LINE2
            GO TO 4
         ELSE IF (LPAR(NLINE) .NE. N2-N1+1) THEN
            WRITE(6,'(//11X,A,I2,A)') 'The numbers of parameters and '//
     &      'ID''s in label #',NLINE,' do not coincide'
            CALL JOBEND(' ')
         END IF
         DO J=N1,N2
            IF (INDEX('0123',LINE(J:J)) .EQ. 0) THEN
               BUFF='Bad ID values: '//LINE(N1:N2)
               CALL JOBEND(BUFF)
            END IF
         END DO
      END IF

*     CONVERT A LINE INTO A'S AND ID'S
      DO I=1,LPAR(NLINE)
         DO J=ISTART(I),IEND(I)
            IF (INDEX('+-0123456789E.',LINE(J:J)) .EQ. 0) THEN
               BUFF = 'A bad parameter '//LINE(ISTART(I):IEND(I))//
     &                ' in label '//LABEL(NLINE)
               CALL JOBEND(BUFF)
            END IF
         END DO

*        GET THE IPAR'TH PARAMETER
         CALL GETPAR(LINE(ISTART(I):IEND(I)),A,IPAR)
         IF (NMODE .NE. 1) THEN
*           GET THE IPAR'TH IDENTIFIER
            IP=ISTART(LPAR(NLINE)+1)+I-1
            READ(LINE(IP:IP),'(I1)') ID(IPAR)
            IF (ID(IPAR) .GT. 3) THEN
               CALL JOBEND('ID(I) must be less or equal to 3')
            END IF
         END IF
         IPAR=IPAR+1
         IF (IPAR .GT. NT) THEN
            CALL JOBEND('The number of parameters has exceeded '//
     &      'the maximum value')
         END IF
      END DO
      END

************************************************************************

      SUBROUTINE GETPAR(LINE,A,IPAR)
*     GET A REAL VALUE, A(IPAR)
      INTEGER IPAR
      REAL A(*)
      CHARACTER LINE*(*),FORM(20)*7
      DATA FORM /'(F1.0)', '(F2.0)', '(F3.0)', '(F4.0)', '(F5.0)',
     &  '(F6.0)', '(F7.0)', '(F8.0)', '(F9.0)', '(F10.0)', '(F11.0)',
     &  '(F12.0)', '(F13.0)', '(F14.0)', '(F15.0)', '(F16.0)',
     &  '(F17.0)', '(F18.0)', '(F19.0)', '(F20.0)'/

      LENGTH=LEN(LINE)
      IF (LENGTH .GT. 20) THEN
         CALL JOBEND('A too long parameter has been input')
      END IF
      READ(LINE,FORM(LENGTH)) A(IPAR)
      END

************************************************************************

      SUBROUTINE PARNUM(NTERMS,ID,LABEL,LPAR,NLINE,IPAR,NPAR,NCHNG)
*     READ REFINABLE-PARAMETER NUMBERS IN INCREMENTAL REFINEMENTS
*     (LABELS ARE USED)
      PARAMETER (NR=400)
*Rev 1.1 2003.05.16 Izumi
      INTEGER NTERMS,ID(*),LPAR(*),NLINE,IPAR(NR,50),NPAR(*),NCHNG
      CHARACTER LABEL(*)*25,LINE*80

*     INPUT
*        NTERMS: NUMBER OF PARAMETERS
*            ID: REFINEMENT IDENTIFIERS
*         LABEL: LABELS FOR LOGICAL PARAMETER LINES
*          LPAR: NUMBERS OF PARAMETERS IN LOGICAL PARAMETER LINES
*         NLINE: NUMBER OF LOGICAL PARAMETER LINES
*     OUTPUT
*          IPAR: NUMBERS FOR RFINABLE PARAMETERS IN EACH CYCLE
*          NPAR: NUMBER OF REFINABLE PARAMETERS IN EACH CYCLE
*         NCHNG: NUMBER OF CYCLES WHERE PARAMETERS ARE REFINED
*                INCREMENTALLY

*     ILOG: LOGICAL LINE NUMBER (NUMBER FOR COMBINED LINES)
      ILOG=1
*     IREF: PRESENT TOTAL NUMBER OF REFINABLE PAREMETERS IN THIS
*           LOGICAL LINE
      IREF=1

      WRITE(6,'(//11X,A)') 'Lines to give variable parameters in '//
     &'iterative cycles'
*     ILINE: PHYSICAL LINE NUMBER
      DO 50 ILINE=1,1000
         READ(4,'(A80)',END=9) LINE

         IF (IREF .EQ. 1 .AND. LINE(1:1) .EQ. '}') THEN
*           END OF REFINABLE PARAMETERS: '}'
            NCHNG=ILOG-1
            RETURN
         END IF
         WRITE(6,'(13X,A)') LINE

         LENGTH=LENSTR(LINE)
*        DETERMINE THE STARTING POINT
         DO I = 1, LENGTH
            IF (LINE(I:I) .NE. ' ') EXIT
         END DO
         ISTART=I

*Rev 1.1 2003.05.16 Izumi {
         IF (ILOG .EQ. 51) CALL JOBEND('Refinable parameters in '//
     &   'each cycle may be specified up to 50 cycles')
* }
         DO J = I, LENGTH
            IF (ISTART .EQ. 0 .AND. LINE(J:J) .EQ. '/') THEN
*              END OF A LOGICAL LINE (' /' HAS BEEN DETECTED)
               EXIT
            ELSE IF (ISTART .GT. 0 .AND. J .EQ. LENGTH .AND.
     &      LINE(J:J) .NE. '/') THEN
*              THE LAST PARAMETER NUMBER OF THIS PHYSICAL LINE
               IPAR(IREF,ILOG)=NUMPAR(NTERMS,LINE(ISTART:J),
     &         LABEL,LPAR,NLINE)
               CALL CHKID(IPAR(IREF,ILOG),ID,NTERMS,IREF)
*              CONTINUE TO THE NEXT PHYSICAL LINE
               GO TO 50
            ELSE IF (J .EQ. LENGTH .AND. ISTART .EQ. 0) THEN
*              ONLY ONE NUMERIC CHARACTER
               IPAR(IREF,ILOG)=NUMPAR(NTERMS,LINE(J:J),LABEL,LPAR,NLINE)
               CALL CHKID(IPAR(IREF,ILOG),ID,NTERMS,IREF)
               EXIT
            ELSE IF (ISTART .GT. 0 .AND. (LINE(J:J) .EQ. ' ' .OR.
     &      LINE(J:J) .EQ. '/')) THEN
*              END OF A PARAMETER NUMBER
               IPAR(IREF,ILOG)=
     &         NUMPAR(NTERMS,LINE(ISTART:J-1),LABEL,LPAR,NLINE)
               CALL CHKID(IPAR(IREF,ILOG),ID,NTERMS,IREF)
               IF (LINE(J:J) .EQ. '/') EXIT
               ISTART=0
            ELSE IF (LINE(J:J) .NE. ' ' .AND. ISTART .EQ. 0) THEN
*              START OF A PARAMETER NUMBER
               ISTART=J
            END IF
         END DO

         NPAR(ILOG)=IREF-1
         ILOG=ILOG+1
         IREF=1
   50 CONTINUE
    9 CALL JOBEND('Missing ''}'' to indicate the end of reference '//
     &'numbers for variable parameters')
      END

************************************************************************

      SUBROUTINE CHKID(IPAR,ID,NTERMS,IREF)
*     CHECK WHETHER OR NOT IPAR IS A NUMBER FOR A REFINABLE PARAMETER
      INTEGER IPAR,ID(*),NTERMS,IREF
      COMMON /FLAG/ NPRE

      IF (IPAR .GT. NTERMS .OR. IPAR .LE. 0) THEN
         WRITE(6,100) IPAR
         CALL JOBEND(' ')
      ELSE IF (ID(IPAR) .NE. 1) THEN
         WRITE(6,200) IPAR
         CALL JOBEND(' ')
      END IF
      IREF=IREF+1
  100 FORMAT(//11X,'Bad refinable parameter number: ',I3)
  200 FORMAT(//11X,'A constant or constrained parameter was regarded ',
     &'as a variable one.'/11X,'Parameter number: ',I3)
      END

************************************************************************

      SUBROUTINE SIMDAT(YFIT,NSTEP,TITLE,DEG,STEP,NREFL)
*     RECORD SIMULATION DATA
      PARAMETER (NB=7000,NPH=8)
      INTEGER H
      REAL YFIT(*),DEG(*)
      CHARACTER*80 TITLE
      CHARACTER*50 FILE20
      INTEGER WIDTH,HEIGHT,YMIN,YMAX,LBGPLOT,LDEL,OFFSETD,MSIZE,FSIZE,
     &  LSIZE,OFFSET,INDREF
      COMMON /IGOR/ WIDTH,HEIGHT,YMIN,YMAX,LBGPLOT,LDEL,OFFSETD,MSIZE,
     &  FSIZE,LSIZE,OFFSET(NPH),INDREF
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
*Rev 1.0f 2000.12.11 Izumi {
      COMMON /PROFI/ PCOT(NB),SHFT(NB),NTSIM(NB),FNORM(NB),TNTSIM(NB),
     &  FWHM(NB),ETA(NB),SQSG(NB),CTH(NB),SIGP(NB),DFWDG(NB),DEDFFG(NB),
     &  DFWDL(NB),DEDFFL(NB),SIGPAP(NB)
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      COMMON /ZSPFN/ PEAKSH(NB),SQH2(NB),SQH2L(NB),SQH2G(NB),ASYM(NB),
     &ETAL(NB),ETAH(NB),FNORM1(2,NB),COMLG(2,NB),COML(2,NB),COMG(2,NB),
     &COEFL(2,NB),COEFG(2,NB),DECAY(2,NB),COMPEAR(2,NB),CDX(2,NB),
     &CDA1(2,NB),CDA2(2,NB),CDW(2,NB),CDRL1(2,NB),CDRL2(2,NB),
     &CDRH1(2,NB),DMPART(2,NB)
* }
      DATA RD2/57.29578/

C     NORMALIZATION OF INTENSITIES
      YMAX=YFIT(1)
      DO 15 J=2,NSTEP
         IF(YFIT(J).GT.YMAX) YMAX=YFIT(J)
   15 CONTINUE
      IF (NPAT .EQ. 2) THEN
         FULLSC=100000.0
      ELSE
         FULLSC=100.0
      END IF
      DO 25 J=1,NSTEP
         YFIT(J)=YFIT(J)*FULLSC/YMAX
   25 CONTINUE
      CALL OPENFILE(4,FILE20)

*     NPAT = 2: Create a Macplot/RietPlot file storing the simulated pattern.
*     NPAT = 3: Not implemented.
*     NPAT = 4: Create a SigmaPlot file storing the simulated pattern.
*     NPAT = 5: Create an Igor Pro file storing the simulated pattern.
      IF (NPAT .EQ. 2) THEN
         REWIND 20
         WRITE(20,'(A80)') TITLE
         WRITE(20,'(I5,F10.3,F10.7,I5)') NSTEP,DEG(1)*RD2,
     &   STEP*RD2,NREFL
         WRITE(20,'(20I7)') (NINT(YFIT(J)),J=1,NSTEP)
         WRITE(20,'(17F8.3)') (PEAK(J)*RD2,J=1,NREFL)
      ELSE IF (NPAT .EQ. 3) THEN
         CALL JOBEND('Simulation is not possible with PLOT')
      ELSE IF (NPAT .EQ. 4) THEN
* Rev 1.02 2001.06.27 Izumi {
*        Gnuplot file (only one phase)
*        Length of tick marks = 100.0/62.0 
*        y coordinate of tick marks = -100/28.0       
         DO J = 1, NSTEP   
            IF (J .LE. NREFL) THEN
*              2-theta, calculated intensity, hkl, peak positions (shifted
*              and non-shifted), and d
               WRITE(20,'(1P,E10.5E1,0P,F9.4,3I4,2F8.3,F8.4)') 
     &         DEG(J)*RD2, YFIT(J), H(J), K(J), L(J), SHPOS(J),
     &         PEAK(J)*RD2, D(J)
            ELSE
*              2-theta and calculated intensity
               WRITE(20,'(1P,E10.5E1,0P,F9.4)') DEG(J)*RD2, YFIT(J)
* }
            END IF
         END DO
      ELSE IF (NPAT .EQ. 5) THEN
*        Igor Pro file (only one phase)
         WRITE(20,'(A)') 'IGOR'
         WRITE(20,'(A)') 'WAVES/O twoth, ycal'
         WRITE(20,'(A)') 'BEGIN'
         DO J=1,NSTEP
            WRITE(20,'(F7.3,F8.3)') DEG(J)*RD2, YFIT(J)
         END DO
         WRITE(20,'(A)') 'END'
* Rev 1.02 2001.06.26 Izumi 
         WRITE(20,'(A)') 'WAVES/O xphase_1, x0phase_1, yphase_1, d_1'
         WRITE(20,'(A)') 'BEGIN'         
* Rev 1.0r 2001.02.15 Izumi {
         DO J = 1, NREFL      
* Rev 1.02 2001.06.27 Izumi {
            WRITE(20,'(F7.3,F8.3,A,F7.4)') SHPOS(J), PEAK(J)*RD2, 
     &      ' 0 ', D(J)
* }
         END DO         
* }
         WRITE(20,'(A)') 'END'
* Rev 1.0v 2001.03.06 Izumi
         WRITE(20,'(3(A,I1))') 'WAVES/B/O h_1, k_1, l_1'
         WRITE(20,'(A)') 'BEGIN'
         DO J=1,NREFL
            WRITE(20,'(I3,2I4)') H(J),K(J),L(J)
         END DO 
         WRITE(20,'(A)') 'END'
*        COMMANDS IN AN IGOR TEXT FILE
         WRITE(20,'(A)') 'X Display ycal vs twoth'
         WRITE(20,'(2(A,F7.2))') 'X SetAxis bottom', DEG(1)*RD2, ',',
     &   DEG(NSTEP)*RD2
         WRITE(20,'(A)') 'X SetAxis left -10,105'
* Rev 1.0v 2001.03.06 Izumi {
         WRITE(20,'(A)') 'X AppendToGraph yphase_1 vs xphase_1'
         WRITE(20,'(2A)') 'X ModifyGraph lsize(ycal)=0.7,rgb(ycal)=',
     &                    '(1,16019,65535)'
*Rev 1.07 2002.08.22 Izumi (indicated by Dr. Fujimori) {
         WRITE(20,'(2A,I1,A)') 
     &   'X ModifyGraph offset(yphase_1)={0,-3.57},',
     &   'mode(yphase_1)=3,marker(yphase_1)=10,msize(yphase_1)=',MSIZE,
     &   ',mrkThick(yphase_1)=0.6,rgb(yphase_1)=(3,52428,1)'
* }
         WRITE(20,'(A,I2,2(A,I4),A)') 
     &   'X ModifyGraph tick=2,mirror=1,font="Times",fSize=', FSIZE, 
     &   ',btLen=5,width=', WIDTH, ',height=', HEIGHT, ',standoff=0'
         WRITE(20,'(2A)') 'X ModifyGraph width=0,height=0,',
     &   'margin(right)=20'
         WRITE(20,'(2A)') 'X ModifyGraph stLen(bottom)=3,',
     &   'nticks(bottom)=6,minor(bottom)=1,sep(bottom)=20'
         WRITE(20,'(A,I2,A)') 'X Label left "\\F''Times''\\Z', LSIZE,
     &   'Intensity"'
         WRITE(20,'(A,I2,A,I2,A)') 
     &   'X Label bottom "\\F''Times''\\Z', LSIZE,
     &   '2\\F''Symbol''\\f02q\\f00\\F''Times''  / \\F''Symbol''°"'
      END IF
      CLOSE(UNIT=20)
      END

************************************************************************

      SUBROUTINE OUT20(DEG,XINT,YFIT,NPTS,STEP,NREFL,TITLE,A,NEXC)

*     RECORD THE RESULTS OF RIETVELD ANALYSIS IN FILE #20

*     THE FOLLOWING VARIABLES AND ARRAYS CAN BE USED TO CREATE FILE #20
*     WITH ANY USER-DESIRED FORMAT

*     (1) VARIABLES
*           NPTS: NUMBER OF INTENSITY DATA (SOME REGIONS MAY BE SKIPPED)
*         NPHASE: NUMBER OF PHASES
*          NREFL: NUMBER OF BRAGG REFLECTIONS
*           STEP: STEP WIDTH IN RADIAN

*     (2) ARRAYS
*            DEG: 2-THETA'S IN RADIAN
*           XINT: OBSERVED INTENSITIES
*           YFIT: CALCULATED INTENSITIES
*             BG: BACKGROUND INTENSITIES
* V1.1e, 1999.02.15, Izumi 
*         DEGNOR: NORMALIZED 2-THETA
*           PEAK: BRAGG POSITIONS IN RADIAN
*           NOPH: PHASE NUMBER FOR PEAKS

      PARAMETER (NB=7000,NPH=8,NAP=150,NT=999,NP=80000)
      REAL DEG(*),XINT(*),YFIT(*),SHIFT(NPH),LENTIC,A(*)
      INTEGER H,ISTART(100),IEND(100),LREFL(NB,NPH),IC(NPH)
      CHARACTER TITLE*80,PARNAM*60,PHNAME*25,FILE20*50
* Rev 1.02 2001.06.26 Izumi {
      CHARACTER LINE*334
      INTEGER WIDTH,HEIGHT,YMIN,YMAX,LBGPLOT,LDEL,OFFSETD,MSIZE,FSIZE,
     &  LSIZE,OFFSET,INDREF
      COMMON /IGOR/ WIDTH,HEIGHT,YMIN,YMAX,LBGPLOT,LDEL,OFFSETD,MSIZE,
     &  FSIZE,LSIZE,OFFSET(NPH),INDREF
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /C/ BG(NP),BGINC(NP),DEGNOR(NP)
      COMMON /DETECT/ VARINV(NP)
      COMMON /EXCREG/ DEGEXC(2,50),NSKIP
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /PAR/ PARNAM(NT),PHNAME(NPH)
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
*Rev 1.0f 2000.12.11 Izumi {
      COMMON /PROFI/ PCOT(NB),SHFT(NB),NTSIM(NB),FNORM(NB),TNTSIM(NB),
     &  FWHM(NB),ETA(NB),SQSG(NB),CTH(NB),SIGP(NB),DFWDG(NB),DEDFFG(NB),
     &  DFWDL(NB),DEDFFL(NB),SIGPAP(NB)
* }
      COMMON /UNCONS/ NUCREF,IPHASE(90),IHKL(3,90),NOPPP(90),LPPP(NB),
     &  RCUT(2,90),NCUT
*Rev 1.0f 2000.12.11 Izumi {
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      COMMON /ZSPFN/ PEAKSH(NB),SQH2(NB),SQH2L(NB),SQH2G(NB),ASYM(NB),
     &ETAL(NB),ETAH(NB),FNORM1(2,NB),COMLG(2,NB),COML(2,NB),COMG(2,NB),
     &COEFL(2,NB),COEFG(2,NB),DECAY(2,NB),COMPEAR(2,NB),CDX(2,NB),
     &CDA1(2,NB),CDA2(2,NB),CDW(2,NB),CDRL1(2,NB),CDRL2(2,NB),
     &CDRH1(2,NB),DMPART(2,NB)
* }
      DATA RD,RD2/0.01745329,57.29578/

      IF (NPAT .LE. 1 .OR. NPAT .GT. 5) THEN
         CALL JOBEND('An invalid NPAT value')
      END IF
      
      IF (NPAT .NE. 5) THEN
         CALL OPENFILE(4,FILE20)
      END IF

      IF (NPAT .EQ. 2) THEN
*        Macplot/RietPlot FILE
*        TWO OR MORE BLOCKS OF INTENSITY DATA CAN BE RECORDED BY EXCLUDING
*        USER-SPECIFIED 2-THETA REGIONS
         NBLOCK=1
         ISTART(1)=1
         DO I = 1, 100
            DO J=ISTART(NBLOCK)+1,NPTS
               IF (J .GE. NPTS) THEN
                  IEND(NBLOCK)=MIN(J,NPTS)
                  GO TO 1
               ELSE IF (DEG(J)-DEG(J-1) .GT. 2.0*STEP) THEN
                  IEND(NBLOCK)=J-1
                  NBLOCK=NBLOCK+1
                  ISTART(NBLOCK)=J
                  EXIT
               END IF
            END DO
         END DO

    1    WRITE(20,'(A)') TITLE
C        NBLOCK: NUMBER OF BLOCKS IN WHICH THE STEP WIDTH IS FIXED
         WRITE(20,'(I6,F10.4,2I6)') NBLOCK,STEP*RD2,NREFL,NPHASE
         DO I=1,NBLOCK
C           WRITE NUMBER OF DATA AND STARTING 2-THETA IN EACH BLOCK
            WRITE(20,'(I5,F10.3)') IEND(I)-ISTART(I)+1,
     &      DEG(ISTART(I))*RD2
            WRITE(20,'(7(2I7,I6:))') (NINT(XINT(J)),NINT(YFIT(J)),
     &      NINT(BG(J)),J=ISTART(I),IEND(I))
         END DO
         WRITE(20,'(17F8.3)') (PEAK(J)*RD2,J=1,NREFL)
         WRITE(20,'(140I1)') (NOPH(J),J=1,NREFL)
      END IF
      
      IF (NPAT .EQ. 3) THEN
*        DMPLOT FILE
*        NO DATA ARE SKIPPED IN THIS CASE; THEREFORE, OBSERVED AND
*        CALCULATED INTENSITIES FOR EXCLUDED REGIONS ARE INCLUDED.
         REWIND 7
         READ(7) NTOTAL,START,STEPW
         READ(7) (XINT(J),J=1,NTOTAL)
         CLOSE(UNIT=7)

         DEG1=DEG(1)
         DEG2=DEG(NPTS)
         FT1=2.0/(DEG2-DEG1)
         FT2=-0.5*(DEG1+DEG2)
*        CALCULATE INTENSITIES BY INCLUDING REGIONS WHICH HAVE BEEN
*        EXCLUDED BY USERS
         JJ = 1
         DO J=1,NTOTAL
            DEG(JJ)=(START+FLOAT(J-1)*STEPW)*RD
            IF (DEG(JJ) .GT. DEG2) EXIT
            IF (DEG(JJ) .GE. DEG1) THEN
               DEGNOR(JJ) = FT1*(DEG(JJ) + FT2)
               YFIT(JJ)=CALINT(DEG,JJ,A)
               XINT(JJ) = XINT(J)
               JJ = JJ + 1
            END IF
         END DO
         NTOTAL = JJ - 1

         WRITE(20,'(A)') TITLE(1:70)
         DO I=1,NPHASE
            IC(I)=0
            DO J=1,NREFL
               IF (NOPH(J) .EQ. I) IC(I)=IC(I)+1
            END DO
         END DO
         WRITE(20,150) NPHASE, (IC(IIPHAS),IIPHAS=1,NPHASE)
  150    FORMAT (1X, 'NO. OF PHASESZ', 2X,I4,/,1X,
     *   ' NO. OF REFLECTIONS IN EACH PHASEQ',2X,8I4)
         WRITE(20,'(1X,A)') 'BRAGG POSITIONSZ'
         DO I=1,NPHASE
            WRITE (20,'(/1X,A,I2,3X,A)') 'PHASE',I,PHNAME(I)
            DO J=1,NREFL
               IF (NOPH(J) .EQ. I) WRITE(20,'(1X,F8.3)') PEAK(J)*RD2
            END DO
         END DO
         WRITE(20,190) NTOTAL, DEG(1)*RD2, STEPW
  190    FORMAT (1X,'NPTSZ',I5,/, 1X, 'THMINZ',F8.3,/,1X,'STEPZ',F8.3,/,
     *   1X,4HYOBS,4X,6HYCALCZ)
         DO J=1,NTOTAL
            WRITE(20,'(1X,2F8.0)') XINT(J),YFIT(J)
         END DO
      END IF
      
      IF (NPAT .EQ. 4 .OR. NPAT .EQ. 5) THEN
*        PREPROCESSING TO CREATE A SigmaPlot/Igor Pro FILE
         DO J = 1, NPHASE
            DO I = 1, NB
               LREFL(I,J)=0
            END DO
            II = 0
            DO I = 1, NREFL
               IF (NOPH(I) .EQ. J) THEN
                  II = II + 1
                  LREFL(II,J) = I
               END IF
            END DO
         END DO
      END IF

* Rev 1.02 2001.06.28 Izumi {
      IF (NPAT .EQ. 4) THEN
*        Gnuplot file
         DO J = 1, NPTS  
            LINE = ' ' 
            IF (J .GE. 2) THEN
               CENDEG = RD2*( DEG(J) + DEG(J-1) )/2.0
               IF (NEXC .EQ. 1) THEN
                  DO I=1,NSKIP
                     IF (CENDEG .GE. DEGEXC(1,I).AND.
     &               CENDEG .LE. DEGEXC(2,I)) THEN
                        GO TO 2
                     END IF
                  END DO
               END IF
            END IF
*           2-theta, Observed, calculated and background intensities
            WRITE(LINE,'(1P,E10.5E1,3E12.5E1)') DEG(J)*RD2,XINT(J),
     &      YFIT(J),BG(J)

    2       IF (J .LE. NREFL) THEN
*              Peak positions, hkl, and d
               DO IPH = 1, NPHASE
                  IS = 47 + (IPH - 1)*36
                  IF (LREFL(J,IPH) .GT. 0) THEN
                     WRITE(LINE(IS:IS+35),'(3I4,2F8.3,F8.4)') 
     &               H(LREFL(J,IPH)), K(LREFL(J,IPH)), L(LREFL(J,IPH)),
     &               SHPOS(LREFL(J,IPH)), PEAK(LREFL(J,IPH))*RD2, 
     &               D(LREFL(J,IPH))
                  END IF
               END DO
               DO IJ = 334, 2, -1
                  IF (LINE(IJ:IJ) .NE. ' ') EXIT
               END DO
            END IF
            WRITE(20,'(A)') LINE(1:IJ)
         END DO
      END IF
* }
      
      IF (NPAT .EQ. 5) THEN
*        Igor Pro file
         WRITE(20,'(A)') 'WAVES/O twoth, yobs, ycal, delta, bg'
         WRITE(20,'(A)') 'BEGIN'
         DO J = 1, NPTS
            IF (J .GE. 2) THEN
               CENDEG = RD2*( DEG(J) + DEG(J-1) )/2.0
               IF (NEXC .EQ. 1) THEN
                  DO I=1,NSKIP
                     IF (CENDEG.GE.DEGEXC(1,I).AND.
     &                   CENDEG.LE.DEGEXC(2,I))
     &                   WRITE(20,'(A)') 'NaN NaN NaN NaN NaN'  
                  END DO
               END IF
            END IF
            SELECT CASE (LDEL)
               CASE (0)
                  DELTAJ = XINT(J)-YFIT(J)
               CASE (1,2)
                  IF (VARINV(1) .EQ. 0.0) THEN
                     WEIGHTJ = 1.0/XINT(J)
                     SIGMAJ = SQRT(XINT(J))
                  ELSE
                     WEIGHTJ = 1.0/VARINV(J)
                     SIGMAJ = SQRT(VARINV(J))
                  END IF
                  IF (LDEL .EQ. 1) THEN
                     DELTAJ = SQRT(WEIGHTJ)*(XINT(J)-YFIT(J))
                  ELSE
*                    Refer to Eq. (1.13) in "The Rietveld Method"
                     DELTAJ = ((XINT(J)-YFIT(J))/XINT(J))/SIGMAJ
                  END IF
            END SELECT
            WRITE(20,'(1P,5E12.5E1)') DEG(J)*RD2,XINT(J),YFIT(J),
     &      DELTAJ,BG(J)
         END DO
         WRITE(20,'(A)') 'END'

*        X AND Y COORDINATES FOR TICK MARKS FOR BRAGG POSITIONS
         DO I = 1, NPHASE
* Rev 1.02 2001.06.26 Izumi {
            WRITE(20,'(4(A,I1))') 'WAVES/O xphase_', I, ', x0phase_', I,
     &      ', yphase_', I, ', d_',I
* }
            WRITE(20,'(A)') 'BEGIN'
            J = 1   
            DO WHILE (LREFL(J,I) .GT. 0)
*Rev 1.0f 2000.12.11 Izumi {
*              Tick marks are shifted (requested by T. Ikeda)
C              WRITE(20,'(F7.3,A,F7.4)')
C    &         PEAK(LREFL(J,I))*RD2,' 0 ',D(LREFL(J,I))
* Rev 1.02 2001.06.26 Izumi {
               WRITE(20,'(F7.3,F8.3,A,F7.4)') SHPOS(LREFL(J,I)), 
     &         PEAK(LREFL(J,I))*RD2, ' 0 ', D(LREFL(J,I))
* }
               J = J + 1
            END DO 
            WRITE(20,'(A)') 'END'
            WRITE(20,'(3(A,I1))') 'WAVES/B/O h_', I, ', k_', I,', l_',I
            WRITE(20,'(A)') 'BEGIN'
            J = 1   
            DO WHILE (LREFL(J,I) .GT. 0)
               WRITE(20,'(I3,2I4)') H(LREFL(J,I)),K(LREFL(J,I)),
     &         L(LREFL(J,I))
               J = J + 1
            END DO 
            WRITE(20,'(A)') 'END'
* Rev 1.02 2001.06.26 Izumi {
            WRITE(20,'(3(A,I1))') 'WAVES/O xppp_', I, ', x0ppp_', I,
     &      ', yppp_', I
* }
            WRITE(20,'(A)') 'BEGIN'
            J = 1   
            DO WHILE (LREFL(J,I) .GT. 0)
               IF (LPPP(LREFL(J,I)).GT.0) THEN
* Rev 1.02 2001.06.27 Izumi {
                  WRITE(20,'(F7.3,F8.3,A)') SHPOS(LREFL(J,I)), 
     &            PEAK(LREFL(J,I))*RD2, ' 0'
* }
               END IF
               J = J + 1
            END DO 
            WRITE(20,'(A)') 'END'
         END DO
*        COMMANDS IN AN IGOR TEXT FILE
*        Deleted on 2000.7.27 to avoid an error in riplx
*        WRITE(20,'(A)')  'X Edit xphase_1,d_1,h_1,k_1,l_1'
         WRITE(20,'(A)') 'X Display ycal,delta vs twoth'
         IF (LBGPLOT .EQ. 1) THEN
            WRITE(20,'(A)') 'X AppendToGraph bg vs twoth'
            WRITE(20,'(2A)') 'X ModifyGraph rgb(bg)=(0,0,0),',
     &      'lsize(bg)=0.7'
         END IF
         WRITE(20,'(A)') 'X AppendToGraph yobs vs twoth'
         WRITE(20,'(2(A,F7.2))') 'X SetAxis bottom', DEG(1)*RD2, ',',
     &   DEG(NPTS)*RD2
         IF (YMIN .NE. 0 .AND. YMAX .NE. 0)
     &   WRITE(20,'(2(A,I8))') 'X SetAxis left', YMIN, ',', YMAX
         DO I = 1, NPHASE
            WRITE(20,'(2(A,I1))') 'X AppendToGraph yphase_', I, 
     &      ' vs xphase_',I
            WRITE(20,'(2(A,I1))') 'X AppendToGraph yppp_', I, 
     &      ' vs xppp_',I
         END DO
         WRITE(20,'(A)') 
     &   'X ModifyGraph lsize(ycal)=0.7,rgb(ycal)=(1,26214,26214)'
         WRITE(20,'(2A)') 'X ModifyGraph mode(yobs)=3,msize(yobs)=1,',
     &   'mrkThick(yobs)=0.4,rgb(yobs)=(39321,1,1)'
         WRITE(20,'(A,I6,A)') 'X ModifyGraph offset(delta)={0,', 
     &   OFFSETD, '},lsize(delta)=0.7,rgb(delta)=(1,16019,65535)'
         DO I = 1, NPHASE
            WRITE(20,'(A,I1,A,I6,6(A,I1),A)')
     &      'X ModifyGraph offset(yphase_', I, ')={0,', OFFSET(I), 
     &      '},mode(yphase_', I, ')=3,marker(yphase_', I, 
     &      ')=10,msize(yphase_', I, ')=', MSIZE,',mrkThick(yphase_', I,
     &      ')=0.6,rgb(yphase_', I, ')=(1,39321,19939)'
            WRITE(20,'(A,I1,A,I6,6(A,I1),A)')
     &      'X ModifyGraph offset(yppp_', I, ')={0,', OFFSET(I), 
     &      '},mode(yppp_', I, ')=3,marker(yppp_', I, 
     &      ')=10,msize(yppp_', I, ')=', MSIZE,',mrkThick(yppp_', I,
     &      ')=0.6,rgb(yppp_', I, ')=(29524,1,58982)'
         END DO
         WRITE(20,'(2A,I2,2(A,I4),A)') 'X ModifyGraph tick=2,mirror=0,'
     &   ,'font="Times",fSize=', FSIZE, ',width=', WIDTH,',height=', 
     &   HEIGHT,',standoff=0,btLen=3,stLen=2'
         WRITE(20,'(A)')  'X ModifyGraph width=0,height=0,margin(top)'//
     &   '=6,margin(right)=14,mirror=2,nticks(bottom)=15,minor=0,'//
     &   'manTick=0'
         WRITE(20,'(2A)') 'X ModifyGraph highTrip(left)=999999,',
     &   'notation(left)=1,nticks(left)=7,tickEnab(left)={-0.1,INF}'
* Rev 1.02 2001.06.27 Izumi
         CALL XATT
         IF (INDREF .EQ. 1) THEN
            WRITE(20,'(A)') 'X AppendToGraph yrefl vs xrefl'
            WRITE(20,'(A)') 'X ModifyGraph rgb(yrefl)=(64000,0,64000)'
         END IF
         WRITE(20,'(A,I2,A)') 'X Label left "\\F''Times''\\Z', LSIZE,
     &   'Intensity"'
         WRITE(20,'(A,I2,A,I2,A)') 
     &   'X Label bottom "\\F''Times''\\Z', LSIZE,
     &   '2\\F''Symbol''\\f02q\\f00\\F''Times''  / \\F''Symbol''°"'
      END IF
      CLOSE(UNIT=20)
      END

************************************************************************

* Rev 1.02 2001.06.27 Izumi {
      REAL FUNCTION SHPOS(J)
*     Peak positions (shifted)
      PARAMETER (NB=7000)
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /PROFI/ PCOT(NB),SHFT(NB),NTSIM(NB),FNORM(NB),TNTSIM(NB),
     &  FWHM(NB),ETA(NB),SQSG(NB),CTH(NB),SIGP(NB),DFWDG(NB),DEDFFG(NB),
     &  DFWDL(NB),DEDFFL(NB),SIGPAP(NB)
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      COMMON /ZSPFN/ PEAKSH(NB),SQH2(NB),SQH2L(NB),SQH2G(NB),ASYM(NB),
     &ETAL(NB),ETAH(NB),FNORM1(2,NB),COMLG(2,NB),COML(2,NB),COMG(2,NB),
     &COEFL(2,NB),COEFG(2,NB),DECAY(2,NB),COMPEAR(2,NB),CDX(2,NB),
     &CDA1(2,NB),CDA2(2,NB),CDW(2,NB),CDRL1(2,NB),CDRL2(2,NB),
     &CDRH1(2,NB),DMPART(2,NB)
      DATA RD2/57.29578/

      SELECT CASE (NPRFN)
         CASE (0)
            SHPOS = PEAK(J) - SHFT(J)
         CASE (1:3)
            SHPOS = PEAK(J) - PEAKSH(J)
      END SELECT
      SHPOS = SHPOS*RD2
      END
* }

************************************************************************

      SUBROUTINE COMPOS(A,SIGMAA,LSIG,ANAME,INDIV,AW)
*     PRINT OUT THE NUMBER OF ATOMS IN THE UNIT CELL AND DISPLACEMENT 
*     PARAMETERS OF ALL THE SITES
      PARAMETER (NB=7000,NP=80000,NT=999,NS=48,NAP=150,NPH=8)
*     2*PI**2 = 19.73921
      PARAMETER (C=19.73921)
      REAL A(*),SIGMAA(*),G1(3,3),GG1(3,3),BETA(3,3),B,AW(*),
     &  SIGBET(6),SZMVP(NPH)
      INTEGER INDIV(*),R,RX,RY,RZ,LCAL(NAP),IG(NAP)
*Rev 1.0w 2001.04.12 Izumi
      CHARACTER ANAME(*)*5,PARNAM*60,PHNAME*25,LINE*136,FORM*6
      COMMON /CC/ APR(NT)
      COMMON /FRDA/ NFR,NMEM,MEED(15),EPSD,SCIO,TSCAT,TSCAT1,TSCAT2,
     &  UCLAG,NDA,NFPX(3)
      COMMON /MULT/ SM(NAP)
      COMMON /H/ NQ(NPH),NCENTR(NPH),LAUEG(NPH)
      COMMON /HALF/ IFAC(NAP,NPH)
      COMMON /LTM/ LISO(NAP,NPH),LMAG(NPH)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /PAR/ PARNAM(NT),PHNAME(NPH)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /RT/ T(3,NS,NPH),HT(NS,NB),NSYM(NPH),NSITE(NPH),
     &  IDSYM(NS,NAP,NPH),R(3,3,NS,NPH),RX(NS,NB),RY(NS,NB),RZ(NS,NB)
      COMMON /SITE/ NOAT(NAP,NPH)
      COMMON /T/ FF(NB),COEF(NPH)
      COMMON /TENS/ G(6,NPH),GG(6,NPH)
      ITB(I,J)=NINT(1.0E6*BETA(I,J))
      ITU(I,J)=NINT(1.0E6*BETA(I,J)/(C*SQRT(GG(I,L)*GG(J,L))))
      ISU(J1,J2,J3)=NINT(SIGBET(J1)/(C*SQRT(GG(J2,L)*GG(J3,L))))

      DATA RD/0.01745329/

*     SZMV IS USED TO CALCULATE WEIGHT FRACTIONS
      SZMV = 0.0
      DO L = 1, NPHASE
         IF (NMODE .GE. 4 .AND. NSITE(L) .EQ. 0) CYCLE
*        WRITE LATTICE PARAMETERS
         WRITE(6,'(/1H0,10X,2A)') 'Lattice parameters (Angstrom or '//
     &   'degree) and unit-cell volume (Angstrom**3) in ',PHNAME(L)
         WRITE(6,'(/1H ,12X,A)')
     &   '    a         b         c       alpha     beta      gamma'//
     &   '        V'
         CALL INVERS(GG(1,L),G(1,L),DET)
*        Unit cell volume
*        Refer to Eq. (3.22) in Y. Saito, "Kagaku Kesshogaku Nyumon," 
*        Kyoritsu, Tokyo (1975), p. 53.
         VOL=SQRT(DET)
         LA=KPHB(L)+NAX
         WRITE(6,'(1H ,10X,3F10.5,3F10.4,F12.4)') (APR(J),J=LA,LA+5),VOL
         IF (LSIG .EQ. 1) THEN
            LINE=' '
            DO J = 21, 71, 10
               LINE(J:J)='-'
            END DO
            DO J = 0, 5
               IF (J .LE. 2) THEN
*                 A, B AND C
                  FORM='(F8.5)'
               ELSE
*                 ALPHA, BETA AND GAMMA
                  FORM='(F8.4)'
               END IF
               IF (ID(LA+J) .EQ. 1) THEN
                  WRITE(LINE(14+10*J:21+10*J),FORM) SIGMAA(LA+J)
               END IF
            END DO

*           STANDARD DEVIATION OF THE UNIT-CELL VOLUME
*           C. Giacovazzo, "Fundamentals of Crystallography," Oxford Univ.
*           Press, Oxford (1992), p. 122.
            A1 = SIGMAA(LA)/APR(LA)
            A2 = SIGMAA(LA+1)/APR(LA+1)
            A3 = SIGMAA(LA+2)/APR(LA+2)
            B1 = SIN(APR(LA+3)*RD)*(COS(APR(LA+3)*RD) - 
     &           COS(APR(LA+4)*RD)*COS(APR(LA+5)*RD))*SIGMAA(LA+3)
            B2 = SIN(APR(LA+4)*RD)*(COS(APR(LA+4)*RD) - 
     &           COS(APR(LA+3)*RD)*COS(APR(LA+5)*RD))*SIGMAA(LA+4)
            B3 = SIN(APR(LA+5)*RD)*(COS(APR(LA+5)*RD) - 
     &           COS(APR(LA+3)*RD)*COS(APR(LA+4)*RD))*SIGMAA(LA+5)
            SIGV = VOL*VOL*(A1*A1 + A2*A2 + A3*A3) +
     &             APR(LA)*APR(LA+1)*APR(LA+2)/(VOL*VOL)*
     &             (B1*B1 + B2*B2 + B3*B3)
            WRITE(LINE(72:83),'(F12.4)') SQRT(SIGV)
            WRITE(6,'(A)') LINE
         END IF

         DO J = 1, NSITE(L)
            LCAL(J)=0
*           DETERMINE THE PARAMETER NUMBERS OF OCCUPATION FACTORS
            IF (J .EQ. 1) THEN
               IG(J)=KPHB(L)+NG1X
            ELSE
               IG(J)=IG(J-1)+NPSITE(J-1,L)
            END IF
         END DO

*        WRITE STRUCTURE PARAMETERS
         IF (INDIV(L) .EQ. 1) THEN
            WRITE(6,'(//1H ,10X,2A/1H ,10X,A//1H ,23X,A)') 
*Rev 1.0w 2001.04.12 Izumi {
     &      'Structure parameters, g, x, y, z, B/Angstrom**2, and '//
     &      'U/Angstrom**2, in ',PHNAME(L),
     &      '                                  100*B/nm**2     '//
     &      '   100*U/nm**2',
     &      'neq  *   g    =   n        x         y         z        B'
     &      //'        U'
* }
         ELSE
            WRITE(6,'(//1H ,10X,2A//1H ,23X,A)') 
     &      'g, x, y, and z in ',PHNAME(L),
     &      'neq  *  g    =    n        x         y         z'
         END IF

         DO J = 1, NSITE(L)
            NPOS = 0
            DO JJ = 1, NSYM(L)
               IF (IDSYM(JJ,J,L) .EQ. 1) NPOS = NPOS + 1
            END DO
            IF (NCENTR(L) .EQ. 1) THEN
               NPOS = IFAC(J,L)*NINT(COEF(L))*NPOS/2
            ELSE
               NPOS = NINT(COEF(L))*NPOS
            END IF
            IF (INDIV(L) .EQ. 1 .AND. LISO(J,L) .NE. 6) THEN
*Rev 1.0w 2001.04.12 Izumi {
*              1/(8*PI**2) = 0.01266515
*              UU and SUU must be calculated before WRITE statements to
*              avoid a bug in Absoft Pro Fortran
               UU = A(IG(J)+4)*0.01266515
               SUU = SIGMAA(IG(J)+4)*0.01266515
               WRITE(6,'(13X,A,I5,2F9.4,3F10.5,F8.3,F10.5)') 
     &         PARNAM(IG(J))(1:9),NPOS,A(IG(J)),NPOS*A(IG(J)),
     &         (A(JJ),JJ=IG(J)+1,IG(J)+4),UU
               IF (LSIG .EQ. 1) 
     &         WRITE(LINE,'(27X,2F9.4,3F10.5,F8.3,F10.5)') 
     &         SIGMAA(IG(J)),SIGMAA(IG(J))*NPOS,
     &         (SIGMAA(JJ),JJ=IG(J)+1,IG(J)+4),SUU
            ELSE
               WRITE(6,'(13X,A,I5,2F9.4,3F10.5,A)') PARNAM(IG(J))(1:9),
     &         NPOS,A(IG(J)),NPOS*A(IG(J)),(A(JJ),JJ=IG(J)+1,IG(J)+3),
     &         '       -         -'
               IF (LSIG .EQ. 1) WRITE(LINE,'(27X,2F9.4,3F10.5,A)')
     &         SIGMAA(IG(J)),SIGMAA(IG(J))*NPOS,
     &         (SIGMAA(JJ),JJ=IG(J)+1,IG(J)+3),'       -         -'
* }
            END IF

            IF (LSIG .EQ. 0) CYCLE
*           PLACE '-' FOR neq, n, AND PARAMETERS WITH ID(I)<>1
*           PLACE '-' FOR neq
            LINE(27:27)='-'
            IF (LINE(29:34) .EQ. ' 0.0 ' .OR.
     &      LINE(30:36) .EQ. ' 0.0000') THEN
*              PLACE '-' IF ID(I)<>1 (g)
               LINE(31:36)='     -'
            END IF
            IF (LINE(39:43) .EQ. ' 0.0 ' .OR.
     &      LINE(39:45) .EQ. ' 0.0000') THEN
*              PLACE '-' IF ID(I)<>1 (n)
               LINE(40:45)='     -'
            END IF
            IF (LINE(48:52) .EQ. ' 0.0 ' .OR.
     &      LINE(48:55) .EQ. ' 0.00000') THEN
*              PLACE '-' IF ID(I)<>1 (x)
               LINE(49:55)='      -'
            END IF
            IF (LINE(58:62) .EQ. ' 0.0 ' .OR.
     &      LINE(58:65) .EQ. ' 0.00000') THEN
*              PLACE '-' IF ID(I)<>1 (y)
               LINE(59:65)='      -'
            END IF
            IF (LINE(68:72) .EQ. ' 0.0 ' .OR.
     &      LINE(68:75) .EQ. ' 0.00000') THEN
*              PLACE '-' IF ID(I)<>1 (z)
               LINE(69:75)='      -'
            END IF
            IF (LINE(78:82) .EQ. ' 0.0 ' .OR.
     &      LINE(78:83) .EQ. ' 0.000') THEN
*              PLACE '-' IF ID(I)<>1 (B)
               LINE(79:83)='    -'
            END IF
            WRITE(6,'(A)') LINE
         END DO
         WRITE(6,'(/13X,A)') 'neq: multiplicity of the Wyckoff position'
     &   //' (number of equivalent points per unit cell)'
         WRITE(6,'(13X,A)') '  n: number of equivalent atoms per '//
     &   'unit cell'

         IF (INDIV(L) .EQ. 1) THEN
*Rev 1.0w 2001.04.12 Izumi {
            WRITE(6,'(//11X,2A)') 'Atomic displacement parameters, '//
     &      '10**6*betaij, 10**6*Uij/Angstrom**2, Beq/Angstrom**2, and '
     &      //'Ueq/Angstrom**2, in ',PHNAME(L)
            WRITE(6,'(11X,A)')    '                                 '//
     &      '             10**8*Uij/nm**2        100*Beq/nm**2        '
     &      //'100*Ueq/nm**2'
            WRITE(6,'(/24X,A)') 'beta11  beta22  beta33  beta12'//
     &      '  beta13  beta23     U11     U22     U33     '//
     &      'U12     U13     U23     Beq       Ueq'
* }
            G1(1,1)=G(1,L)
            G1(2,2)=G(2,L)
            G1(3,3)=G(3,L)
            G1(2,3)=G(4,L)
            G1(1,3)=G(5,L)
            G1(1,2)=G(6,L)
            GG1(1,1)=GG(1,L)
            GG1(2,2)=GG(2,L)
            GG1(3,3)=GG(3,L)
            GG1(2,3)=GG(4,L)
            GG1(1,3)=GG(5,L)
            GG1(1,2)=GG(6,L)
            DO J = 1, NSITE(L)
               IF (LISO(J,L) .EQ. 6) THEN
                  BETA(1,1)=A(IG(J)+4)
                  BETA(2,2)=A(IG(J)+5)
                  BETA(3,3)=A(IG(J)+6)
                  BETA(1,2)=A(IG(J)+7)
                  BETA(1,3)=A(IG(J)+8)
                  BETA(2,3)=A(IG(J)+9)
                  CALL TOB(G1,BETA,BEQ)
                  IF (LSIG .EQ. 1) THEN
                     DO I = 1, 6
                        SIGBET(I)=SIGMAA(IG(J)+I+3)*1.0E6
                     END DO
                  END IF
*Rev 1.0w 2001.04.12 Izumi {
                  WRITE(6,'(13X,A,12I8,F8.3,F10.5)') PARNAM(IG(J))(1:9),
     &            ITB(1,1),ITB(2,2),ITB(3,3),ITB(1,2),ITB(1,3),ITB(2,3),
     &            ITU(1,1),ITU(2,2),ITU(3,3),ITU(1,2),ITU(1,3),ITU(2,3),
     &            BEQ,BEQ*0.01266515
               ELSE
                  B=A(IG(J)+4)
                  CALL TOBETA(GG1,B,BETA)
                  IF (LSIG .EQ. 1) THEN
                     DO I = 1, 6
                        SIGBET(I)=0.0
                     END DO
                  END IF
                  WRITE(6,'(13X,A,12I8,A)') PARNAM(IG(J))(1:9),
     &            ITB(1,1),ITB(2,2),ITB(3,3),ITB(1,2),ITB(1,3),ITB(2,3),
     &            ITU(1,1),ITU(2,2),ITU(3,3),ITU(1,2),ITU(1,3),ITU(2,3),
     &            '       -         -'
               END IF
               IF (LSIG .EQ. 1) THEN
                  WRITE(LINE,'(22X,12I8,A)') (NINT(SIGBET(I)),
     &            I=1,6),(ISU(I,I,I),I=1,3),ISU(4,1,2),ISU(5,1,3),
     &            ISU(6,2,3),'       -         -'
* }
                  DO I = 23, 111, 8
                     IF (LINE(I:I+7) .EQ. '       0') THEN
                        LINE(I+7:I+7)='-'
                     END IF
                  END DO
                  WRITE(6,'(A)') LINE
               END IF
            END DO
      
            WRITE(6,'(6(/13X,A))') 
*Rev 1.0w 2001.04.12 Izumi {
*           Refer to C. Giacovazzo, "Fundamentals of Crystallography," 
*           ed. by C. Giacovazzo, Oxford Univ. Press, Oxford (1992), p. 149.
     &      '  Isotropic Debye-Waller factor = '//
     &      'exp(-B*(sin(theta)/lambda)**2) = '//
     &      'exp(-8*pi**2*U*(sin(theta)/lambda)**2)',
     &      'Anisotropic Debye-Waller factor = exp(-(h**2*beta11 '//
     &      '+ k**2*beta22 + l**2*beta33 + ',
     &     '                                  2*h*k*beta12 + 2*h*l*beta'
     &     //'13 + 2*k*l*beta23))', '                                = '
     &      //'exp(-2*pi**2*(h**2*(a*)**2*U11 + k**2*(b*)**2*U22 + '//
     &      'l**2*(c*)**2*U33 +',
     &      '                                  2*h*k*(a*)*(b*)*U12 + '//
     &      '2*h*l*(a*)*(c*)*U13 + 2*k*l*(b*)*(c*)*U23))',
     &      'Beq, Ueq: equivalent isotropic atomic displacement '//
     &      'parameters calculated from anisotropic atomic displacement'
     &      //' parameters'
* }
         END IF

         WRITE(6,'(//11X,2A)') 'Number and weight of each species in '
     &   //'the unit cell, and density for ',PHNAME(L)
         WRITE(6,'(/13X,A)') 'Atom       N     *    At.wt.  /  '//
     &   '6.02214E23  =     Wt.'
         ZM=0.0
         DO IPOINT = 1, NSITE(L)
            IF (LCAL(IPOINT) .EQ. 1) CYCLE
            X=0.0
            DO J = IPOINT, NSITE(L)
               IF (LCAL(J) .EQ. 0 .AND.
     &         ABS(NOAT(J,L)) .EQ. ABS(NOAT(IPOINT,L))) THEN
                  NPOS=0
                  DO JJ = 1, NSYM(L)
                     IF (IDSYM(JJ,J,L) .EQ. 1) NPOS=NPOS+1
                  END DO
                  IF (NCENTR(L) .EQ. 1) THEN
                     X = X + A(IG(J))*FLOAT(IFAC(J,L)*NPOS)*COEF(L)/2.0
                     IF (L .EQ. NDA)
     &               SM(J) = 0.5*A(IG(J))*FLOAT(IFAC(J,L)*NPOS)/
     &                       FLOAT(NSYM(L))
                  ELSE
                     X = X + A(IG(J))*COEF(L)*FLOAT(NPOS)
                     IF (L .EQ. NDA) SM(J) = A(IG(J))*FLOAT(NPOS)/
     &                                       FLOAT(NSYM(L))
                  END IF
                  LCAL(J)=1
               END IF
            END DO
            WRITE(6,'(13X,A,F11.5,F13.5,15X,1P,E15.6,A)')
     &      ANAME(ABS(NOAT(IPOINT,L))),X,AW(ABS(NOAT(IPOINT,L))),
     &      X*AW(ABS(NOAT(IPOINT,L)))/6.022137E23,' g'
            ZM=ZM+AW(ABS(NOAT(IPOINT,L)))*X
         END DO
         WRITE(6,'(13X,61A)') ('-',JJ=1,61)
         SUMAW=ZM/6.022137E23
         WRITE(6,'(52X,A,1P,E13.6,A)') 'Total =',SUMAW,' g'
*        VOL: UNIT-CELL VOLUMN/CM**3; D: DENSITY/G CM-3
         VOL=VOL*1.0E-24
         WRITE(6,'(13X,A,1P,E13.6,A,1P,E13.6,A,0P,F9.6,A)')
     &   'd = Total/V =',SUMAW,' /',Vol,' = ',SUMAW/VOL,' g/cm**3'
         SZMVP(L) = A(KPHB(L)+NSFX)*ZM*VOL
         SZMV = SZMV + SZMVP(L)  
      END DO
      IF (NPHASE .EQ. 1 .OR. NMODE .GE. 4) RETURN

      WRITE(6,'(//11X,A)') 'Mass fractions of compounds contained'//
     &  ' in the sample'
      DO L = 1, NPHASE
         WRITE(6,'(13X,A,F7.4)') PHNAME(L),SZMVP(L)/SZMV
      END DO
      END

************************************************************************

      SUBROUTINE TOB(G,BETA,B)
*     G: METRIC TENSOR(INPUT)
*     BETA: ANISOTROPIC DISPLACEMENT PARAMETERS (INPUT)
*     B: ISOTROPIC DISPLACEMENT PARAMETER (OUTPUT)
      REAL G(3,3),BETA(3,3),B

      B = 1.333333*(BETA(1,1)*G(1,1) + BETA(2,2)*G(2,2) + 
     &BETA(3,3)*G(3,3) + 2.0*(BETA(1,2)*G(1,2) + BETA(1,3)*G(1,3) +
     &BETA(2,3)*G(2,3)))
      END

************************************************************************

      SUBROUTINE TOBETA(GG,B,BETA)
*     GG: METRIC TENSORS OF THE RECIPROCAL CELL (INPUT)
*     B: ISOTROPIC DISPLACEMENT PARAMETER (INPUT)
*     BETA: ANISOTROPIC DISPLACEMENT PARAMETERS (OUTPUT)
      REAL GG(3,3),B,BETA(3,3)

      C=0.25*B
      BETA(1,1)=C*GG(1,1)
      BETA(2,2)=C*GG(2,2)
      BETA(3,3)=C*GG(3,3)
      BETA(2,3)=C*GG(2,3)
      BETA(1,3)=C*GG(1,3)
      BETA(1,2)=C*GG(1,2)
      END

************************************************************************

      SUBROUTINE RFLCTN(LLIST)
C     PRINT OUT PHASE #, INDICES, INTENSITIES, BRAGG ANGLES, FWHM'S,
C     MULTIPLICITIES, AND STRUCTURE FACTORS. 
      INTEGER H
      PARAMETER (NB=7000,NPH=8,NA=15,NAP=150)
      CHARACTER*9 FWHMCH,OBSINT
      CHARACTER*10 FNUCL,FMAGN
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /B/ U(NB),NREFL
      COMMON /OBSCAL/ OBSI(NB),CALCI(NB),LREFINT(NB)
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /LTM/ LISO(NAP,NPH),LMAG(NPH)
      COMMON /MAGNE/ FMAG(NA,NB),SFRMAG(NB),SFIMAG(NB),Q2(NB)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /ORDER/ IPOINT(NB)
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /PROFI/ PCOT(NB),SHFT(NB),NTSIM(NB),FNORM(NB),TNTSIM(NB),
     &  FWHM(NB),ETA(NB),SQSG(NB),CTH(NB),SIGP(NB),DFWDG(NB),DEDFFG(NB),
     &  DFWDL(NB),DEDFFL(NB),SIGPAP(NB)
*Rev 1.0d 2000.11.09 Kajimoto
      COMMON /SFRI/ SFR(NB),SFI(NB),OTF(NB)
      COMMON /T/ FF(NB),COEF(NPH)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      COMMON /ZSPFN/ PEAKSH(NB),SQH2(NB),SQH2L(NB),SQH2G(NB),ASYM(NB),
     &ETAL(NB),ETAH(NB),FNORM1(2,NB),COMLG(2,NB),COML(2,NB),COMG(2,NB),
     &COEFL(2,NB),COEFG(2,NB),DECAY(2,NB),COMPEAR(2,NB),CDX(2,NB),
     &CDA1(2,NB),CDA2(2,NB),CDW(2,NB),CDRL1(2,NB),CDRL2(2,NB),
     &CDRH1(2,NB),DMPART(2,NB)
      DATA RD2/57.29578/

      IF (LLIST .EQ. 0) THEN
         WRITE(6,100)
      ELSE
         WRITE(6,200)
      ENDIF
  100 FORMAT(/12X,'No.',3X,'Phase',3X,'h',4X,'k',4X,'l',4X,'Code',
     &3X,'2-theta',8X,'d',7X,'Ical',2X,
     &'|F(nucl)| |F(magn)|',4X,'POF',5X,'FWHM',4X,'m',5X,'Dd/d')
  200 FORMAT(/12X,'No.',3X,'Phase',3X,'h',4X,'k',4X,'l',4X,'Code',
     &3X,'2-theta',8X,'d',7X,'Iobs',5X,'Ical',2X,
     &'|F(nucl)| |F(magn)|',4X,'POF',5X,'FWHM',4X,'m',5X,'Dd/d')
     
      YMAX = 0.0
      DO I = 1, NREFL
         IF (LLIST .EQ. 0 .AND. YPEAK(I) .GT. YMAX) YMAX = YPEAK(I)
         IF (LLIST .EQ. 1 .AND. OBSI(I) .GT. YMAX) YMAX = OBSI(I)
      END DO
      YADJ = 100000.0/YMAX
      DO J = 1, NREFL
         I = IPOINT(J)
         IF (L12(I) .EQ. 0) THEN
            LLL = 1
            M = NINT(U(I))
         ELSE
            LLL = 2
            M = NINT(U(I)/R12)
         ENDIF

         FNUCL = '      -'
         FMAGN = '      -'
         IF (NMODE .LE. 3 .OR. NMODE .EQ. 5 .OR. (NMODE .EQ. 4 .AND.
     &   NOPH(I) .GE. 2)) THEN
*           |Fc(nucl)|
*Rev 1.0d 2000.11.09 Kajimoto {
*           WRITE(FNUCL,'(F10.4)') SQRT(FF(I))
            WRITE(FNUCL,'(F10.4)') SQRT(SFR(I)**2 + SFI(I)**2)
* }
*           |Fc(magn)|
            IF (NMODE .LE. 1 .AND. LMAG(NOPH(I)) .EQ. 1) 
     &      WRITE(FMAGN,'(F10.4)') 
     &      SQRT(7.266087*Q2(I)*(SFRMAG(I)**2 + SFIMAG(I)**2))
         END IF
	 
         IF (NMODE .EQ. 1 .AND. NPAT .EQ. 0) THEN
            FWHMCH = '     -'
         ELSE IF (NPRFN .EQ. 0) THEN
            WRITE(FWHMCH,'(F9.4)') FWHM(I)*RD2
            RES = FWHM(I)/TAN(PEAK(I)*0.5)
         ELSE
            WRITE(FWHMCH,'(F9.4)') SQH2(I)*RD2
            RES = SQH2(I)/TAN(PEAK(I)*0.5)
         END IF
      
         IF (LLIST .EQ. 0) THEN
            WRITE(6,300) J,NOPH(I),H(I),K(I),L(I),LLL,PEAK(I)*RD2,D(I),
     &      NINT(YPEAK(I)*YADJ),FNUCL,FMAGN,PROR1(I),FWHMCH,M,RES
  300       FORMAT(10X,I4,I7,I6,2I5,I6,F12.3,F11.5,I9,2A,F8.3,A,I5,
     &      F10.5)
         ELSE
            SELECT CASE (LREFINT(I))
               CASE (0)
                  WRITE(OBSINT,'(I9)') NINT(OBSI(I)*YADJ)
               CASE (1)
                  WRITE(OBSINT,'(A)') '     -'
               CASE (2)
                  WRITE(OBSINT,'(A)') '     W'
               CASE (3)
                  WRITE(OBSINT,'(A)') '     H'
               CASE (4)
                  WRITE(OBSINT,'(A3,I6)') '  G',NINT(OBSI(I)*YADJ)
            END SELECT
	    
            WRITE(6,400) J,NOPH(I),H(I),K(I),L(I),LLL,PEAK(I)*RD2,D(I),
     &      OBSINT,NINT(CALCI(I)*YADJ),FNUCL,FMAGN,PROR1(I),FWHMCH,M,RES
  400       FORMAT(10X,I4,I7,I6,2I5,I6,F12.3,F11.5,A9,I9,2A,F8.3,A,
     &      I5,F10.5)
         END IF
      END DO
      END

************************************************************************

      SUBROUTINE FOURMEM(NSPGR,NSET,VNS,NEXCLD,NWEAK,NHIEND,NGRP)
*     OUTPUT *.hkl AND/OR *.mem FILES

      PARAMETER (NB=7000,NPH=8,NA=15,NAP=150,NT=999,NS=48)
      INTEGER H,NSPGR(*),NSET(*),R,RX,RY,RZ
      DIMENSION IOVLS(NB),IOVLE(NB),IOVLN(NB),JPOINT(NB),DSTORE(NB)
      CHARACTER VNS(NPH)*7
      CHARACTER TITMEM*70,FILE21*50,FILE30*50
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /B/ U(NB),NREFL
*Rev 1.0h 2000.12.13 Izumi
      COMMON /BIJVOET/ IPAIR(NB),LPAIR(NPH)
      COMMON /CC/ APR(NT)
      COMMON /OBSCAL/ OBSI(NB),CALCI(NB),LREFINT(NB)
      COMMON /FMEM/ MH(NB),MK(NB),ML(NB),FMEMR(NB),FMEMI(NB),NHKL
      COMMON /FRDA/ NFR,NMEM,MEED(15),EPSD,SCIO,TSCAT,TSCAT1,TSCAT2,
     &  UCLAG,NDA,NFPX(3)
      COMMON /FRDA2/ TITMEM
      COMMON /H/ NQ(NPH),NCENTR(NPH),LAUEG(NPH)
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /LTM/ LISO(NAP,NPH),LMAG(NPH)
      COMMON /MAGNE/ FMAG(NA,NB),SFRMAG(NB),SFIMAG(NB),Q2(NB)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /ORDER/ IPOINT(NB)
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /PROFI/ PCOT(NB),SHFT(NB),NTSIM(NB),FNORM(NB),TNTSIM(NB),
     &  FWHM(NB),ETA(NB),SQSG(NB),CTH(NB),SIGP(NB),DFWDG(NB),DEDFFG(NB),
     &  DFWDL(NB),DEDFFL(NB),SIGPAP(NB)
      COMMON /RT/ T(3,NS,NPH),HT(NS,NB),NSYM(NPH),NSITE(NPH),
     &  IDSYM(NS,NAP,NPH),R(3,3,NS,NPH),RX(NS,NB),RY(NS,NB),RZ(NS,NB)
      COMMON /SFRI/ SFR(NB),SFI(NB),OTF(NB)
      COMMON /ABPI/ AP(NB),BP(NB)
      COMMON /SIGMAS2/ SIGSCALS(NPH),SCALEAS(NPH)
      COMMON /T/ FF(NB),COEF(NPH)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      COMMON /ZSPFN/ PEAKSH(NB),SQH2(NB),SQH2L(NB),SQH2G(NB),ASYM(NB),
     &ETAL(NB),ETAH(NB),FNORM1(2,NB),COMLG(2,NB),COML(2,NB),COMG(2,NB),
     &COEFL(2,NB),COEFG(2,NB),DECAY(2,NB),COMPEAR(2,NB),CDX(2,NB),
     &CDA1(2,NB),CDA2(2,NB),CDW(2,NB),CDRL1(2,NB),CDRL2(2,NB),
     &CDRH1(2,NB),DMPART(2,NB)
      COMMON /UNCONS/ NUCREF,IPHASE(90),IHKL(3,90),NOPPP(90),LPPP(NB),
     &  RCUT(2,90),NCUT
      P(J) = APR(KPHB(NFR) + J)
      Q(J) = APR(KPHB(NMEM) + J)

      NEXCLD = 0
      NWEAK = 0
      NHIEND = 0
      NGRP = 0

      IF (NFR .GT. 0) THEN
         CALL OPENFILE(5,FILE21)
         REWIND 21
         WRITE(21,'(A/A/A/6F10.5/A/5I5,2A/A/F12.3)')
     &   '# The data were recorded by RIETAN-98 v.1.05a or later.', '#',
     &   '#        a         b         c     alpha      beta     gamma',
     &   (P(J), J = NAX, NGAMX),
     &   '#  NX   NY   NZ  SPG  ORG  VOL',(NFPX(J),J=1,3),
     &   NSPGR(NFR),NSET(NFR),'    ',VNS(NFR)(1:1),
     &   '# Total electrons or coherent scattering lengths',TSCAT
         WRITE(21,250)
  250    FORMAT('#',3X,'h',4X,'k',4X,'l',10X,'FO(R)',10X,'FO(I)',
     &   10X,'FC(R)',10X,'FC(I)')
      END IF
      
      IF (NMEM .GT. 0) THEN
         CALL OPENFILE(7,FILE30)
         REWIND 30
	 
         IF (MEED(7) .NE. 0) THEN
            IBFO = 0
            NN = 0
            DO J = 1, NREFL
               I = IPOINT(J)
               IF (NMEM .NE. NOPH(I) .OR. L12(I) .NE. 0 .OR.
     &         LREFINT(I) .GE. 2) CYCLE
               NN = NN + 1
               DSTORE(NN) = D(I)
               JPOINT(NN) = J
            END DO
            DO I = 1, NN
               J = JPOINT(I)
               II = IPOINT(J)
               IF (I .NE. NN) THEN
                  DLTD = DSTORE(I) - DSTORE(I+1)
                  ITMP = 0
                  IF (DLTD .LE. EPSD) ITMP = 1
                  IF (IBFO .EQ. 0 .AND. ITMP .EQ. 1) THEN
                     NGRP = NGRP + 1
                     IOVLS(NGRP) = J
                     IOVLN(NGRP) = 1
                     LREFINT(II) = 4
                  ELSE IF (IBFO .EQ. 1 .AND. ITMP .EQ. 1) THEN
                     IOVLN(NGRP) = IOVLN(NGRP) + 1
                     LREFINT(II) = 4
                  ELSE IF (IBFO .EQ. 1 .AND. ITMP. EQ. 0) THEN
                     IOVLE(NGRP) =J 
                     IOVLN(NGRP) = IOVLN(NGRP) + 1
                     LREFINT(II) = 4
                  END IF
               ELSE IF (IBFO .EQ. 1) THEN
                  IOVLE(NGRP) = J
                  IOVLN(NGRP) = IOVLN(NGRP) + 1
                  LREFINT(II) = 4
               END IF
               IBFO = ITMP
            END DO
         END IF

*        LREFINT(I) = 0: FULL OBSERVED (POSITIVE) AND CALCULATED 
*                        PROFILES ARE OBTAINED.
*        LREFINT(I) = 1: PART OF THE OBSERVED PROFILE IS LACKING AS A 
*                        RESULT OF EXCLUDING A 2-THETA REGION.
*        LREFINT(I) = 2: FULL OBSERVED (NEGATIVE) AND CALCULATED
*                        PROFILES ARE OBTAINED.
*        LREFINT(I) = 3: PART OF THE PROFILE EXCEEDS 2-THETA(MAX).
*        LREFINT(I) = 4: GROUP OF REFLECTIONS WITH NEARLY EQUAL D-SPACINGS.

         NREFMEED = 0
         LACK = 0
         DO I = 1, NREFL
            IF (NMEM .NE. NOPH(I) .OR. L12(I) .NE. 0) CYCLE
            IF (MEED(8).EQ.0) THEN
               SELECT CASE (LREFINT(I))
                  CASE (0, 1)
                     NREFMEED = NREFMEED + 1
                  CASE (2, 3)
                     LACK = LACK + 1
               END SELECT
            ELSE
               SELECT CASE (LREFINT(I))
                  CASE (0, 1, 2)
                     NREFMEED = NREFMEED + 1
                  CASE (3)
                     LACK = LACK + 1
               END SELECT
            END IF
         END DO

         WRITE(30,'(A/A)') '# These data were recorded by RIETAN-98','#'
         WRITE(30,'(3(A,I1))') '# MEED(1) = ',MEED(1),',  MEED(7) = ',
     &                         MEED(7),',  MEED(8) = ',MEED(8)
         WRITE(30,'(2(A,G13.6))') '# SCIO = ',SCIO,',  EPSD = ',EPSD
         WRITE(30,'(A/A/A/A/6F10.5/A/A/A/5I5,2A)')
     &   '#','# Title (up to 70 characters)', TITMEM,
     &   '#        a         b         c     alpha      beta     gamma',
     &   (Q(J), J = NAX, NGAMX),
     &   '#          POP1           POP2  IPTYP  JPH  JPK  JPL',
     &   '      0.0000000      0.0000000     -1    0    0    0     P',
     &   '#  NX   NY   NZ  SPG  ORG  VOL',(MEED(J), J = 4, 6),
     &   NSPGR(NMEM),NSET(NMEM),'    ',VNS(NMEM)(1:1)
         IF (TSCAT2 .EQ. 0.0) THEN
            NN = NREFMEED
*Rev 1.0r 2001.02.15 Izumi {
C           IF (NCENTR(NMEM) .EQ. 0 .AND. NBEAM .GE. 1) NN = NN/2
            IF (LPAIR(NMEM) .EQ. 1) NN = NN/2
* }
            WRITE(30,'(A/F17.3,F13.7/A/I5)') 
     &      '# Total electrons       lambda',TSCAT1, UCLAG, '# NREF1',NN
         ELSE 
            WRITE(30,'(A,15X,A/2F15.3,F15.7/A/I5)')
     &      '# Total SA b+ & b-', '      lambda', TSCAT1, TSCAT2, UCLAG,
     &      '# NREF1', NREFMEED
         END IF
         WRITE(30,'(2A)') '#   h    k    l     F''o''(R)        ',
     &   ' F''o''(I)        sigma(F)'
      END IF

      DO J = 1, NREFL
         I = IPOINT(J)

         SELECT CASE (LREFINT(I))
            CASE (1)
               NEXCLD = NEXCLD + 1
            CASE (2)
               NWEAK = NWEAK + 1
            CASE (3)
               NHIEND = NHIEND + 1
            CASE DEFAULT
*              DO NOTHING
         END SELECT

*        SKIP -h -k -l REFLECTIONS IN A NONCENTROSYMMETRIC SPACE GROUP (X-RAY)
         IF (IPAIR(I) .NE. 0) CYCLE

         IF (.NOT. ((NFR .EQ. NOPH(I) .OR. NMEM .EQ. NOPH(I)) .AND. 
     &       L12(I) .EQ. 0)) CYCLE
         IF (LREFINT(I) .EQ. 2 .OR. LREFINT(I) .EQ. 3) CYCLE

*        When refining |Fc| for a relaxed reflection, its observed integrated
*        intensity is regarded as better evaluated by numerical integration.
         
         IF ((NMODE .EQ. 3 .OR. NMODE .EQ. 6) .AND. LPPP(I) .GT. 0) 
     &   OBSI(I) = CALCI(I)

         SELECT CASE (NMODE)
            CASE (0,1,5)
               IF (NBEAM .EQ. 0) THEN
                  ACP = SFR(I)
                  BCP = SFI(I)
               ELSE
                  ACP = AP(I)
                  BCP = BP(I)
               END IF
            CASE (2,3)
               ACP = FMEMR(I)
               BCP = FMEMI(I)
         END SELECT
*        |Fo|/|Fc| = SQRT(OBSI(I)/CALCI(I))
         AFO = SQRT(OBSI(I)/CALCI(I))
*        FOS(R,I) = Fc(R,I)*|Fo|/|Fc| - delta(Fc(R,I))
*        delta(Fc(R,I)) gives correction term of anomoulous dispersion.
         FOSR = SFR(I)*(AFO - 1.0) + ACP
         FOSI = SFI(I)*(AFO - 1.0) + BCP
         FOSTAR = SQRT(FOSR*FOSR + FOSI*FOSI)
         FCSTAR = SQRT(ACP*ACP + BCP*BCP)
*        FROSTAR: Real part of Fo* = |Fo*|cos(psi)
         FROSTAR = FOSTAR*ACP/FCSTAR
*        FIOSTAR: Imaginary part of Fo* = |Fo*|sin(psi)
         FIOSTAR = FOSTAR*BCP/FCSTAR

*        IF (NBEAM .EQ. 0) THEN
*           FMOBS = FMAGN*SQRT(OBSI(I)/CALCI(I))
*           FMROBS = FMOBS*SFRMAG(I)/SQRT(SFRMAG(I)**2 + SFIMAG(I)**2)
*           FMIOBS = FMOBS*SFIMAG(I)/
*           SQRT(7.266087*Q2(I)*(SFRMAG(I)**2 + SFIMAG(I)**2))
*           SFMOBS = 0.5*FMOBS*SQRT(1.0/(SCIO*OBSI(I)) +
*    &               (SIGSCALS(NMEM)/SCALEAS(NMEM))**2)
*        END IF
	    
*        Output data for Fourier/D synthesis (*.hkl)

         IF (NFR .EQ. NOPH(I)) WRITE(21,'(3I5,4F15.7)')
     &   H(I),K(I),L(I),FROSTAR,FIOSTAR,ACP,BCP

*        Output data for MEM analysis with MEED (*.mem)

         IF (NMEM .NE. NOPH(I)) CYCLE
         IF (LREFINT(I) .LE. 1) THEN
            SELECT CASE (MEED(8))
               CASE (0)
                  IF (NBEAM.GE.1.AND.MEED(1).EQ.0) THEN
                     TOBSI = CALCI(I)*FOSTAR*FOSTAR/FF(I)
                  ELSE
                     TOBSI = OBSI(I)
                  END IF
                  SFOSTAR = 0.5*FOSTAR*SQRT(1.0/(SCIO*TOBSI) +
     &                      (SIGSCALS(NMEM)/SCALEAS(NMEM))**2)
               CASE (1)
                  IF (NBEAM.GE.1.AND.MEED(1).EQ.0) THEN
                     TCALCI = CALCI(I)*FCSTAR*FCSTAR/FF(I)
                  ELSE
                     TCALCI = CALCI(I)
                  END IF
                  SFOSTAR = 0.5*FCSTAR*SQRT(1.0/(SCIO*TCALCI) +
     &                      (SIGSCALS(NMEM)/SCALEAS(NMEM))**2)
                  FROSTAR = ACP
                  FIOSTAR = BCP
            END SELECT
                  
            IF (LREFINT(I) .EQ. 0) THEN
               WRITE(30,'(3I5,3F15.7)') H(I), K(I), L(I), 
     &         FROSTAR, FIOSTAR, SFOSTAR
            ELSE IF (LREFINT(I) .EQ. 1) THEN
               WRITE(30,'(3I5,3F15.7,A)') H(I), K(I), L(I),
     &         FROSTAR, FIOSTAR, SFOSTAR, '  # = Fc'
            END IF
            DSTORE(I) = 0.0
         ELSE IF (MEED(7) .NE. 0 .AND. LREFINT(I) .EQ. 4) THEN
            SELECT CASE (MEED(8))
               CASE (0)
                  DSTORE(I) = FOSTAR*FOSTAR
               CASE DEFAULT
                  DSTORE(I) = FCSTAR*FCSTAR
            END SELECT
         END IF
      END DO

      IF (NFR .GT. 0) CLOSE(UNIT=21)
      IF (NMEM .EQ. 0) RETURN

      IF (MEED(7) .NE. 0 .AND. NGRP .GT. 0) THEN
         WRITE(30,'(A/I5)') '# NREF2', NGRP
         WRITE(30,'(2A)') '#   h    k    l     G''o''           ',
     &   ' sigma(G)     M'
         DO IGRP = 1, NGRP
            WRITE(30,'(I5)') IOVLN(IGRP)
            DO J = IOVLS(IGRP), IOVLE(IGRP)
               N = IPOINT(J)
*              SKIP -h -k -l REFLECTIONS IN A NONCENTROSYMMETRIC SPACE 
*              GROUP (X-RAY)
               IF (IPAIR(N) .NE. 0) CYCLE
               IF (LREFINT(N) .NE. 4) CYCLE
               M = NINT(U(N))
               IF (NBEAM .GE. 1 .AND. NCENTR(NMEM)  .EQ. 0) M = 2*M
               IF (J .EQ. IOVLS(IGRP)) THEN
                  SUMOBS = 0.0
                  SUMA = 0.0
                  SUMB = 0.0
                  NN = 0
                  DO J2 = IOVLS(IGRP), IOVLE(IGRP)
                     I = IPOINT(J2)
                     IF (LREFINT(I) .NE. 4) CYCLE
                     IF (DSTORE(I) .EQ. 0.0) CALL JOBEND('OOPS !!!')
                     UM = U(I)
                     IF (NBEAM .GE. 1 .AND. NCENTR(NMEM) .EQ. 0) 
     &               UM = 2.0*UM
                     SELECT CASE (MEED(8))
                        CASE (0)
                           IF (NBEAM .GE. 1 .AND. MEED(1) .EQ. 0) THEN
                              SUMOBS = SUMOBS + CALCI(I)*DSTORE(I)/FF(I)
                           ELSE
                              SUMOBS = SUMOBS + OBSI(I)
                           END IF
                        CASE DEFAULT
                           IF (NBEAM.GE.1.AND.MEED(1).EQ.0) THEN
                              SUMOBS = SUMOBS + CALCI(I)*DSTORE(I)/FF(I)
                           ELSE
                              SUMOBS = SUMOBS + CALCI(I)
                           END IF
                     END SELECT
                     SUMA = SUMA + UM*FLP(I)*DSTORE(I)
                     SUMB = SUMB + UM*FLP(I)
                     NN = NN + 1
                  END DO
                  IF (NN .NE. IOVLN(IGRP)) CALL JOBEND('OOPS !!!')
                  GOBS = SQRT(SUMA/SUMB)
                  SIGG = 0.5*GOBS*SQRT(1.0/(SCIO*SUMOBS) +
     &                   (SIGSCALS(NMEM)/SCALEAS(NMEM))**2)
                  WRITE(30,'(3I5,2F15.7,I5)') H(N), K(N), L(N),
     &                                        GOBS, SIGG, M
               ELSE
                  WRITE(30,'(3I5,30X,I5)') H(N), K(N), L(N), M
               END IF
            END DO
         END DO
      ELSE
         WRITE(30,'(A/I5)') '# NREF2', 0
      END IF

      IF (LACK .GT. 0) THEN
         LL = LACK
         IF (NCENTR(NMEM) .EQ. 0 .AND. NBEAM .GE. 1) LL = LACK/2
         WRITE(30,'(A/I5)') '# NREF3', LL
         WRITE(30,'(2A)') '#   h    k    l'
         DO J = 1, NREFL
            I = IPOINT(J)
*           SKIP -h -k -l REFLECTIONS IN A NONCENTROSYMMETRIC SPACE 
*           GROUP (X-RAY)
            IF (IPAIR(I) .NE. 0) CYCLE
            IF (NMEM .EQ. NOPH(I) .AND. L12(I) .EQ. 0 .AND. 
     &      (LREFINT(I) .GE. 2 .AND. LREFINT(I) .LE. 3))
     &      WRITE(30,'(3I5)') H(I),K(I),L(I)
         END DO
      ELSE
         WRITE(30,'(A/I5)') '# NREF3', 0
      END IF
      CLOSE(UNIT=30)
      END

************************************************************************

      SUBROUTINE MKSFF(A,DEG,NPTS)
*     CREATE A FILE FOR STORING INTEGRATED INTENSITIES RESULTING FROM 
*     LE BAIL REFINEMENT
      INTEGER H
      PARAMETER (NB=7000,NPH=8,NAP=150,NP=80000)
      REAL A(*),DEG(*),YMAX(NB)
      CHARACTER FILE23*50
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /B/ U(NB),NREFL
      COMMON /G/ I
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /ORDER/ IPOINT(NB)
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /PROFI/ PCOT(NB),SHFT(NB),NTSIM(NB),FNORM(NB),TNTSIM(NB),
     &  FWHM(NB),ETA(NB),SQSG(NB),CTH(NB),SIGP(NB),DFWDG(NB),DEDFFG(NB),
     &  DFWDL(NB),DEDFFL(NB),SIGPAP(NB)
      COMMON /YOKO/ ABSCISSA(NP),NPOINTS,AWRY
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      COMMON /ZSPFN/ PEAKSH(NB),SQH2(NB),SQH2L(NB),SQH2G(NB),ASYM(NB),
     &ETAL(NB),ETAH(NB),FNORM1(2,NB),COMLG(2,NB),COML(2,NB),COMG(2,NB),
     &COEFL(2,NB),COEFG(2,NB),DECAY(2,NB),COMPEAR(2,NB),CDX(2,NB),
     &CDA1(2,NB),CDA2(2,NB),CDW(2,NB),CDRL1(2,NB),CDRL2(2,NB),
     &CDRH1(2,NB),DMPART(2,NB)
      DATA RD,RD2/0.01745329,57.29578/
      DATA TOL/8.726646E-06/ ! 0.0005 degrees
      P(J)=A(KPHB(NOPH(I))+J)

      CALL OPENFILE(11,FILE23)
      REWIND 23 
*     COPY NPTS INTO NPOINTS AND DEG INTO VARINV, WHICH WILL BE USED 
*     IN SUBROUTINE PRFUNC
      NPOINTS = NPTS
      DO J = 1, NPTS
         ABSCISSA(J) = DEG(J)
      END DO

      Y0 = 0.0
      DO J = 1, NREFL
         I = IPOINT(J)
         IF (.NOT. (NOPH(I) .EQ. 1 .AND. L12(I) .EQ. 0)) CYCLE
*        TTHPEAK: PEAK POSITION FOR THE I'TH REFLECTION
         IF (NPRFN .EQ. 0) THEN
            MODE = 0
*           THE PEAK POSITION IS SANDWICHED BETWEEN X AND XORF
*           +/- 0.8 DEGREES
            X = PEAK(I) - SHFT(I) - 0.8*RD
            XORF = PEAK(I) - SHFT(I) + 0.8*RD
            DO WHILE (.TRUE.)
               CALL SFMIN(X, XORF, MODE, TOL)
               IF (MODE .NE. 1) GO TO 2
               XORF = -PRFUNC(X,A)
            END DO
    2       IF (MODE .EQ. 4) CALL JOBEND('Error termination in SFMIN')
            TTHPEAK = X
            YMAX(I) = -YPEAK(I)*XORF
         ELSE
            TTHPEAK = PEAK(I) - PEAKSH(I)
            YMAX(I) = YPEAK(I)*PRFUNC(TTHPEAK,A)
         END IF
         IF (YMAX(I) .GT. Y0) Y0 = YMAX(I)
      END DO
      YNORM = 100.0/Y0

      WRITE(23,'(A)') '>  h   k   l           FWHM          |F|^2'//
     &'              d          2-th           Area           I/I0'
      DO J = 1, NREFL
         I = IPOINT(J)
         IF (.NOT. (NOPH(I) .EQ. 1 .AND. L12(I) .EQ. 0)) CYCLE
         IF (NPRFN .EQ. 0) THEN
            FWHMOUT = FWHM(I)*RD2
         ELSE
            FWHMOUT = SQH2(I)*RD2
         END IF
C        WRITE(23,'(3I4,3F15.4,F14.3,5X,1P,G14.6,0P,F11.1)') H(I),K(I),
         WRITE(23,'(3I4,F15.4,1X,1P,G14.6,0P,F15.4,F14.3,5X,1P,G14.6,
     &   0P,F11.1)') H(I),K(I),L(I),FWHMOUT,YPEAK(I)/(P(0)*U(I)*FLP(I)*
     &   PROR1(I)),D(I),PEAK(I)*RD2,YPEAK(I),YMAX(I)*YNORM
      END DO
      CLOSE(UNIT=23)
      END
