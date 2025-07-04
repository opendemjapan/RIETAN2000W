      SUBROUTINE REFPAR(LL,NTERMS)
C     DETERMINE WHICH PARAMERTERS TO REFINE IN EACH CYCLE
      PARAMETER(NT=999,NR=400)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /NONLIN/ NLESQ,NAUTO,IDCYCL,NCHNG,NPAR(50),IPAR(NR,50)
      COMMON /RLS/ NRANGE
      
      NRFN=0
      IF(LL.EQ.NCHNG+1) THEN
         DO 10 J=1,NTERMS
            ID(J)=IDSAVE(J)
            IF(ID(J).EQ.1) THEN
               NRFN=NRFN+1
               IR(J)=NRFN
            ENDIF
   10    CONTINUE
      ELSE
         DO 20 J=1,NTERMS
            IF(IDSAVE(J).EQ.2) THEN
               ID(J)=2
            ELSEIF(IDSAVE(J).EQ.0) THEN
               ID(J)=0
            ELSE
               DO 15 JJ=1,NPAR(LL)
                  IF(J.EQ.IPAR(JJ,LL)) THEN
                     ID(J)=1
                     NRFN=NRFN+1
                     IR(J)=NRFN
                     GO TO 20
                  ENDIF
   15          CONTINUE
               ID(J)=0
            ENDIF
   20    CONTINUE
      ENDIF
      END

************************************************************************

      SUBROUTINE CHKCNV(IDCONV,A,ASAVE,CHISQ1,CHISQR,ICONV,NTERMS)
C     IDCONV=1: CONVERGED; IDCONV=0: NOT CONVERGED
      PARAMETER (NT=999,EPSR=1.0E-3)
      DIMENSION A(*),ASAVE(*)
      COMMON /FL0/ NC,TK,CONV,NCONV,FINC,SUMPEN,SUMP1,FACTOR
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN

      IDCONV=0
      RELDEC=(CHISQ1-CHISQR)/CHISQR
      IF(RELDEC.LT.0.0 .OR. RELDEC.GT.CONV) THEN
         ICONV=0
         RETURN
      ELSE
         ICONV=ICONV+1
      ENDIF
      IF(ICONV.LT.NCONV) RETURN
      DX=0.0
      X=0.0
      DO 10 J=1,NTERMS
         IF(ID(J).EQ.1) THEN
            IF(ABS(A(J)-ASAVE(J)).GT.DX) DX=ABS(A(J)-ASAVE(J))
            IF(ABS(ASAVE(J)).GT.X) X=ABS(ASAVE(J))
         ENDIF
   10 CONTINUE
      IF(DX.LE.MAX(1.0,X)*EPSR) IDCONV=1
      END

************************************************************************

      SUBROUTINE IGNORE(NSTEP,NPTS,DEG,XINT)
C     IGNORE ALL DATA POINTS FARTHER THAN PC*FWHM FROM THE NEAREST
C     BRAGG PEAK
      PARAMETER (NB=7000,NP=80000)
      DIMENSION DEG(*),XINT(*)
      COMMON /ABSVR/ ABSORP(NP),VARLEN(NP),ROUGH(NP),NSURFR
      COMMON /B/ U(NB),NREFL
      COMMON /C/ BG(NP),BGINC(NP),DEGNOR(NP)
      COMMON /I/ RDEG(2,NB)
      COMMON /ORDER/ IPOINT(NB)
      COMMON /RLS/ NRANGE

C     NPTS: NUMBER OF INTENSITY DATA USED TO REFINE PARAMETERS
      IF (NRANGE .EQ. 0 .OR. NRANGE .EQ. 3) THEN
         NPTS=NSTEP
         RETURN
      ENDIF
C
C     PACK ARRAYS DEG, XINT, etc.
      NPTS=1
      ISTART=1
      DO K = 1, NSTEP
         DO J=ISTART,NREFL
            I=IPOINT(J)
            IF(DEG(K).GT.RDEG(2,I)) THEN
               ISTART=I+1
            ELSE
               IF(DEG(K).GE.RDEG(1,I)) THEN
                  DEG(NPTS)=DEG(K)
                  XINT(NPTS)=XINT(K)
                  BG(NPTS)=BG(K)
                  DEGNOR(NPTS)=DEGNOR(K)
                  ABSORP(NPTS)=ABSORP(K)
                  VARLEN(NPTS)=VARLEN(K)
                  ROUGH(NPTS)=ROUGH(K)
                  NPTS=NPTS+1
               ENDIF
               EXIT
            ENDIF
         END DO
      END DO
      NPTS=NPTS-1
      END

************************************************************************

      SUBROUTINE INVERS(G1,G2,DET2)
C     G2 === INVERT ===> G1
      DIMENSION G1(6),G2(6)
      DOUBLE PRECISION A11,A22,A33,A23,A13,A12,DET
      A11=G2(1)
      A22=G2(2)
      A33=G2(3)
      A23=G2(4)
      A13=G2(5)
      A12=G2(6)
      DET=(A11*A33-A13*A13)*A22+(A12*A13-A11*A23)*A23+
     &    (A13*A23-A12*A33)*A12
      DET2=DET
      G1(1)=(A22*A33-A23*A23)/DET
      G1(2)=(A11*A33-A13*A13)/DET
      G1(3)=(A11*A22-A12*A12)/DET
      G1(4)=(A12*A13-A11*A23)/DET
      G1(5)=(A12*A23-A13*A22)/DET
      G1(6)=(A13*A23-A12*A33)/DET
      END

************************************************************************

      SUBROUTINE GENREF(AA,NSPGR)
C     GENERATE REFLECTIONS WHICH ARE LOCATED BETWEEN DEG1 AND DEG2
      PARAMETER (NS=48,NA=15,NB=7000,NPH=8,NAP=150)
      INTEGER SECSET,ISET(44),H,R,RX,RY,RZ
      REAL A,B,C,ALPHA,BETA,GAMMA,LAMBDA,ANGMAX,ANGMIN,X,Y,Z
*Rev 1.0i 2000.12.18 Izumi 
      CHARACTER SPGR*10,FORMAT_HKLM*25,PHASENO(NPH)*1
      DIMENSION NSPGR(NPH),AA(*)
*Rev 1.0i 2000.12.18 Izumi 
      LOGICAL LPHASE
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /B/ U(NB),NREFL
*Rev 1.0h 2000.12.13 Izumi
      COMMON /BIJVOET/ IPAIR(NB),LPAIR(NPH)
      COMMON /G/ I
      COMMON /H/ NQ(NPH),NCENTR(NPH),LAUEG(NPH)
      COMMON /INPT/ A,B,C,ALPHA,BETA,GAMMA,LAMBDA,ANGMAX,ANGMIN
      COMMON /IS/ ISET,SECSET
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /ORDER/ IPOINT(NB)
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /PHNO/ IP
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /RT/ T(3,NS,NPH),HT(NS,NB),NSYM(NPH),NSITE(NPH),
     &  IDSYM(NS,NAP,NPH),R(3,3,NS,NPH),RX(NS,NB),RY(NS,NB),RZ(NS,NB)
      COMMON /S/ F(NA,NB),NATOM
      COMMON /S2/ F0(NA,NB)
      COMMON /SG/ SPGR(NPH)
      COMMON /SGN/ NSG,NSC
      COMMON /T/ FF(NB),COEF(NPH)
      COMMON /TENS/ G(6,NPH),GG(6,NPH)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
      COMMON /UNCONS/ NUCREF,IPHASE(90),IHKL(3,90),NOPPP(90),LPPP(NB),
     &  RCUT(2,90),NCUT
      COMMON /FT/ FT1,FT2,PEAKNOR(NB),TANTHNOR(NB) 
*Rev 1.0h 2000.12.13 Izumi
      DATA PHASENO/'1','2','3','4','5','6','7','8'/
      Z(J,X,Y)=ACOS(G(J,IP)/(X*Y))*57.29578
      
      NREFL=0
      DO IP = 1, NPHASE
         A=SQRT(G(1,IP))
         B=SQRT(G(2,IP))
         C=SQRT(G(3,IP))
         ALPHA=Z(4,B,C)
         BETA=Z(5,C,A)
         GAMMA=Z(6,B,A)
         NSG=ABS(NSPGR(IP))
*        IF(NSPGR(IP).LT.0) THEN
*           SECSET=2
*        ELSE
*           SECSET=0
*        ENDIF
         SECSET=0
         IF (LAUEG(IP) .EQ. 3) THEN
*           MONOCLINIC, B-AXIS UNIQUE
            SECSET=2
         ELSE IF (LAUEG(IP) .EQ. 9 .OR. LAUEG(IP) .EQ. 11) THEN
*           TRIGONAL(R LATTICE), HEXAGONAL SETTING
            IF (SPGR(IP)(1:1) .EQ. 'R') SECSET=2
         END IF
         IS = NREFL + 1

*Rev 1.0k 2000.12.28 Izumi {
*        Read hkl and m from a file, hklm.phn (n: phase number) whose first 
*        line gives a format for reading them
         INQUIRE(FILE='hklm.ph'//PHASENO(IP),EXIST=LPHASE)
         IF (LPHASE) THEN
            WRITE(6,'(//11X,2A)') 'Indices and multiplicities have been'
     &      //' input from hklm.ph'//PHASENO(IP)
            OPEN(60,FILE='hklm.ph'//PHASENO(IP),STATUS='OLD')
            READ(60,'(A)') FORMAT_HKLM
            DO I = IS, NB
               READ(60,FORMAT_HKLM,END=3) H(I),K(I),L(I),MM
               IF (H(I) .EQ. 0 .AND. K(I) .EQ. 0 .AND. L(I) .EQ. 0)
     &         GO TO 3
               U(I) = FLOAT(MM)
               NOPH(I) = IP
            END DO
            CALL JOBEND('Too many reflections have been input '//
     &      'from hklm.ph'//CHAR(IP))
    3       NREFL = I - 1
         ELSE
            CALL KDRREF
         END IF
* }
         IE = NREFL
         DO I = IS, IE
            IPAIR(I) = 0
         END DO
         
*Rev 1.0h 2000.12.13 Izumi
*        In the case of X-ray powder diffraction data of a phase with no
*        inversion center, Bijvoet pairs, i.e., hkl and -h -k -l reflections, 
*        are generated because F(hkl) is more or less different from F(-h -k -l).
C        IF (NCENTR(IP) .EQ. 0 .AND. NBEAM .GE. 1 .AND. NMODE .LE. 3) 
C    &   THEN
         IF (LPAIR(IP) .EQ. 1 .AND. NMODE .LE. 3) THEN
* }
            NREFTP = IE - IS + 1
            IF (IE + NREFTP .GT. NB) CALL JOBEND('The total number of '
     &      //'reflections has reached the maximum number')  
            DO IREF = IE + 1, IE + NREFTP
               IPAIR(IREF) = 1
               H(IREF) = -H(IREF-NREFTP)
               K(IREF) = -K(IREF-NREFTP)
               L(IREF) = -L(IREF-NREFTP)
*              EACH MULTIPLICITY IS DIVIED BY 2
               U(IREF-NREFTP) = U(IREF-NREFTP)*0.5
               U(IREF) = U(IREF-NREFTP)
               NOPH(IREF) = IP
            END DO
            NREFL = IE + NREFTP
         END IF
      END DO
      IF (R12 .GT. 0.0 .AND. NREFL*2 .GT. NB) CALL JOBEND
     &('The number of reflections has reached the maximum value')

      DO I = 1, NREFL
         L12(I)=0
         CALL BRANG
      END DO
C     ARRAYS FOR ALPHA-2 PEAKS ARE INSERTED
      IF (R12 .GT. 0.0) THEN
         DO J = NREFL, 1, -1
            DO I = 2*J-1, 2*J
               H(I)=H(J)
               K(I)=K(J)
               L(I)=L(J)
               IPAIR(I)=IPAIR(J)
               NOPH(I)=NOPH(J)
               D(I)=D(J)
               RLV2(I)=RLV2(J)
               IF(MOD(I,2).EQ.0) THEN
                  U(I)=U(J)*R12
                  L12(I)=2
                  CALL BRANG
               ELSE
                  U(I)=U(J)
                  L12(I)=0
                  PEAK(I)=PEAK(J)
                  TH(I)=TH(J)
                  FLP(I)=FLP(J)
               ENDIF
            END DO
         END DO
         NREFL=2*NREFL
      ENDIF

*     INITIAL SORTING
      CALL HEAP(PEAK,IPOINT,NREFL)

C     VARIOUS ARRAYS FOR REFLECTIONS ARE SORTED
      CALL ICOPY(H)
      CALL ICOPY(K)
      CALL ICOPY(L)
      CALL ICOPY(IPAIR)
      CALL ICOPY(NOPH)
      CALL ICOPY(L12)
      CALL RCOPY(U)
      CALL RCOPY(D)
      CALL RCOPY(RLV2)
      CALL RCOPY(TH)
      CALL RCOPY(FLP)

      ISUM = 0
*Rev 1.0q 2001.01.30 Izumi {Correction for individual profile fitting
      IF (NREFL .NE. 1) THEN
         FTPEAK1 = 2.0/(PEAK(NREFL) - PEAK(1))
         FTPEAK2 = -0.5*(PEAK(NREFL) + PEAK(1))
         FTTANTH1 = 2.0/(TAN(TH(NREFL)) - TAN(TH(1)))
         FTTANTH2 = -0.5*(TAN(TH(NREFL)) + TAN(TH(1)))
      END IF
* }

      DO I = 1, NREFL
         IPOINT(I)=I
         IF (L12(I) .EQ. 0) THEN
            CALL HRHT
            CALL ASF
         ELSE
            DO J = I-1, 1, -1
               IF (H(I) .EQ. H(J) .AND. K(I) .EQ. K(J) .AND. 
     &         L(I).EQ.L(J) .AND. NOPH(I) .EQ. NOPH(J)) THEN
C                 L12 FOR A K-ALPHA2 PEAK REPRESENTS THE NUMBER FOR
C                 THE CORRESPONDING K-ALPHA1 PEAK
                  L12(I)=J
                  DO JJ=1,NSYM(NOPH(J))
                     RX(JJ,I)=RX(JJ,J)
                     RY(JJ,I)=RY(JJ,J)
                     RZ(JJ,I)=RZ(JJ,J)
                     HT(JJ,I)=HT(JJ,J)
                  END DO
                  DO JJ = 1, NATOM
                     F(JJ,I)=F(JJ,J)
                     F0(JJ,I)=F0(JJ,J)
                  END DO
                  EXIT
               END IF
            END DO
         ENDIF

*        LPPP(I) = 0: THE PROFILE OF THE I'TH REFLECTION IS CALCULATED
*                     FROM SECONDARY PROFILE PARAMETERS.
*        LPPP(I) = n: THE PROFILE OF THE I'TH REFLECTION IS, AT LEAST,
*                     PARTIALLY CALCULATED FROM THE n'TH GROUP OF 
*                     PRIMARY PROFILE PARAMETERS GIVEN BY THE USER
*                     (AND REFINED) SEPARATELY.  IN THIS CASE,THE 
*                     PRIMARY PROFILE PARAMETERS OF THE I'TH REFLECTION
*                     ARE REFINED INSTEAD OF ITS SECONDARY PROFILE 
*                     PARAMETERS.
         LPPP(I) = 0
         DO J = 1, NUCREF
            IF (NOPH(I) .EQ. IPHASE(J) .AND. H(I) .EQ. IHKL(1,J) .AND. 
     &      K(I) .EQ. IHKL(2,J) .AND. L(I) .EQ. IHKL(3,J)) THEN
               LPPP(I) = J
               IF (L12(I) .EQ. 0) ISUM = ISUM + 1
               EXIT
            END IF
         END DO
*        NORMALIZED 2-THETA AND TAN(THETA)
*Rev 1.0q 2001.01.30 Izumi {Correction for individual profile fitting
         IF (NREFL .GT. 1) THEN
            PEAKNOR(I) = FTPEAK1*(PEAK(I) + FTPEAK2)
            TANTHNOR(I) = FTTANTH1*(TAN(TH(I)) + FTTANTH2)
         ELSE
            PEAKNOR(I) = 0.0
            TANTHNOR(I) = 0.0
         END IF
* }
      END DO
      IF (ISUM .NE. NUCREF) CALL JOBEND('Invalid hkl have been input for
     & a relaxed reflection')
*     ASSIGN FMEM(REAL) AND FMEM(IMAGINARY) TO ALL THE REFLECTIONS
      IF (NMODE .EQ. 2 .OR. NMODE .EQ. 3) THEN
         CALL ASFMEM
      ELSE IF (NMODE .EQ. 4 .AND. NSFF .EQ. 1) THEN
         CALL ASSFF
      END IF
      END

************************************************************************

      SUBROUTINE HEAP(X,IPOINT,NX)
C     HEAP SORT
      DIMENSION X(*),IPOINT(*)
      KC=0
      KX=0
      KW=0
      DO 10 J=1,NX
      IPOINT(J)=J
   10 CONTINUE
      DO 400 NS=2,NX
      N=NS
      NW=N
      W=X(N)
      IW=IPOINT(N)
  200 CONTINUE
      ND=N
      N=N/2
      IF(N.EQ.0) GO TO 350
      KC=KC+1
      IF(W.LE.X(N)) GO TO 350
      KX=KX+1
      X(ND)=X(N)
      IPOINT(ND)=IPOINT(N)
      GO TO 200
  350 CONTINUE
      IF(NW.EQ.ND) GO TO 400
      KW=KW+1
      X(ND)=W
      IPOINT(ND)=IW
  400 CONTINUE
      DO 800 NS=NX,2,-1
      NW=NS
      KW=KW+2
      W=X(NS)
      IW=IPOINT(NS)
      X(NS)=X(1)
      IPOINT(NS)=IPOINT(1)
      N=1
  600 CONTINUE
      NH=N
      N=N*2
      IF(N.GE.NS) GO TO 750
      IF(N+1.GE.NS) GO TO 630
      KC=KC+1
      IF(X(N+1).GT.X(N)) N=N+1
  630 KC=KC+1
      IF(W.GT.X(N)) GO TO 750
      KX=KX+1
      X(NH)=X(N)
      IPOINT(NH)=IPOINT(N)
      GO TO 600
  750 CONTINUE
      IF(NW.EQ.NH) GO TO 800
      KW=KW+1
      X(NH)=W
      IPOINT(NH)=IW
  800 CONTINUE
      END

************************************************************************

      SUBROUTINE ICOPY(IX)
      DIMENSION IX(*)
      PARAMETER (NB=7000)
      COMMON /B/ U(NB),NREFL
      COMMON /BUFFER/ IBUF(NB)
      COMMON /ORDER/ IPOINT(NB)
      DO 10 I=1,NREFL
         IBUF(I)=IX(IPOINT(I))
   10 CONTINUE
      DO 20 I=1,NREFL
         IX(I)=IBUF(I)
   20 CONTINUE
      END

************************************************************************

      SUBROUTINE RCOPY(RX)
      DIMENSION RX(*)
      PARAMETER (NB=7000)
      COMMON /B/ U(NB),NREFL
      COMMON /BUFFER/ RBUF(NB)
      COMMON /ORDER/ IPOINT(NB)
      DO 10 I=1,NREFL
         RBUF(I)=RX(IPOINT(I))
   10 CONTINUE
      DO 20 I=1,NREFL
         RX(I)=RBUF(I)
   20 CONTINUE
      END

************************************************************************

      SUBROUTINE ASFMEM
*     ASSIGN F(MEM) OBTAINED WITH MEED TO ALL THE REFLECTIONS
      PARAMETER (NB=7000,NPH=8)
      INTEGER H
*     ARRAYS AP AND BP ARE USED AS BUFFERS
      COMMON /ABPI/ AP(NB),BP(NB)
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /B/ U(NB),NREFL
*Rev 1.0h 2000.12.13 Izumi
      COMMON /BIJVOET/ IPAIR(NB),LPAIR(NPH)
      COMMON /FMEM/ MH(NB),MK(NB),ML(NB),FMEMR(NB),FMEMI(NB),NHKL
      COMMON /FRDA/ NFR,NMEM,MEED(15),EPSD,SCIO,TSCAT,TSCAT1,TSCAT2,
     &  UCLAG,NDA,NFPX(3)
      COMMON /H/ NQ(NPH),NCENTR(NPH),LAUEG(NPH)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
     
      DO JJ = 1, NHKL
         AP(JJ) = FMEMR(JJ)
         BP(JJ) = FMEMI(JJ)
      END DO

      DO J = 1, NREFL
         IF (IPAIR(J) .EQ. 0) THEN
            DO JJ = 1, NHKL
               IF (H(J) .EQ. MH(JJ) .AND. K(J) .EQ. MK(JJ) .AND.
     &         L(J) .EQ. ML(JJ) .AND. NOPH(J) .EQ. 1) THEN
*              Changing the above line as 
*    &         L(J) .EQ. ML(JJ) .AND. NOPH(J) .EQ. NMEM) THEN
*              failed because NMEM is not read in at this point.
                  FMEMR(J) = AP(JJ)
                  FMEMI(J) = BP(JJ)
                  EXIT
               END IF
            END DO
         ELSE 
            DO JJ = 1, NHKL
               IF (H(J) .EQ. -MH(JJ) .AND. K(J) .EQ. -MK(JJ) .AND.
     &         L(J) .EQ. -ML(JJ) .AND. NOPH(J) .EQ. 1) THEN
*              Changing the above line as 
*    &         L(J) .EQ. -ML(JJ) .AND. NOPH(J) .EQ. NMEM) THEN
*              failed because NMEM is not read in at this point.
                  FMEMR(J) =  AP(JJ)
                  FMEMI(J) = -BP(JJ)
                  EXIT
               END IF
            END DO
         END IF
      END DO
      END

************************************************************************

      SUBROUTINE ASSFF
*     ASSIGN (SCALE FACTOR) X |F|**2 TO ALL THE REFLECTIONS OF THE 1ST
*     PHASE
      PARAMETER (NB=7000,NPH=8)
      INTEGER H
*     ARRAYS AP AND BP ARE USED AS BUFFERS
      COMMON /ABPI/ AP(NB),BP(NB)
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /B/ U(NB),NREFL
      COMMON /FMEM/ MH(NB),MK(NB),ML(NB),FMEMR(NB),FMEMI(NB),NHKL
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /T/ FF(NB),COEF(NPH)
     
      DO J = 1, NREFL
         DO JJ = 1, NHKL
            IF (H(J) .EQ. MH(JJ) .AND. K(J) .EQ. MK(JJ) .AND.
     &      L(J) .EQ. ML(JJ) .AND. NOPH(J) .EQ. 1) THEN
               FF(J) = AP(JJ)
               EXIT
            END IF
         END DO
      END DO
      END

************************************************************************

      SUBROUTINE INITINT(A,DEG,XINT,NSTEP)
*     DETERMINE INITIAL VALUES OF |F|**2
      PARAMETER (NP=80000,NB=7000,NPH=8,NAP=150)
      INTEGER IOPT(3)
      REAL EOPT(1)
      REAL A(*),DEG(*),XINT(*)
      COMMON /A/ IH(NB),IK(NB),IL(NB),L12(NB),NOPH(NB)
      COMMON /B/ U(NB),NREFL
      COMMON /C/ BG(NP),BGINC(NP),DEGNOR(NP)
      COMMON /G/ I
      COMMON /I/ RDEG(2,NB)
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /ORDER/ IPOINT(NB)
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /PROFI/ PCOT(NB),SHFT(NB),NTSIM(NB),FNORM(NB),TNTSIM(NB),
     &  FWHM(NB),ETA(NB),SQSG(NB),CTH(NB),SIGP(NB),DFWDG(NB),DEDFFG(NB),
     &  DFWDL(NB),DEDFFL(NB),SIGPAP(NB)
      COMMON /RLS/ NRANGE
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
      COMMON /YOKO/ ABSCISSA(NP),NPOINTS,AWRY
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      COMMON /ZSPFN/ PEAKSH(NB),SQH2(NB),SQH2L(NB),SQH2G(NB),ASYM(NB),
     &ETAL(NB),ETAH(NB),FNORM1(2,NB),COMLG(2,NB),COML(2,NB),COMG(2,NB),
     &COEFL(2,NB),COEFG(2,NB),DECAY(2,NB),COMPEAR(2,NB),CDX(2,NB),
     &CDA1(2,NB),CDA2(2,NB),CDW(2,NB),CDRL1(2,NB),CDRL2(2,NB),
     &CDRH1(2,NB),DMPART(2,NB)
      DATA TOL/8.726646E-06/ ! 0.0005 degrees
*     TWO ARRAYS USED IN SILUP
      DATA NDEG/2/
      DATA LUP/0/
      DATA IOPT/0,1,0/
      DATA RD/0.01745329/
      P(J)=A(KPHB(NOPH(I))+J)

*     COPY NPTS INTO NPOINTS AND DEG INTO VARINV, WHICH WILL BE USED 
*     IN SUBROUTINE PRFUNC
      NPOINTS = NSTEP
      DO J = 1, NSTEP
         ABSCISSA(J) = DEG(J)
      END DO

      ISTART = 1
      DO J = 1, NREFL
         I = IPOINT(J)
*        SAVE REFLECTION NUMBER I BECAUSE IT MAY BE ALTERED IN
*        LOOP J0
         ISAVE = I
         IF (L12(I) .NE. 0 .OR. NOPH(I) .GT. 1) CYCLE

*        TTHPEAK: PEAK POSITION FOR THE I'TH REFLECTION
         IF (NPRFN .EQ. 0) THEN
            MODE = 0
*           THE PEAK POSITION IS SANDWICHED BETWEEN X AND XORF
*           +/- 0.8 DEGREES
            X = PEAK(I) - SHFT(I) - 0.8*RD
            XORF = PEAK(I) - SHFT(I) + 0.8*RD
            DO WHILE (.TRUE.)
               CALL SFMIN(X, XORF, MODE, TOL)
               IF (MODE .NE. 1) EXIT
               XORF = -PRFUNC(X,A)
            END DO
            IF (MODE .EQ. 4) CALL JOBEND('Error termination in SFMIN')
            TTHPEAK = X
         ELSE
            TTHPEAK = PEAK(I) - PEAKSH(I)
         END IF
         SUMPROF = 0.0
         DO J0 = ISTART, NREFL
            I = IPOINT(J0)
            IF (TTHPEAK .GT. RDEG(2,I)) THEN
               ISTART = J0 + 1
               CYCLE
            ELSE IF (TTHPEAK .LT. RDEG(1,I)) THEN
               EXIT
            END IF
            IF (INCMULT .EQ. 0) THEN
               SUMPROF = SUMPROF + U(I)*PRFUNC(TTHPEAK,A)
            ELSE
               SUMPROF = SUMPROF + PRFUNC(TTHPEAK,A)
            END IF
         END DO

*        RESTORE REFLECTION NUMBER
         I = ISAVE
*        INTERPOLATED OBSERVED INTENSITY AT THE PEAK POSITOIN
*Rev 1.0y 2001.04.28 Izumi {
*        CALL SILUP(TTHPEAK,YOBS,NSTEP,XINT,DEG,NDEG,LUP,IOPT,EOPT)
         CALL SILUP(TTHPEAK,YOBS,NSTEP,DEG,XINT,NDEG,LUP,IOPT,EOPT)
* }
*        INTERPOLATED BACKGROUND INTENSITY AT THE PEAK POSITOIN
         CALL SILUP(TTHPEAK,BACKGD,NSTEP,DEG,BG,NDEG,LUP,IOPT,EOPT)
*        INCMULT = 0: F2 = |Fc|**2
*        INCMULT = 1: F2 = Multiplicity * |Fc|**2
*        Refer to Eq. (3) in A. Altomare et al., J. Appl. Crystallogr., 
*        28 (1995) 842.
         F2 = (YOBS - BACKGD)/(P(0)*FLP(I)*SUMPROF)
         IF (INCMULT .EQ. 0) THEN
            YPEAK(I) = P(0)*FLP(I)*U(I)*F2
         ELSE
            YPEAK(I) = P(0)*FLP(I)*F2
         END IF
         IF (R12 .EQ. 0.0) CYCLE
         DO JJ = I+1, NREFL
            IF (L12(JJ) .EQ. I) EXIT
         END DO
         YPEAK(JJ) = YPEAK(I)*R12*FLP(JJ)/FLP(I)
      END DO
      END

************************************************************************

      REAL FUNCTION PRFUNC(X,A)
*     CALCULATE THE PROFILE FUNCTION AT X
      INTEGER IOPT(2)
      REAL A(*),EOPT(2)
      PARAMETER (NP=80000)
      COMMON /ABSVR/ ABSORP(NP),VARLEN(NP),ROUGH(NP),NSURFR
      COMMON /YOKO/ ABSCISSA(NP),NPOINTS,AWRY

*     INTERPOLATE THE ABSORPTION FACTOR
      LUP = 0      ! Start with a binary search
      IOPT(2) = 0  ! End of the option list
      CALL SILUP(X,VA,NPOINTS,ABSCISSA,ABSORP,3,LUP,IOPT,EOPT)
*     INTERPOLATE THE IRRADIATION WIDTH
      LUP = 0      ! Start with a binary search
      IOPT(2) = 0  ! End of the option list
      CALL SILUP(X,VW,NPOINTS,ABSCISSA,VARLEN,3,LUP,IOPT,EOPT)
      PRFUNC = VA*VW*PRFL(X,A)
      END

************************************************************************

      SUBROUTINE KDRREF
C     GENERATE UNIQUE REFLECTION SETS FOR ALL SPACE GROUPS
      PARAMETER (NB=7000,NGEN=7000,NPH=8)
      INTEGER SYMTAB(679),POINT(2,274),SYMCON(2,109),PTCON(2,48),
     &SCLASS(13),CHK(7),NSG,REFCON(14),NSYST,NTYPE,NSC,DEV,TOT,
     &NCO(32),IEX(32),NCOINC,IC,JND(NGEN),MULT,H,K,L,
     &HMAX,KMAX,LMAX,HMIN,KMIN,NOR,ICODE,IPUNCH,COITAB(5,14),SUB,
     &TEMREF(14),ISET(44),SECSET,REFREQ(14),RRT
      REAL*4 A,B,C,ALPHA,BETA,GAM,LAMBDA,ANGMIN,ANGMAX,DMIN,DMAX,
     &AA,BB,CC,CALPHA,CBETA,CGAMMA,DD(NGEN),
     &ANGLE,PATH,DELAY,TMIN,MIND,MINT,FACTOR
      CHARACTER*8 TITLE(10),CTITLE(2,28),RTITLE(2,14),SYSNAM(2,7)
      INTEGER*2 HH(NGEN),KK(NGEN),LL(NGEN)
      COMMON
     &/INPT/ A,B,C,ALPHA,BETA,GAM,LAMBDA,ANGMAX,ANGMIN
     &/PRT/ TITLE,CTITLE,RTITLE,SYSNAM
     &/SGN/ NSG,NSC
     &/CONDTN/ REFCON,NSYST
     &/TYP/ NTYPE
     &/TUF/ ANGLE,PATH,DELAY,TMIN,MIND
     &/CAL/ DMIN,DMAX,MINT,FACTOR
     &/ER/ TOT
     &/HEAD/ DEV,ICODE,IPUNCH
     &/TAB/ SYMTAB,POINT,SYMCON,PTCON,SCLASS,CHK
     &/COI/ NCO,IEX
     &/IDC/ NCOINC,MULT
     &/HKL/ H,K,L
     &/SORT/ DD,IC
     &/LP/ JND
     &/RR/ HH,KK,LL
     &/CIN/ AA,BB,CC,CALPHA,CBETA,CGAMMA
     &/GPAS/ HMAX,KMAX,LMAX,HMIN,KMIN
     &/CTAB/ COITAB
      COMMON
     &/CNOR/ NOR,SUB
     &/TMP/ TEMREF
     &/IS/ ISET,SECSET
     &/REQ/ REFREQ,RRT
      COMMON /EXC/ LSPSYM(NPH),NSSYM(14,NPH)
      COMMON /PHNO/ IP
*     IF (LSPSYM(IP) .LE. 1) THEN
*        NSC=0
*     ELSE
         NSC=1
*     END IF
      RRT=0
      DEV=2
      NOR=0
      SUB=0
      TOT=0
      CALL INPUT
      IF(TOT.EQ.1) CALL JOBEND(' ')
      CALL COND2
      CALL LIMITS
      CALL GENHKL
      DO 10 J=1,IC
         JND(J)=J
   10 CONTINUE
      CALL OUTPUT(5)
      END

************************************************************************

      SUBROUTINE INPUT
      PARAMETER (NPH=8)
      INTEGER NSG,NSC,DEV,REFCON(14),NSYST,ICODE,IPUNCH,TOT,TEMREF(14),
     &SECSET,RRT,ISET(44),REFREQ(14)
      REAL*4 A,B,C,ALPHA,BETA,GAM,LAMBDA,ANGMAX,ANGMIN,
     &ANGLE,PATH,DELAY,TMIN,MIND
      CHARACTER*8 TITLE(10),CTITLE(2,28),RTITLE(2,14),SYSNAM(2,7)
      COMMON
     &/SGN/ NSG,NSC
     &/INPT/ A,B,C,ALPHA,BETA,GAM,LAMBDA,ANGMAX,ANGMIN
     &/HEAD/ DEV,ICODE,IPUNCH
     &/CONDTN/ REFCON,NSYST
     &/L/ DEG1,DEG2
     &/PRT/ TITLE,CTITLE,RTITLE,SYSNAM
     &/TOF/ ANGLE,PATH,DELAY,TMIN,MIND
     &/ER/ TOT
     &/TMP/ TEMREF
     &/IS/ ISET,SECSET
     &/REQ/ REFREQ,RRT
     &/WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
     &/EXC/ LSPSYM(NPH),NSSYM(14,NPH)
      COMMON /PHNO/ IP
C     CLEAR ARRAYS REFCON, TEMREF, AND REFREQ
      DO I=1,14
         REFCON(I)=0
         TEMREF(I)=0
         REFREQ(I)=0
      END DO
      IF(NSG.GE.1 .AND. NSG.LE.274) GO TO 2
      CALL ERROR(1)
      GO TO 10
    2 IF(NSC.EQ.0) GO TO 1
*     NSSYM HAS BEEN ASSIGNED VALUES REGARDLESS OF LSPSYM
*   2 IF(NSC.EQ.0) GO TO 12
      IF(NSC.EQ.1) GO TO 12
      IF(NSC.EQ.2) GO TO 9
      CALL ERROR(6)
      GO TO 10
   12 DO I=1,14
         REFCON(I)=NSSYM(I,IP)
      END DO
      DO 3 I=1,14
      IF(REFCON(I).GE.0 .AND. REFCON(I).LE.28) GO TO 3
      CALL ERROR(2)
      GO TO 10
    3 CONTINUE
      GO TO 1
    9 CONTINUE
      DO 11 I=1,14
      IF(TEMREF(I).GE.-1 .AND. TEMREF(I).LE.28) GO TO 11
      CALL ERROR(2)
      GO TO 10
   11 CONTINUE
    1 IF(RRT.EQ.0) GO TO 16
      IF(RRT.EQ.1) GO TO 17
      CALL ERROR(8)
      GO TO 10
   17 CONTINUE
   16 CONTINUE
      DEV=1
      IF(DEV.EQ.1) GO TO 7
      IF(DEV.EQ.2) GO TO 8
      CALL ERROR(3)
      GO TO 10
    8 CONTINUE
      GO TO 4
    7 LAMBDA=XLMD
    4 ICODE=1
      IPUNCH=0
      IF(ICODE.GE.1 .AND. ICODE.LE.3) GO TO 5
      CALL ERROR(4)
      GO TO 10
    5 IF(IPUNCH.GE.0 .AND. IPUNCH.LE.3) GO TO 6
      CALL ERROR(5)
      GO TO 10
    6 CONTINUE
      CALL OUTPUT(1)
      IF(SECSET.NE.2) GO TO 15
      DO 13 I=1,44
      IF(ISET(I).EQ.NSG) GO TO 14
   13 CONTINUE
      CALL ERROR(7)
      GO TO 10
   14 NSG=230+I
   15 CONTINUE
   10 END

************************************************************************

      SUBROUTINE COND2
      INTEGER SYMTAB(679),POINT(2,274),SYMCON(2,109),PTCON(2,48),
     &SCLASS(13),CHK(7),NSG,REFCON(14),NSYST,NTYPE,NSGC,CLASS,
     &COND(2,30),TEMP(7),SYMNUM,NSC,NCO(32),IEX(32),NCOINC,MULT,
     &TEMREF(14)
      COMMON
     &/TAB/ SYMTAB,POINT,SYMCON,PTCON,SCLASS,CHK
     &/SGN/ NSG,NSC
     &/CONDTN/ REFCON,NSYST
     &/TYP/ NTYPE
     &/COI/ NCO,IEX
     &/IDC/ NCOINC,MULT
     &/TMP/ TEMREF
      DO 1 I=1,13
      IF(NSG.LT.SCLASS(I)) GO TO 2
    1 CONTINUE
    2 CLASS=I-1
      NSGC=(NSG-SCLASS(I-1))+1
      NSGC=NSGC
      IF(NSC.EQ.1) GO TO 9
    5 DO 20 I=1,30
      DO 21 J=1,2
      COND(J,I)=0
   21 CONTINUE
   20 CONTINUE
      K=POINT(1,NSG)
      L=POINT(2,NSG)
      DO 6 I=1,L
      TEMP(I)=SYMTAB(K)
      K=K+1
    6 CONTINUE
      J=1
      DO 7 I=1,L
      SYMNUM=TEMP(I)
      MM=PTCON(1,SYMNUM)
      N=PTCON(2,SYMNUM)
      DO 8 II=1,N
      COND(1,J)=SYMCON(1,MM)
      COND(2,J)=SYMCON(2,MM)
      MM=MM+1
      J=J+1
    8 CONTINUE
    7 CONTINUE
      DO 16 I=1,14
      REFCON(I)=0
   16 CONTINUE
      DO 15 I=1,30
      IF(COND(2,I).EQ.0) GO TO 10
      REFCON(COND(2,I))=COND(1,I)
   15 CONTINUE
   10 IF(NSC.EQ.0) GO TO 9
    4 DO 3 I=1,14
      IF(TEMREF(I).EQ.-1) GO TO 3
      REFCON(I)=TEMREF(I)
    3 CONTINUE
    9 IF(CLASS.EQ.5) GO TO 105
      IF(CLASS.LE.7) GO TO 205
      IF(CLASS.GE.8 .AND. CLASS.LE.10) GO TO 305
      IF(CLASS.GE.11) GO TO 405
  105 DO 100 I=1,7
      IF(NSG.EQ.CHK(I)) GO TO 110
  100 CONTINUE
      NSYST=6
      GO TO 120
  110 NSYST=5
      GO TO 120
  205 NSYST=CLASS
      GO TO 120
  305 NSYST=CLASS-6
      GO TO 120
  405 NSYST=CLASS-5
  120 CONTINUE
      IF(NSYST.LE.4) GO TO 200
      IF(NSYST.EQ.5) GO TO 210
      IF(NSYST.GE.6) GO TO 220
  200 NTYPE=NSYST
      GO TO 240
  210 NTYPE=1
      GO TO 240
  220 NTYPE=NSYST-2
  240 CONTINUE
   17 CALL OUTPUT(2)
      DO 40 I=1,32
      IF(NSG.LE.NCO(I)) GO TO 41
   40 CONTINUE
   41 NCOINC=IEX(I)
      IF(NCOINC.EQ.0) GO TO 42
      CALL OUTPUT(4)
   42 CONTINUE
   18 CONTINUE
      END

************************************************************************

      SUBROUTINE ERROR(ERR)
      INTEGER ERR,TOT
      COMMON
     &/ER/ TOT
      TOT=1
      GO TO (1,2,3,4,5,6,7,8), ERR
      WRITE(6,100)
      GO TO 20
    1 WRITE(6,11)
      GO TO 20
    2 WRITE(6,12)
      GO TO 20
    3 WRITE(6,13)
      GO TO 20
    4 WRITE(6,14)
      GO TO 20
    5 WRITE(6,15)
      GO TO 20
    6 WRITE(6,16)
      GO TO 20
    7 WRITE(6,17)
      GO TO 20
    8 WRITE(6,18)
   11 FORMAT(//' ',10X,'Illegal space group No.')
   12 FORMAT(//' ',10X,'Illegal conditions')
   13 FORMAT(//' ',10X,'Illegal data type')
   14 FORMAT(//' ',10X,'Illegal nuclear/magnetic code')
   15 FORMAT(//' ',10X,'Illegal output code')
   16 FORMAT(//' ',10X,'Illegal test for non-standard conditions')
   17 FORMAT(//' ',10X,'There is no alternative setting for this space g
     &roup.'/' ',10X,'Job abandoned.')
   18 FORMAT(//' ',10X,'Illegal reflection type condition')
  100 FORMAT(//' ',10X,'Error in system')
   20 END

************************************************************************

      SUBROUTINE LIMITS
      INTEGER NTYPE,HMAX,KMAX,LMAX,HMIN,KMIN,H,K,L,DEV,ICODE,IPUNCH
      REAL*4 A,B,C,ALPHA,BETA,GAM,LAMBDA,ANGMAX,ANGMIN,DMIN,DMAX,RAD,S,
     &VOL,AA,BB,CC,CALPHA,CBETA,CGAMMA,FACTOR,PATH,ANGLE,TMIN,DELAY,
     &MIND,MINT
      CHARACTER*8 TITLE(10),CTITLE(2,28),RTITLE(2,14),SYSNAM(2,7)
      COMMON
     &/INPT/ A,B,C,ALPHA,BETA,GAM,LAMBDA,ANGMAX,ANGMIN
     &/TYP/ NTYPE
     &/TOF/ ANGLE,PATH,DELAY,TMIN,MIND
     &/CAL/ DMIN,DMAX,MINT,FACTOR
     &/HEAD/ DEV,ICODE,IPUNCH
     &/PRT/ TITLE,CTITLE,RTITLE,SYSNAM
     &/HKL/ H,K,L
     &/CIN/ AA,BB,CC,CALPHA,CBETA,CGAMMA
     &/GPAS/ HMAX,KMAX,LMAX,HMIN,KMIN
      DATA RAD/0.0174532925/
      ALPHA=ALPHA*RAD
      BETA=BETA*RAD
      GAM=GAM*RAD
      S=(ALPHA+BETA+GAM)/2
      VOL=2*A*B*C*SQRT(SIN(S)*SIN(S-ALPHA)*SIN(S-BETA)*SIN(S-GAM))
      AA=B*C*SIN(ALPHA)/VOL
      BB=C*A*SIN(BETA)/VOL
      CC=A*B*SIN(GAM)/VOL
      CALPHA=(COS(BETA)*COS(GAM)-COS(ALPHA))/(SIN(BETA)*SIN(GAM))
      CBETA=(COS(GAM)*COS(ALPHA)-COS(BETA))/(SIN(GAM)*SIN(ALPHA))
      CGAMMA=(COS(ALPHA)*COS(BETA)-COS(GAM))/(SIN(ALPHA)*SIN(BETA))
      IF(ABS(CALPHA).LT.0.001) CALPHA=0.
      IF(ABS(CBETA).LT.0.001) CBETA=0.
      IF(ABS(CGAMMA).LT.0.001) CGAMMA=0.
      IF(DEV.EQ.1) GO TO 1
      FACTOR=2*252.77*PATH*SIN(0.5*ANGLE*RAD)
      DMIN=(TMIN+DELAY)/FACTOR
      IF(DMIN.LT.MIND) DMIN=MIND
      MINT=DMIN*FACTOR-DELAY
      DMAX=1000
      GO TO 2
    1 DMIN=LAMBDA/(2.0*SIN(0.5*ANGMAX*RAD))
      DMAX=LAMBDA/(2.0*SIN(0.5*ANGMIN*RAD))
    2 CALL OUTPUT(3)
      HMAX=INT(A/DMIN)
      KMAX=INT(B/DMIN)
      LMAX=INT(C/DMIN)
      HMIN=0
      KMIN=0
      H=1
      K=0
      L=0
      IF(NTYPE-2) 5,6,10
    5 KMIN=-KMAX
    6 HMIN=-HMAX
   10 CONTINUE
      END

************************************************************************

      SUBROUTINE GENHKL
      PARAMETER (NB=7000,NGEN=7000)
      INTEGER IC,HMIN,MAXH,MINH,KMAX,KMIN,LMAX,NTYPE,REFCON(14),NSYST,
     &H,K,L,HMAX,IH,IK,IL,RT,NSG,NSC,REFREQ(14),RRT
      INTEGER*2 HH(NGEN),KK(NGEN),LL(NGEN)
      REAL*4 QMAX,DMIN,DMAX,AA,BB,CC,CALPHA,CBETA,CGAMMA,HA,HB,HC,ZAP,
     &AH,BK,CL,Q,D,DD(NGEN),MINT,FACTOR,MINM,MAXM
      COMMON
     &/SGN/ NSG,NSC
     &/CAL/ DMIN,DMAX,MINT,FACTOR
     &/TYP/ NTYPE
     &/CONDTN/ REFCON,NSYST
     &/HKL/ H,K,L
     &/CIN/ AA,BB,CC,CALPHA,CBETA,CGAMMA
     &/SORT/ DD,IC
     &/RR/ HH,KK,LL
     &/GPAS/ HMAX,KMAX,LMAX,HMIN,KMIN
     &/REQ/ REFREQ,RRT
      SAVE NTOTAL
      DATA NTOTAL/0/

*     A change to adapt for Hitachi's case    
      init = 1
      IC=1
  194 QMAX=1/DMIN
      HA=AA**2
      HB=2*AA*(K*BB*CGAMMA+L*CC*CBETA)
      HC=(K*BB)**2+(L*CC)**2+2*K*L*BB*CC*CALPHA-QMAX**2
      ZAP=(HB**2-4*HA*HC)
      IF(ZAP.LT.0.0) GO TO 116
 2750 MAXM=(-HB+SQRT(ZAP))/(2*HA)
      MINM=(-HB-SQRT(ZAP))/(2*HA)
      IF(MAXM.LT.0) MAXM=MAXM-1
      IF(MINM.GT.0) MINM=MINM+1
      MAXH=MAXM
      MINH=MINM
*     A change to adapt for Hitachi's case    
*     IF(H.GT.MAXH) GO TO 117
      if (h .gt. maxh .and. init .eq. 0) then
         go to 117
      end if
      init = 0
 2950 IF(NTYPE.GT.2) GO TO 206
      IF(K.EQ.0 .AND. L.EQ.0) GO TO 206
  205 H=MINH
  206 CONTINUE
      IF(NSYST.NE.2) GO TO 208
      IF(NSG.GT.230) GO TO 207
      IF(K.EQ.0 .AND. H.LT.0) GO TO 122
      GO TO 208
  207 IF(L.EQ.0 .AND. H.LT.0) GO TO 122
  208 IF(NSYST.NE.7) GO TO 40
      IH=H
      IK=K
      IL=L
      H=IABS(IH)
      K=IABS(IK)
      L=IABS(IL)

*     CHECK ALL THE SPECIAL SUBSETS ------------------------------------
   40 DO 120 LOOP=1,14
      RT=LOOP
      IF (LOOP .EQ. 1 .AND. H .EQ. 0 .AND. K .EQ. 0) THEN
      ELSE IF (LOOP .EQ. 2 .AND. H .EQ. 0 .AND. L .EQ. 0) THEN
      ELSE IF (LOOP .EQ. 3 .AND. H .EQ. 0) THEN
      ELSE IF (LOOP .EQ. 4 .AND. K .EQ. 0 .AND. L .EQ. 0) THEN
      ELSE IF (LOOP .EQ. 5 .AND. K .EQ. 0) THEN
      ELSE IF (LOOP .EQ. 6 .AND. L .EQ. 0) THEN
      ELSE IF (LOOP .EQ. 7) THEN
      ELSE IF (LOOP .EQ. 8 .AND. H .EQ. 0 .AND. K .EQ. L) THEN
      ELSE IF (LOOP .EQ. 9 .AND. H .EQ. K .AND. L .EQ. 0) THEN
      ELSE IF (LOOP .EQ. 10 .AND. H .EQ. K) THEN
      ELSE IF (LOOP .EQ. 11 .AND. H .EQ. L .AND. K .EQ. 0) THEN
      ELSE IF (LOOP .EQ. 12 .AND. K .EQ. L) THEN
      ELSE IF (LOOP .EQ. 13 .AND. H .EQ. L) THEN
      ELSE IF (LOOP .EQ. 14 .AND. H .EQ. K .AND. K .EQ. L) THEN
      ELSE
         GO TO 120
      END IF
C     IF(H.EQ.K .AND. H.EQ.L) GO TO 41
C     RT=0
C     IF(L.NE.0) RT=RT+1
C     IF(K.NE.0) RT=RT+2
C     IF(H.NE.0) RT=RT+4
C     IF(H.EQ.K .AND. H.NE.0) RT=RT+3
C     IF(H.EQ.L .AND. H.NE.0) RT=RT+6
C     IF(L.EQ.K .AND. L.NE.0) RT=RT+5
C     GO TO 42
C  41 RT=14
   42 CONTINUE
      IF(NSYST.NE.7) GO TO 43
      H=IH
      K=IK
      L=IL
   43 CONTINUE
      J=RT
      IF(REFREQ(J).NE.0) GO TO 122
      J=7
   44 II=REFCON(J)
      IF(II.EQ.0) GO TO 51
   31 GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     &23,24,25,26,27,28),  II
    1 IF(MOD(H,2).EQ.0) GO TO 51
      GO TO 122
    2 IF(MOD(K,2).EQ.0) GO TO 51
      GO TO 122
    3 IF(MOD(L,2).EQ.0) GO TO 51
      GO TO 122
    4 IF(MOD(K+L,2).EQ.0) GO TO 51
      GO TO 122
    5 IF(MOD(H+L,2).EQ.0) GO TO 51
      GO TO 122
    6 IF(MOD(H+K,2).EQ.0) GO TO 51
      GO TO 122
    7 IF(MOD(H,2).EQ.0 .AND. MOD(K,2).EQ.0 .AND. MOD(L,2).EQ.0 .OR.
     &MOD(H,2).NE.0 .AND. MOD(K,2).NE.0 .AND. MOD(L,2).NE.0) GO TO 51
      GO TO 122
    8 IF(MOD(K+L,4).EQ.0) GO TO 51
      GO TO 122
    9 IF(MOD(H+L,4).EQ.0) GO TO 51
      GO TO 122
   10 IF(MOD(H+K,4).EQ.0) GO TO 51
      GO TO 122
   11 IF(MOD(2*H+L,2).EQ.0) GO TO 51
      GO TO 122
   12 IF(MOD(2*H+L,4).EQ.0) GO TO 51
      GO TO 122
   13 IF(MOD(H+K+L,2).EQ.0) GO TO 51
      GO TO 122
   14 IF(MOD(-H+K+L,3).EQ.0) GO TO 51
      IF(NSYST.EQ.6) GO TO 51
      GO TO 122
   15 IF(MOD(H-K+L,3).EQ.0) GO TO 51
      IF(NSYST.EQ.6) GO TO 51
      GO TO 122
   16 IF(MOD(H,4).EQ.0) GO TO 51
      GO TO 122
   17 IF(MOD(K,4).EQ.0) GO TO 51
      GO TO 122
   18 IF(MOD(L,3).EQ.0) GO TO 51
      GO TO 122
   19 IF(MOD(L,4).EQ.0) GO TO 51
      GO TO 122
   20 IF(MOD(L,6).EQ.0) GO TO 51
      GO TO 122
   21 IF(IABS(H).GE.IABS(K) .AND. IABS(K).GE.IABS(L)) GO TO 29
      GO TO 122
   22 IF(MOD(2*H+K,2).EQ.0) GO TO 51
      GO TO 122
   23 IF(MOD(2*H+K,4).EQ.0) GO TO 51
      GO TO 122
   24 IF(MOD(H+2*K,2).EQ.0) GO TO 51
      GO TO 122
   25 IF(MOD(H+2*K,4).EQ.0) GO TO 51
      GO TO 122
   26 IF(MOD(H,2).EQ.0 .AND. MOD(K,2).EQ.0) GO TO 51
      GO TO 122
   27 IF(MOD(K,2).EQ.0 .AND. MOD(L,2).EQ.0) GO TO 51
      GO TO 122
   28 IF(MOD(H,2).EQ.0 .AND. MOD(L,2).EQ.0) GO TO 51
      GO TO 122
   29 IF(K+L.EQ.0 .AND. H.LT.0) GO TO 122
      IF(H+K.EQ.0 .AND. K.LT.0) GO TO 122
   51 IF(J.NE.7) GO TO 120
      IF(RT.EQ.7) GO TO 120
      J=RT
      GO TO 44
  120 CONTINUE
*     ------------------------------------------------------------------
  119 AH=AA*H
      BK=BB*K
      CL=CC*L
      Q=SQRT(AH*AH+BK*BK+CL*CL+2*AH*BK*CGAMMA+2*BK*CL*CALPHA
     &+2*CL*AH*CBETA)
      D=1./Q
      IF (D .GT. DMAX .OR. D .LT. DMIN) GO TO 122
  121 IF(NTOTAL+IC .GT. NGEN) THEN
         CALL JOBEND('The value of NGEN (NB) is too small')
      ENDIF
      HH(IC)=H
      KK(IC)=K
      LL(IC)=L
      DD(IC)=D
      IC=IC+1
  122 H=H+1
      IF(MAXH-H) 116,206,206
  116 K=K+1
      IF(KMAX-K) 117,114,114
  117 L=L+1
      IF(LMAX-L) 196,113,113
  113 IF(NTYPE.EQ.5) GO TO 108
  107 K=KMIN
  114 IF(NTYPE.EQ.4 .OR. NTYPE.EQ.5) GO TO 109
      H=HMIN
      GO TO 110
  108 K=L
  109 H=K
  110 CONTINUE
      GO TO 194
  196 CONTINUE
      IC=IC-1
      NTOTAL=NTOTAL+IC
      END

************************************************************************

      SUBROUTINE OUTPUT(IND)
      PARAMETER (NB=7000,NGEN=7000,NT=999,NPH=8)
      INTEGER NSG,REFCON(14),NSYST,IND,NSC,DEV,NCOINC,IC,
     &JND(NGEN),MULT,NOR,ICODE,IPUNCH,SUB,H,K,L,SECSET,ISET(44),
     &RRT,REFREQ(14)
      REAL*4 A,B,C,ALPHA,BETA,GAM,LAMBDA,ANGMAX,ANGMIN,DMIN,DMAX,
     &DD(NGEN),X,STHE,ATHE2,RAD,RH,RK,RL,ANGLE,PATH,MINT,TIME,
     &DELAY,TMIN,MIND,FACTOR,AA,BB,CC,CALPHA,CBETA,CGAMMA
      CHARACTER*8 TITLE(10),CTITLE(2,28),RTITLE(2,14),SYSNAM(2,7)
      CHARACTER PARNAM*48,PHNAME*25
      INTEGER*2 HH(NGEN),KK(NGEN),LL(NGEN)
      LOGICAL SET(28)
      COMMON
     &/A/ IH(NB),IK(NB),IL(NB),L12(NB),NOPH(NB)
     &/B/ U(NB),NREFL
     &/INPT/ A,B,C,ALPHA,BETA,GAM,LAMBDA,ANGMAX,ANGMIN
     &/EXC/ LSPSYM(NPH),NSSYM(14,NPH)
     &/SGN/ NSG,NSC
     &/CONDTN/ REFCON,NSYST
     &/TOF/ ANGLE,PATH,DELAY,TMIN,MIND
     &/CAL/ DMIN,DMAX,MINT,FACTOR
     &/PRT/ TITLE,CTITLE,RTITLE,SYSNAM
     &/HEAD/ DEV,ICODE,IPUNCH
      COMMON
     &/IDC/ NCOINC,MULT
     &/SORT/ DD,IC
     &/LP/ JND
     &/RR/ HH,KK,LL
     &/CNOR/ NOR,SUB
     &/CIN/ AA,BB,CC,CALPHA,CBETA,CGAMMA
     &/IS/ ISET,SECSET
     &/REFNO/ IREFL
     &/REQ/ REFREQ,RRT
     &/PRLV/ NPRINT
     &/PHNO/ IP
     &/PAR/ PARNAM(NT),PHNAME(NPH)
      COMMON /CONHKL/ NHKL(NPH),ICH(5,3,7,NPH),NCON(7,NPH),NLINE(NPH)
      DATA RAD/0.0174532925/

      IREFL=NREFL+1
      IF (IREFL .GT. NB) CALL JOBEND('The total number of reflections'//
     &' has reached the maximum number')  
      GO TO (21,22,23,28,29), IND
   21 CONTINUE
   56 CONTINUE
      GO TO 24
   22 IF(NPRINT.EQ.2) WRITE(6,3007) PHNAME(IP)
 3007 FORMAT(//11X,'Symmetry conditions in ',A)
      DO 1 I=1,28
      SET(I)=.FALSE.
    1 CONTINUE
      IF(REFCON(7).EQ.0 .OR. REFCON(7).EQ.21) GO TO 53
      GO TO 52
   53 IF(NPRINT.EQ.2) WRITE(6,205)
      GO TO 55
   52 DO 50 J=1,28
      IF(REFCON(7).EQ.J) GO TO 51
      GO TO 50
   51 IF(NPRINT.EQ.2) WRITE(6,203) CTITLE(1,J),CTITLE(2,J)
   50 CONTINUE
      IF(NPRINT.EQ.2) WRITE(6,206) RTITLE(1,7),RTITLE(2,7)
   55 DO 100 J=1,28
      DO 200 I=1,14
      IF(I.EQ.7) GO TO 200
      IF(REFCON(I).EQ.J) GO TO 201
      GO TO 200
  201 IF(SET(J)) GO TO 202
      SET(J)=.TRUE.
      IF(NPRINT.EQ.2) THEN
         WRITE(6,203) CTITLE(1,J),CTITLE(2,J)
         WRITE(6,206) RTITLE(1,I),RTITLE(2,I)
      ENDIF
      GO TO 200
  202 IF(NPRINT.EQ.2) WRITE(6,204) RTITLE(1,I),RTITLE(2,I)
  200 CONTINUE
  100 CONTINUE
  203 FORMAT(' ',12X,2A8)
  204 FORMAT(' ',30X,2A8)
  205 FORMAT(' ',12X,'No lattice absences')
  206 FORMAT('+',30X,2A8)
      GO TO 24
   23 CONTINUE
      GO TO (25,26), DEV
   25 CONTINUE
      GO TO 27
   26 CONTINUE
   27 CONTINUE
      GO TO 24
   28 CONTINUE
      GO TO 24
   29 IF(DEV.EQ.1) GO TO 426
      DO 424 I=1,IC
      J=JND(I)
      TIME=DD(I)*FACTOR-DELAY
      TIME=TIME
      RH=HH(J)
      RK=KK(J)
      RL=LL(J)
      CALL SEMULT(RH,RK,RL,NSG,MULT)
      H=HH(J)
      K=KK(J)
      L=LL(J)
      IF(NSYST.EQ.6 .AND. (REFCON(7).EQ.14 .OR. REFCON(7).EQ.15)) GOTO 5
      GO TO 3
    5 IF(REFCON(7).EQ.14) GO TO 2
    7 IF(MOD(H-K+L,3).EQ.0) GO TO 3
      SUB=SUB+1
      GO TO 4
    2 IF(MOD(-H+K+L,3).EQ.0) GO TO 3
      SUB=SUB+1
      GO TO 4
    3 IF(LSPSYM(IP).NE.1 .OR. IPOS(H,K,L,NHKL(IP),ICH(1,1,1,IP),
     &NCON(1,IP),NLINE(IP)) .EQ. 1) THEN
         IH(IREFL)=H
         IK(IREFL)=K
         IL(IREFL)=L
         U(IREFL)=MULT
         NOPH(IREFL)=IP
         IREFL=IREFL+1
         IF (IREFL .GT. NB) CALL JOBEND('The total number of '//
     &   'reflections has reached the maximum number')  
      ENDIF
      IF(IPUNCH.EQ.0) GO TO 4
      IF(IPUNCH.EQ.2) GO TO 425
      IF(IPUNCH.EQ.1) GO TO 4
  425 CONTINUE
    4 IF(H+K.EQ.0 .AND. NCOINC.NE.2) GO TO 31
      IF(NCOINC.EQ.0) GO TO 31
      CALL COINC(J)
   31 CONTINUE
  424 CONTINUE
      GO TO 427
  426 CONTINUE
      DO 421 I=1,IC
      J=JND(I)
      X=1./(2.*DD(I))
      STHE=X*LAMBDA
      ATHE2=2.0*ASIN(STHE)/RAD
      ATHE2=ATHE2
      RH=HH(J)
      RK=KK(J)
      RL=LL(J)
      CALL SEMULT(RH,RK,RL,NSG,MULT)
      H=HH(J)
      K=KK(J)
      L=LL(J)
      IF(NSYST.EQ.6 .AND. (REFCON(7).EQ.14 .OR. REFCON(7).EQ.15)) GO TO
     &16
      GO TO 14
   16 IF(REFCON(7).EQ.14) GO TO 13
   18 IF(MOD(H-K+L,3).EQ.0) GO TO 14
      SUB=SUB+1
      GO TO 15
   13 IF(MOD(-H+K+L,3).EQ.0) GO TO 14
      SUB=SUB+1
      GO TO 15
   14 IF(LSPSYM(IP).NE.1 .OR. IPOS(H,K,L,NHKL(IP),ICH(1,1,1,IP),
     &NCON(1,IP),NLINE(IP)) .EQ. 1) THEN
         IH(IREFL)=H
         IK(IREFL)=K
         IL(IREFL)=L
         U(IREFL)=MULT
         NOPH(IREFL)=IP
         IREFL=IREFL+1
         IF (IREFL .GT. NB) CALL JOBEND('The total number of '//
     &   'reflections has reached the maximum number')  
      ENDIF
      IF(IPUNCH.EQ.0) GO TO 15
      IF(IPUNCH.EQ.2) GO TO 422
      IF(IPUNCH.EQ.1) GO TO 15
  422 CONTINUE
   15 IF(H+K.EQ.0 .AND. NCOINC.NE.2) GO TO 30
      IF(NCOINC.EQ.0) GO TO 30
      CALL COINC(J)
   30 CONTINUE
  421 CONTINUE
  427 IF(IPUNCH.EQ.0) GO TO 61
      IF(IPUNCH.EQ.2) GO TO 62
      IF(IPUNCH.EQ.1) GO TO 61
   62 CONTINUE
   61 NREFL=IREFL-1
   24 CONTINUE
   57 END

************************************************************************

      SUBROUTINE SEMULT(H,K,L,ISG,MULT)
      REAL*4 H,K,L,HH,KK,LL
      DIMENSION MSG(19),MULTM(14,14)
      DATA MSG/2,15,74,88,142,148,167,176,194,206,230,243,248,251,261,
     &263,268,270,274/
      DATA MULTM/2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6, 2, 6, 6,
     &           2, 2, 2, 4, 4, 6, 6, 6, 6, 6, 6, 2, 6, 6,
     &           2, 4, 4, 8, 8, 6, 6,12,12,12,24, 4, 6,12,
     &           2, 2, 2, 4, 4, 6, 6, 6, 6, 6, 6, 2, 6, 6,
     &           2, 4, 4, 8, 8, 6, 6,12,12,12,24, 2, 6,12,
     &           2, 2, 4, 4, 8, 6,12, 6,12,12,24, 4, 6,12,
     &           2, 4, 8, 8,16, 6,12,12,24,24,48, 4, 6,12,
     &           2, 4, 4, 8, 8, 6, 6,12,12,12,12, 4, 6, 6,
     &           2, 2, 4, 4, 4, 6, 6, 6, 6,12,12, 4, 6, 6,
     &           2, 4, 8, 8, 8, 6,12,12,12,24,24, 4, 6, 6,
     &           2, 4, 4, 8, 8, 6, 6,12,12,12,12, 2, 6, 6,
     &           2, 4, 8, 8,16, 6,12,12,24,24,24, 4, 6, 6,
     &           2, 4, 8, 8,16, 6,12,12,24,24,24, 4, 6, 6,
     &           2, 4, 8, 8, 8, 6,12,12,12, 8, 8, 4, 2, 2/
      IF(ISG.LT.195) GO TO 50
      IF(ISG.GT.230 .AND. ISG.LT.269) GO TO 50
      HH=H
      KK=K
      LL=L
      H=ABS(HH)
      K=ABS(KK)
      L=ABS(LL)
   50 IF(H.EQ.K .AND. H.EQ.L) GO TO 100
      J=0
      IF(L.NE.0.0) J=J+1
      IF(K.NE.0.0) J=J+2
      IF(H.NE.0.0) J=J+4
      IF(H.EQ.K .AND. H.NE.0.0) J=J+3
      IF(H.EQ.L .AND. H.NE.0.0) J=J+6
      IF(L.EQ.K .AND. L.NE.0.0) J=J+5
      GO TO 90
  100 J=14
   90 IF(ISG.LT.195) GO TO 200
      IF(ISG.GT.230 .AND. ISG.LT.269) GO TO 200
      H=HH
      K=KK
      L=LL
      IF(ISG.GT.270) GO TO 40
      IF(ISG.LT.206 .OR. ISG.GT.268) GO TO 70
   40 I=11
      GO TO 600
   70 I=10
      GO TO 600
  200 IF(ISG.LT.146 .OR. ISG.GT.167) GO TO 500
      NISG=ISG-145
      GO TO (2,4,2,5,1,5,1,5,1,3,1,5,1,5,3,3,5,5,1,1,3,3), NISG
    2 I=13
      GO TO 600
    3 I=14
      IF(H.NE.0.0 .AND. K.NE.0.0 .AND. L.NE.0.0) GO TO 600
      IF(ABS(H+K+L).LT.0.0005) I=13
      GO TO 600
    5 IF(J.EQ.3 .OR. J.EQ.5) GO TO 6
      IF(J.EQ.8 .OR. J.EQ.11) GO TO 6
      IF(J.EQ.10 .OR. J.EQ.14) J=5
      GO TO 1
    6 J=10
    1 I=7
      GO TO 800
    4 CONTINUE
  500 DO 300 I=1,19
      IF(MSG(I).GE.ISG) GO TO 400
  300 CONTINUE
  400 IF(I.GE.18) I=I-8
      IF(I.GE.13) I=I-10
  800 IF(I.NE.7) GO TO 700
      IF(ABS(H+K).LT.0.0005) GO TO 60
      IF(L.NE.0.0) GO TO 600
      IF(ABS(H+2*K).LT.0.0005) I=6
      IF(ABS(2*H+K).LT.0.0005) I=6
      GO TO 600
  700 IF(I.NE.9) GO TO 600
      IF(ABS(H+K).LT.0.0005) GO TO 80
      IF(ABS(H+2*K).LT.0.0005) I=8
      IF(ABS(2*H+K).LT.0.0005) I=8
      GO TO 600
   60 I=6
      GO TO 600
   80 I=8
  600 MULT=MULTM(I,J)
      END

************************************************************************

      SUBROUTINE COINC(J)
      PARAMETER (NB=7000,NGEN=7000)
      INTEGER NCOINC,IC,RT,COITAB(5,14),INDEX,NSG,H,K,L,IH,IK,IL,NSC,
     &NOR,MULT,REFCON(14),NSYST,DEV,ICODE,IPUNCH
      REAL*4 DD(NGEN)
      INTEGER*2 HH(NGEN),KK(NGEN),LL(NGEN)
      CHARACTER*8 TITLE(10),CTITLE(2,28),RTITLE(2,14),SYSNAM(2,7)
      COMMON
     &/IDC/ NCOINC,MULT
     &/RR/ HH,KK,LL
     &/SORT/ DD,IC
     &/CTAB/ COITAB
     &/SGN/ NSG,NSC
     &/CNOR/ NOR,SUB
     &/CONDTN/ REFCON,NSYST
     &/PRT/ TITLE,CTITLE,RTITLE,SYSNAM
     &/HEAD/ DEV,ICODE,IPUNCH
      H=HH(J)
      K=KK(J)
      L=LL(J)
      IF(NSYST.NE.7) GO TO 40
      IH=H
      IK=K
      IL=L
      H=IABS(IH)
      K=IABS(IK)
      L=IABS(IL)
   40 IF(H.EQ.K .AND. H.EQ.L) GO TO 41
      RT=0
      IF(L.NE.0) RT=RT+1
      IF(K.NE.0) RT=RT+2
      IF(H.NE.0) RT=RT+4
      IF(H.EQ.K .AND. H.NE.0) RT=RT+3
      IF(H.EQ.L .AND. H.NE.0) RT=RT+6
      IF(L.EQ.K .AND. L.NE.0) RT=RT+5
      GO TO 42
   41 RT=14
   42 CONTINUE
      IF(NSYST.NE.7) GO TO 43
      H=IH
      K=IK
      L=IL
   43 CONTINUE
      INDEX=COITAB(NCOINC,RT)
      IF(INDEX.EQ.0) GO TO 4
      GO TO (1,2,3), INDEX
    1 IH=K
      IK=H
      CALL COINC2(IH,IK,L)
      GO TO 4
    2 IH=-H
      IK=-K
      CALL COINC2(IH,IK,L)
      GO TO 4
    3 IH=K
      IK=H
      CALL COINC2(IH,IK,L)
      IH=-H
      IK=-K
      CALL COINC2(IH,IK,L)
      IH=-K
      IK=-H
      CALL COINC2(IH,IK,L)
    4 CONTINUE
      END

************************************************************************

      SUBROUTINE COINC2(H,K,L)
      PARAMETER (NB=7000,NGEN=7000,NPH=8)
      INTEGER H,K,L,NCOINC,MULT,DEV,ICODE,IPUNCH,NOR,IC,REFCON(14),
     &NSYST
      REAL*4 DD(NGEN)
      CHARACTER*8 TITLE(10),CTITLE(2,28),RTITLE(2,14),SYSNAM(2,7)
      COMMON
     &/A/ IH(NB),IK(NB),IL(NB),L12(NB),NOPH(NB)
     &/B/ U(NB),NREFL
     &/EXC/ LSPSYM(NPH),NSSYM(14,NPH)
     &/IDC/ NCOINC,MULT
     &/PRT/ TITLE,CTITLE,RTITLE,SYSNAM
     &/HEAD/ DEV,ICODE,IPUNCH
     &/CNOR/ NOR,SUB
     &/SORT/ DD,IC
     &/CONDTN/ REFCON,NSYST
     &/PHNO/ IP
     &/REFNO/ IREFL
      COMMON /CONHKL/ NHKL(NPH),ICH(5,3,7,NPH),NCON(7,NPH),NLINE(NPH)
      IF(NSYST.EQ.6 .AND. (REFCON(7).EQ.14 .OR. REFCON(7).EQ.15)) GO TO
     &20
      GO TO 3
   20 IF(REFCON(7).EQ.14) GO TO 2
      GO TO 7
    2 IF(MOD(-H+K+L,3).EQ.0) GO TO 3
      GO TO 24
    7 IF(MOD(H-K+L,3).EQ.0) GO TO 3
      GO TO 24
    3 IF(LSPSYM(IP).NE.1 .OR. IPOS(H,K,L,NHKL(IP),ICH(1,1,1,IP),
     &NCON(1,IP),NLINE(IP)) .EQ. 1) THEN
         IF (IREFL .GT. NB) CALL JOBEND
     &   ('The number of reflections has reached the maximum value')
         IH(IREFL)=H
         IK(IREFL)=K
         IL(IREFL)=L
         U(IREFL)=MULT
         NOPH(IREFL)=IP
         IREFL=IREFL+1
         IF (IREFL .GT. NB) CALL JOBEND('The total number of '//
     &   'reflections has reached the maximum number')  
         NOR=NOR+1
      ENDIF
      IF(IPUNCH.EQ.0) GO TO 24
      IF(IPUNCH.EQ.2) GO TO 25
      IF(IPUNCH.EQ.1) GO TO 24
   25 CONTINUE
   24 END

************************************************************************

      SUBROUTINE SETNUM
*     CONSTANTS TO DETERMINE PARAMETER NUMBERS
      PARAMETER (NR=400)
      CHARACTER LINE*80
      COMMON /NONLIN/ NLESQ,NAUTO,IDCYCL,NCHNG,NPAR(50),IPAR(NR,50)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X

      DO J = 1, 10
         DO I = 1, NR
            IPAR(I,J) = 0.0
         END DO
      END DO
      
*    NROUGH: FIRST SURFACE-ROUGHNESS PARAMETER
*    NBCKGR: FIRST BACKGROUND PARAMETER
*     NUMBG: NUMBER OF BACKGROUND PARAMETERS
*     NPPPX: ABSOLUTE NUMBER FOR THE FIRST PRIMARY PROFILE PARAMETER 
*            (WHEN REFLECTIONS BROADENED ANISOTROPICALLY ARE SPECIFIED)
*    NRELAX: NUMBER OF RELAXED REFLECTIONS
*     NPRIM: NUMBER OF PRIMARY PROFILE PARAMETERS PER REFLECTION
*      NSFX: SCALE FACTOR
*     NPO1X: PREFERRED-ORIENTATION PARAMETER 1
*     NPO2X: PREFERRED-ORIENTATION PARAMETER 2
*       NAX: A
*       NBX: B
*       NCX: C
*     NALPX: ALPHA
*     NBETX: BETA
*     NGAMX: GAMMA
*      NOTX: OVERALL ISOTROPIC DISPLACEMENT PARAMETER
*      NG1X: OCCUPATION FACTOR OF THE 1ST SITE

      REWIND 4
      NRELAX = 0
      DO I = 1, 100000
         READ(4,'(A)',END=1) LINE
         DO J = 1, 80
            IF (LINE(J:J) .NE. ' ') EXIT
         END DO
         IDOT = INDEX(LINE,'.')
         IF (LINE(J:J+2) .EQ. 'PPP' .AND. LINE(J+4:J+4) .EQ. '_' .AND.
     &   IDOT .GE. J+6) THEN
            NRELAX = NRELAX + 1
         ELSE IF (INDEX('0123456789.-',LINE(J:J)) .GT. 0) THEN
            CYCLE
         ELSE IF (NRELAX .GT. 0) THEN
            EXIT
         END IF
      END DO
    1 REWIND 4

      NROUGH = 5
      NBCKGR = 9
      NUMBG = 12
      NPPPX = NBCKGR + NUMBG
*     NPRIM is given in SUBROUTINE FLSITE.
      NSFX = 0
      NPO1X = 15
      NPO2X = NPO1X + 1
      NAX = NPO2X + 1
      NBX = NPO2X + 2
      NCX = NPO2X + 3
      NALPX = NPO2X + 4
      NBETX = NPO2X + 5
      NGAMX = NPO2X + 6
      NOTX = NGAMX + 1
      NG1X = NOTX + 1
      END

************************************************************************

      BLOCK DATA VINIT
      INTEGER TRIC(2),MONO1(21),ORTHO1(146),TET1(178),TRIG1(43),
     &HEX(53),CUBIC1(108),MONO2(21),ORTHO2(18),TET2(48),TRIG2(11),
     &CUBIC2(30),SYMTAB(679),POINT1(256),POINT2(256),POINT3(36),
     &POINT(2,274),SYMCON(2,109),PTCON(2,48),SCLASS(13),CHK(7),
     &NCO(32),IEX(32),COITAB(5,14),SECSET,ISET(44)
      CHARACTER*8 TITLE(10),CTITLE(2,28),RTITLE(2,14),SYSNAM(2,7)
      COMMON
     &/TAB/ SYMTAB,POINT,SYMCON,PTCON,SCLASS,CHK
     &/PRT/ TITLE,CTITLE,RTITLE,SYSNAM
     &/COI/ NCO,IEX
     &/CTAB/ COITAB
     &/IS/ ISET,SECSET
      EQUIVALENCE (TRIC(1),  SYMTAB(1)),
     &            (MONO1(1), SYMTAB(3)),
     &            (ORTHO1(1),SYMTAB(24)),
     &            (TET1(1),  SYMTAB(170)),
     &            (TRIG1(1), SYMTAB(348)),
     &            (HEX(1),   SYMTAB(391)),
     &            (CUBIC1(1),SYMTAB(444)),
     &            (MONO2(1), SYMTAB(552)),
     &            (ORTHO2(1),SYMTAB(573)),
     &            (TET2(1),  SYMTAB(591)),
     &            (TRIG2(1), SYMTAB(639)),
     &            (CUBIC2(1),SYMTAB(650)),
     &            (POINT1(1),POINT(1,1)),
     &            (POINT2(1),POINT(1,129)),
     &            (POINT3(1),POINT(1,257))
      DATA TRIC/
     & 1, 1/
      DATA MONO1/
     & 1, 1,44, 3, 1, 1,20, 3, 3,20, 1, 1,44, 3, 1,20, 1,20,
     &44, 3,20/
      DATA ORTHO1/
     & 1, 1,44, 1,41,43, 1,41,43,44, 4,44, 4, 5, 6, 6, 1, 1,
     &16, 1,12,16, 1,15, 1,12,15, 1,13,16, 1,17, 1,11,15, 1,
     &13,15, 1,13,17, 4, 4,16, 4,12,16, 2, 2,11, 2,15, 2,11,
     &15, 5, 5,14,18, 6, 6,11,15, 6,15, 1, 1,13,17,21, 1,12,
     &16, 1,11,15,21, 1,19, 1,13,17,19, 1,17,19, 1,12,16,19,
     & 1,11,15, 1,12,16,21, 1,11,16, 1,13,17, 1,21, 1,11,16,
     &21, 1,11,16,19, 1,13,19, 4,16, 4,16,19, 4, 4,12,16, 4,
     &19, 4,12,16,19, 5, 5,14,18,22, 6, 6,11,15, 6,11,16,19,
     & 6,19/
      DATA TET1/
     & 1, 1,46, 1,44, 1,46, 6, 6,46, 1, 6, 1, 1,44, 1,21, 1,
     &21,44, 6, 6,34,46, 1, 1,41,43, 1,46, 1,41,43,46, 1,44,
     & 1,41,43,44, 1,46, 1,41,43,46, 6, 6,46, 1, 1,11,15, 1,
     &12,16, 1,13,17, 1,12,16,23, 1,13,17,23, 1,23, 1,11,15,
     &23, 6, 6,12,16, 6,25, 6,12,16,25, 1, 1,23, 1,41,43, 1,
     &23,41,43, 1, 1,12,16, 1,11,15, 1,13,17, 6, 6,12,16, 6,
     & 6,25, 1, 1,12,16,23, 1,11,15,21, 1,13,17,23,21, 1,11,
     &15, 1,13,17,23, 1,21, 1,12,16,23,21, 1,23, 1,12,16, 1,
     &11,15,23,21, 1,13,17,21, 1,11,15,23, 1,13,17, 1,23,21,
     & 1,12,16,21, 6, 6,12,16, 6,25,34, 6,12,16,25,34/
      DATA TRIG1/
     & 1, 1,45, 1,45, 9, 1, 9, 1, 1, 1,45, 1,45, 1,45, 1,45,
     & 9, 1, 1, 1,12,16, 1,23, 9, 9,23,28,31, 1, 1,23, 1, 1,
     &12,16, 9, 9,23,28,31/
      DATA HEX/
     & 1, 1,47, 1,47, 1,45, 1,45, 1,44, 1, 1, 1,44, 1, 1,47,
     & 1,47, 1,45, 1,45, 1,44, 1, 1,12,16,23, 1,12,16, 1,23,
     & 1, 1,12,16, 1, 1,23, 1, 1,12,16,23, 1,12,16, 1,23/
      DATA CUBIC1/
     & 1, 5, 6, 1,41,43,44, 6,41,43,44, 1, 1,13,17,21, 5, 5,
     &14,18,22, 6, 1,11,16,19, 6,34,35,36, 1, 1,41,43,44, 5,
     & 5,42,48,46, 6, 1,42,48,46, 1,42,48,46, 6,42,48,46, 1,
     & 5, 6, 1,24,29,32, 5,23,28,31, 6,25,30,33, 1, 1,13,17,
     &21,24,29,32, 1,24,29,32, 1,13,17,21, 5, 5,23,28,31, 5,
     &14,18,22, 5,23,28,31,14,18,22, 6, 6,34,35,36,25,30,33/
      DATA MONO2/
     & 1, 1,43, 4, 1, 1,16, 4, 4,16, 1, 1,43, 4, 1,16, 1,16,
     &43, 4,16/
      DATA ORTHO2/
     & 1,13,17,21, 1,11,15,21, 1,21, 4,12,16,19, 5,14,18,22/
      DATA TET2/
     & 1,21, 1,21,44, 6,34,46, 1,11,15,21, 1,13,17,23,21, 1,
     &21, 1,12,16,23,21, 1,11,15,23,21, 1,13,17,21, 1,23,21,
     & 1,12,16,21, 6,25,19, 6,12,16,25,19/
      DATA TRIG2/
     & 7, 7, 7, 7, 7,12,16, 7, 7,12,16/
      DATA CUBIC2/
     & 1,13,17,21, 5,14,18,22, 1,13,17,21,24,29,32, 1,13,17,
     &21, 5,14,18,22, 5,23,28,31,14,18,22/
      DATA POINT1/
     &  1,  1,  2,  1,  3,  1,  4,  2,  6,  1,  7,  1,  8,  2, 10,  1,
     & 11,  2, 13,  1, 14,  2, 16,  1, 17,  2, 19,  3, 22,  2, 24,  1,
     & 25,  2, 27,  3, 30,  4, 34,  2, 36,  1, 37,  1, 38,  1, 39,  1,
     & 40,  1, 41,  2, 43,  3, 46,  2, 48,  3, 51,  3, 54,  2, 56,  3,
     & 59,  3, 62,  3, 65,  1, 66,  2, 68,  3, 71,  1, 72,  2, 74,  2,
     & 76,  3, 79,  1, 80,  3, 83,  1, 84,  3, 87,  2, 89,  1, 90,  4,
     & 94,  3, 97,  4,101,  2,103,  4,107,  3,110,  4,114,  3,117,  4,
     &121,  3,124,  3,127,  2,129,  4,133,  4,137,  3,140,  2,142,  3,
     &145,  1,146,  3,149,  2,151,  4,155,  1,156,  4,160,  1,161,  3,
     &164,  4,168,  2,170,  1,171,  2,173,  2,175,  2,177,  1,178,  2,
     &180,  1,181,  1,182,  1,183,  2,185,  2,187,  3,190,  1,191,  3,
     &194,  1,195,  3,198,  2,200,  4,204,  2,206,  4,210,  2,212,  4,
     &216,  1,217,  2,219,  1,220,  3,223,  3,226,  3,229,  4,233,  4,
     &237,  2,239,  4,243,  1,244,  3,247,  2,249,  4,253,  1,254,  2,
     &256,  3,259,  4,263,  1,264,  3,267,  3,270,  3,273,  1,274,  3,
     &277,  1,278,  2,280,  1,281,  4,285,  4,289,  5,294,  3,297,  4/
      DATA POINT2/
     &301,  2,303,  5,308,  2,310,  3,313,  5,318,  4,322,  4,326,  3,
     &329,  3,332,  4,336,  1,337,  3,340,  3,343,  5,348,  1,349,  2,
     &351,  2,353,  1,354,  1,355,  1,356,  1,357,  1,358,  2,360,  2,
     &362,  2,364,  2,366,  1,367,  1,368,  1,369,  3,372,  2,374,  1,
     &375,  4,379,  1,380,  2,382,  1,383,  3,386,  1,387,  4,391,  1,
     &392,  2,394,  2,396,  2,398,  2,400,  2,402,  1,403,  1,404,  2,
     &406,  1,407,  2,409,  2,411,  2,413,  2,415,  2,417,  1,418,  4,
     &422,  3,425,  2,427,  1,428,  3,431,  1,432,  2,434,  1,435,  4,
     &439,  3,442,  2,444,  1,445,  1,446,  1,447,  4,451,  4,455,  1,
     &456,  4,460,  1,461,  4,465,  1,466,  4,470,  4,474,  1,475,  4,
     &479,  1,480,  4,484,  1,485,  4,489,  4,493,  4,497,  1,498,  1,
     &499,  1,500,  4,504,  4,508,  4,512,  1,513,  7,520,  4,524,  4,
     &528,  1,529,  4,533,  4,537,  7,544,  1,545,  7,552,  1,553,  2,
     &555,  1,556,  1,557,  2,559,  1,560,  2,562,  1,563,  2,565,  1,
     &566,  2,568,  3,571,  2,573,  4,577,  4,581,  2,583,  4,587,  4,
     &591,  2,593,  3,596,  3,599,  4,603,  5,608,  2,610,  5,615,  5/
      DATA POINT3/
     &620,  4,624,  3,627,  4,631,  3,634,  5,639,  1,640,  1,641,  1,
     &642,  1,643,  3,646,  1,647,  3,650,  4,654,  4,658,  7,665,  4,
     &669,  4,673,  7/
      DATA SYMCON/
     & 0, 7, 4, 7, 5, 7, 6, 7, 7, 7,13, 7,14, 7,15, 7,21, 7,
     & 0, 0, 2, 3, 2, 8, 2, 2, 3, 3, 2, 8, 3, 1, 4, 3, 2, 2,
     & 3, 1, 8, 3,17, 2,19, 1, 2, 8, 1, 5, 1,11, 1, 4, 3, 5,
     & 1,11, 3, 1, 5, 5, 1, 4, 3, 1, 9, 5,16, 4,19, 1, 1,11,
     & 1, 6, 1, 9, 1, 4, 2, 6, 1, 9, 2, 2, 6, 6, 1, 4, 2, 2,
     &10, 6,16, 4,17, 2, 1, 9, 3,10, 3, 1, 1,14, 3,10, 3, 1,
     & 1,14,12,10,16,14,19, 1, 1, 9, 3,10, 3, 1, 1,14, 3, 5,
     & 1,11, 3, 1, 2,13, 2, 2, 1,14, 2,13, 2, 2, 1,14,23,13,
     &17, 2,16,14, 1,11, 1,12, 1, 4, 1,14, 1,12, 1, 4, 1,14,
     &25,12,16, 4,16,14, 2, 8,26, 6, 1, 9, 1, 4, 2, 2,27, 3,
     & 2, 8, 2, 2, 3, 1,28, 5, 1,11, 1, 4, 3, 1, 0, 0, 0, 0,
     & 0, 0, 0, 0, 1, 4,16, 4, 2, 2, 3, 1,18, 1,19, 1,20, 1,
     &17, 2/
      DATA PTCON/
     &  1,  1,  2,  1,  3,  1,  4,  1,  5,  1,  6,  1,  7,  1,  8,  1,
     &  9,  1, 10,  1, 11,  3, 14,  3, 17,  3, 20,  4, 24,  3, 27,  3,
     & 30,  3, 33,  4, 37,  3, 40,  3, 43,  3, 46,  4, 50,  3, 53,  3,
     & 56,  4, 60,  3, 63,  3, 66,  3, 69,  3, 72,  4, 76,  3, 79,  3,
     & 82,  4, 86,  4, 90,  4, 94,  4, 98,  1, 99,  1,100,  1,101,  1,
     &102,  1,103,  1,104,  1,105,  1,106,  1,107,  1,108,  1,109,  1/
      DATA SCLASS/1,3,16,75,143,168,195,231,244,249,262,269,275/
      DATA CHK/146,148,155,160,161,166,167/
      DATA CTITLE/
     &'h=2n;   ','        ','k=2n;   ','        ','l=2n;   ','        ',
     &'k+l=2n; ','        ','h+l=2n; ','        ','h+k=2n; ','        ',
     &'h,k,l od','d/even; ','k+l=4n; ','        ','h+l=4n; ','        ',
     &'h+k=4n; ','        ','2h+l=2n;','        ','2h+l=4n;','        ',
     &'h+k+l=2n',';       ','-h+k+l=3','n;      ','h-k+l=3n',';       ',
     &'h=4n;   ','        ','k=4n;   ','        ','l=3n;   ','        ',
     &'l=4n;   ','        ','l=6n;   ','        ','|h|>=|k|','>=|l|;  ',
     &'2h+k=2n;','        ','2h+k=4n;','        ','h+2k=2n;','        ',
     &'h+2k=4n;','        ','h=2n,k=2','n;      ','k=2n,l=2','n;      ',
     &'h=2n,l=2','n;      '/
      DATA RTITLE/
     &'00l;    ','        ','0k0;    ','        ','0kl;    ','        ',
     &'h00;    ','        ','h0l;    ','        ','hk0;    ','        ',
     &'All refl','ections;','0kk;    ','        ','hh0;    ','        ',
     &'hhl;    ','        ','h0h;    ','        ','hkk;    ','        ',
     &'hkh;    ','        ','hhh;    ','        '/
      DATA SYSNAM/
     &'triclini','c       ','monoclin','ic      ','orthorho','mbic    ',
     &'tetragon','al      ','rhombohe','dral    ','hexagona','l       ',
     &'cubic   ','        '/
      DATA NCO/74,88,142,145,146,147,148,149,150,151,152,153,154,
     &155,156,157,158,159,161,163,165,167,176,194,206,248,251,261,
     &263,268,270,274/
      DATA IEX/0,1,0,2,4,2,4,5,3,5,3,5,3,0,3,5,3,5,0,5,3,0,1,0,4,
     &0,1,0,2,3,4,0/
      DATA COITAB/0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,1,1,
*     BUG POINTED OUT BY M. ONODA
C    &1,0,1,1,0,1,0,1,3,1,1,2,0,0,1,0,0,0,0,0,0,0,0,2,0,0,2,0,0,
     &1,0,1,1,0,1,0,1,3,1,1,2,0,1,1,0,0,0,0,0,0,0,0,2,0,0,2,0,1,
     &1,0,0,1,3,1,0,2,1,3,1,0,2,0,2,0,0,2/
      DATA ISET/3,4,5,6,7,8,9,10,11,12,13,14,15,48,50,59,68,70,85,86,88,
     &125,126,129,130,133,134,137,138,141,142,146,148,155,160,161,166,
     &167,201,203,222,224,227,228/
      END

************************************************************************

      SUBROUTINE POSREF(COL,IC,NCON)
C     DETERMINE THE COEFFICIENTS IN EQUATIONS
C     (C1)H+(C2)K+(C3)L=(C4)N+(C5); .....
C     UP TO THREE EQUATIONS MAY BE INPUT.  THEY ARE SEPARATED BY ','
C     OR ';'
C     IC(I,J): I SCANS COEFFICIENTS C1 - C5
C              J SCANS EQUATIONS (RELATED BY 'AND')
      CHARACTER COL(80)*1,CH(14)*1
      INTEGER IC(5,3),NCON,NUM,SIGN
      DATA CH/'h','k','l','n','0','1','2','3','4','5','6','7','8','9'/
C
      DO 20 J=1,3
         DO 10 I=1,5
            IC(I,J)=0
   10    CONTINUE
   20 CONTINUE
      J=1
C     LOUT=0:  AN INTEGER HAS BEEN ENCOUNTERED
C     LOUT=1:  AN INTEGER HAS NOT YET BEEN ENCOUNTERED OR HAS BEEN
C              STORED IN ARRAY IC
      LOUT=1
      SIGN=1
      NUM=1
      DO 50 ICOL=1,80
         IF (COL(ICOL) .EQ. '-') THEN
            SIGN=-1
         ELSE IF (COL(ICOL) .NE. ' ') THEN
            DO 25 I=1,14
               IF (COL(ICOL) .EQ. CH(I)) THEN
                  LP=I
                  GO TO 1
               END IF
   25       CONTINUE
            LP=0
    1       IF (LP .GT. 0 .AND. LP .LE. 4) THEN
               IC(LP,J)=SIGN*NUM
               LOUT=1
               SIGN=1
               NUM=1
            ELSE IF (LP .GE. 5) THEN
               NUM=LP-5
               LOUT=0
            ELSE IF (COL(ICOL) .EQ. ';' .OR. COL(ICOL) .EQ. ',') THEN
               IF (LOUT .EQ. 0) IC(5,J)=SIGN*NUM
               LOUT=1
               SIGN=1
               NUM=1
               J=J+1
            END IF
         END IF
   50 CONTINUE
      IF (LOUT .EQ. 0) IC(5,J)=SIGN*NUM
      NCON=J
      END

************************************************************************

      INTEGER FUNCTION IPOS(H,K,L,NHKL,ICH,NCON,NLINE)
*     CHECK WHETHER THIS REFLECTION IS POSSIBLE
      INTEGER H,K,L,ICH(5,3,7),NCON(7)

      IPOS=1
      IF (NHKL .EQ. 2 .AND. H .NE. 0) THEN
         RETURN
      ELSE IF (NHKL .EQ. 3 .AND. K .NE. 0) THEN
         RETURN
      ELSE IF (NHKL .EQ. 4 .AND. L .NE. 0) THEN
         RETURN
      ELSE IF (NHKL .EQ. 5 .AND. H .NE. K) THEN
         RETURN
*     ELSE IF (NHKL .EQ. 6 .AND. H .NE. -H) THEN
*        RETURN
      END IF

      DO 50 J=1,NLINE
         DO 40 I=1,NCON(J)
            IF (ICH(4,I,J) .NE. 0 .AND. MOD(ICH(1,I,J)*H + ICH(2,I,J)
     &      *K + ICH(3,I,J)*L - ICH(5,I,J), ICH(4,I,J)) .EQ. 0) THEN
C              DO NOTHING
            ELSE IF (ICH(4,I,J) .EQ. 0 .AND. ICH(1,I,J)*H + ICH(2,I,J)
     &      *K + ICH(3,I,J)*L .EQ. ICH(5,I,J)) THEN
C              DO NOTHING
            ELSE
               GO TO 50
            END IF
   40    CONTINUE
         RETURN
   50 CONTINUE
      IPOS=0
      END

************************************************************************

      SUBROUTINE HRHT
C     HS AND H*TS
C     'JIKKEN KAGAKU KOHZA,' VOL. 6, MARUZEN, TOKYO (1977), P. 60.
      INTEGER H,R,RX,RY,RZ
      PARAMETER (NB=7000,NS=48,NAP=150,NPH=8)
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /G/ I
      COMMON /RT/ T(3,NS,NPH),HT(NS,NB),NSYM(NPH),NSITE(NPH),
     &  IDSYM(NS,NAP,NPH),R(3,3,NS,NPH),RX(NS,NB),RY(NS,NB),RZ(NS,NB)
     
      IP=NOPH(I)
      DO J = 1, NSYM(IP)
         RX(J,I) = H(I)*R(1,1,J,IP) + K(I)*R(2,1,J,IP) + 
     &             L(I)*R(3,1,J,IP)
         RY(J,I) = H(I)*R(1,2,J,IP) + K(I)*R(2,2,J,IP) +
     &             L(I)*R(3,2,J,IP)
         RZ(J,I) = H(I)*R(1,3,J,IP) + K(I)*R(2,3,J,IP) +
     &             L(I)*R(3,3,J,IP)
         HT(J,I) = FLOAT(H(I))*T(1,J,IP) + FLOAT(K(I))*T(2,J,IP) +
     &             FLOAT(L(I))*T(3,J,IP)
      END DO
      END

************************************************************************

      SUBROUTINE UPDATE(A,IDPOS,ID1ST)
C     UPDATE VARIOUS ARRAYS SUCH AS FWHM AND PEAK POSITIONS
      PARAMETER (NB=7000,NS=48,NAP=150,NT=999,NCS=400,NA=15,NPH=8)
      INTEGER H,HP,R,RX,RY,RZ,HANIS
      DIMENSION A(*)
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /ANISBR/ HANIS(NPH),KANIS(NPH),LANIS(NPH),
     &  CANIS(NPH),PCOS(NB)
      COMMON /B/ U(NB),NREFL
      COMMON /BUFFER/ TEMP(NB)
      COMMON /DERPR/ POP2(NB),COMPO(NB),CS(NB),CCC(6,NB)
      COMMON /G/ I
      COMMON /I/ RDEG(2,NB)
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /LTM/ LISO(NAP,NPH),LMAG(NPH)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /ORDER/ IPOINT(NB)
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /PROFI/ PCOT(NB),SHFT(NB),NTSIM(NB),FNORM(NB),TNTSIM(NB),
     &  FWHM(NB),ETA(NB),SQSG(NB),CTH(NB),SIGP(NB),DFWDG(NB),DEDFFG(NB),
     &  DFWDL(NB),DEDFFL(NB),SIGPAP(NB)
      COMMON /RT/ T(3,NS,NPH),HT(NS,NB),NSYM(NPH),NSITE(NPH),
     &  IDSYM(NS,NAP,NPH),R(3,3,NS,NPH),RX(NS,NB),RY(NS,NB),RZ(NS,NB)
      COMMON /S/ F(NA,NB),NATOM
*     F0: Nondispersive atomic scattering factor
      COMMON /S2/ F0(NA,NB)
      COMMON /SFRI/ SFR(NB),SFI(NB),OTF(NB)
      COMMON /ABPI/ AP(NB),BP(NB)
      COMMON /T/ FF(NB),COEF(NPH)
      COMMON /TENS/ G(6,NPH),GG(6,NPH)
      COMMON /UCC/ IDPEAK(NPH),IDSF(NPH),IDPRF(NPH),IDPO(NPH)
      COMMON /V/ NPROR(NPH),HP(NPH),KP(NPH),LP(NPH),LSUM(NPH)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
      COMMON /Z/ PC,CHGPC
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      COMMON /ZSPFN/ PEAKSH(NB),SQH2(NB),SQH2L(NB),SQH2G(NB),ASYM(NB),
     &ETAL(NB),ETAH(NB),FNORM1(2,NB),COMLG(2,NB),COML(2,NB),COMG(2,NB),
     &COEFL(2,NB),COEFG(2,NB),DECAY(2,NB),COMPEAR(2,NB),CDX(2,NB),
     &CDA1(2,NB),CDA2(2,NB),CDW(2,NB),CDRL1(2,NB),CDRL2(2,NB),
     &CDRH1(2,NB),DMPART(2,NB)
      COMMON /UNCONS/ NUCREF,IPHASE(90),IHKL(3,90),NOPPP(90),LPPP(NB),
     &  RCUT(2,90),NCUT
*Rev 1.1 2003.05.21
      SAVE PCORG,ICHGPC

      P(J)=A(KPHB(NOPH(I))+J)

C     ID1ST=0: THIS SUBROUTINE IS CALLED FOR THE FIRST TIME
*Rev 1.1 2003.05.25 {
      IF (NMODE .EQ. 4 .AND. NSFF .NE. 1 .AND. ID1ST .EQ. 0) THEN
         PCORG = PC
*        CHANGE THE PROFILE CUTOFF IN THE FIRST CYCLE
         PC = CHGPC*PC
         ICHGPC = 1
         IF (ABS(CHGPC-1.0) .GT. 1.0E-4) WRITE(6,'(//1H ,10X,A,F8.4)')
     &   'PC has been changed to',PC
      END IF
* }

C     IDPEAK=1: CHANGE BRAGG POSITIONS
C     IDPRF=1: CHANGE PROFILE PARAMETERS
C     IDSF=1: CHANGE STRUCTURE-FACTOR PARAMETER(S)
C     IDPOS=0: BRAGG POSITIONS CANNOT BE CALCULATED
C     IDPOS=1: BRAGG POSITIONS AND PEAK SHAPE FUNCTION CAN BE CALCULATED
C     IDPOS=-1: A PART OF THE PROFILE FUNCTION IS UNREASONABLE
      IDPOS=1
      CALL SETPAR(A)
      CALL TENSOR(A)
      IDBR=0
      DO IP = 1, NPHASE
         IF (ID1ST .EQ. 0) THEN
            IDPEAK(IP)=0
            IDPRF(IP)=1
            SELECT CASE (NMODE)
               CASE (0:3,5)
                  IDSF(IP) = 1
               CASE (4)
*                 CONVENTIONAL LE BAIL REFINEMENT
                  IF (IP .EQ. 1) THEN
*                    NO STRUCTURE FACTORS ARE REQUIED FOR THE FIRST PHASE 
                     IDSF(1) = 0
                  ELSE
*                    STRUCTURE FACTORS ARE CALCULATED EXCEPT THE FIRST PHASE 
                     IDSF(IP) = 1
                  END IF
               CASE (6)
*                 INDIVIDUAL PROFILE FITTING
*                 DO NOTHING BECAUSE INTEGRATED INTENSITIES ARE
*                 REFINED
            END SELECT
            IDPO(IP)=1
         END IF
         IDBR=IDBR+IDPEAK(IP)
      
*        CPROR IS USED IN SUBROUTINE PROR
         IF (NPROR(IP) .GE. 1)
     &   CPROR(IP) = FLOAT(HP(IP)*HP(IP))*GG(1,IP) + 
     &   FLOAT(KP(IP)*KP(IP))*GG(2,IP) + FLOAT(LP(IP)*LP(IP))*GG(3,IP)
     &   + 2.0*(FLOAT(KP(IP)*LP(IP))*GG(4,IP) + 
     &   FLOAT(LP(IP)*HP(IP))*GG(5,IP) + FLOAT(HP(IP)*KP(IP))*GG(6,IP))

*        CANIS IS USED IN SUBROUTINE LBROAD
         CANIS(IP) = FLOAT(HANIS(IP)*HANIS(IP))*GG(1,IP) + 
     &   FLOAT(KANIS(IP)*KANIS(IP))*GG(2,IP) +
     &   FLOAT(LANIS(IP)*LANIS(IP))*GG(3,IP) + 
     &   2.0*(FLOAT(KANIS(IP)*LANIS(IP))*GG(4,IP) +
     &   FLOAT(LANIS(IP)*HANIS(IP))*GG(5,IP) + 
     &   FLOAT(HANIS(IP)*KANIS(IP))*GG(6,IP))
      END DO

*Rev 1.1 2003.05.25 {
      IF (NMODE .EQ. 4 .AND. NSFF .NE. 1 .AND. ID1ST .NE. 0 .AND.
     &ICHGPC .EQ. 1 .AND. IDPRF(1) .EQ. 1) THEN
*        PC is restored because IDPRF for the 1st phase is now 1.
         PC = PCORG
         ICHGPC = 0
         IF (ABS(CHGPC-1.0) .GT. 1.0E-4) WRITE(6,'(//1H ,10X,A,F8.4)')
     &   'PC has been restored to',PC
      END IF
* }

C     CALCULATE VARIOUS ARRAYS USING NEW INDICES AND PARAMETERS
      DO I = 1, NREFL
         IF (IDPEAK(NOPH(I)) .EQ. 1) THEN
            CALL BRANG
            IF (D(I) .EQ. 0.0) THEN
               IDPOS=0
               RETURN
            ENDIF
         END IF

*        A PEAK POSITION FOR A RELAXED REFLECTION IS RESTORED
         IF (NMODE .EQ. 6 .AND. LPPP(I) .GT. 0) CALL BRANGPF(A)
         IF ((NMODE .NE. 1 .OR. NPAT .GT. 0) .AND. IDPRF(NOPH(I)) 
     &   .EQ. 1) THEN
C           RDEG(1,I)-RDEG(2,I): REGION IN WHICH REFLECTION I
C           CONTRIBUTES TO THE INTENSITY
            
            SELECT CASE (NPRFN)
               CASE (0)
                  CALL PREDE0(A,IDPOS)
                  IF (IDPOS .EQ. -1) RETURN
               CASE (1:3)
                  SELECT CASE (NPRFN)
                     CASE (1,2)
                        CALL PREDE12(A,IDPOS)
                     CASE (3)
                        CALL PREDE3(A,IDPOS)
                  END SELECT
                  IF (IDPOS .EQ. -1) RETURN

                  IF (PC .GT. 1.0) THEN
                     IF (NCUT .EQ. 1 .AND. LPPP(I) .NE. 0) THEN
                        RDEG(1,I) = RCUT(1,LPPP(I))
                        RDEG(2,I) = RCUT(2,LPPP(I))
                     ELSE 
*                       FWHM = 2.0/(1.0+A)*W*A (LOWER ANGLE)
                        RDEG(1,I) = PEAK(I) - PEAKSH(I) - 
     &                              PC*2.0/(1.0+ASYM(I))*SQH2(I)*ASYM(I)
*                       FWHM = 2.0/(1.0+A)*W (HIGHER ANGLE)
                        RDEG(2,I) = PEAK(I) - PEAKSH(I) + 
     &                              PC*2.0/(1.0+ASYM(I))*SQH2(I)
                     END IF
                  END IF

*                 IF NMODE = 3, 6 AND REFLECTION I IS RELAXED, |Fc| MAY BE 
*                 REFINED.  THEREFORE, FF(I) MUST BE UPDATED HERE THOUGH IT
*                 IS NOT A PROFILE PARAMETER.
                  IF ((NMODE .EQ. 3 .OR. NMODE .EQ. 6) .AND. 
     &            LPPP(I) .GT. 0) THEN
                     SELECT CASE (NPRFN)
                        CASE (1,3)
                           FF(I) = A(NOPPP(LPPP(I)) + 4)**2
                        CASE (2)
                           FF(I) = A(NOPPP(LPPP(I)) + 5)**2
                     END SELECT
                  END IF
            END SELECT
         END IF

         IF (L12(I) .NE. 0) CYCLE
         IF (IDPO(NOPH(I)) .EQ. 1) CALL PROR(A)
         IF (IDSF(NOPH(I)) .EQ. 1) CALL SF(A,RX(1,I),RY(1,I),RZ(1,I),
     &   HT(1,I),NSITE(NOPH(I)),NSYM(NOPH(I)),IDSYM(1,1,NOPH(I)))
      END DO

      DO I = 1, NREFL
         II=L12(I)
         IF (II .NE. 0) THEN
C           ARRAYS FOR ALPHA-2 PEAKS
            PROR1(I)=PROR1(II)
            DPROR(I)=DPROR(II)
            IF (NPROR(NOPH(I)) .LE. 2) POP2(I)=POP2(II)
            IF (IDSF(NOPH(I)) .EQ. 1) THEN
               FF(I)=FF(II)
               IF (NMODE .LE. 1) THEN
                  SFR(I)=SFR(II)
                  SFI(I)=SFI(II)
               END IF
               OTF(I)=OTF(II)
            END IF
         END IF
         IF ((NMODE .EQ. 4 .OR. NMODE .EQ. 5) .AND. NOPH(I) .EQ. 1 .AND.
     &   ID1ST .NE. 0) CYCLE
*        YPEAK WILL BE DETERMINED IN SUBROUTINE INITINT WHEN NSFF = 0
         IF (NMODE .EQ. 4 .AND. NOPH(I) .EQ. 1 .AND. NSFF .EQ. 2) THEN
*           |F|**2 or m*|F|**2 is set at 100.0
            SELECT CASE (INCMULT)
               CASE (0)
*                 NO MULTIPLICITY IS INCLUDED IN 100.0
                  YPEAK(I) = P(0)*U(I)*FLP(I)*100.0
               CASE (1)
*                 A MULTIPLICITY IS INCLUDED IN 100.0
                  YPEAK(I) = P(0)*FLP(I)*100.0
                  IF (II .NE. 0) YPEAK(I) = YPEAK(I)*R12
            END SELECT
         ELSE
            YPEAK(I) = P(0)*U(I)*FLP(I)*PROR1(I)*FF(I)
         END IF
      END DO
      IF (PC .LT. 1.0 .AND. (IDPRF(1) .EQ. 1 .OR. ID1ST .EQ. 0)) 
     &  CALL CUTOFF(A)

      DO I = 1, NREFL
         TEMP(I)=PEAK(I)
      END DO
      CALL EXCTES(TEMP,IPOINT,NREFL)
      END

************************************************************************

      SUBROUTINE CUTOFF(A)
*     DETERMINE THE PROFILE CUT-OFF FOR EACH REFLECTION (SPECIFIC
*     FOR THE PSEUDO-VOIGT FUNCTION OF THOMPSON, COX, AND HASTINGS)
      EXTERNAL PRFL0
      PARAMETER (NB=7000)
      INTEGER IOPT(3)
      REAL EOPT(1)
      REAL A(*),TTHPEAK(100),TTHMAX(100),TTHMIN(100)
      COMMON /B/ U(NB),NREFL
      COMMON /G/ I
      COMMON /I/ RDEG(2,NB)
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /UNCONS/ NUCREF,IPHASE(90),IHKL(3,90),NOPPP(90),LPPP(NB),
     &  RCUT(2,90),NCUT
      COMMON /PROFI/ PCOT(NB),SHFT(NB),NTSIM(NB),FNORM(NB),TNTSIM(NB),
     &  FWHM(NB),ETA(NB),SQSG(NB),CTH(NB),SIGP(NB),DFWDG(NB),DEDFFG(NB),
     &  DFWDL(NB),DEDFFL(NB),SIGPAP(NB)
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      COMMON /ZSPFN/ PEAKSH(NB),SQH2(NB),SQH2L(NB),SQH2G(NB),ASYM(NB),
     &ETAL(NB),ETAH(NB),FNORM1(2,NB),COMLG(2,NB),COML(2,NB),COMG(2,NB),
     &COEFL(2,NB),COEFG(2,NB),DECAY(2,NB),COMPEAR(2,NB),CDX(2,NB),
     &CDA1(2,NB),CDA2(2,NB),CDW(2,NB),CDRL1(2,NB),CDRL2(2,NB),
     &CDRH1(2,NB),DMPART(2,NB)
      DATA RD/0.01745329/
*     TWO ARRAYS USED IN SILUP
      DATA NDEG/2/
      DATA LUP/0/
      DATA IOPT/0,1,0/
      DATA TOL/8.726646E-06/ ! 0.0005 degrees
      DATA PRFRNG/10.0/

      I = 1
      DO IPROF = 1, 100
*        TTHPEAK: PEAK POSITION FOR THE I'TH REFLECTION
*        YMAX: INTENSITY AT THE PEAK POSITION
         IF (NPRFN .EQ. 0) THEN
            MODE = 0
*           THE PEAK POSITION IS SANDWICHED BETWEEN X AND XORF
*           +/- 0.8 DEGREES
            X = PEAK(I) - SHFT(I) - 0.8*RD
            XORF = PEAK(I) - SHFT(I) + 0.8*RD
            DO WHILE (.TRUE.)
               CALL SFMIN(X, XORF, MODE, TOL)
               IF (MODE .NE. 1) EXIT
               XORF = -PRFL(X,A)
            END DO
            IF (MODE .EQ. 4) CALL JOBEND('Error termination in SFMIN')
            TTHPEAK(IPROF) = X
            YMAX = -XORF
         ELSE
            TTHPEAK(IPROF) = PEAK(I) - PEAKSH(I)
            YMAX = PRFL(TTHPEAK(IPROF),A)
         END IF

         XMIN = TTHPEAK(IPROF) - PRFRNG*RD
         IF (XMIN .LT. 0.0) XMIN = 0.001*RD
         TTHMIN(IPROF) = ZBRENT(PRFL0, XMIN, TTHPEAK(IPROF), 0.001*RD,
     &   A, YMAX)
         IF (TTHMIN(IPROF) .EQ. 0.0) TTHMIN(IPROF) = XMIN
         XMAX = TTHPEAK(IPROF) + PRFRNG*RD
         IF (XMAX .GT. 178.0*RD) XMAX = 178.0*RD
         TTHMAX(IPROF) = ZBRENT(PRFL0, TTHPEAK(IPROF), XMAX, 0.001*RD,
     &   A, YMAX)
         IF (TTHMAX(IPROF) .EQ. 0.0) TTHMAX(IPROF) = XMAX
         IF (NREFL .LE. 70) THEN
            RDEG(1,I) = TTHMIN(IPROF)
            RDEG(2,I) = TTHMAX(IPROF)
            I = I + 1
            IF (I .GT. NREFL) RETURN
         ELSE
            DO J = I + 1, NREFL
               IF (PEAK(J) .GT. PEAK(I) + 2.0*RD) GO TO 4
            END DO
*           IN THE CASE OF THE NREFL'TH REFLECTION, RTIME IS ALWAYS
*           CALCULATED BECAUSE INTERPOLATION IS NOT ACCURATE FOR IT.
            IF (I .EQ. NREFL) EXIT
            J = NREFL
    4       I = J
         END IF
      END DO
      NPRFL = IPROF - 1

*     DETERMINE THE PROFILE CUT-OFF BY INTERPOLATION
      DO I = 1, NREFL
         IF (NCUT .EQ. 1 .AND. LPPP(I) .NE. 0) THEN
            RDEG(1,I) = RCUT(1,LPPP(I))
            RDEG(2,I) = RCUT(2,LPPP(I))
            CYCLE
         END IF
         IF (NPRFN .EQ. 0) THEN
            PEAKCOR = PEAK(I) - SHFT(I)
         ELSE
            PEAKCOR = PEAK(I) - PEAKSH(I)
         END IF
         CALL SILUP
     &   (PEAKCOR,RDEG(2,I),NPRFL,TTHPEAK,TTHMAX,NDEG,LUP,IOPT,EOPT)
         CALL SILUP
     &   (PEAKCOR,RDEG(1,I),NPRFL,TTHPEAK,TTHMIN,NDEG,LUP,IOPT,EOPT)
      END DO
      END

************************************************************************

      FUNCTION ZBRENT(FUNC,X1,X2,TOL,AA,YMAX)
*     
      INTEGER ITMAX
      REAL ZBRENT,TOL,X1,X2,FUNC,EPS,AA(*),YMAX
      EXTERNAL FUNC
      PARAMETER (ITMAX=100,EPS=3.E-8)
      INTEGER ITER
      REAL A,B,C,D,E,FA,FB,FC,P,Q,R,S,TOL1,XM
      
      ZBRENT = 0.0
      A=X1
      B=X2
      FA=FUNC(A,AA,YMAX)
      FB=FUNC(B,AA,YMAX)
*     A root is not sandwiched between the two extremes
      IF ((FA .GT. 0. .AND. FB .GT. 0.) .OR. 
     &(FA .LT. 0. .AND. FB .LT. 0.)) RETURN
      C=B
      FC=FB
      DO ITER=1,ITMAX
        IF((FB.GT.0..AND.FC.GT.0.).OR.(FB.LT.0..AND.FC.LT.0.))THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          ZBRENT=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B,AA,YMAX)
      END DO
      CALL JOBEND('ZBRENT exceeding maximum iterations')
      ZBRENT=B
      END

************************************************************************

      FUNCTION PRFL0(TTH,A,YMAX)
      REAL A(*),YMAX
      COMMON /Z/ PC,CHGPC

      PRFL0 = PRFL(TTH,A) - YMAX*PC
      END

************************************************************************

      SUBROUTINE PREDE0(A,IDPOS)
*     UPDATE ARRAYS FOR THE CALCULATION OF THE PSEUDO-VOIGT FUNCTION
*     OF THOMPSON, COX, AND HASTING MADE ASYMMETRIC BY HOWARD'S
*     PROCEDURE

      INTEGER H
      REAL A(*),COFT(6),COFN(3)
      PARAMETER (NB=7000,NAP=150,NPH=8)
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /ANISBR/ HANIS(NPH),KANIS(NPH),LANIS(NPH),
     &  CANIS(NPH),PCOS(NB)
      COMMON /G/ I
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /L/ DEG1,DEG2
      COMMON /MC/ CTHM,NTRAN,PCOR
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /PROFI/ PCOT(NB),SHFT(NB),NTSIM(NB),FNORM(NB),TNTSIM(NB),
     &  FWHM(NB),ETA(NB),SQSG(NB),CTH(NB),SIGP(NB),DFWDG(NB),DEDFFG(NB),
     &  DFWDL(NB),DEDFFL(NB),SIGPAP(NB)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM

      DATA COFT /1.0, 2.69269, 2.42843, 4.47163, 0.07842, 1.0/
      DATA COFN /1.36603, -0.47719, 0.11116/
      P(J)=A(KPHB(NOPH(I))+J)

*     STORE ARRAYS TO CALCULATE PROFILE FUNCTIONS
      CALL LBROAD
      CTH(I) = 0.0
      IF (ABS(TH(I) * 57.29578 - 45.0) .GT. 1.0E-3) 
     &  CTH(I) = 1.0/TAN(PEAK(I))
      CSTH = COS(TH(I))
      STTH = SIN(PEAK(I))
      TTH = TAN(TH(I))

*     PEAK-SHIFT TERMS
*     Suggestion by Dr. Bordet on 1995.12.12
*Rev 1.07 2002.08.22 Izumi {
      IF (NTRAN .EQ. 2) THEN
*        TRANSMISSION GEOMETRY
         SHFT(I) = A(1) + A(2)*SIN(TH(I)) + A(3)*STTH
      ELSE
         SHFT(I) = A(1) + A(2)*CSTH + A(3)*STTH
      END IF
C     IF (NTRAN .EQ. 0) THEN
C        SHFT(I) = A(1) + A(2)*CSTH + A(3)*STTH
C     ELSE
C        SHFT(I) = A(1) + A(2)*SIN(TH(I)) + A(3)*STTH
C     END IF
* }

*     SIG: GAUSSIAN VARIANCE (NOT Hk**2 BUT sigma**2)

*     ORIGINAL EXPRESSION OF SIG
*     THERE MAY BE A SERIOUS PROBLEM WHEN 
*     (P(1)*TTH + P(2))*TTH + P(3) + P(4)/CSTH**2 < 0
C     SIG = ABS((P(1)*TTH + P(2))*TTH + P(3) + P(4)/CSTH**2)

*     MY EXPRESSION OF SIG
      SIG = (P(1)*TTH + P(2))*TTH + P(3) + P(4)/CSTH**2
      IF (SIG .LE. 0.0) THEN
         IDPOS = -1
         RETURN
      END IF

*     GAM: LORENTZIAN FWHM.  PCOS(I) HAS BEEN CALCULATED IN LBROAD.

*     ORIGINAL EXPRESSION OF GAM
C     GAM = MAX(0.0, (P(5) + P(6)*PCOS(I))/CSTH + 
C    &      (P(7) + P(8)*PCOS(I))*TTH)

*     MY EXPRESSION FOR H.k(Lorentz)
      GAM = (P(5) + P(6)*PCOS(I))/CSTH + (P(7) + P(8)*PCOS(I))*TTH
      IF (GAM .LT. 0.0) THEN
         IDPOS = -1
         RETURN
      END IF

      FWHMG = 2.354820 * SQRT(SIG) ! H.k(Gauss)
*     FWHMGL = FWHM(GAUSS) + FWHM(LORENTZ)
      FWHMGL =  FWHMG + GAM ! H.k(Gauss) + H.k(Lorentz)

      IF (NASYM .EQ. 1) THEN
**       As * cot(2-theta)
         PCOT(I) = P(9) * CTH(I)
*        CF. THE CAPTION OF FIG. 1 IN HOWARD'S PAPER
         KT = MAX(1, IFIX(2.5 * ABS(PCOT(I) / FWHMGL))) ! Original in GSAS
*        KT = MAX(1, NINT(2.5 * ABS(PCOT(I) / FWHMGL))) ! Modified by me
*        In the case of the high-resolution X-ray powder diffractometer 
*        at NIRIM, fixing KT at 1 gave the lowest Rwp.
*        IF (NBEAM .EQ. 1) THEN
*           KT = 1
*        ELSE
*           KT = MAX(1, IFIX(2.5 * ABS(PCOT(I) / FWHMGL)))
*        END IF

**       NTSIM: NT IN GSAS (SIM = SIMPSON)
         NTSIM(I) = 2 * KT + 1
*        FNORM: NORMALIZATION FACTOR
         FNORM(I) = 1.0/FLOAT(6*KT)
*        TNTSIM: TNT IN GSAS (SIM = SIMPSON)
         TNTSIM(I) = FLOAT((NTSIM(I) - 1)**2)
      END IF

*     SQSG: sigma
*     ORIGINAL EXPRESSION OF sigma
*     SQSG(I) = MAX(SQRT(SIG), 0.001)
*     SIMILAR TO THE ORIGINAL EXPRESSION, BUT 0.001 IS REDUCED 
      SQSG(I) = MAX(SQRT(SIG), 0.0000001)

*     MY EXPRESSION OF sigma     
*     SQSG(I) = SQRT(SIG)
 
      FWHG = 2.354820 * SQSG(I)
      PGL = FWHG**5
*     SUMHM = GAMMA**5
      SUMHM = PGL
*     DSDL: d(SUMHM)/d(GAM) = d(SUMHM)/d(H_kL). Derivative wrt H.k(Lorentz)
      DSDL = 0.0
*     DSDG: d(SUMHM)/d(FWHG) = d(SUMHM)/d(H_kG). Derivative wrt H.k(Gauss)
      DSDG = 0.0
      DO ITRM = 1, 5
         PGL = PGL / FWHG
         DSDL = DSDL + FLOAT(ITRM)*COFT(ITRM+1)*PGL
         DSDG = DSDG + FLOAT(6-ITRM)*COFT(ITRM)*PGL
         PGL = PGL * GAM
         SUMHM = SUMHM + COFT(ITRM+1)*PGL
      END DO
**    FWHM(I) = SUMHM**0.2, FWHM: GAMMA = H.k
      FWHM(I) = EXP(0.2*LOG(SUMHM))
*     FRAC: gamma/GAMMA
      FRAC = GAM / FWHM(I) ! H.k(Lorentz) / H.k

*     DEDF: d(eta)/d(gamma/GAMMA)
      DEDF = 0.0
      PF = 1.0
**    ETA: MIXING PARAMETER
      ETA(I) = 0.0
      DO ITRM = 1,3
         DEDF = DEDF + FLOAT(ITRM)*COFN(ITRM)*PF
         PF = PF * FRAC
         ETA(I) = ETA(I) + COFN(ITRM)*PF
      END DO
      IF (ETA(I) .LT. 0.0 .OR. ETA(I) .GT. 1.0) THEN
         IDPOS = -1
         RETURN
      END IF

*     8*LN(2) = 5.545177
      SIGP(I) = FWHM(I)*FWHM(I)/5.545177

*     DFWDG: d(GAMMA)/d(GAMMA.g)
      DFWDG(I) = 0.2 * DSDG * FWHM(I) / SUMHM
*     4*ln(2) = 2.7725887
      SIGPAP(I) = 2.7725887*DFWDG(I)/FWHMG

*     DFRDG: d(gamma/GAMMA)/d(GAMMA.g)
      DFRDG = -FRAC * DFWDG(I) / FWHM(I)
*     DEDFFG IS USED IN SUBROUTINE PVOIGT TO CALCULATE DFDS
      DEDFFG(I) = DEDF * DFRDG

*     DFWDL: d(GAMMA)/d(gamma)
      DFWDL(I) = 0.2 * DSDL * FWHM(I) / SUMHM

*     DFRDL: d(gamma/GAMMA)/d(gamma)
      DFRDL = (1.0 - FRAC*DFWDL(I))/FWHM(I)
*     DEDFFL IS USED IN SUBROUTINE PVOIGT TO CALCULATE DFDG
      DEDFFL(I) = DEDF * DFRDL 
      END

************************************************************************

      SUBROUTINE PREDE12(A,IDPOS)
*     UPDATE ARRAYS FOR THE CALCULATION OF THE SPLIT-TYPE PSEUDO-VOIGT
*     FUNCTION
      PARAMETER (NB=7000,NAP=150,NPH=8)
      INTEGER H
      REAL A(*)
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /FT/ FT1,FT2,PEAKNOR(NB),TANTHNOR(NB) 
      COMMON /G/ I
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /PRLV/ NPRINT
      COMMON /UNCONS/ NUCREF,IPHASE(90),IHKL(3,90),NOPPP(90),LPPP(NB),
     &  RCUT(2,90),NCUT
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      COMMON /ZSPFN/ PEAKSH(NB),SQH2(NB),SQH2L(NB),SQH2G(NB),ASYM(NB),
     &ETAL(NB),ETAH(NB),FNORM1(2,NB),COMLG(2,NB),COML(2,NB),COMG(2,NB),
     &COEFL(2,NB),COEFG(2,NB),DECAY(2,NB),COMPEAR(2,NB),CDX(2,NB),
     &CDA1(2,NB),CDA2(2,NB),CDW(2,NB),CDRL1(2,NB),CDRL2(2,NB),
     &CDRH1(2,NB),DMPART(2,NB)
      
      P(J) = A(KPHB(NOPH(I))+J)
      DATA TPI/0.6366198/   ! 2.0/PI
      DATA CGAU/0.9394373/  ! 2.0*SQRT(LN(2.0)/PI)
      
      CALL SPLITCOM(A,H2,IDPOS) 
      IF (IDPOS .EQ. -1) RETURN
      
*     MIXING PARAMETER (LOWER ANGLE)
      ETAL(I) = P(9) + P(10)*PEAK(I)
*     MIXING PARAMETER (HIGHER ANGLE)
      ETAH(I) = P(11) + P(12)*PEAK(I)
      
      IF (LPPP(I) .GT. 0) CALL MARGOT(A)
      
      IF (H2 .LT. 0.0 .OR. ASYM(I) .LE. 0.0 .OR. ASYM(I) .GT. 40.0 .OR. 
     &ETAL(I) .LT. 0.0 .OR. ETAL(I) .GT. 1.0 .OR. ETAH(I) .LT. 0.0 .OR. 
     &ETAH(I) .GT. 1.0) THEN
         IDPOS = -1
         IF (NPRINT .GE. 1) THEN
            WRITE(6,'(/11X,2(A,I5))')  'I = ',I,',  PHASE = ',NOPH(I)
            WRITE(6,'(11X,A,G12.6)') 'ASYM = ',ASYM(I)
            WRITE(6,'(11X,A,G12.6)') 'ETAL = ',ETAL(I)
            WRITE(6,'(11X,A,G12.6)') 'ETAH = ',ETAH(I)
         END IF
         RETURN
      END IF

*     The normalization factor is common to the original and modified
*     split-type pseudo-Voigt functions
*     NORMALIZATION FACTOR (LOWER ANGLE)
      FNORM1(1,I) = FACNOR12(ASYM(I),ETAL(I),ETAH(I))
*     NORMALIZATION FACTOR (HIGHER ANGLE)
      FNORM1(2,I) = FACNOR12(1.0/ASYM(I),ETAH(I),ETAL(I))
      
      IF (NPRFN .EQ. 1 .OR. (NPRFN .EQ. 2 .AND. LPPP(I) .EQ. 0)) THEN
*        Split-type pseudo-Voigt function of Toraya 

*        PART COMMON TO THE LORENTZ AND GAUSS COMPONENTS (LOWER ANGLE)
         COMLG(1,I) = ((1.0 + ASYM(I))/(ASYM(I)*SQH2(I)))**2
*        PART COMMON TO THE LORENTZ AND GAUSS COMPONENTS (HIGHER ANGLE)
         COMLG(2,I) = ((1.0 + ASYM(I))/SQH2(I))**2
      
*        COEFFICIENT FOR THE LORENTZ COMPONENT (LOWER ANGLE)
         COEFL(1,I) = ETAL(I)*TPI/SQH2(I)
*        COEFFICIENT FOR THE LORENTZ COMPONENT (HIGHER ANGLE)
         COEFL(2,I) = ETAH(I)*TPI/SQH2(I)
      
*        COEFFICIENT FOR THE GAUSS COMPONENT (LOWER ANGLE)
         COEFG(1,I) = (1.0 - ETAL(I))*CGAU/SQH2(I)
*        COEFFICIENT FOR THE GAUSS COMPONENT (HIGHER ANGLE)
         COEFG(2,I) = (1.0 - ETAH(I))*CGAU/SQH2(I)
      ELSE
*        Modified split-type pseudo-Voigt function

*        PART OF THE LORENTZ AND GAUSS COMPONENTS (LOWER ANGLE)
         COML(1,I) = ((1.0 + ASYM(I))/(ASYM(I)*SQH2L(I)))**2
         COMG(1,I) = -0.6931472*((1.0 + ASYM(I))/(ASYM(I)*SQH2G(I)))**2
*        PART COMMON TO THE LORENTZ AND GAUSS COMPONENTS (HIGHER ANGLE)
         COML(2,I) = ((1.0 + ASYM(I))/SQH2L(I))**2
         COMG(2,I) = -0.6931472*((1.0 + ASYM(I))/SQH2G(I))**2
         
*        COEFFICIENT FOR THE LORENTZ COMPONENT (LOWER ANGLE)
         COEFL(1,I) = ETAL(I)*TPI/SQH2L(I)
*        COEFFICIENT FOR THE LORENTZ COMPONENT (HIGHER ANGLE)
         COEFL(2,I) = ETAH(I)*TPI/SQH2L(I)
      
*        COEFFICIENT FOR THE GAUSS COMPONENT (LOWER ANGLE)
         COEFG(1,I) = (1.0 - ETAL(I))*CGAU/SQH2G(I)
*        COEFFICIENT FOR THE GAUSS COMPONENT (HIGHER ANGLE)
         COEFG(2,I) = (1.0 - ETAH(I))*CGAU/SQH2G(I)
      END IF
      END

************************************************************************

      SUBROUTINE PREDE3(A,IDPOS)
*     UPDATE ARRAYS FOR THE CALCULATION OF THE SPLIT-TYPE PEARSON VII
*     FUNCTION
      PARAMETER (NB=7000,NAP=150,NPH=8)
      INTEGER H
      REAL A(*),PBUNBO(2)
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /FT/ FT1,FT2,PEAKNOR(NB),TANTHNOR(NB) 
      COMMON /G/ I
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /PRLV/ NPRINT
      COMMON /UNCONS/ NUCREF,IPHASE(90),IHKL(3,90),NOPPP(90),LPPP(NB),
     &  RCUT(2,90),NCUT
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      COMMON /ZSPFN/ PEAKSH(NB),SQH2(NB),SQH2L(NB),SQH2G(NB),ASYM(NB),
     &ETAL(NB),ETAH(NB),FNORM1(2,NB),COMLG(2,NB),COML(2,NB),COMG(2,NB),
     &COEFL(2,NB),COEFG(2,NB),DECAY(2,NB),COMPEAR(2,NB),CDX(2,NB),
     &CDA1(2,NB),CDA2(2,NB),CDW(2,NB),CDRL1(2,NB),CDRL2(2,NB),
     &CDRH1(2,NB),DMPART(2,NB)
      
      P(J) = A(KPHB(NOPH(I))+J)
      
      CALL SPLITCOM(A,H2,IDPOS)  
      
*     PART OF THE DERIVATIVE, dm/dm0 and dm/dm1 (lower angle)
      DMPART(1,I) = -1.578/(P(9) + P(10)*PEAK(I))**2 + 0.98
*     DECAY PARAMETER (LOWER ANGLE)
      DECAY(1,I) = 1.578/(P(9) + P(10)*PEAK(I)) - 1.517 + 
     &            0.980 * (P(9) + P(10)*PEAK(I))
*     PART OF THE DERIVATIVE, dm/dm0 and dm/dm1 (higher angle)
      DMPART(2,I) = -1.578/(P(11) + P(12)*PEAK(I))**2 + 0.98
*     DECAY PARAMETER (HIGHER ANGLE)
      DECAY(2,I) = 1.578/(P(11) + P(12)*PEAK(I)) - 1.517 + 
     &            0.980 * (P(11) + P(12)*PEAK(I))
      
      IF (LPPP(I) .GT. 0) CALL MARGOT(A)
      
      IF (H2 .LT. 0.0 .OR. ASYM(I) .LE. 0.0 .OR. ASYM(I) .GT. 40.0 .OR. 
     &DECAY(1,I) .LT. 0.0 .OR. DECAY(1,I) .GT. 35.0 .OR. DECAY(2,I) 
     &.LT. 0.0 .OR. DECAY(2,I) .GT. 35.0) THEN
         IDPOS = -1
         IF (NPRINT .GE. 1) THEN
            WRITE(6,'(/11X,2(A,I5))')  'I = ',I,',  PHASE = ',NOPH(I)
            WRITE(6,'(11X,A,G12.6)') 'H**2 = ',H2
            WRITE(6,'(11X,A,G12.6)') 'ASYM = ',ASYM(I)
            WRITE(6,'(11X,A,G12.6)') 'ML   = ',DECAY(1,I)
            WRITE(6,'(11X,A,G12.6)') 'MH   = ',DECAY(2,I)
         END IF
         RETURN
      END IF

*     PARTS OF THE PEARSON VII FUNCTION (LOWER ANGLE)
      COMLG(1,I) = ((1.0 + ASYM(I))/(ASYM(I)*SQH2(I)))**2
      COMPEAR(1,I) = COMLG(1,I)*TR(DECAY(1,I))
      
*     PARTS OF THE PEARSON VII FUNCTION (HIGHER ANGLE)
      COMLG(2,I) = ((1.0 + ASYM(I))/SQH2(I))**2
      COMPEAR(2,I) = COMLG(2,I)*TR(DECAY(2,I))
      
*     NORMALIZATION FACTOR (LOWER ANGLE)
      FNORM1(1,I) = FACNOR3(ASYM(I),SQH2(I),DECAY(1,I),DECAY(2,I),
     &              PBUNBO(1))
*     NORMALIZATION FACTOR (HIGHER ANGLE)
      FNORM1(2,I) = FACNOR3(1.0/ASYM(I),SQH2(I),DECAY(2,I),DECAY(1,I),
     &              PBUNBO(2))
     
      CDX(1,I) = FNCDX(DECAY(1,I),ASYM(I),SQH2(I))
      CDX(2,I) = FNCDX(DECAY(2,I),1.0/ASYM(I),SQH2(I))
      
      CDA1(1,I) = FNCDA1(DECAY(1,I),ASYM(I),PBUNBO(1))
      CDA1(2,I) = FNCDA1(DECAY(2,I),1.0/ASYM(I),PBUNBO(2))
      
      CDA2(1,I) = FNCDA2(DECAY(1,I),ASYM(I))
      CDA2(2,I) = FNCDA2(DECAY(2,I),1.0/ASYM(I))
      
      CDW(1,I) = FNCDW(DECAY(1,I),ASYM(I),SQH2(I))
      CDW(2,I) = FNCDW(DECAY(2,I),1.0/ASYM(I),SQH2(I))
      
      CDRL1(1,I) = FNCDRL1(DECAY(1,I),ASYM(I),PBUNBO(1))
      CDRL1(2,I) = FNCDRL1(DECAY(2,I),1.0/ASYM(I),PBUNBO(2))
      
      CDRL2(1,I) = FNCDRL2(DECAY(1,I),ASYM(I),SQH2(I))
      CDRL2(2,I) = FNCDRL2(DECAY(2,I),1.0/ASYM(I),SQH2(I))
      
      CDRH1(1,I) = FNCDRH(DECAY(2,I),ASYM(I),PBUNBO(1))
      CDRH1(2,I) = FNCDRH(DECAY(1,I),1.0/ASYM(I),PBUNBO(2))
      END
      
************************************************************************
      
      SUBROUTINE SPLITCOM(A,H2,IDPOS)
*     CALCULATE ARRAYS COMMON TO THE SPLIT-TYPE FUNCTIONS
      PARAMETER (NB=7000,NPH=8,NAP=150)
      INTEGER H, HANIS
      REAL A(*)
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /ANISBR/ HANIS(NPH),KANIS(NPH),LANIS(NPH),
     &  CANIS(NPH),PCOS(NB)
      COMMON /FT/ FT1,FT2,PEAKNOR(NB),TANTHNOR(NB) 
      COMMON /G/ I
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      COMMON /ZSPFN/ PEAKSH(NB),SQH2(NB),SQH2L(NB),SQH2G(NB),ASYM(NB),
     &ETAL(NB),ETAH(NB),FNORM1(2,NB),COMLG(2,NB),COML(2,NB),COMG(2,NB),
     &COEFL(2,NB),COEFG(2,NB),DECAY(2,NB),COMPEAR(2,NB),CDX(2,NB),
     &CDA1(2,NB),CDA2(2,NB),CDW(2,NB),CDRL1(2,NB),CDRL2(2,NB),
     &CDRH1(2,NB),DMPART(2,NB)
     
      P(J) = A(KPHB(NOPH(I))+J)
      
      CALL LBROAD
      
      TTH = TAN(TH(I))
      
      SELECT CASE (NSHIFT)
         CASE (0)
            PEAKSH(I) = A(1)
         CASE (1)
            PEAKSH(I) = A(1) + A(2)*COS(PEAK(I)) + A(3)*SIN(PEAK(I)) + 
     &                  A(4)*TTH
         CASE (2)
            PEAKSH(I) = ((A(4)*PEAK(I) + A(3))*PEAK(I) + A(2))*PEAK(I) +
     &                  A(1)
         CASE (3) 
            PEAKSH(I) = ((A(4)*TTH + A(3))*TTH + A(2))*TTH + A(1)
         CASE (4)
            CALL SLESUM(PEAKNOR(I),3,A(1),PEAKSH(I))
         CASE (5)
            CALL SLESUM(TANTHNOR(I),3,A(1),PEAKSH(I))
      END SELECT
      
*     FWHM**2
*     H2 = (P(1)*TTH + P(2))*TTH + P(3)
*     PCOS(I) should be squared
      H2 = ((P(1) + P(13)*PCOS(I)**2)*TTH + P(2))*TTH + P(3) + 
     &     P(14)*(PCOS(I)/COS(TH(I)))**2
      IF (H2 .LT. 0.0) THEN
         WRITE(6,'(/11X,2(A,I5))')  'I = ',I,',  PHASE = ',NOPH(I)
         WRITE(6,'(11X,A,G12.6)') 'H**2 = ',H2
         IDPOS = -1
         RETURN
      END IF
*     FWHM
      SQH2(I) = SQRT(H2)
      
*     DEGREE OF PROFILE ASYMMETRY
      COSCTH = 1.0/SIN(TH(I))
      ASYM(I) = P(5) + P(6)*(1.4142136 - COSCTH) + 
     &          P(7)*(2.0 - COSCTH**2)
     
      END

************************************************************************

      REAL FUNCTION FNCDX(RL,A,W)
*     COMMON PART IN dF/dx

      FNCDX = 2.0*RL*((1.0 + A)/(A*W))**2*TR(RL)
      END

************************************************************************

      REAL FUNCTION FNCDA1(RL,A,PBUNBO)
*     COMMON PART 1 IN dF/dA

      FNCDA1 = 1.0/(1.0 + A) - SGAMMA(RL - 0.5)/(SGAMMA(RL)*
     &  SQRT(TR(RL))*PBUNBO)
      END

************************************************************************

      REAL FUNCTION FNCDA2(RL,A)
*     COMMON PART 2 IN dF/dA

      FNCDA2 = 2.0*RL*TR(RL)*(1.0 + A)/A**3
      END

************************************************************************

      REAL FUNCTION FNCDW(RL,A,W)
*     COMMON PART IN dF/dW

      FNCDW = 2.0*RL*((1.0 + A)/A)**2*TR(RL)/W**3
      END

************************************************************************

      REAL FUNCTION FNCDRL1(RL,A,PBUNBO)
*     COMMON PART IN dF/dRL

      FNCDRL1 = A*SGAMMA(RL - 0.5)/(PBUNBO *SGAMMA(RL)*SQRT(TR(RL)))*
     &  (2.0**(1.0/RL)*0.6931472/(2.0*TR(RL)*RL**2) - 
     &  DGDR(RL)/SGAMMA(RL))
      END

************************************************************************

      REAL FUNCTION FNCDRL2(RL,A,W)
*     COMMON PART IN dF/dRL

      FNCDRL2 = ((1.0 + A)/(A*W))**2/RL*2.0**(1.0/RL)*0.6931472
      END

************************************************************************

      REAL FUNCTION FNCDRH(RH,A,PBUNBO)
*     COMMON PART IN dF/dRH

      FNCDRH = SGAMMA(RH - 0.5)/(PBUNBO *SGAMMA(RH)*SQRT(TR(RH)))*
     &(2.0**(1.0/RH)*0.6931472/(2.0*TR(RH)*RH**2) - DGDR(RH)/SGAMMA(RH))
      END

************************************************************************

      REAL FUNCTION TR(R)
      TR = 2.0**(1.0/R) - 1.0
      END

************************************************************************

      REAL FUNCTION DGDR(X)
*     DERIVATIVE OF THE GAMMA FUNCTION (ACCURATE BETWEEN 0.1 < X < 35.0)
      PARAMETER (NDATA=13961)
      REAL XT(NDATA),YT(NDATA),EOPT(2)
      INTEGER IOPT(5)
      DATA IENTER/0/
      SAVE XT,YT,IOPT,IENTER
 
*     Initialization
      IF (IENTER .EQ. 0) THEN
         I = 0
         DO I = 1, NDATA
            XX = 0.1 + FLOAT(I-1)*0.0025
            YT(I) = SGAMMA(XX)
         END DO
         XT(1) = 0.1
         XT(2) = 0.0025
         IOPT(2) = 3
         IOPT(3) = 2
         IOPT(4) = 1
         IOPT(5) = 0
         IENTER = 1
      END IF
      LUP = 3
      CALL SILUP(X,Y,13961,XT,YT,3,LUP,IOPT,EOPT)
      DGDR = EOPT(2)
      END

************************************************************************

      SUBROUTINE MARGOT(A)
*     SET UP PRIMARY PROFILE PARAMETERS FOR A RELAXED REFLECTION
      PARAMETER (NT=999,NB=7000)
      REAL A(*)
      INTEGER H
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /G/ I
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /UNCONS/ NUCREF,IPHASE(90),IHKL(3,90),NOPPP(90),LPPP(NB),
     &  RCUT(2,90),NCUT
      COMMON /ZSPFN/ PEAKSH(NB),SQH2(NB),SQH2L(NB),SQH2G(NB),ASYM(NB),
     &ETAL(NB),ETAH(NB),FNORM1(2,NB),COMLG(2,NB),COML(2,NB),COMG(2,NB),
     &COEFL(2,NB),COEFG(2,NB),DECAY(2,NB),COMPEAR(2,NB),CDX(2,NB),
     &CDA1(2,NB),CDA2(2,NB),CDW(2,NB),CDRL1(2,NB),CDRL2(2,NB),
     &CDRH1(2,NB),DMPART(2,NB)
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM

      IS = NOPPP(LPPP(I))
      SELECT CASE (NPRFN)
         CASE (1)
*           Split-type pseudo-Voigt function
            IF (L12(I) .EQ. 0) THEN
*              Monochromatic beam
*              Initial PPP's are set at values calculated from secondary 
*              profile parematers when they are zero, ID(IS) = 0, or 
*              NMODE = 1 (simulation; ID's are set at zero in 
*              SUBROUTINE RDPARA).
               IF (A(IS) .EQ. 0.0 .OR. IDSAVE(IS) .EQ. 3) 
     &         A(IS) = SQH2(I)
               IF (A(IS+1) .EQ. 0.0 .OR. IDSAVE(IS+1) .EQ. 3) 
     &         A(IS+1) = ASYM(I)
               IF (A(IS+2) .EQ. 0.0 .OR. IDSAVE(IS+2) .EQ. 3) 
     &         A(IS+2) = ETAL(I)
               IF (A(IS+3) .EQ. 0.0 .OR. IDSAVE(IS+3) .EQ. 3) 
     &         A(IS+3) = ETAH(I)
*              UPDATE PPP'S
               SQH2(I) = A(IS)
               ASYM(I) = A(IS+1)
               ETAL(I) = A(IS+2)
               ETAH(I) = A(IS+3)
            ELSE
*              The primary profile parameters for a K_Alpha-2 peak is 
*              set equal to those for the corresponding K_Alpha-1 peak
               SQH2(I) = SQH2(L12(I))
               ASYM(I) = ASYM(L12(I))
               ETAL(I) = ETAL(L12(I))
               ETAH(I) = ETAH(L12(I))
            END IF
         CASE (2)
*           Modified split-type pseudo-Voigt function
            IF (L12(I) .EQ. 0) THEN
*              K-alpha1 or neutron
*              Initial PPP's are set at values calculated from secondary 
*              profile parematers when they are zero, ID(IS) = 0, or 
*              NMODE = 1 (simulation; ID's are set at zero in SUBROUTINE
*              RDPARA).
               IF (A(IS) .EQ. 0.0 .OR. IDSAVE(IS) .EQ. 3) 
     &         A(IS) = SQH2(I)
*              ASSUME THAT W1 = W2
               IF (A(IS+1) .EQ. 0.0 .OR. IDSAVE(IS+1) .EQ. 3) 
     &         A(IS+1) = SQH2(I)
               IF (A(IS+2) .EQ. 0.0 .OR. IDSAVE(IS+2) .EQ. 3) 
     &         A(IS+2) = ASYM(I)
               IF (A(IS+3) .EQ. 0.0 .OR. IDSAVE(IS+3) .EQ. 3) 
     &         A(IS+3) = ETAL(I)
               IF (A(IS+4) .EQ. 0.0 .OR. IDSAVE(IS+4) .EQ. 3) 
     &         A(IS+4) = ETAH(I)
*              UPDATE PPP'S
               SQH2L(I) = A(IS)
               SQH2G(I) = A(IS+1)
               ASYM(I) = A(IS+2)
               ETAL(I) = A(IS+3)
               ETAH(I) = A(IS+4)
            ELSE
*              The primary profile parameters for a K_Alpha-2 peak is 
*              set equal to those for the corresponding K_Alpha-1 peak
               SQH2L(I) = SQH2L(L12(I))
               SQH2G(I) = SQH2G(L12(I))
               ASYM(I) = ASYM(L12(I))
               ETAL(I) = ETAL(L12(I))
               ETAH(I) = ETAH(L12(I))
            END IF
         CASE (3)
*           Split-type Pearson VII function
            IF (L12(I) .EQ. 0) THEN
*              K-alpha1 or neutron
*              Initial PPP's are set at values calculated from secondary 
*              profile parematers when they are zero, ID(IS) = 0, or 
*              NMODE = 1 (simulation; ID's are set at zero in SUBROUTINE
*              RDPARA).
               IF (A(IS) .EQ. 0.0 .OR. IDSAVE(IS) .EQ. 3) 
     &         A(IS) = SQH2(I)
               IF (A(IS+1) .EQ. 0.0 .OR. IDSAVE(IS+1) .EQ. 3) 
     &         A(IS+1) = ASYM(I)
               IF (A(IS+2) .EQ. 0.0 .OR. IDSAVE(IS+2) .EQ. 3) 
     &         A(IS+2) = DECAY(1,I)
               IF (A(IS+3) .EQ. 0.0 .OR. IDSAVE(IS+3) .EQ. 3) 
     &         A(IS+3) = DECAY(2,I)
*              UPDATE PPP'S
               SQH2(I) = A(IS)
               ASYM(I) = A(IS+1)
               DECAY(1,I) = A(IS+2)
               DECAY(2,I) = A(IS+3)
            ELSE
*              The primary profile parameters for a K_Alpha-2 peak is 
*              set equal to those for the corresponding K_Alpha-1 peak
               SQH2(I) = SQH2(L12(I))
               ASYM(I) = ASYM(L12(I))
               DECAY(1,I) = DECAY(1,L12(I))
               DECAY(2,I) = DECAY(2,L12(I))
            END IF
      END SELECT
      END

************************************************************************

      FUNCTION FACNOR12(A,ETAL,ETAH)
*     NORMALIZATION FACTOR FOR THE SPLIT-TYPE PSEUDO-VOIGHT FUNCTION
      DATA RPILN2/1.4756646/ ! SQRT(PI*LN(2.0))

      COM = ETAH + (1.0 - ETAH)*RPILN2
      FACNOR12 = (1.0 + A)*COM/(ETAL + (1.0 - ETAL)*RPILN2 + A*COM)
      END

************************************************************************

      FUNCTION FACNOR3(A,W,RL,RH,PBUNBO)
*     NORMALIZATION FACTOR FOR THE SPLIT-TYPE PEARSON VII FUNCTION

*     PARTS OF THE NORMALIZATION FACTOR
      PBUNBO = A*SGAMMA(RL - 0.5)/(SGAMMA(RL)*SQRT(TR(RL))) +
     &           SGAMMA(RH - 0.5)/(SGAMMA(RH)*SQRT(TR(RH)))
*     SQRT(PI) = 1.7724539
      FACNOR3 = 2.0*(1.0 + A)/(1.7724539*W*PBUNBO)
      END

************************************************************************

      SUBROUTINE EXCTES(A,IPOINT,K)
C     TESTED EXCHANGE
      DIMENSION A(*),IPOINT(*)
      DO 10 I=1,K
         IPOINT(I)=I
   10 CONTINUE
    1 ISW=0
      DO 20 I=2,K
         IF(A(I-1).LE.A(I)) GO TO 20
         W=A(I)
         A(I)=A(I-1)
         A(I-1)=W
         IW=IPOINT(I)
         IPOINT(I)=IPOINT(I-1)
         IPOINT(I-1)=IW
         ISW=1
   20 CONTINUE
      IF(ISW.NE.0) GO TO 1
      END

************************************************************************

      SUBROUTINE MLAT
*     RETURN THE MULTIPLICITY OF THE CELL
      PARAMETER (NPH=8,NAP=150,NB=7000)
      COMMON /EXC/ LSPSYM(NPH),NSSYM(14,NPH)
      COMMON /H/ NQ(NPH),NCENTR(NPH),LAUEG(NPH)
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /T/ FF(NB),COEF(NPH)

      DO 10 IP=1,NPHASE
         NCON=NSSYM(7,IP)
         IF (NCON .EQ. 0 .OR. NCON .EQ. 21) THEN
*           P OR R LATTICE
            IF(NCENTR(IP).EQ.1) THEN
               COEF(IP)=2.0
            ELSE
               COEF(IP)=1.0
            ENDIF
         ELSE IF (NCON .EQ. 4 .OR. NCON .EQ. 5 .OR. NCON .EQ. 6 .OR.
     &   NCON .EQ. 13) THEN
*           A, B, C, OR I LATTICE
            IF(NCENTR(IP).EQ.1) THEN
               COEF(IP)=4.0
            ELSE
               COEF(IP)=2.0
            ENDIF
         ELSE IF (NCON .EQ. 7) THEN
*           F LATTICE
            IF(NCENTR(IP).EQ.1) THEN
               COEF(IP)=8.0
            ELSE
               COEF(IP)=4.0
            ENDIF
         ELSE IF (NCON .EQ. 14 .OR. NCON .EQ. 15) THEN
*           R LATTICE (HEXAGONAL SETTING)
            IF(NCENTR(IP).EQ.1) THEN
               COEF(IP)=6.0
            ELSE
               COEF(IP)=3.0
            ENDIF
         ENDIF
   10 CONTINUE
      END

************************************************************************

      SUBROUTINE LIST2(G,ID,DEG,XINT,YFIT,NTERMS,NRFN,NPTS,
     &  ISCALE,IDIF,NPRINT)
C     (1) LIST OF CORRELATION MATRIX
C     (2) LIST OF 2*THETA, Y(O), AND Y(C)
C     (3) LINE PRINTER PLOT OF Y(O), Y(C), AND Y(O)-Y(C)
      DOUBLE PRECISION G(*)
      CHARACTER*1 ISPACE,IDOT,ISTAR,IPLUS,IMINUS,IBAR
      PARAMETER (NT=999)
      DIMENSION ITEMP(NT),ID(*),DEG(*),XINT(*),YFIT(*),KS(19)
      COMMON /L/ DEG1,DEG2
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /X/ STEP
      DATA RD2/57.29578/
      DATA ISPACE,IDOT,ISTAR,IPLUS,IMINUS,IBAR/' ','.','*','+','0','-'/

      IF (NPRINT .EQ. 2) THEN
C        CORRELATION MATRIX
         WRITE(6,'(//,11X,A)') 'Correlation matrix'
         L=1
         DO 12 J=1,NTERMS
            IF(ID(J).EQ.1) THEN
               ITEMP(L)=J
               L=L+1
            ENDIF
   12    CONTINUE
         IA=1
         LIM=19
    1    LIM=MIN0(NRFN,LIM)
         WRITE(6,100) (ITEMP(I),I=IA,LIM)
  100    FORMAT(//13X,19I6)
         DO 20 I=1,NRFN
            L=0
            DO 15 J=IA,LIM
               L=L+1
               IF(I.GE.J) THEN
                  KS(L)=NINT(100.0D0*G((I-1)*I/2+J)/SQRT(G((I-1)*I/2+I)*
     &            G((J-1)*J/2+J)))
               ELSE
                  KS(L)=NINT(100.0D0*G((J-1)*J/2+I)/SQRT(G((J-1)*J/2+J)*
     &            G((I-1)*I/2+I)))
               ENDIF
   15       CONTINUE
            WRITE(6,110) ITEMP(I),(KS(J),J=1,L)
   20    CONTINUE
  110    FORMAT(8X,I5,19I6)
         IF(LIM.LT.NRFN) THEN
            IA=LIM+1
            LIM=LIM+19
            GO TO 1
         ENDIF
      END IF
C
      IF (NPRINT .EQ. 2) THEN
C        LIST 2*THETA, Y(O), AND Y(C)
         NP=NPTS/300
         WRITE(6,120) (((DEG(300*N-360+60*I+J)*RD2,
     &   NINT(XINT(300*N-360+60*I+J)),NINT(YFIT(300*N-360+60*I+J)),
     &   I=1,5),J=1,60),N=1,NP)
  120    FORMAT(('1'/5X,5('  2O',6X,'Yobs',3X,'Ycalc',4X)/'+',4X,
     &   5('   -',22X)/,60(3X,5(1X,F7.3,1X,2(I7,1X),1X),/)))
         NPTS2=NPTS-NP*300
         IF(NPTS2.EQ.0) RETURN
         NCOL=NPTS2/60
         WRITE(6,130)
  130    FORMAT('1'/5X,5('  2O',6X,'Yobs',3X,'Ycalc',4X)/'+',4X,
     &   5('   -',22X))
         NLINES=NPTS2-NCOL*60
         IF(NLINES.NE.0) THEN
            NCOL1=NCOL+1
            NP=300*NP-60
            DO 30 J=1,NLINES
               WRITE(6,140) (DEG(NP+60*I+J)*RD2,NINT(XINT(NP+60*I+J)),
     &         NINT(YFIT(NP+60*I+J)),I=1,NCOL1)
   30       CONTINUE
  140       FORMAT(3X,5(1X,F7.3,1X,2(I7,1X),1X))
            NP=NP+NLINES
         ENDIF
         DO 35 J=1,60-NLINES
            WRITE(6,140) (DEG(NP+60*I+J)*RD2,NINT(XINT(NP+60*I+J)),
     &      NINT(YFIT(NP+60*I+J)),I=1,NCOL)
   35    CONTINUE
      END IF
      END

************************************************************************

      FUNCTION ABSCOR(SIGA,SIGI,ATOMNO,NABS,SIN2TH,RADIUS,MUR)
C     ABSORPTION CORRECTION FOR A CYLINDRICAL CONTAINER AT THETA
C     ALL THE QUANTITIES ARE EXPRESSED IN SI DERIVED UNITS
C     SIN2TH: (SIN(THETA))**2
      PARAMETER (NA=15)
      REAL SIGA(*),SIGI(*),ATOMNO(*),MU,MUR
      COMMON /E/ AA(4,NA),BB(4,NA),C(NA)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD

      IF (NBEAM .EQ. 0) THEN
*        NEUTRON DIFFRACTION
         MU=0.0
         DO 10 J=1,NABS
C           SIGA(J): ABSORPTION CROSS SECTION/M**2 AT 0.179817E-9 M
C           SIGI(J): INCOHERENT SCATTERING CROSS SECTION/M**2
C           C(J): COHERENT SCATTERING LENGTH/FM
C           4.0*PI*(C(J)*1.0E-15)**2: COHERENT SCATTERING CROSS SECTION
C           THE SUM OF INCOHERENT AND COHERENT SCATTERING IS THE SELF-
C           SHIELDING.
C           THE COHERENT SCATTERING IS ASSUMED TO BE ISOTROPIC.
*           SIGABS: ABSORPTION CROSS SECTION AT A WAVELENGTH OF XLMD
            SIGABS = SIGA(J)*(XLMD*1E-10)/0.179817E-9
            MU = MU + ATOMNO(J)*(SIGABS + SIGI(J) +
     &      4.0*3.141593*(C(J)*1.0E-15)**2)
   10    CONTINUE
         MUR = MU * RADIUS
      END IF

*     ANALYTICAL EXPRESSION FOR ABSROPTION CORRECTION
*     A. W. HEWAT, ACTA CRYSTALLOGR., SECT. A, 35, 248 (1979).
      ABSCOR = 
     & EXP((-1.7133 + 0.0368*SIN2TH + (0.0927 + 0.3750*SIN2TH)*MUR)*MUR)
      END

************************************************************************

      SUBROUTINE BRANG
C     CALCULATE THE D SPACING FOR THE I'TH REFLECTION AND FUNCTIONS
C     AND THE FUNCTIONS OF D
C     PEAK: TWO-THETA
C       TH: THETA
C        D: D SPACING
C     RLV2: (SIN(THETA)/LAMBDA)**2
C      FLP: LORENTZ AND POLARIZATION FACTORS
      INTEGER H
      PARAMETER (NB=7000,NPH=8)
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /G/ I
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /MC/ CTHM,NTRAN,PCOR
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /TENS/ G(6,NPH),GG(6,NPH)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD

      N=NOPH(I)
      Z = FLOAT(H(I)*H(I))*GG(1,N) + FLOAT(K(I)*K(I))*GG(2,N)+
     &  FLOAT(L(I)*L(I))*GG(3,N) + 2.0*(FLOAT(K(I)*L(I))*GG(4,N) +
     &  FLOAT(L(I)*H(I))*GG(5,N) + FLOAT(H(I)*K(I))*GG(6,N))
      IF (Z .LE. 0.0) THEN
         D(I)=0.0
         RETURN
      ENDIF
      D(I)=1.0/SQRT(Z)
      RLV2(I)=0.25*Z
C
*     CALCULATE WAVELENGTH-DEPENDENT QUANTITIES
      SELECT CASE (L12(I))
         CASE (0)
            IF (XLMDH/D(I) .LT. 1.0) THEN
               TH(I) = ASIN(XLMDH/D(I))
            ELSE
               D(I) = 0.0
               RETURN
            END IF
         CASE DEFAULT
            IF (XLMD2H/D(I) .LT. 1.0) THEN
               TH(I) = ASIN(XLMD2H/D(I))
            ELSE
               D(I) = 0.0
               RETURN
            END IF
      END SELECT
      PEAK(I) = TH(I) + TH(I)

*     LORENTZ AND POLARIZATION FACTORS
*     Refer to EQ. (3.25) in the manual of Fullprof
      SELECT CASE (NBEAM)
         CASE (0)
*           PCOR = 0.0
            FLP(I) = 0.5/(SIN(TH(I))**2*COS(TH(I)))
         CASE (1)
            SELECT CASE (NTRAN)
               CASE (0,1,3)
*                 Reflection or Debye-Scherrer geometry
*                 PCOR = 0.5
                  FLP(I) = (0.25 + 0.25*CTHM*COS(PEAK(I))**2)/
     &                     (SIN(TH(I))**2*COS(TH(I)))
               CASE (2)
*                 Transmission flat-plate geometry
*                 Refer to Eqs. (3.26) and (3.27) in the manual of Fullprof
*                 Polarization factor
                  FP = PCOR*(1.0+SQRT(CTHM)*COS(PEAK(I))**2)/
     &            (1.0+SQRT(CTHM)) + (1.0-PCOR)*
     &            (1.0+CTHM*COS(PEAK(I))**2)/(1.0+CTHM)
                  FLP(I) = FP/SIN(PEAK(I))
            END SELECT
         CASE (2)
*           PCOR and CTHM have been read in the main program.
            FLP(I) = (0.5 - 0.5*PCOR*(1.0 - CTHM*COS(PEAK(I))**2))/
     &               (SIN(TH(I))**2*COS(TH(I)))
      END SELECT
      END

************************************************************************

      SUBROUTINE BRANGPF(A)
C     CALCULATE THE D SPACING FOR THE I'TH REFLECTION AND FUNCTIONS
C     AND THE FUNCTIONS OF D
C     PEAK: TWO-THETA
C       TH: THETA
C        D: D SPACING
C      FLP: LORENTZ AND POLARIZATION FACTORS
      REAL A(*)
      INTEGER H
      PARAMETER (NB=7000,NPH=8)
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /G/ I
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /MC/ CTHM,NTRAN,PCOR
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
      COMMON /UNCONS/ NUCREF,IPHASE(90),IHKL(3,90),NOPPP(90),LPPP(NB),
     &  RCUT(2,90),NCUT
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM

      IF (L12(I) .EQ. 0) THEN
*        MONOCHROMATIC BEAM
         SELECT CASE (NPRFN)
            CASE (1,3)
               PEAK(I) = A(NOPPP(LPPP(I)) + 5) 
            CASE (2)
               PEAK(I) = A(NOPPP(LPPP(I)) + 6) 
         END SELECT
         TH(I) = 0.5*PEAK(I)
         D(I) = XLMDH/SIN(TH(I))
         RLV2(I) = 0.25/D(I)**2
      ELSE
*        K-ALPHA2 RADIATION
*        Cf. H. Toraya, Jikken Kagaku Kohza 10, 4th ed., Maruzen, p. 313.
*        RLAMBD = lambda(K-alpha2)/lambda(K-alpha1)
         TH(I) = ASIN(RLAMBD*SIN(TH(L12(I))))
         PEAK(I) = TH(I) + TH(I)
         D(I) = D(L12(I))
         RLV2(I) = RLV2(L12(I))
      END IF

*     LORENTZ AND POLARIZATION FACTORS
      SELECT CASE (NBEAM)
         CASE (1)
            FLP(I) = (1.0 + CTHM*COS(PEAK(I))**2)/
     &               (SIN(TH(I))**2*COS(TH(I)))
         CASE (0,2)
            FLP(I) = 1.0/(SIN(TH(I))**2*COS(TH(I)))
      END SELECT
      END

************************************************************************

      SUBROUTINE PROR(A)
C     PREFERRED-ORIENTATION CORRECTION OF THE I'TH REFLECTION
      PARAMETER (NB=7000,NPH=8,NAP=150,NS=48)
      INTEGER H,HP,R,RX,RY,RZ,IH(49),IK(49),IL(49)
      REAL A(*)
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /B/ U(NB),NREFL
*Rev 1.0h 2000.12.13 Izumi
      COMMON /BIJVOET/ IPAIR(NB),LPAIR(NPH)
      COMMON /DERPR/ POP2(NB),COMPO(NB),CS(NB),CCC(6,NB)
      COMMON /G/ I
      COMMON /H/ NQ(NPH),NCENTR(NPH),LAUEG(NPH)
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /RT/ T(3,NS,NPH),HT(NS,NB),NSYM(NPH),NSITE(NPH),
     &  IDSYM(NS,NAP,NPH),R(3,3,NS,NPH),RX(NS,NB),RY(NS,NB),RZ(NS,NB)
      COMMON /TENS/ G(6,NPH),GG(6,NPH)
      COMMON /V/ NPROR(NPH),HP(NPH),KP(NPH),LP(NPH),LSUM(NPH)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
      DATA DEGR,PI/1.570796,3.141593/
      P(J)=A(KPHB(N)+J)

      N = NOPH(I)
      IF (NPROR(N) .EQ. 0) THEN
         PROR1(I)=1.0
         RETURN
      ENDIF

*     GENERATE THE INDICES OF EQUIVALENT REFLECTIONS
      IF (LSUM(N) .EQ. 1) THEN
         IMULT=1
         DO 30 J = 1, NSYM(N)
*           CF. EQ. (5.5-18), P. 243 IN SAKURAI'S BOOK
            IBEFOR = IMULT
            IH(IMULT)=H(I)*R(1,1,J,N)+K(I)*R(2,1,J,N)+L(I)*R(3,1,J,N)
            IK(IMULT)=H(I)*R(1,2,J,N)+K(I)*R(2,2,J,N)+L(I)*R(3,2,J,N)
            IL(IMULT)=H(I)*R(1,3,J,N)+K(I)*R(2,3,J,N)+L(I)*R(3,3,J,N)
            DO JJ = 1, IMULT-1
               IF (IH(JJ) .EQ. IH(IMULT) .AND. IK(JJ) .EQ. IK(IMULT) 
*Rev 1.0h 2000.12.13 Izumi {
C    &         .AND. IL(JJ) .EQ. IL(IMULT)) GO TO 2
     &         .AND. IL(JJ) .EQ. IL(IMULT)) GO TO 30
            END DO                        
            IMULT = IMULT + 1
*           The following part was deleted because hkl and -h-k-l give
*           the same VV value (see below) because only (cos(phi))**2 = 
*           1 - (sin(phi))**2 is contained in VV.
C   2       IF (NCENTR(N) .EQ. 1 .OR. LPAIR(N) .EQ. 1) THEN
*              Indices derived from inverted equivalent positions
*              in a space group with an inversion center at the origin
C              IH(IMULT) = -IH(IBEFOR)
C              IK(IMULT) = -IK(IBEFOR)
C              IL(IMULT) = -IL(IBEFOR)
C              DO JJ = 1, IMULT-1
C                 IF (IH(JJ) .EQ. IH(IMULT) .AND. IK(JJ) .EQ. IK(IMULT) 
C    &            .AND. IL(JJ) .EQ. IL(IMULT)) GO TO 30
C              END DO                       
C              IMULT = IMULT + 1
C           END IF
   30    CONTINUE
*        MULT: Multiplicity for averaging the preferred-orientation
*              function.  This may differ from the multiplicity of a
*              reflection
         MULT=IMULT-1
* }
      END IF

*     MARCH-DOLLASE FUNCTION
      IF (NPROR(N) .EQ. 3) THEN
         IF (LSUM(N) .EQ. 0) THEN
            TT = (COSPHI(N,H(I),K(I),L(I)))**2
            Q = P(NPO1X)**2*TT + (1.0 - TT)/P(NPO1X)
            PROR1(I) = 1.0/Q**1.5
            DPROR(I) = - 1.5 * PROR1(I) / Q * 
     &                 (2.0*P(NPO1X)*TT - (1.0 - TT)/P(NPO1X)**2)
         ELSE
            PROR1(I) = 0.0
            DPROR(I) = 0.0
            DO 35 J=1,MULT
               TT = (COSPHI(N,IH(J),IK(J),IL(J)))**2
               Q = P(NPO1X)**2*TT + (1.0 - TT)/P(NPO1X)
               VV = 1.0/Q**1.5
               PROR1(I) = PROR1(I) + VV
               DPROR(I) = DPROR(I) - 1.5 * VV / Q * 
     &         (2.0*P(NPO1X)*TT - (1.0 - TT)/P(NPO1X)**2)
   35       CONTINUE
            PROR1(I) = PROR1(I)/FLOAT(MULT)
            DPROR(I) = DPROR(I)/FLOAT(MULT)
            IF (NCENTR(N) .EQ. 0 .AND. NBEAM .GE. 1)
     &      DPROR(I) = 0.5*DPROR(I)
         END IF
         RETURN
      END IF

*     TORAYA-MARUMO FUNCTION 
      IF (LSUM(N) .EQ. 0) THEN
         TT = COSPHI(N,H(I),K(I),L(I))
         IF (TT .GE. 1.0) THEN
            PHI=0.0
         ELSE IF (TT .LE. -1.0) THEN
            PHI=PI
         ELSE
            PHI=ACOS(TT)
         END IF
         IF (PHI .GT. DEGR) PHI=PI-PHI
         IF (NPROR(N) .EQ. 2) PHI=DEGR-PHI
         GAUSS = EXP(-P(NPO2X)*PHI**2)
         PROR1(I) = P(NPO1X) + (1.0-P(NPO1X))*GAUSS
         DPROR(I) = 1.0 - GAUSS
         POP2(I) = -(1.0-P(NPO1X))*GAUSS*PHI**2
      ELSE
         PROR1(I)=0.0
         DPROR(I)=0.0
         POP2(I)=0.0
         DO 50 J=1,MULT
            TT = COSPHI(N,IH(J),IK(J),IL(J))
            IF (TT .GE. 1.0) THEN
               PHI=0.0
            ELSE IF (TT .LE. -1.0) THEN
               PHI=PI
            ELSE
               PHI=ACOS(TT)
            END IF
            IF (PHI .GT. DEGR) PHI=PI-PHI
            IF (NPROR(N) .EQ. 2) PHI=DEGR-PHI
            GAUSS = EXP(-P(NPO2X)*PHI**2)
            PROR1(I) = PROR1(I) + P(NPO1X) + (1.0-P(NPO1X))*GAUSS
            DPROR(I) = DPROR(I) + 1.0 - GAUSS
            POP2(I) = POP2(I) - (1.0-P(NPO1X))*GAUSS*PHI**2
   50    CONTINUE
         PROR1(I)=PROR1(I)/FLOAT(MULT)
         DPROR(I)=DPROR(I)/FLOAT(MULT)
         POP2(I)=POP2(I)/FLOAT(MULT)
      END IF
      END

************************************************************************

      REAL FUNCTION COSPHI(N,H,K,L)
*     COSPHI: COS(PHI), WHERE PHI IS THE ACCUTE ANGLE BETWEEN 
*             (HP)A* + (KP)B* +(LP)C* AND HA* + KB* +LC*
      INTEGER H,HP
      PARAMETER (NB=7000,NPH=8)
      COMMON /G/ I
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /TENS/ G(6,NPH),GG(6,NPH)
      COMMON /V/ NPROR(NPH),HP(NPH),KP(NPH),LP(NPH),LSUM(NPH)

      APROR = FLOAT(H*HP(N))*GG(1,N) + FLOAT(K*KP(N))*GG(2,N) +
     &        FLOAT(L*LP(N))*GG(3,N) + FLOAT(K*LP(N)+L*KP(N))*
     &        GG(4,N) + FLOAT(H*LP(N)+L*HP(N))*GG(5,N) +
     &        FLOAT(H*KP(N)+K*HP(N))*GG(6,N)
      COSPHI = APROR*D(I)/SQRT(CPROR(N))
      END

************************************************************************

      SUBROUTINE LBROAD
*     CALCULATE PCOS WHICH IS THE COSINE OF THE ACCUTE ANGLE BETWEEN 
*     (HANIS)A* + (KANIS)B* + (LANSI)C* AND HA* + KB* +LC*
      INTEGER H,HANIS
      PARAMETER (NB=7000,NPH=8,NAP=150)
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /ANISBR/ HANIS(NPH),KANIS(NPH),LANIS(NPH),
     &  CANIS(NPH),PCOS(NB)
      COMMON /G/ I
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /TENS/ G(6,NPH),GG(6,NPH)

      N=NOPH(I)
      ANIS = FLOAT(H(I)*HANIS(N))*GG(1,N) + FLOAT(K(I)*KANIS(N))*GG(2,N)
     &  + FLOAT(L(I)*LANIS(N))*GG(3,N) + 
     &  FLOAT(K(I)*LANIS(N) + L(I)*KANIS(N))*GG(4,N) +
     &  FLOAT(H(I)*LANIS(N) + L(I)*HANIS(N))*GG(5,N) +
     &  FLOAT(H(I)*KANIS(N) + K(I)*HANIS(N))*GG(6,N)

      PCOS(I) = ABS(ANIS/SQRT(CANIS(N)))*D(I)
      END

************************************************************************

      SUBROUTINE PACK(A,C,NTERMS)
C     A: UNPACKED ARRAY ---> PACKED ARRAY
C     C: RESERVE UNPACKED ARRAY
      PARAMETER (NT=999)
      REAL A(*),C(*)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      DO 10 J=1,NTERMS
         C(J)=A(J)
         IF(ID(J).EQ.1) A(IR(J))=A(J)
   10 CONTINUE
      END

************************************************************************

      SUBROUTINE UNPACK(B,C,NTERMS)
C     B: PACKED ARRAY ---> UNPACKED ARRAY
C     C: RESERVE UNPACKED ARRAY
      PARAMETER (NT=999)
      REAL B(*),C(*)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      DO 10 J=1,NTERMS
         IF(ID(J).EQ.1) C(J)=B(IR(J))
   10 CONTINUE
      DO 20 J=1,NTERMS
         B(J)=C(J)
   20 CONTINUE
      END

************************************************************************

      SUBROUTINE PENFN(ALPHA,BETA,A,NTERMS)
C     CONTRIBUTION OF THE PENALTY FUNCTION IS ADDED TO ARRAYS ALPHA AND
C     BETA
      PARAMETER (NR=400,NT=999,NCS=400)
      DOUBLE PRECISION ALPHA,BETA
      DIMENSION ALPHA(*),BETA(*),A(*),DCNSTR(NR),DCPM(NT)
      COMMON /CNSTR/ C1(NCS),C2(10,NCS),IN(10,NCS),NCNPAR(NCS),
     &  NACNS(NCS),NCNSTR
      COMMON /FL0/ NC,TK,CONV,NCONV,FINC,SUMPEN,SUMP1,FACTOR
      COMMON /NONLIN/ NLESQ,NAUTO,IDCYCL,NCHNG,NPAR(50),IPAR(NR,50)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN

      DO 80 I=1,NC
         FCNSTR=CON(I,A)
         IF(FCNSTR.EQ.0.0) GO TO 80
C        DCPM: DERIVATIVES OF CONSTRAINTS WITH RESPECT TO PARAMETERS
C        FOR WHICH ID(J)<>0.  CENTERED DIVIDED DIFFERENCES ARE USED.
         DO 10 J=1,NTERMS
            IF(ID(J).NE.0) THEN
               ASAVE=A(J)
               IF(ASAVE.EQ.0.0) THEN
                  DELTAA=0.002
               ELSE
                  DELTAA=ABS(ASAVE)*0.002
               ENDIF
               A(J)=ASAVE-DELTAA
               F1=CON(I,A)
               A(J)=ASAVE+DELTAA
               F2=CON(I,A)
               DCPM(J)=(F2-F1)/(DELTAA+DELTAA)
               A(J)=ASAVE
            ENDIF
   10    CONTINUE
C
C        DCNSTR: DERIVATIVES WITH RESPECT PARAMETERS FOR WHICH ID(J)=1
         DO 40 J=1,NTERMS
            IF(ID(J).EQ.1) THEN
               DCNSTR(IR(J))=DCPM(J)
               DO 30 J2=1,NCNSTR
                  DO 20 J3=1,NCNPAR(J2)
                     IF(IN(J3,J2).EQ.J) DCNSTR(IR(J))=DCNSTR(IR(J))+
     &               C2(J3,J2)*DCPM(NACNS(J2))
   20             CONTINUE
   30          CONTINUE
            ENDIF
   40    CONTINUE
         DO 60 J=1,NTERMS
            IF(ID(J).EQ.1) THEN
               IF(IDCYCL.EQ.1) BETA(IR(J))=BETA(IR(J))-TK*FCNSTR*
     &         DCNSTR(IR(J))
               DO 50 K=1,J
                  JK=(IR(J)-1)*IR(J)/2+IR(K)
                  IF(ID(K).EQ.1) ALPHA(JK)=ALPHA(JK)+
     &            TK*DCNSTR(IR(J))*DCNSTR(IR(K))
   50          CONTINUE
            ENDIF
   60    CONTINUE
   80 CONTINUE
      END

************************************************************************

      SUBROUTINE INPCON(NC,NSCONS,EXPCTD,DEVDA,XYZC,NPHCON,IXCON)
*     READ NONLINEAR CONSTRAINTS IMPOSED ON INTERATOMIC DISTANCES
*     AND BOND ANGLES
      PARAMETER (MAXNL=1000)
      INTEGER NSCONS(4,MAXNL),NPHCON,IXCON(*)
      REAL EXPCTD(*),DEVDA(*),XYZC(4,3,3,MAXNL)
      CHARACTER LINE*80

      WRITE(6,'(//11X,A)') 'Nonlinear constraints imposed on '//
     &'interatomic diatances and/or bond angles'
      WRITE(6,'(15X,A)') 'I    NSCONS     EXPCTD    DEVDA'
      READ(4,'(A)') LINE
      I=1
      DO WHILE (LINE(1:1) .NE. '}')
         IF (I .GT. MAXNL) 
     &   CALL JOBEND('Too large number of nonlinear constraints')
         READ(LINE,*) NSCONS(1,I),EXPCTD(I),DEVDA(I)
         WRITE(6,'(11X,I5,I10,F11.4,F9.4)') 
     &   I,NSCONS(1,I),EXPCTD(I),DEVDA(I)
         READ(4,'(A)') LINE
         I=I+1
      END DO
      NC=I-1
      CALL ELEONORA(NC,NSCONS,XYZC,NPHCON,IXCON)
      END

************************************************************************

      SUBROUTINE ELEONORA(NSELEC,NSCONS,XYZC,NPHCON,IXCON)
*     DETERMINE ATOM NUMBERS AND COEFFICIENTS FROM WHICH FRACTIONAL 
*     COORDINATES ARE CALCULATED

*     1) INPUT
*     NSELC: TOTAL NUMBER OF USER-SELECTED INTERATOMIC DISTANCES AND 
*            BOND ANGLES 
*     NSCONS(1,I): SERIES NUMBER OF INTERATOMIC DISTANCES AND BOND
*                  ANGLES ON WHICH THE I'TH NONLINEAR CONSTRAINT IS 
*                  IMPOSED

*     2) OUTPUT
*     NSCONS(2,I)-NSCONS(J,I): SERIES NUMBERS OF ATOMS RELATED TO THE 
*                              INTERATOMIC DISTANCE (J=3) OR THE BOND 
*                              ANGLE (J=4) ON WHICH THE I'TH NONLINEAR 
*                              CONSTRAINT IS IMPOSED.  
*                              NSCONS(4,I)=0: CONSTRAINT IMPOSED ON AN
*                                             INTERATOMIC DISTANCE
*     XYZC(L,K,J,I): COEFFICIENTS TO CALCULATE FRACTIONAL COORDINATES 
*                    FROM THOSE OF ATOMS IN THE ASYMMETRIC UNIT
*                    I: CONSTRAINT NUMBER
*                    J: ATOM PAIR/TRIO NUMBER
*                    K: X', Y', AND Z' COORDINATES
*                    L: FOUR COEFFICIENTS TO CALCULATE FRACTIONAL
*                       COORDINATES
*     NPHCON: PHASE NUMBER (READ IN FROM *.ffe)
*     IXCON: PARAMETER NUMBERS FOR ATOMS IN THE ASYMMETRIC UNIT

*     MAXNL: MAXIMUM NUMBER OF NONLINEAR CONSTRAINTS INPUT BY THE USER
      PARAMETER (NB=7000,NPH=8,NAP=150,NS=48,MAXNL=1000)
      INTEGER NSELEC,NSCONS(4,MAXNL),NPHCON,IXCON(*),R,RX,RY,RZ
      REAL CONVFC(4,3,0:1000),XYZC(4,3,3,MAXNL)
      CHARACTER FILE10*50,LINE*136
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /RT/ T(3,NS,NPH),HT(NS,NB),NSYM(NPH),NSITE(NPH),
     &  IDSYM(NS,NAP,NPH),R(3,3,NS,NPH),RX(NS,NB),RY(NS,NB),RZ(NS,NB)
*     S=0: X'=X; Y'=Y; Z'=Z (NO TRANSFORMATION)
      DATA ((CONVFC(K,J,0),K=1,4),J=1,3)/1.0,0.0,0.0,0.0, 
     &  0.0,1.0,0.0,0.0, 0.0,0.0,1.0,0.0/

*     INPUT FILE: FILE #10 (FILE OUTPUT BY ORFFE: *.ffe)
      CALL OPENFILE(8,FILE10)
      REWIND 10
      READ(10,'(A)') LINE
*Rev 1.0u 2001.02.26 Izumi {
C     DO WHILE (LINE(1:30) .NE. '0ATOMIC PARAMETERS FOR PHASE #')
      DO WHILE (LINE(1:30) .NE. ' ATOMIC PARAMETERS FOR PHASE #')
* }
         READ(10,'(A)') LINE
      END DO
*     READ THE PHASE NUMBER
      READ(LINE,'(30X,I1)') NPHCON

*     READ THE NUMBER OF ATOMS IN THE ASYMMETRIC UNIT
      READ(10,'(A)') LINE
*Rev 1.0u 2001.02.26 Izumi {
C     DO WHILE (LINE(1:21) .NE. '0SYMMETRY INFORMATION')
      DO WHILE (LINE(1:21) .NE. ' SYMMETRY INFORMATION')
* }
         READ(10,'(A)') LINE
      END DO
*     SKIP TWO LINE
      READ(10,'()') 
      READ(10,'()') 

*     READ SYMMETRY OPERATIONS
      READ(10,'(A)') LINE
      ISYM=1
      DO WHILE (LINE(1:76) .NE. ' ')
         IF (ISYM .GT. 1000) CALL JOBEND('Array CONVFC overflowed')
         CALL L2CON(LINE,CONVFC(1,1,ISYM))
         ISYM=ISYM+1
         READ(10,'(A)') LINE
      END DO

      READ(10,'(A)') LINE
*Rev 1.0u 2001.02.26 Izumi {
C     DO WHILE (LINE(1:34) .NE. '0INTERATOMIC DISTANCE IN ANGSTROMS')
      DO WHILE (LINE(1:34) .NE. ' INTERATOMIC DISTANCE IN ANGSTROMS')
* }
         READ(10,'(A)',END=7) LINE
      END DO

*     DETERMINE ATOM NUMBERS AND COEFFICIENTS FOR INTERATOMIC DISTANCES
      DO WHILE (LINE .NE. ' ')
         READ(10,'(A)',END=7) LINE
*Rev 1.0u 2001.02.26 Izumi {
C        IF (LINE(1:2) .NE. '0-') THEN
         IF (LINE(1:2) .NE. ' -') THEN
* }
            CALL ATMNOC(NSELEC,CONVFC,LINE,NSCONS,XYZC,2)
            READ(10,'()',END=7)
         END IF
      END DO

      READ(10,'(A)',END=7) LINE
*Rev 1.0u 2001.02.26 Izumi {
C     DO WHILE (LINE(1:23) .NE. '0BOND ANGLE IN DEGREES.')
      DO WHILE (LINE(1:23) .NE. ' BOND ANGLE IN DEGREES.')
* }
         READ(10,'(A)') LINE
      END DO

*     DETERMINE ATOM NUMBERS AND COEFFICIENTS FOR BOND ANGLES
      DO WHILE (LINE .NE. ' ')
         READ(10,'(A)',END=7) LINE
         CALL ATMNOC(NSELEC,CONVFC,LINE,NSCONS,XYZC,3)
         READ(10,'()',END=7)
      END DO
    7 CLOSE(UNIT=10)

*     PARAMETER NUMBERS OF X COORDINATES OF ATOMS FOR NPHCON'TH PHASE
      DO J=1,NSITE(NPHCON)
         IF (J .EQ. 1) THEN
            IX=KPHB(NPHCON)+NG1X+1
         ELSE
            IX=IX+NPSITE(J-1,NPHCON)
         END IF
*        IXCON: PARAMETER NUMBERS OF X COORDINATES
         IXCON(J)=IX
      END DO
      RETURN
      END

************************************************************************

      SUBROUTINE L2CON(LINE,CONVFC)
*     CONVERT A LINE INTO FOUR COEFFICIENTS TO CALCULATE ATOM POSITIONS
      CHARACTER LINE*136,XYZ(2,3)*2
      REAL CONVFC(4,3)
      DATA ((XYZ(J,K),J=1,2),K=1,3)/'+X','-X','+Y','-Y','+Z','-Z'/

*     CONVFC: COEFFICIENTS USED TO CONVERT FRACTIONAL COORDINATES INTO
*             THOSE OF OTHER EQUIVALENT POSITIONS
*     X'=CONVFC(1,1)*X+CONVFC(2,1)*Y+CONVFC(3,1)*Z+CONVFC(4,1)
*     Y'=CONVFC(1,2)*X+CONVFC(2,2)*Y+CONVFC(3,2)*Z+CONVFC(4,2)
*     Z'=CONVFC(1,3)*X+CONVFC(2,3)*Y+CONVFC(3,3)*Z+CONVFC(4,3)
*       X, Y, AND Z: FRACTIONAL COORDINATES OF AN ATOM IN THE ASYMMETRIC 
*                    UNIT
*     X', Y, AND Z': FRACTIONAL COORDINATES OF AN ATOM OUTSIDE THE 
*                    ASYMMETRIC UNIT
*     J SCANS X', Y', AND Z'
      J=1
      DO I=17,65,24
         DO K=1,3
            IF (INDEX(LINE(I+8:I+11),XYZ(1,K)) .GT. 0) THEN
               CONVFC(K,J)=1.0
            ELSE IF (INDEX(LINE(I+8:I+11),XYZ(2,K)) .GT. 0) THEN
               CONVFC(K,J)=-1.0
            ELSE
               CONVFC(K,J)=0.0
            END IF
         END DO
         READ(LINE(I:I+7),'(F8.6)') CONVFC(4,J)
         J=J+1
      END DO
      END

************************************************************************

      SUBROUTINE ATMNOC(NSELEC,CONVFC,LINE,NSCONS,XYZC,NATOM)
*     DETERMINE ATOM NUMBERS AND FOUR COEFFICIENTS TO CALCULATE
*     FRACTIONAL COORDINATES
*     NATOM: NUMBER OF ATOMS RELATED TO THE INTERATOMIC DISTANCE 
*            (NATOM = 2) OR THE BOND ANGLE (NATOM = 3)
*
      PARAMETER (MAXNL=1000)
      INTEGER NSELEC,NSCONS(4,MAXNL),NATOM,IP(3)
      REAL CONVFC(4,3,0:1000),XYZC(4,3,3,MAXNL),C2XYZ(3,0:124)
      CHARACTER LINE*136
*     CELL TRANSLATIONS
      DATA ((C2XYZ(I,J),I=1,3),J=0,124)/
     &   0.0, 0.0, 0.0, 
     &  -1.0, 0.0, 0.0, -1.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, 1.0, 1.0, 
     &   0.0,-1.0, 0.0,  1.0,-1.0, 0.0,  0.0,-1.0, 1.0, 
     &   1.0,-1.0, 1.0, -1.0,-1.0, 0.0, -1.0,-1.0, 1.0,  0.0, 0.0,-1.0, 
     &   1.0, 0.0,-1.0,  0.0, 1.0,-1.0,  1.0, 1.0,-1.0, 
     &  -1.0, 0.0,-1.0, -1.0, 1.0,-1.0,  0.0,-1.0,-1.0,  1.0,-1.0,-1.0, 
     &  -1.0,-1.0,-1.0,  1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  1.0, 1.0, 0.0,  
     &   0.0, 0.0, 1.0,  1.0, 0.0, 1.0,  0.0, 1.0, 1.0,  1.0, 1.0, 1.0,
     &  -2.0, 0.0, 0.0,  2.0, 0.0, 0.0,  0.0,-2.0, 0.0,  0.0, 2.0, 0.0,
     &   0.0, 0.0,-2.0,  0.0, 0.0, 2.0,
     &  -2.0,-2.0, 0.0, -2.0,-1.0, 0.0, -2.0, 1.0, 0.0, -2.0, 2.0, 0.0,
     &  -2.0, 0.0,-2.0, -2.0, 0.0,-1.0, -2.0, 0.0, 1.0, -2.0, 0.0, 2.0,
     &  -1.0,-2.0, 0.0, -1.0, 2.0, 0.0, -1.0, 0.0,-2.0, -1.0, 0.0, 2.0,
     &   0.0,-2.0,-2.0,  0.0,-2.0,-1.0,  0.0,-2.0, 1.0,  0.0,-2.0, 2.0,
     &   0.0,-1.0,-2.0,  0.0,-1.0, 2.0,  0.0, 1.0,-2.0,  0.0, 1.0, 2.0,
     &   0.0, 2.0,-2.0,  0.0, 2.0,-1.0,  0.0, 2.0, 1.0,  0.0, 2.0, 2.0,
     &   1.0,-2.0, 0.0,  1.0, 2.0, 0.0,  1.0, 0.0,-2.0,  1.0, 0.0, 2.0,
     &   2.0,-2.0, 0.0,  2.0,-1.0, 0.0,  2.0, 1.0, 0.0,  2.0, 2.0, 0.0,  
     &   2.0, 0.0,-2.0,  2.0, 0.0,-1.0,  2.0, 0.0, 1.0,  2.0, 0.0, 2.0,
     &  -2.0,-2.0,-2.0, -2.0,-2.0,-1.0, -2.0,-2.0, 1.0, -2.0,-2.0, 2.0, 
     &  -2.0,-1.0,-2.0, -2.0,-1.0,-1.0, -2.0,-1.0, 1.0, -2.0,-1.0, 2.0, 
     &  -2.0, 1.0,-2.0, -2.0, 1.0,-1.0, -2.0, 1.0, 1.0, -2.0, 1.0, 2.0, 
     &  -2.0, 2.0,-2.0, -2.0, 2.0,-1.0, -2.0, 2.0, 1.0, -2.0, 2.0, 2.0, 
     &  -1.0,-2.0,-2.0, -1.0,-2.0,-1.0, -1.0,-2.0, 1.0, -1.0,-2.0, 2.0, 
     &  -1.0,-1.0,-2.0, -1.0,-1.0, 2.0, -1.0, 1.0,-2.0, -1.0, 1.0, 2.0, 
     &  -1.0, 2.0,-2.0, -1.0, 2.0,-1.0, -1.0, 2.0, 1.0, -1.0, 2.0, 2.0, 
     &   1.0,-2.0,-2.0,  1.0,-2.0,-1.0,  1.0,-2.0, 1.0,  1.0,-2.0, 2.0, 
     &   1.0,-1.0,-2.0,  1.0,-1.0, 2.0,  1.0, 1.0,-2.0,  1.0, 1.0, 2.0, 
     &   1.0, 2.0,-2.0,  1.0, 2.0,-1.0,  1.0, 2.0, 1.0,  1.0, 2.0, 2.0, 
     &   2.0,-2.0,-2.0,  2.0,-2.0,-1.0,  2.0,-2.0, 1.0,  2.0,-2.0, 2.0, 
     &   2.0,-1.0,-2.0,  2.0,-1.0,-1.0,  2.0,-1.0, 1.0,  2.0,-1.0, 2.0, 
     &   2.0, 1.0,-2.0,  2.0, 1.0,-1.0,  2.0, 1.0, 1.0,  2.0, 1.0, 2.0, 
     &   2.0, 2.0,-2.0,  2.0, 2.0,-1.0,  2.0, 2.0, 1.0,  2.0, 2.0, 2.0
     &   / 
     
      READ(LINE,'(1X,I5)') ISER
*     I SCANS NONLINEAR CONSTRAINTS
      DO I=1,NSELEC
         IF (ISER .EQ. NSCONS(1,I)) THEN
*           NSCONS(4,I)=0: THIS CONSTRAINT IS IMPOSED ON AN INTERATOMIC 
*                          DISTANCE (OVERWRITTEN WHEN ATOM = 3)
            NSCONS(4,I)=0
            READ(LINE,'(12X,3(I2,9X))') (NSCONS(J,I),J=2,NATOM+1)
            READ(LINE,'(12X,3(3X,I5,3X))') (IP(J),J=1,NATOM)
*           J SCANS ATOM PAIR/TRIO NUMBERS
            DO J=1,NATOM
*              IC: C IN THE OUTPUT OF ORFFE
               IC=IP(J)/1000
*              IS: S IN THE OUTPUT OF ORFFE
               IS=IP(J)-IC*1000
*              K SCANS X', Y', AND Z' COORDINATES
               DO K=1,3
*                 L SCANS FOUR COEFFICIENTS
                  DO L=1,3 
                     XYZC(L,K,J,I)=CONVFC(L,K,IS)
                  END DO
*                 TRANSLATE USING IC
                  XYZC(4,K,J,I)=CONVFC(4,K,IS)+C2XYZ(K,IC)
               END DO
            END DO
            RETURN
         END IF
      END DO
      END

************************************************************************

      FUNCTION CON(I,A)
*     CONSTRAINTS IMPOSED ON INTERATOMIC DISTANCES AND BOND ANGLES
*     THIS FUNCTION RETURNS THE I'TH VIOLATED CONSTRAINTS.  EACH
*     CONSTRAINT IS FORMULATED BY USING THE CURRENT VALUES OF
*     PARAMETERS (A).  CELL CONSTANTS SHOULD NOT BE CONTAINED IN THESE
*     CONSTRAINTS.

*     FUNCTIONS DIS AND ANG ARE USED TO OBTAIN THE VALUES OF INTERATOMIC
*     DISTANCES AND BOND ANGLES, RESPECTIVELY.

*     IN A CONSTRAINT IMPOSED ON AN INTERATOMIC DIATANCE, THE DISTANCE
*     ARE EXPECTED TO BE EXPCTD(I) +/- DEVDA(I).  IF THE DISTANCE IS
*     WITHIN THIS RANGE, CON IS SET EQUAL TO ZERO.  A CONSTRAINT FOR A
*     BOND ANGLE CAN BE FORMULATED IN A SIMILAR WAY.

      PARAMETER (MAXNL=1000,NAP=150)
      REAL A(*),XYZ(3,3)
      COMMON /CONSDA/ NPHCON,NSCONS(4,MAXNL),XYZC(4,3,3,MAXNL),
     &  EXPCTD(MAXNL),DEVDA(MAXNL),IXCON(NAP)

      IF (NSCONS(4,I) .EQ. 0) THEN
*        A CONSTRAINT IS IMPOSED ON AN INTERATOMIC DISTANCE 
         NATOM = 2
      ELSE 
*        A CONSTRAINT IS IMPOSED ON A BOND ANGLE
         NATOM = 3
      END IF
      DO J=1,NATOM
         NX=IXCON(NSCONS(J+1,I))
*        K SCANS X', Y', AND 'Z COORDINATES
         DO K=1,3
            XYZ(K,J) = XYZC(1,K,J,I)*A(NX) + XYZC(2,K,J,I)*A(NX+1) +
     &                 XYZC(3,K,J,I)*A(NX+2) + XYZC(4,K,J,I)
         END DO
      END DO
      IF (NATOM .EQ. 2) THEN
         DIS12 = DIS(XYZ(1,1),XYZ(2,1),XYZ(3,1),
     &           XYZ(1,2),XYZ(2,2),XYZ(3,2),NPHCON)
         CON = MIN(0.0, DEVDA(I) - ABS(DIS12 - EXPCTD(I))) / 
     &         SQRT(DIS12)
      ELSE
         ANG123 = ANG(XYZ(1,1),XYZ(2,1),XYZ(3,1),XYZ(1,2),XYZ(2,2),
     &            XYZ(3,2),XYZ(1,3),XYZ(2,3),XYZ(3,3),NPHCON)
         CON = MIN(0.0, DEVDA(I) - ABS(ANG123 - EXPCTD(I))) / 
     &         SQRT(ANG123)
      END IF
      END

************************************************************************

      FUNCTION DIS(X1,Y1,Z1,X2,Y2,Z2,IP)
C     DIS: DISTANCE BETWEEN TWO ATOMS WHICH ARE LOCATED AT (X1,Y1,Z1)
C          AND (X2,Y2,Z2)
C     IP: PHASE NUMBER
C     THIS FUNCTION IS USED TO GIVE CONSTRAINTS IN SUBROUTINE CON
      PARAMETER (NPH=8)
      COMMON /TENS/ G(6,NPH),GG(6,NPH)

      DIS=G(1,IP)*(X1-X2)**2+G(2,IP)*(Y1-Y2)**2+G(3,IP)*(Z1-Z2)**2+
     &  2.0*(G(4,IP)*(Y1-Y2)*(Z1-Z2)+G(5,IP)*(Z1-Z2)*(X1-X2)+
     &  G(6,IP)*(X1-X2)*(Y1-Y2))
      DIS=SQRT(DIS)
      END

************************************************************************

      FUNCTION ANG(X2,Y2,Z2,X1,Y1,Z1,X3,Y3,Z3,IP)
C     ANG: BOND ANGLE (X2,Y2,Z2) - (X1,Y1,Z1) - (X3,Y3,Z3)
C     IP: PHASE NUMBER

      AA=DIS(X2,Y2,Z2,X3,Y3,Z3,IP)
      CC=DIS(X1,Y1,Z1,X2,Y2,Z2,IP)
      BB=DIS(X1,Y1,Z1,X3,Y3,Z3,IP)
      COSAA=(BB*BB+CC*CC-AA*AA)/(2.0*BB*CC)
      ANG=57.29578*ACOS(COSAA)
      END
