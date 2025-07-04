      SUBROUTINE CURFIT(X,Y,NSTEP,NPTS,NTERMS,A,AMIN,FLAMDA,FLAMDC,
     &  FLAMDP,DAMP,IFAIL,YFIT,CHISQ1,CHISQR,CSMIN,ICYCLE,LDM,NFNEV)
C     NONLINEAR LEAST-SQUARES FITTING BY MARQUARDT'S METHOD AND
C     GAUSS-NEWTON METHOD
      PARAMETER (NB=7000,NR=400,NSF=400,NT=999,NPH=8,NOD=NR*(NR+1)/2)
      PARAMETER (NP=80000,NAP=150)
      DOUBLE PRECISION G(NOD),ALPHA(NOD),BETA(NR),DELTAB(NR),AD(NR)
      SAVE ALPHA,AD
      DIMENSION X(*),Y(*),A(*),AMIN(*),YFIT(*),
     &  B(NT),D(NT),DERIV(NR),DSTFAC(NSF,NB),ISFP(NT)
      INTEGER JSAVE(NR)
      COMMON /ALP/ LALPHA
      COMMON /DETECT/ VARINV(NP)
      COMMON /FL0/ NC,TK,CONV,NCONV,FINC,SUMPEN,SUMP1,FACTOR
      COMMON /NONLIN/ NLESQ,NAUTO,IDCYCL,NCHNG,NPAR(50),IPAR(NR,50)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /TWO/ G,FIRST(NR)
      EXTERNAL PACK
C
      IF (IDCYCL .EQ. 0 .AND. LALPHA .EQ. 1) THEN
         DO J = 1, NRFN
            ALPHA(J*(J+1)/2)=AD(J)
         END DO
         GO TO 3
      ENDIF
C     NFNEV: NUMBER OF FUNCTION EVALUATIONS IN THIS CYCLE
      NFNEV=0
      CALL DERSP(A,DSTFAC,ISFP)
C
C     ALPHA: LOWER TRIANGLE OF THE COEFFICIENT MATRIX
C     BETA: RIGHT-HAND SIDE OF THE NORMAL EQUATION
C     1/Y(I):  WEIGHT
      DO J = 1, NRFN
         BETA(J)=0.0D0
         JSAVE(J)=(J-1)*J/2
      END DO
      DO J = 1, NRFN*(NRFN+1)/2
         ALPHA(J)=0.0D0
      END DO
C     NREST: NUMBER FOR THE FIRST REFINABLE PARAMETER WHICH FOLLOWS
C            BACKGROUND PARAMETERS
      NREST=1
      DO J = 1, NBCKGR + NUMBG -1
         IF (ID(J) .EQ. 1) NREST=NREST+1
      END DO
   
      IF (VARINV(1) .EQ. 0.0) THEN
         DO I = 1, NPTS
            CALL FDERIV(I,A,X,DERIV,DSTFAC,ISFP,NREST)
            DELTAY = Y(I) - YFIT(I)
            DO J = 1, NRFN
               BETA(J) = BETA(J) + DERIV(J)/Y(I)*DELTAY
               DO K = 1, J
                  ALPHA(K+JSAVE(J)) = ALPHA(K+JSAVE(J)) +
     &            DERIV(J)/Y(I)*DERIV(K)
               END DO
            END DO
         END DO
      ELSE
         DO I = 1, NPTS
            CALL FDERIV(I,A,X,DERIV,DSTFAC,ISFP,NREST)
            DELTAY = Y(I) - YFIT(I)
            DO J = 1, NRFN
               BETA(J) = BETA(J) + DERIV(J)*VARINV(I)*DELTAY
               DO K = 1, J
                  ALPHA(K+JSAVE(J)) = ALPHA(K+JSAVE(J)) +
     &            DERIV(J)*VARINV(I)*DERIV(K)
               END DO
            END DO
         END DO
      END IF
 
C     INVERT ARRAY ALPHA TO CALCULATE STANDARD DEVIATIONS OF PARAMETERS
    3 IF (IDCYCL .EQ. 0) THEN
C        INVERSION OF THE POSITIVE DEFINITE MATRIX
         CALL CHLSK(NRFN,ALPHA,G,0.0D0)
         CALL MATINV(G,NRFN)
         RETURN
      END IF
C     CONTRIBUTION OF THE PENALTY FUNCTION TO ARRAYS ALPHA AND BETA
      IF (NC .GT. 0) CALL PENFN(ALPHA,BETA,A,NTERMS)
C
C     PACK ARRAY A
      CALL PACK(A,D,NTERMS)
C     STORE COEFFICIENTS OF FLAMDA, THAT IS, DIAGONAL TERMS OF THE
C     COEFFICIENT MATRIX
      DO J = 1, NRFN
         AD(J)=ALPHA(J*(J+1)/2)
      END DO
C
    4 IF (NLESQ .EQ. 0) THEN
C        GRADIENT-EXPANSION ALGORITHM.
         DO J = 1, NRFN
            ALPHA(J*(J+1)/2)=AD(J)*(1.0+FLAMDA)
         END DO
      ENDIF
C     CHOLESKY DECOMPOSITION
      CALL CHLSK(NRFN,ALPHA,G,0.0D0)
C     SOLVE SIMULTANEOUS LINEAR EQUATIONS
      CALL SOLVE(NRFN,G,DELTAB,BETA)
C
      IF (NLESQ .EQ. 0) THEN
         CALL RSSE(A,B,DELTAB,D,X,Y,YFIT,AMIN,NPTS,NSTEP,
     &   NTERMS,CSMIN,CHISQR,1.0,IDPOS,IFAIL,NFNEV,FLAMDA,FLAMDP)
         IF (NAUTO .EQ. 0 .OR. CHISQ1 .LT. CHISQR .OR. (NAUTO .GE. 1 
     &   .AND. ICYCLE .GE. NCHNG+1)) CALL LAMBDA(DELTAB,BETA,AD,G,
     &   FLAMDA,FLAMDC,IDPOS,LDM,CHISQ1,CHISQR)
         IF(LDM.EQ.1) THEN
            WRITE(6,'(//11X,A,I2/11X,A)') 'Cycle #',
     &      ICYCLE,'Too large Marquardt parameter'
         ELSEIF(CHISQ1.LT.CHISQR) THEN
            GO TO 4
         ENDIF
      ELSE
C        DAMP: VARIABLE DAMPING (RELAXATION) FACTOR
         DAMP=MIN(1.0,2.0*DAMP)
    6    CALL RSSE(A,B,DELTAB,D,X,Y,YFIT,AMIN,NPTS,NSTEP,
     &   NTERMS,CSMIN,CHISQR,DAMP,IDPOS,IFAIL,NFNEV,FLAMDA,FLAMDP)
         IF(IDPOS.NE.1 .OR. CHISQ1.LT.CHISQR) THEN
            IF(DAMP.LE.2.0**(-5)) THEN
               WRITE(6,150) ICYCLE
               IFAIL=IFAIL+1
               LDM=1
               IF(IFAIL.GE.2) WRITE(6,200) ICYCLE,IFAIL
            ELSE
C              REDUCE THE DAMPING FACTOR AND CREATE A NEW SET OF
C              PARAMETERS
               DAMP=DAMP*0.5
               GO TO 6
            ENDIF
         ENDIF
      ENDIF
  150 FORMAT(//11X,'Cycle #',I2/11X,'Failure.  The objective function do
     &es not decrease in this cycle.')
  200 FORMAT(//11X,'Cccle #',I2/11X,'The number of failures has reached
     &the maximum value (',I1,').')
C
C     COPY THE NEW VALUES OF PARAMETERS
      DO 90 J=1,NTERMS
         A(J)=B(J)
   90 CONTINUE
      END

************************************************************************

      FUNCTION AOFFN(X,Y,YFIT,NPTS,NSTEP,A,IDHKL,IDPOS)
C     UPDATE INDICES AND FWHM AND CALCULATE THE VALUE OF THE AUGMENTED
C     OBJECTIVE FUNCTION USING NEW PARAMETERS
      PARAMETER (NP=80000,NPH=8,NAP=150,NB=7000,NR=400)
      DOUBLE PRECISION CHISQ
      DIMENSION X(*),Y(*),YFIT(*),A(*)
      COMMON /ORDER/ IPOINT(NB)
      COMMON /B/ U(NB),NREFL
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /ABSVR/ ABSORP(NP),VARLEN(NP),ROUGH(NP),NSURFR
      COMMON /C/ BG(NP),BGINC(NP),DEGNOR(NP)
      COMMON /DETECT/ VARINV(NP)
      COMMON /FL0/ NC,TK,CONV,NCONV,FINC,SUMPEN,SUMP1,FACTOR
      COMMON /NONLIN/ NLESQ,NAUTO,IDCYCL,NCHNG,NPAR(50),IPAR(NR,50)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /RLS/ NRANGE
      PARAMETER (NT=999)
      COMMON /SAVEYP/ SYPEAK(NB),SAVEA(NT),NTOTPAR
      COMMON /Z/ PC,CHGPC
      SAVE CSMIN
      DATA CSMIN/1.0E20/
	  
*     RESTORE YPEAK FOR PARAMETERS GIVING THE MINIMUM SUM OF SQUARES
*Rev 1.2 2003.05.28 {
      IF ((NMODE .EQ. 4 .AND. NCONST .EQ. 0) .OR. NMODE .EQ. 5) 
     &  CALL SCYPEAK(1)
* }
      IDPOS=1
      IF(NRANGE.EQ.1) THEN
         REWIND 40
         IF (VARINV(1) .EQ. 0.0) THEN
            READ(40) (X(J),Y(J),BG(J),DEGNOR(J),ABSORP(J),VARLEN(J),
     &      ROUGH(J), J=1,NSTEP)
         ELSE
            READ(40) (X(J),Y(J),VARINV(J),BG(J),DEGNOR(J),ABSORP(J),
     &      VARLEN(J),ROUGH(J),J=1,NSTEP)
         END IF
      ENDIF
      IF(IDHKL.EQ.1) CALL UPDATE(A,IDPOS,1)

      IF (IDPOS .NE. 1) THEN
*        This may cause errors when using Digital Fortran
*        AOFFN=1.7E38
         AOFFN=1.0E20
         RETURN
      ENDIF
      CALL IGNORE(NSTEP,NPTS,X,Y)

*     UPDATE YPEAK, i.e., INTEGRATED INTENSITIES
*     INTEGINT IS CALLED EXCEPT FOR THE FIRST TIME
*Rev 1.2 2003.05.28 {
      IF (((NMODE .EQ. 4 .AND. NCONST .EQ. 0) .OR. NMODE .EQ. 5) 
     &  .AND. IDHKL .EQ. 1) CALL INTEGINT(A,X,YFIT,Y,NPTS)
* }
      CHISQ=0.0D0
      IF (VARINV(1) .EQ. 0.0) THEN
         DO J=1,NPTS
            YFIT(J)=CALINT(X,J,A)
            CHISQ=CHISQ+(Y(J)-YFIT(J))**2/Y(J)
         END DO
      ELSE
         DO J=1,NPTS
            YFIT(J)=CALINT(X,J,A)
            CHISQ=CHISQ+(Y(J)-YFIT(J))**2*VARINV(J)
         END DO
      END IF
C
C     PENALTY FUNCTION
      SUMPEN=0.0
      DO J = 1, NC
         SUMPEN=SUMPEN+TK*(CON(J,A))**2
      END DO
      AOFFN=CHISQ+SUMPEN

*Rev 1.2 2003.05.28 {
      IF (((NMODE .EQ. 4 .AND. NCONST .EQ. 0) .OR. NMODE .EQ. 5) 
     &  .AND. AOFFN .LT. CSMIN) THEN
* }
*        SAVE THE CURRENT PARAMETER AND YPEAK GIVING THE MINIMUM 
*        SUM OF SQUARES
         CSMIN = AOFFN
*        SAVE CURRENT A
         DO J = 1, NTOTPAR
            SAVEA(J) = A(J)
         END DO
*        SAVE CURRENT YPEAK
         DO J = 1, NREFL
            I = IPOINT(J)
            SYPEAK(J) = YPEAK(I)
         END DO
      END IF
      END

************************************************************************

      SUBROUTINE SCYPEAK(LSC)
*     SAVE OR COPY ARRAY YPEAK
      PARAMETER (NB=7000,NPH=8)
      COMMON /B/ U(NB),NREFL
      COMMON /ORDER/ IPOINT(NB)
      PARAMETER (NT=999)
      COMMON /SAVEYP/ SYPEAK(NB),SAVEA(NT),NTOTPAR
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)

      SELECT CASE (LSC)
         CASE (0)
*           SAVE YPEAK IN SYPEAK
            DO J = 1, NREFL
               I = IPOINT(J)
               SYPEAK(J) = YPEAK(I)
            END DO
         CASE (1)
*           COPY SYPEAK INTO YPEAK
            DO J = 1, NREFL
               I = IPOINT(J)
               YPEAK(J) = SYPEAK(I)
            END DO
      END SELECT
      END

************************************************************************

      SUBROUTINE FDERIV(K,A,DEG,DERIV,DSTFAC,ISFP,NREST)
*     PARTIAL DERIVATIVES THE MODEL FUNCTION WITH RESPECT TO REFINABLE
*     PARAMETERS
      PARAMETER (NP=80000,NB=7000,NT=999,NR=400,NSF=400,NPH=8,NAP=150)
      INTEGER H
      REAL A(*)
      DIMENSION DEG(*),DERIV(*),DSTFAC(NSF,NB),ISFP(*)
      COMMON /A/ H(NB),IK(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /ABSVR/ ABSORP(NP),VARLEN(NP),ROUGH(NP),NSURFR
      COMMON /B/ U(NB),NREFL
      COMMON /C/ BG(NP),BGINC(NP),DEGNOR(NP)
      COMMON /DERPR/ POP2(NB),COMPO(NB),CS(NB),CCC(6,NB)
      COMMON /DPRFL/ DPRDT,SIGPART,GAMPART,DPRDC7,DPRDS,DPRDD
      COMMON /DPRSPV/ DFDX,DFDA,DFDW,DFDW1,DFDW2,DFDEL,DFDEH,DFDRL,DFDRH
      COMMON /G/ I
      COMMON /I/ RDEG(2,NB)
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /LEGEN/ MAXDEG
      COMMON /MC/ CTHM,NTRAN,PCOR
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /ORDER/ IPOINT(NB)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /RLS/ NRANGE
      COMMON /UCC/ IDPEAK(NPH),IDSF(NPH),IDPRF(NPH),IDPO(NPH)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      COMMON /UNCONS/ NUCREF,IPHASE(90),IHKL(3,90),NOPPP(90),LPPP(NB),
     &  RCUT(2,90),NCUT
      SAVE ISTART
      LP(J)=KPHB(NOPH(I))+J
      LD(J)=ID(LP(J))
      LR(J)=IR(LP(J))

C     RESET ARRAY DERIV
      DO J = 1, NRFN
         DERIV(J)=0.0
      END DO
      
*     NUMERICAL PARTIAL DERIVATIVES OF BACKGROUND PARAMETERS
*     LEGENDRE POLYNOMIALS ARE USED.
*     WE ASSUME THAT b0 IS ALWAYS REFINED WHEN REFINING BACKGROUND
*     PARAMETERS AND THAT NO BACKGROUND PARAMETERS ARE SKIPPED UP
*     TO THE TERM WITH THE MAXDEG DEGREE.
      IF (ID(NBCKGR) .EQ. 1) THEN
         CALL FLEG(DEGNOR(K), DERIV(IR(NBCKGR)), MAXDEG+1)
         DO J = IR(NBCKGR), IR(NBCKGR)+MAXDEG
            SELECT CASE(NRANGE)
               CASE(0)
                  DERIV(J) = VARLEN(K)*DERIV(J)
               CASE(3)
                  DERIV(J) = VARLEN(K)*DERIV(J)*BGINC(K)
            END SELECT
         END DO
      END IF
      IF (K .EQ. 1) ISTART=1
      CALINT = 0.0
      DO J0 = ISTART, NREFL
         I=IPOINT(J0)
         IF (DEG(K) .GT. RDEG(2,I)) THEN
            ISTART=J0+1
            CYCLE
         ELSE IF (DEG(K) .LT. RDEG(1,I)) THEN
            EXIT
         ENDIF
         
C        PROFILE FUNCTION OF THE I'TH REFLECTION AT DEG(K)
         PRFIK = PRFL(DEG(K),A)
         CALINT = CALINT + YPEAK(I)*PRFIK
C
C        SCALE FACTOR
         IF(LD(NSFX).EQ.1) DERIV(LR(NSFX))=DERIV(LR(NSFX))+PRFIK*CS(I)
C
C        PARAMETERS CONTAINED IN THE PROFILE-SHAPE FUNCTION
         IF (IDPRF(NOPH(I)) .EQ. 1) CALL DERPF(A,DERIV,PRFIK)
C
C        PREFERRED-ORIENTATION PARAMETERS
         IF (LD(NPO1X) .EQ. 1) 
     &   DERIV(LR(NPO1X)) = DERIV(LR(NPO1X)) + PRFIK*COMPO(I)
         IF (LD(NPO2X) .EQ. 1)
     &   DERIV(LR(NPO2X)) = DERIV(LR(NPO2X)) + PRFIK*POP2(I)
C
C        ELEMENTS OF THE METRIC TENSOR
         IF (IDPEAK(NOPH(I)) .EQ. 1) THEN
            JJ = 1
*           DPRDT, DFDX: DERIVATIVE OF THE PROFILE FUNCTION WITH
*           RESPECT TO DELTA_2-THETA
            SELECT CASE (NPRFN)
               CASE (0)
                  DTEMP = DPRDT
               CASE (1:3)
                  DTEMP = DFDX
            END SELECT
            DO J = NAX, NGAMX
               IF (LD(J) .EQ. 1) 
     &         DERIV(LR(J)) = DERIV(LR(J)) + DTEMP*CCC(JJ,I)
               JJ = JJ + 1
            END DO
         END IF
C
         IF (IDSF(NOPH(I)) .EQ. 1 .AND. NMODE .EQ. 0) THEN
C           OVERALL ISOTROPIC DISPLACEMENT PARAMETER AND STRUCTURE 
C           PARAMETERS
            DO J = LP(NOTX), KPHE(NOPH(I))
               IF (ID(J) .EQ. 1) DERIV(IR(J))=DERIV(IR(J))+PRFIK*
     &         DSTFAC(ISFP(J),I)
            END DO
         END IF
      END DO

*     DERIVATIVES WITH RESPECT TO PEAK-SHIFT PARAMETERS
      DO J = 1, 4
         DERIV(IR(J)) = ABSORP(K)*VARLEN(K)*ROUGH(K)*DERIV(IR(J))
      END DO
      IF (NSURFR .GT. 0) 
     &  CALL DERSR(0.5*DEG(K),A,DERIV,ABSORP(K)*VARLEN(K)*CALINT)

*     DERIVATIVES WITH RESPECT TO BACKGROUND PARAMETERS ARE SKIPPED
      DO J = NREST, NRFN
         DERIV(J) = ABSORP(K)*VARLEN(K)*ROUGH(K)*DERIV(J)
      END DO
      END

************************************************************************

      SUBROUTINE DERPF(A,DERIV,PRFIK)
*     CALCULATE DERIVATIVES WITH RESPECT TO PARAMETERS CONTAINED IN
*     THE PROFILE-SHAPE FUNCTION
      PARAMETER (NT=999,NB=7000,NPH=8,NAP=150)
      INTEGER H,HANIS
      REAL A(*),DERIV(*)
      LOGICAL LCONS,LINTEG,LPEAKPOS
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /ANISBR/ HANIS(NPH),KANIS(NPH),LANIS(NPH),
     &  CANIS(NPH),PCOS(NB)
      COMMON /DPRFL/ DPRDT,SIGPART,GAMPART,DPRDC7,DPRDS,DPRDD
      COMMON /DPRSPV/ DFDX,DFDA,DFDW,DFDW1,DFDW2,DFDEL,DFDEH,DFDRL,DFDRH
      COMMON /FT/ FT1,FT2,PEAKNOR(NB),TANTHNOR(NB) 
      COMMON /G/ I
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /MC/ CTHM,NTRAN,PCOR
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /PROFI/ PCOT(NB),SHFT(NB),NTSIM(NB),FNORM(NB),TNTSIM(NB),
     &  FWHM(NB),ETA(NB),SQSG(NB),CTH(NB),SIGP(NB),DFWDG(NB),DEDFFG(NB),
     &  DFWDL(NB),DEDFFL(NB),SIGPAP(NB)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /T/ FF(NB),COEF(NPH)
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      COMMON /ZSPFN/ PEAKSH(NB),SQH2(NB),SQH2L(NB),SQH2G(NB),ASYM(NB),
     &ETAL(NB),ETAH(NB),FNORM1(2,NB),COMLG(2,NB),COML(2,NB),COMG(2,NB),
     &COEFL(2,NB),COEFG(2,NB),DECAY(2,NB),COMPEAR(2,NB),CDX(2,NB),
     &CDA1(2,NB),CDA2(2,NB),CDW(2,NB),CDRL1(2,NB),CDRL2(2,NB),
     &CDRH1(2,NB),DMPART(2,NB)
      COMMON /UNCONS/ NUCREF,IPHASE(90),IHKL(3,90),NOPPP(90),LPPP(NB),
     &  RCUT(2,90),NCUT

*     NOPPP(LPPP(I)) + J: SERIAL PARAMETER NUMBER FOR THE J'TH 
*     PRIMARY PROFILE PARAMETER (J = 0 - (3 OR 4)) FOR REFLECTION I
*     SPECIFIED BY THE USER
      LD2(J) = ID(NOPPP(LPPP(I)) + J)
      LR2(J) = IR(NOPPP(LPPP(I)) + J)

*     LCONS(J) = .TRUE.: The (J+1)-th primary profile parameter is
*     conventionally calculated from secondary profile parameters.
*     Reflection I is either constrained (LPPP(I) = 0) or relaxed 
*     (LPPP(I) > 0).  In the latter case, LD2(J) = 3.

*     LCONS(J) = .FALSE.: The (J+1)-th primary profile parameter is
*     input and refined/fixed.

      LCONS(J) = LPPP(I) .EQ. 0 .OR. 
     &           (LPPP(I) .GT. 0 .AND. LD2(J) .EQ. 3)
     
*     LINTEG(J) = .TRUE.: The (J+1)-th primary profile parameter, which
*     is |Fc|, for relaxed reflection I is refined in MEM-based whole-pattern 
*     fitting and individual profile fitting.

      LINTEG(J) = (NMODE .EQ. 3 .OR. NMODE .EQ. 6) .AND. LPPP(I) .GT. 0 
     &            .AND. LD2(J) .EQ. 1

*     LPEAKPOS(J) = .TRUE.: THE (J+1)-th primary profile parameter, which
*     is 2-theta(peak), for relaxed reflection I is refined in individual
*     profile fitting.

      LPEAKPOS(J) = NMODE .EQ. 6 .AND. LPPP(I) .GT. 0 .AND. 
     &              LD2(J) .EQ. 1

      TTH = TAN(TH(I))
      SELECT CASE (NPRFN)
         CASE (0)
*           Pseudo-Voigt function of Thompson, Cox, and Hastings
*           Profile relaxtion cannot be applied to this profile function
            CSTH = COS(TH(I))
            YPD = YPEAK(I)*DPRDT
         
*           (1) DERIVATIVES WITH RESPECT TO THE THREE PEAK-SHIFT 
*               PARAMETERS
*           ZERO-POINT SHIFT, Z
            IF (ID(1) .EQ. 1) DERIV(1) = DERIV(1) + YPD
*           SPECIMEN-DISPLACEMENT PARAMETER, Ds
*Rev 1.07 2002.08.22 Izumi {
            IF (ID(2) .EQ. 1 .AND. NTRAN .EQ. 2) THEN
*              TRANSMISSION GEOMETRY
               DERIV(IR(2)) = DERIV(IR(2)) + YPD*SIN(TH(I))
            ELSE IF (ID(2) .EQ. 1) THEN                              
               DERIV(IR(2)) = DERIV(IR(2)) + YPD*CSTH
            END IF
C           IF (ID(2) .EQ. 1 .AND. NTRAN .EQ. 0) THEN
C              DERIV(IR(2)) = DERIV(IR(2)) + YPD*CSTH
C           ELSE IF (ID(2) .EQ. 1 .AND. NTRAN .EQ. 1) THEN                              
C              DERIV(IR(2)) = DERIV(IR(2)) + YPD*SIN(TH(I))
C           END IF
* }
*           SPECIMEN-TRANSPARENCY PARAMETER, Ts
            IF (ID(3) .EQ. 1) DERIV(IR(3)) = DERIV(IR(3)) + 
     &                                       YPD*SIN(PEAK(I))
     
*           (2) DERIVATIVES WITH RESPECT TO THE PROFILE PARAMETERS
*           GAUSSIAN FWHM PARAMETER, U
            CALL DERCOM(DERIV, 1, SIGPART*TTH*TTH)
*           GAUSSIAN FWHM PARAMETER, V
            CALL DERCOM(DERIV, 2, SIGPART*TTH)
*           GAUSSIAN FWHM PARAMETER, W
            CALL DERCOM(DERIV, 3, SIGPART)
*           SCHERRER COEFFICIENT FOR GAUSSIAN BROADENING, P
            CALL DERCOM(DERIV, 4, SIGPART/(CSTH*CSTH))
*           LORENTZIAN SCHERRER BROADENING, X
            CALL DERCOM(DERIV, 5, GAMPART/CSTH)
*           ANISOTROPY COEFFICIENT FOR SCHERRER BROADENING, Xe
            CALL DERCOM(DERIV, 6, GAMPART*PCOS(I)/CSTH)
*           STRAIN BROADENING, Y
            CALL DERCOM(DERIV, 7, GAMPART*TTH)
*           ANISOTROPY COEFFICIENT FOR STRAIN BROADENING, Ye
            CALL DERCOM(DERIV, 8, GAMPART*PCOS(I)*TTH)
            IF (NASYM .EQ. 0) THEN
*              r.s
               CALL DERCOM(DERIV, 9, DPRDS)
*              r.d      
               CALL DERCOM(DERIV, 10, DPRDD)
            ELSE
*              ASYMMETRY PARAMETER, As
               CALL DERCOM(DERIV, 9, FNORM(I)*DPRDC7*CTH(I))
            END IF
         CASE (1,2)
*           Split-type pseudo-Voigt function
            CALL DERSHIFT(TTH,DERIV)
      
            IF (NPRFN .EQ. 1 .OR. (NPRFN .EQ. 2 .AND. LPPP(I) .EQ. 0)) 
     &      THEN

*              Derivatives of the model function, f, with respect to 
*              secondary profile parameters, y, are obtained with
*              df/dy = YPEAK(I)*(dG/dx)*(dx/dy), 
*              where G is the profile function, and x is a primary 
*              profile parameter.

*              Derivatives of the model function with respect to primary 
*              profile parameters for a relaxed reflection are given by 
*              df/dx = YPEAK(I)*(dG/dx).
               CALL DERWA(TTH,DERIV)

*              ETAL
               IF (LCONS(2)) THEN
*                 eta0 FOR THE LORENTZIAN FRACTION (LOWER ANGLE)
                  CALL DERCOM(DERIV, 9, DFDEL)
*                 eta1 FOR THE LORENTZIAN FRACTION (LOWER ANGLE)
                  CALL DERCOM(DERIV, 10, DFDEL*PEAK(I))
               ELSE IF (LD2(2) .EQ. 1) THEN
                  DERIV(LR2(2)) = DERIV(LR2(2)) + YPEAK(I)*DFDEL
               END IF
         
*              ETAH
               IF (LCONS(3)) THEN
*                 eta0 FOR THE LORENTZIAN FRACTION (HIGHER ANGLE)
                  CALL DERCOM(DERIV, 11, DFDEH)
*                 eta1 FOR THE LORENTZIAN FRACTION (HIGHER ANGLE)
                  CALL DERCOM(DERIV, 12, DFDEH*PEAK(I))
               ELSE IF (LD2(3) .EQ. 1) THEN
                  DERIV(LR2(3)) = DERIV(LR2(3)) + YPEAK(I)*DFDEH
               END IF
	       
*              |Fc|
*              The derivative of the model function, f, with respect to |Fc|
*              df/d|Fc| = (df/d|Fc|**2)*(d|Fc|**2/d|Fc|)
*                       = (df/d|Fc|**2)*2*|Fc|

               IF (LINTEG(4)) DERIV(LR2(4)) =
*    &         DERIV(LR2(4)) + PRFIK*YPEAK(I)/FF(I)*2.0*SQRT(FF(I))
     &         DERIV(LR2(4)) + 2.0*PRFIK*YPEAK(I)/SQRT(FF(I))

*              2-theta_k: 2-theta of reflection k at the peak position
*              Delta.2-theta_{ik}: 2-theta_i - 2-theta_k
*              i: step number
*              dG/d(2-theta_k) 
*              = dG/d(Delta.2-theta_{ik}) * d(Delta.2-theta_{ik})/d(2-theta_k)
*              = -dG/d(Delta.2-theta_{ik})
*              = -DFDX

*              2-theta(peak)
               IF (LPEAKPOS(5)) DERIV(LR2(5)) = DERIV(LR2(5)) -
     &         YPEAK(I)*DFDX
         
            ELSE
      
*              The modified split-type pseudo-Voigt function is used for
*              relaxed reflections.
*              ID should not be set at 3 in this case.

*              SQH2L
               IF (LD2(0) .EQ. 1) 
     &         DERIV(LR2(0)) = DERIV(LR2(0)) + YPEAK(I)*DFDW1
     
*              SQH2G
               IF (LD2(1) .EQ. 1) 
     &         DERIV(LR2(1)) = DERIV(LR2(1)) + YPEAK(I)*DFDW2

*              ASYM
               IF (LD2(2) .EQ. 1)
     &         DERIV(LR2(2)) = DERIV(LR2(2)) + YPEAK(I)*DFDA

*              ETAL
               IF (LD2(3) .EQ. 1) 
     &         DERIV(LR2(3)) = DERIV(LR2(3)) + YPEAK(I)*DFDEL
         
*              ETAH
               IF (LD2(4) .EQ. 1)
     &         DERIV(LR2(4)) = DERIV(LR2(4)) + YPEAK(I)*DFDEH

*              |Fc|
               IF (LINTEG(5)) DERIV(LR2(5)) = 
     &         DERIV(LR2(5)) + 2.0*PRFIK*YPEAK(I)/SQRT(FF(I))
         
*              2-theta(peak)
               IF (LPEAKPOS(6)) DERIV(LR2(6)) = DERIV(LR2(6)) -
     &         YPEAK(I)*DFDX
            END IF

         CASE (3)
            CALL DERSHIFT(TTH,DERIV)
            
            CALL DERWA(TTH,DERIV)

*           m_L
            IF (LCONS(2)) THEN
*              m_L0 (LOWER ANGLE)
               CALL DERCOM(DERIV, 9, DFDRL*DMPART(1,I))
*              m_L1 (LOWER ANGLE)
               CALL DERCOM(DERIV, 10, DFDRL*DMPART(1,I)*PEAK(I))
            ELSE IF (LD2(2) .EQ. 1) THEN
               DERIV(LR2(2)) = DERIV(LR2(2)) + YPEAK(I)*DFDRL
            END IF
         
*           m_H
            IF (LCONS(3)) THEN
*              m_H0 (HIGHER ANGLE)
               CALL DERCOM(DERIV, 11, DFDRH*DMPART(2,I))
*              m_H1 (HIGHER ANGLE)
               CALL DERCOM(DERIV, 12, DFDRH*DMPART(2,I)*PEAK(I))
            ELSE IF (LD2(3) .EQ. 1) THEN
               DERIV(LR2(3)) = DERIV(LR2(3)) + YPEAK(I)*DFDRH
            END IF
	    
*           |Fc|
            IF (LINTEG(4)) DERIV(LR2(4)) =
     &      DERIV(LR2(4)) + 2.0*PRFIK*YPEAK(I)/SQRT(FF(I))

*           2-theta(peak)
            IF (LPEAKPOS(5)) DERIV(LR2(5)) = DERIV(LR2(5)) -
     &      YPEAK(I)*DFDX
      END SELECT
      END

************************************************************************

      SUBROUTINE DERWA(TTH,DERIV)
*     DERIVATIVE WITH RESPECT TO FWHM AND ASYMMETRY PARAMETERS
      PARAMETER (NB=7000,NT=999,NPH=8)
      REAL DERIV(*)
      LOGICAL LCONS
      COMMON /ANISBR/ HANIS(NPH),KANIS(NPH),LANIS(NPH),
     &  CANIS(NPH),PCOS(NB)
      COMMON /DPRSPV/ DFDX,DFDA,DFDW,DFDW1,DFDW2,DFDEL,DFDEH,DFDRL,DFDRH
      COMMON /G/ I
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /ZSPFN/ PEAKSH(NB),SQH2(NB),SQH2L(NB),SQH2G(NB),ASYM(NB),
     &ETAL(NB),ETAH(NB),FNORM1(2,NB),COMLG(2,NB),COML(2,NB),COMG(2,NB),
     &COEFL(2,NB),COEFG(2,NB),DECAY(2,NB),COMPEAR(2,NB),CDX(2,NB),
     &CDA1(2,NB),CDA2(2,NB),CDW(2,NB),CDRL1(2,NB),CDRL2(2,NB),
     &CDRH1(2,NB),DMPART(2,NB)
      COMMON /UNCONS/ NUCREF,IPHASE(90),IHKL(3,90),NOPPP(90),LPPP(NB),
     &  RCUT(2,90),NCUT
*     NOPPP(LPPP(I)) + J: SERIAL PARAMETER NUMBER FOR THE J'TH 
*     PRIMARY PROFILE PARAMETER (J = 0 - (3 OR 4)) FOR THIS REFLECTION 
*     SPECIFIED BY THE USER
      LD2(J) = ID(NOPPP(LPPP(I)) + J)
      LR2(J) = IR(NOPPP(LPPP(I)) + J)
*     LCONS(J) = .TRUE.: The (J+1)-th primary profile parameter is
*     conventionally calculated from secondary profile parameters.
*     Reflection I is either constrained (LPPP(I) = 0) or relaxed 
*     (LPPP(I) > 0).  In the latter case, LD2(J+1) = 3.

*     LCONS(J) = .FALSE.: The (J+1)-th primary profile parameter is
*     input and refined/fixed.

      LCONS(J) = LPPP(I) .EQ. 0 .OR. 
     &           (LPPP(I) .GT. 0 .AND. LD2(J) .EQ. 3)

*     SQH2
*     0.5/SQRT(U*TTH**2 + V*TTH + W) MUST BE MULTIPLIED BECAUSE 
*     FWHM = SQRT((U + Ue*cos(phi))*TTH**2 + V*TTH + W + 
*            Pe*cos(phi)/(cos(theta))**2
      IF (LCONS(0)) THEN
*        FWHM PARAMETER, U
         CALL DERCOM(DERIV, 1, 0.5/SQH2(I)*DFDW*TTH*TTH)
*        FWHM PARAMETER, V
         CALL DERCOM(DERIV, 2, 0.5/SQH2(I)*DFDW*TTH)
*        FWHM PARAMETER, W
         CALL DERCOM(DERIV, 3, 0.5/SQH2(I)*DFDW)
*        PARAMETERS RELATED TO ANISOTROPIC BROADENING
*        PCOS(I) should be squared
*        ANISOTROPIC STRAIN BROADNING, Ue
*        CALL DERCOM(DERIV, 13, 0.5/SQH2(I)*DFDW*PCOS(I)*TTH*TTH)
         CALL DERCOM(DERIV, 13, 0.5/SQH2(I)*DFDW*PCOS(I)**2*TTH*TTH)
*        ANISOTROPIC SCHERRER BROADNING, Pe
         CALL DERCOM(DERIV, 14, 
*    &   0.5/SQH2(I)*DFDW*PCOS(I)/(COS(TH(I)))**2)
     &   0.5/SQH2(I)*DFDW*(PCOS(I)/COS(TH(I)))**2)
      ELSE IF (LD2(0) .EQ. 1) THEN
         DERIV(LR2(0)) = DERIV(LR2(0)) + YPEAK(I)*DFDW
      END IF
     
*     PROFILE ASYMMETRY
      IF (LCONS(1)) THEN
*        a0
         CALL DERCOM(DERIV, 5, DFDA)
*        a1
         COSCTH = 1.0/SIN(TH(I))
         CALL DERCOM(DERIV, 6, DFDA*(1.4142136 - COSCTH))
*        a2
         CALL DERCOM(DERIV, 7, DFDA*(2.0 - COSCTH**2))
      ELSE IF (LD2(1) .EQ. 1) THEN
         DERIV(LR2(1)) = DERIV(LR2(1)) + YPEAK(I)*DFDA
      END IF
      END

************************************************************************

      SUBROUTINE DERSHIFT(TTH,DERIV)
*     CALCULATE DERIVATIVES WITH RESPECT TO PEAK-SHIFT PARAMETERS
      PARAMETER (NB=7000,NT=999,NPH=8)
      REAL DERIV(*),DERSH(4)
      COMMON /DPRSPV/ DFDX,DFDA,DFDW,DFDW1,DFDW2,DFDEL,DFDEH,DFDRL,DFDRH
      COMMON /FT/ FT1,FT2,PEAKNOR(NB),TANTHNOR(NB) 
      COMMON /G/ I
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
            
      YPD = YPEAK(I)*DFDX
      IF (NSHIFT .LE. 3 .AND. ID(1) .EQ. 1) DERIV(1) = DERIV(1) + YPD
      SELECT CASE (NSHIFT)
         CASE (1)
            IF (ID(2) .EQ. 1) DERIV(IR(2)) = 
     &                        DERIV(IR(2)) + YPD*COS(PEAK(I))
            IF (ID(3) .EQ. 1) DERIV(IR(3)) = 
     &                        DERIV(IR(3)) + YPD*SIN(PEAK(I))
            IF (ID(4) .EQ. 1) DERIV(IR(4)) = DERIV(IR(4)) + YPD*TTH
         CASE (2)
            IF (ID(2) .EQ. 1) DERIV(IR(2)) = DERIV(IR(2)) + YPD*PEAK(I)
            IF (ID(3) .EQ. 1) DERIV(IR(3)) = DERIV(IR(3)) + 
     &                                       YPD*PEAK(I)**2
            IF (ID(4) .EQ. 1) DERIV(IR(4)) = DERIV(IR(4)) + 
     &                                       YPD*PEAK(I)**3
         CASE (3)
            IF (ID(2) .EQ. 1) DERIV(IR(2)) = DERIV(IR(2)) + YPD*TTH
            IF (ID(3) .EQ. 1) DERIV(IR(3)) = DERIV(IR(3)) + YPD*TTH**2
            IF (ID(4) .EQ. 1) DERIV(IR(4)) = DERIV(IR(4)) + YPD*TTH**3
         CASE (4)
            CALL FLEG(PEAKNOR(I), DERSH, 4) 
            DO J = 1, 4
               IF (ID(J) .EQ. 1) 
     &         DERIV(IR(J)) = DERIV(IR(J)) + YPD*DERSH(J)
            END DO
         CASE (5) 
            CALL FLEG(TANTHNOR(I), DERSH, 4) 
            DO J = 1, 4
               IF (ID(J) .EQ. 1) DERIV(IR(J)) = DERIV(IR(J)) + 
     &                                          YPD*DERSH(J)
            END DO
      END SELECT
      END

************************************************************************

      SUBROUTINE DERCOM(DERIV, J, X)
*     CALCULATE THE DERIVATIVE WITH RESPECT TO A PROFILE PARAMETER
      PARAMETER (NB=7000,NT=999,NPH=8,NAP=150)
      INTEGER H
      REAL DERIV(*)
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /G/ I
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      LP(J) = KPHB(NOPH(I))+J
      LD(J) = ID(LP(J))
      LR(J) = IR(LP(J))
      LR2(J) = IR(KPHB(1)+J)

      IF (LD(J) .EQ. 0) THEN
         RETURN
      ELSE IF (LD(J) .EQ. 1) THEN
         DERIV(LR(J)) = DERIV(LR(J)) + YPEAK(I)*X
      ELSE IF (ID(KPHB(1)+J) .EQ. 1) THEN
*        THIS PROFILE PARAMETER IS CONSTRAINED TO BE EQUAL TO THAT FOR
*        THE FIRST PHASE (ID = 2)
         DERIV(LR2(J)) = DERIV(LR2(J)) + YPEAK(I)*X
      END IF
      END

************************************************************************

      SUBROUTINE DERSR(THETA,A,DERIV,ABSBI)
*     DERIVATIVES WITH RESPECT TO SURFACE-ROUGHNESS PARAMETERS
      PARAMETER (NP=80000,NT=999)
      REAL A(*),DERIV(*)
      COMMON /ABSVR/ ABSORP(NP),VARLEN(NP),ROUGH(NP),NSURFR
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      
      P = A(NROUGH)
      Q = A(NROUGH+1)
      R = A(NROUGH+2)
      T = A(NROUGH+3)
      SELECT CASE (NSURFR)
         CASE (1)
            SINTH = SIN(THETA)
            IF (ID(NROUGH) .EQ. 1) DERIV(IR(NROUGH)) = ABSBI*
     &      (R*(-EXP(-Q) + EXP(-Q/SINTH)))
            IF (ID(NROUGH+1) .EQ. 1) DERIV(IR(NROUGH+1)) = ABSBI*
     &      (R*P*(EXP(-Q) - EXP(-Q/SINTH)/SINTH))
            IF (ID(NROUGH+2) .EQ. 1) DERIV(IR(NROUGH+2)) = ABSBI*
     &      (1.0 - P*(EXP(-Q) - EXP(-Q/SINTH)) - 1.0 - 
     &      T*(THETA - 1.570796))
            IF (ID(NROUGH+3) .EQ. 1) DERIV(IR(NROUGH+3)) = ABSBI*
     &      (1.0 - R)*(THETA - 1.570796)
         CASE (2)
            IF (ID(NROUGH+3) .EQ. 1) DERIV(IR(NROUGH+3)) = 
     &      ABSBI*(-THETA + 1.570796)
         CASE (3)
            SINTH = SIN(THETA)
            IF (ID(NROUGH) .EQ. 1) DERIV(IR(NROUGH)) = ABSBI*
     &      (-EXP(-Q) + EXP(-Q/SINTH))
            IF (ID(NROUGH+1) .EQ. 1) DERIV(IR(NROUGH+1)) = ABSBI*
     &      (P*(EXP(-Q) - EXP(-Q/SINTH)/SINTH))
         CASE (4)
            SINTH = SIN(THETA)
            IF (ID(NROUGH) .EQ. 1) DERIV(IR(NROUGH)) = ABSBI*Q*
     &      (-1.0 + Q - (1.0 - Q/SINTH)/SINTH)   
            IF (ID(NROUGH+1) .EQ. 1) DERIV(IR(NROUGH+1)) = ABSBI*
     &      Q*(-1.0 + 2.0*Q + (-1.0 + 2.0*Q/SINTH)/SINTH)
      END SELECT
      END
            
************************************************************************

      SUBROUTINE FLEG(X,PL,NL)
*     FIT LEGENDRE POLYNOMIALS UP TO SOME ORDER N1-1 THROUGH A DATA SET
      INTEGER NL
      REAL X,PL(NL)
      INTEGER J
      REAL D,F1,F2,TWOX

      PL(1)=1.
      PL(2)=X
      IF (NL .GT. 2) THEN
        TWOX=2.*X
        F2=X
        D=1.
        DO J=3,NL
          F1=D
          F2=F2+TWOX
          D=D+1.
          PL(J)=(F2*PL(J-1)-F1*PL(J-2))/D
        END DO
      END IF
      END

************************************************************************

      SUBROUTINE CHLSK(N,A,G,EPS)
C     CHOLESKI DECOMPOSITION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*),G(*)

      IF(A(1).LE.EPS) GO TO 100
      G(1)=SQRT(A(1))
      DO 10 I=2,N
         I1=I*(I-1)/2+1
         G(I1)=A(I1)/G(1)
   10 CONTINUE
      DO 60 J=2,N
         JJ=J*(J-1)/2
         SUM=A(JJ+J)
         JMINUS=J-1
         DO 20 K=1,JMINUS
            SUM=SUM-G(JJ+K)**2
   20    CONTINUE
         IF(SUM.LE.EPS) GO TO 100
         G(JJ+J)=SQRT(SUM)
         IF(J.EQ.N) GO TO 60
         JPLUS=J+1
         DO 50 I=JPLUS,N
            II=I*(I-1)/2
            SUM=A(II+J)
            DO 40 K=1,JMINUS
               SUM=SUM-G(II+K)*G(JJ+K)
   40       CONTINUE
            G(II+J)=SUM/G(JJ+J)
   50    CONTINUE
   60 CONTINUE
      RETURN
  100 CALL JOBEND('The coefficient matrix is not positive definite')
      END

************************************************************************

      SUBROUTINE RSSE(A,B,DELTAB,D,X,Y,YFIT,AMIN,NPTS,NSTEP,
     &NTERMS,CSMIN,CHISQR,DAMP,IDPOS,IFAIL,NFNEV,FLAMDA,FLAMDP)
C     UPDATE PARAMETERS AND CALCULATE NEW CHI-SQUARE
      PARAMETER (NT=999)
      DOUBLE PRECISION DELTAB(*)
      DIMENSION A(*),B(*),D(*),X(*),Y(*),YFIT(*),AMIN(*)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      EXTERNAL UNPACK

      DO J = 1, NRFN
C        COPY PACKED ARRAY A INTO ARRAY B
         B(J)=A(J)
C        ALTER PARAMETERS
         B(J)=B(J)+DAMP*DELTAB(J)
      END DO
      CALL UNPACK(B,D,NTERMS)
C
C     CALCULATE CHISQR USING NEW PARAMETERS
      CHISQR=AOFFN(X,Y,YFIT,NPTS,NSTEP,B,1,IDPOS)
      NFNEV=NFNEV+1
C     FLAMDP: MARQUARDT PARAMETER ADOPTED IN THE PRESENT CYCLE
      FLAMDP=FLAMDA
      IF(CHISQR.LT.CSMIN) THEN
         CSMIN=CHISQR
         DO J = 1, NTERMS
            AMIN(J)=B(J)
         END DO
         IFAIL=0
      ENDIF
      END

************************************************************************

      SUBROUTINE BACKSP(LINE)
*     RETURN TO THE LINE WHERE THE FIRST DATUM IS WRITTEN
      CHARACTER*80 LINE,LINE2
      PARAMETER (N=4)

    1 CONTINUE
         BACKSPACE N
         READ(N,'(A)') LINE2
         IF (LINE .EQ. LINE2) THEN
            BACKSPACE N
            RETURN
         END IF
      GO TO 1
      END

************************************************************************

      INTEGER FUNCTION LAPO(LINE)
*     LAPO=1: THE FIRST CHARACTER IS AN APOSTROPHE
      CHARACTER*80 LINE

      DO I = 1, 80
         IF (LINE(I:I) .NE. ' ') EXIT
      END DO
      IF (LINE(I:I) .EQ. '''') THEN
         LAPO=1
      ELSE
         LAPO=0
      END IF
      END

************************************************************************

      SUBROUTINE SOLVE(N,G,X,B)
C     SOLVE CONSISTENT SETS OF EQUATIONS WHOSE COEFFICIENT MATRICES
C     ARE SYMMETRIC AND NON-NEGATIVE DEFINITE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(*),X(*),B(*)
      X(1)=B(1)/G(1)
      DO 20 I=2,N
         S=B(I)
         II=I*(I-1)/2
         IMINUS=I-1
         DO 10 J=1,IMINUS
            S=S-G(II+J)*X(J)
   10    CONTINUE
         X(I)=S/G(II+I)
   20 CONTINUE
      X(N)=X(N)/G(N*(N+1)/2)
      DO 40 IBACK=2,N
         I=N+1-IBACK
         S=X(I)
         IPLUS=I+1
         JI=I*(I+3)/2
         DO 30 J=IPLUS,N
            S=S-G(JI)*X(J)
            JI=JI+J
   30    CONTINUE
         X(I)=S/G(I*(I+1)/2)
   40 CONTINUE
      END

************************************************************************

      SUBROUTINE LAMBDA(DELTAB,BETA,AD,G,FLAMDA,FLAMDC,IDPOS,LDM,
     &CHISQ1,CHISQR)
C     FLETCHER'S ALGORITHM FOR DETERMINING FLAMDA
      PARAMETER (NT=999)
      DOUBLE PRECISION DELTAB(*),BETA(*),AD(*),G(*)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
C     CALCULATE PREDICTED DECREASE IN CHI-SQUARE ASSUMING LINEARITY
      PRED=0.0
      DO 10 J=1,NRFN
         PRED=PRED+DELTAB(J)*(BETA(J)+FLAMDA*DELTAB(J)*AD(J))
   10 CONTINUE
      IF(PRED.LE.0.0) THEN
         R1=0.5
      ELSE
         R1=(CHISQ1-CHISQR)/PRED
      ENDIF
      IF(R1.GT.0.75 .AND. FLAMDA.NE.0.0) THEN
         FLAMDA=0.5*FLAMDA
         IF(FLAMDA.LT.FLAMDC) FLAMDA=0.0
      ELSEIF(R1.LT.0.25) THEN
         IF(IDPOS.EQ.1) THEN
            U=MAX(2.0,MIN(2.0-R1,10.0))
         ELSE
            U=10.0
         ENDIF
         IF(FLAMDA.EQ.0.0) THEN
C           CALCULATE THE SUM OF THE DIAGONAL TERMS OF THE INVERSE
C           OF THE COEFFICIENT MATRIX
            CALL MATINV(G,NRFN)
            TRACE=0.0
            DO 20 J=1,NRFN
               TRACE=TRACE+G(J*(J+1)/2)*AD(J)
   20       CONTINUE
            FLAMDC=1.0/TRACE
            FLAMDA=FLAMDC
            U=U*0.5
         ENDIF
         FLAMDA=FLAMDA*U
         IF(FLAMDA.GT.1E4) LDM=1
      ENDIF
      END

************************************************************************

      SUBROUTINE MATINV(M,N)
C     INVERSION OF MATRIX M WHICH HAS BEEN DECOMPOSED BY THE
C     CHOLESKY METHOD
      DOUBLE PRECISION M(*),AM
      I0=0
      DO 30 I=1,N
      I0=I0+I-1
      J0=I0-I+1
      DO 31 J=I,N
      AM=0.0D0
      IF(J.EQ.I) AM=1.0D0
      J0=J0+J-1
      KI=I0+1
      DO 31 K=I,J
      JK=J0+K
      KI=KI+K-1
      IF(K.NE.J) AM=AM-M(JK)*M(KI)
 31   IF(K.EQ.J) M(KI)=AM/M(JK)
 30   CONTINUE
C     INVERSION OF MATRIX M
      IJ=0
      I0=0
      DO 40 I=1,N
      I0=I0+I-1
      DO 40 J=1,I
      AM=0.0D0
      K0=I0-I+1
      DO 41 K=I,N
      K0=K0+K-1
      KI=K0+I
      KJ=K0+J
 41   AM=AM+M(KI)*M(KJ)
      IJ=IJ+1
 40   M(IJ)=AM
      END

************************************************************************

      SUBROUTINE DERSP(A,DSTFAC,ISFP)
C     CALCULATE ARRAYS USED FOR THE CALCULATION OF PARTIAL DERIVATIVES
C     OF PARAMETERS.
      INTEGER R,RX,RY,RZ,H,HP
      PARAMETER (NB=7000,NR=400,NSF=400,NT=999,NCS=400,NPH=8,NAP=150,
     &  NS=48)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      DIMENSION A(*),DSTFAC(NSF,NB),ISFP(*),DPARA(NT)
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /B/ U(NB),NREFL
      COMMON /CNSTR/ C1(NCS),C2(10,NCS),IN(10,NCS),NCNPAR(NCS),
     &  NACNS(NCS),NCNSTR
      COMMON /DERPR/ POP2(NB),COMPO(NB),CS(NB),CCC(6,NB)
      COMMON /G/ I
      COMMON /H/ NQ(NPH),NCENTR(NPH),LAUEG(NPH)
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /LTM/ LISO(NAP,NPH),LMAG(NPH)
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /RT/ T(3,NS,NPH),HT(NS,NB),NSYM(NPH),NSITE(NPH),
     &  IDSYM(NS,NAP,NPH),R(3,3,NS,NPH),RX(NS,NB),RY(NS,NB),RZ(NS,NB)
      COMMON /SFRI/ SFR(NB),SFI(NB),OTF(NB)
      COMMON /ABPI/ AP(NB),BP(NB)
      COMMON /T/ FF(NB),COEF(NPH)
      COMMON /UCC/ IDPEAK(NPH),IDSF(NPH),IDPRF(NPH),IDPO(NPH)
*     LP ==> LP2 TO PREVENT THE OVERLAP OF THE NAME
      COMMON /V/ NPROR(NPH),HP(NPH),KP(NPH),LP2(NPH),LSUM(NPH)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
      LP(J)=KPHB(NOPH(I))+J
      P(J)=A(LP(J))
      LD(J)=ID(LP(J))

C     ASSIGN NUMBERS FOR REFINABLE STRUCTURE-FACTOR PARAMETERS TO
C     ARRAY ISFP
      IS=0
      DO 20 J=1,NPHASE
         DO 10 JJ=KPHB(J)+NOTX,KPHE(J)
            IF(ID(JJ).EQ.1) THEN
               IS=IS+1
               IF(IS.GT.NSF) CALL JOBEND('NSF is too small')
               ISFP(JJ)=IS
            ENDIF
   10    CONTINUE
   20 CONTINUE
C
      DO I = 1, NREFL

*        SCALE FACTOR
         IF(LD(NSFX).EQ.1) CS(I)=YPEAK(I)/P(NSFX)

*        PREFERRED ORIENTATION
         IF (LD(NPO1X) .EQ. 1) COMPO(I)=YPEAK(I)/PROR1(I)*DPROR(I)
         IF (LD(NPO2X) .EQ. 1) POP2(I)=YPEAK(I)/PROR1(I)*POP2(I)         

*        METRIC TENSORS
         IF(IDPEAK(NOPH(I)).EQ.1) THEN
*           -D(I)**2*TAN(TH(I)): COEFFICIENT IN d(TH(I))/dX, WHERE X IS
*                                AN ELEMENT OF THE METRIC TENSOR
            CT = -D(I)**2*TAN(TH(I))*YPEAK(I)
            IL=LAUEG(NOPH(I))
            IF (LD(NAX) .EQ. 1) THEN
               IF (IL .GE. 14 .OR. IL .EQ. 8 .OR. IL .EQ. 10) THEN
                  CCC(1,I)=CT*FLOAT(H(I)**2+K(I)**2+L(I)**2)
               ELSE IF (IL .EQ. 6 .OR. IL .EQ. 7) THEN
                  CCC(1,I)=CT*FLOAT(H(I)**2+K(I)**2)
               ELSE IF (IL .EQ. 9 .OR. IL .EQ. 11 .OR. IL .EQ. 12 .OR.
     &         IL .EQ. 13) THEN
                  CCC(1,I)=CT*FLOAT(H(I)**2+K(I)**2+H(I)*K(I))
               ELSE
                  CCC(1,I)=CT*FLOAT(H(I)**2)
               ENDIF
            END IF
            IF (LD(NBX) .EQ. 1) CCC(2,I)=CT*FLOAT(K(I)**2)
            IF (LD(NCX) .EQ. 1) CCC(3,I)=CT*FLOAT(L(I)**2)
            IF (LD(NALPX) .EQ. 1) THEN
               IF (IL .EQ. 8 .OR. IL .EQ. 10) THEN
                  CCC(4,I)=
     &            (CT+CT)*FLOAT(K(I)*L(I)+L(I)*H(I)+H(I)*K(I))
               ELSE
                  CCC(4,I)=(CT+CT)*FLOAT(K(I)*L(I))
               END IF
            END IF
            IF (LD(NBETX) .EQ. 1) CCC(5,I)=(CT+CT)*FLOAT(H(I)*L(I))
            IF (LD(NGAMX) .EQ. 1) CCC(6,I)=(CT+CT)*FLOAT(H(I)*K(I))
         END IF
C
         IF(IDSF(NOPH(I)).EQ.0) CYCLE
	 
C        DSTFAC: PARTIAL DERIVATIVES OF YPEAK(I) WITH RESPECT TO
C                STRUCTURE PARAMETERS
	 
         IF(LD(NOTX).EQ.1)
     &      DSTFAC(ISFP(LP(NOTX)),I) = -2.0*RLV2(I)*YPEAK(I)
         IF (L12(I) .EQ. 0) THEN
            CALL ANALD(A,DPARA,RX(1,I),RY(1,I),RZ(1,I),HT(1,I),
     &      NSITE(NOPH(I)),NSYM(NOPH(I)),IDSYM(1,1,NOPH(I)))
         END IF
	 
         DO J = LP(NG1X), KPHE(NOPH(I))
            IF (ID(J) .NE. 1) CYCLE
            IF (L12(I) .EQ. 0) THEN
               TEMP = DPARA(J)
C              CONTRIBUTION OF THE CONSTRAINED PARAMETERS
               DO J2 = 1, NCNSTR
                  DO J3 = 1, NCNPAR(J2)
                     IF (IN(J3,J2) .EQ. J)
     &               TEMP = TEMP + C2(J3,J2)*DPARA(NACNS(J2))
                  END DO
               END DO
               IF (R12 .GT. 0.0) THEN
                  DSTFAC(ISFP(J),I) = TEMP/FF(I)
               ELSE
                  DSTFAC(ISFP(J),I) = TEMP/FF(I)*YPEAK(I)
               END IF
               IF (NQ(NOPH(I)) .GE. 3) DSTFAC(ISFP(J),I)=
     &         DSTFAC(ISFP(J),I)*OTF(I)
            ELSE
               DSTFAC(ISFP(J),I) = DSTFAC(ISFP(J),L12(I))*YPEAK(I)
C              DSTFAC FOR THE K-ALPHA1 PEAK IS ALSO CALCULATED
               DSTFAC(ISFP(J),L12(I)) = DSTFAC(ISFP(J),L12(I))*
     &         YPEAK(L12(I))
            ENDIF
         END DO
      END DO
      END

************************************************************************

      SUBROUTINE INTTBL(VNS,L,NSPGR,NSET,LAUEG,NCENTR,NSYM,SPGR,LSPSYM,
     &  NSSYM,SYM,LINNSS,NLNSS)
*     READ INFORMATION ABOUT SPACE GROUP
      PARAMETER (NS=48,NPH=8)
      INTEGER NSPGR(*),NSET(*),LAUEG(*),NCENTR(*),NSYM(*),LSPSYM(*),
     &  NSSYM(14,NPH)
      CHARACTER VNS*7,FORM(3)*4,SPGR(*)*10,COORD*20,
     &  SYM(3,NS,NPH)*10,LINNSS(*)*40,LINE*80
      DATA FORM/'(I1)','(I2)','(I3)'/

C        LAUEG: LAUE GROUP NUMBER
C              LAUEG= 1  TRICLINIC (-1)   1/2
C                   = 2  MONOCLINIC-1 (2/M) A-AXIS   1/4
C                   = 3  MONOCLINIC-2 (2/M) B-AXIS   1/4
C                   = 4  MONOCLINIC-3 (2/M) C-AXIS   1/4
C                   = 5  ORTHORHOMBIC (MMM)   1/8
C                   = 6  TETRAGONAL-1 (4/M)   1/8
C                   = 7  TETRAGONAL-2 (4/MMM)   1/16
C                   = 8  TRIGONAL-1  RHOMBOHEDRAL LATTICE (-3)   1/6
C                   = 9  TRIGONAL-2  HEXAGONAL LATTICE (-3)   1/6
C                   =10  TRIGONAL-3  RHOMBOHEDRAL LATTICE (-3M)   1/12
C                   =11  TRIGONAL-4  HEXAGONAL LATTICE (-3M)   1/12
C                   =12  HEXAGONAL-1 (6/M)   1/12
C                   =13  HEXAGONAL-2 (6/MMM)   1/24
C                   =14  CUBIC-1 (M3)   1/24
C                   =15  CUBIC-2 (M3M)   1/48
C     NCENTR=0: LACKING A CENTER OF SYMMETRY
C     NCENTR=1: POSSESSING A CENTER OF SYMMETRY
C         NSYM: NUMBER OF EQUIVALENT COORDINATES
C         SPGR: NAME OF THE SPACE GROUP

*     VNS: VOL.#-SPACE GROUP #-SETTING # (E.G., 'A-123-2', 'I-56')
      IF (VNS(1:1) .NE. 'I' .AND. VNS(1:1) .NE. 'A') THEN
         CALL JOBEND('Wrong volume name in: '//VNS)
      ELSE IF (VNS(2:2) .NE. '-') THEN
         CALL JOBEND('Invalid delimiter in: '//VNS)
      END IF
      MINUS=INDEX(VNS(3:),'-')
      IF (MINUS .NE. 0) THEN
         READ(VNS(3:),FORM(MINUS-1)) NSPGR(L)
         READ(VNS(MINUS+3:),'(I1)') NSET(L)
      ELSE
         DO 20 I=3,7
            IF (VNS(I:I) .EQ. ' ') GO TO 2
   20    CONTINUE
    2    READ(VNS(3:I-1),FORM(I-3)) NSPGR(L)
         NSET(L)=1
      END IF
      IF(NSPGR(L) .LT. 1 .OR. NSPGR(L) .GT. 230)
     &CALL JOBEND('Number of space group out of range')
      CALL OPENSPGR(VNS(1:1),NFILE)
      REWIND NFILE
    4 READ(NFILE,'(5I3,3X,A)',IOSTAT=IOS)
     &NSPGR1,NSET1,LAUEG(L),NCENTR(L),NSYM(L),SPGR(L)
      IF (IOS .LT. 0) THEN
         CLOSE(UNIT=NFILE,STATUS='KEEP')
         CALL JOBEND('Bad space group specification: '//VNS)
      ELSE IF (NSPGR(L) .NE. NSPGR1 .OR. NSET(L) .NE. NSET1) THEN
         DO 40 J=1,NSYM(L)+1
            READ(NFILE,'()')
   40    CONTINUE
         GO TO 4
      END IF
      IF (NSYM(L) .GT. NS) THEN
         CALL JOBEND('Too small NS value')
      ELSE IF (LSPSYM(L) .LE. 1) THEN
         READ(NFILE,'(14I2)') (NSSYM(I,L),I=1,14)
         DO 50 J=1,NSYM(L)
            READ(NFILE,'(A20)') COORD
            CALL DECODE(COORD,SYM(1,1,L),J,L)
   50    CONTINUE
      ELSE IF (LSPSYM(L) .EQ. 2) THEN
         READ(4,*) NCENTR(L),SPGR(L)
         DO 55 J=1,14
            NSSYM(J,L)=0
   55    CONTINUE
         DO 60 J=1,14
            READ(4,'(A)') LINNSS(J)
            IF (LINNSS(J)(1:1) .EQ. '}') GO TO 6
            CALL NSCONV(LINNSS(J),NSSYM(1,L))
   60    CONTINUE
         READ(4,'(A)') LINE
         IF (LINE(1:1) .NE. '}') CALL JOBEND
     &   ('Too many non-standard symmetry lines')
    6    NLNSS=J-1
         DO 70 J=1,NS
            READ(4,'(A)') LINE
            IF (LINE(1:1) .EQ. '}') GO TO 7
            READ(LINE,'(A20)') COORD
            CALL DECODE(COORD,SYM(1,1,L),J,L)
   70    CONTINUE
    7    NSYM(L)=J-1
      ELSE
         CALL JOBEND('Bad LSPSYM value')
      END IF
      IF (LAUEG(L) .EQ. 8 .OR. LAUEG(L) .EQ. 10) NSSYM(7,L)=21
      END

************************************************************************

      SUBROUTINE PARAM(ANAME,ISONUM,A,ID,NMODE,NTERMS,LABEL,LPAR,NLINE)
C     ASSIGN REMARKS ON PARAMETERS
      PARAMETER (NT=999,NR=400,NA=15,NAP=150,NPH=8,NS=48,NB=7000)
      CHARACTER ANAME(*)*5,PARNAM*60,PHNAME*25,LABEL(*)*25
*Rev 1.0t 2001.02.24 Izumi {
C     CHARACTER TOPCOL(NT)*16
      CHARACTER TOPCOL(NT)*26
* }
      INTEGER ISONUM(NAP,NPH),ID(*),R,RX,RY,RZ,LPAR(*),NLINE,HP
      REAL A(*)
      COMMON /H/ NQ(NPH),NCENTR(NPH),LAUEG(NPH)
      COMMON /LTM/ LISO(NAP,NPH),LMAG(NPH)
      COMMON /MAGSC/ MAGATM(NA),LCMFF(NA),CMFF(7,NA)
      COMMON /PAR/ PARNAM(NT),PHNAME(NPH)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /RT/ T(3,NS,NPH),HT(NS,NB),NSYM(NPH),NSITE(NPH),
     &  IDSYM(NS,NAP,NPH),R(3,3,NS,NPH),RX(NS,NB),RY(NS,NB),RZ(NS,NB)
      COMMON /SITE/ NOAT(NAP,NPH)
      COMMON /V/ NPROR(NPH),HP(NPH),KP(NPH),LP2(NPH),LSUM(NPH)
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      LP(J)=KPHB(L)+J
 
      IF (NPRFN .EQ. 0) THEN
         PARNAM(1)='Zero-point shift, Z'
         PARNAM(2)='Specimen-displacement parameter, Ds'
         PARNAM(3)='Specimen-transparency parameter, Ts'
         PARNAM(4)='Dummy'
      ELSE
         PARNAM(1)='Peak-shift parameter, t0'
         PARNAM(2)='Peak-shift parameter, t1'
         PARNAM(3)='Peak-shift parameter, t2'
         PARNAM(4)='Peak-shift parameter, t3'
      END IF
      PARNAM(NROUGH)='Surface-roughness parameter, q0'
      PARNAM(NROUGH+1)='Surface-roughness parameter, q1'
      PARNAM(NROUGH+2)='Surface-roughness parameter, q2'
      PARNAM(NROUGH+3)='Surface-roughness parameter, q3'
      PARNAM(NBCKGR)='Background parameter, b0'
      PARNAM(NBCKGR+1)='Background parameter, b1'
      PARNAM(NBCKGR+2)='Background parameter, b2'
      PARNAM(NBCKGR+3)='Background parameter, b3'
      PARNAM(NBCKGR+4)='Bacgground parameter, b4'
      PARNAM(NBCKGR+5)='Background parameter, b5'
      PARNAM(NBCKGR+6)='Background parameter, b6'
      PARNAM(NBCKGR+7)='Background parameter, b7'
      PARNAM(NBCKGR+8)='Background parameter, b8'
      PARNAM(NBCKGR+9)='Background parameter, b9'
      PARNAM(NBCKGR+10)='Background parameter, b10'
      PARNAM(NBCKGR+11)='Background parameter, b11'
      
*     PRIMARY PROFILE PARAMETERS
      IF (NPRIM .GT. 0) THEN
         DO IS = NPPPX, NPPPX + (NRELAX - 1)*NPRIM, NPRIM
            SELECT CASE (NPRFN)
               CASE (1)
*                 SPLIT-TYPE PSEUDO-VOIGT FUNCTION
                  IFOR = 6
                  WRITE(PARNAM(IS),'(A)')   'FWHM'
                  WRITE(PARNAM(IS+1),'(A)') 'As'
                  WRITE(PARNAM(IS+2),'(A)') 'eta_L'
                  WRITE(PARNAM(IS+3),'(A)') 'eta_H'
                  SELECT CASE (NMODE)
                     CASE (3)
                        WRITE(PARNAM(IS+4),'(A)') '|Fc|'	 
                     CASE (6)
                        WRITE(PARNAM(IS+4),'(A)') 'const*|Fc|'
                        WRITE(PARNAM(IS+5),'(A)') '2-theta(peak)' 
                  END SELECT
               CASE (2)
*                 MODIFIED SPLIT-TYPE PSEUDO-VOIGT FUCNTION FOR RELAXED 
*                 REFLECTIONS
                  IFOR = 15
                  WRITE(PARNAM(IS),'(A)')   'FWHM(Lorentz)'
                  WRITE(PARNAM(IS+1),'(A)') 'FWHM(Gauss)'
                  WRITE(PARNAM(IS+2),'(A)') 'As'
                  WRITE(PARNAM(IS+3),'(A)') 'eta_L'
                  WRITE(PARNAM(IS+4),'(A)') 'eta_H'
                  SELECT CASE (NMODE)
                     CASE (3)
                        WRITE(PARNAM(IS+5),'(A)') '|Fc|'	 
                     CASE (6)
                        WRITE(PARNAM(IS+5),'(A)') 'const*|Fc|'
                        WRITE(PARNAM(IS+6),'(A)') '2-theta(peak)' 
                  END SELECT
               CASE (3)
*                 SPLIT-TYPE PEARSON VII FUNCTION
                  IFOR = 6
                  WRITE(PARNAM(IS),'(A)')   'FWHM'
                  WRITE(PARNAM(IS+1),'(A)') 'As'
                  WRITE(PARNAM(IS+2),'(A)') 'm_L'
                  WRITE(PARNAM(IS+3),'(A)') 'm_H'
                  SELECT CASE (NMODE)
                     CASE (3)
                        WRITE(PARNAM(IS+4),'(A)') '|Fc|'	 
                     CASE (6)
                        WRITE(PARNAM(IS+4),'(A)') 'const*|Fc|'
                        WRITE(PARNAM(IS+5),'(A)') '2-theta(peak)' 
                  END SELECT
            END SELECT
         END DO
      END IF

      DO L = 1, NPHASE
         PARNAM(LP(NSFX))='Scale factor, s'
         SELECT CASE (NPRFN)
            CASE (0)
               PARNAM(LP(1))='Gaussian FWHM parameter, U'
               PARNAM(LP(2))='Gaussian FWHM parameter, V'
               PARNAM(LP(3))='Gaussian FWHM parameter, W'
               PARNAM(LP(4))=
     &         'Scherrer coefficient for Gaussian broadening, P'
               PARNAM(LP(5))='Lorentzian Scherrer broadening X'
               PARNAM(LP(6))='Anisotropic Scherrer broadening, Xe'
               PARNAM(LP(7))='Strain broadening, Y'
               PARNAM(LP(8))='Anisotropic strain broadening, Ye'
               SELECT CASE (NASYM)
                  CASE (0)
                     PARNAM(LP(9))=
     &               '(Source width)/(detector distance), rs'
                     PARNAM(LP(10))=
     &               '(Detector width)/(detector distance), rd'
                  CASE (1)
                     PARNAM(LP(9))='Asymmetry parameter, As'
                     PARNAM(LP(10))='Dummy'
               END SELECT
               PARNAM(LP(11))='Dummy'
               PARNAM(LP(12))='Dummy'
               PARNAM(LP(13))='Dummy'
               PARNAM(LP(14))='Dummy'
            CASE (1:3)
               PARNAM(LP(1))='FWHM parameter, U'
               PARNAM(LP(2))='FWHM parameter, V'
               PARNAM(LP(3))='FWHM parameter, W'
               PARNAM(LP(4))='Dummy'
               PARNAM(LP(5))='Asymmetry parameter, a0'
               PARNAM(LP(6))='Asymmetry parameter, a1'
               PARNAM(LP(7))='Asymmetry parameter, a2'
               PARNAM(LP(8))='Dummy'
               SELECT CASE (NPRFN)
                  CASE (1,2)
                     PARNAM(LP(9))='Mixing parameter, eta_L0'
                     PARNAM(LP(10))='Mixing parameter, eta_L1'
                     PARNAM(LP(11))='Mixing parameter, eta_H0'
                     PARNAM(LP(12))='Mixing parameter, eta_H1'
                  CASE (3)
                     PARNAM(LP(9))='Decay parameter, m_L0'
                     PARNAM(LP(10))='Decay parameter, m_L1'
                     PARNAM(LP(11))='Decay parameter, m_H0'
                     PARNAM(LP(12))='Decay parameter, m_H1'
               END SELECT
               PARNAM(LP(13))='Anisotropic strain broadening, Ue'
               PARNAM(LP(14))='Anisotropic Scherrer broadening, Pe'
         END SELECT
         IF (NPROR(L) .LE. 2) THEN
            PARNAM(LP(NPO1X))='Preferred-orientation parameter, p1'
            PARNAM(LP(NPO2X))='Preferred-orientation parameter, p2'
         ELSE
            PARNAM(LP(NPO1X))='Preferred-orientation parameter, r'
            PARNAM(LP(NPO2X))='Dummy'
         END IF
         PARNAM(LP(NAX))='Lattice parameter, a'
         PARNAM(LP(NBX))='Lattice parameter, b'
         PARNAM(LP(NCX))='Lattice parameter, c'
         PARNAM(LP(NALPX))='Lattice parameter, alpha'
         PARNAM(LP(NBETX))='Lattice parameter, beta'
         PARNAM(LP(NGAMX))='Lattice parameter, gamma'
         PARNAM(LP(NOTX))='Overall isotropic displacement'//
     &   ' parameter, Q'
         LL=LP(NG1X)
         DO KK = 1, NSITE(L)
            PARNAM(LL)='           Occupation factor, g'
            PARNAM(LL+1)='           Fractional coordinate, x'
            PARNAM(LL+2)='           Fractional coordinate, y'
            PARNAM(LL+3)='           Fractional coordinate, z'
            IF (LISO(KK,L) .EQ. 1) THEN
               PARNAM(LL+4)='           Isotropic displacement'//
     &         ' parameter, B'
               NEXT=LL+5
            ELSE IF (LISO(KK,L) .EQ. 6) THEN
               PARNAM(LL+4)='           Anisotropic'//
     &         ' displacement parameter, beta11'
               PARNAM(LL+5)='           Anisotropic'//
     &         ' displacement parameter, beta22'
               PARNAM(LL+6)='           Anisotropic'//
     &         ' displacement parameter, beta33'
               PARNAM(LL+7)='           Anisotropic'//
     &         ' displacement parameter, beta12'
               PARNAM(LL+8)='           Anisotropic'//
     &         ' displacement parameter, beta13'
               PARNAM(LL+9)='           Anisotropic'//
     &         ' displacement parameter, beta23'
               CALL B2BETA(A(LL+4),L)
               NEXT=LL+10
            ELSE
               NEXT=LL+4
            END IF
            IF (MAGATM(NOAT(KK,L)) .EQ. 1) THEN
               PARNAM(NEXT)='           Magnetic moment, mu'
            END IF
            LL=LL+NPSITE(KK,L)
         END DO
         IF (LMAG(L) .NE. 1) CYCLE
*Rev 1.0c 2000.11.08 Izumi {
         IF (LAUEG(L) .LE. 5) THEN
            PARNAM(KPHE(L)-2)='           phi(a)'
            PARNAM(KPHE(L)-1)='           phi(b)'
            PARNAM(KPHE(L))='           phi(c)'
* }
         ELSE IF (LAUEG(L) .LE. 13) THEN
            PARNAM(KPHE(L))='           phi'
         END IF
      END DO
      WRITE(6,320)
  320 FORMAT(//11X,'Information on profile and structure parameters')

*     PARAMETER INPUT WITH LABELS
*     TOPCOL: LABELS TO BE PRINTED IN FRONT OF PARAMETER NUMBERS
      DO I = 1, NTERMS
         TOPCOL(I)=' '
      END DO
      IPAR=1
      DO I = 1,NLINE
*Rev 1.0t 2001.02.24 Izumi {
C        LENLBL=INDEX(LABEL(I),'/')-1
C        IF (LENLBL .EQ. -1) LENLBL=LENSTR(LABEL(I))
C        IF (LENLBL .GE. 15) THEN
C           TOPCOL(IPAR)=LABEL(I)
C        ELSE
C           TOPCOL(IPAR)(16-LENLBL:15)=LABEL(I)
C        END IF
         DO J = 25, 1, -1
            IF (LABEL(I)(J:J) .NE. ' ') EXIT
         END DO
         TOPCOL(IPAR)(26-J:25) = LABEL(I)

C        TOPCOL(IPAR)(16:16)=':'
         TOPCOL(IPAR)(26:26)=':'
* }

         IF (LABEL(I)(1:3) .EQ. 'PPP' .AND. LABEL(I)(5:5) .EQ. '_'
     &   .AND. INDEX(LABEL(I)(7:),'.') .GT. 0) THEN
*           n_h.k.l IS ATTATCHED TO SIGMA0 FOR A SPECIFIED REFLECTION
            WRITE(PARNAM(IPAR)(IFOR:),'(4A)') 'for phase ',
     &      LABEL(I)(4:4),': ',LABEL(I)(6:)
            DO J = 21, 60
               IF (PARNAM(IPAR)(J:J) .EQ. '.') PARNAM(IPAR)(J:J) = ' '
            END DO
         END IF
         
*        ADD A SITE NAME IN A OCCUPATION FACTOR LINE
*Rev 1.0t 2001.02.24 Izumi {
C        IF (INDEX(PARNAM(IPAR),'Occupation factor, g') .EQ. 12) 
C    &   PARNAM(IPAR)(1:9) = LABEL(I)(1:LENLBL)
         IF (INDEX(PARNAM(IPAR),'Occupation factor, g') .EQ. 12) THEN
            ISLASH = INDEX(LABEL(I),'/')
            IF (ISLASH .GE. 11) THEN
               PARNAM(IPAR)(1:9) = LABEL(I)(1:9)
            ELSE
               PARNAM(IPAR)(1:9) = LABEL(I)(1:ISLASH-1)
            END IF
         END IF
* }
         IPAR = IPAR + LPAR(I)
      END DO
      IF (NMODE .EQ. 1) THEN
         WRITE(6,550) (TOPCOL(J),J,A(J),PARNAM(J),J=1,NTERMS)
*Rev 1.0t 2001.02.24 Izumi {
C 550    FORMAT(/30X,'No.',6X,'A'/(11X,A,I5,2X,1P,G13.6,6X,A))
  550    FORMAT(/40X,'No.',6X,'A'/(11X,A,I5,2X,1P,G13.6,6X,A))
* }
      ELSE
         WRITE(6,600) (TOPCOL(J),J,A(J),ID(J),PARNAM(J),J=1,NTERMS)
*Rev 1.0t 2001.02.24 Izumi {
C 600    FORMAT(/30X,'No.',6X,'A',12X,'ID'/
  600    FORMAT(/40X,'No.',6X,'A',12X,'ID'/
* }
     &   (11X,A,I5,2X,1P,G13.6,5X,I2,6X,A))
      END IF
      END

************************************************************************

      SUBROUTINE B2BETA(A,L)
*     CONVERT AN EQUIVALENT ISOTROPIC DISPLACEMENT PARAMETER INTO BETA'S
      PARAMETER (NPH=8)
      REAL A(*),GG1(3,3),BETA(3,3)
      COMMON /TENS/ G(6,NPH),GG(6,NPH)

      IF (A(2) .NE. 0.0 .OR. A(3) .NE. 0.0 .OR. A(4) .NE. 0.0 .OR.
     &  A(5) .NE. 0.0 .OR. A(6) .NE. 0.0) RETURN
      GG1(1,1)=GG(1,L)
      GG1(2,2)=GG(2,L)
      GG1(3,3)=GG(3,L)
      GG1(2,3)=GG(4,L)
      GG1(1,3)=GG(5,L)
      GG1(1,2)=GG(6,L)
      CALL TOBETA(GG1,A(1),BETA)
      A(1)=BETA(1,1)
      A(2)=BETA(2,2)
      A(3)=BETA(3,3)
      A(4)=BETA(1,2)
      A(5)=BETA(1,3)
      A(6)=BETA(2,3)
      DO 20 I=4,6
         IF (ABS(A(I)) .LT. 4.9999E-7) A(I)=0.0
   20 CONTINUE
      END

************************************************************************

      SUBROUTINE NSCONV(LINE,NSSYM)
C     CONVERT A REFLECTION TYPE AND A REFLECTION CONDITION INTO AN
C     ELEMENT OF ARRAY NSSYM
      CHARACTER LINE*40,REFTYP*3,CONDTN*14,RTITLE(14)*3,CTITLE(28)*14
      INTEGER NSSYM(*)
      DATA RTITLE/
     &'00l','0k0','0kl','h00','h0l','hk0','hkl','0kk','hh0',
     &'hhl','h0h','hkk','hkh','hhh'/
*     DATA CTITLE/
      DATA (CTITLE(J),J=1,28)/
     &'h=2n','k=2n','l=2n','k+l=2n','h+l=2n','h+k=2n',
     &'h,k,l odd/even','k+l=4n','h+l=4n','h+k=4n','2h+l=2n',
     &'2h+l=4n','h+k+l=2n','-h+k+l=3n','h-k+l=3n',
     &'h=4n','k=4n','l=3n','l=4n','l=6n','|h|>=|k|>=|l|',
     &'2h+k=2n','2h+k=4n','h+2k=2n',
     &'h+2k=4n','h=2n,k=2n','k=2n,l=2n','h=2n,l=2n'/

      DO 10 I=40,1,-1
         IF (LINE(I:I) .NE. ' ') GO TO 1
   10 CONTINUE
    1 NEND=I
      REFTYP=' '
      IP=0
      DO 20 I=1,NEND
         IF (LINE(I:I) .EQ. ':') THEN
            GO TO 2
         ELSE IF (LINE(I:I) .NE. ' ') THEN
            IP=IP+1
            IF (IP .GE. 4) CALL JOBEND('Bad reflection type')
            REFTYP(IP:IP)=LINE(I:I)
         END IF
   20 CONTINUE
      CALL JOBEND
     &('Missing a colon in a non-standard symmetry condition')
    2 CONDTN=' '
      IP=1
      DO 30 J=I+1,NEND
         IF (LINE(J:J) .NE. ' ') THEN
            CONDTN(IP:IP)=LINE(J:J)
            IP=IP+1
         END IF
   30 CONTINUE
      DO 40 I=1,14
         IF (REFTYP .EQ. RTITLE(I)) GO TO 5
   40 CONTINUE
      CALL JOBEND('Bad reflection type: '//REFTYP)
    5 DO 50 J=1,28
         IF (CONDTN .EQ. CTITLE(J)) GO TO 6
   50 CONTINUE
      CALL JOBEND('Bad reflection condition: '//CONDTN)
    6 NSSYM(I)=J
      END

************************************************************************

      SUBROUTINE EQPOS(NPHASE,PHNAME,NSITE,ISONUM,ANAME,NOAT,NSYM,IDSYM,
     &  NCENTR,IFAC,SYM,NSSYM)
*     PRINT OUT EQUIVALENT POSITIONS FOR EACH SITE
      PARAMETER(NPH=8,NAP=150,NS=48,NB=7000)
      INTEGER NSITE(*),ISONUM(NAP,NPH),NOAT(NAP,NPH),NSYM(*),
     &  IDSYM(NS,NAP,NPH),NCENTR(*),IFAC(NAP,NPH),NSSYM(14,NPH)
      CHARACTER PHNAME(*)*25,ANAME(*)*5,SYM(3,NS,NPH)*10,CHAR4*4,
     &  INTCHR*4,XYZ(3)*10
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /T/ FF(NB),COEF(NPH)

      DO L = 1, NPHASE
         IF (NMODE .GE. 4 .AND. NSITE(L) .EQ. 0) CYCLE
         WRITE(6,'(//11X,A,A)')
     &   'Equivalent positions for each site in ',PHNAME(L)
         NCON=NSSYM(7,L)
         IF (NCON .EQ. 4) THEN
*           A LATTICE
            WRITE(6,'(13X,A/)') '(0,0,0)+  (0,1/2,1/2)+'
         ELSE IF (NCON .EQ. 5) THEN
*           B LATTICE
            WRITE(6,'(13X,A/)') '(0,0,0)+  (1/2,0,1/2)+'
         ELSE IF (NCON .EQ. 6) THEN
*           C LATTICE
            WRITE(6,'(13X,A/)') '(0,0,0)+  (1/2,1/2,0)+'
         ELSE IF (NCON .EQ. 13) THEN
*           I LATTICE
            WRITE(6,'(13X,A/)') '(0,0,0)+  (1/2,1/2,1/2)+'
         ELSE IF (NCON .EQ. 7) THEN
*           F LATTICE
            WRITE(6,'(13X,A/)') '(0,0,0)+  (0,1/2,1/2)+  (1/2,0,1/2)+'//
     &      '  (1/2,1/2,0)'
         ELSE IF (NCON .EQ. 14) THEN
*           RHOMBOHEDRALLY CENTERED (OBVERSE SETTING)
            WRITE(6,'(13X,A/)') '(0,0,0)+  (2/3,1/3,1/3)+  '//
     &      '(1/3,2/3,2/3)+'
         ELSE IF (NCON .EQ. 15) THEN
*           RHOMBOHEDRALLY CENTERED (REVERSE SETTING)
            WRITE(6,'(13X,A/)') '(0,0,0)+  (1/3,2/3,1/3)+  '//
     &      '(2/3,1/3,2/3)+'
         END IF

         DO J = 1, NSITE(L)
            IF(ISONUM(J,L).NE.0) THEN
               CHAR4=INTCHR(ISONUM(J,L))
               WRITE(6,350) ANAME(NOAT(J,L)),CHAR4
            ELSE
               WRITE(6,350) ANAME(NOAT(J,L))
            ENDIF
            NPOS=0
            DO JJ = 1, NSYM(L)
               IF(IDSYM(JJ,J,L).EQ.1) NPOS=NPOS+1
            END DO
            IF (NCENTR(L) .EQ. 1) THEN
               WRITE(6,370) IFAC(J,L),IFAC(J,L)*NINT(COEF(L))*NPOS/2
            ELSE
               WRITE(6,375) NINT(COEF(L))*NPOS
            ENDIF
            DO JJ = 1, NSYM(L)
               IF(IDSYM(JJ,J,L).EQ.1) WRITE(6,'(21X,3(A10,3X))')
     &         (SYM(LL,JJ,L),LL=1,3)
            END DO
            IF (NCENTR(L) .EQ. 1 .AND. IFAC(J,L) .EQ. 2) THEN
               DO JJ = 1, NSYM(L)
                  IF (IDSYM(JJ,J,L) .EQ. 1) THEN
                     XYZ(1)=SYM(1,JJ,L)
                     XYZ(2)=SYM(2,JJ,L)
                     XYZ(3)=SYM(3,JJ,L)
                     CALL MINUS(XYZ)
                  END IF
               END DO
            END IF
         END DO
      END DO
  350 FORMAT(' ',12X,A5,A4)
  370 FORMAT('+',26X,'Factor =',I2,5X,'neq =',I4)
  375 FORMAT('+',26X,'neq =',I4)
      END

************************************************************************

      SUBROUTINE MINUS(XYZ)
*     PRINT OUT INVERTED COORDINATES IN A CETROSYMMETRIC SPACE GROUP
      CHARACTER*10 XYZ(3),TEMP(3)

      DO 50 I=1,3
         TEMP(I)=' '
         IP=1
         J=1
    1    CONTINUE
            IF (XYZ(I)(J:J) .EQ. ' ') THEN
*              DO NOTHING
            ELSE IF (XYZ(I)(J:J) .EQ. '-' .AND. IP .EQ. 1) THEN
               K=J+1
               DO 20 J=K,10
                  IF (XYZ(I)(J:J) .NE. ' ') GO TO 3
   20          CONTINUE
    3          TEMP(I)(1:1)=XYZ(I)(J:J)
               IP=2
            ELSE IF (XYZ(I)(J:J) .NE. '+' .AND. IP .EQ. 1) THEN
               TEMP(I)(1:1)='-'
               TEMP(I)(2:2)=XYZ(I)(J:J)
               IP=3
            ELSE
               IF (XYZ(I)(J:J) .EQ. '-') THEN
                  TEMP(I)(IP:IP)='+'
               ELSE IF (XYZ(I)(J:J) .EQ. '+') THEN
                  TEMP(I)(IP:IP)='-'
               ELSE
                  TEMP(I)(IP:IP)=XYZ(I)(J:J)
               END IF
               IP=IP+1
            END IF
            J=J+1
         IF (J .LE. 10 .AND. IP .LE. 10) GO TO 1
   50 CONTINUE
      WRITE(6,'(21X,3(A,3X))') (TEMP(J),J=1,3)
      END

************************************************************************

      SUBROUTINE CONSTR(NCNSTR,NTERMS,LABEL,LPAR,NLINE)
*     READ CONSTRAINT(S) AND DIVIDE THEM INTO ONE OR MORE LINES
      INTEGER LPAR(*)
      CHARACTER LABEL(*)*25,LINE*80,LINE2*80,COL(80)*1
      EQUIVALENCE (LINE,COL(1))

      IF (NCNSTR .EQ. 0) THEN
*        NO LINEAR CONSTRAINTS ARE IMPOSED.  SKIP CONSTRAINTS IF ANY.
         DO 25 J=1,10000
            READ(4,'(A)') LINE 
            DO 20 I=1,80
               IF (LINE(I:I) .NE. ' ') GO TO 2
   20       CONTINUE
    2       IF (LINE(I:I) .NE. 'A' .AND. LINE(I:I) .NE. '}') THEN
               BACKSPACE 4
               RETURN
            END IF
   25    CONTINUE
      END IF

*     LINEAR CONSTRAINTS ARE IMPOSED.  READ CONSTRAINTS.
    4 WRITE(6,'(//11X,A)') 'Linear constraints'
      ICNSTR=0
      DO 80 ILINE=1,1000
         READ(4,'(A)') LINE2
         IF (LINE2(1:1) .EQ. '}') THEN
            CALL JOBEND('The number of linear constraints is smaller'//
     &      ' than that expected from the ID values of parameters')
         END IF
         WRITE(6,'(13X,A)') LINE2
         IPOINT=1
         DO 40 J=1,1000
            ISEM=INDEX(LINE2,';')
            IF (ISEM .EQ. 0) THEN
               LINE=LINE2(IPOINT:80)
               IF (LINE .EQ. ' ') GO TO 80
            ELSE
               LINE2(ISEM:ISEM)=' '
               LINE=LINE2(IPOINT:ISEM-1)
               IPOINT=ISEM+1
            END IF
            CALL LINCON(COL,NTERMS,LABEL,LPAR,NLINE)
            ICNSTR=ICNSTR+1
            IF (ICNSTR .EQ. NCNSTR .AND. ISEM .EQ. 0) THEN
               GO TO 3
            ELSE IF (ICNSTR .EQ. NCNSTR .AND. ISEM .NE. 0) THEN
               CALL JOBEND('Too large number of linear constraints')
            ELSE IF (ISEM .EQ. 0 .OR. IPOINT .GE. 81) THEN
               GO TO 80
            END IF
   40    CONTINUE
   80 CONTINUE
    3 READ(4,'(A)') LINE2

C*    END OF CONSTRAINTS: '}'
      IF (LINE2(1:1) .NE. '}') THEN
         CALL JOBEND('Missing ''}'' to indicate the end of linear '//
     &   'constraints')
      END IF
      END

************************************************************************

      SUBROUTINE LINCON(COL1,NTERMS,LABEL,LPAR,NLINE)
C     DECODE A LINE IN WHICH A LINEAR CONSTRAINT IS DESCRIBED
      PARAMETER (NT=999,NCS=400)
      INTEGER LPAR(*)
      CHARACTER LABEL(*)*25
      CHARACTER COL1(80)*1,LINE*80,COL(80)*1,FORM1*11,CH(39)*1
      CHARACTER BUFF*120
      EQUIVALENCE (LINE,COL(1))
      COMMON /CNSTR/ C1(NCS),C2(10,NCS),IN(10,NCS),NCNPAR(NCS),
     &  NACNS(NCS),NCNSTR
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      DATA CH/'0','1','2','3','4','5','6','7','8','9','.','-','+',
     &'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',
     &'P','Q','R','S','T','U','V','W','X','Y','Z'/
C
C     COPY COL1 INTO COL, REMOVING ' ' IN A LINE
      IP=1
      DO 12 I=1,80
         IF (COL1(I) .NE. ' ') THEN
            COL(IP)=COL1(I)
            IP=IP+1
         END IF
   12 CONTINUE
      LENGTH=IP-1
      IF (COL(1) .NE. 'A' .OR. COL(2) .NE. '(') THEN
         CALL JOBEND('Invalid linear constraint')
      END IF
C
C     LEFT-HAND SIDE OF THE EQUATION
      DO 20 ICOL=1,LENGTH
         IF (COL(ICOL) .EQ. '(') THEN
            IB=ICOL
         ELSE IF (COL(ICOL) .EQ. ')') THEN
            IE=ICOL
            GO TO 3
         END IF
   20 CONTINUE
    3 NLEFT=NUMPAR(NTERMS,LINE(IB+1:IE-1),LABEL,LPAR,NLINE)
      
C     ICON: CONSTRAINT NUMBER
      ICON=0
      DO 25 J=1,NTERMS
         IF (ID(J) .EQ. 2) THEN
            ICON=ICON+1
            IF (J .EQ. NLEFT) GO TO 4
         END IF
   25 CONTINUE
      BUFF = 'Check the left hand side of the linear constraint: '//
     &LINE(IB-1:IE)
      CALL JOBEND(BUFF)
    4 NACNS(ICON)=NLEFT
      C1(ICON)=0.0
      DO 27 I=1,10
         C2(I,ICON)=1.0
   27 CONTINUE
C
C     LOUT=0: THE TOP OF A NUMBER HAS BEEN ENCOUNTERED
C     LOUT=1: A NUMBER HAS NOT YET BEEN ENCOUNTERED OR HAS BEEN
C             STORED IN ARRAY C1, C2, OR IN
      LOUT=1
      IA=1
      IB=IE+2

C     RIGHT-HAND SIDE OF THE EQUATION
C     LINEAR CONSTRAINT: A(NACNS(J))=C1(J)+C2(1,J)*A(IN(1,J))+C2(2,J)*
C                        A(IN(2,J))+ .....
C     PARAMETERS CONTAINED IN THE ABOVE FORMULA MUST BE STRUCTURE-FACTOR
C     PARAMETERS.
C     THE ABOVE EQUALITY CONSTRAINT IS USED IN SUBROUTINE SETPAR
C     NCNSTR: NUMBER OF CONSTRAINTS
C     NCNPAR: NUMBER OF PARAMETERS CONTAINED IN THE RIGHT SIDE OF THE
C             EQUALITY CONSTRAINT
      DO 70 ICOL=IB,LENGTH
         IF (COL(ICOL) .EQ. ')') THEN
            IE=ICOL
            IN(IA,ICON)=NUMPAR(NTERMS,LINE(IB:IE-1),LABEL,LPAR,NLINE)
            IA=IA+1
            LOUT=1
         ELSE IF (COL(ICOL) .EQ. '*') THEN
            IE=ICOL
            WRITE(FORM1,110) IB-1,IE-IB
  110       FORMAT('(',I2,'X,F',I2,'.0)')
            READ(LINE,FORM1) C2(IA,ICON)
            LOUT=1
         ELSE IF ((COL(ICOL) .EQ. '+' .OR. (COL(IB-1) .NE. '(' .AND.
     &   COL(ICOL) .EQ. '-')) .AND. COL(ICOL-1) .NE. 'E') THEN
            IF (LOUT .EQ. 1) THEN
               IB=ICOL
               LOUT=0
            ELSE
               IE=ICOL
               WRITE(FORM1,110) IB-1,IE-IB
               READ(LINE,FORM1) C1(ICON)
               IB=ICOL
            END IF
         ELSE IF (LINE(ICOL:ICOL+1) .EQ. 'A(' .AND.
     &   LOUT .EQ. 0 .AND. COL(IB) .EQ. '+') THEN
            C2(IA,ICON)=+1.0
            LOUT=1
         ELSE IF (LINE(ICOL:ICOL+1) .EQ. 'A(' .AND.
     &   LOUT .EQ. 0 .AND. COL(IB) .EQ. '-') THEN
            C2(IA,ICON)=-1.0
            LOUT=1
         ELSE IF (ICOL .EQ. LENGTH) THEN
            IF (LOUT .EQ. 0) THEN
               WRITE(FORM1,110) IB-1,ICOL-IB+1
               READ(LINE,FORM1) C1(ICON)
            END IF
            GO TO 9
         ELSE IF (LINE(ICOL:ICOL+1) .EQ. 'A(') THEN
*           DO NOTHING
         ELSE
            DO 30 J=1,39
               IF (COL(ICOL) .EQ. CH(J)) THEN
                  LP=J
                  GO TO 5
               END IF
   30       CONTINUE
            LP=0
    5       IF (LP .GT. 0 .AND. LOUT .EQ. 1) THEN
               IB=ICOL
               LOUT=0
            END IF
         END IF
   70 CONTINUE
    9 NCNPAR(ICON)=IA-1
      END

************************************************************************

      SUBROUTINE APRNTD(A,NTERMS)
C     PRINTED VALUES OF ARRAY A
      PARAMETER (NT=999,NPH=8,NAP=150)
      DIMENSION A(*)
      COMMON /CC/ APR(NT)
      COMMON /H/ NQ(NPH),NCENTR(NPH),LAUEG(NPH)
      COMMON /LTM/ LISO(NAP,NPH),LMAG(NPH)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /TENS/ G(6,NPH),GG(6,NPH)
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      DATA RD,RD2,RDSQR/0.0174532925,57.29578,0.0003046174/
      LP(J)=KPHB(L)+J
      Z(J1,J2,J3)=ACOS(G(J1,L)/(APR(LP(J2))*APR(LP(J3))))*RD2
C
C     RAD ==> DEG
      DO J = 1, NBCKGR - 1
         APR(J)=A(J)*RD2
      END DO
      DO J = NBCKGR, NTERMS
         APR(J)=A(J)
      END DO
      IF (NPRIM .GT. 0) THEN
         DO J = NPPPX, NPPPX + (NRELAX - 1)*NPRIM, NPRIM
*           RAD ==> DEG (FWHM FOR A RELAXED REFLECTION)
            APR(J) = A(J)*RD2
            IF (NPRFN .EQ. 2) APR(J+1) = A(J+1)*RD2
            IF (NMODE .EQ. 6) THEN
*              2-theta(peak) in degree
               SELECT CASE (NPRFN)
                  CASE (1,3)
                     APR(J+5) = A(J+5)*RD2
                  CASE (2)
                     APR(J+6) = A(J+6)*RD2
               END SELECT
            END IF
         END DO
      END IF

***   IF A NEW PROFILE FUNCTION IS IMPLEMENTED, THE FOLLOWING PART
***   SHOULD BE MODIFIED.
      DO L = 1, NPHASE
*        RAD**2 ==> DEG**2 (U, V, W, AND P)
         DO J=LP(1),LP(4)
            APR(J)=A(J)/RDSQR
         END DO
         IF (NPRFN .EQ. 0) THEN
*           RAD ==> DEG (X, Xe, Y, Ye, and As)
            DO J = LP(5), LP(9)
               IF (NASYM .EQ. 0 .AND. J .EQ. LP(9)) EXIT
               APR(J) = A(J)*RD2
            END DO
         END IF
         
C        LATTICE CONSTANTS
         APR(LP(NAX))=SQRT(G(1,L))
         APR(LP(NBX))=SQRT(G(2,L))
         APR(LP(NCX))=SQRT(G(3,L))
         APR(LP(NALPX))=Z(4,NBX,NCX)
         APR(LP(NBETX))=Z(5,NCX,NAX)
         APR(LP(NGAMX))=Z(6,NAX,NBX)
*        GIVE STRICT VALUES OF ALPHA, BETA, OR GAMMA
         LG=LAUEG(L)
         IF ((LG .GE. 5 .AND. LG .LE. 7) .OR. LG .GE. 14) THEN
            APR(LP(NALPX))=90.0
            APR(LP(NBETX))=90.0
            APR(LP(NGAMX))=90.0
         ELSE IF (LG .EQ. 2) THEN
            APR(LP(NBETX))=90.0
            APR(LP(NGAMX))=90.0
         ELSE IF (LG .EQ. 3) THEN
            APR(LP(NALPX))=90.0
            APR(LP(NGAMX))=90.0
         ELSE IF (LG .EQ. 4) THEN
            APR(LP(NALPX))=90.0
            APR(LP(NBETX))=90.0
         ELSE IF (LG .EQ. 9 .OR. LG .EQ. 11 .OR. LG .EQ. 12 .OR.
     &   LG .EQ. 13) THEN
            APR(LP(NALPX))=90.0
            APR(LP(NBETX))=90.0
            APR(LP(NGAMX))=120.0
         END IF
         IF (LMAG(L) .NE. 1) CYCLE
C        RAD ==> DEG (PHI)
*Rev 1.0c 2000.11.08 Izumi
         IF (LAUEG(L) .LE. 5) THEN
            APR(KPHE(L)-2)=A(KPHE(L)-2)*RD2
            APR(KPHE(L)-1)=A(KPHE(L)-1)*RD2
         END IF
         APR(KPHE(L))=A(KPHE(L))*RD2
      END DO
      END

************************************************************************

      SUBROUTINE CCTENS(A)
C     CONVERT THE DIRECT UNIT CELL CONSTANTS INTO THE METRIC TENSORS
      PARAMETER (NPH=8,NAP=150)
      DIMENSION A(*)
      COMMON /H/ NQ(NPH),NCENTR(NPH),LAUEG(NPH)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /TENS/ G(6,NPH),GG(6,NPH)
      DATA RD/0.0174532925/
      LP(J)=KPHB(L)+J
      P(J)=A(LP(J))

      DO 10 L=1,NPHASE
         G(1,L)=P(NAX)**2
         G(2,L)=P(NBX)**2
         G(3,L)=P(NCX)**2
         G(4,L)=P(NBX)*P(NCX)*COS(P(NALPX)*RD)
         G(5,L)=P(NCX)*P(NAX)*COS(P(NBETX)*RD)
         G(6,L)=P(NAX)*P(NBX)*COS(P(NGAMX)*RD)
         CALL INVERS(GG(1,L),G(1,L),DET)
   10 CONTINUE
      END

************************************************************************

      SUBROUTINE SYMOP(A)
C     EXAMINE WHETHER OR NOT EACH SYMMETRY OPERATION IS NECESSARY
      INTEGER R,RX,RY,RZ
      LOGICAL OVR
      PARAMETER (NB=7000,NS=48,NAP=150,NT=999,NPH=8)
      DIMENSION A(*),XYZ(3,NS)
      COMMON /EXC/ LSPSYM(NPH),NSSYM(14,NPH)
      COMMON /H/ NQ(NPH),NCENTR(NPH),LAUEG(NPH)
      COMMON /HALF/ IFAC(NAP,NPH)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /RT/ TV(3,NS,NPH),HT(NS,NB),NSYM(NPH),NSITE(NPH),
     &  IDSYM(NS,NAP,NPH),R(3,3,NS,NPH),RX(NS,NB),RY(NS,NB),RZ(NS,NB)
      DELTA(I)=ABS(XYZ(I,L1)-XYZ(I,L2))-0.001
*     IF AN INVERTED POSITION OVERLAPS WITH ANOTHER EQUIVALENT POSITION,
*     OVR(I) BECOMES .TRUE.
*     SPECIAL CARE IS TAKEN FOR 1/2 BECAUSE 1/2 IS EQUIVALENT TO -1/2.
      OVR(I) = ABS(XYZ(I,L1)+XYZ(I,L2)) .LT. 0.001 .OR.
     &         (ABS(XYZ(I,L1)-0.5) .LT. 0.0001 .AND.
     &         ABS(XYZ(I,L2)-0.5) .LT. 0.0001)

      DO L = 1, NPHASE
         IF (NMODE .EQ. 4 .AND. L .EQ. 1) CYCLE
         IF (NMODE .EQ. 6) CYCLE
         NCON=NSSYM(7,L)
C        IDSYM=0: THIS SYMMETRY OPERATION IS SKIPPED
C        IDSYM=1: THIS SYMMETRY OPERATION IS NOT SKIPPED
C        A(IX):   X
C        A(IX+1): Y
C        A(IX+2): Z
C        EACH COORDINATE, THAT IS, XYZ(1,JJ), XYZ(2,JJ), XYZ(3,JJ) ARE
C        NORMALIZED AS -1/2<XYZ<1/2
         IX=KPHB(L)+NG1X+1
         DO J = 1, NSITE(L)
            IFAC(J,L)=2
            DO JJ = 1, NSYM(L)
               IDSYM(JJ,J,L)=1
               DO K = 1, 3
                  XYZ(K,JJ) = FLOAT(R(K,1,JJ,L))*A(IX) + 
     &            FLOAT(R(K,2,JJ,L))*A(IX+1) +
     &            FLOAT(R(K,3,JJ,L))*A(IX+2) + TV(K,JJ,L)
                  XYZ(K,JJ) = XNORM(XYZ(K,JJ))
               END DO
            END DO
C           WHEN TWO POINTS HAVE APPROXIMATELY THE SAME FRACTIONAL
C           COORDINATES, THE POINTS WITH HIGHER SYMMETRY # IS SKIPPED.
            DO L1 = 1, NSYM(L)
               DO L2 = L1+1, NSYM(L)
                  IF(IDSYM(L2,J,L).EQ.0) CYCLE
                  IF(DELTA(1).LT.0.0 .AND.
     &            DELTA(2).LT.0.0 .AND. DELTA(3).LT.0.0) THEN
                     IDSYM(L2,J,L)=0
                  ELSE IF (NCENTR(L) .EQ. 1 .AND. OVR(1) .AND.
     &            OVR(2) .AND. OVR(3)) THEN
                     IDSYM(L2,J,L)=0
                  ELSE IF(NCON .EQ. 0 .OR. NCON .EQ. 21) THEN
C                    DO NOTHING FOR THE PRIMITIVE LATTICE
                  ELSE
C                    COMPLEX LATTICE
                     IDSYM(L2,J,L)=ICOINC(XYZ(1,L1),XYZ(2,L1),XYZ(3,L1),
     &               XYZ(1,L2),XYZ(2,L2),XYZ(3,L2),NCON)
*                    CHECK WHETHER OR NOT INVERTED AND TRANSLATED
*                    POSITIONS OVERLAP WITH (X1,Y1,Z1)
                     IF (NCENTR(L) .EQ. 1 .AND. IDSYM(L2,J,L) .EQ. 1)
     &               THEN
                        IDSYM(L2,J,L)=ICOINC(XYZ(1,L1),XYZ(2,L1),
     &                  XYZ(3,L1),-XYZ(1,L2),-XYZ(2,L2),-XYZ(3,L2),NCON)
                     END IF
                  ENDIF
               END DO
            END DO
            IF(NCENTR(L).EQ.1) CALL INVRSN(XYZ,L,J,NCON)
C           NEXT SITES
            IX=IX+NPSITE(J,L)
         END DO
      END DO
      END

************************************************************************

      INTEGER FUNCTION ICOINC(X1,Y1,Z1,X2,Y2,Z2,NCON)
C     ICOINC=0: TRANSLATED POSITIONS OVERLAPS WITH (X1,Y1,Z1)
      ICOINC=1
      IF(NCON .EQ. 13) THEN
C        BODY CENTERED
         ICOINC=IOVRLP(X1,Y1,Z1,X2+0.5,Y2+0.5,Z2+0.5)
      ELSEIF(NCON .EQ. 4) THEN
C        A-FACE CENTERED
         ICOINC=IOVRLP(X1,Y1,Z1,X2,Y2+0.5,Z2+0.5)
      ELSEIF(NCON .EQ. 5) THEN
C        B-FACE CENTERED
         ICOINC=IOVRLP(X1,Y1,Z1,X2+0.5,Y2,Z2+0.5)
      ELSEIF(NCON .EQ. 6) THEN
C        C-FACE CENTERED
         ICOINC=IOVRLP(X1,Y1,Z1,X2+0.5,Y2+0.5,Z2)
      ELSEIF(NCON .EQ. 7) THEN
C        ALL-FACE CENTERED
         ICOINC=IOVRLP(X1,Y1,Z1,X2,Y2+0.5,Z2+0.5)
         IF(ICOINC.EQ.0) RETURN
         ICOINC=IOVRLP(X1,Y1,Z1,X2+0.5,Y2,Z2+0.5)
         IF(ICOINC.EQ.0) RETURN
         ICOINC=IOVRLP(X1,Y1,Z1,X2+0.5,Y2+0.5,Z2)
      ELSEIF(NCON .EQ. 14 .OR. NCON .EQ. 15) THEN
C        RHOMBOHEDRALLY CENTERED (HEXAGONAL AXES)
         ICOINC=IOVRLP(X1,Y1,Z1,X2+0.6666667,Y2+0.3333333,Z2+0.3333333)
         IF(ICOINC.EQ.0) RETURN
         ICOINC=IOVRLP(X1,Y1,Z1,X2+0.3333333,Y2+0.6666667,Z2+0.6666667)
      ENDIF
      END

************************************************************************

      INTEGER FUNCTION IOVRLP(X1,Y1,Z1,X,Y,Z)
C     IOVRLP=0: THE TWO POINTS OVERLAP WITH EACH OTHER
      DELTA(A,B)=ABS(A-XNORM(B))-0.001
      IOVRLP=1
      IF(DELTA(X1,X).LT.0.0 .AND. DELTA(Y1,Y).LT.0.0 .AND.
     &  DELTA(Z1,Z).LT.0.0) IOVRLP=0
      END

************************************************************************

      FUNCTION XNORM(X)
*     NORMALIZE X BETWEEN -0.5 AND 0.5

      XNORM = X
    1 CONTINUE
         IF (XNORM .LT. -0.49999) THEN
            XNORM = XNORM + 1.0
         ELSE IF (XNORM .GT. 0.50001) THEN
            XNORM = XNORM - 1.0
         ELSE
            RETURN
         ENDIF
      GO TO 1
      END

************************************************************************

      SUBROUTINE INVRSN(XYZ,L,J,NCON)
C     EXAMINE WHETHER OR NOT (-X,-Y,-Z) COINCIDES WITH OTHER EQUIVALENT
C     POSITIONS
      INTEGER R,RX,RY,RZ
      PARAMETER (NB=7000,NS=48,NAP=150,NT=999,NPH=8)
      DIMENSION XYZ(3,NS)
      COMMON /HALF/ IFAC(NAP,NPH)
      COMMON /RT/ TV(3,NS,NPH),HT(NS,NB),NSYM(NPH),NSITE(NPH),
     &  IDSYM(NS,NAP,NPH),R(3,3,NS,NPH),RX(NS,NB),RY(NS,NB),RZ(NS,NB)
      DELTA(I)=ABS(XYZ(I,K)-XNORM(-XYZ(I,1)))-0.001

C     IFAC=1: THE OCCUPATION FACTOR IS MULTIPLIED BY 0.5 IF X,Y,Z
C             COINCIDES WITH -X,-Y,-Z.
      DO 20 K=1,NSYM(L)
         IF(IDSYM(K,J,L).EQ.0) GO TO 20
         IF(DELTA(1).LT.0.0 .AND.
     &   DELTA(2).LT.0.0 .AND. DELTA(3).LT.0.0) THEN
            IFAC(J,L)=1
            RETURN
         ELSE IF (NCON .EQ. 0 .OR. NCON .EQ. 21) THEN
C           DO NOTHING FOR THE PRIMITIVE LATTICE
         ELSE IF (ICOINC(XYZ(1,1),XYZ(2,1),XYZ(3,1),-XYZ(1,K),-XYZ(2,K),
     &   -XYZ(3,K),NCON).EQ.0) THEN
            IFAC(J,L)=1
            RETURN
         ENDIF
   20 CONTINUE
      END

************************************************************************

      SUBROUTINE RDBKGD(NSTEP,NEXC)
*     READ BACKGROUND INTENSITIES
      PARAMETER (NP=80000)
*     RELEASE THIS ARRAY IF FORTRAN 90 IS AVAILABLE
      REAL DEG(NP)
      CHARACTER FILE8*50
      INTEGER NTOTAL,NEXC,NSKIP,NSTEP
      COMMON /C/ BG(NP),BGINC(NP),DEGNOR(NP)
      COMMON /EXCREG/ DEGEXC(2,50),NSKIP

      CALL OPENFILE(3,FILE8)
      REWIND 8
      DO J = 1, NP
         READ(8,*,END=9) DEG(J),BG(J)
      END DO
    9 NTOTAL = J - 1

*     DELETE INTENSITY DATA TO BE SKIPPED
      NPOINTS = 0
      DO LL = 1, NTOTAL
         X1 = DEG(LL)
         IF (NEXC .EQ. 1) THEN
            DO J = 1, NSKIP
               IF (X1 .GT. DEGEXC(1,J) - 0.0005 .AND. X1 .LT.
     &         DEGEXC(2,J) + 0.0005) CYCLE
            END DO
         END IF
         NPOINTS = NPOINTS + 1
         BG(NPOINTS) = BG(LL)
         BGINC(NPOINTS) = BG(LL)
      END DO
      IF (NSTEP .NE. NPOINTS) CALL JOBEND('The total number of backgroun
     &ds is not equal to that of diffraction intensities')
      END

************************************************************************

      SUBROUTINE SMINT(X,Y,DEG,XINT,NBG,NSTEP)
C     SMOOTHING USING FIVE POINTS
      PARAMETER (NP=80000)
      COMMON /C/ BG(NP),BGINC(NP),DEGNOR(NP)
      DIMENSION X(*),Y(*),DEG(*),XINT(*)
      REAL EOPT(1)
      INTEGER IOPT(3)
      DATA IOPT /0, 1, 0/
      DATA NDEG, LUP /3, 0/

      WRITE(6,100) (X(J), Y(J), J = 1, NBG)
  100 FORMAT(//11X,'Two-theta versus background (zero for smoothing)'/
     &  (13X,F7.3,F10.1))
      DO J = 1, NBG
         X(J) = X(J)*0.01745329
         IF (Y(J) .GT. 0.0) CYCLE
         DO MM = 3, NSTEP-2
            IF (ABS(DEG(MM)-X(J)) .LT. 2.0E-5) GO TO 1
         END DO
         WRITE(6,110) X(J)/0.01745329
         CALL JOBEND(' ')
    1    Y(J) = (-3.0*XINT(MM-2)+12.0*XINT(MM-1)+17.0*XINT(MM)+
     &          12.0*XINT(MM+1)-3.0*XINT(MM+2))/35.0
      END DO
  110 FORMAT(//11X,'A bad background point or zero background in an excl
     &uded region: ',F7.3)

C     LAGRANGE INTERPOLATION
      DO J = 1, NSTEP
         CALL SILUP(DEG(J),G,NBG,X,Y,NDEG,LUP,IOPT,EOPT)
         BG(J) = G
      END DO
      END

************************************************************************

      CHARACTER*4 FUNCTION INTCHR(J)
C     CONVERT AN INTEGER INTO FOUR CHARACTERS

      INTCHR='    '
      IF(J.GE.0 .AND. J.LE.9) THEN
         WRITE(INTCHR,'(A,I1,A)') '(',J,')'
      ELSEIF(J.GE.10) THEN
         WRITE(INTCHR,'(A,I2,A)') '(',J,')'
      ENDIF
      END

************************************************************************

      SUBROUTINE TENSOR(A)
C     STORE THE METRIC TENSORS
      PARAMETER (NB=7000,NPH=8,NAP=150)
      DIMENSION A(*)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /TENS/ G(6,NPH),GG(6,NPH)
      P(J)=A(KPHB(L)+J)

      DO 10 L=1,NPHASE
         GG(1,L)=P(NAX)
         GG(2,L)=P(NBX)
         GG(3,L)=P(NCX)
         GG(4,L)=P(NALPX)
         GG(5,L)=P(NBETX)
         GG(6,L)=P(NGAMX)
         CALL INVERS(G(1,L),GG(1,L),DET)
   10 CONTINUE
      END

************************************************************************

      SUBROUTINE PRUCC
C     PRINT OUT RECIPROCAL UNIT CELL CONSTANTS
      CHARACTER PARNAM*60,PHNAME*25
      PARAMETER (NT=999,NPH=8,NAP=150)
      COMMON /PAR/ PARNAM(NT),PHNAME(NPH)
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /TENS/ G(6,NPH),GG(6,NPH)
      Z2(J,X,Y)=ACOS(GG(J,L)/(X*Y))*57.29578

      DO 10 L=1,NPHASE
         CALL INVERS(GG(1,L),G(1,L),DET)
         A2=SQRT(GG(1,L))
         B2=SQRT(GG(2,L))
         C2=SQRT(GG(3,L))
         ALPHA2=Z2(4,B2,C2)
         BETA2=Z2(5,C2,A2)
         GAMMA2=Z2(6,B2,A2)
         WRITE(6,100) PHNAME(L),A2,B2,C2,ALPHA2,BETA2,GAMMA2,
     &   1.0/SQRT(DET),SQRT(DET)
   10 CONTINUE
  100 FORMAT(//11X,'Dimensions of the reciprocal cell for '
     &,A/13X,'a* =',F8.5,5X,'b* =',F8.5,5X,'c* =',F8.5,5X,
     &'alpha* = ',F7.3,5X,'beta* = ',F7.3,5X,'gamma* = ',F7.3/
     &13X,'V* =',1P,G14.6,5X,'(V =',1P,G14.6,')')
      END

************************************************************************

      SUBROUTINE DECODE(COORD,SYM,N,L)
C     CONVERT THREE CHARACTER DATA INTO A ROTATION MATRIX (R) AND
C     A TRANSLATION VECTOR (T)
      PARAMETER (NB=7000,NS=48,NAP=150,NPH=8)
      INTEGER R,RX,RY,RZ
      CHARACTER CHAR*1,SYM(3,NS)*10,COORD*20,XYZ(3)*1
      DIMENSION IX(3)
      COMMON /RT/ T(3,NS,NPH),HT(NS,NB),NSYM(NPH),NSITE(NPH),
     &  IDSYM(NS,NAP,NPH),R(3,3,NS,NPH),RX(NS,NB),RY(NS,NB),RZ(NS,NB)
      DATA (XYZ(J),J=1,3)/'x','y','z'/

      I1=INDEX(COORD,',')
      DO 10 J=I1+2,19
         READ(COORD(J:J),'(A1)') CHAR
         IF(CHAR.EQ.',') THEN
            I2=J
            GO TO 1
         ENDIF
   10 CONTINUE
    1 SYM(1,N)=COORD(1:I1-1)
      SYM(2,N)=COORD(I1+1:I2-1)
      SYM(3,N)=COORD(I2+1:20)
C
      DO 30 J=1,3
         T(J,N,L)=0.0
         DO 20 K=1,3
            R(J,K,N,L)=0
   20    CONTINUE
C
C        ROTATION MATRIX
         DO 25 JJ=1,3
            IX(JJ)=INDEX(SYM(J,N),XYZ(JJ))
            IF(IX(JJ).GT.0) THEN
               IF(IX(JJ).EQ.1) THEN
                  R(J,JJ,N,L)=1
               ELSE
                  IF(SYM(J,N)(IX(JJ)-1:IX(JJ)-1).EQ.'+') THEN
                     R(J,JJ,N,L)=1
                  ELSE
                     R(J,JJ,N,L)=-1
                  ENDIF
               ENDIF
            ENDIF
   25    CONTINUE
C
C        TRANSLATION VECTOR
         IT=INDEX(SYM(J,N),'/')
         IF(IT.GT.0) THEN
            READ(SYM(J,N)(IT-1:IT-1),'(I1)') N1
            READ(SYM(J,N)(IT+1:IT+1),'(I1)') N2
            T(J,N,L)=FLOAT(N1)/FLOAT(N2)
            IF(IT.NE.2 .AND. SYM(J,N)(IT-2:IT-2).EQ.'-')
     &      T(J,N,L)=-T(J,N,L)
         ENDIF
   30 CONTINUE
      END

************************************************************************

      SUBROUTINE CONDIR(NTERMS,NSTEP,A,DEG,XINT,YFIT)
C     NONLINEAR LEAST-SQUARES FITTING BY POWELL'S CONJUGATE DIRECTION
C     METHOD
      CHARACTER PARNAM*60,PHNAME*25
      PARAMETER (NB=7000,NT=999,NCS=400,NR=400,NPH=8)
      DIMENSION A(*),DEG(*),XINT(*),YFIT(*),B(NT)
      COMMON /CNSTR/ C1(NCS),C2(10,NCS),IN(10,NCS),NCNPAR(NCS),
     &  NACNS(NCS),NCNSTR
      COMMON /FL0/ NC,TK,CONV,NCONV,FINC,SUMPEN,SUMP1,FACTOR
      COMMON /ONE/ Y(NR),S(NR),FX,FY
      COMMON /PAR/ PARNAM(NT),PHNAME(NPH)
      COMMON /POW/ MITER,STEP,ACC
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /PRLV/ NPRINT
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /RLS/ NRANGE
      COMMON /THREE/ N,NFUNCT,NDRV,ITER,INDIC,IPRINT
      COMMON /TWO/ DIRECT(NR,NR),BEFORE(NR),FIRST(NR)
      EXTERNAL PACK,UNPACK
      
      WRITE(6,100)
  100 FORMAT(///7X,'*** Nonlinear least-squares fitting by the conjugate
     & direction method (NLESQ = 2) ***'/)
      WRITE(6,130) MITER,STEP,ACC
  130 FORMAT(/11X,'MITER =',I2,5X,'STEP = ',F5.3,5X,
     &'ACC = ',1P,E9.3)
      IF(NC.NE.0) WRITE(6,150) NC,TK
  150 FORMAT(11X,'NCONS =',I3,5X,'T(0) =',1P,G13.6/)
      IF (MITER .EQ. 0) RETURN
      LCOUNT=0
      DO 10 J=1,NTERMS
         IF(ID(J).EQ.1 .AND. FIRST(IR(J)).NE.0.0) THEN
            IF(LCOUNT.EQ.0) THEN
               WRITE(6,'(/11X,A)') 'Step(s) specified by user'
               LCOUNT=1
            ENDIF
            WRITE(6,'(13X,A,I3,A,F7.4)') 'STEP(',J,') =',FIRST(IR(J))
         ENDIF
   10 CONTINUE
      CALL PACK(A,B,NTERMS)
      N=NRFN
      NFUNCT=0
      NDRV=0
      INDIC=2
      IPRINT=NPRINT
      CALL MINI(YFIT,B,XINT,DEG,NTERMS,NSTEP,A)
      CALL UNPACK(A,B,NTERMS)
      CALL UPDATE(A,IDPOS,1)
      END

************************************************************************

      SUBROUTINE MARGAU(NCYCL,NTERMS,NSTEP,NPTS,A,DEG,XINT,YFIT,
     &  SIGMAA,NESD)
C     NONLINEAR LEAST-SQUARES FITTING BY THE GAUSS-NEWTON METHOD OR
C     THE MODIFIED MARQUARDT METHOD
      CHARACTER PARNAM*60,PHNAME*25,INTCHR*4,CHAR4*4
      DOUBLE PRECISION G
      PARAMETER (NA=15,NB=7000,NT=999,NR=400,NPH=8,NAP=150,
     &  NOD=NR*(NR+1)/2)
      DIMENSION A(*),DEG(*),XINT(*),SIGMAA(NT),YFIT(*),
     &  AMIN(NT),ASAVE(NT),DELTA(NT)
      COMMON /A/ IH(NB),IK(NB),IL(NB),L12(NB),NOPH(NB)
      COMMON /B/ U(NB),NREFL
      COMMON /CC/ APR(NT)

      COMMON /FL0/ NC,TK,CONV,NCONV,FINC,SUMPEN,SUMP1,FACTOR
      COMMON /FRDA/ NFR,NMEM,MEED(15),EPSD,SCIO,TSCAT,TSCAT1,TSCAT2,
     &  UCLAG,NDA,NFPX(3)
      COMMON /G/ I
      COMMON /NONLIN/ NLESQ,NAUTO,IDCYCL,NCHNG,NPAR(50),IPAR(NR,50)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /OLD/ AOLD(NT),AOUT(NT)
      COMMON /ORDER/ IPOINT(NB)
      COMMON /PAR/ PARNAM(NT),PHNAME(NPH)
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /PRLV/ NPRINT
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /RLS/ NRANGE
      COMMON /S/ F(NA,NB),NATOM
      COMMON /TWO/ G(NOD),AD(NR)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
      COMMON /S2/ F0(NA,NB)
*     SIGSCAL is used in REFLCTN.
      COMMON /SIGMAS2/ SIGSCALS(NPH),SCALEAS(NPH)
      COMMON /SAVEYP/ SYPEAK(NB),SAVEA(NT),NTOTPAR
      DATA DELTA/NT*0.0/

      IF (IDCYCL .EQ. 0) THEN
         WRITE(6,100)
         WRITE(6,'(11X,A,I1)') 'NESD = ',NESD
      ELSE IF (NLESQ .EQ. 0) THEN
         WRITE(6,110)
      ELSE IF (NLESQ .EQ. 1) THEN
         WRITE(6,120)
      END IF
  100 FORMAT(1H1/7X,'*** R factors, final parameters, and their estimate
     &d standard deviations ***'//)
  110 FORMAT(///7X,'*** Nonlinear least-squares fitting by the Marquardt
     & method (NLESQ = 0) ***'//)
  120 FORMAT(///7X,'*** Nonlinear least-squares fitting by the Gauss-New
     &ton method (NLESQ = 1) ***'//)
      IF (IDCYCL .EQ. 1) THEN
         WRITE(6,'(11X,A,I2,5X,A,I3,5X,A,1P,E9.3,5X,A,I2)')
     &   'NAUTO =',NAUTO,'NCYCL =',NCYCL,'CONV = ',CONV,
     &   'NCONV =',NCONV
         IF (NC .NE. 0) WRITE(6,135) NC,TK,FINC
  135    FORMAT(11X,'NCONS =',I3,5X,'T(0) =',1P,G13.6,5X,'FINC =',G13.6)
         IF (NAUTO .NE. 0) THEN
            WRITE(6,137)
            DO JJ = 1, NCHNG
               WRITE(6,138) JJ,(IPAR(J,JJ),J=1,NPAR(JJ))
            END DO
         END IF
C        AMIN: ARRAY STORING PARAMETERS WHICH GIVE THE LOWEST VALUE OF
C              CHISQR
         DO J = 1, NTERMS
            AMIN(J)=A(J)
         END DO
      END IF
  137 FORMAT(//11X,'Numbers for variable parameters in preparatory cycle
     &s')
  138 FORMAT(' ',12X,'Cycle #',I2,':',20I4,:,/5(' ',22X,20I4:/))
C
C        LTK: FLAG TO INDICATE THE ALTERATION OF THE PENALTY PARAMETER
C         IT: COUNTER FOR THE ALTERATION OF THE PENALTY PARAMETER
C        LDM: FLAG TO INDICATE THAT THE DAMPING FACTOR OR MARQUARDT
C             PARAMETER HAS REACHED THE LIMIT VALUE
C      IFAIL: NUMBER OF FAILURE IN THE GAUSS-NEWTON METHOD
C     FLAMDA: MARQUARDT PARAMETER
C       DAMP: DAMPING FACTOR
C      ICONV: INCREASED BY 1 IF RELDEC <= CONV
C     IDCONV: =0, NOT CONVERGED;  =1, CONVERGED
      LTK=0
      IT=0
      LDM=0
      IFAIL=0
      FLAMDA=0.0
      DAMP=1.0
      ICONV=0
      IDCONV=0
*Rev 1.2 2003.05.28 {
      IF (IDCYCL .EQ. 0 .AND. ((NMODE .EQ. 4 .AND. NCONST .EQ. 0)
     &  .OR. NMODE .EQ. 5)) THEN
*        RESOTRE THE PARAMETER GIVING THE MINIMUM SUM OF SQUARES
         DO J = 1, NTOTPAR
            A(J) = SAVEA(J)
         END DO
         CALL SCYPEAK(1)
      END IF
C
C     IREFUC: COUNTER TO INDICATE THE NUMBER OF CYCLES IN WHICH LATTICE
C             PARAMETERS ARE REFINED
      IREFUC=0
      DO ICYCLE = 1, NCYCL
         IF (ICYCLE .EQ. 1 .AND. (NMODE .EQ. 4 .OR. NMODE .EQ. 5))
     &   CALL SCYPEAK(0)
         IF (NAUTO .NE. 0 .AND. ICYCLE .LE. NCHNG+1)
     &   CALL REFPAR(ICYCLE,NTERMS)
C        SAVE THE CURRENT VALUES OF ARRAY A (USED IN SUBROUTINE CHKCNV)
         DO J = 1, NTERMS
            ASAVE(J)=A(J)
         END DO
C
         IF (LTK .EQ. 1) THEN
            SUMP1=FACTOR*SUMPEN
            CHISQ1=CHISQR-SUMPEN+SUMP1
            CSMIN=CHISQ1
C           RESET LTK
            LTK=0
         ELSE IF (ICYCLE .EQ. 1) THEN
            CHISQ1=AOFFN(DEG,XINT,YFIT,NPTS,NSTEP,A,0,IDPOS)
            CHISQR=CHISQ1
            SUMP1=SUMPEN
            CSMIN=CHISQ1
            FLAMDC=0.0
            IF (IDCYCL .EQ. 1) WRITE(6,'(//11X,A)')
     &      'Calculation using initial parameters'
            IF (NC .GT. 0) THEN
               CHAR4=INTCHR(IT)
               WRITE(6,140) CHISQ1,CHISQ1-SUMP1,SUMP1,CHAR4,TK
  140          FORMAT(' ',10X,'AOF =',1P,G13.6,5X,'OF =',
     &         1P,G13.6,5X,'PF =',1P,G13.6,5X,'T',A4,' =',G13.6)
            ELSE
               WRITE(6,'(11X,A,1P,G13.6)') 'OF =',CHISQ1
            END IF
         ELSE
            CHISQ1=CHISQR
            SUMP1=SUMPEN
         END IF
C
         CALL REFMG(ID,IREFUC)
         IF (NBEAM .NE. 0 .AND. IREFUC .NE. 0 .AND. 
     &   MOD(IREFUC,4) .EQ. 0) THEN
            DO I = 1, NREFL
               IF (L12(I) .EQ. 0) THEN
                  CALL ASF
               ELSE
                  DO JJ = 1, NATOM
                     F(JJ,I)=F(JJ,L12(I))
                  END DO
               END IF
            END DO
         ENDIF
         CALL CURFIT(DEG,XINT,NSTEP,NPTS,NTERMS,A,AMIN,
     &   FLAMDA,FLAMDC,FLAMDP,DAMP,IFAIL,YFIT,CHISQ1,CHISQR,CSMIN,
     &   ICYCLE,LDM,NFNEV)
         IF (LDM .EQ. 1) GO TO 4
         IF (IDCYCL .EQ. 1) THEN
            WRITE(6,'(//11X,A,I2)') 'Cycle #',ICYCLE
            IF(NC.EQ.0) THEN
               WRITE(6,'(11X,A,I5,5X,A,1P,G13.6)')
     &         'NPTS =',NPTS,'OF =',CHISQR
            ELSE
               CHAR4=INTCHR(IT)
               WRITE(6,150) NPTS,CHISQR,CHISQR-SUMPEN,SUMPEN,CHAR4,TK
  150          FORMAT(11X,'NPTS =',I6,5X,'AOF =',1P,G13.6,5X,'OF =',
     &         G13.6,5X,'PF =',G13.6,5X,'T',A4,' =',G13.6)
            END IF
            IF (NLESQ .EQ. 0) THEN
               WRITE(6,'(11X,A,I2,5X,A,1P,G9.3)')
     &         'NFNEV =',NFNEV,'LAMBDA = ',FLAMDP
            ELSE
               WRITE(6,'(11X,A,F7.4,5X,A,I2)') 'DAMP =',DAMP,
     &         'NFNEV =',NFNEV
            END IF
         END IF
         IF (NPRINT .GE. 1 .OR. IDCYCL .EQ. 0) THEN
            CALL APRNTD(A,NTERMS)
            IF (IDCYCL .EQ. 0) THEN
*              CALL RFACTR(A,DEG,YFIT,XINT,NPTS,1)
               CALL STDEV(APR,SIGMAA,NTERMS,NPTS,CHISQR-SUMPEN,NESD)
               IF (NLESQ .EQ. 2 .OR. NCYCL .EQ. 0 .OR. NPRINT .EQ. 0)
     &         THEN
                  WRITE(6,'(/12X,A,8X,A,12X,A)') 'No.','A','SIGMA'
               ELSE
                  WRITE(6,'(/12X,A,8X,A,12X,A,6X,A)') 'No.','A','SIGMA',
     &            'DELTA.A/SIGMA'
               END IF
            ELSE
               CALL RFACTR(A,DEG,YFIT,XINT,NPTS,0)
               WRITE(6,'(/12X,A)') 'No.    A(old)    +   Delta.A   '//
     &         '=  A(Refined)    ID'
            END IF
C
            DO JJ = 1, NTERMS
               SELECT CASE (IDCYCL)
                  CASE (0)
                     IF (NLESQ .EQ. 2 .OR. NCYCL .EQ. 0 .OR. 
     &               NPRINT .EQ. 0) THEN
                        SELECT CASE (ID(JJ))
                           CASE (1) 
                              WRITE(6,200)
     &                        JJ,APR(JJ),SIGMAA(JJ),PARNAM(JJ)
                           CASE DEFAULT
                              WRITE(6,210) JJ,APR(JJ),PARNAM(JJ)
                        END SELECT
                     ELSE
                        SELECT CASE (ID(JJ))
                           CASE (1)
                              WRITE(6,205) JJ,APR(JJ),SIGMAA(JJ),
     &                        ABS(DELTA(JJ))/SIGMAA(JJ),PARNAM(JJ)
                           CASE DEFAULT
                              WRITE(6,215) JJ,APR(JJ),PARNAM(JJ)
                        END SELECT
                     END IF
                  CASE DEFAULT
                     IF (ID(JJ) .NE. 0 .OR. LATVAR(JJ) .EQ. 1) THEN
                        DELTA(JJ)=APR(JJ)-AOUT(JJ)
                     ELSE
                        DELTA(JJ)=0.0
                     ENDIF
C                    AOUT: THE VALUES OF ARRAY A IN THE PREVIOUS CYCLE
                     WRITE(6,220) JJ,AOUT(JJ),DELTA(JJ),
     &               APR(JJ),ID(JJ),PARNAM(JJ)
                     AOUT(JJ)=APR(JJ)
               END SELECT
            END DO
         ENDIF
*        SCALE FACTORS AND ITS E.S.D.'S (WILL BE USED IN SUBROUTINE RFLCTN)
         DO IPH = 1, NPHASE
            SCALEAS(IPH) = APR(KPHB(IPH)+0)
            SIGSCALS(IPH) = SIGMAA(KPHB(IPH)+0)
         END DO

         IF (IDCYCL .EQ. 0 .AND. NPRINT .GE. 1) THEN
            CALL TENSOR(A)
            CALL PRUCC
         ENDIF
C
    4    IF (IDCYCL .EQ. 1 .AND. ICYCLE .GT. NCHNG+1) THEN
            CALL CHKCNV(IDCONV,A,ASAVE,CHISQ1,CHISQR,ICONV,NTERMS)
         END IF
         IF ((NC .EQ. 0 .OR. SUMP1 .EQ. 0.0 .OR. IFAIL .GE. 2) .AND.
     &   ((LDM .EQ. 1 .AND. ICYCLE .GT. NCHNG) .OR. IDCONV .EQ. 1)) THEN
            DO J = 1, NTERMS
               A(J)=AMIN(J)
            END DO
            RETURN
         ELSE IF((NAUTO .NE. 0 .AND. ICYCLE .LE. NCHNG) .OR. 
     &   (IDCYCL .EQ. 1 .AND. SUMPEN .GT. 0.0 .AND. 
     &   (LDM .EQ. 1 .OR. IDCONV .EQ. 1))) THEN
C           RESET THE MARQUARDT PARAMETER, DAMPING FACTOR, COUNTER OF
C           FAILURE
            FLAMDA=0.0
            DAMP=1.0
            IFAIL=0
            IF (LDM .EQ. 1) THEN
               DO JJ = 1, NTERMS
                  A(JJ)=AMIN(JJ)
               END DO  
               LDM=0
            ENDIF
            IF (NAUTO .NE. 0 .AND. ICYCLE .LE. NCHNG) CYCLE
C           INCREASE THE PENALTY PARAMETER
            IF (SUMPEN .GT. 0.2*CHISQR) THEN
               FACTOR=MIN(FINC,1.5)
               IF (ICYCLE .EQ. 1 .OR. FACTOR*SUMPEN .GT. 0.8*CHISQR) 
     &         THEN
                  WRITE(6,250)
                  RETURN
               ENDIF
            ELSE
               FACTOR=MIN(FINC,0.25*CHISQR/SUMPEN)
            END IF
            TK=FACTOR*TK
            IT=IT+1
            IF (TK .GT. 10E25) THEN
               WRITE(6,250)
               RETURN
            END IF
            LTK=1
         END IF
      END DO
  200 FORMAT(' ',10X,I3,2X,2(1P,G13.6,2X),3X,A)
  205 FORMAT(' ',10X,I3,2X,3(1P,G13.6,2X),3X,A)
  210 FORMAT(' ',10X,I3,2X,1P,G13.6,20X,A)
  215 FORMAT(' ',10X,I3,2X,1P,G13.6,35X,A)
  220 FORMAT(' ',10X,I3,2X,1P,3(G13.6,1X),I4,3X,A)
  250 FORMAT(//11X,'Too large penalty term.  Stop refining.')
      END

************************************************************************

      INTEGER FUNCTION LATVAR(JJ)
C     EXAMINE WHETHER OR NOT AT LEAST ONE OF THE LATTICE PARAMETERS
C     HAS BEEN REFINED
      PARAMETER (NPH=8,NT=999,NAP=150)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN

      DO 50 L=1,NPHASE
         IF (JJ .GE. KPHB(L)+NAX .AND. JJ .LE. KPHB(L)+NGAMX) THEN
            DO 10 J=KPHB(L)+NAX,KPHB(L)+NGAMX
               IF (ID(J) .EQ. 1) THEN
                  LATVAR=1
                  RETURN
               END IF
   10       CONTINUE
         END IF
   50 CONTINUE
      LATVAR=0
      END

************************************************************************

      SUBROUTINE STDEV(APR,SIGMAA,NTERMS,NPTS,CHISQR,NESD)
C     CALCULATE STANDARD DEVIATIONS OF REFINABLE PARAMETERS
      PARAMETER (NB=7000,NT=999,NPH=8,NR=400,NAP=150,NOD=NR*(NR+1)/2)
      DOUBLE PRECISION G
      DIMENSION APR(*),SIGMAA(*)
      COMMON /B/ U(NB),NREFL
      COMMON /ERMAT/ SI2,SP2
      COMMON /H/ NQ(NPH),NCENTR(NPH),LAUEG(NPH)
      COMMON /LTM/ LISO(NAP,NPH),LMAG(NPH)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /TWO/ G(NOD),AD(NR)
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      DATA RD,RDSQR/0.0174532925,0.0003046174/
      LP(J)=KPHB(L)+J
      P(J)=APR(LP(J))
      LD(J)=ID(LP(J))
      R(J)=SIGMAA(LP(J))
C
      DO 10 J=1,NTERMS
         IF(ID(J).EQ.1) THEN
C           DIAGONAL ELEMENT OF THE INVERSE OF THE MAXTIX OF THE
C           NORMAL EQUATION
            SIGMAA(J)=G(IR(J)*(IR(J)+1)/2)
         ELSE
            SIGMAA(J)=0.0
         ENDIF
   10 CONTINUE
C
C     ADJUST RIETVELD E.S.D. TO PROVIDE COMPARABILITY WITH INTEGRATED
C     INTENSITY REFINEMENT
C     CF. H. G. SCOTT, J. APPL. CRYSTALLOGR., 16 159 (1983)
C     NPP: NUMBER OF PARAMETERS DESCRIBING THE PROFILE
C     NPC: NUMBER OF CRYSTAL STRUCTURE PARAMETERS
      SP2=CHISQR/(NPTS-NRFN)
      NPP=0
      DO 20 J=1,KPHB(1)-1
         IF(ID(J).EQ.1) THEN
            SIGMAA(J)=SQRT(SIGMAA(J)*SP2)
            NPP=NPP+1
         ENDIF
   20 CONTINUE
   
*     RAD ==> DEG (PEAK-SHIFT PARAMETERS)
      DO J = 1, NROUGH-1
         IF(ID(J).EQ.1) SIGMAA(J)=SIGMAA(J)/RD
      END DO
      
      IF (NPRIM .GT. 0) THEN
         DO J = NPPPX, NPPPX + (NRELAX - 1)*NPRIM, NPRIM
*           RAD ==> DEG (FWHM FOR A RELAXED REFLECTION)
            IF (ID(J) .EQ. 1) SIGMAA(J) = SIGMAA(J)/RD
            IF (ID(J+1) .EQ. 1 .AND. NPRFN .EQ. 2) 
     &      SIGMAA(J+1) = SIGMAA(J+1)/RD
         END DO
      END IF

      DO L=1,NPHASE
      
***      IF A NEW PROFILE FUNCTION IS IMPLEMENTED, THE FOLLOWING PART MUST BE
***      MODIFIED
*        RAD**2 ==> DEG**2 (U, V, W, AND P)
         DO JJ=LP(1),LP(4)
            SIGMAA(JJ)=SIGMAA(JJ)/RDSQR
         END DO
         IF (NPRFN .EQ. 0) THEN
*           RAD ==> DEG (X, Xe, Y, Ye, and As)
            DO JJ = LP(5), LP(9)
* V2.0, 1999.04.02, Izumi
               IF (NASYM .EQ. 0 .AND. JJ .EQ. LP(9)) EXIT
               SIGMAA(JJ)=SIGMAA(JJ)/RD
            END DO
         END IF
         
         DO J=LP(NSFX),LP(NGAMX)
            IF(ID(J).EQ.1) THEN
               SIGMAA(J)=SQRT(SIGMAA(J)*SP2)
               NPP=NPP+1
            ENDIF
         END DO
      END DO

*     THE STANDARD DEVIATIONS OF STRUCTURE PARAMETERS ARE CALCULATED
*     ACCORDING TO THE PROCEDURE OF SCOTT (1983).
      IF (NESD .EQ. 1) THEN
         NPC=NRFN-NPP
         IF (NREFL .NE. NPC) THEN
            SI2=1.0+(CHISQR+FLOAT(NPP-NPTS))/FLOAT(NREFL-NPC)
         ELSE
            SI2=-1.0
         END IF
         IF (SI2 .LT. 0.0) THEN
            WRITE(6,'(//11X,A/11X,A/11X,A//)')
     &      'The value of SI**2 is negative or its denominator is zero.'
     &      ,'The standard deviations of structural parameters are '
     &      //'calculated according to the conventional method.'
         END IF
      END IF
	  
      DO L = 1, NPHASE
         DO J = LP(NOTX), KPHE(L)
            IF (NESD .EQ. 0) THEN
*              CONVENTIONAL METHOD OF CALCULATING E.S.D.'S
               SIGMAA(J) = SQRT(SIGMAA(J)*SP2)
            ELSE IF (ID(J) .EQ. 1 .AND. SI2 .GE. 0.0) THEN
*              SCOTT'S METHOD OF CALCULATING E.S.D.'S
               SIGMAA(J) = SQRT(SIGMAA(J)*SI2)
            ELSE IF (ID(J) .EQ. 1 .AND. SI2 .LT. 0.0) THEN
*              CONVENTIONAL METHOD OF CALCULATING E.S.D.'S
               SIGMAA(J) = SQRT(SIGMAA(J)*SP2)
            END IF
         END DO
         IF (LMAG(L) .EQ. 1) THEN
*Rev 1.0c 2000.11.08 Izumi
            IF (LAUEG(L) .LE. 5) THEN
               DO J = KPHE(L)-2, KPHE(L)
                  IF (ID(J) .EQ. 1) SIGMAA(J) = SIGMAA(J)/RD
               END DO
            ELSE
               IF (ID(KPHE(L)) .EQ. 1) SIGMAA(KPHE(L)) =
     &         SIGMAA(KPHE(L))/RD
            END IF
         END IF
      END DO
C
      DO 90 L=1,NPHASE
         CALL CPSIG(L,LAUEG,SIGMAA,KPHB)
C        COMPUTE STANDARD DEVIATIONS OF THE DIRECT UNIT CELL DIMENSIONS
         A=P(NAX)
         B=P(NBX)
         C=P(NCX)
         ALPHA=P(NALPX)*RD
         BETA=P(NBETX)*RD
         GAMMA=P(NGAMX)*RD
         S1=R(NAX)**2
         S2=R(NBX)**2
         S3=R(NCX)**2
         S4=R(NALPX)**2
         S5=R(NBETX)**2
         S6=R(NGAMX)**2
C        PROPAGATION OF THE ERROR
C        CF. C. H. KELSEY, MINERAL. MAG., 33, 809 (1964)
C        AND EQ. (5.3-82) IN T. SAKURAI, "XRAY-ANALYSIS OF CRYSTAL
C        CRYSTAL STRUCTURES," P. 233.
C
         IF(LD(NAX).EQ.1)
     &   SIGMAA(LP(NAX))=(0.5*A**3)**2*S1+(0.5*A*B**2*COS(GAMMA)
     &   **2)**2*S2+(0.5*A*C**2*COS(BETA)**2)**2*S3+(A*B*C*COS(BETA)
     &   *COS(GAMMA))**2*S4+(A**2*C*COS(BETA))**2*S5+(A**2*B*COS(GAMMA))
     &   **2*S6
         IF(LD(NBX).EQ.1)
     &   SIGMAA(LP(NBX))=(0.5*A**2*B*COS(GAMMA)**2)**2*S1+(0.5*B**3
     &   )**2*S2+(0.5*B*C**2*COS(ALPHA)**2)**2*S3+(B**2*C*COS(ALPHA)
     &   )**2*S4+(A*B*C*COS(ALPHA)*COS(GAMMA))**2*S5+(A*B**2*COS(GAMMA))
     &   **2*S6
         IF(LD(NCX).EQ.1)
     &   SIGMAA(LP(NCX))=(0.5*A**2*C*COS(BETA)**2)**2*S1+(0.5*B**2
     &   *C*COS(ALPHA)**2)**2*S2+(0.5*C**3)**2*S3+(B*C**2*COS(ALPHA)
     &   )**2*S4+(A*C**2*COS(BETA))**2*S5+(A*B*C*COS(ALPHA)*COS(BETA))**
     &   2*S6
         IF(LD(NALPX).EQ.1)
     &   SIGMAA(LP(NALPX))=(A**2*(2.0*COS(BETA)*COS(GAMMA)-COS(ALPHA)*
     &   (COS(BETA)**2+COS(GAMMA)**2))/(2.0*SIN(ALPHA)))**2*S1+
     &   (0.25*B**2*SIN(2.0*ALPHA))**2*S2+(0.25*C**2*SIN(2.0
     &   *ALPHA))**2*S3+(B*C*SIN(ALPHA))**2*S4+(A*C*SIN(ALPHA)
     &   *COS(GAMMA))**2*S5+(A*B*SIN(ALPHA)*COS(BETA))**2*S6
         IF(LD(NBETX).EQ.1)
     &   SIGMAA(LP(NBETX))=(0.25*A**2*SIN(2.0*BETA))**2*S1+(B**2*(2.0
     &   *COS(ALPHA)*COS(GAMMA)-COS(BETA)*(COS(ALPHA)**2+COS(GAMMA)
     &   **2))/(2.0*SIN(BETA)))**2*S2+(0.25*C**2*SIN(2.0*BETA))**2
     &   *S3+(B*C*SIN(BETA)*COS(GAMMA))**2*S4+(A*C*SIN(BETA))**2*S5+
     &   (A*B*SIN(BETA)*COS(ALPHA))**2*S6
         IF(LD(NGAMX).EQ.1)
     &   SIGMAA(LP(NGAMX))=(0.25*A**2*SIN(2.0*GAMMA))**2*S1+(0.25
     &   *B**2*SIN(2.0*BETA))**2*S2+(C**2*(2.0*COS(ALPHA)*COS(BETA)-
     &   COS(GAMMA)*(COS(ALPHA)**2+COS(BETA)**2))/(2.0*SIN(GAMMA))
     &   )**2*S3+(B*C*SIN(GAMMA)*COS(BETA))**2*S4+(A*C*SIN(GAMMA)
     &   *COS(ALPHA))**2*S5+(A*B*SIN(GAMMA))**2*S6
         DO 70 J=NAX,NCX
            IF(LD(J).EQ.1) SIGMAA(LP(J))=SQRT(SIGMAA(LP(J)))
   70    CONTINUE
C
C        RAD ==> DEG
         DO 80 J=NALPX,NGAMX
            IF(LD(J).EQ.1) SIGMAA(LP(J))=SQRT(SIGMAA(LP(J)))/RD
   80    CONTINUE
         CALL CPSIG(L,LAUEG,SIGMAA,KPHB)
   90 CONTINUE
      END

************************************************************************

      SUBROUTINE CPSIG(L,LAUEG,SIGMAA,KPHB)
*     COPY STANDARD DEVIAIONS OF LATTICE PARAMETERS OR ELEMENTS OF
*     THE METRIC TENSOR TO CONSTRAINED PARAMETERS
      INTEGER L,LAUEG(*),KPHB(*)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      REAL SIGMAA(*)
      LP(J)=KPHB(L)+J
      R(J)=SIGMAA(LP(J))

      IF(LAUEG(L).EQ.6 .OR. LAUEG(L).EQ.7 .OR. LAUEG(L).EQ.9 .OR.
     &LAUEG(L).EQ.11 .OR. LAUEG(L).EQ.12 .OR. LAUEG(L).EQ.13) THEN
         SIGMAA(LP(NBX))=R(NAX)
      ELSEIF(LAUEG(L).EQ.8 .OR. LAUEG(L).EQ.10) THEN
         SIGMAA(LP(NBX))=R(NAX)
         SIGMAA(LP(NCX))=R(NAX)
         SIGMAA(LP(NBETX))=R(NALPX)
         SIGMAA(LP(NGAMX))=R(NALPX)
      ELSEIF(LAUEG(L).GE.14) THEN
         SIGMAA(LP(NBX))=R(NAX)
         SIGMAA(LP(NCX))=R(NAX)
      END IF
      END

************************************************************************

      SUBROUTINE REFPOW(A,NTERMS)
C     EXAMINE WHETHER OR NOT A GROUP OF PARAMETERS ARE REFINED IN THE
C     CONJUGATE DIRECTION METHOD
      PARAMETER (NPH=8,NT=999,NAP=150)
      DIMENSION A(*)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /OLD/ AOLD(NT),AOUT(NT)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /UCC/ IDPEAK(NPH),IDSF(NPH),IDPRF(NPH),IDPO(NPH)
      LP(J)=KPHB(L)+J

      DO L = 1, NPHASE
         IDPRF(L) = 0
      END DO
      
*     CHECK PARAMETERS RELATED TO PEAK SHIFTS
      LSKIP = 0
      DO J = 1, NBCKGR - 1        
         IF (ID(J) .EQ. 1 .AND. A(J) .NE. AOLD(J)) THEN
            DO L = 1, NPHASE
               IDPRF(L)=1
            END DO
            LSKIP = 1
            EXIT
         END IF
      END DO

      IF (LSKIP .EQ. 0) THEN
*        CHECK PRIMARY PROFILE PARAMETERS
         DO J = NPPPX, NPPPX + NRELAX*NPRIM -1
            IF (ID(J) .EQ. 1) THEN
               DO L = 1, NPHASE
                  IDPRF(L)=1
               END DO
               EXIT
            END IF
         END DO
      END IF

      DO 40 L = 1, NPHASE
         IDPEAK(L)=0
         IDSF(L)=0
         IDPO(L)=0

*        CHECK PREFERRED-ORIENTATION PARAMETERS
         DO J = LP(NPO1X), LP(NPO2X)
            IF (ID(J) .EQ. 1 .AND. A(J) .NE. AOLD(J)) IDPO(L)=1
         END DO

*        CHECK LATTICE PARAMETERS    
         DO J = LP(NAX), LP(NGAMX)
            IF (ID(J) .EQ. 1 .AND. A(J) .NE. AOLD(J)) THEN
               IDPEAK(L)=1
               IDPRF(L)=1
               IF (NMODE .EQ. 0) IDSF(L)=1
               GO TO 40
            END IF
         END DO

*        CHECK SECONDARY PROFILE PARAMETERS
         IF (IDPRF(L) .EQ. 0) THEN         
            DO J = LP(1), LP(NPO1X-1)
               IF (ID(J) .EQ. 1) THEN
                  IDPRF(L)=1
                  EXIT
               END IF
            END DO
         END IF

*        CHECK STRUCTURE PARAMETERS
         DO J = LP(NOTX), KPHE(L)
            IF (ID(J) .EQ. 1 .AND. A(J) .NE. AOLD(J)) THEN
               IDSF(L)=1
               EXIT
            END IF
         END DO
   40 CONTINUE

      DO J = 1, NTERMS
         AOLD(J)=A(J)
      END DO
      END

************************************************************************

      SUBROUTINE REFMG(ID,IREFUC)
C     EXAMINE WHETHER OR NOT A GROUP OF PARAMETERS ARE REFINED IN THE
C     GAUSS-NEWTON OR MARQUARDT METHOD
      PARAMETER (NPH=8,NAP=150)
      DIMENSION ID(*)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /UCC/ IDPEAK(NPH),IDSF(NPH),IDPRF(NPH),IDPO(NPH)
      LP(J)=KPHB(L)+J

      DO L = 1, NPHASE
         IDPRF(L) = 0
      END DO
      
*     CHECK PARAMETERS RELATED TO PEAK SHIFTS
      LSKIP = 0
      DO J = 1, NBCKGR -1         
         IF (ID(J) .EQ. 1) THEN
            DO L = 1, NPHASE
               IDPRF(L)=1
            END DO
            LSKIP = 1
            EXIT
         END IF
      END DO

      IF (LSKIP .EQ. 0) THEN
*        CHECK PRIMARY PROFILE PARAMETERS
         DO J = NPPPX, NPPPX + NRELAX*NPRIM -1
            IF (ID(J) .EQ. 1) THEN
               DO L = 1, NPHASE
                  IDPRF(L)=1
               END DO
               EXIT
            END IF
         END DO
      END IF

      IFIRST=IREFUC
      DO 40 L = 1, NPHASE
         IDPEAK(L)=0
         IDSF(L)=0
         IDPO(L)=0

*        CHECK PREFERRED-ORIENTATION PARAMETERS
         DO J = LP(NPO1X), LP(NPO2X)
            IF (ID(J) .EQ. 1) THEN
               IDPO(L)=1
            END IF
         END DO

*        CHECK LATTICE PARAMETERS    
         DO J = LP(NAX), LP(NGAMX)
            IF (ID(J) .EQ. 1) THEN
               IDPEAK(L)=1
               IDPRF(L)=1
               IF (NMODE .EQ. 0) IDSF(L)=1
               IF(IFIRST.EQ.IREFUC) IREFUC=IREFUC+1
               GO TO 40
            ENDIF
         END DO

*        CHECK SECONDARY PROFILE PARAMETERS
         IF (IDPRF(L) .EQ. 0) THEN         
            DO J = LP(1), LP(NPO1X-1)
               IF(ID(J) .NE. 0) THEN
                  IDPRF(L)=1
                  EXIT
               ENDIF
            END DO
         END IF

*        CHECK STRUCTURE PARAMETERS
         DO J = LP(NOTX), KPHE(L)
            IF(ID(J).EQ.1) THEN
               IDSF(L)=1
               EXIT
            ENDIF
         END DO
   40 CONTINUE
      END

************************************************************************

      SUBROUTINE ABSVARL(NSTEP,SIGA,SIGI,ATOMNO,NREAL,DEG,RADIUS,MUR,
     &  SABS,DSANG,RGON,SWIDTH)
*     CALCULATE ABSORPTION FACTORS AND ILLUMINATION LENGTHS
      PARAMETER (NP=80000)
      REAL SIGA(*),SIGI(*),ATOMNO(*),DEG(*),RADIUS,MUR
      COMMON /ABSVR/ ABSORP(NP),VARLEN(NP),ROUGH(NP),NSURFR
      COMMON /MC/ CTHM,NTRAN,PCOR
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD

      DO J = 1,NSTEP     
*        ABSORP: ABSORPTION FACTOR
*        VARLEN: ILLUMINATION LENGTH FROM A VARIABLE DIVERGENCE SLIT
         SELECT CASE (NBEAM)
            CASE (0)
*              NEUTRON DIFFRACTION
               ABSORP(J) = ABSCOR(SIGA,SIGI,ATOMNO,NREAL,
     &                     (SIN(0.5*DEG(J)))**2,RADIUS,MUR)
               VARLEN(J) = 1.0
            CASE (1)
*              CONVENTIONAL X-RAY DIFFRACTION
               SELECT CASE (NTRAN)
                  CASE (0)
*                    BRAGG-BRENTANO & FIXED DS
                     ABSORP(J) = 1.0
                     VARLEN(J) = 1.0
                  CASE (1)
*                    BRAGG-BRENTANO & AUTO DS
*                    C. Weidenthaler et al., Acta Crystallogr., Sect. B, 
*                    53 (1997) 429.
                     ABSORP(J) = 1.0
                     TERM1 = 0.5*DSANG*COS(0.5*DEG(J))/SIN(0.5*DEG(J))
                     TERM2 = RGON/(SWIDTH*COS(0.5*DEG(J)))
                     VARLEN(J) = (SQRT(TERM2**2 + 1.0) - TERM2)/TERM1
                  CASE (2)
*                    TRANSMISSION GEOMETRY
*                    Refer TO Eq. (2.3.1.20) in International Tables, Vol. C
                     ABSORP(J) = 
     &               EXP(-SABS/COS(0.5*DEG(J)))/COS(0.5*DEG(J))
                     VARLEN(J) = 1.0
                  CASE (3)
*                    DEBYE-SCHERRER GIOMETRY
*                    REFER TO EQ. (3.25) IN THE MANUAL OF Fullprof
                     ABSORP(J) = ABSCOR(SIGA,SIGI,ATOMNO,NREAL,
     &                           (SIN(0.5*DEG(J)))**2,RADIUS,MUR)
                     VARLEN(J) = 1.0
               END SELECT
            CASE (2)
*              SYNCHROTRON X-RAY DIFFRACTION
               IF (MUR .EQ. 0.0) THEN
*                 FLAT-PLATE SAMPLE
                  ABSORP(J) = 1.0
               ELSE
*                 DEBYE-SCHERRER GEOMETRY
                  ABSORP(J) = ABSCOR(SIGA,SIGI,ATOMNO,NREAL,
     &                        (SIN(0.5*DEG(J)))**2,RADIUS,MUR)
               END IF
               VARLEN(J) = 1.0
         END SELECT
      END DO
      END

************************************************************************

      SUBROUTINE INCDEG(A,ANGMAX)
*     INCREASE DEG2 IN SUCH A WAY THAT SOME REFLECTIONS WHOSE 
*     DIFFRACTION ANGLES ARE SLIGHTLY HIGHER THAN DEG2 ARE INCLUDED
      PARAMETER (NPH=8,NAP=150)
      REAL A(*)
      COMMON /L/ DEG1,DEG2
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      COMMON /Z/ PC,CHGPC
      DATA RD2/57.29578/
      P(J)=A(KPHB(L)+J)

      TTH = TAN(0.5*DEG2)
*     Include reflections whose parts of profiles overpass 2-theta(max)
C     This value should be changed appropriately.
      PCEXT = 2.0
*     PROFILE PARAMETERS FOR THE FIRST PHASE WILL BE USED
      L = 1
      IF (NPRFN .EQ. 0) THEN
*        USE U, V, AND W FOR THE 1ST PHASE
         CSTH = COS(0.5*DEG2)
         SIG = (P(1)*TTH + P(2))*TTH + P(3) + P(4)/CSTH**2
*        Neglect anisotropic broadening
         GAM = P(5) / CSTH + P(7) * TTH
         IF (SIG .LT. 0.0 .OR. GAM .LT. 0.0) THEN
            CALL JOBEND('Negative FWHM(Gauss)**2 and/or '//
     &      'FWHM(Lorentz) resulted from the current profile'//
     &      'parameters at the highest angle')
         END IF
         FWHG = 2.354820 * SQRT(SIG) 
         FWHMGL =  FWHG + GAM
         ANGMAX = RD2*(DEG2 + PCEXT*FWHMGL)
      ELSE
         FWHMGL = (P(1)*TTH + P(2))*TTH + P(3) 
         IF (FWHMGL .LT. 0.0) THEN
            CALL JOBEND('A Negative FWHM resulted from the current '
     &      //'profile parameters at the highest angle')
         END IF
         COSCTH = 1.0/SIN(0.5*DEG2)
         ASYM2 = P(5) + P(6)*(1.4142136 - COSCTH) + 
     &          P(7)*(2.0 - COSCTH**2)
         ANGMAX = RD2*(DEG2 + PCEXT*2.0*ASYM2/(1.0+ASYM2)*SQRT(FWHMGL))
      END IF
      END

