* RIETAN-2000, Rev 1.2, 2003.05.28

C    ***************************************************************
C    ***************************************************************
C    **                                                           **
C    **                       RIETAN-2000                         **
C    **                                                           **
C    **          A FORTRAN PROGRAM FOR PATTERN-FITTING            ** 
C    **   STRUCTURE REFINEMENT AND PATTERN DECOMPOSITION WITH     **
C    **   ANGLE-DISPERSIVE X-RAY AND NEUTRON DIFFRACTION DATA     **
C    **                                                           **
C    **            COPYRIGHT (C) 2000, BY F. IZUMI                **
C    **                  ALL RIGHTS RESERVED                      **
C    **                                                           **
C    **  NATIONAL INSTITUTE FOR RESEARCH IN INORGANIC MATERIALS   **
C    **                                                           **
C    ***************************************************************
C    ***************************************************************

* FILES USED IN THIS PROGRAM
*  #1: spgri
* #11: spgra
*  #2: asfdc
*  #3: INTENSITY DATA OF THE SAMPLE, *.int
*  #4: INPUT DATA CREATED BY A PREPROCESSOR (RECORD LENGTH = 80)
*  #5: ORIGINAL INPUT DATA INCLUDING COMMENTS (RECORD LENGTH = 80), *.ins
*  #6: PRINTER OUTPUT, *.lst
*  #7: A TEMPORAY FILE STORING 2-THETA AND OBSERVED INTENSITIES
*  #8: BACKGROUND INTENSITIES, *.bkg
*  #9: OUTPUT FOR ORFFE, *.xyz
* #10: *.ffe FILE CREATED BY ORFFE
* #20: OUTPUT FOR IGOR PRO, *.pat
* #21: OUTPUT FOR FOUSYN, *.hkl
* #22: INITIAL INTEGRATED INTENSITIES FOR LE BAIL REFINEMENT, *.ffi
* #23: INTEGRATED INTENSITIES RESULTING FROM LE BAIL REFINEMENT, *.ffo
* #30: OUTPUT FOR MEED, *.mem
* #32: OUTPUT FROM MEED, *.fba
* #55: COPY OF FILE #5 USED DURING THE UPDATE OF REFINED PARAMETERS
* #60: FILE STORING hkl and m, hklm.phn

      PROGRAM RIETAN
C     NB:  MAXIMUM NUMBER OF PEAKS
C     NP:  MAXIMUM NUMBER OF INTENSITY DATA
C     NA:  MAXIMUM NUMBER OF CHEMICAL SPECIES INCLUDING VIRTUAL ONES
C     NT:  MAXIMUM NUMBER OF PARAMETERS
C     NR:  MAXIMUM NUMBER OF REFINABLE PARAMETERS
C     NSF: MAXIMUM NUMBER OF REFINABLE STRUCTURE-FACTOR PARAMETERS
C          (Q, G, X, Y, Z, DISPLACEMENT PARAMETER(S), MU, AND PHI)
C     NS:  MAXIMUM NUMBER OF SYMMETRY OPERATIONS (<=48)
C     NAP: MAXIMUM NUMBER OF ATOMS IN THE ASYMMETRIC UNIT
C     NCS: MAXIMUM NUMBER OF LINEAR CONSTRAINTS
C     NPH: MAXIMUM NUMBER OF PHASES
C     NGEN: MAXIMUM NUMBER OF REFLECTIONS GENERATED TEMPORARILY IN THE
C           PROGRAM (NGEN >= NB)
C     MAXLAB: MAXIMUM NUMBER OF LABELS
 
C     CHANGE THE FOLLOWING PARAMETER VALUES IF NECESSARY
      PARAMETER (NB=7000,NP=80000,NA=15,NT=999,NR=400,NSF=400,
C    &NS=48,NAP=150,NCS=400,NPH=8,NGEN=7000,MAXNL=1000,NCOMPAR=400,
     &NS=48,NAP=150,NCS=400,NPH=8,NGEN=7000,MAXNL=1000,MAXLAB=800)
      PARAMETER (NOD=NR*(NR+1)/2)

      INTEGER H,HP,R,RX,RY,RZ,SECSET,HANIS
      INTEGER NSPLBL(NT)
      REAL MUR
      DOUBLE PRECISION G
      CHARACTER TNAME*2,XTARG(7)*2,CRYSYS(15)*36,
     &  ANAME*5,SYM*10,SPGR*10,PARNAM*60,
     &  TITLE*80,LINE*80,FSSET*14,PHNAME*25,
     &  LINNSS(14)*40,REFTYP(NPH)*3,LCON(7,NPH)*80,ANAME2(NA)*5,
     &  VNS(NPH)*7,LABEL(NT)*25,TITMEM*70
      CHARACTER*50 FILE8*50
      DIMENSION DEG(NP),XINT(NP),YFIT(NP),A(NT),ANAME(NA),
     &  X(100),Y(100),DC(2,6),IDREAL(NA),SYM(3,NS,NPH),
     &  ISONUM(NAP,NPH),WL(4,7),NSPGR(NPH),AW(NA),INDIV(NPH),
     &  NSET(NPH),SIGMAA(NT),LPAR(MAXLAB),SIGA(NA),ATOMNO(NA),
     &  CATOM(NA),SIGI(NA)
      COMMON /EXCREG/ DEGEXC(2,50),NSKIP
      COMMON /A/ H(NB),K(NB),IL(NB),L12(NB),NOPH(NB)
      COMMON /ABSVR/ ABSORP(NP),VARLEN(NP),ROUGH(NP),NSURFR
      COMMON /ALP/ LALPHA
      COMMON /ANISBR/ HANIS(NPH),KANIS(NPH),LANIS(NPH),
     &  CANIS(NPH),PCOS(NB)
      COMMON /B/ U(NB),NREFL
*Rev 1.0h 2000.12.13 Izumi
      COMMON /BIJVOET/ IPAIR(NB),LPAIR(NPH)
      COMMON /C/ BG(NP),BGINC(NP),DEGNOR(NP)
      COMMON /CC/ APR(NT)
      COMMON /CNSTR/ C1(NCS),C2(10,NCS),IN(10,NCS),NCNPAR(NCS),
     &  NACNS(NCS),NCNSTR
      COMMON /CONHKL/ NHKL(NPH),ICH(5,3,7,NPH),NCON(7,NPH),MLINE(NPH)
      COMMON /CONSDA/ NPHCON,NSCONS(4,MAXNL),XYZC(4,3,3,MAXNL),
     &  EXPCTD(MAXNL),DEVDA(MAXNL),IXCON(NAP)
      COMMON /DERPR/ POP2(NB),COMPO(NB),CS(NB),CCC(6,NB)
      COMMON /DETECT/ VARINV(NP)
      COMMON /E/ AA(4,NA),BB(4,NA),C(NA)
      COMMON /EXC/ LSPSYM(NPH),NSSYM(14,NPH)
      COMMON /FL0/ NC,TK,CONV,NCONV,FINC,SUMPEN,SUMP1,FACTOR
      COMMON /FRDA/ NFR,NMEM,MEED(15),EPSD,SCIO,TSCAT,TSCAT1,TSCAT2,
     &  UCLAG,NDA,NFPX(3)
      COMMON /FRDA2/ TITMEM
      COMMON /FT/ FT1,FT2,PEAKNOR(NB),TANTHNOR(NB) 
      COMMON /G/ I
      COMMON /H/ NQ(NPH),NCENTR(NPH),LAUEG(NPH)
      COMMON /HALF/ IFAC(NAP,NPH)
      COMMON /I/ RDEG(2,NB)
      COMMON /INPT/ DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,DUM7,ANGMAX,ANGMIN
      COMMON /IS/ ISET(44),SECSET
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /L/ DEG1,DEG2
      COMMON /LTM/ LISO(NAP,NPH),LMAG(NPH)
      COMMON /MAGNE/ FMAG(NA,NB),SFRMAG(NB),SFIMAG(NB),Q2(NB)
      COMMON /MAGSC/ MAGATM(NA),LCMFF(NA),CMFF(7,NA)
      COMMON /MC/ CTHM,NTRAN,PCOR
      COMMON /MIX/ FRCTN(5,10),MFRAC(5,10),NOFRAC(10),NMIX
*Rev 1.1 2003.05.16 Izumi (also changed in other subroutines)
      COMMON /NONLIN/ NLESQ,NAUTO,IDCYCL,NCHNG,NPAR(50),IPAR(NR,50)
      COMMON /NUML/ NUPDT
*Rev 1.2 2003.05.28 Izumi (also changed in other subroutines)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /OLD/ AOLD(NT),AOUT(NT)
      COMMON /ORDER/ IPOINT(NB)
      COMMON /PAR/ PARNAM(NT),PHNAME(NPH)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
*     'IPHASE'IS REPLACED WITH 'JPHASE' HERE BECAUSE IT IS ALSO USED IN  
*     COMMON AREA /UNCONS/
      COMMON /PHNO/ JPHASE
      COMMON /POW/ MITER,STEP2,ACC
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /PRLV/ NPRINT
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /RLS/ NRANGE
      COMMON /RT/ T(3,NS,NPH),HT(NS,NB),NSYM(NPH),NSITE(NPH),
     &  IDSYM(NS,NAP,NPH),R(3,3,NS,NPH),RX(NS,NB),RY(NS,NB),RZ(NS,NB)
      COMMON /S/ F(NA,NB),NATOM
      COMMON /S2/ F0(NA,NB)
      COMMON /SFRI/ SFR(NB),SFI(NB),OTF(NB)
      COMMON /ABPI/ AP(NB),BP(NB)
      COMMON /SG/ SPGR(NPH)
      COMMON /SITE/ NOAT(NAP,NPH)
      COMMON /SGN/ NSG9,NSC9
      COMMON /T/ FF(NB),COEF(NPH)
      COMMON /TENS/ DUMG(6,NPH),GG(6,NPH)
      COMMON /TWO/ G(NOD),FIRST(NR)
      COMMON /U/ DELTF1(NA),DELTF2(NA)
      COMMON /UCC/ IDPEAK(NPH),IDSF(NPH),IDPRF(NPH),IDPO(NPH)
      COMMON /V/ NPROR(NPH),HP(NPH),KP(NPH),LP(NPH),LSUM(NPH)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
      COMMON /X/ STEP
*Rev 1.1 2003.05.19 Izumi
      COMMON /Z/ PC,CHGPC
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      COMMON /ZSPFN/ PEAKSH(NB),SQH2(NB),SQH2L(NB),SQH2G(NB),ASYM(NB),
     &ETAL(NB),ETAH(NB),FNORM1(2,NB),COMLG(2,NB),COML(2,NB),COMG(2,NB),
     &COEFL(2,NB),COEFG(2,NB),DECAY(2,NB),COMPEAR(2,NB),CDX(2,NB),
     &CDA1(2,NB),CDA2(2,NB),CDW(2,NB),CDRL1(2,NB),CDRL2(2,NB),
     &CDRH1(2,NB),DMPART(2,NB)
      COMMON /UNCONS/ NUCREF,IPHASE(90),IHKL(3,90),NOPPP(90),LPPP(NB),
     &  RCUT(2,90),NCUT
      INTEGER WIDTH,HEIGHT,YMIN,YMAX,LBGPLOT,LDEL,OFFSETD,MSIZE,FSIZE,
     &  LSIZE,OFFSET,INDREF
      COMMON /IGOR/ WIDTH,HEIGHT,YMIN,YMAX,LBGPLOT,LDEL,OFFSETD,MSIZE,
     &  FSIZE,LSIZE,OFFSET(NPH),INDREF
     
*     Kb: CuK-beta
      DATA (XTARG(I),I=1,7)/'Cr','Fe','Co','Cu','Mo','Ag','Kb'/
      DATA ((WL(J,I),J=1,4),I=1,7)/
     &  2.28970,   2.293606,  1.14485,   1.146803,
     &  1.936042,  1.939980,  0.968021,  0.96999,
     &  1.788965,  1.792850,  0.8944825, 0.896425,  
     &  1.540562,  1.544390,  0.770281,  0.772195,
     &  0.7093,    0.713590,  0.35465,   0.356795,
     &  0.5594075, 0.563798,  0.2797038, 0.281899,
     &  1.39222,   1.39222,   0.69611,   0.69611/
      DATA (CRYSYS(I),I=1,15)/'triclinic, -1','monoclinic, 2/m (a-axis)'
     &,'monoclinic, 2/m (b-axis)','monoclinic, 2/m (c-axis)',
     &'orthorhombic, mmm','tetragonal, 4/m','tetragonal, 4/mmm',
     &'trigonal, -3 (rhombohedral lattice)','trigonal, -3 (hexagonal lat
     &tice)','trigonal, -3m (rhombohedral lattice)','trigonal, -3m (hexa
     &goanl lattice)','hexagonal,6/m','hexagonal, 6/mmm',
     &'cubic, m3','cubic, m3m'/
      DATA RD,RD2/0.01745329,57.29578/
C
      LD(J)=ID(KPHB(L)+J)
      P(J)=A(KPHB(L)+J)
	  
      CALL INITIAL
      WRITE(6,'(54X,A///(20X,A))') 'Welcome to the',
     &'RRRRRRRRRRR     IIIIIIII    EEEEEEEEEEEE  TTTTTTTTTTTT  '//
     &' AAAAAAAAAA    NN        NN',
     &'RRRRRRRRRRRR    IIIIIIII    EEEEEEEEEEEE  TTTTTTTTTTTT  '//
     &'AAAAAAAAAAAA   NNN       NN',
     &'RR        RR       II       EE                 TT       '//
     &'AA        AA   NNNN      NN',
     &'RR        RR       II       EE                 TT       '//
     &'AA        AA   NN NN     NN',
     &'RR        RR       II       EE                 TT       '//
     &'AA        AA   NN  NN    NN',
     &'RRRRRRRRRRRR       II       EEEEEEEEEEE        TT       '//
     &'AAAAAAAAAAAA   NN   NN   NN'
      WRITE(6,'(20X,A)')
     &'RRRRRRRRRRR        II       EEEEEEEEEEE        TT       '//
     &'AAAAAAAAAAAA   NN    NN  NN',
     &'RR    RRR          II       EE                 TT       '//
     &'AA        AA   NN     NN NN',
     &'RR     RRR         II       EE                 TT       '//
     &'AA        AA   NN      NNNN',
     &'RR      RRR        II       EE                 TT       '//
     &'AA        AA   NN       NNN',
     &'RR       RRR    IIIIIIII    EEEEEEEEEEEE       TT       '//
     &'AA        AA   NN        NN',
     &'RR        RR    IIIIIIII    EEEEEEEEEEEE       TT       '//
     &'AA        AA   NN        NN'
      WRITE(6,'(//17X,A///)')
     &'system for pattern-fitting structure refinement with '//
     &'X-ray and neutron diffraction data'
      CALL TINK
      CALL SETNUM

C**   TITLE
      READ(4,'(A80)') TITLE
      WRITE(6,'(7X,A///11X,A,A/)')
     &'*** Data given by user and system ***','Title: ',TITLE

C**   CONTROL DATA

*     NBEAM = 0: NEUTRON DIFFRACTION.
*     NBEAM = 1: X-RAY DIFFRACTION.
*     NBEAM = 2: SYNCHROTRON X-RAY DIFFRACTION.
 
*     NMODE = 0: RIETVELD ANALYSIS OF POWDER DIFFRACTION DATA.
*     NMODE = 1: CALCULATION OF POWDER DIFFRACTION INTENSITIES
*                PLUS SIMULATION OF POWDER PATTERNS FOR NPAT <> 0.
*     NMODE = 2: PROFILE FITTING OF THE WHOLE PATTERN USING STRUCTURE
*                FACTORS OBTAINED BY MEED.
*     NMODE = 3: THE SAME AS NMODE = 2 EXCEPT THAT |Fc|'S OF RELAXED 
*                REFLECTIONS ARE NOT FIXED BUT REFINED.
*     NMODE = 4: CONVENTIONAL LE BAIL REFINEMENT.
*     NMODE = 5: LE BAIL REFINEMENT WHERE INITIAL STRUCTURE FACTORS ARE
*                CALCULATED FROM STRUCTURE PARAMETERS.
*     NMODE = 6: INDIVIDUAL PROFILE FITTING.

C     NPRINT = 0: A MINIMUM PRINTOUT IS DESIRED.
C     NPRINT = 1: A STANDARD PRINTOUT IS DESIRED.
C     NPRINT = 2: A COMPLETE PRINTOUT IS DESIRED.

      READ(4,*) NBEAM,NMODE,NPRINT
      WRITE(6,120) NBEAM,NMODE,NPRINT
  120 FORMAT(/11X,'NBEAM =',I2,5X,'NMODE =',I2,5X,'NPRINT =',I2)
      IF (NBEAM .LT. 0 .OR. NBEAM .GT. 2 .OR. NMODE .LT. 0 .OR.
     &NMODE .GT. 6 .OR. NPRINT .LT. 0 .OR. NPRINT .GT. 2) THEN
         CALL JOBEND('Control parameter out of range')
      END IF
      
      SELECT CASE (NBEAM)

         CASE (0)
*           Neutron diffraction
*             XLMD: Wavelength/Angstrom.
*           RADIUS: Radius of the cylindrical container/cm.
*           DENSTY: Density/g.cm**(-3). 0 for no correction.
            READ(4,*) XLMD,RADIUS,DENSTY
            WRITE(6,'(//11X,A,F8.5,A,5X,A,F6.2,A,5X,A,F8.5,A)') 
     &      'Wavelength =',XLMD,' Angstroms','RADIUS =',RADIUS,' cm',
     &      'DENSTY =',DENSTY,' g/cm**3'
            XLMDH = 0.5*XLMD
            RADIUS = RADIUS*1.0E-2
            R12 = 0.0

         CASE (1)
*           Characteristic X rays    
*           TNAME: Target element.
*             R12: Intensity(alpha2)/Intensity(alpha1).
*            CTHM: cos(2*alpha)**2, and 2*alpha = diffraction angle in the
*                  crystal monochromator.
**          Lorentz-polarization factor = 
*           (1 - PCOR + PCOR*CTHM*cos(2*theta)**2)/(2*sin(theta)**2*cos(theta)),
*           where PCOR = 0.5 in the case of characteristic X rays.
*           Refer to Eq. (3.25) in the manual of Fullprof.
**          Flag for correction of surface roughness
*           NSURFR = 0: No correction.
*                    1: Model combining NSURFR = 2 and 3.
*                    2: Model of Sparks et al.
*                    3: Model of Suortti et al. 
*                    4: Model of Pitschke et al.
**          NTRAN = 0: Bragg-Brentano geometry (conventional divergence slit).
*                 = 1: Bragg-Brentano geometry (automatic divergence slit).
*                 = 2: Transmission geometry (e.g., Guinier diffractometer).
*                 = 3: Debye-Scherrer geometry.
            READ(4,'(A80)') LINE
            CALL BACKSP(LINE)
            IF(LAPO(LINE) .EQ. 0) CALL JOBEND('The name of target '//
     &      'element should be enclosed by a pair of apostrophes')
            READ(4,*) TNAME,R12,CTHM,NSURFR,NTRAN
            DO J = 1, 7
               IF (TNAME .EQ. XTARG(J)) THEN
                  NTARG = J
                  XLMD = WL(1,J)
                  XLMD2 = WL(2,J)
                  XLMDH = WL(3,J)
                  XLMD2H = WL(4,J)
                  RLAMBD = XLMD2/XLMD
                  EXIT
               END IF
            END DO
            IF (J. EQ. 8) CALL JOBEND('Bad target name: '//TNAME)
            WRITE(6,150) TNAME,XLMD,XLMD2,R12,CTHM,NSURFR
  150       FORMAT(//11X,'Radiation: ',A2/11X,'Wavelength(alpha1) =',
     &      F9.6,5X,'Wavelength(alpha2) =',F9.6,5X,
     &      'I(alpha2)/I(alpha1) =',F9.6,5X,'CTHM =',F9.6/
     &      11X,'NSURFR =',I2)
            IF (NSURFR .LT. 0 .OR. NSURFR .GT. 4)
     &      CALL JOBEND('Invalid NSURFR value')
            SELECT CASE (NTRAN)
               CASE (0)
                  WRITE(6,'(11X,A,I2)') 'NTRAN =',NTRAN
               CASE (1)
*                  DSANG: Divergence angle/degree (0 for a conventional
*                         divergence slit.
*                   RGON: Goniomator radius/mm.
*                 SWIDTH: Illumination length/mm.
                  READ(4,*) DSANG,RGON,SWIDTH
                  WRITE(6,'(11X,A,I2,5X,A,F8.4,A,6X,A,F6.2,A,5X,A,F5.2,
     &            A)') 'NTRAN =',NTRAN,'DSANG =',DSANG,' degrees',
     &            'RGON =',RGON,' mm','SWIDTH =',SWIDTH,' mm'
                  DSANG = DSANG*RD
               CASE (2) 
*                 PCOR: Fraction of the perfect crystal contribution.
*                 SABS: (Linear absorption coefficient)*(effective thickness).
                  READ(4,*) PCOR,SABS
                  WRITE(6,'(11X,A,I2,5X,A,F7.4,5X,A,F8.5)') 'NTRAN =',
     &            NTRAN,'PCOR =',PCOR,'SABS =',SABS
               CASE (3)
*                 MUR: (Linear absorption coefficient)*(radius).
                  READ(4,*) MUR
                  WRITE(6,'(11X,A,I2,5X,A,F9.6)') 'NTRAN =',NTRAN,
     &            'MUR =',MUR
               CASE DEFAULT
                  CALL JOBEND('0 <= NTRAN <= 3')
            END SELECT

         CASE (2)
*           Synchrotron X-ray diffraction
*           XLMD: Wavelength/Angstrom.
*           PCOR: I0(perpendicular)/I0(parallel). I0: incident intensity.
*           CTHM: cos(2*alpha)**2 for the crystal monochromator (see above).
*            MUR: (Linear absorption coefficient)*(radius)
            READ(4,*) XLMD,PCOR,CTHM,MUR
            WRITE(6,'(//11X,A,F8.5,A,5X,A,F7.4,5X,A,F9.6,5X,A,F9.6)') 
     &      'Wavelength =',XLMD,' Angstroms','PCOR =',PCOR,'CTHM =',
     &      CTHM,'MUR =',MUR
            XLMDH = 0.5*XLMD
            R12 = 0.0

         CASE DEFAULT            
            CALL JOBEND('Invalid NBEAM value')

      END SELECT

C**   CHEMICAL SYMBOLS
C     ANAME(KATOM): CHEMICAL SYMBOL (+VALENCE)
C     '*' IS ADDED IMMEDIATELY AFTER THE CHEMICAL SYMBOL OF THE MAGNETIC ATOM
      READ(4,'(A)') LINE
      CALL BACKSP(LINE)
      IF (LAPO(LINE) .EQ. 0) THEN
         CALL JOBEND('The name of pure species should be enclosed by '//
     &   'a pair of apostrophes')
      END IF
      DO J=1,NA
         ANAME(J)=' '
      END DO
      IF (NBEAM .EQ. 0) THEN
         READ(4,*) (ANAME(J),CATOM(J),J=1,NA)
      ELSE
         READ(4,*) (ANAME(J),J=1,NA)
      END IF

C     NREAL: NUMBER OF PURE CHEMICAL SPECIES
      DO J = NA, 1, -1
         IF (ANAME(J) .NE. ' ') EXIT
      END DO
      NREAL = J
      IF (NREAL .EQ. NA) CALL JOBEND
     &('The number of pure chemical species reached the maximum number')

*      CMFF: COEFFICIENTS IN ANALYTICAL APPROXIMATIOS TO MAGNETIC FORM FACTORS
*     LCMFF: L (= 0, 2, 4, OR 6) IN EQUATIONS (4.4.5.2) AND (4.4.5.3)
*     CF. INTERNATIONAL TABLES, VOL. C, PP. 391-399
*     IDREAL=0: THE DATA FOR THIS CHEMICAL SPECIES HAS NOT BEEN READ
*     IDREAL=1: THE DATA FOR THIS CHEMICAL SPECIES HAS BEEN READ
      DO J=1,NREAL
         IF(NBEAM.EQ.0 .AND. (INDEX(ANAME(J),'+').NE.0 .OR.
     &   INDEX(ANAME(J),'-').NE.0)) CALL JOBEND
     &   ('Ionic species should not be input in neutron diffraction')
         IDREAL(J)=0
*        MAGATM=0: NONMAGNETIC ATOM
*        MAGATM=1: MAGNETIC ATOM
         MAGATM(J)=0
         ANAME2(J)=ANAME(J)
         ISTAR=INDEX(ANAME(J),'*')
         IF (ISTAR .NE. 0) THEN
            ANAME2(J)(ISTAR:ISTAR)=' '
            MAGATM(J)=1
            READ(4,*) LCMFF(J),(CMFF(JJ,J),JJ=1,7)
            IF (LCMFF(J) .NE. 0 .AND. LCMFF(J) .NE. 2 .AND. 
     &      LCMFF(J) .NE. 4 .AND. LCMFF(J) .NE. 6)
     &      CALL JOBEND('Invalid LCMFF value')
         END IF
      END DO

C*    DISPERSION CORRECTION (NBEAM=2)
C     READ DISPERSION TERMS FOR SR
      IF (NBEAM .EQ. 2) READ(4,*) (DELTF1(J),DELTF2(J),J=1,NREAL)

      CALL ASFDC(NREAL,IDREAL,ANAME,ANAME2,AW,NBEAM,AA,BB,C,DC,
     &  DELTF1,DELTF2,NTARG,SIGA,SIGI,CATOM,ATOMNO,DENSTY)

*     CHECK WHETER NMIX > 0
      CALL CHKMIX(MIX)
      IF (MIX .EQ. 0) THEN
         NMIX=0
      ELSE
         SELECT CASE (NBEAM)
            CASE (0)
               WRITE(6,'(11X,91A1)') ('-',J=1,91)
            CASE DEFAULT
               WRITE(6,'(10X,122A1)') ('-',J=1,122)
         END SELECT
         WRITE(6,247)
         IMIX=1
         DO LL = NREAL+1, NA+1
            READ(4,'(A80)') LINE
            IF (LINE(1:1) .EQ. '}') EXIT
            IF (LL .EQ. NA+1) THEN
               CALL JOBEND('Too large number of chemical species')
            END IF
            CALL RDFRCT(LINE,ANAME,LL,IMIX,NREAL,AW)
            WRITE(6,255) LL,ANAME(LL),
     &      (ANAME(MFRAC(J,IMIX)),FRCTN(J,IMIX),J=1,NOFRAC(IMIX))
            IMIX=IMIX+1
         END DO
         NMIX=LL-NREAL-1
      END IF
  247 FORMAT(' ',10X,'Virtual atom, constituents, and their fractions')
  255 FORMAT(10X,I5,3X,A5,5X,A5,F11.6,3(6X,A5,F11.6))
C     NATOM: NUMBER OF ATOMIC SPECIES (INCLUDING VIRTUAL SPECIES)
      NATOM=NREAL+NMIX
      IF (NATOM .GT. NA) THEN
         CALL JOBEND('Number of atomic species out of range')
      END IF
 
      LL=0
      DO J = 1, NREAL
         IF (MAGATM(J) .NE. 1) CYCLE
         IF (LL .EQ. 0) WRITE(6,'(//11X,A/)') 
     &   'Coefficients for analytical approximations to magnetic '//
     &   'form factors'
         WRITE(6,'(13X,A,A,I2,A,7F8.4)') ANAME(J),'  L =',LCMFF(J),
     &   '  CMFF: ',(CMFF(JJ,J),JJ=1,7)
         LL=1
      END DO

C**   INFORMATION ON EACH PHASE
C       PHNAME: NAME OF A PHASE
C          VNS: VOLUME NAME-SPACE GROUP #-SETTING # OF INTERNATIONAL
C               TABLES
C     LSPSYM=0: STANDARD SYMMETRY IS REQUIRED
C     LSPSYM=1: SOME SPECIAL SYMMETRY IS REQUIRED
C     LSPSYM=2: NON-STANDARD SYMMETRY IS REQUIRED
*Rev 1.0h 2000.12.13 Izumi {
C     LPAIR=0: No Bijvoet pairs (hkl & -h-k-l) are generated.
C     LPAIR=1: Bijvoet pairs (hkl & -h-k-l) are generated.
* }
C      INDIV=0: NEITHER ISOTROPIC NOR ANISOTROPIC DISPLACEMENT 
C               PARAMETERS ARE INPUT FOR INDIVIDUAL SITES (THE OVERALL 
C               ISOTROPIC DISPLACEMENT PARAMETER IS USED)
C      INDIV=1: ISOTROPIC AND/OR ANISOTROPIC DISPLACEMENT PARAMETERS 
C               ARE INPUT
C*    PREFERRED ORIENTATION
C     NPROR=0: DO NOT CORRECT PREFERRED ORIENTATION (ENTER ONLY '/')
C     NPROR=1: PLATY CRYSTALS (TORAYA-MARUMO FUNCTION)
C     NPROR=2: NEEDLE CRYSTALS (TORAYA-MARUMO FUNCTION)
C     NPROR=3: MARCH-DOLLASE FUNCTION
C     (HP, KP, LP): PREFERRED-ORIENTATION VECTOR
C     WHEN PREFERRED ORIENTATION IS NEGLIGIBLE, NPROR IS SET TO 0.
C     THE PREFERRED-ORIENTATION VECTOR, (HP)A* + (KP)B* + (LP)C*, IS
C     PERPENDICULAR TO THE (HP KP LP) PLANE.
C     THE DIRECTION OF THE VECTOR CORRESPONDS TO THE NORMAL TO THE
C     PLATE CRYSTALLITES OR TO THE NEEDLE AXIS IN NEEDLE CRYSTALS.
C     ANISOTROPIC LORENTZIAN BROADENING
C     (HANIS, KANIS, LANIS): ANISOTROPIC BROADENING AXIS.
 
      DO L = 1, NPH+1
         READ(4,'(A)') LINE
         IF (LINE(1:1) .EQ. '}') EXIT
         IF (L .EQ. NPH+1) THEN
            CALL JOBEND('Too large number of phases')
         END IF
         CALL BACKSP(LINE)
*Rev 1.0h 2000.12.13 Izumi {
         IF (NBEAM .EQ. 0) THEN
*           No Bijvoet pairs are generated in neutron diffraction
            LPAIR(L) = 0
            READ(4,*) PHNAME(L),VNS(L),LSPSYM(L),INDIV(L),NPROR(L),
     &      HP(L),KP(L),LP(L),LSUM(L),HANIS(L),KANIS(L),LANIS(L)
         ELSE
            READ(4,*) PHNAME(L),VNS(L),LSPSYM(L),LPAIR(L),INDIV(L),
     &      NPROR(L),HP(L),KP(L),LP(L),LSUM(L),HANIS(L),KANIS(L),
     &      LANIS(L)
         END IF
* }
         IF (HP(L) .EQ. 0 .AND. KP(L) .EQ. 0 .AND. LP(L) .EQ. 0)
     &   CALL JOBEND('Invalid preferred-orientation vector: (0,0,0)')
         IF (HANIS(L) .EQ. 0 .AND. KANIS(L) .EQ. 0 .AND. 
     &   LANIS(L) .EQ. 0) CALL JOBEND
     &   ('Invalid anisotropic-broadening axis: (0,0,0)')

C*       TYPE OF POSSIBLE REFLECTIONS
         IF (LSPSYM(L) .EQ. 1) CALL TYPEREF(L,REFTYP,LCON)

         CALL INTTBL(VNS(L),L,NSPGR,NSET,LAUEG,NCENTR,NSYM,SPGR,LSPSYM,
     &   NSSYM,SYM,LINNSS,NLNSS)

         WRITE(6,'(//11X,A,I1,A,A)') 'Phase #',L,': ',PHNAME(L)
         IF(NSET(L) .EQ. 1) THEN
            FSSET=')'
         ELSE
            WRITE(FSSET,'(A,I1,A)') ', Setting #',NSET(L),')'
         END IF
         WRITE(6,'(11X,5A,I4,A)')
     &   'Space group: ',SPGR(L),' (VOL. ',VNS(L)(1:1),',',NSPGR(L),
     &   FSSET
         WRITE(6,277) CRYSYS(LAUEG(L))
  277    FORMAT(11X,'Crystal system and Laue-symmetry class: ',A)
*Rev 1.0h 2000.12.13 Izumi {
         IF (NBEAM .GE. 1) THEN
            WRITE(6,'(11X,A,I2)') 'LPAIR =',LPAIR(L)
            IF (LPAIR(L) .EQ. 1 .AND. NCENTR(L) .EQ. 1) CALL JOBEND
     &      ('No Bijvoet pairs should be generated for a space group '//
     &      'with an inversion center at the origin')
         END IF
* }
         IF (NPROR(L) .NE. 0) THEN
            WRITE(6,278) NPROR(L),HP(L),KP(L),LP(L),LSUM(l)
         ELSE
            WRITE(6,278) 0
         END IF
         IF (NPROR(L) .LT. 0 .OR. NPROR(L) .GT. 3) THEN
            WRITE(6,'(//11X,A,I2)') 'NPROR out of range: ',NPROR(L)
            CALL JOBEND(' ')
         ELSE IF (NPROR(L) .EQ. 1 .AND. NBEAM .EQ. 0) THEN
*           NEUTRON DIFFRACTION, CYLINDRICAL CONTAINER
            NPROR(L)=2
         END IF
         WRITE(6,279) HANIS(L),KANIS(L),LANIS(L)
         IF (LSPSYM(L) .EQ. 1) THEN
            WRITE(6,'(/11X,A)') 'Special symmetry conditions'
            DO I = 1, MLINE(L)
               IF (I .EQ. 1) THEN
                  WRITE(6,'(13X,3A)') REFTYP(L),': ',LCON(I,L)
               ELSE
                  WRITE(6,'(19X,A)') LCON(I,L)
               END IF
            END DO
         ELSE IF (LSPSYM(L) .EQ. 2) THEN
            WRITE(6,'(/11X,A,I3/11X,A,I2/11X,A)') 'LAUEG  =',LAUEG(L),
     &      'NCENTR =',NCENTR(L),'Non-standard symmetry conditions'
            WRITE(6,'(13X,A)') (LINNSS(J),J=1,NLNSS)
         END IF
*Rev 1.0j 2000.12.21 Izumi {
         SELECT CASE (NCENTR(L))
            CASE (1) 
               WRITE(6,280)
            CASE DEFAULT
               WRITE(6,285)
         END SELECT
         DO J=1,NSYM(L)
            WRITE(6,290) (SYM(KK,J,L),KK=1,3)
         END DO
* }
      END DO

      NPHASE=L-1
  278 FORMAT(/11X,'Preferred-orientation correction'/13X,'NPROR =',I2,:,
     &5X,'Preferred-orientation vector: (',I2,',',I2,',',I2,')',5X,
     &'LSUM =',I2)
  279 FORMAT(/11X,'Anisotropic-broadening axis: (',I2,',',I2,',',I2,')')
  280 FORMAT(/11X,'Coordinates of equivalent positions exclusing those g
     &enerated by lattice centering')
  285 FORMAT(/11X,'Coordinates of equivalent positions')
  290 FORMAT(' ',12X,3(A10,3X))
 
C     NPRFN = 0: PSEUDO-VOIGT FUNCTION OF THOMPSON, COX, AND HASTINGS
C     NPRFN = 1: SPLIT-TYPE PSEUDO-VOIGT FUNCTION DESCRIBED IN:
C                H. Toraya, J. Appl. Crystallogr. 23 (1990) 485.
C     NPRFN = 2: SPLIT-TYPE PSEUDO-VOIGT FUNCTION PLUS MODIFIED
C                SPLIT-TYPE PSEUDO-VOIGT FUCNTION (W1 <> W2) FOR
C                RELAXED REFLECTIONS
C     NPRFN = 3: SPLIT-TYPE PEARSON VII FUNCTION

*     NASYM = 0: The pseudo-Voigt function of Thompson, Cox, and
*                Hastings are made asymmetric according to a procudure 
*                of Finger, Cox, and Jephcoat (1994).
*     NASYM = 1: The pseudo-Voigt function of Thompson, Cox, and
*                Hastings are made asymmetric according to a procudure  
*                of Howard (1982).

*     NSHIFT = 0: Peak shift = t0
*     NSHIFT = 1: Peak shift = t0 + t1*cos(2th) + t2*sin(2th) + t3*tan(th)
*     NSHIFT = 2: Peak shift = t0 + t1*(2th) + t2*(2th)**2 + t3*(2th)**3
*     NSHIFT = 3: Peak shift = t0 + t1*tan(th) + t2*(tan(th))**2 +
*                              t3*(tan(th))**3
*     NSHIFT = 4: Peak shift = Legendre polynomial of normalized 2th 
*     NSHIFT = 5: Peak shift = Legendre polynomial of normalized tan(th)

      READ(4,*) NPRFN
      WRITE(6,'(//11X,A,I2)') 'NPRFN =',NPRFN
      IF (NPRFN .LT. 0 .OR. NPRFN .GT. 3) THEN
         CALL JOBEND('NPRFN out of range')
      ELSE IF (NPRFN .EQ. 0) THEN
         READ(4,*) NASYM
         WRITE(6,'(11X,A,I2)') 'NASYM =',NASYM
         IF (NASYM .LT. 0 .OR. NASYM .GT. 1) 
     &   CALL JOBEND('Invalid NASYM value')
         IF (NMODE .EQ. 6) CALL JOBEND('NPRFN should be greater than '//
     &   'zero when NMODE = 6')
      ELSE IF (NPRFN .GE. 1) THEN
         READ(4,*) NSHIFT
         WRITE(6,'(11X,A,I2)') 'NSHIFT =',NSHIFT
         IF (NSHIFT .LT. 0 .OR. NSHIFT .GT. 5) 
     &   CALL JOBEND('Invalid NSHIFT value')
      END IF

      CALL MLAT
      CALL RDPARA(NRFN,NCNSTR,NTERMS,A,ID,AOUT,IDSAVE,IR,
     &  NMODE,FIRST,NSPLBL,LABEL,LPAR,NLINE)
C
      IP=1
      DO L=1,NPHASE
         CALL FLSITE(A,L,NOAT,INDIV,LISO,MAGATM,NSITE,LMAG,
     &   NATOM,ISONUM,NLINE,LPAR,LABEL,IP,ANAME)
      END DO
      IF (NMODE .GE. 2) CALL CHGID(NTERMS,NCNSTR)
      
      IF (NTERMS .LT. KPHE(NPHASE)) THEN
         CALL JOBEND('The number of parameters input is too small')
      ELSE IF (NTERMS .GT. KPHE(NPHASE)) THEN
         CALL JOBEND('The number of parameters input is too large')
      END IF

      CALL CCTENS(A)
      CALL PARAM(ANAME,ISONUM,A,ID,NMODE,NTERMS,LABEL,LPAR,NLINE)
      CALL ACONV(NPHASE,LMAG,A,KPHB,KPHE,LAUEG)
 
C     PRINT OUT DIRECT AND RECIPROCAL CELL CONSTANTS
      IF (NPRINT .GE. 1) CALL PRUCC
C
      CALL CONSTR(NCNSTR,NTERMS,LABEL,LPAR,NLINE)
      CALL SETPAR(A)
C
      CALL SYMOP(A)
*     NDA MUST BE DEFINED BECAUSE IT IS USED IN SUBROUTINE COMPOS
      NDA=0
      IF (NPRINT .NE. 0) THEN
         CALL EQPOS(NPHASE,PHNAME,NSITE,ISONUM,ANAME,NOAT,NSYM,IDSYM,
     &   NCENTR,IFAC,SYM,NSSYM)
         CALL APRNTD(A,NTERMS)
         CALL COMPOS(A,SIGMAA,0,ANAME,INDIV,AW)
      END IF
C
      DO L = 1, NPHASE
         IF (P(NPO1X) .NE. 0.0 .OR. P(NPO2X) .NE. 0.0 .OR. 
     &   (NMODE .NE. 1 .AND. LD(NPO1X) .NE. 0) .OR. 
     &   (NMODE .NE. 1 .AND. LD(NPO2X) .NE. 0)) THEN
*           PREFERRED ORIENTATION SHOULD BE CORRECTED
            SELECT CASE (NPROR(L))
               CASE (0)
                  CALL JOBEND
     &            ('NPROR for this phase should be 1, 2, or 3')
               CASE (3)
                  IF (P(NPO2X) .NE. 0) CALL JOBEND
     &            ('Parameter p2 should be 0 when using the '//
     &            'March-Dollase function')
            END SELECT
         END IF
         IF (NMODE .GE. 4 .AND. (LD(NPO1X) .NE. 0 .OR. LD(NPO2X)  
     &   .NE. 0)) THEN
            CALL JOBEND('Neither p1 nor p2 should be refined when '//
     &      'NMODE >= 4')
         END IF
      END DO
C
*     NCUT = 0: 2-Theta ranges for relaxed reflections are calculated
*               by RIETAN
*     NCUT = 1: 2-Theta ranges for relaxed reflections are input by
*               the user
      READ(4,*) NCUT
      WRITE(6,'(//11X,A,I2)') 'NCUT =',NCUT
      IF (NCUT .LT. 0 .OR. NCUT .GT. 1) THEN
         CALL JOBEND('Invalid NCUT value')
      ELSE IF (NPRFN .EQ. 0 .AND. NCUT .EQ. 1) THEN
         CALL JOBEND('NCUT should be zero when NPRFN = 0') 
      ELSE IF (NCUT .EQ. 1) THEN
*        READ 2-THETA RANGES FOR RELAXED REFLECTIONS
         WRITE(6,'(//11X,A)') '2-Theta ranges for relaxed reflections'
         DO J = 1, NUCREF
            READ(4,*) RCUT(1,J),RCUT(2,J)
            WRITE(6,'(13X,F7.3,A,F7.3)') RCUT(1,J),'-',RCUT(2,J)
            RCUT(1,J) = RD*RCUT(1,J)
            RCUT(2,J) = RD*RCUT(2,J)
         END DO
      END IF

      IF (NMODE .NE. 1) THEN
C        NEXC = 0: USE ALL THE INTENSITY DATA.
C        NEXC = 1: SKIP SOME INTENSITY DATA.
         READ(4,*) NEXC
C
C*       EXCLUDED REGIONS (NEXC=1)
C        DEGEXC(1,J): LOWER LIMIT, DEGEXC(2,J): UPPER LIMIT
C*       END OF REGION: '}'
C        NSKIP: NUMBER OF SKIPPED REGIONS
         IF (NEXC .EQ. 1) THEN
            DO J = 1, 50
               READ(4,'(A80)') LINE
               IF (LINE(1:1) .EQ. '}') EXIT
               CALL BACKSP(LINE)
               READ(4,*) (DEGEXC(MM,J),MM=1,2)
            END DO
            IF (J .EQ. 51) CALL JOBEND('The number of excluded regions '
     &      //'has reached the maximum value')
            NSKIP=J-1
            WRITE(6,480) ((DEGEXC(MM,J),MM=1,2),J=1,NSKIP)
         ENDIF
C
C        READ STEP-SCANNED INTENSITY DATA RECORDED ON DISK
         CALL RDINT(THINIT,STEP,NEXC,NSTEP,DEG,XINT)
         IF (NMODE .EQ. 2 .OR. NMODE .EQ. 3) CALL RDMEED
      ENDIF
  480 FORMAT(//11X,'Excluded regions in the diffraction data'/
     &(13X,F7.3,'-',F7.3))
C
      IF (NMODE .NE. 1) THEN
C        NRANGE=0: REFINE THE BACKGROUND.
C        NRANGE=1: FIX BACKGROUND INTENSITIES AT INTERPOLATED VALUES.
C        NRANGE=2: FIX BACKGROUND INTENSITIES AT INPUT VALUES.
C        NRANGE=3: BACKGROUND FUNCTION = (BACKGROUND INTENSITIES INPUT 
C                  FROM *.bkg) * (LEGENDRE POLYNOMIALS).
         READ(4,*) NRANGE
         IF (NRANGE .EQ. 1) THEN
C*          BACKGROUND SPECIFICATION (NRANGE=1)
C           X: ARRAY OF TWO-THETA (PLACE '/' AFTER THE LAST VALUE)
            DO J=1,100
               X(J)=-1.0
            END DO
            CALL OPENFILE(3,FILE8)
            READ(8,*) (X(J), Y(J), J = 1, 100)
            IF (X(100) .GT. 0.0) THEN
               CALL JOBEND('The number of background specifications has'
     &         //'reached the maximum value')
            END IF
            DO J = 1, 100
               IF (X(J) .LT. 0.0) EXIT
            END DO
            NBG=J-1
            CALL SMINT(X,Y,DEG,XINT,NBG,NSTEP)
         ELSE IF (NRANGE .GE. 2) THEN
            CALL RDBKGD(NSTEP,NEXC)
         ENDIF
      END IF

C
C*    RANGE OF TWO-THETA
C     REFLECTIONS BETWEEN DEG1 AND DEG2 ARE CALCULATED

C     USTP: STEP WIDTH

C     NPAT = 0: DO NOT CREATE ANY FILE STORING DIFFRACTION INTENSITIES.
C     NPAT = 1: CREATE A PostScript FILE STORING STRUCTURE-REFINEMENT PATTERNS
C               (NMODE <> 1) OR SIMULATED PATTERNS (NMODE = 1).
C     NPAT = 2: CREATE A Macplot/RietPlot-FORMAT FILE STORING STRUCTURE-REFINEMENT 
C               PATTERNS (NMODE <> 1) OR SIMULATED PATTERNS (NMODE = 1).
C     NPAT = 3: CREATE A PLOT-FORMAT FILE STORING STRUCTURE-REFINEMENT PATTERNS.
C     NPAT = 4: CREATE A SigmaPlot-FORMAT FILE STORING STRUCTURE-REFINEMENT 
C               PATTERNS.
C     NPAT = 5: CREATE An Igor-FORMAT FILE STORING STRUCTURE-REFINEMENT PATTERNS.
 
      IF (NMODE .EQ. 1) THEN
*Rev 1.03 2000.12.26 Izumi
         NRANGE = 0
         READ(4,*) DEG1,DEG2,USTP,NPAT
         DEG1=DEG1*RD
         DEG2=DEG2*RD
*          WIDTH: Width of the graph.
*         HEIGHT: Height of the graph.
*           YMIN: Minimum value for the y axis (zero for the default value).
*           YMAX: Maximum value for the y axis (zero for the default value).
*        LBGPLOT: = 0, Do not plot the background.
*                 = 1, Plot the background.
*           LDEL: = 0, Plot delta_y (observed minus calculated intensities).
*                 = 1, Plot SQRT(wi)*delta_y.
*                 = 2, Plot [delta_y/(observed intensity)]/sigma_i.
*                   For LDEL = 2, refer to R. A. Young, "The Rietveld Method,"
*                   Oxford University Press, Oxford (1993), p. 24.
*        OFFSETD: Offset of delta_y or SQRT(wi)*delta_y.
*          MSIZE: Length of tick marks.
*          FSIZE: Size of numerical values attached to axes.
*          LSIZE: Size of axis labels.
*     INDREF = 0: Output no profiles for individual reflections
*            = 1: Output profiles for individual reflections
*        OFFSETn: Offset of tick marks along the y axis for the n'th phase.
         SELECT CASE (NPAT)
            CASE (1) 
               CALL JOBEND('NPAT = 1 is no longer supported')
            CASE (2:)
C              STEP: STEP INTERVAL
               IF (USTP .EQ. 0.0) USTP = 0.01
               STEP = USTP*RD
               NSTEP=(DEG2-DEG1)/STEP+1
               IF (NSTEP .GT. NP-2) THEN
                  NSTEP=NP-2
                  STEP=(DEG2-DEG1)/FLOAT(NSTEP-1)
               END IF
               IF (NPAT .EQ. 5)
     &         READ(4,*) WIDTH,HEIGHT,LBGPLOT,MSIZE,FSIZE,LSIZE
         END SELECT
      ELSE
         DEG1=DEG(1)
         DEG2=DEG(NSTEP)
         READ(4,*) NPAT
         IF (NPAT .EQ. 5) THEN
            DO J = 1, NPH
               OFFSET(J) = 0
            END DO
            READ(4,*) WIDTH,HEIGHT,YMIN,YMAX,LBGPLOT,LDEL,OFFSETD,MSIZE,
     &      FSIZE,LSIZE,INDREF,(OFFSET(J),J=1,NPH)
         END IF
      END IF

C     COEFFICIENTS FOR THE NORMALIZATION OF TWO-THETA
      FT1=2.0/(DEG2-DEG1)
      FT2=-0.5*(DEG1+DEG2)
      WRITE(6,570) DEG1*RD2,DEG2*RD2
  570 FORMAT(//11X,'2-theta(min) =',F8.3,5X,'2-theta(max) =',F8.3)
      IF (DEG1 .LT. 0.1*RD .OR. DEG2 .LT. 0.1*RD .OR. DEG1 .GT. 180.0*RD
     &  .OR. DEG2 .GT. 180.0*RD .OR. DEG2 .LE. DEG1)
     &   CALL JOBEND('Bad diffraction angle(s)')
      WRITE(6,'(//11X,A,I2)') 'NPAT =',NPAT
      IF (NMODE .NE. 1 .OR. NPAT .GT. 0) 
     &WRITE(6,590) STEP*RD2,NSTEP
  590 FORMAT(' ',10X,'STEP =',F6.3,5X,'NSTEP =',I5)
      IF ((NMODE .NE. 1 .OR. NPAT .GT. 0) .AND. (NRANGE .EQ. 0 .OR. 
     &NRANGE .EQ. 3)) THEN
         DO J=1,NSTEP
            IF (NMODE .EQ. 1 .AND. NPAT .GT. 0) 
     &      DEG(J)=DEG1+STEP*FLOAT(J-1)
            DEGNOR(J) = FT1*(DEG(J) + FT2)
         END DO
      END IF

      CALL ABSVARL(NSTEP,SIGA,SIGI,ATOMNO,NREAL,DEG,RADIUS,MUR,SABS,
     &  DSANG,RGON,SWIDTH)

      IF (NRANGE .EQ. 1) THEN
         OPEN(UNIT=40,ACCESS='SEQUENTIAL',
     &   FORM='UNFORMATTED',STATUS='SCRATCH')
         IF (VARINV(1) .EQ. 0.0) THEN
            WRITE(40) (DEG(J),XINT(J),BG(J),DEGNOR(J),ABSORP(J),
     &      VARLEN(J),ROUGH(J),J = 1, NSTEP)
         ELSE
            WRITE(40) (DEG(J),XINT(J),VARINV(J),BG(J),DEGNOR(J),
     &      ABSORP(J),VARLEN(J),ROUGH(J), J = 1, NSTEP)  
         END IF
         REWIND 40
      END IF
C
      IF (NMODE .NE. 1 .OR. NPAT .GT. 0) THEN
C*       PROFILE PARAMETERS 
C        PROFILE CUTOFF
         READ(4,*) PC
         IF (PC .LT. 1.0) THEN
            WRITE(6,'(11X,A,F7.5,A)') 'Cut-off = ',PC,
     &      '*(PEAK INTENSITY)'
         ELSE
            WRITE(6,'(11X,A,F5.2,A)') 'Cut-off = +/-',PC,'*FWHM'
            IF (NPRFN .EQ. 0) CALL JOBEND('PC should be less than 1 '//
     &      'for TCH''s pseudo-Voigt function')
         END IF 
         IF (NMODE .NE. 6) THEN
            CALL INCDEG(A,ANGMAX)
            ANGMIN=RD2*DEG1
         ELSE
            CALL ANGRANGE(A,ANGMAX,ANGMIN)
         END IF
      ELSE
         ANGMAX=RD2*DEG2
         ANGMIN=RD2*DEG1
      ENDIF

*     Initial values of multiplicity X |Fc|**2 for the 1st phase are
*     NSFF = 0: estimated according to the Wilson statistics.
*     NSFF = 1: read in from *.ffi.
*     NSFF = 2: all set at 100.0.

*     INCMULT = 0: Multiplicities are not included in integrated intensities.
*     INCMULT = 1: Multiplicities are included in integrated intensities.

*Rev 1.2 2003.05.28 {
*     NCONST = 0: Integrated intensities are varied during least-squares fitting.
*     NCONST = 1: Integrated intensities remain constant during least-squares fitting.
* }

*Rev 1.1 2003.05.17
      CHGPC = 1.0        
*Rev 1.2 2003.05.28
      NCONST = 0
      IF (NMODE .EQ. 4) THEN
         READ(4,*) NSFF
         WRITE(6,'(//11X,A,I2)') 'NSFF =',NSFF
         INCMULT = 0
*Rev 1.2 2003.05.28 {
         SELECT CASE (NSFF)
            CASE (0,2)
               READ(4,*) INCMULT, CHGPC
               WRITE(6,'(11X,A,I2,5X,A,F7.3)') 'INCMULT =',INCMULT,
     &         'CHGPC =',CHGPC
               IF (INCMULT .LT. 0 .OR. INCMULT .GT. 1)
     &         CALL JOBEND('Invalid INCMULT value')
            CASE (1)
               READ(4,*) NCONST
               WRITE(6,'(11X,A,I2)') 'NCONST =',NCONST
               IF (NCONST .LT. 0 .OR. NCONST .GT. 1)
     &         CALL JOBEND('Invalid NCONST value')
*              File *.ffi should be given for phase #1
               CALL RDSFF
            CASE DEFAULT
               CALL JOBEND('Invalid NSFF value')
         END SELECT
      END IF
* }
      CALL GENREF(A,NSPGR)

      CALL UPDATE(A,IDPOS,0)

      IF (IDPOS .EQ. -1) THEN
         CALL JOBEND('An invalid FWHM, A, or eta resulted from the '//
     &   'current profile parameters')
      ELSE IF (IDPOS .EQ. 0) THEN
         CALL JOBEND('d spacings cannot be calculated from the lattice '
     &   //'constants')
      END IF
      IF (NMODE .EQ. 4 .AND. NSFF .EQ. 0) CALL INITINT(A,DEG,XINT,NSTEP)

      IF ((NMODE .EQ. 1 .AND. NPAT .EQ. 0) .OR. 
     &(NMODE .EQ. 1 .AND. NPAT .GT. 0 .AND. NPRINT .GE. 1) .OR.
     &(NMODE .NE. 1 .AND. NPRINT .EQ. 2)) THEN
         CALL UPDATE(A,IDPOS,1)
         WRITE(6,630)
         CALL RFLCTN(0)
      ENDIF
  630 FORMAT(///7X,'*** Summary of possible reflections (based on the us
     &er-supplied parameters) ***'/)
      IF (NMODE .EQ. 1 .AND. NPAT .GT. 0) THEN
C        CALCULATE INTENSITIES FROM PROFILE AND STRUCTURE PARAMETERS
         DO J = 1, NSTEP
            YFIT(J)=CALINT(DEG,J,A)
         END DO
         IF (NPAT .GE. 2) THEN
            CALL SIMDAT(YFIT,NSTEP,TITLE,DEG,STEP,NREFL)
         END IF
      ENDIF
      IF (NMODE .EQ. 1) THEN
         CLOSE(UNIT=4)
         CALL JOBEND(' ')
      ENDIF
 
*     NLESQ=0: MARQUARDT'S METHOD (USING DERIVATIVES)
*     NLESQ=1: GAUSS-NEWTON METHOD (USING DERIVATIVES)
*     NLESQ=2: CONJUGATE DIRECTION METHOD (WITHOUT USING DERIVATIVES)
*      NESD=0: STANDARD DEVIATIONS ARE ESTIMATED BY THE CONVENTIONAL
*              METHOD
*      NSED=1: STANDARD DEVIATIONS ARE ESTIMATED BY SCOTT'S METHOD
      READ(4,*) NLESQ, NESD
      IF (NLESQ .LT. 0 .OR. NLESQ .GE. 3) THEN
         CALL JOBEND('Bad NLESQ value')
      ELSE IF (NESD .LT. 0 .OR. NESD .GT. 1) THEN
         CALL JOBEND('Bad NESD value')
      END IF
	  
*     CONTROL DATA FOR NONLINEAR LEAST-SQUARES FITTING USING PARTIAL
*     DERIVATIVES (NLESQ=0 OR 1)
*      NAUTO=0: REFINE VARIABLE PARAMETERS SIMULTANEOUSLY
*      NAUTO=1: REFINE VARIABLE PARAMATERS INCREMENTALLY
*      NAUTO=2: PARAMETERS TO BE REFINED IN EACH CYCLE ARE DETERMINED
*               AUTOMATICALLY
*      NAUTO=3: NAUTO=2 PLUS CHECK THE RESULTING PARAMETERS BY
*               POWELL'S METHOD
*        NCYCL: MAXIMUM NUMBER OF ITERATIONS
* CONV & NCONV: IF (S(N-1)-S(N))/S(N) <= CONV NCONV TIMES IN SUCCESSION,
*               FINISH THE REFINEMENT (NC=0) OR GO TO THE NEXT REFINEMENT
*               STAGE IN WHICH A LARGER VALUE OF TK IS USED (NC>0).  IN
*               THE ABOVE FORMULA, S DENOTES THE OBJECTIVE FUNCTION.

*     CONTROL DATA FOR THE CONJUGATE DIRECTION METHOD
*     MITER: MAXIMUM NUMBER OF ITERATIONS
*     STEP2: THE INITIAL STEP SIZE
*     ACC:   THE REQUIRED ACCURACY IN THE FUNCTION AND VECTOR VALUES
*     TK:    PENALTY PARAMETER (WEIGHT IMPOSED ON THE SUM OF VIOLATED
*            CONSTRAINTS
 
*      NC = 0: NO NONLINEAR CONSTRAINS ARE IMPOSED 
*      NC = 1: NONLINEAR CONSTRAINTS ARE IMPOSED 
*          TK: INITIAL VALUE OF THE PENALTY PARAMETER (WEIGHT IMPOSED 
*              ON THE SUM OF VIOLATED CONSTRAINTS)
*        FINC: TK(K) = FINC*TK(K-1)
 
*     LALPHA=0: THE COEFFICIENT MATRIX HAS NOT BEEN CALCULATED
*     LALPHA=1: THE COEFFICIENT MATRIX HAS ALREADY BEEN CALCULATED
*     IDCYCL=0: SUBROUTINE CURFIT IS CALLED TO CALCULATE STANDARD
*               DEVIAIONS OF PARAMETERS
*     IDCYCL=1: SUBROUTINE CURFIT IS CALLED TO REFINE THE VALUES OF
*               PARAMETERS

      MITER = 5
      INC=0

      LALPHA=0
      IF (NLESQ .LE. 1) THEN
         READ(4,*) NAUTO,NCYCL,CONV,NCONV,NC,TK,FINC
         IF (NAUTO .EQ. 3 .AND. (NMODE .EQ. 4 .OR. NMODE .EQ. 5)) THEN
            CALL JOBEND('The Le Bail method should not be used when '//
     &      'NAUTO = 3')
         ELSE IF (NAUTO .EQ. 1) THEN
            CALL PARNUM(NTERMS,ID,LABEL,LPAR,NLINE,IPAR,NPAR,NCHNG)
         ELSE IF (NAUTO .GE. 2) THEN
            CALL DETREF(IPAR,NPAR,NCHNG,ID)
C*          CONDITIONS OF THE CONJUGATE DIRECTION METHOD
            IF (NAUTO .EQ. 3) READ(4,*) MITER,STEP2,ACC
         END IF
         IF (NC .EQ. 1) THEN
            CALL INPCON(NC,NSCONS,EXPCTD,DEVDA,XYZC,NPHCON,IXCON)
         ELSE IF (NC .NE. 0) THEN
            CALL JOBEND('NC must be 0 or 1')
         END IF
         IF (NCYCL .GT. 0) THEN
            IDCYCL=1
            CALL MARGAU(NCYCL,NTERMS,NSTEP,NPTS,A,DEG,XINT,YFIT,
     &      SIGMAA,NESD)
            LALPHA=1
         ELSE
            LALPHA=0
         END IF
         IF (NAUTO .EQ. 3) THEN
C           CHECK THE CONVERGENCE TO THE GLOBAL MINIMUM
            NLESQ=2
            CALL CONDIR(NTERMS,NSTEP,A,DEG,XINT,YFIT)
            LALPHA=0
         ENDIF
C        STANDARD DEVIATIONS OF PARAMETERS ARE CALCULATED
         NLESQ=1
         NCYCL=1
         NC=0
         NAUTO=0
         IDCYCL=0
         CALL MARGAU(NCYCL,NTERMS,NSTEP,NPTS,A,DEG,XINT,YFIT,SIGMAA,
     &   NESD)   
      ELSE
         READ(4,*) MITER,STEP2,ACC,NC,TK
         IF (NC .EQ. 1) CALL INPCON(NC,NSCONS,EXPCTD,DEVDA,XYZC,NPHCON,
     &   IXCON)
         CALL CONDIR(NTERMS,NSTEP,A,DEG,XINT,YFIT)
C        STANDARD DEVIATIONS OF PARAMETERS ARE CALCULATED
         NLESQ=1
         NAUTO=0
         NCYCL=1
         NC=0
         IDCYCL=0
         CALL MARGAU(NCYCL,NTERMS,NSTEP,NPTS,A,DEG,XINT,YFIT,
     &   SIGMAA,NESD)
      END IF
      CLOSE(UNIT=10,STATUS='DELETE')

C     NUPDT = 0: DO NOT UPDATE PARAMETER VALUES IN THE INPUT FILE.
C     NUPDT = 1: UPDATE VARIABLE PARAMETERS IN THE PACKING MODE.

C       NFR = 0: NO FOURIER FILE IS CREATED.
C       NFR > 0: OUTPUT H, K, L, SIN(THETA/LAMBDA), FO, FC, COS(PHI),
C                SIN(PHI), FO-FC FOR THE NFR'TH PHASE.

C      NMEM = 0: NO FILE FOR MEM IS CREATED.
C      NMEM > 0: A FILE FOR MEM IS CREATED.

C       NDA = 0: NEITHER BOND DISTANCES NOR ANGLES ARE CALCULATED.
C       NDA > 0: BOND DISTANCES AND/OR ANGLES FOR THE NDA'TH PHASE
C              ARE CALCULATED.
C
      READ(4,*) NUPDT, NFR, NMEM, NDA
      WRITE(6,'(//11X,4(A,I2,5X))') 'NUPDT =', NUPDT, 'NFR =', NFR,
     &'NMEM =', NMEM, 'NDA =', NDA
      IF (NUPDT .LT. 0) THEN
         CALL JOBEND('NUPDT is negative')
      ELSE IF (NUPDT .EQ. 2) THEN
         CALL JOBEND('NUPDT = 0 or 1 in the current version')
      ELSE IF (NUPDT .GT. 2) THEN
         CALL JOBEND('NUPDT is too large')
      ELSE IF (NUPDT .NE. 0 .AND. (NCYCL .EQ. 0 .OR. MITER .EQ. 0)) THEN
         CALL JOBEND('Parameters should not be updated if NCYCL = 0 '//
     &   'or MITER = 0')
      ELSE IF (NFR .LT. 0) THEN
         CALL JOBEND('NFR is negative')
      ELSE IF (NMEM .LT. 0) THEN
         CALL JOBEND('NMEM is negative')
      ELSE IF (NDA .LT. 0) THEN
         CALL JOBEND('NDA is negative')
      ELSE IF (NFR .GT. NPHASE) THEN
         CALL JOBEND('NFR is not compatible with the number of phases')
      ELSE IF (NMEM .GT. NPHASE) THEN
         CALL JOBEND('NMEM is not compatible with the number of phases')
      ELSE IF (NDA .GT. NPHASE) THEN
         CALL JOBEND('NDA is not compatible with the number of phases')
      ELSE IF ((NMODE .EQ. 2 .OR. NMODE .EQ. 3) .AND. NMEM .NE. 1) THEN
         CALL JOBEND('MEM analysis can be applied to only the 1st '//
     &   'phase')
      ELSE IF (NMODE .GE. 4 .AND. (NFR .NE. 0 .OR. NDA .NE. 0)) THEN
         CALL JOBEND('NFR/NDA should be zero when NMODE >= 4')
      END IF
      
*      TITMEM: TITLE OUTPUT FOR *.mem
*     MEED(1) = 0: THE CONTRIBUTION OF ANOMOLOUS DISPERSION IS CONSIDERED
*                  IN THE CALCULATION OF E.S.D.'S OF INTEGRATED INTENSITIES.
*             = 1: DO NOTHING.
*     MEED(2): SPACE GROUP NUMBER (NOT USED HERE)
*     MEED(3): SETTING NUMBER (NOT USED HERE)
*     NFPX(1), MEED(4): PIXEL NUMBER ALONG THE X AXIS
*     NFPX(2), MEED(5): PIXEL NUMBER ALONG THE Y AXIS
*     NFPX(3), MEED(6): PIXEL NUMBER ALONG THE Z AXIS
*     MEED(7) = 0: DECOMPOSE ALL THE RELECTIONS BY RIETAN
*             = 1: PART OF OVERLAPPED REFLECTIONS ARE GROUPED
*     MEED(8) = 0: 'OBSERVED' Fo'S ARE OUTPUT
*             = 1: Fc'S IN RIETVELD ANALYSIS ARE OUTPUT
*        EPSD: REFLECTIONS WHOSE DIFFERENCE IN THE D SPACING ARE SMALLER THAN EPSD
*              WILL BE GROUPED.
*              EFFECTIVE WHEN ONLY MEED(7) = 1
*        SCIO: SCALE FACTOR FOR INTEGRATED INTENSITIES TO ADJUST 
*              SIGMA(Io).
*       TSCAT: TOTAL ELECTRON NUMBERS (X-RAY DIFFRACTION) OR TOTAL
*              COHERENT SCATTERING LENGTHS (NEUTRON DIFFRACTION).
*      TSCAT1: TOTAL ELECTRON NUMBERS (X-RAY DIFFRACTION) OR TOTAL
*              POSITIVE COHERENT SCATTERING LENGTHS (NEUTRON DIFFRACTION).
*      TSCAT2: ZERO (X-RAY DIFFRACTION) OR TOTAL NEGATIVE COHERENT SCATTERING
*              LENGTHS (NEUTRON DIFFRACTION).
*       UCLAG: LAGRANGE'S UNDETERMINED COEFFICIENT.

      IF (NFR .GT. 0) THEN
         READ(4,*) (NFPX(J),J=1,3), TSCAT
         WRITE(6,635) ('NFPX(',J,') =',NFPX(J),J=1,3),'TSCAT =',TSCAT
      END IF
  635 FORMAT(//11X,'Parameters related to the file for Fourier',
     &  '/D synthesis'//11X,3(A,I1,A,I4,5X),A,1P,G13.6)

      IF (NMEM .GT. 0) THEN
         READ(4,*) TITMEM, MEED(1),(MEED(J),J=4,8), EPSD, SCIO, TSCAT1, 
     &   TSCAT2, UCLAG
         WRITE(6,640) TITMEM,'MEED(1) =',MEED(1),
     &   ('MEED(',J,') =',MEED(J),J=4,8),'EPSD =',EPSD,'SCIO =',
     &    SCIO,'TSCAT1 =',TSCAT1,'TSCAT2 =',TSCAT2,'UCLAG =',UCLAG
      END IF
  640 FORMAT(//11X,'Parameters related to the file for MEM analysis'
     &  //11X,'TITMEM : ',A,//11X,A,I4,5(5X,A,I1,A,I4),
     &  //9X,5(2X,A,1P,G13.6))

      CALL RFACTR(A,DEG,YFIT,XINT,NPTS,1)
C     PRINT OUT A TABLE OF INDICES, RELATIVE INTENSITIES, PEAK
C     POSITIONS, AND SO ON
      IF (NMODE .LT. 4) CALL COMPOS(A,SIGMAA,1,ANAME,INDIV,AW)

      IF (NFR .GE. 1 .OR. NMEM .GE. 1)
     &CALL FOURMEM(NSPGR,NSET,VNS,NEXCLD,NWEAK,NHIEND,NGRP)
      IF (NPRINT .GE. 1 .OR. NMEM .GT. 0) THEN
         WRITE(6,650)
         CALL RFLCTN(1)
         IF (NFR .GE. 1 .OR. NMEM .GE. 1)
     &   WRITE(6,660) NEXCLD, NWEAK, NHIEND, NGRP
      END IF
  650 FORMAT('1'/7X,'*** Summary of possible reflections (based on the r
     &efined parameters) ***'/)

  660 FORMAT(/11X,'NEXCLD =',I4,5X,'NWEAK =',I4,5X,'NHIEND =',I4,
     &  5X,'NGRP =',I4)
     
      IF (NPRINT .EQ. 2) THEN
*        THE PRINTER OUTPUT OF THE RIETVELD-REFINEMENT PATTERN IS NO
*        LONGER SUPPORTED
         ISCALE=0
         IDIF=0
         CALL LIST2(G,ID,DEG,XINT,YFIT,NTERMS,NRFN,NPTS,ISCALE,IDIF,
     &   NPRINT)
      END IF

C     RECORD RESULTS OF RIETVELD REFINEMENT
      IF (NMODE .NE. 1 .AND. (NPAT .GE. 2 .AND. NPAT .LE. 5)) THEN
         CALL OUT20(DEG,XINT,YFIT,NPTS,STEP,NREFL,TITLE,A,NEXC)
      END IF
      IF (NDA .GT. 0) THEN
         CALL MKFLOR(A,G,PHNAME(NDA),SIGMAA,NTERMS,NESD)
      END IF

*Rev 1.2 2003.5.28 {
*     Integrated intensities, YPEAK, kept constant during least-squares
*     fitting are finally calculated from refined parameters.
      IF (NMODE .EQ. 4 .AND. NCONST .EQ. 1) 
     &  CALL INTEGINT(A,DEG,YFIT,XINT,NPTS)
* }
*     Output *.ffo
      IF (NMODE .EQ. 4 .OR. NMODE .EQ. 5) CALL MKSFF(A,DEG,NPTS)

      IF ((NAUTO .EQ. 0 .AND. (NUPDT .EQ. 1 .OR. NUPDT .EQ. 2) .AND.
     &  IDCYCL .EQ. 0) .OR. (NLESQ .EQ. 2 .AND. (NUPDT .EQ. 1 .OR.
     &  NUPDT .EQ. 2)))  THEN   
         CALL OVERWR(A,NTERMS,NSPLBL,LABEL,LPAR,NLINE)
      END IF
*Rev 1.07 2002.08.22 Izumi
*     Update a CrystalMaker text file or a VICS file
      CALL UPDATE_LPP(A,SIGMAA)
*Rev 1.0a 2000.10.28 Izumi
      CLOSE(UNIT=5)
*Rev 1.0e 2000.11.28 Izumi
      CLOSE(UNIT=9)
      IF (NPRINT .EQ. 2) THEN
      WRITE(6,'(///11X,A/11X,A/11x,A//61X,A)') '"We should note that sto
     &pping and finishing are not the same thing.'
     &,'The refined model must make physical and chemical sense,',
     &'or it is not finished ... and even then it might be wrong!"',
     &'R. A. Young, 1993'
      WRITE(6,'(//11X,A/11X,A/11X,A//64X,A)')
     &'"If the fit of the assumed model is not adequate, the precision a
     &nd',
     &'accuracy of the parameters cannot be validly assessed by statisti
     &cal',
     &'methods."','E. Prince, 1981'
      END IF
      CALL JOBEND(' ')
      END

************************************************************************

      SUBROUTINE TINK
*     NAMED AFTER 'TINKER BELL' IN 'PETER PAN AND WENDY'
*     (1) DELETE COMMENTS AND VACANT LINES IN AN INPUT FILE.
*     (2) DECODE EACH LINE, EXTRACT VARIABLE(S), AND WRITE THEM TO UNIT
*         IOUT
*     (3) DETERMINE WHETHER OR NOT LINES ARE READ IN WHICH ARE SANDWICH-
*         ED BETWEEN 'If ..... then', 'else if ..... then', 'else', AND
*         'end if'.
      PARAMETER (MFLAG=200)
      CHARACTER LINE*80,FLAG(MFLAG)*10,LONGLINE*90
      INTEGER IVALUE(MFLAG)

*     FILE NUMBER FOR *.ins 
      INS = 5
*     SCRATCH FILE WHERE 'else' AND 'else if' HAVE BEEN CONVERTED
      INP = 2
*     SCRATCH FILE WHERE COMMNETS HAVE BEEN DELETED
      IOUT = 4

      REWIND(UNIT=INS)
      DO WHILE (.TRUE.)
*        CHECK WHETHER THE LENGTH OF THIS LINE EXCEEDS 80.
         READ(INS,'(A)',END=2) LONGLINE
         IF (LONGLINE(81:) .NE. ' ') CALL JOBEND
     &   ('Too long line: '//LONGLINE)
      END DO
    2 CALL ELSECONV(INS,INP)

      OPEN(UNIT=IOUT,STATUS='SCRATCH',ACCESS='SEQUENTIAL')
*     LABEL: =0, A LABEL DOES NOT APPEAR AS YET.
*            =1, LINES FOR LABELS, PARAMETERS, AND IDENTIFIERS ARE NOW
*                BEING INPUT.
*            =2  LINES FOR LABELS, PARAMETERS, AND IDENTIFIERS HAVE
*                ALREADY BEEN PASSED THROUGH.
      LABEL=0
      IA=0
*     NFLAG: NUMBER OF FLAGS
      NFLAG=0
*     LSKIP=1: LINES SANDWICHED BETWEEN THIS 'IF ..... THEN' AND
*              'END IF' PAIR ARE SKIPPED.
      LSKIP=0
*     IFENTR=1: THE PRESENT LINE IS SANDWICHED BETWEEN 'IF ..... THEN'
*               AND 'END IF' LINES.
      IFENTR=0
      DO WHILE (.TRUE.)
         READ(INP,'(A)',END=9) LINE
         CALL CHKLBL(LINE,LABEL,ID,IA)
         ICOM=INDEX(LINE,'#')
         IBEGIN=INDEX(LINE,'{')
         IEND=INDEX(LINE,'}')
         IF (LSKIP .EQ. 1 .AND. ID .EQ. 0) THEN
*           SKIP THIS LINE
         ELSE IF (INDEX(LINE,'Go to *') .GT. 0) THEN
            CALL SKIP(LINE(INDEX(LINE,'Go to *')+6:),INP)
            IF (IFENTR .EQ. 1) IFENTR = 0
         ELSE IF (LINE(1:1) .EQ. '*' .AND. LINE(2:2) .NE. ' ') THEN
*           A LABEL (DESTINATION OF 'Go to')
*           SKIP THIS LINE
         ELSE IF (ICOM .EQ. 1 .OR. LINE(:80) .EQ. ' ') THEN
*           SKIP THIS LINE
         ELSE IF (ICOM .NE. 0 .AND. LINE(:ICOM-1) .EQ. ' ') THEN
*           SKIP THIS LINE
         ELSE IF (IBEGIN .GT. 0 .AND. LINE(IBEGIN+1:) .EQ. ' ') THEN
*           (COMMENT) + '{'
*           SKIP THIS LINE
         ELSE IF (IEND .GT. 1 .AND. LINE(IEND+1:) .NE. ' ' .AND. 
     &   LINE(1:IEND-1) .EQ. ' ' .AND. LABEL .NE. 1) THEN
*           SPACE(S) + '}' + COMMENT (NOT PARAMETER LINES)
            WRITE(IOUT,'(A)') '}'
         ELSE IF (ICOM .NE. 0 .AND. IEND .EQ. 0 .AND. 
     &   IRSIN(LINE) .EQ. 0) THEN
            WRITE(IOUT,'(A)') LINE(:ICOM-1)
         ELSE IF (IEND .EQ. 1 .AND. LABEL .NE. 1) THEN
            WRITE(IOUT,'(A)') '}'
         ELSE IF (IEND .GT. 1 .AND. LINE(1:IEND-1) .NE. ' ' .AND.
     &   LABEL .NE. 1 .AND. IRSIN(LINE) .EQ. 0) THEN
*           INPUT DATA OTHER THAN PARAMETER LINES + '}'
            WRITE(IOUT,'(A/A)') LINE(:IEND-1),'}'
         ELSE IF (IEND .GT. 1 .AND. LINE(1:IEND-1) .NE. ' ' .AND.
     &   LABEL .NE. 1) THEN
*           Extract a variable and put '}'             
            CALL REGIST(LINE,IRSIN(LINE),NFLAG,FLAG,IVALUE,IOUT)
            WRITE(IOUT,'(A)') '}'
         ELSE IF (IEND .EQ. 1 .AND. LABEL .EQ. 1) THEN
            LABEL=2
            WRITE(IOUT,'(A)') 'ENDPARA'
         ELSE IF (IEND .GT. 1 .AND. LINE(:IEND-1) .EQ. ' ' .AND.
     &   LABEL .NE. 1) THEN
*           END OF DATA EXCEPT PARAMETERS TO CALCULATE INTENSITIES
            WRITE(IOUT,'(A)') '}'
         ELSE IF (IEND .GT. 1 .AND. LINE(:IEND-1) .EQ. ' ' .AND.
     &   LABEL .EQ. 1) THEN
            LABEL=2
            WRITE(IOUT,'(A)') 'ENDPARA'
         ELSE IF (IEND .GT. 1 .AND. ID .EQ. 0 .AND. LABEL .EQ. 1 .AND.
     &   (IEND .EQ. 80 .OR. LINE(IEND+1:) .EQ. ' ')) THEN
*           IDENTIFIERS + '}' IN THIS LINE (LACKING A COMMNET)
            CALL JOBEND
     &      ('"}" should be placed at the top of a separate line')
         ELSE IF (IEND .GT. 1 .AND. LINE(IEND+1:) .NE. ' ' .AND. 
     &   LABEL .EQ. 1 .AND. LINE(1:IEND-1) .NE. ' ') THEN
*           IDENTIFIERS + '}' + COMMENT
            CALL JOBEND
     &      ('"}" should be placed at the top of a separate line')
         ELSE IF (IFENTR .EQ. 1 .AND. (ID .EQ. 7 .OR. ID .EQ. 8 .OR.
     &   ID .EQ. 9)) THEN
*           INVALID 'IF ..... THEN' HAS BEEN DETECTED
            CALL JOBEND('''If ..... then'' nested')
         ELSE IF (IFENTR .EQ. 0 .AND. (ID .EQ. 7 .OR. ID .EQ. 8 .OR.
     &   ID .EQ. 9)) THEN
*           VALID 'IF ..... THEN' HAS BEEN DETECTED
            IFENTR=1
            CALL IFCHK(LINE,NFLAG,FLAG,IVALUE,LSKIP)
         ELSE IF (IFENTR .EQ. 0 .AND. ID .GE. 1 .AND. ID .LE. 6) THEN
*           INVALID 'END IF' HAS BEEN DETECTED
            CALL JOBEND('''end if'' without ''If ..... then''')
         ELSE IF (IFENTR .EQ. 1 .AND. ID .GE. 1 .AND. ID .LE. 6) THEN
*           VALID 'END IF' HAS BEEN DETECTED
            LSKIP=0
            IFENTR=0
         ELSE IF (IRSIN(LINE) .GT. 0 .AND. INDEX(LINE,'!') .GT. 0 .AND.
     &   IRSIN(LINE) .LT. INDEX(LINE,'!')) THEN
*           SKIP THIS LINE: (VARIABLE NAME) = (ITS VALUE)! COMMENT
         ELSE IF (IRSIN(LINE) .NE. 0) THEN
            CALL REGIST(LINE,IRSIN(LINE),NFLAG,FLAG,IVALUE,IOUT)
         ELSE 
            CALL PRLINE(LINE,IOUT)
         END IF
      END DO
    9 END
	  
************************************************************************

      SUBROUTINE ELSECONV(INS,IOUT)
*     CONVERT LINES RELATED TO 'else if' and 'else'
*     At present, only 'If ... then', 'else', and 'else if'
*     are allowed; If blocks such as 'IF ... THEN', 'Else', 'ELSE',
*     and  'Else if' cannot be accepted. 
*     IFELSE = 1: 'If' has been passed through.
*     IFELSE = 2: 'else' has been passed through.
*     IFELSE = 3: 'else if' has been passed through.
      INTEGER IFELSE
      CHARACTER LINE*80,STARNUM*3
      SAVE LINE, STARNUM

      REWIND(UNIT=INS)
*     Scratch file where 'else' and 'else if' are converted.
      OPEN(UNIT=IOUT,STATUS='SCRATCH',ACCESS='SEQUENTIAL')
      N = 1
      IFELSE = 0
      DO I = 1, 100000
         READ(INS,'(A)',END=8) LINE
         IF (LINE .EQ. ' ') THEN
            WRITE(IOUT,'(A)') ' '
            CYCLE
         END IF

         DO IS = 1, 80
*Rev 1.0m 2001.01.09 Izumi {
            IF (LINE(IS:IS) .NE. ' ') EXIT
*           IF (LINE(IS:IS) .NE. ' ') THEN
*              EXIT
*           ELSE IF (LINE(IS:IS) .EQ. '#') THEN
*              WRITE(IOUT,'(A)') LINE
*              CYCLE
*           END IF
* }
         END DO

         IF (LINE(IS:IS+2) .EQ. 'If ' .AND. INDEX(LINE(IS+3:),'then')
     &   .GT. 0) THEN
            WRITE(IOUT,'(A)') LINE
            IF (IFELSE .GT. 0) THEN
               WRITE(6,'(//11X,3A)') 'Nesting of if blocks has been ',
     &         'detected: ', LINE
               CALL JOBEND(' ')
            END IF
            IFELSE = 1
         ELSE IF (LINE(IS:IS+6) .EQ. 'else if') THEN
            CALL GLABEL(N,STARNUM)
            WRITE(IOUT,'(2A/A/2A)') 'Go to ',STARNUM,'end if','I',
     &      LINE(IS+6:)
            IF (IFELSE .EQ. 0) THEN
*Rev 1.0n 2001.01.10 Izumi
               WRITE(6,'(//11X,2A,I4)') '''else if'' without a ',
     &         'corresponding ''If'' or ''else if'' at Line #',I
               CALL JOBEND(' ')
            END IF
            IFELSE = 3
         ELSE IF (LINE(IS:IS+4) .EQ. 'else ') THEN
            CALL GLABEL(N,STARNUM)
            WRITE(IOUT,'(2A/A)') 'Go to ',STARNUM,'end if'
            IF (IFELSE .EQ. 0) THEN
*Rev 1.0n 2001.01.10 Izumi {
               WRITE(6,'(//11X,2A,I4)') '''else'' without a ',
     &         'corresponding ''If'' or ''else if'' at Line #',I
* }
               CALL JOBEND(' ')
            END IF
            IFELSE = 2
         ELSE IF (LINE(IS:IS+5) .EQ. 'end if') THEN
            SELECT CASE (IFELSE)
               CASE (0)
*Rev 1.0n 2001.01.10 Izumi
                  WRITE(6,'(//11X,2A,I4)') '''end if'' without a ',
     &            'corresponding ''If'' or ''else'' at Line #',I
                  CALL JOBEND(' ')
               CASE (1)
                  WRITE(IOUT,'(A)') LINE
               CASE (2)
                  CALL GLABEL(N,STARNUM)
                  WRITE(IOUT,'(A)') STARNUM
                  N = N + 1
               CASE (3)
                  CALL GLABEL(N,STARNUM)
                  WRITE(IOUT,'(A/A)') LINE,STARNUM
                  N = N + 1
            END SELECT
            IFELSE = 0
         ELSE
            WRITE(IOUT,'(A)') LINE
         END IF
      END DO
    8 REWIND IOUT
      END
	  
************************************************************************

      SUBROUTINE GLABEL(N,STARNUM)
      CHARACTER*3 STARNUM
	  
      IF (N .LT. 10) THEN
         WRITE(STARNUM,'(A,I1)') '*',N
      ELSE IF (N .LT. 100) THEN
         WRITE(STARNUM,'(A,I2)') '*',N
      ELSE
         CALL JOBEND('Too many else''s have appeared')
      END IF
      END

************************************************************************

      SUBROUTINE CHKLBL(LINE,LABEL,ID,IA)
*     SET LABEL AT 1 WHEN A LABEL FOR PARAMETERS APPEARED FOR THE FIRST
*     TIME
      CHARACTER LINE*80

      CALL KWCHK(LINE,ID)
      IF (LABEL .NE. 0 .OR. ID .NE. 0 .OR. INDEX(LINE,'#') .NE. 0 .OR.
     &  LINE .EQ. ' ') RETURN
*     IA: =0, A PHASE LINE HAS NOT APPEARED AS YET.
*         =1, A PHASE LINE HAS ALREADY BEEN ENCOUNTERED.
      IF (INDEX(LINE,'''A-') .GT. 0 .OR. INDEX(LINE,'''I-') .GT. 0) IA=1
      IF (IA .EQ. 0) RETURN

      IF (IRSIN(LINE) .NE. 0) RETURN
      DO I = 1, 80
         IF (LINE(I:I) .NE. ' ') EXIT
      END DO
      IF (INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ',LINE(I:I)) .EQ. 0) RETURN
      DO J = I + 1, 80
         IF (LINE(J:J) .EQ. ' ' .AND. J .NE. 4) THEN
            LABEL=1
            RETURN
         ELSE IF (LINE(J:J) .EQ. ' ' .AND. J .EQ. 4 .AND. LINE(1:3) .NE.
     &   'HKL' .AND. LINE(1:3) .NE. 'H0L' .AND. LINE(1:3) .NE. 'HK0'
     &   .AND. LINE(1:3) .NE. 'HHL') THEN
            LABEL=1
            RETURN
         ELSE IF (INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789',
     &   LINE(J:J)) .EQ. 0) THEN
*           A NON-ALPHANUMERIC CHARACTER HAS APPEARED
            RETURN
         END IF
      END DO
      CALL JOBEND('A too long label appeared')
      END

************************************************************************

      SUBROUTINE SKIP(STLBL,INP)
*     SKIP THE PART FROM 'Go to *Destination' TO '*Destination'
      CHARACTER STLBL*(*),LINE*80

      DO WHILE (.TRUE.)
         READ(INP,'(A)',END=2) LINE
         IF (LINE .EQ. STLBL) RETURN
      END DO
    2 WRITE(6,100) STLBL
  100 FORMAT(//11X,'Label not found: ',A)
      CALL JOBEND(' ')
      END

************************************************************************

      SUBROUTINE KWCHK(LINE,ID)
*     CHECK WHETHER OR NOT A KEYWORD IS PRESENT IN THE PRESENT LINE
      PARAMETER (NKW=9)
      CHARACTER LINE*80,KEYWRD(NKW)*10
      INTEGER LENGTH(NKW)
      DATA (KEYWRD(I),I=1,NKW) /'END IF ','End if ','end if ','ENDIF ',
     &  'Endif ','endif ','IF ','If ','if '/
      DATA (LENGTH(I),I=1,NKW) /7,7,7,6,6,6,3,3,3/

      ID=0
*Rev 1.0o 2001.01.29 Izumi {
C     ICOM = INDEX(LINE,'#')
C     IF (ICOM .EQ. 1 .OR. (ICOM .GT. 1 .AND. LINE(:ICOM-1) .EQ. ' '))
C    &  RETURN
      DO J = 1, 80
	   IF (LINE(J:J) .NE. ' ') EXIT
	END DO

      DO I = 1, NKW
C        IPOS=INDEX(LINE,' '//KEYWRD(I)(:LENGTH(I)-1))
C        IF (LINE(:LENGTH(I)) .EQ. KEYWRD(I)(:LENGTH(I)) .OR.
C    &   INDEX(LINE,' '//KEYWRD(I)) .NE. 0 .OR. IPOS+LENGTH(I) .EQ. 81)
C    &   THEN
         IF (LINE(J:J+LENGTH(I)-1) .EQ. KEYWRD(I)(:LENGTH(I))) THEN
            ID=I
            IF (ID .GE. 7 .AND. INDEX(LINE,' THEN') .EQ. 0 .AND.
     &      INDEX(LINE,' then') .EQ. 0)
     &      CALL JOBEND('''then'' is missing after ''If''')
            EXIT
         END IF
* }
      END DO
      END

************************************************************************

      INTEGER FUNCTION IRSIN(LINE)
*     CHECK WHETHER AN INTEGER (I), A REAL (R), OR A STRING (S) IS GIVEN
*     A VALUE IN THIS LINE
      CHARACTER LINE*80

      IRSIN=0
*     LINE(I1:I2): FIRST WORD
      DO I = 1, 80
         IF (LINE(I:I) .NE. ' ') EXIT
      END DO
      I1=I
      IF (INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ',LINE(I1:I1)) .EQ. 0) RETURN

      DO I = I1, 80
         IF (LINE(I:I) .EQ. ' ' .OR. LINE(I:I) .EQ. '=') GO TO 2
         IF (INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789@',LINE(I:I)) 
     &   .EQ. 0) RETURN
      END DO
      RETURN
    2 I2=I

      IF (I1 .EQ. 1 .AND. I2-1 .EQ. 1 .AND. 
     &   INDEX('HKL',LINE(1:1)) .NE. 0) THEN
*        A LINE SUCH AS 'H=2N+1' IS INPUT (A VALUE IS NOT GIVEN)
         RETURN
      END IF
         
      DO I = I2, 80
         IF (LINE(I:I) .NE. ' ') GO TO 3
      END DO
      RETURN
    3 I3=I
      IF (LINE(I3:I3) .NE. '=') RETURN

      DO I = I3+1, 80
         IF (LINE(I:I) .NE. ' ') GO TO 4
      END DO
      RETURN
    4 I4=I
*     CHECK THE FIRST CHARACTER AFTER '='
      IF (INDEX('0123456789+-.''',LINE(I4:I4)) .EQ. 0) RETURN
      
*     THE FIRST CHARACTER OF A VARIABLE NAME
      IRSIN=I1
      END      

************************************************************************

      SUBROUTINE REGIST(LINE,N1,NFLAG,FLAG,IVALUE,IOUT)
*     GET THE FLAG NAME OF AN INTEGER VARIABLE AND ITS VALUE, THEN WRITE
*     ONLY THE VALUE OF THE INTEGER VARIABLE OR A REAL VARIABLE TO UNIT
*     IOUT (=4)
      PARAMETER (MFLAG=200)
      CHARACTER LINE*80,FLAG(*)*10
      INTEGER IVALUE(*)
      CHARACTER*4 FMTINT

      DO I = N1+1, 80
         IF (LINE(I:I) .EQ. ' ' .OR. LINE(I:I) .EQ. '=') EXIT
      END DO
      N2=I-1
      IF (N2-N1+1 .GT. 10) THEN
         CALL JOBEND('Too long variable name')
      END IF

*     LINE(IV1:IV2): INTEGER (e.g. 2), REAL (e.g. 3.14159), OR STRING
*                    (e.g. 'CU')
      DO I = N2, 80
         IF (LINE(I:I) .EQ. '=') GO TO 3
      END DO
      CALL JOBEND('''='' is missing after a variable name')
    3 IEQ=I
      DO I = IEQ+1, 80
         IF (LINE(I:I) .NE. ' ') GO TO 4
      END DO
      CALL JOBEND('A value is missing after ''=''')
    4 IV1=I

      IF (LINE(IV1:IV1) .EQ. '''') THEN
*        CHARACTER VARIABLE
         DO I = IV1+1, 80
            IF (LINE(I:I) .EQ. '''') GO TO 6
         END DO
         CALL JOBEND('A string should be enclosed by a pair of '//
     &   'apostrophes')
    6    IV2=I
         WRITE(IOUT,'(A)') LINE(IV1:IV2)    
         RETURN
      END IF

*     REAL OR INTEGER
      DO I = IV1, 80
         IF (LINE(I:I) .EQ. ' ' .OR. LINE(I:I) .EQ. ':') EXIT
         IF (INDEX('0123456789EeDd+-.',LINE(I:I)) .EQ. 0) THEN
*Rev 1.01 2001.06.12 Izumi {
            CALL JOBEND('Invalid real or integer value, '//LINE(I:I)//
     &      ', in: '//LINE)
* }
         END IF
      END DO
      IV2=I-1
*     An integer variable whose tail is '@' is not written
*Rev 1.0k 2001.01.02 Izumi {
C     IF (LINE(IV2:IV2) .NE. '@') WRITE(IOUT,'(A)') LINE(IV1:IV2)
      IF (LINE(N2:N2) .NE. '@') WRITE(IOUT,'(A)') LINE(IV1:IV2)
* }

*     REAL VARIABLE (DOES NOT REGISTER)
      IF (INDEX('IJKLMN',LINE(N1:N1)) .EQ. 0) RETURN

      IF (INDEX(LINE(IV1:IV2),'.') .GT. 0  .OR. INDEX(LINE(IV1:IV2),'E')
     &  .GT. 0 .OR. INDEX(LINE(IV1:IV2),'e') .GT. 0 .OR. 
     &  INDEX(LINE(IV1:IV2),'.') .GT. 0) THEN
           CALL JOBEND('An real value is input for an integer '//
     &     'variable: '//LINE(IV1:IV2))
      END IF      
*     INTEGER VARIABLE
      NFLAG=NFLAG+1
      IF (NFLAG .GT. MFLAG) CALL JOBEND('NFLAG exceeded MFLAG')
      FLAG(NFLAG)=LINE(N1:N2)
      DO I = 1, NFLAG-1
         IF (FLAG(NFLAG) .EQ. FLAG(I)) CALL JOBEND
     &   ('The same name of an integer variable has appeared again: '//
     &   FLAG(I))
      END DO
      IF (IV2-IV1+1 .GT. 9) THEN
         CALL JOBEND('Integer value out of range')
      ELSE
         WRITE(FMTINT,'(A,I1,A)') '(I',IV2-IV1+1,')'
         READ(LINE(IV1:IV2),FMTINT) IVALUE(NFLAG)	 
      END IF
      END      

************************************************************************

      SUBROUTINE IFCHK(LINE,NFLAG,FLAG,IVALUE,LSKIP)
*     DECODE AN 'IF ..... THEN' LINE
      CHARACTER LINE*80,FLAG(*)*10
      INTEGER IVALUE(*)

      N1=INDEX(LINE,'IF ')
      IF (N1 .EQ. 0) N1=INDEX(LINE,'If ')
      IF (N1 .EQ. 0) N1=INDEX(LINE,'if ')
*     CHECK THE FIRST LOGICAL EXPRESSION
      CALL CHKVAL(LINE,NFLAG,FLAG,IVALUE,N1+3,LSKIP,N6)
      
*     LINE(NA:NB): 'THEN', 'OR', OR 'AND'
      DO I = N6+1, 80
         IF (LINE(I:I) .NE. ' ') GO TO 1
      END DO
      CALL JOBEND('''then'', ''or'', or ''and'' is missing')
    1 NA=I
      DO I = NA+1, 80
         IF (LINE(I:I) .EQ. ' ') EXIT
      END DO
      NB=I-1
      IF (LINE(NA:NB) .NE. 'THEN' .AND. LINE(NA:NB) .NE. 'then'
     &  .AND. LINE(NA:NB) .NE. 'OR' .AND. LINE(NA:NB) .NE. 'or' .AND.
     &  LINE(NA:NB) .NE. 'AND' .AND. LINE(NA:NB) .NE. 'and') THEN
         CALL JOBEND('''then'', ''or'', or ''and'' is missing')
      END IF
     
*     LSKIP=0: THE LOGICAL EXPRESSION IS TRUE (READ LINES BETWEEN
*              'IF ..... THEN' AND 'END IF')
*     LSKIP=1: THE LOGICAL EXPRESSION IS FALSE (SKIP LINES BETWEEN
*              'IF ..... THEN' AND 'END IF')
      IF (NB-NA+1 .EQ. 4 .AND. (LINE(NA:NB) .EQ. 'THEN' .OR. 
     &  LINE(NA:NB) .EQ. 'then')) THEN
*        THE LOGICAL EXPRESION HAS ALREADY BEEN CHECKED
         RETURN
      ELSE IF (NB-NA+1 .EQ. 2 .AND. (LINE(NA:NB) .EQ. 'OR' .OR.
     &  LINE(NA:NB) .EQ. 'or') .AND. LSKIP .EQ. 0) THEN
          RETURN
      ELSE IF (NB-NA+1 .EQ. 3 .AND. (LINE(NA:NB) .EQ. 'AND' .OR.
     &  LINE(NA:NB) .EQ. 'and') .AND. LSKIP .EQ. 1) THEN
          RETURN
      END IF

*     CHECK THE SECOND LOGICAL EXPRESSION
      CALL CHKVAL(LINE,NFLAG,FLAG,IVALUE,NB+1,LSKIP,N6)
      END

************************************************************************

      SUBROUTINE CHKVAL(LINE,NFLAG,FLAG,IVALUE,NSTART,LSKIP,N6)
*     CHECK EACH LOGICAL EXPRESSION
      CHARACTER LINE*80,FLAG(*)*10,BUFF*120
      INTEGER IVALUE(*),LSKIP

*     LINE(N1:N2): VARIABLE NAME
      DO I = NSTART, 80
         IF (LINE(I:I) .NE. ' ') EXIT
      END DO
      N1=I
      DO I = N1+1, 80
         IF (INDEX(' >=<',LINE(I:I)) .NE. 0) EXIT
      END DO
      N2=I-1
      DO I = 1, NFLAG
         IF (FLAG(I) .EQ. LINE(N1:N2)) GO TO 3
      END DO
      BUFF = 'Invalid flagname: '//LINE(N1:N2)
      CALL JOBEND(BUFF)
*     NN: FLAG NUMBER FOR THIS VARIABLE
    3 NN=I

*     LINE(N3:N4): '>', '=', '<', '>=', '<=', OR '<>' 
      DO I = N2+1, 80
         IF (LINE(I:I) .NE. ' ') EXIT
      END DO
      N3=I
      DO I = N3+1, 80
         IF (INDEX('>=<',LINE(I:I)) .EQ. 0) EXIT
      END DO
      N4=I-1

*     LINE(N5:N6): INTEGER VALUE (= IFIG)
      DO I = N4+1, 80
         IF (LINE(I:I) .NE. ' ') EXIT
      END DO
      N5=I
      DO I = N5+1, 80
         IF (LINE(I:I) .EQ. ' ') EXIT
      END DO
      N6=I-1
      IF (N6-N5+1 .EQ. 1) THEN
         READ(LINE(N5:N6),'(I1)') IFIG
      ELSE IF (N6-N5+1 .EQ. 2) THEN
         READ(LINE(N5:N6),'(I2)') IFIG
      ELSE
         BUFF = 'Integer value out of range: '//LINE(N5:N6)
         CALL JOBEND(BUFF)
      END IF

*     LSKIP=0: THE LOGICAL EXPRESSION IS TRUE (READ LINES BETWEEN 'IF
*              ..... THEN' AND 'END IF')
*     LSKIP=1: THE LOGICAL EXPRESSION IS FALSE (SKIP LINES BETWEEN 'IF
*              ..... THEN' AND 'END IF')
      LSKIP=1
      IF (N4-N3+1 .EQ. 1 .AND. LINE(N3:N4) .EQ. '>' .AND. 
     &  IVALUE(NN) .GT. IFIG) THEN
         LSKIP=0
      ELSE IF (N4-N3+1 .EQ. 1 .AND. LINE(N3:N4) .EQ. '=' .AND. 
     &  IVALUE(NN) .EQ. IFIG) THEN
         LSKIP=0
      ELSE IF (N4-N3+1 .EQ. 1 .AND. LINE(N3:N4) .EQ. '<' .AND. 
     &  IVALUE(NN) .LT. IFIG) THEN
         LSKIP=0
      ELSE IF (N4-N3+1 .EQ. 2 .AND. LINE(N3:N4) .EQ. '>=' .AND. 
     &  IVALUE(NN) .GE. IFIG) THEN
         LSKIP=0
      ELSE IF (N4-N3+1 .EQ. 2 .AND. LINE(N3:N4) .EQ. '<=' .AND. 
     &  IVALUE(NN) .LE. IFIG) THEN
         LSKIP=0
      ELSE IF (N4-N3+1 .EQ. 2 .AND. LINE(N3:N4) .EQ. '<>' .AND. 
     &  IVALUE(NN) .NE. IFIG) THEN
         LSKIP=0
      END IF
      END

************************************************************************

      SUBROUTINE CHKMIX(MIX)
*     MIX=0: NMIX = 0; MIX=1: NMIX > 0
      CHARACTER*80 LINE

      MIX=1
      READ(4,'(A)') LINE
      ISQ=0
      DO I = 1, 80
         IF (LINE(I:I) .EQ. '''') ISQ=ISQ+1
         IF (ISQ .EQ. 3) GO TO 1
      END DO

      IF (ISQ .EQ. 2) THEN
C        ONLY ONE STRING IS INPUT IN THIS LINE, WHICH IS THEREFORE 
C        REGARDED AS A PHASE LINE 
         MIX=0
         BACKSPACE 4
         RETURN
      ELSE
         CALL JOBEND('Bad imaginary species line or phase line')
      END IF
    1 IF (LINE(I+1:I+2) .EQ. 'I-' .OR. LINE(I+1:I+2) .EQ. 'A-') THEN
*        THIS IS NOT AN IMAGINARY SPECIES LINE BUT A PHASE LINE
         MIX=0
      END IF
      BACKSPACE 4
      END
      
************************************************************************

      SUBROUTINE RDFRCT(LINE,ANAME,LL,IMIX,NREAL,AW)
*     READ THE CONSTITUENTS AND FRACTIONS FOR A SITE
      PARAMETER (NA=15)
      REAL AW(*)
      CHARACTER*5 ANAME(*),NAMEA(5)
      CHARACTER*80 LINE
      COMMON /MAGSC/ MAGATM(NA),LCMFF(NA),CMFF(7,NA)
      COMMON /MIX/ FRCTN(5,10),MFRAC(5,10),NOFRAC(10),NMIX
      COMMON /U/ DELTF1(NA),DELTF2(NA)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD

      IF (LAPO(LINE) .EQ. 0) THEN
         CALL JOBEND('The name of virtual atom'
     &   //' should be enclosed by apostrophes')
      END IF
      DO I = 80, 1, -1
         IF (LINE(I:I) .EQ. '/') THEN
            EXIT
         ELSE IF (LINE(I:I) .NE. ' ') THEN
            CALL JOBEND('The last character in a virtual species line '
     &      //'should be ''/''')
         END IF
      END DO
      CALL BACKSP(LINE)

      DO J = 1, 5
         NAMEA(J)=' '
      END DO

C*    CONSTITUENTS OF SITES
C     NOFRAC: NUMBER OF ATOMIC SPECIES OCCUPYING THE SITE
C     MFRAC:  NUMBER WHICH REPRESENTS THE KIND OF THE ATOM
C     FRCTN:  FRACTION OF THE MFRAC'TH ATOM
      READ(4,*) ANAME(LL),(NAMEA(J),FRCTN(J,IMIX),J=1,5)
      DO J = 5, 1, -1
         IF (NAMEA(J) .NE. ' ') THEN
            NOFRAC(IMIX)=J
            EXIT
         END IF
      END DO
      DO J = 1, NOFRAC(IMIX)
         MFRAC(J,IMIX)=NUMATM(NAMEA(J),ANAME,NREAL)
      END DO

      MAGATM(LL)=0
      DELTF2(LL)=0.0
      AW(LL)=0.0
      DO J = 1, NOFRAC(IMIX)
         IF (NBEAM .GE.1) THEN
            DELTF2(LL)=DELTF2(LL)+DELTF2(MFRAC(J,IMIX))*FRCTN(J,IMIX)
         ELSE IF (MAGATM(MFRAC(J,IMIX)) .EQ. 1) THEN
C           ONE OF THE CONSTITUENTS OF THE SITE IS OF MAGNETIC MOMENT
            MAGATM(LL)=1
         END IF
         AW(LL)=AW(LL)+AW(MFRAC(J,IMIX))*FRCTN(J,IMIX)
      END DO
      END

************************************************************************

      SUBROUTINE ASFDC(NREAL,IDREAL,ANAME,ANAME2,AW,NBEAM,AA,BB,C,DC,
     &  DELTF1,DELTF2,NTARG,SIGA,SIGI,CATOM,ATOMNO,DENSTY)
*     SEARCH THE ASFDC DATABASE
C     AA(N,KATOM), BB(N,KATOM), AND C(KATOM): COEFFICIENTS FOR ANALYTIC
C     APPROXIMATIONS OF SCATTERING FACTORS
C     DELTF1: REAL DISPERSION TERM
C     DELTF2: IMAGINARY DISPERSION TERM
C     'INTERNATIONAL TABLES FOR X-RAY CRYSTALLOGRAPHY,' VOL. IV (1974).
      PARAMETER (NA=15)
      INTEGER IDREAL(*)
      REAL AW(*),AA(4,NA),BB(4,NA),C(*),DC(2,6),DELTF1(*),DELTF2(*),
     &  SIGA(*),SIGI(*),CATOM(*),ATOMNO(*),MW
      CHARACTER*5 DNAME,ANAME(*),ANAME2(*)

      IREAL=1
      CALL OPENASFDC
      REWIND 2
* The total number of chemical species is greater than 211
*     DO 50 J = 1, 211
      DO 50 J = 1, 300
         READ(2,'(A5)',END=2) DNAME
         DO JJ = 1, NREAL
            IF (IDREAL(JJ) .EQ. 0 .AND. DNAME .EQ. ANAME2(JJ)) THEN
               BACKSPACE(UNIT=2)
               READ(2,'(9X,F10.0)') AW(JJ)
               IF (NBEAM .GE. 1) THEN
                  READ(2,*) (AA(N,JJ),BB(N,JJ),N=1,4),C(JJ)
                  READ(2,*) ((DC(MM,NN),MM=1,2),NN=1,6)
                  IF (NBEAM .EQ. 1) THEN
                     DELTF1(JJ)=DC(1,NTARG)
                     DELTF2(JJ)=DC(2,NTARG)
                  END IF
               ELSE
                  READ(2,*) (AA(N,JJ),BB(N,JJ),N=1,4),DUM1,C(JJ)
                  READ(2,*) (SIGA(JJ),II=1,12),SIGA(JJ),SIGI(JJ)
C                 BARN ==> M**2
                  SIGA(JJ)=SIGA(JJ)*1.0E-28
                  SIGI(JJ)=SIGI(JJ)*1.0E-28
                  DELTF1(JJ) = 0.0
		  DELTF2(JJ) = 0.0
               END IF
               IDREAL(JJ)=1
               IREAL=IREAL+1
               IF(IREAL.GT.NREAL) GO TO 3
               GO TO 50
            ENDIF
         END DO
C        SKIP TWO LINES
         READ(2,'()')
         READ(2,'()')
   50 CONTINUE
    2 WRITE(6,'(//11X,A)') 'Wrong chemical symbol(s)'
      DO J = 1, NREAL
         IF(IDREAL(J).EQ.0) WRITE(6,'(13X,A5)') ANAME(J)
      END DO
      CALL JOBEND(' ')

    3 IF (NBEAM .EQ. 0) THEN
C        MW: MOLECULAR WEIGHT
         MW=0.0
         DO J = 1, NREAL
            MW=MW+CATOM(J)*AW(J)
         END DO
C        ATOMNO: NUMBER OF ATOMS PER UNIT VOLUME/M**(-3)
         DO J = 1, NREAL
            ATOMNO(J)=DENSTY/MW*6.0221367E23*CATOM(J)*1.0E6
         END DO
         WRITE(6,100) (J,ANAME(J),C(J),SIGA(J),SIGI(J),CATOM(J),AW(J),
     &   ATOMNO(J),J=1,NREAL)
         RETURN
      END IF
  100 FORMAT(//11X,
     &'Scattering lengths, cross sections, and amounts of atoms'//13X,
     &'No.  Atom    b/fm      SIGMAA/m**2   SIGMAI/m**2   CATOM/mol  At.
     &wt.       ATOMNO/m**(-3)'/
     &(13X,I2,3X,A5,0P,F8.3,2(5X,1PE9.3),4X,0P,F7.3,F13.4,5X,1PE10.4))

      WRITE(6,150) (JJ,ANAME(JJ),(AA(N,JJ),
     &BB(N,JJ),N=1,4),C(JJ),DELTF1(JJ),DELTF2(JJ),JJ=1,NREAL)
  150 FORMAT(//' ',10X,'Coefficients for analytic approximations of scat
     &tering factors and dispersion terms'//10X,'No. Atom',4X,
     &'a1',9X,'b1',9X,'a2',9X,'b2',9X,'a3',9X,'b3',9X,'a4',9X,'b4',9X,
     &'c',7X,'DELTF1  DELTF2'/(10X,I2,1X,A5,F10.6,8F11.6,2F8.3))
      END

************************************************************************

      SUBROUTINE FLSITE(A,L,NOAT,INDIV,LISO,MAGATM,NSITE,LMAG,NATOM,
     &  ISONUM,NLINE,LPAR,LABEL,IP,ANAME)
*     DETERMINE LISO, LMAG, AND ISONUM
      PARAMETER (NB=7000,NAP=150,NPH=8)
      REAL A(*)
      INTEGER NOAT(NAP,NPH),INDIV(*),LISO(NAP,NPH),LPAR(*),
     &  MAGATM(*),NSITE(*),LMAG(*),ISONUM(NAP,NPH)
      CHARACTER LABEL(*)*25,FORM(4)*4,ANAME(*)*5
      COMMON /H/ NQ(NPH),NCENTR(NPH),LAUEG(NPH)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /UNCONS/ NUCREF,IPHASE(90),IHKL(3,90),NOPPP(90),LPPP(NB),
     &  RCUT(2,90),NCUT
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      P(J) = A(KPHB(LL)+J)
      DATA FORM/'(I1)', '(I2)', '(I3)', '(I4)'/

*     DETERMINE THE VALUES OF NSITE AND NOAT
C     NSITE: NUMBER OF ATOMS IN ASYMMETRIC UNIT
C     NOAT: TYPES OF ATOMS IN ONE ASYMMETRIC UNIT, REPRESENTED BY
C           CODE NUMBERS

      ISITE=0
      IHIT=0
      IF (L .EQ. 1 .AND. (NMODE .EQ. 4 .OR. NMODE .EQ. 6)) THEN
*        INTEGER I BECOMES INDEFINITE WHEN THE NEXT LOOP IS SKIPPED.
         I = IP
         GO TO 1
      END IF

      DO I = IP, NLINE
         ISLA=INDEX(LABEL(I),'/')
         IF (ISLA .EQ. 0 .AND. IHIT .EQ. 1) THEN
            EXIT
         ELSE IF (ISLA .GT. 0) THEN
            ISITE=ISITE+1
            IF (ISITE .GT. NAP) CALL JOBEND('NAP is too small')
            IHIT=1
            LENGTH=LENSTR(LABEL(I))
            IF (INDEX('123456789',LABEL(I)(ISLA+1:ISLA+1)) .GT. 0) THEN
               READ(LABEL(I)(ISLA+1:LENGTH),FORM(LENGTH-ISLA))
     &         NOAT(ISITE,L)
            ELSE
               NOAT(ISITE,L)=NUMATM(LABEL(I)(ISLA+1:LENGTH),ANAME,
     &         NATOM)
            END IF
C           LISO=0: OVERALL TEMPERATURE FACTOR
C           LISO=1: INDIVIDUAL ISOTROPIC TEMPERATURE FACTOR
C           LISO=6: INDIVIDUAL ANISOTROPIC TEMPERATURE FACTOR
            IF (LPAR(I) .LE. 6) THEN
               IF (INDIV(L) .EQ. 1) THEN
                  LISO(ISITE,L)=1
               ELSE
                  LISO(ISITE,L)=0
               END IF
            ELSE
               IF (INDIV(L) .EQ. 0) THEN
                  CALL JOBEND('Check the value(s) of INDIV')
               END IF
               LISO(ISITE,L)=6
            END IF
C           NPSITE: STRUCTURE-FACTOR PARAMETERS PER SITE
            NPSITE(ISITE,L)=4+LISO(ISITE,L)+MAGATM(NOAT(ISITE,L))
            IF (NPSITE(ISITE,L) .LT. LPAR(I)) THEN
               WRITE(6,'(//11X,A,I3,A,I1,A)')
     &         'The number of parameters of site #',ISITE,
     &         ' for phase #',L,' is too large'
            ELSE IF (NPSITE(ISITE,L) .GT. LPAR(I)) THEN
               WRITE(6,'(//11X,A,I3,A,I1,A)')
     &         'The number of parameters of site #',ISITE,
     &         ' for phase #',L,' is too small'
               CALL JOBEND(' ')
            END IF
         END IF
      END DO
    1 NSITE(L)=ISITE

      IF (NMODE .EQ. 4 .AND. L .EQ. 1 .AND. NSITE(1) .NE. 0)
     &CALL JOBEND('No structure parameters should be given for '//
     &'the first phase when NMODE = 4')
      IF (NMODE .EQ. 4 .AND. L .GT. 1 .AND. NSITE(L) .EQ. 0)
     &CALL JOBEND('Structure parameters should be given except for '//
     &'the first phase when NMODE = 4')
      IF ((NMODE .LE. 3 .OR. NMODE .EQ. 5) .AND. NSITE(L) .EQ. 0)
     &CALL JOBEND('Structure parameters should be given for all the '//
     &'phases when NMODE <= 3 or NMODE = 5')
      IF (NMODE .EQ. 6 .AND. NSITE(L) .GT. 0)
     &CALL JOBEND('No structure parameters should be given when '//
     &'NMODE = 6')
      IP=I+1

C     LMAG=0: NONMAGNETIC MATERIAL
C     LMAG=1: MAGNETIC MATERIAL
      LMAG(L)=0
      DO J=1,NSITE(L)
         IF (MAGATM(NOAT(J,L)) .EQ. 1) THEN
            LMAG(L)=1
            EXIT
         END IF
      END DO
C     ISONUM: ARRAY TO ASSIGN A SERIES OF NEMBERS FOR EACH CHEMICAL
C             SPECIES
      DO J=1,NSITE(L)
         IF(NOAT(J,L).LT.1 .OR. NOAT(J,L).GT.NATOM) THEN
            WRITE(6,'(//11X,A,I2)') 'Atomic number out of range: ',
     &      NOAT(J,L)
            CALL JOBEND(' ')
         ENDIF
         KK=1
         DO JJ=1,J-1
            IF (NOAT(J,L) .EQ. NOAT(JJ,L)) KK=KK+1
         END DO
         ISONUM(J,L)=KK
         IF(KK.EQ.1) THEN
            DO JJ=J+1,NSITE(L)
               IF(NOAT(J,L).EQ.NOAT(JJ,L)) CYCLE
            END DO
C           A SERIES NUMBER IS NOT ASSIGNED
            ISONUM(J,L)=0
         ENDIF
      END DO
      
*     NPRIM: NUMBER OF PRIMARY PROFILE PARAMETERS PER REFLECTION
      SELECT CASE (NPRFN)
         CASE (0)
            NPRIM = 0
         CASE (1,3)
            SELECT CASE (NMODE)
               CASE (3)
*                 |Fc| is added
                  NPRIM = 5
               CASE (6)
*                 const*|Fc| and 2-theta(peak) are added
                  NPRIM = 6
               CASE DEFAULT
                  NPRIM = 4
            END SELECT
         CASE (2)
            SELECT CASE (NMODE)
               CASE (3)
*                 |Fc| is added
                  NPRIM = 6
               CASE (6)
*                 const*|FC| and 2-theta(peak) are added
                  NPRIM = 7
               CASE DEFAULT
                  NPRIM = 5
            END SELECT
      END SELECT
      
C     KPHB(L): PARAMETER NUMBER FOR THE SCALE FACTOR OF THE L'TH PHASE
C     KPHE(L): PARAMETER NUMBER FOR THE LAST STRUCTURE-FACTOR
C              PARAMETER OF THE L'TH PHASE
      IF (L .EQ. 1) THEN
         KPHB(L) = NPPPX + NRELAX*NPRIM
      ELSE
         KPHB(L)=KPHE(L-1)+1
      END IF
      NPARL=0
      DO J = 1, NSITE(L)
         NPARL = NPARL + NPSITE(J,L)
      END DO
      KPHE(L) = KPHB(L) + NPARL + NOTX
      IF (LMAG(L) .EQ. 1) THEN
         SELECT CASE (LAUEG(L))
            CASE (2,4)
*              MONOCLINIC (A- AND C-AXIS UNIQUE)
               CALL JOBEND('Use the b-axis as the principal axis in '//
     &         'magnetic materials')
            CASE (1,3,5)
*              TRICLINIC, MONOCLINIC (B-AXIS UNIQUE), AND ORTHORHOMBIC
               KPHE(L)=KPHE(L)+3
            CASE (6:13)
*              TETRAGONAL, HEXAGONAL, AND RHOMBOHEDRAL
               KPHE(L)=KPHE(L)+1
         END SELECT
      END IF
      
*     DECODE LABELS AND READ IN VARIOUS VALUES FOR REFLECTIONS 
*     BROADENED ANISOTROPICALLY
      IPAR = 1
*     SERIAL NUMBER FOR REFLECTIONS BROADENED ANISOTROPICALLY
      IANIS = 1
      DO I = 1, NLINE
         LDOTH = INDEX(LABEL(I),'.')
         IF (LABEL(I)(1:3) .EQ. 'PPP' .AND. LABEL(I)(5:5) .EQ. '_'
     &   .AND. LDOTH .GE. 7) THEN
*           A LABEL FOR A REFLECTION BROADENED ANISOTROPICALLY 
*           FORMAT: PPPn_h.k.l. where n is the phase number, and
*           h, k, l are reflection indices.
            IF (IANIS .GT. 90) CALL JOBEND
     &      ('Too many reflections broadened anisotropically')

*           DETERMINE THE PHASE NUMBER FOR THIS REFLECTION
            READ(LABEL(I)(4:4),'(I1)') IPHASE(IANIS)
*           INDEX h
            READ(LABEL(I)(6:LDOTH-1),*) IHKL(1,IANIS)
            LDOTK = LDOTH + INDEX(LABEL(I)(LDOTH+1:),'.')
*           INDEX k
            READ(LABEL(I)(LDOTH+1:LDOTK-1),*) IHKL(2,IANIS)
            LDOTL = LDOTH + INDEX(LABEL(I)(LDOTH+1:),' ')
*           INDEX l
            READ(LABEL(I)(LDOTK+1:LDOTL-1),*) IHKL(3,IANIS)

*           THE SERIAL PARAMETER NUMBER OF THE FIRST PRIMARY PROFILE 
*           PARAMETER FOR THIS REFLECTION SPECIFIED BY THE USER
            NOPPP(IANIS) = IPAR
            IANIS = IANIS + 1
         END IF
         IPAR = IPAR + LPAR(I)
      END DO
      NUCREF = IANIS -1

*     THE VALUE OF NQ IS REFERRED TO IN SUBROUTINE SF
      IF (INDIV(L) .EQ. 1 .AND. NCENTR(L) .EQ. 0) THEN
         NQ(L)=1
      ELSE IF(INDIV(L) .EQ. 1 .AND. NCENTR(L) .EQ. 1) THEN
         NQ(L)=2
      ELSE IF(INDIV(L) .EQ. 0 .AND. NCENTR(L) .EQ. 0) THEN
         NQ(L)=3
      ELSE IF(INDIV(L) .EQ. 0 .AND. NCENTR(L) .EQ. 1) THEN
         NQ(L)=4
      END IF
*     THE FOLLOWING LINE IS REQUIRED BECAUSE SUBROUTINE FLSITE IS
*     INCLUDED IN A LOOP SCANNING L IN THE MAIN PROGRAM
      IF (L .NE. NPHASE) RETURN

*     CHECK PARAMETERS RELATED TO PROFILE ASYMMETRY
*     'L' is not used because it is an argument of this subroutine.
      DO LL = 1, NPHASE
         SELECT CASE (NPRFN)
            CASE (0)
               IF (NASYM .EQ. 0 .AND. (P(10) .EQ. 0.0 .OR. 
     &         P(11) .NE. 0.0)) THEN
                  CALL JOBEND('Check rd and the next parameter')
               ELSE IF (NASYM .EQ. 1 .AND. P(10) .NE. 0.0) THEN
                  CALL JOBEND('Check the parameter following As')
               END IF
            CASE (1:3)
               IF (P(11) .EQ. 0.0) CALL JOBEND
     &         ('Check decay parameters')
         END SELECT
      END DO
      END

************************************************************************

      SUBROUTINE CHGID(NTERMS,NCNSTR)
*     SET ID VALUES AT 0 FOR STRUCTURE PARAMETERS BECAUSE THEY ARE
*     NEITHER REFINED NOR CONSTRAINED WHEN NMODE >= 2
      PARAMETER (NAP=150,NA=15,NT=999,NPH=8)
      COMMON /R/ ID(NT),IR(NT),IDSAVE(NT),NRFN
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)

      DO L = 1, NPHASE
         DO J = KPHB(L)+NOTX, KPHE(L)
            IF (ID(J) .NE. 0) THEN
               ID(J) = 0
               IDSAVE(J) = 0
            END IF
         END DO
      END DO
      
*     REDETERMINE NCNSTR
      NCNSTR = 0
      DO J = 1,NTERMS
         IF (ID(J) .EQ. 2) NCNSTR = NCNSTR + 1
      END DO
      END

************************************************************************

      INTEGER FUNCTION NUMATM(ATOM,ANAME,N)
*     RETURN AN NOAT VALUE OF A SITE
      PARAMETER(NA=15,NB=7000)
      CHARACTER ATOM*(*),ANAME(*)*(*),TMP*5

*     N: NREAL OR NATOM
      DO I = 1, N
         IF (ATOM .EQ. ANAME(I)) THEN
            NUMATM=I
            RETURN
         END IF
      END DO
*     Dummy statement to avoid a warning message in Digital Fortran
      NUMATM = 0
      TMP=ATOM
      CALL JOBEND('A bad chemical species name was input: '//TMP)
      END

************************************************************************

      SUBROUTINE DETREF(IPAR,NPAR,NCHNG,ID)
*     DETERMINE PARAMETERS TO BE REFINED IN EACH CYCLE
      PARAMETER (NR=400,NAP=150,NPH=8)
*Rev 1.1 2003.05.16 Izumi
      INTEGER IPAR(NR,50),NPAR(*),ID(*)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      LP(I)=KPHB(L)+I
      LD(I)=ID(LP(I))

*     1ST CYCLE: LINEAR PARAMETERS
      IP=1
*     BACKGROUND PARAMETERS
      DO I = NBCKGR, NBCKGR + NUMBG - 1
         IF (ID(I) .EQ. 1) THEN
            IPAR(IP,1)=I
            IP=IP+1
         END IF
      END DO
*     SCALE FACTOR(S)
      DO L = 1, NPHASE
         IF (LD(NSFX) .EQ. 1) THEN
            IPAR(IP,1)=LP(NSFX)
            IP=IP+1
         END IF
      END DO
      NPAR(1)=IP-1

*     2ND CYCLE: SECONDARY PROFILE PARAMETERS #1
*     U, V, W, X, Y, a0, a1, and a2
      IF (NPAR(1) .EQ. 0) THEN
         ICYC=1
      ELSE
         ICYC=2
      END IF
      IP=1
      DO L = 1, NPHASE
         IF (NMODE .NE. 1 .AND. NMODE .NE. 6) THEN
            SELECT CASE (NPRFN)
               CASE (0)
*                 Pseudo-Voigt function of Thompson, Cox, and Hastings
*                 U, V, AND W
                  DO I = NSFX+1, NSFX+3
                     IF (LD(I) .EQ. 1) THEN
                        IPAR(IP,ICYC)=LP(I)
                        IP=IP+1
                     END IF
                  END DO
*                 X (LORENTZIAN SCHERRER BROADENING)
                  IF (LD(NSFX+5) .EQ. 1) THEN
                     IPAR(IP,ICYC)=LP(NSFX+5)
                     IP=IP+1
                  END IF
*                 Y (STRAIN BROADENING)
                  IF (LD(NSFX+7) .EQ. 1) THEN
                     IPAR(IP,ICYC)=LP(NSFX+7)
                     IP=IP+1
                  END IF
               CASE (1:3)
*                 Split-type pseudo-Voigt and Pearson VII functions
*                 U, V, AND W
                  DO I = NSFX+1, NSFX+3
                     IF (LD(I) .EQ. 1) THEN
                        IPAR(IP,ICYC)=LP(I)
                        IP=IP+1
                     END IF
                  END DO
*                 a0, a1, and a2
                  DO I = NSFX+5, NSFX+7
                     IF (LD(I) .EQ. 1) THEN
                        IPAR(IP,ICYC)=LP(I)
                        IP=IP+1
                     END IF
                  END DO
            END SELECT
         END IF
      END DO
      NPAR(ICYC)=IP-1

*     3RD CYCLE: SECONDARY PROFILE PARAMETERS #2
*     Xe, Ye, eL0, mL0, eH0, mH0
      IF (NPAR(ICYC) .GT. 0) ICYC=ICYC+1
      IP=1
      IS=0
      DO L = 1, NPHASE
         IF (NMODE .NE. 1 .AND. NMODE .NE. 6) THEN
            SELECT CASE (NPRFN)
               CASE (0)
*                 Anisotropic Scherrer broadening (Xe AND Ye)
*                 Xe (ANISOTROPIC SCHERRER BROADENING)
                  IF (LD(NSFX+6) .EQ. 1) THEN
                     IPAR(IP,ICYC)=LP(NSFX+6)
                     IP=IP+1
                  END IF
*                 Ye (ANISOTROPIC STRAIN BROADENING)
                  IF (LD(NSFX+8) .EQ. 1) THEN
                     IPAR(IP,ICYC)=LP(NSFX+8)
                     IP=IP+1
                  END IF
               CASE (1:3)
*                 Split-type pseudo-Voigt and Pearson VII functions 
*                 Dumping factors eta_L0, or m_L0
                  IF (LD(NSFX+9) .EQ. 1) THEN
                     IPAR(IP,ICYC)=LP(9)
                     IP=IP+1
                     IS=IS+1
                  ELSE
                     IS=IS
                  END IF
*                 Dumping factors eta_H0 or m_H0
                  IF (LD(NSFX+11) .EQ. 1) THEN
                     IPAR(IP,ICYC)=LP(11)
                     IP=IP+1
                     IS=IS+1
                  ELSE
                     IS=IS
                  END IF
            END SELECT
         END IF
      END DO
      NPAR(ICYC)=IP-1

*     4TH CYCLE: SECONDARY PROFILE PARAMETERS #3
*     rs, As, eL1, mL1, eH1, mH1
      IF (NPAR(ICYC) .GT. 0) ICYC=ICYC+1
      IP=1
      IS=0
      DO L = 1, NPHASE
         IF (NMODE .NE. 1 .AND. NMODE .NE. 6) THEN
            SELECT CASE (NPRFN)
               CASE (0)
*                 FINGER'S ASYMMETRY PARAMETER rs OR 
*                 HOWARD'S ASYMMETRY PARAMETER As
*                 (Source width)/(detector distance) rs OR As
                  IF (LD(NSFX+9) .EQ. 1) THEN
                     IPAR(IP,ICYC)=LP(NSFX+9)
                     IP=IP+1
                     IS=-1
                  END IF
               CASE (1:3)
*                 Split-type profile function
*                 Dumping factors etaL_1 or mL_1
                  IF (LD(NSFX+10) .EQ. 1) THEN
                     IPAR(IP,ICYC)=LP(10)
                     IP=IP+1
                     IS=IS+1
                  ELSE
                     IS=IS
                  END IF
*                 Dumping factors etaH_1 or mH_1
                  IF (LD(NSFX+12) .EQ. 1) THEN
                     IPAR(IP,ICYC)=LP(12)
                     IP=IP+1
                     IS=IS+1
                  ELSE
                     IS=IS
                  END IF
            END SELECT
         END IF
      END DO
      NPAR(ICYC)=IP-1

*     5TH CYCLE: SECONDARY PROFILE PARAMETERS #4
      IF (NPAR(ICYC) .GT. 0) ICYC=ICYC+1
      IP=1
      DO L = 1, NPHASE
         IF (NMODE .NE. 1 .AND. NMODE .NE. 6 .AND. IS .NE. 0) THEN
            SELECT CASE (NPRFN)
               CASE (0) 
                  IF (NASYM .EQ. 0 .AND. LD(NSFX+10) .EQ. 1 .AND. IS 
     &            .EQ. -1) THEN
*                    FINGER'S ASYMMETRY PARAMETER rd
*                    (Detector width)/(detector distance), rd
                     IPAR(IP,ICYC)=LP(NSFX+10)
                     IP=IP+1
                  ELSE IF (NASYM .EQ. 1) THEN
*                    Pseudo-Voigt function of Thompson, Cox, and Hastings
*                    U, V, AND W
                     DO I = NSFX+1,NSFX+3
                        IF (LD(I) .EQ. 1) THEN
                           IPAR(IP,ICYC)=LP(I)
                           IP=IP+1
                        END IF
                     END DO
*                    HOWARD'S ASYMMETRY PARAMETER As
                     IF (LD(NSFX+9) .EQ. 1) THEN
                        IPAR(IP,ICYC)=LP(NSFX+9)
                        IP=IP+1
                     END IF
                  END IF
               CASE (1:3)
*                 Split-type pseudo-Voigt and Pearson VII functions
*                 U, V, AND W
                  DO I = NSFX+1,NSFX+3
                     IF (LD(I) .EQ. 1) THEN
                        IPAR(IP,ICYC)=LP(I)
                        IP=IP+1
                     END IF
                  END DO
*                 a0, a1, and a2
                  DO I = NSFX+5, NSFX+7
                     IF (LD(I) .EQ. 1) THEN
                        IPAR(IP,ICYC)=LP(I)
                        IP=IP+1
                     END IF
                  END DO
*                 Dumping factors eta or m
                  DO I = NSFX+9, NSFX+12
                     IF (LD(I) .EQ. 1) THEN
                        IPAR(IP,ICYC)=LP(I)
                        IP=IP+1
                     END IF
                  END DO
            END SELECT
         END IF
      END DO
      NPAR(ICYC)=IP-1

*     6TH CYCLE: LATTICE CONSTANTS
      IF (NPAR(ICYC) .GT. 0) ICYC=ICYC+1
      IP=1
      DO L = 1, NPHASE
         DO I = NAX, NGAMX
            IF (LD(I) .EQ. 1) THEN
               IPAR(IP,ICYC)=LP(I)
               IP=IP+1
            END IF
         END DO
      END DO
      NPAR(ICYC)=IP-1

*     7TH CYCLE: PREFERRED-ORIENTATION PARAMETER 2 PLUS SCALE FACTOR
      IF (NPAR(ICYC) .GT. 0) ICYC=ICYC+1
      IP=1
      DO L = 1, NPHASE
         IF (LD(NPO1X) .EQ. 0 .AND. LD(NPO2X) .EQ. 0) CYCLE
*        SCALE FACTOR
         IF (LD(NSFX) .EQ. 1) THEN
            IPAR(IP,ICYC)=LP(NSFX)
            IP=IP+1
         END IF
         IF (LD(NPO1X) .EQ. 1) THEN
            IPAR(IP,ICYC)=LP(NPO1X)
            IP=IP+1
         END IF
         IF (LD(NPO2X) .EQ. 1) THEN
            IPAR(IP,ICYC)=LP(NPO2X)
            IP=IP+1
         END IF
      END DO
      NPAR(ICYC)=IP-1

*     8TH CYCLE: ATOMIC AND MAGNETIC STRUCTURAL PARAMETERS
*     g, x, y, z, B(or beta(ij)), mu phi
      IF (NPAR(ICYC) .GT. 0) ICYC=ICYC+1
      IP=1
      DO L = 1, NPHASE
         DO I = NOTX, KPHE(L)
            IF (LD(I) .EQ. 1) THEN
               IPAR(IP,ICYC)=LP(I)
               IP=IP+1
            END IF
         END DO
      END DO
      NPAR(ICYC)=IP-1

*     9TH CYCLE: PRIMARY PROFILE PARAMETERS #1
*     FWHM (OR FWHM(G), FWHM(L)), As
      IF (NPAR(ICYC) .GT. 0) ICYC=ICYC+1
      IP=1
      DO L = 1, NPHASE
         IF (NMODE .NE. 1 .AND. NPRIM .GT. 0) THEN
            SELECT CASE (NPRFN)
               CASE (1,3)
                  DO I = NPPPX, NPPPX + (NRELAX - 1)*NPRIM, NPRIM 
*                    FWHM
                     IF (ID(I) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I
                        IP=IP+1
                     END IF
*                    As
                     IF (ID(I+1) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+1
                        IP=IP+1
                     END IF
                  END DO
               CASE (2)
                  DO I = NPPPX, NPPPX + (NRELAX - 1)*NPRIM, NPRIM 
*                    FWHM(L)
                     IF (ID(I) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I
                        IP=IP+1
                     END IF
*                    FWHM(G)
                     IF (ID(I+1) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+1
                        IP=IP+1
                     END IF
*                    As
                     IF (ID(I+2) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+2
                        IP=IP+1
                     END IF
                  END DO
            END SELECT
         END IF
      END DO
      NPAR(ICYC)=IP-1

*     10TH CYCLE: PRIMARY PROFILE PARAMETERS #2
*     etaL (or mL) and etaH (or mH)
      IF (NPAR(ICYC) .GT. 0) ICYC=ICYC+1
      IP=1
      IS=0
      DO L = 1, NPHASE
         IF (NMODE .NE. 1 .AND. NPRIM .GT. 0) THEN
            SELECT CASE (NPRFN)
               CASE (1,3)
                  DO I = NPPPX, NPPPX + (NRELAX - 1)*NPRIM, NPRIM 
*                    etaL or mL
                     IF (ID(I+2) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+2
                        IP=IP+1
                        IS=IS+1
                     END IF
*                    etaH or mH
                     IF (ID(I+3) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+3
                        IP=IP+1
                        IS=IS+1
                     END IF
                     IF (ID(I+2) .NE. 1 .AND. ID(I+3) .NE. 1) THEN
                        IS=IS
                     END IF
                  END DO
               CASE (2)
                  DO I = NPPPX, NPPPX + (NRELAX - 1)*NPRIM, NPRIM 
*                    etaL
                     IF (ID(I+3) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+3
                        IP=IP+1
                        IS=IS+1
                     END IF
*                    etaH
                     IF (ID(I+4) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+4
                        IP=IP+1
                        IS=IS+1
                     END IF
                     IF (ID(I+3) .NE. 1 .AND. ID(I+4) .NE. 1) THEN
                        IS=IS
                     END IF
                  END DO
            END SELECT
         END IF
      END DO
      NPAR(ICYC)=IP-1

*     11TH CYCLE: All PRIMARY PROFILE PARAMETERS, |Fc| AND 2-THETA(PEAK) 
*     FWHM(OR FWHM(G), FWHM(L)), As, etaL(or mL) and etaH(or mH)
      IF (NPAR(ICYC) .GT. 0) ICYC=ICYC+1
      IP=1
      DO L = 1, NPHASE
         IF (NMODE .NE. 1 .AND. NPRIM .GT. 0 .AND. IS .NE. 0) THEN
            SELECT CASE (NPRFN)
               CASE (1,3)
                  DO I = NPPPX, NPPPX + (NRELAX - 1)*NPRIM, NPRIM 
*                    FWHM
                     IF (ID(I) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I
                        IP=IP+1
                     END IF
*                    As
                     IF (ID(I+1) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+1
                        IP=IP+1
                     END IF
*                    etaL or mL
                     IF (ID(I+2) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+2
                        IP=IP+1
                     END IF
*                    etaH or mH
                     IF (ID(I+3) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+3
                        IP=IP+1
                     END IF
*                    |Fc|
                     IF (NMODE .EQ. 3 .AND. ID(I+4) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+4
                        IP=IP+1
                     ELSE IF (NMODE .NE. 6) THEN
                        CYCLE
                     END IF
*                    const*|Fc|, 2-THETA(PEAK)
                     IF (ID(I+4) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+4
                        IP=IP+1
                     END IF
                     IF (ID(I+5) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+5
                        IP=IP+1
                     END IF
                  END DO
               CASE (2)
                  DO I = NPPPX, NPPPX + (NRELAX - 1)*NPRIM, NPRIM 
*                    FWHM(L)
                     IF (ID(I) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I
                        IP=IP+1
                     END IF
*                    FWHM(G)
                     IF (ID(I+1) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+1
                        IP=IP+1
                     END IF
*                    As
                     IF (ID(I+2) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+2
                        IP=IP+1
                     END IF
*                    etaL
                     IF (ID(I+3) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+3
                        IP=IP+1
                     END IF
*                    etaH
                     IF (ID(I+4) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+4
                        IP=IP+1
                     END IF
                     IF (NMODE .EQ. 3 .AND. ID(I+5) .EQ. 1) THEN
*                       |Fc|
                        IPAR(IP,ICYC)=I+5
                        IP=IP+1
                     ELSE IF (NMODE .NE. 6) THEN
                        CYCLE
                     END IF
*                    const*|Fc|, 2-THETA(PEAK)
                     IF (ID(I+5) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+5
                        IP=IP+1
                     END IF
                     IF (ID(I+6) .EQ. 1) THEN
                        IPAR(IP,ICYC)=I+6
                        IP=IP+1
                     END IF
* }
                  END DO
            END SELECT
         END IF
      END DO
      NPAR(ICYC)=IP-1

      IF (NPAR(ICYC) .GT. 0) THEN
         NCHNG=ICYC
      ELSE
         NCHNG=ICYC-1
      END IF
      END

************************************************************************

      SUBROUTINE RDPARA(NRFN,NCNSTR,NTERMS,A,ID,AOUT,IDSAVE,IR,
     &NMODE,FIRST,NSPLBL,LABEL,LPAR,NLINE)
*     READ LABELS, PARAMETERS, AND REFINEMENT INDICATORS
      PARAMETER (NT=999,NP=80000,NR=400,NCS=400,NPH=8,NAP=150)
      REAL A(*),FIRST(*),AOUT(*)
      INTEGER ID(*),IR(*),IDSAVE(*)
      INTEGER LPAR(*),NLINE
      INTEGER NSPLBL(*)
      CHARACTER LABEL(*)*25,LAB1*25,LAB2*25
      COMMON /ABSVR/ ABSORP(NP),VARLEN(NP),ROUGH(NP),NSURFR
      COMMON /LEGEN/ MAXDEG
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      PARAMETER (NB=7000)
      COMMON /SAVEYP/ SYPEAK(NB),SAVEA(NT),NTOTPAR

C**   PARAMETER AND CONDITIONS OF REFINEMENT

C     (1) GLOBAL PARAMETERS
C              A(1): Z (NPRFN = 0) or t0 (NPRFN = 1-3).
C              A(2): Ds (NPRFN = 0) or t1 (NPRFN = 1-3).
C              A(3): Ts (NPRFN = 0) or t2 (NPRFN = 1-3).
C              A(4): Dummy (NPRFN = 0) or t3 (NPRFN = 1-3).
C              A(NROUGH)-A(NROUGH+3): Surface-roughness parameter
C              A(NBCKGR)-A(NBCKGR+NUMBG-1): Background parameters.

C     (2) PRIMARY PROFILE PARAMETERS
C         A(NPPPX)-A(NPPPX+NPRIM*NRELAX-1)

C     (2) PHASE-DEPENDENT PARAMETERS
C         J IS THE NUMBER USED AS THE DUMMY ARGUMENT IN STATEMENT
C         FUNCTIONS P
C         KPHB(1) = NPPPX + NPRIM*NRELAX
C         P(J)=A(KPHB(L)+J)
C         WHERE L DENOTES THE PHASE NUMBER.

C      J   NAME
C      0   Scale factor, s.
C      1   Gaussian FWHM parameter (NPRFN = 0) or 
C          FWHM parameter (NPRFN = 1-3), U.
C      2   Gaussian FWHM parameter (NPRFN = 0) or 
C          FWHM parameter (NPRFN = 1-3), V.
C      3   Gaussian FWHM parameter (NPRFN = 0) or 
C          FWHM parameter (NPRFN = 1-3), W.
C      4   P (NPRFN = 0) or dummy (NPRFN = 1-3).
C      5   X (NPRFN = 0) or a0 (NPRFN = 1-3).
C      6   Xe (NPRFN = 0) or a1 (NPRFN = 1-3).
C      7   Y (NPRFN = 0) or a2 (NPRFN = 1-3).
C      8   Ye (NPRFN = 0) or dummy (NPRFN = 1-3).
C      9   r_s (NPRFN = 0 and NASYM = 0), As (NPRFN = 0 and NASYM =1), 
C          eta_L0 (NPRFN = 1 and 2), or m_L0 (NPRFN =3).
C     10   r_D (NPRFN = 0 and NASYM = 0), Dummy (NPRFN = 0), 
C          eta_L1 (NPRFN = 1 and 2), or m_L1 (NPRFN =3).
C     11   Dummy (NPRFN = 0), eta_H0 (NPRFN = 1 and 2), or 
C          m_H0 (NPRFN =3).
C     12   Dummy (NPRFN = 0), eta_H1 (NPRFN = 1 and 2), or 
C          m_H1 (NPRFN =3).
C     13   Dummy (NPRFN = 0) or anisotropic strain broadening, Ue
C          (NPRFN = 1-3)
C     14   Dummy (NPRFN = 0) or anisotropic Scherrer broadening, Pe
C          (NPRFN = 1-3).
C     15   Preferred-orientation parameter, p1.
C     16   Preferred-orientation parameter, p2.
C     17   a --> (a*)**2.
C     18   b --> (b*)**2.
C     19   c --> (c*)**2.
C     20   alpha --> (b*)(c*)cos(alpha*).
C     21   beta  --> (c*)(a*)cos(beta*).
C     22   gamma --> (a*)(b*)cos(gamma*).
C     23   Overall isotropic displacement parameter, Q.
C     24-  Structure-factor parameters.
C
C     ORDER OF STRUCTURE-FACTOR PARAMETERS
C        Occupation factor (g).
C        Fractional coordinates (x, y, z).
C        Isotropic or anisotropic displacement parameter(s)
C        Magnetic moment in Bohr magneton, mu

C     PHI(S) ARE ADDED AT THE TAIL OF STRUCTURE-FACTOR PARAMETERS
C     IF MAGNETIC SCATTERING IS OBSERVED
C     COORDINATES ARE GIVEN FOR ONE POSITION ONLY FOR EACH SET OF
C     SYMMETRY-RELATED POSITIONS
C
C     ID(LL)=0: THIS PARAMETER IS FIXED.
C     ID(LL)=1: THIS PARAMETER IS REFINED.
C     ID(LL)=2: THIS PARAMETER IS CONSTRAINED.
C     IR: REFINABLE-PARAMETER NUMBER.
C     NRFN: NUMBER OF REFINABLE PARAMETERS.
C
      NRFN=0
      NCNSTR=0
*     INITIALIZE ID FOR SUBROUTINE MARGOT (IN THE CASE OF SIMULATION)
      DO I = 1, NT
         IDSAVE(I) = 0
      END DO
      CALL DECPAR(NTERMS,A,ID,NSPLBL,LABEL,LPAR,NLINE,NPHL)
      NTOTPAR = NTERMS

      DO I = 1, NLINE - 1
         ISL = INDEX(LABEL(I),'/')
         IF (ISL .GT. 0) THEN
            LAB1 = LABEL(I)(1:ISL-1)
         ELSE
            LAB1 = LABEL(I)
         END IF
         DO J = I + 1, NLINE
            ISL = INDEX(LABEL(J),'/')
            IF (ISL .GT. 0) THEN
               LAB2 = LABEL(J)(1:ISL-1)
            ELSE
               LAB2 = LABEL(J)
            END IF
            IF (LAB1 .EQ. LAB2) CALL JOBEND
     &      ('The same label has been input: '//LAB1)
         END DO
      END DO
      IF (NMODE .NE. 1) THEN
         DO LL = 1, NTERMS
            IF (ID(LL) .EQ. 2) THEN
               NCNSTR=NCNSTR+1
               IF (NCNSTR .GT. NCS) CALL JOBEND('The number of'//
     &         ' linear constraints has reached the maximum value')
            ELSE IF (ID(LL) .EQ. 1) THEN
               NRFN=NRFN+1
               IF (NRFN .EQ. NR) CALL JOBEND('The number of refinable'
     &         //' parameters has reached the maximum value')
               IR(LL)=NRFN
               FIRST(NRFN)=0.0
            ENDIF
            IDSAVE(LL)=ID(LL)
            AOUT(LL)=A(LL)
         END DO
      END IF
*     CHECK THE MAXIMUM NUMBER OF TERMS 
*     SCAN BACKGROUND PARAMETERS
*     ENDING NUMBER FOR THE BACKGROUND PARAMETERS = NBCKGR + NUMBG -1
*     STARTING NUMBER FOR THE BACKGROUND PARAMETERS = NBCKGR
      DO J = NBCKGR + NUMBG -1, NBCKGR, -1
         IF (NMODE .NE. 1 .AND. (A(J) .NE. 0.0 .OR. ID(J) .EQ. 1)) EXIT
         IF (NMODE .EQ. 1 .AND. A(J) .NE. 0.0) EXIT
      END DO
      MAXDEG = J - NBCKGR
      
*     CHECK PEAK-SHIFT PARAMETERS IN INDIVIDUAL PROFILE FITTING
      IF (NMODE .EQ. 6) THEN
         DO J = 1, 4
            IF (A(J) .NE. 0.0) CALL JOBEND('Peak-shift parameters, '//
     &      'A(1)-A(4), should be zero when NMODE = 6')
         END DO
      END IF

*     CHECK A AND ID FOR SURFACE-ROUGHNESS PARAMETER
      IF (NMODE .NE. 1) THEN
         SELECT CASE (NSURFR)
            CASE (0)
               DO J = NROUGH, NROUGH + 3
                  IF (ID(J) .NE. 0) CALL JOBEND('Surface-roughness '//
     &            'parameters should not be refined when NSURFR = 0')
               END DO
            CASE (2)
               DO J = NROUGH, NROUGH + 2
                  IF (ID(J) .NE. 0) CALL JOBEND('q0, q1, and q2 '//
     &            'should not be refined when NSURFR = 2')
               END DO
            CASE (3,4)
               DO J = NROUGH+2, NROUGH + 3
                  IF (ID(J) .NE. 0) CALL JOBEND('q2 and q3 '//
     &            'should not be refined when NSURFR = 3')
               END DO
         END SELECT
      END IF
      END

************************************************************************

      SUBROUTINE ACONV(NPHASE,LMAG,A,KPHB,KPHE,LAUEG)
*     CONVERT PARAMETER VALUES
      PARAMETER (NPH=8)
      INTEGER LMAG(*),KPHB(*),KPHE(*),LAUEG(*)
      REAL A(*)
      COMMON /O/ NMODE,NPAT,NSFF,INCMULT,NCONST
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /TENS/ G(6,NPH),GG(6,NPH)
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      DATA RD,RDSQR/0.01745329,0.0003046174/
      LP(J)=KPHB(L)+J

*     DEG ==> RAD (PEAK-SHIFT PARAMETERS)
      DO J = 1, NBCKGR - 1
         A(J) = A(J)*RD
      END DO
      
      IF (NPRIM .GT. 0) THEN
         DO J = NPPPX, NPPPX + (NRELAX - 1)*NPRIM, NPRIM
*           DEG ==> RAD (FWHM FOR A RELAXED REFLECTION)
            A(J) = A(J)*RD
            IF (NPRFN .EQ. 2) A(J+1) = A(J+1)*RD
            IF (NMODE .EQ. 6) THEN
*              2-theta(peak) in radian
               SELECT CASE (NPRFN)
                  CASE (1,3)
                     A(J+5) = A(J+5)*RD
                  CASE (2)
                     A(J+6) = A(J+6)*RD
               END SELECT
            END IF
         END DO
      END IF
      
***   IF A NEW PROFILE FUNCTION IS IMPLEMENTED, THE FOLLOWING PART
***   SHOULD BE MODIFIED.
      DO L = 1, NPHASE
*        DEG**2 ==> RAD**2 (U, V, W, AND P)
         DO JJ=LP(1),LP(4)
            A(JJ)=A(JJ)*RDSQR
         END DO
         IF (NPRFN .EQ. 0) THEN
*           DEG ==> RAD (X, Xe, Y, Ye, and As)
            DO JJ = LP(5), LP(9)
               IF (NASYM .EQ. 0 .AND. JJ .EQ. LP(9)) EXIT
               A(JJ) = A(JJ)*RD
            END DO
         END IF
         A(LP(NAX))=GG(1,L)
         A(LP(NBX))=GG(2,L)
         A(LP(NCX))=GG(3,L)
         A(LP(NALPX))=GG(4,L)
         A(LP(NBETX))=GG(5,L)
         A(LP(NGAMX))=GG(6,L)
         IF (LMAG(L) .NE. 1) CYCLE
C        DEG ==> RAD (PHI)
         IF (LAUEG(L) .LE. 5) THEN
            A(KPHE(L)-2)=A(KPHE(L)-2)*RD
            A(KPHE(L)-1)=A(KPHE(L)-1)*RD
         END IF
         A(KPHE(L))=A(KPHE(L))*RD
      END DO
      END

************************************************************************

      FUNCTION CALINT(DEG,K,A)
C     SUMMATION OF THE CONTRIBUTIONS FROM NEIGHBORING BRAGG REFLECTIONS
      PARAMETER (NP=80000,NB=7000,NPH=8,NAP=150)
      INTEGER H
      REAL A(*),DEG(*)
      COMMON /A/ H(NB),IK(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /ABSVR/ ABSORP(NP),VARLEN(NP),ROUGH(NP),NSURFR
      COMMON /B/ U(NB),NREFL
      COMMON /C/ BG(NP),BGINC(NP),DEGNOR(NP)
      COMMON /G/ I
      COMMON /I/ RDEG(2,NB)
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /L/ DEG1,DEG2
      COMMON /MC/ CTHM,NTRAN,PCOR
      COMMON /ORDER/ IPOINT(NB)
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /PR/ PROR1(NB),CPROR(NPH),FLP(NB),DPROR(NB),YPEAK(NB)
      COMMON /RLS/ NRANGE
      COMMON /T/ FF(NB),COEF(NPH)
      COMMON /WL/ NBEAM,XLMD,XLMDH,XLMD2H,R12,RLAMBD
      COMMON /Z/ PC,CHGPC
C     THE NUMBER OF THE STARTING REFLECTION IS SAVED IN VARIABLE ISTART.
      SAVE ISTART

      IF (K .EQ. 1) ISTART=1
      CALINT = 0.0
      DO J = ISTART, NREFL
         I=IPOINT(J)
         IF (DEG(K) .GT. RDEG(2,I)) THEN
            ISTART=J+1
         ELSE IF (DEG(K) .LT. RDEG(1,I)) THEN
            EXIT
         ELSE
            CALINT=CALINT+YPEAK(I)*PRFL(DEG(K),A)
         END IF
      END DO

      THETA = 0.5*DEG(K)
*     CALCULATE SURFACE-ROUGHNESS CORRECTIONS
      IF (NSURFR .EQ. 0) THEN
         ROUGH(K) = 1.0
      ELSE
         P = A(NROUGH)
         Q = A(NROUGH+1)
         R = A(NROUGH+2)
         T = A(NROUGH+3)
         SELECT CASE (NSURFR)
            CASE (1)
               ROUGH(K) = R*(1.0 - P*(EXP(-Q) - EXP(-Q/SIN(THETA)))) +
     &                    (1.0 - R)*(1.0 + T*(THETA - 1.570796))
            CASE (2)
               ROUGH(K) = 1.0 - T*(THETA - 1.570796)
            CASE (3)
               ROUGH(K) = 1.0 - P*(EXP(-Q) - EXP(-Q/SIN(THETA)))
            CASE (4)
               ROUGH(K) = 1.0 - P*Q*((1.0 - Q) +
     &                    (1.0 - Q/SIN(THETA))/SIN(THETA))
         END SELECT
      END IF
      
      IF (NRANGE .EQ. 0) THEN
         BG(K) = FLEGEN(A,DEGNOR(K))
      ELSE IF (NRANGE .EQ. 3) THEN
         BG(K) = BGINC(K)+FLEGEN(A,DEGNOR(K))
      END IF
      CALINT = ABSORP(K)*VARLEN(K)*ROUGH(K)*CALINT + BG(K)
      END

************************************************************************

      REAL FUNCTION FLEGEN(A,X)
*     CALCULATE THE BACKGROUND CORRECTION USING A LEGENDRE POLYNOMIAL
      REAL A(*)
      COMMON /LEGEN/ MAXDEG
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X

      CALL SLESUM(X,MAXDEG,A(NBCKGR),FLEGEN)
      END

************************************************************************

      FUNCTION PRFL(TWOTH,A)
C     TWO PROFILE FUNCTIONS
*     DTTH: DELTA(2-THETA)
      PARAMETER (NB=7000,NPH=8,NAP=150)
      REAL A(*)
      COMMON /A/ H(NB),K(NB),L(NB),L12(NB),NOPH(NB)
      COMMON /DPRFL/ DPRDT,SIGPART,GAMPART,DPRDC7,DPRDS,DPRDD
      COMMON /DPRSPV/ DFDX,DFDA,DFDW,DFDW1,DFDW2,DFDEL,DFDEH,DFDRL,DFDRH
      COMMON /G/ I
      COMMON /J/ PEAK(NB),TH(NB),RLV2(NB),D(NB)
      COMMON /PHASE/ NPHASE,KPHB(NPH),KPHE(NPH),NPSITE(NAP,NPH)
      COMMON /PROFI/ PCOT(NB),SHFT(NB),NTSIM(NB),FNORM(NB),TNTSIM(NB),
     &  FWHM(NB),ETA(NB),SQSG(NB),CTH(NB),SIGP(NB),DFWDG(NB),DEDFFG(NB),
     &  DFWDL(NB),DEDFFL(NB),SIGPAP(NB)
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      COMMON /ZSPFN/ PEAKSH(NB),SQH2(NB),SQH2L(NB),SQH2G(NB),ASYM(NB),
     &ETAL(NB),ETAH(NB),FNORM1(2,NB),COMLG(2,NB),COML(2,NB),COMG(2,NB),
     &COEFL(2,NB),COEFG(2,NB),DECAY(2,NB),COMPEAR(2,NB),CDX(2,NB),
     &CDA1(2,NB),CDA2(2,NB),CDW(2,NB),CDRL1(2,NB),CDRL2(2,NB),
     &CDRH1(2,NB),DMPART(2,NB)
      COMMON /UNCONS/ NUCREF,IPHASE(90),IHKL(3,90),NOPPP(90),LPPP(NB),
     &  RCUT(2,90),NCUT
      LOGICAL USE_ASYM
*     RD22 = (180/PI)**2
      DATA RD2, RD22/57.29578, 3282.806/
      P(J) = A(KPHB(NOPH(I))+J)

      SELECT CASE (NPRFN)
         CASE (0)
*           PSEUDO-VOIGT FUNCTION OF THOMPSON, COX, AND HASTINGS
            IF (NASYM .EQ. 0) THEN
               USE_ASYM = .TRUE.
*              Unit: degree for Gamma, TwoTH, and TwoTH0, 
*              degree**(-1) for dPRdT and dPRdG. 
*              PIRORIN(Eta, Gamma, S_L, D_L, TwoTH, TwoTH0, 
*    &         dPRdT, dPRdG, dPRdE, dPRdS, dPRdD, Use_Asym)
*              With PIRORIN, the integrated intensity is normalized with 
*              the unit of the abscissa in degrees.
               PRFL = RD2*PIRORIN(ETA(I), FWHM(I)*RD2, P(9), P(10),
     &         TWOTH*RD2, (PEAK(I)-SHFT(I))*RD2, DPRDT, DPRDG, DPRDE, 
     &         DPRDS, DPRDD, USE_ASYM)
               DPRDT = -DPRDT*RD22
               DPRDG = DPRDG*RD22
               DPRDE = DPRDE*RD2
               DPRDS = DPRDS*RD2
               DPRDD = DPRDD*RD2
               SIGPART = DPRDG*SIGPAP(I)
               GAMPART = DPRDG*DFWDL(I)
            ELSE
               DTT = TWOTH - PEAK(I) + SHFT(I)
*              dH/d(sigma**2).  H is the summation of pseudo-Voigt 
*              functions.
               SIGPART = 0.0
*              dH/d(gamma)
               GAMPART = 0.0
*              dH/d(DELTA.T')
               DPRDT = 0.0
               DPRDC7 = 0.0
               PRFL = 0.0
*              MULTI-TERM SIMPSON'S RULE INTEGRATION DESCRIBED IN
*              C.J. HOWARD, J. APPL. CRYSTALLOGR. 15 (1982) 615.      
               DO IT = 1, NTSIM(I)
*                 SIMK: gi in GSAS (p. 122)
                  IF (MOD(IT,2) .EQ. 0) THEN
                     SIMK = 4.0
                  ELSE IF (IT .EQ. 1 .OR. IT .EQ. NTSIM(I)) THEN
                     SIMK = 1.0
                  ELSE
                     SIMK = 2.0
                  END IF
*                 SIMC: fi in GSAS (p. 122)
                  SIMC = FLOAT((IT-1)**2) / TNTSIM(I)
                  DT = DTT + PCOT(I)*SIMC
                  CALL PVOIGT(DT,FUNC,DFDT,DFDS,DFDG)
                  PRFL = PRFL + SIMK*FUNC
                  SIGPART = SIGPART + SIMK*DFDS
                  GAMPART = GAMPART + SIMK*DFDG
                  DFDT = SIMK * DFDT
                  DPRDT = DPRDT + DFDT
                  DPRDC7 = DPRDC7 + SIMC*DFDT
               END DO
               PRFL = FNORM(I)*PRFL
*              DERIVATIVE WITH RESPECT TO SIGMA**2
               SIGPART = FNORM(I) * SIGPART
               GAMPART = FNORM(I) * GAMPART
*              DPRDT: dH/d(DELTA.T')
               DPRDT = FNORM(I) * DPRDT
            END IF
         CASE (1,2)
*           (MODIFIED) SPLIT-TYPE PSEUDO-VOIGT FUNCTIONS
            CALL DETXJ(TWOTH - PEAK(I),PEAKSH(I),DTT,DTT2,J)
            IF (NPRFN .EQ. 1 .OR. (NPRFN .EQ. 2 .AND. LPPP(I) .EQ. 0)) 
     &      THEN
*              SPLIT-TYPE PEUDO-VOIGT FUNCTION
*              H. Toraya, J. Appl. Crystallogr. 23 (1990) 485.
               FLOR = 1.0/(1.0 + COMLG(J,I)*DTT2)
               FGAU = EXP(-0.6931472*COMLG(J,I)*DTT2)
               SUMLG = COEFL(J,I)*FLOR + COEFG(J,I)*FGAU
               SELECT CASE (J)
                  CASE (1)
                     CALL DERSPV(FNORM1(J,I),COMLG(J,I),SQH2(I),ASYM(I),
     &               ETAL(I),ETAH(I),DTT,FLOR,FGAU,SUMLG,DFDX,DFDA,DFDW,
     &               DFDEL,DFDEH)
                  CASE (2)
                     CALL DERSPV(FNORM1(J,I),COMLG(J,I),SQH2(I),
     &               1.0/ASYM(I),ETAH(I),ETAL(I),DTT,FLOR,FGAU,SUMLG,
     &               DFDX,DFDA,DFDW,DFDEH,DFDEL)
*                    A' = 1/A ==> DF/DA = (DF/DA')*(DA'/DA) = 
*                    -(DF/DA')/A**2
                     DFDA = -DFDA/ASYM(I)**2
               END SELECT
            ELSE
*              MODIFIED SPLIT-TYPE PSEUDO-VOIGT FUNCTION
               FLOR = 1.0/(1.0 + COML(J,I)*DTT2)
               FGAU = EXP(COMG(J,I)*DTT2)
               SUMLG = COEFL(J,I)*FLOR + COEFG(J,I)*FGAU
               SELECT CASE (J)
                  CASE (1)
                     CALL DERMSPV(FNORM1(J,I),COML(J,I),SQH2L(I),
     &               SQH2G(I),ASYM(I),ETAL(I),ETAH(I),DTT,FLOR,FGAU,
     &               SUMLG,DFDX,DFDA,DFDW1,DFDW2,DFDEL,DFDEH)
                  CASE (2)
                     CALL DERMSPV(FNORM1(J,I),COML(J,I),SQH2L(I),
     &               SQH2G(I),1.0/ASYM(I),ETAH(I),ETAL(I),DTT,FLOR,FGAU,
     &               SUMLG,DFDX,DFDA,DFDW1,DFDW2,DFDEH,DFDEL)
*                    A' = 1/A ==> DF/DA = (DF/DA')*(DA'/DA) = 
*                    -(DF/DA')/A**2
                     DFDA = -DFDA/ASYM(I)**2
               END SELECT
            END IF
            PRFL = FNORM1(J,I)*SUMLG
         CASE (3)
*           SPLIT-TYPE PEARSON VII FUNCTION
            CALL DETXJ(TWOTH - PEAK(I),PEAKSH(I),DTT,DTT2,J)
            PRFL = FNORM1(J,I)/(1.0 + COMPEAR(J,I)*DTT2)**DECAY(J,I)
            TEMP = PRFL
            CALL DERSPVII(DTT,DTT2,TEMP,CDX(J,I),CDA1(J,I),CDA2(J,I),
     &      CDW(J,I),SQH2(I),CDRL1(J,I),CDRL2(J,I),CDRH1(J,I),
     &      COMPEAR(J,I))
*           -(DF/DA')/A**2
            IF (J .EQ. 2) DFDA = -DFDA/ASYM(I)**2
      END SELECT
      END

************************************************************************

      SUBROUTINE DETXJ(DTTH,PEAKSH,DTT,DTT2,J)
*     DETERMINE DELTA_2-THETA AND J

      DTT = DTTH + PEAKSH
      DTT2 = DTT*DTT
      IF (DTT .LE. 0.0) THEN
*        LOWER ANGLE REGION
         J = 1
      ELSE
*        HIGHER ANGLE REGION
         J = 2
      END IF
      END

************************************************************************

      SUBROUTINE DERSPVII(X,X2,PRFL,CDX,CDA1,CDA2,CDW,W,CDRL1,CDRL2,
     &  CDRH1,COMPEAR)
*     CALCULATE THE DERIVATIVES OF THE SPLIT PEARSON VII FUNCTION
*     WITH RESPECT TO x, A, W, m_l, and m_h
      COMMON /DPRSPV/ DFDX,DFDA,DFDW,DFDW1,DFDW2,DFDEL,DFDEH,DFDRL,DFDRH
      
          T = 1.0 + COMPEAR*X2
       DFDX = -PRFL*CDX*X/T
       DFDA =  PRFL*(CDA1 + CDA2/T)
       DFDW =  PRFL*(CDW*X2/T - 1.0/W)
      DFDRL = -PRFL*(CDRL1 - LOG(T) - CDRL2*X2/T)
      DFDRH = -PRFL*CDRH1
      END

************************************************************************

      SUBROUTINE DERSPV(FNORM1,COMLG,W,A,ETAL,ETAH,X,FLOR,FGAU,SUMLG,
     &  DFDX,DFDA,DFDW,DFDEL,DFDEH)
*     CALCULATE THE DERIVATIVES OF THE SPLIT PSEUDO-VOIGT FUNCTION
*     WITH RESPECT TO x, A, W, eta_l, and eta_h

      DATA TPI/0.6366198/    ! 2.0/PI
      DATA RPILN2/1.4756646/ ! SQRT(PI*LN(2.0))
      DATA CGAU/0.9394373/   ! 2.0*SQRT(LN(2.0)/PI)
      DATA CDFN/-0.4756646/  ! 1.0 - SQRT(PI*LN(2.0))
      
*     4.0/PI = 1.2732395
*     4.0*LN(2.0)*SQRT(LN(2.0)/PI) = 1.3023366
      SUMCOM = X*(1.2732395*ETAL*FLOR**2 + 1.3023366*(1.0-ETAL)*FGAU)
      
      DFDX = -FNORM1*COMLG/W*SUMCOM
     
      DFDA = FNORM1 * ((-FNORM1 + 1.0)/(1.0 + A)*SUMLG +
     &                X*(1.0 + A)/(A*W)**3 * SUMCOM)
     
      DFDW = FNORM1/W**2* 
     &       (COMLG*X*SUMCOM - TPI*ETAL*FLOR - CGAU*(1.0 - ETAL)*FGAU)
     
*     'BUNBO' FOR THE NORMALIZATION FACTOR
      BUNBO = ETAL + (1.0 - ETAL)*RPILN2 + 
     &        A*(ETAH + (1.0 - ETAH)*RPILN2)
     
      DFDEL = FNORM1*(-CDFN/BUNBO*SUMLG + (TPI*FLOR - CGAU*FGAU)/W)
      
      DFDEH = CDFN/BUNBO*((1.0 + A) - A*FNORM1) * SUMLG
      END

************************************************************************

      SUBROUTINE DERMSPV(FNORM1,COML,W1,W2,A,ETAL,ETAH,X,FLOR,FGAU,
     &  SUMLG,DFDX,DFDA,DFDW1,DFDW2,DFDEL,DFDEH)
*     CALCULATE THE DERIVATIVES OF THE MODIFIED SPLIT PSEUDO-VOIGT 
*     FUNCTION WITH RESPECT TO x, A, W, eta_l, and eta_h

      DATA TPI/0.6366198/    ! 2.0/PI
      DATA PI/3.14159265/    ! PI
      DATA RPILN2/1.4756646/ ! SQRT(PI*LN(2.0))
      DATA CGAU/0.9394373/   ! 2.0*SQRT(LN(2.0)/PI)
      DATA CDFN/-0.4756646/  ! 1.0 - SQRT(PI*LN(2.0))
      
*     LN(2.0)*SQRT(LN(2.0)/PI) = 0.32558415
      COMXA = ETAL/(PI*W1**3)*FLOR**2 + (1.0-ETAL)*0.32558415/W2**3*FGAU
      
      DFDX = -FNORM1*4.0*((1.0+A)/A)**2*X*COMXA
     
      DFDA = FNORM1 * ((-FNORM1 + 1.0)/(1.0 + A)*SUMLG + 
     &       4.0*X**2*(1.0 + A)/A**3*COMXA)
      
      DFDW1 = FNORM1*ETAL/(PI*W1**2)*((2.0*FLOR*X)**2*COML - 2.0*FLOR)
     
*     SQRT(LN(2.0)/PI) = 0.46971864
*     4.0*LN(2.0) = 2.772589
      DFDW2 = FNORM1*(1.0-ETAL)*0.46971864/W2**2*
     &        (2.772589*(X*(1.0+A)/(A*W2))**2 - 2.0)*FGAU
     
*     'BUNBO' FOR THE NORMALIZATION FACTOR
      BUNBO = ETAL + (1.0 - ETAL)*RPILN2 + 
     &        A*(ETAH + (1.0 - ETAH)*RPILN2)
     
      DFDEL = FNORM1*(-CDFN/BUNBO*SUMLG + TPI/W1*FLOR - CGAU/W2*FGAU)
     
      DFDEH = CDFN/BUNBO*((1.0 + A) - A*FNORM1) * SUMLG
      END

************************************************************************

      SUBROUTINE PVOIGT(DT,FUNC,DFDT,DFDS,DFDG)
*     CALCULATE THE PSEUDO-VOIGT FUNCTION MADE ASYMMETRIC ACCORDING TO
*     THE PROCEDURE OF HOWARD (1982)
*     THE PSEUDO-VOIGT FUNCTION IS DESCRIBED IN P. THOMPSON, D.E. COX, 
*     AND J.B. HASTINGS, J. APPL. CRYSTALLOGR. 20 (1987) 79. 
      PARAMETER (NB=7000)
      COMMON /G/ I
      COMMON /PROFI/ PCOT(NB),SHFT(NB),NTSIM(NB),FNORM(NB),TNTSIM(NB),
     &  FWHM(NB),ETA(NB),SQSG(NB),CTH(NB),SIGP(NB),DFWDG(NB),DEDFFG(NB),
     &  DFWDL(NB),DEDFFL(NB),SIGPAP(NB)

      BUNBO = (0.5*FWHM(I))**2 + DT*DT
*     LORENTZIAN COMPONENT
      TL = FWHM(I) / (6.283185*BUNBO)

*     DTLDT: dL/d(DELTA.T')
      DTLDT = -2.0*DT*TL/BUNBO

*     DTLDFW: dL/d(GAMMA)
      DTLDFW = TL/FWHM(I) - 3.141593*TL*TL

      EX = MAX(-20.0, -0.5*DT*DT/SIGP(I))

*     GAUSSIAN COMPONENT
*     SQRT(4*LN(2)/PI) = 0.9394373
      TG = 0.9394373 * EXP(EX) / FWHM(I) 

*     FUNC: F(DELTA.T')
      FUNC = ETA(I)*TL + (1.0 - ETA(I))*TG

*     TS: (1 - eta) * dG/d(GAMMA)
      TS = -2.0*(1.0 - ETA(I))*TG*(EX + 0.5)/FWHM(I)

*     DFDT: dF/d(DELTA.T').  In SUBROUTINE PSVOIGT, DFDX      
      DFDT = ETA(I)*DTLDT - (1.0 - ETA(I))*TG*DT/SIGP(I)

*     dF/d(GAMMA.g)
      DFDS = DEDFFG(I)*(TL - TG) + (ETA(I)*DTLDFW + TS)*DFWDG(I)

*     DFDS: dF/d(sigma**2).  1.17741 = 0.5*SQRT(8*LN(2))
      DFDS = 1.17741 * DFDS / SQSG(I)

*     DFDG: dF/d(gamma) 
      DFDG = DEDFFL(I)*(TL - TG) + (ETA(I)*DTLDFW + TS)*DFWDL(I)
      END

************************************************************************

      SUBROUTINE TYPEREF(L,REFTYP,LCON)
*     TYPE OF POSSIBLE REFLECTIONS
      PARAMETER (NPH=8)
      CHARACTER REFTYP(NPH)*3,LCON(7,NPH)*80,COL(80),LINE*80
      COMMON /CONHKL/ NHKL(NPH),ICH(5,3,7,NPH),NCON(7,NPH),MLINE(NPH)
      EQUIVALENCE (LINE,COL(1))

      REFTYP(L)='hkl'
      READ(4,'(A)') LINE
      DO J = 1, 80
         IF (LINE(J:J) .NE. ' ') EXIT
      END DO
      REFTYP(L) = LINE(J:)

*     SELECT CASE using a character type is not supported in g77
C     SELECT CASE (REFTYP(L))
C        CASE ('hkl')
C           NHKL(L)=1
C        CASE ('0kl')
C           NHKL(L)=2
C        CASE ('h0l') 
C           NHKL(L)=3
C        CASE ('hk0')
C           NHKL(L)=4
C        CASE ('hhl') 
C           NHKL(L)=5
*        CASE ('h-hl')
*           NHKL(L)=6
C        CASE DEFAULT
C           CALL JOBEND('Bad reflection type: '//REFTYP(L))
C     END SELECT
      IF (REFTYP(L) .EQ. 'hkl') THEN
         NHKL(L)=1
      ELSE IF (REFTYP(L) .EQ. '0kl') THEN
         NHKL(L)=2
      ELSE IF (REFTYP(L) .EQ. 'h0l') THEN
         NHKL(L)=3
      ELSE IF (REFTYP(L) .EQ. 'hk0') THEN
         NHKL(L)=4
      ELSE IF (REFTYP(L) .EQ. 'hhl') THEN
         NHKL(L)=5
*     ELSE IF (REFTYP(L) .EQ. ('h-hl') THEN
*        NHKL(L)=6
      ELSE
         CALL JOBEND('Bad reflection type: '//REFTYP(L))
      END IF

      DO ICOND = 1, 7
C*       CONDITIONS FOR POSSIBLE REFLECTIONS
C*       END OF CONDITIONS: '}'
         READ(4,'(A)') LINE
         IF (LINE(1:1) .EQ. '}') EXIT
         LCON(ICOND,L)=LINE
         CALL POSREF(COL,ICH(1,1,ICOND,L),NCON(ICOND,L))
      END DO
      MLINE(L) = ICOND -1
      END

************************************************************************

      SUBROUTINE ANGRANGE(A,ANGMAX,ANGMIN)
*     DETERMINE ANGULAR RANGE FOR KDRREF IN INDIVIDUAL PROFILE FITTING
      REAL A(*)
      COMMON /L/ DEG1,DEG2
      COMMON /PARNO/ NROUGH,NBCKGR,NUMBG,NPPPX,NRELAX,NPRIM,NSFX,NPO1X,
     &  NPO2X,NAX,NBX,NCX,NALPX,NBETX,NGAMX,NOTX,NG1X
      COMMON /ZPRFN/ NPRFN,NSHIFT,NASYM
      DATA RD, RD2/0.01745329, 57.29578/
     
      ANGMIN = DEG1
      ANGMAX = DEG2
      DO IS = NPPPX, NPPPX + (NRELAX - 1)*NPRIM, NPRIM
         SELECT CASE (NPRFN)
            CASE (1,3)
*              SPLIT-TYPE PSEUDO-VOIGT OR PEARSON VII FUNCTION
               TWOTH = A(IS+5)
            CASE (2)
*              MODIFIED SPLIT-TYPE PSEUDO-VOIGT FUCNTION
               TWOTH = A(IS+6)
         END SELECT
         IF (TWOTH .LT. DEG1) THEN
            ANGMIN = TWOTH - 0.03*RD
         ELSE IF (TWOTH .GT. DEG2) THEN
            ANGMAX = TWOTH + 0.03*RD
         END IF
      END DO
      ANGMIN = ANGMIN*RD2
      ANGMAX = ANGMAX*RD2
      END
