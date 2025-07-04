      SUBROUTINE WRITE_BLOCK5(OUT,TMP)
C
      INTEGER OUT,TMP
      CHARACTER CMD*300
      CHARACTER DATE_TIME*16

	LOGICAL RESULT
C
C--------------------------------------------------------------------
C     POWDER SPECIMEN AND EXPERIMENTAL DATA
C--------------------------------------------------------------------
      CALL GET_DATE_TIME(DATE_TIME)
      CMD= 'GENERAL_BEGIN'
      CALL FINDSTRING(0,TMP,CMD,13,RESULT)
      CALL SKIPLINE(TMP,1)
      READ(TMP,*) METHOD
      IF (METHOD.EQ.1) THEN
         READ(TMP,*) NBEAM
         READ(TMP,*) WL
      END IF
      READ(TMP,*) RADIUS
      READ(TMP,*) DENSTY 
      CMD= '###'
      CALL FINDSTRING(1,TMP,CMD,13,RESULT)
      CALL SKIPLINE(TMP,1)

	IF (METHOD.EQ.1) THEN
        READ(TMP,*) THMIN
        READ(TMP,*) THMAX
        READ(TMP,*) STEP
        READ(TMP,*) NSTEP
      ENDIF
      WRITE(OUT,900)
      WRITE(OUT,902)
      WRITE(OUT,900)
      WRITE(OUT,500)
      WRITE(OUT,900)
      WRITE(OUT,501)
      WRITE(OUT,900)
      WRITE(OUT,502)
* Corrected by F. Izumi on 2001.02.21 
C     WRITE(OUT,503) DATE_TIME
      WRITE(OUT,503) DATE_TIME(1:10)
* Correction end
      WRITE(OUT,504)
      WRITE(OUT,505)
      WRITE(OUT,506)
      WRITE(OUT,507)
      WRITE(OUT,901)
      WRITE(OUT,901)
      WRITE(OUT,900)
      WRITE(OUT,529)
      WRITE(OUT,508)
      WRITE(OUT,509)
      WRITE(OUT,510)
      WRITE(OUT,511)
      WRITE(OUT,512)
      WRITE(OUT,513)
      WRITE(OUT,514)
      WRITE(OUT,515)
      IF (METHOD.EQ.0) THEN

         WRITE(OUT,516) 'tof'

	ELSE
         WRITE(OUT,516) 'step'

	ENDIF

      WRITE(OUT,517)
      WRITE(OUT,901)
      WRITE(OUT,901)
C      WRITE(OUT,900)
      IF (METHOD.EQ.1) THEN
C         WRITE(OUT,903)
         WRITE(OUT,518) "'Cu K\\a~1~'"
         WRITE(OUT,519) WL
C         WRITE(OUT,518)
C         WRITE(OUT,519)
C         WRITE(OUT,520) 'Cu K\\a~1~',WL
C         WRITE(OUT,900)
         WRITE(OUT,521)
         WRITE(OUT,900)
         WRITE(CMD,'(F7.3)') THMIN
         DO I=1,7
            IF (CMD(I:I).NE.' ') GOTO 1 
         END DO 
 1       I1=I
         WRITE(OUT,525) CMD(I1:I1+8)
         WRITE(CMD,'(F7.3)') THMAX
         DO I=1,7
            IF (CMD(I:I).NE.' ') GOTO 2 
         END DO 
 2       I1=I
         WRITE(OUT,526) CMD(I1:I1+8)
         WRITE(CMD,'(F7.3)') STEP
         DO I=1,7
            IF (CMD(I:I).NE.' ') GOTO 3 
         END DO 
 3       I1=I
         WRITE(OUT,527) CMD(I1:I1+8)
         WRITE(CMD,*) NSTEP
         DO I=1,10
            IF (CMD(I:I).NE.' ') GOTO 4 
         END DO 
 4       I1=I
         WRITE(OUT,528) CMD(I1:20)
      ELSE
         WRITE(OUT,522)
         WRITE(OUT,523)
         WRITE(OUT,524)
      END IF
      WRITE(OUT,900)


C--------------------------------------------------------------------
C     FORMATS FOR  POWDER SPECIMEN AND EXPERIMENTAL DATA
C--------------------------------------------------------------------
 500  FORMAT('# POWDER SPECIMEN AND EXPERIMENTAL DATA') 
 501  FORMAT('data_RIETAN_p_01') 
 502  FORMAT('_pd_block_id') 
 503  FORMAT(6X,"'",A,'|','POWSET_01','|','..creator_name..','|',
     &       '..instr_name..',"'")
 504  FORMAT('_pd_meas_datetime_initialed',12X,'?')
 505  FORMAT('_pd_meas_info_author_name',14X,"""?""")
 506  FORMAT('_pd_meas_info_author_email',13X,'?')
 507  FORMAT('_pd_meas_info_author_address')

 529  FORMAT('_pd_calc_method',23X,"""Rietveld Refinement""")
 508  FORMAT('_diffrn_ambient_temperature',11X,'?')
 509  FORMAT('_diffrn_ambient_environment',11X,'?')
 510  FORMAT('_diffrn_source',24X,"'",'?',"'")
 511  FORMAT('_diffrn_source_target',17X,'?')
 512  FORMAT('_diffrn_source_type',19X,'?')
 513  FORMAT('_diffrn_measurement_device_type',7X,"'",'?',"'")
 514  FORMAT('_diffrn_detector',22X,"'",'?',"'")
 515  FORMAT('_diffrn_detector_type',17X,'?')
 516  FORMAT('_pd_meas_scan_method',18X,A)
 517  FORMAT('_pd_meas_special_details')
 518  FORMAT('_diffrn_radiation_type',16X,A)
 519  FORMAT('_diffrn_radiation_wavelength',9X,F7.5)
 520  FORMAT(6X,A,/6X,F8.5) 
 521  FORMAT('_diffrn_radiation_monochromator',6X,'none')
 522  FORMAT('_pd_instr_dist_src/spec',12X,'?')
 523  FORMAT('_pd_instr_dist_src/detc',12X,'?')
 524  FORMAT('_pd_meas_2theta_fixed',14X,'?')
 525  FORMAT('_pd_meas_2theta_range_min',8X,A)
 526  FORMAT('_pd_meas_2theta_range_max',8X,A)
 527  FORMAT('_pd_meas_2theta_range_inc',8X,A)
 528  FORMAT('_pd_meas_number_of_points',8X,A)
C--------------------------------------------------------------------
C
 900  FORMAT(' ')
 901  FORMAT(';')
 902  FORMAT('#',70('='))
 903  FORMAT('loop_')
 904  FORMAT('#',70('-'))
C
      END



C--------------------------------------------------------------------
C--------------------------------------------------------------------
      SUBROUTINE WRITE_BLOCK6(OUT,TMP)
C
      CHARACTER CMD*300
      INTEGER OUT,TMP

	LOGICAL RESULT
C--------------------------------------------------------------------
C     REFINEMENT DATA
C--------------------------------------------------------------------
      CMD= 'RFACTOR_BEGIN'
      CALL FINDSTRING(0,TMP,CMD,13,RESULT)
      CALL SKIPLINE(TMP,1)
* Corrected by F. Izumi on 2001.02.20 
* In the case of multi-phase analysis, R factors are not written correctly.
* For convenience, this program is stopped here.
C     READ(TMP,'(4F6.2)') WRF,RF,QQ,WRE
      READ(TMP,'(4F6.2)',ERR=9) WRF,RF,QQ,WRE
	GO TO 8
    9 STOP
    8 CONTINUE
      GFE=WRF/WRE
      CMD= 'CONSTR_BEGIN'
      CALL FINDSTRING(1,TMP,CMD,12,RESULT)
      CALL SKIPLINE(TMP,1)
      READ(TMP,*) NCNSTR

      WRITE(OUT,900)
      WRITE(OUT,902)
      WRITE(OUT,900)
      WRITE(OUT,600)
      WRITE(OUT,900)
      WRITE(OUT,601)
      WRITE(OUT,901)
      WRITE(OUT,901)
      WRITE(OUT,900)
      WRITE(OUT,602)
      WRITE(OUT,603)
      WRITE(OUT,604)
      WRITE(OUT,901)
      WRITE(OUT,901)
      WRITE(OUT,900)
      WRITE(CMD,'(F6.3)') RF/100.0
      DO I=1,6
         IF (CMD(I:I).NE.' ') GOTO 1
      END DO
 1    I1=I
      WRITE(OUT,605) CMD(I1:I1+10)
      WRITE(CMD,'(F6.3)') WRF/100.0
      DO I=1,6
         IF (CMD(I:I).NE.' ') GOTO 2
      END DO
 2    I1=I
      WRITE(OUT,606) CMD(I1:I1+10)
      WRITE(CMD,'(F6.3)') WRE/100.0
      DO I=1,6
         IF (CMD(I:I).NE.' ') GOTO 3
      END DO
 3    I1=I
      WRITE(OUT,607) CMD(I1:I1+10)
      WRITE(OUT,608)
      WRITE(OUT,901)
      WRITE(OUT,901)
      WRITE(OUT,900)
      WRITE(OUT,609)
      WRITE(OUT,610)
      WRITE(OUT,611)
      WRITE(OUT,612)
      WRITE(OUT,613)
      WRITE(OUT,614)
      WRITE(OUT,615)
      WRITE(CMD,*) NCNSTR
      DO I=1,10
         IF (CMD(I:I).NE.' ') GOTO 4
      END DO
 4    I1=I
      WRITE(OUT,616) CMD(I1:I1+10)
      WRITE(CMD,'(F5.2)') GFE
      DO I=1,10
         IF (CMD(I:I).NE.' ') GOTO 5
      END DO
 5    I1=I
      WRITE(OUT,617) CMD(I1:I1+10)

      WRITE(OUT,900)
      

C
C--------------------------------------------------------------------
C     FORMATS FOR  POWDER SPECIMEN AND EXPERIMENTAL DATA
C--------------------------------------------------------------------
 600  FORMAT('# REFINEMENT DATA') 
 601  FORMAT('_pd_proc_ls_special_details') 
 602  FORMAT('_pd_proc_ls_profile_function',6X,'?') 
 603  FORMAT('_pd_proc_ls_background_function',3X,"'",'?',"'") 
 604  FORMAT('_pd_proc_ls_pref_orient_corr') 

 605  FORMAT('_pd_proc_ls_prof_R_factor',9X,A) 
 606  FORMAT('_pd_proc_ls_prof_wR_factor',8X,A) 
 607  FORMAT('_pd_proc_ls_prof_wR_expected',6X,A) 
 608  FORMAT('_refine_special_details') 

 609  FORMAT('_refine_ls_structure_factor_coef',4X,'Inet') 
 610  FORMAT('_refine_ls_matrix_type',14X,'?') 
 611  FORMAT('_refine_ls_weighting_sheme',10X,"'",'?',"'") 
 612  FORMAT('_refine_ls_hydrogen_treatment',7X,'noref') 
 613  FORMAT('_refine_ls_extinction_method',8X,'none') 
 614  FORMAT('_refine_ls_extinction_coef',10X,'?') 

 615  FORMAT('_refine_ls_number_parameters',8X,'?') 
 616  FORMAT('_refine_ls_number_constraints',7X,A) 
 617  FORMAT('_refine_ls_goodness_of_fit_all',6X,A) 
C--------------------------------------------------------------------
C
 900  FORMAT(' ')
 901  FORMAT(';')
 902  FORMAT('#',70('='))
 903  FORMAT('loop_')
 904  FORMAT('#',70('-'))
C
C
      END






