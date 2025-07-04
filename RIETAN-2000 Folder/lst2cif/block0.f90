      SUBROUTINE WRITE_BLOCK0(OUT)
      INTEGER OUT
      CHARACTER DATE_TIME*16

C --------------------------------------------------------------------
C     MAIN BLOCK
C--------------------------------------------------------------------
      CALL GET_DATE_TIME(DATE_TIME)
      WRITE(OUT,902)
      WRITE(OUT,900)
      WRITE(OUT,001)
      WRITE(OUT,900)
      WRITE(OUT,902)
      WRITE(OUT,900)
* Corrected by F. Izumi on 2001.02.21 
C     WRITE(OUT,002)
C     WRITE(OUT,003) DATE_TIME
      WRITE(OUT,004) DATE_TIME(1:10)
* Correction end
      WRITE(OUT,005) 
* Corrected by F. Izumi on 2001.02.20 
C     WRITE(OUT,006)
C     WRITE(OUT,007) DATE_TIME
C     WRITE(OUT,901)
C     WRITE(OUT,900)
* Correction end
C
C--------------------------------------------------------------------
C     FORMATS FOR MAIN BLOCK
C--------------------------------------------------------------------
 001  FORMAT('data_RIETAN_publ')
 002  FORMAT('_pd_block_id')
 003  FORMAT(6X,A,'|','RIETAN_overall','|','..creator_name..','|',
     &       '..instr_name..')
 004  FORMAT('_audit_creation_date',2X,A)
* Corrected by F. Izumi on 2001.02.20 
C005  FORMAT('_audit_creation_method',2X,"""from *.lst file using RTOCIF
C    &""")
 005  FORMAT('_audit_creation_method',2X,'''Converted from *.lst using l
     &st2cif''')
* Correction end
 006  FORMAT('_audit_update_record')
* Corrected by F. Izumi on 2001.02.20 
C007  FORMAT(';',1X,A,16X,'Initial CIF as created by RTOCIF')
 007  FORMAT(';',1X,A,16X,'Created with lst2cif')
* Correction end
C--------------------------------------------------------------------
C
 900  FORMAT(' ')
 901  FORMAT(';')
 902  FORMAT('#',70('='))
 903  FORMAT('loop_')
C
      END
