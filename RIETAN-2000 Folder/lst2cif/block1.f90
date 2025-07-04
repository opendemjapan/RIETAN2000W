      SUBROUTINE WRITE_BLOCK1(OUT)
C
      INTEGER OUT
C--------------------------------------------------------------------
C     SUBMISSION DETAILS
C--------------------------------------------------------------------
      WRITE(OUT,900)
      WRITE(OUT,902)
      WRITE(OUT,900)
      WRITE(OUT,100)
      WRITE(OUT,900)
      WRITE(OUT,900)
      WRITE(OUT,101)
      WRITE(OUT,102)
      WRITE(OUT,901)
      WRITE(OUT,901)
      WRITE(OUT,103)
      WRITE(OUT,104)
      WRITE(OUT,105)
      WRITE(OUT,900)
      WRITE(OUT,106)
      WRITE(OUT,107)
      WRITE(OUT,900)
      WRITE(OUT,108)
      WRITE(OUT,901)
      WRITE(OUT,901)
      WRITE(OUT,900)
C--------------------------------------------------------------------
C
C--------------------------------------------------------------------
C     FORMATS FOR  SUBMISSION DETAILS
C--------------------------------------------------------------------
 100  FORMAT('# SUBMISSION DETAILS')
 101  FORMAT('_publ_contact_author_name',10X,"'",'?',"'")
 102  FORMAT('_publ_contact_author_address')
 103  FORMAT('_publ_contact_author_email',9X,'?')
 104  FORMAT('_publ_contact_author_fax',11X,"'",'?',"'")
 105  FORMAT('_publ_contact_author_phone',9X,"'",'?',"'")
 106  FORMAT('_publ_requested_journal',12X,"'",'?',"'")
 107  FORMAT('_publ_requested_category',11X,'?')
 108  FORMAT('_publ_contact_letter')
C--------------------------------------------------------------------
C
C
 900  FORMAT(' ')
 901  FORMAT(';')
 902  FORMAT('#',70('='))
 903  FORMAT('loop_')
C
      END
