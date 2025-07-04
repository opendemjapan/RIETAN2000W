      SUBROUTINE WRITE_BLOCK2(OUT)
C
      INTEGER OUT
C
C--------------------------------------------------------------------
C     PROCESSING SUMMARY
C--------------------------------------------------------------------
      WRITE(OUT,900)
      WRITE(OUT,902)
      WRITE(OUT,900)
      WRITE(OUT,200)
      WRITE(OUT,900)
      WRITE(OUT,201)
      WRITE(OUT,202)
      WRITE(OUT,203)
      WRITE(OUT,204)
      WRITE(OUT,205)
      WRITE(OUT,206)
      WRITE(OUT,207)
      WRITE(OUT,208)
      WRITE(OUT,209)
      WRITE(OUT,210)
      WRITE(OUT,211)
      WRITE(OUT,212)
      WRITE(OUT,213)
      WRITE(OUT,214)
      WRITE(OUT,215)
      WRITE(OUT,216)
      WRITE(OUT,217)
      WRITE(OUT,218)
      WRITE(OUT,219)
      WRITE(OUT,220)
      WRITE(OUT,221)
      WRITE(OUT,222)
      WRITE(OUT,223)
      WRITE(OUT,224)
      WRITE(OUT,900)
C--------------------------------------------------------------------
C
C--------------------------------------------------------------------
C     FORMATS FOR PROCESSING SUMMARY
C--------------------------------------------------------------------
 200  FORMAT('# PROCESSING SUMMARY (IUCr Office Use Only)')
 201  FORMAT('#_journal_date_recd_electronic')
 202  FORMAT('#_journal_date_to_coeditor')
 203  FORMAT('#_journal_date_from_coeditor')
 204  FORMAT('#_journal_date_accepted')
 205  FORMAT('#_journal_date_printers_first')
 206  FORMAT('#_journal_date_printers_final')
 207  FORMAT('#_journal_date_proofs_out')
 208  FORMAT('#_journal_date_proofs_in')
 209  FORMAT('#_journal_coeditor_name')
 210  FORMAT('#_journal_coeditor_code')
 211  FORMAT('#_journal_coeditor_notes')
 212  FORMAT('#_journal_techeditor_code')
 213  FORMAT('#_journal_paper_category')
 214  FORMAT('#_journal_compatibility_tag')
 215  FORMAT('#_journal_techeditor_notes')
 216  FORMAT('#_journal_coden_ASTM')
 217  FORMAT('#_journal_name_full')
 218  FORMAT('#_journal_year')
 219  FORMAT('#_journal_volume')
 220  FORMAT('#_journal_issue')
 221  FORMAT('#_journal_page_first')
 222  FORMAT('#_journal_page_last')
 223  FORMAT('#_journal_suppl_publ_number')
 224  FORMAT('#_journal_suppl_publ_pages')
C--------------------------------------------------------------------
C
 900  FORMAT(' ')
 901  FORMAT(';')
 902  FORMAT('#',70('='))
 903  FORMAT('loop_')
C
      END
