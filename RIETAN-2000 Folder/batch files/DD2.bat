REM Batch file to run RIETAN-2000, browse the output, calculate interatomic distances and bond angles by ORFFE, plot the result of Rietveld refinement or simulation, and open *.ins on Windows 2000/XP.

ECHO OFF
CLS
REM If you enter 'set' in the Command Prompt, all environment variables including 'OS' are displayed.
IF "%OS%" == "Windows_NT" GOTO NT2kXP
Echo This batch file cannot be used on Windows 95/98/Me.
PAUSE
GOTO END

:NT2kXP
REM A directory where three setting files are stored
SET SETTINGS=%USERPROFILE%\Settings_RIETAN

IF NOT EXIST "%SETTINGS%\abs_path_rietan" (
   Echo Execute setupDD.bat before using DD.bat for the first time.
   PAUSE
   GOTO END
)

REM Environment variable RIETAN: location of rietan.exe, orffe.exe, etc.
REM Note that the tail of RIETAN is '\'
PUSHD "%SETTINGS%"
FOR /F "DELIMS=""" %%I IN (abs_path_rietan) DO SET RIETAN=%%I
POPD

REM This batch file has been clicked
REM Placing "%1" in the next line causes an error when *.ins was dragged and dropped onto DD.bat.
REM The error probably occurs because %1 contains a pair of double quotation marks.
IF "%~1" == "" (
   IF EXIST "%SETTINGS%\base_name_ins" (
      REM Read the base name and absolue path for *.ins
      PUSHD "%SETTINGS%"
      FOR /F "DELIMS=""" %%I IN (base_name_ins) DO SET SAMPLE=%%I
      FOR /F "DELIMS=""" %%I IN (abs_path_ins) DO (
         POPD
         CHDIR /D %%I
      )
      GOTO Path_pat_ins
   ) ELSE (
      GOTO Exec_path_ins
   )
)

REM Drag & drop *.ins
REM The full (drive+path+file) name of %1 may be changed into an abbreviated name.
REM The full name can be obtained by 'DIR %1 /B /S'.
REM SAMPLE: Base name for the file that has been dragged & dropped
REM File name = (base name).*, where '*' denotes an extension
REM %1 contains a pair of double quotation marks in this case.
IF EXIST %1 (
   FOR /F "DELIMS=""" %%I IN ('DIR %1 /B /S') DO SET BASE_NAME=%%~nI
) ELSE (
   REM Only the base name has been input in the command prompt: Exriet (base name)
   FOR /F "DELIMS=""" %%I IN ('DIR %1.ins /B /S') DO SET BASE_NAME=%%~nI
)

REM Register *.ins
REM No space should be placed just before '>' to avoid insertion of ' ' before CR and LF
ECHO %BASE_NAME%>"%SETTINGS%\base_name_ins"
SET SAMPLE=%BASE_NAME%
IF EXIST %1 (
   FOR /F "DELIMS=""" %%I IN ('DIR %1 /B /S') DO ECHO %%~dpI>"%SETTINGS%\abs_path_ins"
) ELSE (
   REM Only the base name has been input
   FOR /F "DELIMS=""" %%I IN ('DIR %1.ins /B /S') DO ECHO %%~dpI>"%SETTINGS%\abs_path_ins"
)
CHDIR /D %~dp1
GOTO Path_pat_ins

:Exec_path_ins
"%RIETAN%path_ins.exe"
REM Specification of *.ins has been canceled
PUSHD "%HOMEDRIVE%\
IF NOT EXIST temp_base GOTO Quit

FOR /F "DELIMS=""" %%I IN (temp_base) DO (
   ECHO %%I>"%SETTINGS%\base_name_ins"
   SET SAMPLE=%%I
   DEL temp_base
)

FOR /F "DELIMS=""" %%I IN (temp_path) DO (
   ECHO %%I>"%SETTINGS%\abs_path_ins"
   DEL temp_path
   POPD
   CHDIR /D %%I
)

:Path_pat_ins
REM Add '.pat' and '.ins' to PATHEXT for opening *.pat and *.ins automatically
SET PATHEXT=%PATHEXT%; .pat; .ins

REM Move *.pat and *.xyz temporarily
IF EXIST "%SAMPLE%.pat" MOVE "%SAMPLE%.pat" "%SETTINGS%"
IF EXIST "%SAMPLE%.xyz" MOVE "%SAMPLE%.xyz" "%SETTINGS%"

ECHO RIETAN-2000 is now running ...
REM  Input *.ins: Standard input (record length of 80).
REM        *.int: X-ray/neutron diffraction data.
REM        *.bkg: Background intensities.
REM        *.ffe: Input data created by ORFFE for imposing constraints on interatomic distances and/or bond angles.
REM        *.fba: Data created by MEED for MEM-based whole-pattern fitting.
REM        *.ffi: Initial integrated intensities for Le Bail refinement.
REM Output *.pat: Data for plotting Rietveld-refinement patterns or a simulated pattern.
REM        *.hkl: Data for Fourier/D synthesis by FOUSYN.
REM        *.xyz: Data for calculating interatomic distances and bond angles by ORFFE.
REM        *.mem: Data for MEM analysis by MEED.
REM        *.ffo: Integrated intensities resulting from Le Bail refinement.
REM        *.vcs: VICS (VIsualization of Crystal Structures) text file.
REM        *.lst: Standard output.
"%RIETAN%rietan.exe" "%SAMPLE%.ins" "%SAMPLE%.int" "%SAMPLE%.bkg" "%SAMPLE%.pat" "%SAMPLE%.hkl" "%SAMPLE%.xyz" "%SAMPLE%.mem" "%SAMPLE%.ffe" "%SAMPLE%.fba" "%SAMPLE%.ffi" "%SAMPLE%.ffo" "%SAMPLE%.vcs" | "%RIETAN%tee.exe" "%SAMPLE%.lst"

REM Execute orffe.exe, which should be located in the same folder as rietan.exe
IF EXIST "%SAMPLE%.xyz" (
   "%RIETAN%orffe.exe" "%SAMPLE%.xyz" | "%RIETAN%tee.exe" "%SAMPLE%.dst"
   IF EXIST "%SETTINGS%\%SAMPLE%.xyz" DEL "%SETTINGS%\%SAMPLE%.xyz"
) ELSE (
   IF EXIST "%SETTINGS%\%SAMPLE%.xyz" MOVE "%SETTINGS%\%SAMPLE%.xyz" .
)

REM Display *.pat
REM In what follows, file names after "START" are titles displayed in window title bars.
REM Enter "help START" in the Command Prompt to learn arguments of START.
REM The base name of *.ins should not contain a space because %SAMPLE%.plt cannot be enclosed by " ".
IF EXIST "%SAMPLE%.pat" (
   REM If *.plt exist, not Igor.exe but wgnuplot.exe is launched.
   IF EXIST "%SAMPLE%.plt" (
      REM %SAMPLE%.plt should not be enclosed by " ".
      START "%SAMPLE%.plt" /B "%RIETAN%Gnuplot\wgnuplot.exe" /noend %SAMPLE%.plt
      START "%SAMPLE%.plt" /B "%SAMPLE%.plt"
   ) ELSE (
      START "%SAMPLE%.pat" /B "%SAMPLE%.pat"
   )
   IF EXIST "%SETTINGS%\%SAMPLE%.pat" DEL "%SETTINGS%\%SAMPLE%.pat"
) ELSE (
   IF EXIST "%SETTINGS%\%SAMPLE%.pat" MOVE %SETTINGS%\%SAMPLE%.pat" .
)

REM Display *.ins for editing it
START "%SAMPLE%.ins" /B "%SAMPLE%.ins"

REM Browse the output file(s) using less.exe
IF EXIST "%SAMPLE%.xyz" (
   "%RIETAN%less.exe" -MIS "%SAMPLE%.lst" "%SAMPLE%.dst"
) ELSE (
   "%RIETAN%less.exe" -MIS "%SAMPLE%.lst"
)
GOTO END

:Quit
Wscript "%RIETAN%message2.vbs"
:END
