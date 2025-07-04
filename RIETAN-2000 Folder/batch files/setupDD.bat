REM Record the abosolute path for rietan.exe
ECHO OFF
CLS

SET SETTINGS=%HOMEDRIVE%%HOMEPATH%Settings_RIETAN
REM Refer to p. 322 in Akiyama (2000) for a technique using NUL
IF NOT EXIST %SETTINGS%\NUL (
   MKDIR %SETTINGS%
) ELSE (
  DEL /Q %SETTINGS%\*
)

FOR /F "DELIMS=""" %%I IN ('DIR ..\programs\rietan.exe /B /S') DO (
ECHO %%~dpI>%SETTINGS%\abs_path_rietan
SET RIETAN=%%~dpI
)
ECHO The absolute path for rietan.exe is %RIETAN%%.
ECHO It has been recorded in %SETTINGS%\abs_path_rietan.
Wscript "%RIETAN%message2.vbs"
