REM Record the abosolute path for rietan.exe
ECHO OFF
CLS

SET SETTINGS=%USERPROFILE%\Settings_RIETAN
IF NOT EXIST "%SETTINGS%\" (
   MKDIR "%SETTINGS%"
) ELSE (
  DEL /Q "%SETTINGS%\*"
)

FOR /F "DELIMS=""" %%I IN ('DIR ..\programs\rietan.exe /B /S') DO (
ECHO %%~dpI>"%SETTINGS%\abs_path_rietan"
SET RIETAN=%%~dpI
)
ECHO The absolute path for rietan.exe is %RIETAN%
ECHO It has been recorded in %SETTINGS%\abs_path_rietan
Wscript "%RIETAN%message2.vbs"
