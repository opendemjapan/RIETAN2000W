ECHO OFF
CLS
REM Batch file to run ORFFE and browse the resulting output file
REM Change three environment variables in lines starting with 'SET' appropriately.
REM Double quotation marks are unnecessary even when spaces are included in folder/file names.

REM Absolute path for the location of *.xyz
REM 'LOCXYZ=.' if *.xyz is contained in the same directory as this batch file
SET LOCXYZ=C:\RIETAN-2000 Folder\templates
REM Base name.  Name of the input file = (base name).xyz.
SET SAMPLE=FapatiteE
ECHO The sample name is '%SAMPLE%'.
REM Absolute path for the location of orffe.exe (the same as that of rietan.exe)
SET ORFFE=C:\RIETAN-2000 Folder\programs

REM Change the current directory
IF "%OS%" == "Windows_NT" CHDIR /D "%LOCXYZ%"
IF NOT "%OS%" == "Windows_NT" "%ORFFE%\PUSHD" "%LOCXYZ%"

ECHO ORFFE is now running ...
REM  Input *.xyz: ORFFE instructions.
REM Output *.dst: Sandard output.
REM        *.ffe: Data for imposing restraints on interatomic distances and/or bond angles (a copy of *.dst).
"%ORFFE%\orffe.exe" "%SAMPLE%.xyz" | "%ORFFE%\tee.exe" "%SAMPLE%.dst"

IF "%OS%" == "Windows_NT" Wscript "%ORFFE%\message2.vbs"
IF NOT "%OS%" == "Windows_NT" Wscript "%ORFFE%\message1.vbs"
