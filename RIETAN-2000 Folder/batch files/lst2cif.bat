ECHO OFF
CLS
REM Batch file to run lst2cif
REM Change environment variables LOCLST, SAMPLE, and LST2CIF in lines starting with 'SET' appropriately.
REM Double quotation marks are unnecessary even when spaces are included in folder/file names.

REM Absolute path for the location of *.ins and *.int
REM 'LOCLST=.' if *.ins is contained in the same directory as this batch file
SET LOCLST=C:\RIETAN-2000 Folder\templates
SET SAMPLE=FapatiteE
REM Absolute path for the location of lst2cif
SET LST2CIF=C:\RIETAN-2000 Folder\programs
REM Absolute path for the location of spgri, spgra, and tee.exe
SET DATABASE=%LST2CIF%

ECHO Data files are located in '%LOCLST%'.
ECHO The base name is '%SAMPLE%'.
REM Change the current directory
IF "%OS%" == "Windows_NT" CHDIR /D "%LOCLST%"
IF NOT "%OS%" == "Windows_NT" "%LST2CIF%\PUSHD" "%LOCLST%"

REM  Input *.lst: Output file of RIETAN.
REM Output *.cif: CIF.
REM        *.l2c: Printer output.
"%LST2CIF%\lst2cif.exe" "%SAMPLE%.lst" | "%DATABASE%\tee.exe" "%SAMPLE%.l2c"
DEL "%SAMPLE%.l2c"
PAUSE
