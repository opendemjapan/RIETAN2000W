REM Batch file to run RIETAN-2000 and display the resulting *.pat if any
ECHO OFF
CLS

REM Change environment variables, LOCINS, SAMPLE, RIETAN, and PATVIEWER in lines starting with 'SET' appropriately.
REM No folder/file/base names need to be enclosed by a pair of double quotation marks even when they include spaces.

REM When the total length of environment variables is large in Windows 95/98/Me, the initial size for environment variables must be increased.
REM Richt-click the icon of sample.bat, select 'Properties', click the 'Memory' tab, and change the size for environment variables from 'Auto' to 512 or 768.

REM Absolute path for the location of *.ins and *.int
REM 'LOCINS=.' if *.ins is contained in the same directory as this batch file
SET LOCINS=C:\RIETAN-2000 Folder\templates
REM Base name.  File name = (base name).*, where '*' denotes an extension such as 'ins' and 'int'.
SET SAMPLE=FapatiteE
REM Absolute path for the location of rietan.exe, orffe.exe, etc.
SET RIETAN=C:\RIETAN-2000 Folder\programs
REM An application (other than gnuplot) to plot Rietveld-refinement and simulated patterns
SET PATVIEWER=E:\Applications\Igo Pro\Igor.exe

ECHO Data files are located in '%LOCINS%'.
ECHO The base name is '%SAMPLE%'.
REM Change the current directory
IF "%OS%" == "Windows_NT" CHDIR /D "%LOCINS%"
IF NOT "%OS%" == "Windows_NT" "%RIETAN%\pushd" "%LOCINS%"

IF EXIST "%SAMPLE%.pat" RENAME "%SAMPLE%.pat" "%SAMPLE%.sav"

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
"%RIETAN%\rietan.exe" "%SAMPLE%.ins" "%SAMPLE%.int" "%SAMPLE%.bkg" "%SAMPLE%.pat" "%SAMPLE%.hkl" "%SAMPLE%.xyz" "%SAMPLE%.mem" "%SAMPLE%.ffe" "%SAMPLE%.fba" "%SAMPLE%.ffi" "%SAMPLE%.ffo" "%SAMPLE%.vcs" | "%RIETAN%\tee.exe" "%SAMPLE%.lst"

IF NOT EXIST "%SAMPLE%.pat" GOTO NPAT0
REM The base name of *.ins should not contain a space because %SAMPLE%.plt cannot be enclosed by " ".
IF EXIST "%SAMPLE%.plt" "%RIETAN%\Gnuplot\wgnuplot.exe" /noend %SAMPLE%.plt
IF NOT EXIST "%SAMPLE%.plt" "%PATVIEWER%" "%SAMPLE%.pat"
IF EXIST "%SAMPLE%.sav" DEL "%SAMPLE%.sav"

:NPAT0
IF "%OS%" == "Windows_NT" Wscript "%RIETAN%\message2.vbs"
IF NOT "%OS%" == "Windows_NT" Wscript "%RIETAN%\message1.vbs"
