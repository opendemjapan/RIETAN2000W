Installing and Executing RIETAN-2000, ORFFE, and Other Related Programs on Microsoft Windows

RIETAN-2000 is a multi-purpose pattern-fitting system.

ORFFE is a program for calculating geometric structure parameters such as interatomic distances and bond angles:

W. R. Busing, K. O. Martin, and H. A. Levy, Report ORNL-TM-306, Oak Ridge National Laboratory, Tennessee (1964).

Its original program was modified in such a way that five columns are assigned to each atom designation (refer to the tail of each user input file, *.ins).

The current versions of RIETAN-2000 and ORFFE for Microsoft Windows (hereafter abbreviated to Windows) have been compiled and linked with Compaq Visual Fortran 6.6B.

They will be able to be run on IBM-PC compatible machines (CPU: Pentium II/Pro or later) where Windows 95/98/Me/NT/2000/XP have been installed but cannot be run on 80486 machines (obsolete!).

The use of mission-critical operating systems, Windows NT/2000/XP, is strongly recommended because of the instability of Windows 95/98/Me and the ability of the MS-DOS Prompt much lower than that of the Command Prompt.  For example, the concurrent execution of RIETAN-2000 and tee.exe is practically impossible with Windows 95/98/Me.

In what follows, 'folder' (in Windows) may be replaced with 'directory' (in the DOS/Command Prompt), and vice versa.

1. How to install RIETAN-2000

At first, select a mode where extensions for files, e.g., 'ins', 'lst', and 'pat', are explicitly displayed, which is set by turning off 'Hide file extensions for known file types' at 'View' in 'Folder Options' under 'Tools' in the case of Windows 2000.

RIETAN-2000 can be installed according to the following procedure:

1.1 Archive file

If you have installed a previous version of RIETAN-2000, deleting the whole content is recommended to avoid problems arising the mixing of old and new files.  In addition, please use new template files because the order and number of input data may be changed.

Extracting the archive file, rietan2000w.tbz, in an appropriate folder (e.g., C:\Program Files), you will get a folder named 'RIETAN-2000 Folder' there.  Six subfolders (programs, rietan2000, orffe, lst2cif, templates, batch files), a text file (readme_win.txt, i.e., this document), reflection.ipf, and wgnuplot.ini are included in 'RIETAN-2000 Folder'.

Beware lest 'RIETAN-2000 Folder' is moved into a folder with a very long absolute path name on the use of a batch file, DD2.bat, described later, DD2.bat causes troubles under such a situation. For example, DD2.bat works well if 'RIETAN-2000 Folder' is put in a typical folder of C:\Program Files\.

The second file, reflection.ipf, stores a macro named reflection for Igor Pro.  Move it to folder Igor Procedures under folder Igor Pro.  Then, an item, reflection, appears when clicking the Macros pull-down menu in Igor Pro.  Refer to Section 5 for its operation.

The third file, wgnuplot.ini, is a preferences file for gnuplot 3.7.1+1.2.0 described below.  Move it to the Windows directory, e.g., C:\Winnt40 in Windows NT and C:\Winnt in Windows 2000.  Refer to Section 7 to get detailed information about gnuplot.

Contents of the five folders:

(a) programs
 1) rietan.exe: Executable binary program of RIETAN-2000 for Intel Pentium II/Pro or later.
 2) rietan_pn4.exe: Executable binary program of RIETAN-2000 optimized for Intel Pentium 4.
 3) asfdc: File storing atomic scattering factors, scattering lengths, etc.
 4) spgra: File storing information on 230 space groups described in International Tables, Vol. A.
 5) spgri: File storing information on 230 space groups described in International Tables, Vol. I.
 6) orffe.exe: Executable binary program of ORFFE.
 7) lst2cif.exe: Executable binary program of lst2cif.
 8) tee.exe: UNIX-compatible command.
 9) less.exe: Command contained in the GNU utility.
10) pushd.exe: UNIX-compatible command.
11) path_ins.exe: Executable binary program used in DD.bat or DD2.bat described in (e).
12) messagei.vbs (i = 1-4): One-line scripts written in VBScript.  They are called in batch files described in (e).
13) Gnuplot: A folder containing gnuplot 3.7.1+1.2.0 for Windows and files associated with it.

Never transfer the above files and folder 'Gnuplot' from folder 'programs'.  In adition, never change the file and folder name except for rietan.exe.

The four executable binary programs, rietan.exe, rietan_pn4.exe, orffe.exe, and lst2cif.exe, were obtained with Compaq Visual Fortran, Ver. 6.6B.

If you have a PC whose CPU is Pentium 4, rename the two executable files of RIETAN0-2000 in the following manner:
rietan.exe --> rietan_pn2.exe
rietan_pn4.exe --> rietan.exe

To my regret, the current version of Compaq Visual Fortran does not yield a stable executable program of RIETAN-2000 optimized for AMD Athlon.

Free software less.exe (v3.81) and tee.exe were obtained from "GNU utilities for Win32" at http://unxutils.sourceforge.net/.  The GNU utility command, less.exe, is used to browse the output file, *.lst as well as *.dst if any.  The UNIX-compatible command, tee.exe, enables us to run an application such as RIETAN-2000 while browsing its output scrollng in the Command Prompt window.

Pushd.exe (Ver. 1.1) is free software distributed by Kenji Arai (http://www2k.biglobe.ne.jp/~araken/pushd.htm).  It is required to change both the current directory and drive simultaneously in the DOS prompt on Windows 95/98/Me (pushd.exe is included in Windows NT/2000/XP).  It is used not in DD.bat or DD2.bat but in RIETAN-2000.bat described in (c).

Path_ins.exe was obtained with Absoft Pro Fortran for Windows, Ver. 7.0.  No source program is included in this distribution.  Contact the author by to get the source program.

Gnuplot 3.7.1+1.2.0 enhanced by Masahito Yamaga is available for download from http://www.yama-ga.com/gnuplot/.  Because he permits its redistribution, it was included in the 'programs' folder for the convenience of users of RIETAN-2000 for Windows.  If you want to redistribute it further, please read Copyright_plus.txt in the Gnuplot folder carefully.

(b) rietan2000
Source program of RIETAN-2000 (*.f90) and a makefile for obtaining rietan.exe (placed in folder programs) from *.f90 (orffe.mak).

The source program is essentially the same as that for UNIX/Linux except that
USE DFLIB
USE DFPORT
are inserted in windos.f90 to call UNIX-compatible subprograms.

(c) orffe
Source program of ORFFE (*.f90) and a makefile for obtaining orffe.exe (placed in folder programs) from *.f90 (orffe.mak).  Refer to Section 3.

When using Windows NT/2000/XP, you can change various settings including the size and number of the buffer, window size, layout, font, colors of the background, and font for the Command Prompt window by right-clicking its title bar and selecting 'Properties...'.  The buffer size must be sufficiently large on the use of DD.bat (DD2.bat) because both outputs from RIETAN-2000 and ORFFE are displayed in the same window.  The width of 140 would be appropriate for the buffer and window.  Raster Fonts with a size of 8x12 is highly recommended.

(d) lst2cif
Source program of lst2cif (*.f90) and a makefile for obtaining lst2cif.exe (placed in folder programs) from *.f90 (lst2cif.mak).  Refer to Sections 4.

(d) templates
Template files (*.ins) and intensity data files (*.int) for neutron, conventional X-ray, and synchrotron X-ray powder diffraction are contained together with a gnuplot scripts folder.

Template files whose names contain 'E' and 'J' just before periods are those where comments are written in English and Japanese, respectively.  Needless to say, contents of *E.int and *J.int are just the same with each other.  Three gnuplot script files for only *E.ins and *E.int are contained in the gnuplot scripts folder

Many 'hankaku-kana' characters are contained in the template files written in Japanese to decrease the number of characters in each line as many as possible.  It should be pointed out that they are converted into 'zenkaku-kana' if these files are attached to mails without compressing them.  This causes exceeding the limit of 80 columns for each line, which leads to errors.

For details in the gnuplot script files, refer to Section 7.

1) Neutron powder diffraction
Tl2223E.ins and Tl2223J.ins for Rietveld analysis with the corresponding intensity files Tl2223E.int and Tl2223J.int, respectively.

2) Powder diffraction with characteristic X rays (Cu K-alpha radiation)
FapatiteE.ins and FapatiteJ.ins for Rietveld analysis with the corresponding intensity files FapatiteE.int and FapatiteJ.int, respectively.

3) Powder diffraction with characteristic X rays (Cu K-alpha1 radiation)
Cu3Fe4P6E.ins and Cu3Fe4P6E.int for three-phase Rietveld analysis.  No template file written in Japanese is provided at present.  Only profile parameters for the pseudo-Voigt function of Thompson et al. and asymmetry correction of Howard are effective in Cu3Fe4P6E.ins.

BaSO4E.ins and BaSO4E.int for Le Bail analysis.  Neither template file written in Japanese nor script file for gnuplot is provided at present.  Only NPRFN = 2 is  effective in BaSO4E.ins.

4) Synchrotron X-ray powder diffraction
CimetidineE.ins and CimetidineJ.ins for Rietveld analysis with the corresponding intensity files CimetidineE.int and CimetidineJ.int, respectively.

(e) batch files
1) RIETAN-2000.bat: Batch file to execute RIETAN-2000 and Igor Pro successively (for Windows 95/98/Me).
2) setupDD.bat: Batch file to execute before using DD.bat at first (for Windows NT).
3) setupDD2.bat: Batch file to execute before using DD2.bat at first (for Windows 2000/XP).
4) DD.bat: Batch file to execute RIETAN-2000 and other programs successively (for Windows NT).
5) DD2.bat: Batch file to execute RIETAN-2000 and other programs successively (for Windows 2000/XP).
6) ORFFE.bat: Batch file to execute ORFFE.
7) lst2cif.bat: Batch file to execute lst2cif.

These batch files may be moved to anywhere including the desktop.

1.2 Modifying the batch files

(a) RIETAN-2000.bat

Copy RIETAN-2000.bat to any folder whithin the same drive as *.ins.  Of course, the folder containing sample.ins and sample.int is fine, where 'sample' is a sample name.  Rename RIETAN-2000.bat appropriately, e.g., sample.bat (this name will hereafter be used).

Modify sample.bat using a text editor.  This template batch file was prepared for FapatiteE.*.

The following four environment variables have to be changed:
1) LOCINS: (absolute) path for the location of sample.ins ('.' if sample.bat is located in the same folder as sample.ins).
2) SAMPLE: Base name without a period and an extension in *.ins and *.int (in the present case "sample").
3) RIETAN: Absolute path for the location of rietan.exe, asfdc, spgra, spgri, and tee.exe.
4) PATVIEWER: Absolute path + file name of an application (other than gnuplot) for plotting Rietveld-refinement and simulated patterns (sample.pat).

Edit only lines containing SET commands as required.

When the total length of environment variables is large in Windows 95/98/Me, the initial size for environment variables must be increased.  Richt-click the icon of sample.bat, select 'Properties', click the 'Memory' tab, and change the size for environment variables from 'Auto' to 512 or 768.

The following line is related to a graph-plotting application such as Igor Pro:

IF NOT EXIST "%SAMPLE%.plt" "%PATVIEWER%" "%SAMPLE%.pat"

If you have not installed any graphing software by yourself, comment out it by placing 'REM ' at the top of each line.

If this line is not commented out and PATVIEWER is correctly set, Igor Pro will be launched to read in sample.pat, if any, in the absence of sample.plt after execution of RIETAN-2000.  Even if a scientific graphing software other than Igor Pro is used, you can change the above line appropriately.

RIETAN-2000.bat is a mere template file; you may change it as you like.  For example, a batch file for the Command Prompt window is also fine for those who like a character user interface.

You may feel that such a manner of editing a batch file is too old-fashioned.  However, this makes it possible to run RIETAN without specifying an input file, *.ins, using an 'open file' dialog-box, which is rather tedious.  In addition, I hate this dialog-box because it is not very convenient and user-friendly in contrast with the corresponding one in the Mac OS.

Despite the use of tee.exe, the output from RIETAN-2000 is displayed not during its execution but after the end of the calculatoin on Windows 95/98/Me, which is a reason why Windows NT/2000/XP is preferred to Windows 95/98/Me to run RIETAN-2000.

(b) DD.bat and DD2.bat

DD.bat for Windows NT and DD2.bat for Windows 2000/XP are elaborate and powerful batch files to run RIETAN-2000, ORFFE, Igor Pro, and an editor successively.  If you hope to skip running ORFFE, Igor Pro, or the editor, change corresponding parts to comments.

Combination of rietan.exe and DD.bat/DD2.bat is a good guide to transportation of a Fortran program for a character user interface to Windows in a cool manner.  Because it fully utilizes features supported only in the Command Prompt, it cannot be used on Windows 95/98/Me providing us with the MS-DOS Prompt of poor features.

You must double-click the icon of setupDD.bat (setupDD2.bat) before using DD.bat (DD2.bat) for the first time.  Then, the absolute path for rietan.exe is recorded in %HOMEDRIVE%%HOMEPATH%Settings_RIETAN\abs_path_rietan on the use of setupDD.bat, and in %USERPROFILE%\Settings_RIETAN\abs_path_rietan on the use of setupDD2.bat. Here, %HOMEDRIVE% (usually 'C:') is the drive of the home folder, and %HOMEPATH% (usually '\') is the home folder, and %USERPROFILE% is the path of a folder (under %HOMEDRIVE%\Documents and Settings) where the profile of a logon user is stored.  The content of abs_path_rietan will later be input by DD.bat (DD2.bat).

In Windows, it is essential to associate file extensions (three characters after a period in each file name) with applications.  All the file extensions other than 'pat' for input/output files read/created by RIETAN-2000 have to be associated with your editor.

However, file extension 'ins' is associated with 'Internet Communication Settings' by default.  Then, delete this file association and associate 'ins' with your editor, which can be carried out at 'File Types' in 'Folder Options...' under the 'Tools' menu in the case of Windows 2000.  File associations are also possible by combining commands ASSOC and FTYPE in the Command Prompt (type "HELP ASSOC" and "HELP FTYPE" in the Command Prompt).  Of course, the original association should be resumed if you encounter some troubles on the use of the Internet.

In addition, file extension 'pat' should be associated with Igor Pro when using it to browse and print out Rietveld-refinement and simulated patterns with it.  In the case of Windows NT, select text/plain for the content type (MIME) further.

(c) ORFFE.bat
Skip this section unless RIETAN-2000.bat is preferred to DD.bat/DD2.bat.

Copy ORFFE.bat to anywhere as you like.  Rename ORFFE.bat appropriately, e.g., sample_FFE.bat (this name will hereafter be used).

Modify sample_FFE.bat using a text editor.  This template batch file was prepared for FapatiteE.*.

You must change the following three environment variables:

(a) LOCXYZ: (absolute) path for the location of sample.xyz (current folder: '.').
(b) SAMPLE: Sample name (in the present case 'sample').
(c) ORFFE: Absolute path for the location of orffe.exe (rietan.exe and tee.exe)

ORFFE.bat is a mere template file; you may modify it as you like.  For example, a batch file for the Command Prompt window is also fine for those who like a character user interface.

(d) lst2cif.bat
Copy lst2cif.bat to anywhere as you like.  Rename lst2cif.bat appropriately, e.g., sample_lst2cif.bat (this name will hereafter be used).

Modify sample_lst2cif.bat using a text editor.

You must change the following three environment variables:

(a) LOCLST: (absolute) path for the location of sample.lst (current folder: '.').
(b) SAMPLE: Sample name (in the present case 'sample').
(c) LST2CIF: Absolute path for the location of lst2cif.exe (rietan.exe, orffe.exe, and tee.exe).

lst2cif.bat is a mere template file; you may modify it as you like.

2. How to execute RIETAN-2000

(a) Windows 95/98/Me

Double-click the icon of sample.bat to run RIETAN-2000 (rietan.exe).

After the execution of RIETAN-2000, a DOS Prompt window is opened, and the output from RIETAN-2000 is displayed.  The same content is stored in sample.lst thanks to tee.exe.  Then, Igor Pro will be launched automatically to open sample.pat provided that sample.pat exists in the current folder.  You must quit Igor Pro before quiting the output window created by RIETAN-2000.

To quit the output window, click 'OK' in a VBScirpt window and then the X button at the top right-hand corner in the DOS Prompt window.  Use a text editor to browse sample.lst later.

Refer to comment lines in sample.bat for further details.

(b) Windows NT/2000/XP

Drag and drop *.ins onto DD.bat (Windows NT) or DD2.bat (Windows 2000/XP).

On the use of DD.bat, the location of *.ins is stored in %HOMEDRIVE%%HOMEPATH%Settings_RIETAN\abs_path_ins, and its base name in %HOMEDRIVE%%HOMEPATH%Settings_RIETAN\base_name_ins.

On the use of DD2.bat, the location of *.ins is stored in %USERPROFILE%\Settings_RIETAN\abs_path_ins, and its base name in %USERPROFILE%\Settings_RIETAN\base_name_ins.

Then, rietan.exe is lauched to read in *.ins.

If you double-click DD.bat (DD2.bat) in the absence of the above two files, an open dialog-box prompts the specification of *.ins.

If you drag and drop another *.ins file onto DD.bat (DD2.bat), rietan.exe is launched to read in it, and its location and basename in the above two files are updated.

Just after launching RIETAN-2000, a Command Prompt window, which works as a kind of a console, is opened immediately; the tee command makes it possible to display the output from RIETAN-2000 simultaneously while saving the same content in sample.lst.

After the calculation, ORFFE (orffe.exe) calculates geometric structure parameters (if *.xyz is new), Igor Pro (NPAT = 5 and *.plt does not exist) or gnuplot (NPAT = 4 and *.plt exists) displays the content of *.pat (if it is new), and *.ins is opened for the next refinement.  Finally, the output from RIETAN-2000 (*.lst) plus ORFFE (*.dst, if any) is opened by the GNU command, less.exe.  With less.exe, you can use various key operations, e.g.,

        Downward arrow --> Forward one line.
          Upward arrow --> Backward one line.
f, SPACE, or PAGE DOWN --> Scroll forward one window.
            b, PAGE UP --> Scroll backward one window.
                  HOME --> Move the cursor to the beginning of the line.
                   END --> Move the cursor to the end of the line.
              /pattern --> Search forward in the file containing the pattern.
                    :p --> Move to the previous file, i.e., *.lst
                           (used only when both *.lst and *.dst are being browsed).
                    :n --> Move to the next file, i.e., *.dst
                           (used only when both *.lst and *.dst are being browsed).
                     q --> Exit.

Clicking the left button of the mouse within the part displaying the text makes it impossible to operate less.exe further (use the left button to select the window for less.exe).  If you have happened to click the left button, click the right button, which will undo your mistake.

For details in less.exe, press 'h' (help) while browsing the output using 'less'.  You will be surprised to see a number of commands that can be used in less.exe.

Of course, extensions pat and ins must be related to a scientific graph program (e.g., Igor Pro) and a text editor you use, respectively.  If *.ins has been updated, previous *.ins file is stored as *.bak.

After the above operations, to double-click DD.bat (DD2.bat) will launch rietan.exe to read in registered *.ins and *.int.

In an alternative method of launching rietan.exe with DD.bat or DD2.bat, 'DD (base name)' or 'DD2 (base name)' is input under a directory where *.ins and *.int are located in the Command Prompt. For example, a line
DD sample
or
DD2 sample
is entered to deal with a pair of sample.ins and sample.int.  Of course, DD.bat (DD2.bat) has to be placed in a folder whose path has been set up.  In this case, the base name (sample) should contain no space in it.

Before using DD.bat (DD2.bat) again, quit the output window created by RIETAN-2000 and close *.ins in the text editor.  You need not quit Igor Pro or close the current graph window (waves for older graph windows will be overwritten).

To quit the output window, click 'OK' in a VBScirpt window.

(c) Files created by RIETAN-2000

In addition to sample.lst, sample.pat (Igor Pro file), sample.hklifile for fousyn), sample.xyz (ORFFE file), sample.mem (file for MEM analysis), sample.ffe (file for imposing restraints), and sample.ffo (file to input integrated intensities for Le Bail refinement) are created according to values of related flags.

A file named sample.ffi, which stores initial integrated intensities for Le Bail refinement, is obtained by renaming sample.ffo as much as is desired.  The file named sample.fba is created by MEED or PRIMA for REMEDY cycles.

The use of a scientific graphing software, Igor Pro, is strongly recommended to display and print out sample.pat (NPAT=5).  Igor Pro is commercial software produced by WaveMetrics, Inc.; refer to http://www.wavemetrics.com/ to obtain detailed information about it.  Its Japanese version (more expensive than its English version) is available from HULINKS (http://www.hulinks.co.jp/).

Please check 'Make a table for new experiments' in 'Experiment Settings' under menu Misc before using Igor Pro.  Other options (NPAT = 2 and 4) may work without problem, but they are no longer supported.

(d) On UNIX/GNU tools

By the way, I love my idea of running rietan.exe and tee.exe simultaneously.  How do you feel on it?  In a similar manner, I have been using other UNIX/GNU utility commands in the Command Prompt window, i.e., cat, cp, grep, bzip2, head, less, ls, mv, pwd, ps, pwd, rm, tail, and tar.  They can be downloaded from http://unxutils.sourceforge.net/ and/or http://athani.pair.com/msteed/software/unix_on_nt/.

There are so many cases where such utilities for the CUI are more convenient than operations in the GUI of Windows.

3. How to execute ORFFE

Double-click the icon of sample_FFE.bat.  A Command Prompt window is opened immediately; the tee command makes it possible to display the output from ORFFE simultaneously while saving the same content in sample.dst. After the execution of ORFFE and browsing its output, press any key to quit from the Command Prompt window.  Use a text editor for browsing sample.dst later.

Refer to comment lines in ORFFE.bat for further details.

In addition to sample.dst, sample.ffe for imposing restraints upon interatomic distances and bond angles is created if NDA > 0 and sample.ffe does not exist in the current folder.

When no instructions 2 are included in sample.xyz, ORFFE generates a series of instructions 2 to calculate A-B-C bond angles (B: atoms in the asymmetric unit) on the basis of all the B-A and B-C interatomic distances that have been evaluated with instructions 201.  These instructions 2 are appended to the tail of *.xyz.

Running ORFFE again to input the resulting *.xyz with ORFFE yields the bond angles.  In this case, *xyz remains unchanged.

4. How to execute lst2cif

This program developed by Ruben A. Dilanian during his stay at the National Institute for Research in Inorganic Materials is used to convert an output file of RIETAN-2000 (sample.lst) into a CIF (Crystallographic Information File).  For details in the CIF, refer to http://www.iucr.org/iucr-top/cif/.

Nearly all structure-drawing programs, e.g., VENUS (http://homepage.mac.com/fujioizumi/visualization/VENUS.html), ATOMS (http://www.shapesoftware.com/), Diamond (http://www.crystalimpact.com/diamond/index.html), Crystallographica (http://www.oxfordcryosystems.co.uk/software/), Balls & Sticks (http://www.ccp14.ac.uk/ccp/web-mirrors/toycrate/), and DRAWxtl (http://www.lwfinger.net/drawxtl/) in the case of Windows, can read in CIFs.  I highly recommend VENUS and Balls & Sticks because they are available free of charge and fully utilizes the OpenGL technology, the premier environment for developing portable, interactive 2D and 3D graphics applications.

Double-click the batch file for lst2cif to launch it.  Then, lst2cif convertes sample.lst (output file of RIETAN-2000) to a CIF, i.e., sample.cif.

Note that part of structure-drawing programs can import only CIFs for single phases.  In such cases, please delete lines for other phases.

Connections with some powerful programs such as ORFFE, Igor Pro, meed/fousyn/mevius/dsurf/riplx, and EXPO are important advantages of RIETAN-2000 over other Rietveld-analysis programs.  With lst2cif, RIETAN-2000 and all of the above structure-drawing programs can be associated in constructing, modifying, and visualizing structural models.  Structures under investigation can be explored through rotation, transformation, and manipulation.  It is convenient particularly when dealing with structures containing many atoms in their asymmetric units.  Have fun!

5. How to use the reflection macro

This pretty macro for Igor Pro displays hkl, lattice-plane spacing, and 2-theta for a selected reflection.

After plotting Rietveld-refinement patterns or a simulation pattern, press Ctrl-I (or Graph --> Show Info).  Drag and drop cursor A (Never use cursor B!) onto the tick mark (between observed/calculated patterns and a difference pattern) of a reflection and then select Macro --> Reflection.  Information about the specified reflection will appear in the history area as follows:

hkl = 1 1 1,  d = 3.8728 angstroms,  2-theta =  22.893 degrees

6. How to plot Rietveld-refinement patterns with gnuplot

Gnuplot is a command-driven interactive plotting program for all the major operating systems and graphic file formats.  Information about gnuplot is available from its official WWW site: http://www.gnuplot.info/.  The manual of gnuplot 3.7 (compressed as gnuplot.pdf.gz) can be downloaded from ftp://ftp.gnuplot.info/pub/gnuplot/.  A detailed primer is presented at http://t16web.lanl.gov/Kawano/gnuplot/index-e.html (in English) and http://t16web.lanl.gov/Kawano/gnuplot/ (in Japanese).

A Windows version of gnuplot 3.7.1 (compressed as gp371w32.zip) is distributed free of charge at ftp://ftp.gnuplot.info/pub/gnuplot/.  I highly recommend the enhanced version 3.7.1+1.2.0 included in the 'programs' folder.

Procedures when using gnuplot Ver. 3.7.1+1.2.0 (English menu) and one of the gnuplot script file, FapatiteE.plt, will be described below to instruct how to manipulate it.  A series of commands to plot Rietveld-refinement patterns is recorded in this file.  In what follows, let's assume that RIETAN-2000 is executed as a stand-alone program according to the procedure described in Section 2(a).

FapatiteE.plt should be transferred from the gnuplot scripts folder to the same folder as FapatiteE.ins and FapatiteE.int when NPAT is set at 4.

Prior to plotting Rietveld-refinement patterns, you have to change various settings (xmin, xmax, ymin, ymax, offset_delta, len_bar, y_phase1, y_phase2, etc.) in FapatiteE.plt.  To give proper settings, it is helpful to see only observed and calculated intensities.

Opening a text file, FapatiteE.plt, using a text editor, you will find two file names in this file: FapatiteE.eps and FapatiteE.pat.  The Encapsulated PostScript (EPS) file, FapatiteE.eps, is output when changing the terminal to 'postscript', and FapatiteE.pat is input to plot Rietveld-refinement patterns.  Of course, these file names must be changed when dealing with other *.pat files.

You will also find a line "iout = 1" (display observed, calculated, and difference patterns on the screen) in FapatiteE.plt.  Change this into "iout = 0" (display observed and calculated intensities on the screen) and save FapatiteE.plt (usually by pressing Ctrl+S).  Never forget to save it since interapplication communication is not possible between gnuplot and the text editor, in constrast with gnuplot for Mac OS.

Launch RIETAN-2000 with RIETAN-2000.bat.  After Rietveld analysis of fluorapatite with FapatiteE.ins and FapatiteE.int, RIETAN-2000 creates FapatiteE.pat for plotting Rietveld-refinement patterns with gnuplot when NPAT = 4.  Then, gnuplot is run to input in turn FapatiteE.plt and FapatiteE.pat.  A gnuplot graph window displaying observed and calculated intensities will appear immediately.  The working directory in gnuplot is now changed to a directory containing FapatiteE.plt, which can be confirmed by entering "pwd" in the gnuplot window.  Browsing this graph, you can easily determine appropriate setting values in FapatiteE.plt.

FapatiteE.plt can be opened (loaded) in a gnuplot window by either by pressing the Open button or by selecting Open ... in the File menu.

In FapatiteE.plt, characters after '#' are regarded as comments in the same manner as in *.ins.  Now, iout should be set at 1.  Be careful not to input a space after '\' indicating the continuity of a line.

After changing various settings in FapatiteE.plt appropriately with the text editor, save it.  Reopen it in a gnuplot window by either by pressing the Open button or by selecting Open ... in the File menu.  Then, observed, calculated, and difference patterns for fluorapatite are are immediately plotted together with short vertical bars indicating peak positions in the gnuplot graph window.

In this way, opening FapatiteE.plt by gnuplot after editing it by the text editor yields a new graph reflecting the changes in settings.  Modifying the script file with the editor and plotting graphs gnuplot are repeated until a satisfactory graph is obtained.

An alternative way of replotting FapatiteE.plt is as follows.  Press the up arrow key or the Prev button until a command "load 'FapatiteE.plt'" is displayed after the prompt "gnuplot> " (UNIX-like history feature).  Then, pressing the Enter key to input this command yields a new graph drawn on the basis of the new settings.

Any commands of gnuplot can be entered in the gnuplot window; for example,

pwd: print the name of the working directory.
cd 'directory_name': change the working directory.  Note that directory_name has to be enclosed by the pair of the single quotation marks.
load 'file_name': Open a gnuplot script file to execute a series of commands contained in it.
save <option> 'file_name': save a file with <option>.
show all: show all the current settings.

In the case of Cu K-alpha1 radiation, a red line indicating a 2-theta of 50 degrees can be plotted with 2-theta and d values printed when you copy and paste the following line just after 'gnuplot> ' in the Console window:

x = 50; set noarrow; set arrow from x,ymin to x,ymax nohead; print "x = ",x,", d = ",1.540562*0.5/sin(0.0087266463*x); replot

In this way, two or more commands are allowed to be input in one line by placing ';' between two commands.  Another line indicating a d spacing of 2.25 angstroms is plotted by entering the following line:

d = 2.25; x = 114.591559*asin(1.540562*0.5/d); set noarrow; set arrow from x,ymin to x,ymax nohead; print "x = ",x,", d = ",d; replot

The graphic device is initially set at 'windows' (iout = 0 and 1) on the use of Windows 95/98/Me/NT/2000/XP, which is the case in FapatiteE.plt.  If iout is set at 2, an EPS file named FapatiteE.eps is created by opening FapatiteE.plt.

The resultant EPS file, FapatiteE.eps, can be imported by graphic software like Adobe Illustrator and CorelDRAW to add characters, lines, arrows, etc.  Further, GSview (http://www.cs.wisc.edu/~ghost/gsview/index.htm) enables you to browse/print it.

Various other graphic divices are supported, e.g., hpgl, latex, tgif, pict, and x11.

After changing various settings for the Rietveld-refinement patterns, launching RIETAN-2000 through DD.bat (DD2.bat) is more convenient on Windows NT/2000/XP; with DD.bat (DD2.bat), FapatiteE.plt can be opened automatically by gnuplot.  If FapatiteE.plt exists in the same folder as *.ins and *.int on the use of DD.bat (DD2.bat), DD.bat (DD2.bat) judges that gnuplot will read in FapatiteE.pat to output Rietveld-refinement patterns.  Hence, FapatiteE.plt should never be included in that folder if NPAT is set at 5  to use Igor Pro.

Line styles (style, color, and width), fonts, background, etc. in the gnuplot graph window can be changed by clicking the right button of the mouse in this window and selecting items you like to change.  Finally, the last item, 'Update wgnuplot.init', should be selected.  Note that the current size of the gnuplot graph window is also recorded in this case.  You can alternatively change wgnuplot.ini using a text editor.

Refer to the manual of gnuplot to learn a number of commands in gnuplot. A Help menu can also be utilized.

7. On the VICS file

VICS is a program for VIsualization of Crystal Structures and contained in the VENUS package.

VICS can import user input files, *.ins, for RIETAN-2000 as well as CIFs, *.cif, created by lst2cif (File --> Import --> CIF or Command-2).

With the Save Dialog Box, VICS outputs a VICS text file, *.vcs, where crystal data and a variety of settings for drawing crystal structures are recorded.

Just after Rietveld analysis, RIETAN-2000 for Windows automatically updates lattice and structure parameters in *.vcs on the basis of refinement results.  The extension of this file needs to be 'vcs'.  For example, if the name of the user input file is sample.ins, the corresponding VICS text file should be sample.vcs.

If the resulting *.vcs is opened by VICS with the Open Dialog Box, a crystal structure is displayed from the refined lattice/structure parameters and the various settings recorded in *.vcs. Then, you can rotate/magnify/reduce objects, display interatomic distances, bond angles, and torsion angles, delete unnecessary atoms, and so forth, as you like.  To drag and drop *.vcs onto the icon of VICS.EXE (executable binary file for VICS) is an alernative to inptting *.vcs.

8. Request

Please send bug reports, comments, and constructive suggestions to the author:

Fujio Izumi
Advanced Materials Laboratory
National Institute for Materials Science
1-1 Namiki, Tsukuba, Ibaraki 305-0044, Japan
E-mail: IZUMI.Fujio@nims.go.jp
Web site: http://homepage.mac.com/fujioizumi/index.html
