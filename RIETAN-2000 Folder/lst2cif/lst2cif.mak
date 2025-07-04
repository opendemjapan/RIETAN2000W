# Microsoft Developer Studio Generated NMAKE File, Based on lst2cif.dsp
!IF "$(CFG)" == ""
CFG=lst2cif - Win32 Debug
!MESSAGE No configuration specified. Defaulting to lst2cif - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "lst2cif - Win32 Release" && "$(CFG)" != "lst2cif - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "lst2cif.mak" CFG="lst2cif - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "lst2cif - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "lst2cif - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "lst2cif - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\lst2cif.exe"


CLEAN :
	-@erase "$(INTDIR)\block0.obj"
	-@erase "$(INTDIR)\block1.obj"
	-@erase "$(INTDIR)\block2.obj"
	-@erase "$(INTDIR)\block3.obj"
	-@erase "$(INTDIR)\block4.obj"
	-@erase "$(INTDIR)\block5.obj"
	-@erase "$(INTDIR)\date_time.obj"
	-@erase "$(INTDIR)\readlst.obj"
	-@erase "$(INTDIR)\tools.obj"
	-@erase "$(INTDIR)\windos.obj"
	-@erase "$(INTDIR)\writecif.obj"
	-@erase "$(OUTDIR)\lst2cif.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90=df.exe
F90_PROJ=/alignment:dcommons /alignment:sequence /architecture:p5 /compile_only /fixed /math_library:fast /nologo /optimize:5 /tune:p5 /warn:nofileopt /module:"Release/" /object:"Release/" 
F90_OBJS=.\Release/

.SUFFIXES: .fpp

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.fpp{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

CPP=cl.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /Fp"$(INTDIR)\lst2cif.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

RSC=rc.exe
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\lst2cif.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\lst2cif.pdb" /machine:I386 /out:"$(OUTDIR)\lst2cif.exe" 
LINK32_OBJS= \
	"$(INTDIR)\block0.obj" \
	"$(INTDIR)\block1.obj" \
	"$(INTDIR)\block2.obj" \
	"$(INTDIR)\block3.obj" \
	"$(INTDIR)\block4.obj" \
	"$(INTDIR)\block5.obj" \
	"$(INTDIR)\date_time.obj" \
	"$(INTDIR)\readlst.obj" \
	"$(INTDIR)\tools.obj" \
	"$(INTDIR)\windos.obj" \
	"$(INTDIR)\writecif.obj"

"$(OUTDIR)\lst2cif.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "lst2cif - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

ALL : "$(OUTDIR)\lst2cif.exe"


CLEAN :
	-@erase "$(INTDIR)\block0.obj"
	-@erase "$(INTDIR)\block1.obj"
	-@erase "$(INTDIR)\block2.obj"
	-@erase "$(INTDIR)\block3.obj"
	-@erase "$(INTDIR)\block4.obj"
	-@erase "$(INTDIR)\block5.obj"
	-@erase "$(INTDIR)\date_time.obj"
	-@erase "$(INTDIR)\DF60.PDB"
	-@erase "$(INTDIR)\readlst.obj"
	-@erase "$(INTDIR)\tools.obj"
	-@erase "$(INTDIR)\windos.obj"
	-@erase "$(INTDIR)\writecif.obj"
	-@erase "$(OUTDIR)\lst2cif.exe"
	-@erase "$(OUTDIR)\lst2cif.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90=df.exe
F90_PROJ=/check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt /module:"Debug/" /object:"Debug/" /pdbfile:"Debug/DF60.PDB" 
F90_OBJS=.\Debug/

.SUFFIXES: .fpp

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.fpp{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /Fp"$(INTDIR)\lst2cif.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

RSC=rc.exe
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\lst2cif.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\lst2cif.pdb" /debug /machine:I386 /out:"$(OUTDIR)\lst2cif.exe" /pdbtype:sept 
LINK32_OBJS= \
	"$(INTDIR)\block0.obj" \
	"$(INTDIR)\block1.obj" \
	"$(INTDIR)\block2.obj" \
	"$(INTDIR)\block3.obj" \
	"$(INTDIR)\block4.obj" \
	"$(INTDIR)\block5.obj" \
	"$(INTDIR)\date_time.obj" \
	"$(INTDIR)\readlst.obj" \
	"$(INTDIR)\tools.obj" \
	"$(INTDIR)\windos.obj" \
	"$(INTDIR)\writecif.obj"

"$(OUTDIR)\lst2cif.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("lst2cif.dep")
!INCLUDE "lst2cif.dep"
!ELSE 
!MESSAGE Warning: cannot find "lst2cif.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "lst2cif - Win32 Release" || "$(CFG)" == "lst2cif - Win32 Debug"
SOURCE=.\Source\block0.f90

"$(INTDIR)\block0.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\Source\block1.f90

"$(INTDIR)\block1.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\Source\block2.f90

"$(INTDIR)\block2.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\Source\block3.f90

"$(INTDIR)\block3.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\Source\block4.f90

"$(INTDIR)\block4.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\Source\block5.f90

"$(INTDIR)\block5.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\Source\date_time.f90

"$(INTDIR)\date_time.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\Source\readlst.f90

"$(INTDIR)\readlst.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\Source\tools.f90

"$(INTDIR)\tools.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\Source\windos.f90

"$(INTDIR)\windos.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\Source\writecif.f90

"$(INTDIR)\writecif.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)



!ENDIF 

