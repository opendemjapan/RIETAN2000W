# Microsoft Developer Studio Generated NMAKE File, Based on rietan.dsp
!IF "$(CFG)" == ""
CFG=rietan - Win32 Debug
!MESSAGE No configuration specified. Defaulting to rietan - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "rietan - Win32 Release" && "$(CFG)" != "rietan - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "rietan.mak" CFG="rietan - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "rietan - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "rietan - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "rietan - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\rietan.exe"


CLEAN :
	-@erase "$(INTDIR)\jobend.obj"
	-@erase "$(INTDIR)\math77lib.obj"
	-@erase "$(INTDIR)\rietan1.obj"
	-@erase "$(INTDIR)\rietan2.obj"
	-@erase "$(INTDIR)\rietan3.obj"
	-@erase "$(INTDIR)\rietan4.obj"
	-@erase "$(INTDIR)\rietan5.obj"
	-@erase "$(INTDIR)\windos.obj"
	-@erase "$(OUTDIR)\rietan.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90_PROJ=/alignment:dcommons /alignment:sequence /architecture:pn2 /compile_only /debug:none /fixed /math_library:fast /nologo /optimize:5 /tune:pn2 /warn:nofileopt /module:"Release/" /object:"Release/" 
F90_OBJS=.\Release/
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /Fp"$(INTDIR)\rietan.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\rietan.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\rietan.pdb" /machine:I386 /out:"$(OUTDIR)\rietan.exe" 
LINK32_OBJS= \
	"$(INTDIR)\jobend.obj" \
	"$(INTDIR)\math77lib.obj" \
	"$(INTDIR)\rietan1.obj" \
	"$(INTDIR)\rietan2.obj" \
	"$(INTDIR)\rietan3.obj" \
	"$(INTDIR)\rietan4.obj" \
	"$(INTDIR)\rietan5.obj" \
	"$(INTDIR)\windos.obj"

"$(OUTDIR)\rietan.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "rietan - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

ALL : "$(OUTDIR)\rietan.exe"


CLEAN :
	-@erase "$(INTDIR)\DF60.PDB"
	-@erase "$(INTDIR)\jobend.obj"
	-@erase "$(INTDIR)\math77lib.obj"
	-@erase "$(INTDIR)\rietan1.obj"
	-@erase "$(INTDIR)\rietan2.obj"
	-@erase "$(INTDIR)\rietan3.obj"
	-@erase "$(INTDIR)\rietan4.obj"
	-@erase "$(INTDIR)\rietan5.obj"
	-@erase "$(INTDIR)\windos.obj"
	-@erase "$(OUTDIR)\rietan.exe"
	-@erase "$(OUTDIR)\rietan.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90_PROJ=/check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt /module:"Debug/" /object:"Debug/" /pdbfile:"Debug/DF60.PDB" 
F90_OBJS=.\Debug/
CPP_PROJ=/nologo /MLd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /Fp"$(INTDIR)\rietan.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\rietan.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\rietan.pdb" /debug /machine:I386 /out:"$(OUTDIR)\rietan.exe" /pdbtype:sept 
LINK32_OBJS= \
	"$(INTDIR)\jobend.obj" \
	"$(INTDIR)\math77lib.obj" \
	"$(INTDIR)\rietan1.obj" \
	"$(INTDIR)\rietan2.obj" \
	"$(INTDIR)\rietan3.obj" \
	"$(INTDIR)\rietan4.obj" \
	"$(INTDIR)\rietan5.obj" \
	"$(INTDIR)\windos.obj"

"$(OUTDIR)\rietan.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

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

.SUFFIXES: .fpp

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.fpp{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("rietan.dep")
!INCLUDE "rietan.dep"
!ELSE 
!MESSAGE Warning: cannot find "rietan.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "rietan - Win32 Release" || "$(CFG)" == "rietan - Win32 Debug"
SOURCE=.\Source\jobend.f90

"$(INTDIR)\jobend.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\Source\math77lib.f90

"$(INTDIR)\math77lib.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\Source\rietan1.f90

"$(INTDIR)\rietan1.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\Source\rietan2.f90

"$(INTDIR)\rietan2.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\Source\rietan3.f90

"$(INTDIR)\rietan3.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\Source\rietan4.f90

"$(INTDIR)\rietan4.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\Source\rietan5.f90

"$(INTDIR)\rietan5.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\Source\windos.f90

"$(INTDIR)\windos.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)



!ENDIF 

