# Microsoft Developer Studio Project File - Name="isart" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) QuickWin Application" 0x0107

CFG=isart - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "isart.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "isart.mak" CFG="isart - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "isart - Win32 Release" (based on "Win32 (x86) QuickWin Application")
!MESSAGE "isart - Win32 Debug" (based on "Win32 (x86) QuickWin Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "isart - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "c:\bin"
# PROP Intermediate_Dir "c:\tmp"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /include:"Release/" /libs:qwin /nologo /warn:nofileopt
# ADD F90 /alignment:dcommons /compile_only /fpscomp:nolibs /include:"c:\tmp/" /libs:qwin /nologo /warn:nofileopt /warn:unused
# SUBTRACT F90 /fpscomp:general /fpscomp:ioformat
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /entry:"WinMainCRTStartup" /subsystem:windows /machine:I386 /nodefaultlib:"dfconsol.lib"
# ADD LINK32 kernel32.lib /nologo /stack:0x3d0900 /entry:"WinMainCRTStartup" /subsystem:windows /incremental:yes /machine:I386
# SUBTRACT LINK32 /pdb:none

!ELSEIF  "$(CFG)" == "isart - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "c:\bin\debug"
# PROP Intermediate_Dir "c:\tmp\debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /include:"Debug/" /libs:qwin /nologo /warn:argument_checking /warn:nofileopt
# ADD F90 /alignment:dcommons /check:bounds /compile_only /debug:full /fpscomp:nolibs /include:"c:\tmp\debug/" /libs:qwin /nologo /warn:argument_checking /warn:nofileopt /warn:unused
# SUBTRACT F90 /fpscomp:general /fpscomp:ioformat
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /entry:"WinMainCRTStartup" /subsystem:windows /debug /machine:I386 /nodefaultlib:"dfconsol.lib" /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /entry:"WinMainCRTStartup" /subsystem:windows /incremental:no /debug /machine:I386 /nodefaultlib:"dfconsol.lib" /pdbtype:sept

!ENDIF 

# Begin Target

# Name "isart - Win32 Release"
# Name "isart - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=..\common\addtgt.f
DEP_F90_ADDTG=\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"tglist.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=.\adjust.f
DEP_F90_ADJUS=\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"tglist.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=..\common\cancel.f
# End Source File
# Begin Source File

SOURCE=..\common\clight_avg.f
# End Source File
# Begin Source File

SOURCE=..\common\cmlib.f
# End Source File
# Begin Source File

SOURCE=..\common\defaults.f
DEP_F90_DEFAU=\
	{$(INCLUDE)}"kalman.h"\
	{$(INCLUDE)}"realtime.h"\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"tglist.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=..\common\ellipse_param.f
DEP_F90_ELLIP=\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"tglist.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=..\common\fastfit_rt.f
# End Source File
# Begin Source File

SOURCE=..\common\fillbuff.f
DEP_F90_FILLB=\
	{$(INCLUDE)}"realtime.h"\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=..\common\fixem.f
# End Source File
# Begin Source File

SOURCE=..\common\fourt.f
# End Source File
# Begin Source File

SOURCE=..\common\getacc.f
DEP_F90_GETAC=\
	{$(INCLUDE)}"kalman.h"\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"tglist.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=..\common\getacc_rt.f
DEP_F90_GETACC=\
	{$(INCLUDE)}"realtime.h"\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"tglist.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=..\common\getraw.f
DEP_F90_GETRA=\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=.\gettgt_s.f
DEP_F90_GETTG=\
	{$(INCLUDE)}"realtime.h"\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"tglist.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=..\common\imgen.f
# End Source File
# Begin Source File

SOURCE=..\common\imgen_rt.f
DEP_F90_IMGEN=\
	{$(INCLUDE)}"realtime.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=.\imgenu.f
# End Source File
# Begin Source File

SOURCE=..\common\iqgen.f
DEP_F90_IQGEN=\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=.\isarlogo.f
# End Source File
# Begin Source File

SOURCE=.\isart.f
DEP_F90_ISART=\
	{$(INCLUDE)}"kalman.h"\
	{$(INCLUDE)}"motime.h"\
	{$(INCLUDE)}"realtime.h"\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=..\common\kalman.f
DEP_F90_KALMA=\
	{$(INCLUDE)}"kalman.h"\
	{$(INCLUDE)}"tglist.h"\
	
# End Source File
# Begin Source File

SOURCE=.\lamps.f
DEP_F90_LAMPS=\
	{$(INCLUDE)}"realtime.h"\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=.\maritime.f
DEP_F90_MARIT=\
	{$(INCLUDE)}"realtime.h"\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=..\common\mathsubs.f
# End Source File
# Begin Source File

SOURCE=.\mdiags.f
DEP_F90_MDIAG=\
	{$(INCLUDE)}"motime.h"\
	{$(INCLUDE)}"realtime.h"\
	
# End Source File
# Begin Source File

SOURCE=.\mocomp.f
DEP_F90_MOCOM=\
	{$(INCLUDE)}"sarprm.h"\
	
# End Source File
# Begin Source File

SOURCE=.\movie.f
DEP_F90_MOVIE=\
	{$(INCLUDE)}"kalman.h"\
	{$(INCLUDE)}"motime.h"\
	{$(INCLUDE)}"realtime.h"\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"tglist.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=.\movie0.f
DEP_F90_MOVIE0=\
	{$(INCLUDE)}"motime.h"\
	{$(INCLUDE)}"realtime.h"\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"tglist.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=.\norden.f
DEP_F90_NORDE=\
	{$(INCLUDE)}"realtime.h"\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=..\common\pga.f
DEP_F90_PGA_F=\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=..\common\polar.f
# End Source File
# Begin Source File

SOURCE=..\common\quant.f
# End Source File
# Begin Source File

SOURCE=.\radars.f
DEP_F90_RADAR=\
	{$(INCLUDE)}"realtime.h"\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=..\common\rawdat.f
DEP_F90_RAWDA=\
	{$(INCLUDE)}"aps137hd.h"\
	
# End Source File
# Begin Source File

SOURCE=..\common\realtlib.f
# End Source File
# Begin Source File

SOURCE=..\common\remove_pg.f
# End Source File
# Begin Source File

SOURCE=.\rfcomb0.f
# End Source File
# Begin Source File

SOURCE=..\common\sarsubs.f
# End Source File
# Begin Source File

SOURCE=..\common\scpsubs.f
# End Source File
# Begin Source File

SOURCE=..\common\shipsim.f
# End Source File
# Begin Source File

SOURCE=.\subimg.f
DEP_F90_SUBIM=\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=.\subimg0.f
DEP_F90_SUBIMG=\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=.\subimg_rt.f
DEP_F90_SUBIMG_=\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=..\common\syssun.f
# End Source File
# Begin Source File

SOURCE=.\target.f
DEP_F90_TARGE=\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=.\target_rt.f
DEP_F90_TARGET=\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=..\common\tcolor.f
DEP_F90_TCOLO=\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"tglist.h"\
	{$(INCLUDE)}"updates.h"\
	

!IF  "$(CFG)" == "isart - Win32 Release"

# ADD F90 /optimize:4

!ELSEIF  "$(CFG)" == "isart - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\common\tgtbox.f
DEP_F90_TGTBO=\
	{$(INCLUDE)}"sarprm.h"\
	{$(INCLUDE)}"tglist.h"\
	{$(INCLUDE)}"updates.h"\
	
# End Source File
# Begin Source File

SOURCE=..\common\unirnor.f
# End Source File
# Begin Source File

SOURCE=.\uwbsub.f
# End Source File
# Begin Source File

SOURCE=..\common\wininit.f
DEP_F90_WININ=\
	"..\common\scr_stuf.f"\
	{$(INCLUDE)}"rmovie.h"\
	
# End Source File
# Begin Source File

SOURCE=..\common\wtffti.f
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
