// MinSurfaceCommands.cpp : Defines the initialization routines for the DLL.
//

#include "stdafx.h"
#include "MinSurfaceCommandsApp.h"

//
//	Note!
//
//    A Rhino plug-in is an MFC DLL.
//
//		If this DLL is dynamically linked against the MFC
//		DLLs, any functions exported from this DLL which
//		call into MFC must have the AFX_MANAGE_STATE macro
//		added at the very beginning of the function.
//
//		For example:
//
//		extern "C" BOOL PASCAL EXPORT ExportedFunction()
//		{
//			AFX_MANAGE_STATE(AfxGetStaticModuleState());
//			// normal function body here
//		}
//
//		It is very important that this macro appear in each
//		function, prior to any calls into MFC.  This means that
//		it must appear as the first statement within the 
//		function, even before any object variable declarations
//		as their constructors may generate calls into the MFC
//		DLL.
//
//		Please see MFC Technical Notes 33 and 58 for additional
//		details.
//

// CMinSurfaceCommandsApp

BEGIN_MESSAGE_MAP(CMinSurfaceCommandsApp, CWinApp)
END_MESSAGE_MAP()

// The one and only CMinSurfaceCommandsApp object
CMinSurfaceCommandsApp theApp;

// CMinSurfaceCommandsApp initialization

BOOL CMinSurfaceCommandsApp::InitInstance()
{
  // CRITICAL: DO NOT CALL RHINO SDK FUNCTIONS HERE!
  // Only standard MFC DLL instance initialization belongs here. 
  // All other significant initialization should take place in
  // CMinSurfaceCommandsPlugIn::OnLoadPlugIn().

	CWinApp::InitInstance();

	// https://gist.github.com/hi2p-perim/7855506
	int cpuinfo[4];
	__cpuid(cpuinfo, 1);
	bool avxSupportted = cpuinfo[2] & (1 << 28) || false;
	bool osxsaveSupported = cpuinfo[2] & (1 << 27) || false;
	if (osxsaveSupported && avxSupportted)
	{
		// _XCR_XFEATURE_ENABLED_MASK = 0
		unsigned long long xcrFeatureMask = _xgetbv(0);
		avxSupportted = (xcrFeatureMask & 0x6) == 0x6;
	}

	if (!avxSupportted) {
		RhinoApp().Print(L"Can't load the MinSurface command line interface since this processor does not seem to support AVX.");
		return FALSE;
	}
	return TRUE;
}

int CMinSurfaceCommandsApp::ExitInstance()
{
  // CRITICAL: DO NOT CALL RHINO SDK FUNCTIONS HERE!
  // Only standard MFC DLL instance clean up belongs here. 
  // All other significant cleanup should take place in either
  // CMinSurfaceCommandsPlugIn::OnSaveAllSettings() or
  // CMinSurfaceCommandsPlugIn::OnUnloadPlugIn().
  return CWinApp::ExitInstance();
}
