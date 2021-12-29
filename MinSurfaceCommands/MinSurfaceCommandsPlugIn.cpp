// MinSurfaceCommandsPlugIn.cpp : defines the initialization routines for the plug-in.
//

#include "StdAfx.h"
#include "rhinoSdkPlugInDeclare.h"
#include "MinSurfaceCommandsPlugIn.h"
#include "Resource.h"

// The plug-in object must be constructed before any plug-in classes derived
// from CRhinoCommand. The #pragma init_seg(lib) ensures that this happens.
#pragma warning(push)
#pragma warning(disable : 4073)
#pragma init_seg(lib)
#pragma warning(pop)

// Rhino plug-in declaration
RHINO_PLUG_IN_DECLARE

// Rhino plug-in name
// Provide a short, friendly name for this plug-in.
RHINO_PLUG_IN_NAME(L"MinSurface");

// Rhino plug-in id
// Provide a unique uuid for this plug-in.
RHINO_PLUG_IN_ID(L"219FFEA1-4ED5-45E9-96F7-42FB81EEF80A");

// Rhino plug-in version
// Provide a version number string for this plug-in.
RHINO_PLUG_IN_VERSION("0.1.3.8")

// Rhino plug-in description
// Provide a description of this plug-in.
RHINO_PLUG_IN_DESCRIPTION(L"MinSurface plug-in for Rhinoceros®");

// Rhino plug-in icon resource id
// Provide an icon resource this plug-in.
// Icon resource should contain 16, 24, 32, 48, and 256-pixel image sizes.
RHINO_PLUG_IN_ICON_RESOURCE_ID(IDI_ICON);

// Rhino plug-in developer declarations
// TODO: fill in the following developer declarations with
// your company information. Note, all of these declarations
// must be present or your plug-in will not load.
//
// When completed, delete the following #error directive.
RHINO_PLUG_IN_DEVELOPER_ORGANIZATION(L"Dr. Mathias Fuchs");
RHINO_PLUG_IN_DEVELOPER_ADDRESS(L"Hallgartenstr. 6");
RHINO_PLUG_IN_DEVELOPER_COUNTRY(L"Germany");
RHINO_PLUG_IN_DEVELOPER_PHONE(L"+49/174 584 6212");
RHINO_PLUG_IN_DEVELOPER_FAX(L"-");
RHINO_PLUG_IN_DEVELOPER_EMAIL(L"mathias@mathiasfuchs.com");
RHINO_PLUG_IN_DEVELOPER_WEBSITE(L"https://mathiasfuchs.com");
RHINO_PLUG_IN_UPDATE_URL(L"https://mathiasfuchs.com");

// The one and only CMinSurfaceCommandsPlugIn object
static class CMinSurfaceCommandsPlugIn thePlugIn;

/////////////////////////////////////////////////////////////////////////////
// CMinSurfaceCommandsPlugIn definition

CMinSurfaceCommandsPlugIn& MinSurfaceCommandsPlugIn()
{
	// Return a reference to the one and only CMinSurfaceCommandsPlugIn object
	return thePlugIn;
}

CMinSurfaceCommandsPlugIn::CMinSurfaceCommandsPlugIn()
{
	// Description:
	//   CMinSurfaceCommandsPlugIn constructor. The constructor is called when the
	//   plug-in is loaded and "thePlugIn" is constructed. Once the plug-in
	//   is loaded, CMinSurfaceCommandsPlugIn::OnLoadPlugIn() is called. The
	//   constructor should be simple and solid. Do anything that might fail in
	//   CMinSurfaceCommandsPlugIn::OnLoadPlugIn().

	// TODO: Add construction code here
	m_plugin_version = RhinoPlugInVersion();
}

/////////////////////////////////////////////////////////////////////////////
// Required overrides

const wchar_t* CMinSurfaceCommandsPlugIn::PlugInName() const
{
	// Description:
	//   Plug-in name display string.  This name is displayed by Rhino when
	//   loading the plug-in, in the plug-in help menu, and in the Rhino
	//   interface for managing plug-ins.

	// TODO: Return a short, friendly name for the plug-in.
	return RhinoPlugInName();
}

const wchar_t* CMinSurfaceCommandsPlugIn::PlugInVersion() const
{
	// Description:
	//   Plug-in version display string. This name is displayed by Rhino
	//   when loading the plug-in and in the Rhino interface for managing
	//   plug-ins.

	// TODO: Return the version number of the plug-in.
	return m_plugin_version;
}

GUID CMinSurfaceCommandsPlugIn::PlugInID() const
{
	// Description:
	//   Plug-in unique identifier. The identifier is used by Rhino to
	//   manage the plug-ins.

	// TODO: Return a unique identifier for the plug-in.
	// {219FFEA1-4ED5-45E9-96F7-42FB81EEF80A}
	return ON_UuidFromString(RhinoPlugInId());
}

/////////////////////////////////////////////////////////////////////////////
// Additional overrides

BOOL CMinSurfaceCommandsPlugIn::OnLoadPlugIn()
{
	// Description:
	//   Called after the plug-in is loaded and the constructor has been
	//   run. This is a good place to perform any significant initialization,
	//   license checking, and so on.  This function must return TRUE for
	//   the plug-in to continue to load.

	// Remarks:
	//    Plug-ins are not loaded until after Rhino is started and a default document
	//    is created.  Because the default document already exists
	//    CRhinoEventWatcher::On????Document() functions are not called for the default
	//    document.  If you need to do any document initialization/synchronization then
	//    override this function and do it here.  It is not necessary to call
	//    CPlugIn::OnLoadPlugIn() from your derived class.

	// TODO: Add plug-in initialization code here.

	return TRUE;
}

void CMinSurfaceCommandsPlugIn::OnUnloadPlugIn()
{
	// Description:
	//    Called one time when plug-in is about to be unloaded. By this time,
	//    Rhino's mainframe window has been destroyed, and some of the SDK
	//    managers have been deleted. There is also no active document or active
	//    view at this time. Thus, you should only be manipulating your own objects.
	//    or tools here.

	// TODO: Add plug-in cleanup code here.
}
