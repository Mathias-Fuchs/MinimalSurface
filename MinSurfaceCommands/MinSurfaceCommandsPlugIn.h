
#pragma once

class CMinSurfaceCommandsPlugIn : public CRhinoUtilityPlugIn
{
public:
  CMinSurfaceCommandsPlugIn();
  ~CMinSurfaceCommandsPlugIn() = default;

  // Required overrides
  const wchar_t* PlugInName() const override;
  const wchar_t* PlugInVersion() const override;
  GUID PlugInID() const override;
  BOOL OnLoadPlugIn() override;
  void OnUnloadPlugIn() override;

private:
  ON_wString m_plugin_version;
};
CMinSurfaceCommandsPlugIn& MinSurfaceCommandsPlugIn();



