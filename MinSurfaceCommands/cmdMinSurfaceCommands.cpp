#include "StdAfx.h"
#include "LaplaceData.h"
#include "MinSurfaceCommandsPlugIn.h"

// int n - the root of the number of desired vertices
ON_Mesh CreateDiskMesh(int n) {
	if (n < 2) return ON_Mesh();
	ON_Mesh m;
	int quadcount = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			double x = 2 * (double)i / (n - 1) - 1;
			double y = 2 * (double)j / (n - 1) - 1;

			double u = x * sqrt(1 - y * y / 2.0);
			double v = y * sqrt(1 - x * x / 2.0);

			m.SetVertex(j + n * i, ON_3dPoint(u, v, 0));
			if (i == 0 || j == 0) continue;
			m.SetQuad(quadcount++,
				j - 1 + n * (i - 1),
				j - 0 + n * (i - 1),
				j - 0 + n * (i - 0),
				j - 1 + n * (i - 0));
		}
	return m;
}




// Do NOT put the definition of class CCommandMinSurfaceCommands in a header
// file. There is only ONE instance of a CCommandMinSurfaceCommands class
// and that instance is the static theMinSurfaceCommandsCommand that appears
// immediately below the class definition.

class MiniOneCommand : public CRhinoCommand
{
public:
	// The one and only instance of CCommandMinSurfaceCommands is created below.
	// No copy constructor or operator= is required.
	// Values of member variables persist for the duration of the application.

	// CCommandMinSurfaceCommands::CCommandMinSurfaceCommands()
	// is called exactly once when static theMinSurfaceCommandsCommand is created.
	MiniOneCommand() = default;

	// CCommandMinSurfaceCommands::~CCommandMinSurfaceCommands()
	// is called exactly once when static theMinSurfaceCommandsCommand is destroyed.
	// The destructor should not make any calls to the Rhino SDK. 
	// If your command has persistent settings, then override 
	// CRhinoCommand::SaveProfile and CRhinoCommand::LoadProfile.
	~MiniOneCommand() = default;

	// Returns a unique UUID for this command.
	// If you try to use an id that is already being used, then
	// your command will not work. Use GUIDGEN.EXE to make unique UUID.
	UUID CommandUUID() override
	{
		// {F1C2449E-7E92-4DE0-A5F1-5A1281E81D62}
		static const GUID MinSurfaceCommandsCommand_UUID =
		{ 0xF1C2449E, 0x7E92, 0x4DE0, { 0xA5, 0xF1, 0x5A, 0x12, 0x81, 0xE8, 0x1D, 0x62 } };
		return MinSurfaceCommandsCommand_UUID;
	}

	// Returns the English command name.
	// If you want to provide a localized command name, then override 
	// CRhinoCommand::LocalCommandName.
	const wchar_t* EnglishCommandName() override { return L"MinSurfFromCircularBoundary"; }

	// Rhino calls RunCommand to run the command.
	CRhinoCommand::result RunCommand(const CRhinoCommandContext& context) override;
};

// The one and only CCommandMinSurfaceCommands object
// Do NOT create any other instance of a CCommandMinSurfaceCommands class.
static class MiniOneCommand theMinSurfaceCommandsCommand;

CRhinoCommand::result MiniOneCommand::RunCommand(const CRhinoCommandContext& context) {

	if (!context.IsInteractive()) {
		RhinoApp().Print(L"The minimal surface command only works in interactive context.");
		return CRhinoCommand::failure;
	}

	CRhinoGetObject go;
	int degree = 100;
	int nrBoundaryVertices = 300;
	go.AddCommandOptionInteger(RHCMDOPTNAME(L"Degree"), &degree, L"Degree of output mesh as surface", 0, 10000);
	go.AddCommandOptionInteger(RHCMDOPTNAME(L"MeshDensity"), &nrBoundaryVertices, L"Number of vertices along the boundary of the output mesh", 100, 100000);
	go.SetCommandPrompt(L"Select a closed curve for the minimal surface boundary.");
	go.SetGeometryFilter(CRhinoGetObject::curve_object);
	go.SetGeometryAttributeFilter(CRhinoGetObject::closed_curve);
	while (true) {
		CRhinoGet::result res = go.GetObjects(1, 1);
		if (res == CRhinoGet::option) continue;
		if (res != CRhinoGet::object) return CRhinoCommand::cancel;
		break;
	}
	const ON_Curve* targetCurve = go.Object(0).Curve();
	if (!targetCurve) return CRhinoCommand::failure;

	//// number of control points, tells about the complexity of the curve
	//int nrCont = _targetCurve->NurbsCurve()->CVCount();
	//int crDeg = _targetCurve->Degree();

	//  number of boundary subdivisions for computation of the polynomials
	// one needs to fit about 2*degree real parameters, each subdivision gives one real observation
	int n = 100 * degree;
	double t0, t1;
	targetCurve->GetDomain(&t0, &t1);

	std::vector<double> targetPointsx, targetPointsy, targetPointsz;

	targetPointsx.reserve(n);
	targetPointsy.reserve(n);
	targetPointsz.reserve(n);

	std::vector<double> s(n);
	std::vector<double> t(n);

	for (int i = 0; i < n; i++) s[i] = (double) i / n;
	
	targetCurve->GetNormalizedArcLengthPoints(n, s.data(), t.data());

	for (int i = 0; i < n; i++) {
		auto p = targetCurve->PointAt(t[i]);
		targetPointsx.push_back(p.x);
		targetPointsy.push_back(p.y);
		targetPointsz.push_back(p.z);
	}

	LaplaceData kx(targetPointsx, degree);
	LaplaceData ky(targetPointsy, degree);
	LaplaceData kz(targetPointsz, degree);

	ON_Mesh disk = CreateDiskMesh(nrBoundaryVertices / 4);

	int nv = disk.VertexCount();
	for (int i = 0; i < nv; i++) {
		const ON_3dPoint& p = disk.Vertex(i);
		auto goal = ON_3dPoint(
			kx.eval(p.x, p.y),
			ky.eval(p.x, p.y),
			kz.eval(p.x, p.y)
		);
		disk.SetVertex(i, goal);
	}

	RhinoApp().ActiveDoc()->AddMeshObject(disk);
	return CRhinoCommand::success;
}
