#include "StdAfx.h"
#include "LaplaceData.h"
#include "MinSurfaceCommandsPlugIn.h"

// int n - the root of the number of desired vertices
ON_Mesh CreateDiskMesh(int n) {
	if (n < 2) return ON_Mesh();
	ON_Mesh m;
	m.ReserveVertexCapacity(n * n);
	m.m_F.Reserve(n * n - 2 * n + 1);
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


ON_Mesh CreateUncappedCylinder(int axis_divisions, int height_divisions) {
	ON_Mesh M;
	if (axis_divisions <= 0 || height_divisions <= 0) return M;
	M.ReserveVertexCapacity(axis_divisions * height_divisions);
	M.m_F.Reserve(2 * (axis_divisions * (height_divisions - 1)));

	int vc = 0;
	for (int th = 0; th < axis_divisions; th++)
	{
		double x = cos(2.0 * M_PI * th / (double)axis_divisions);
		double y = sin(2.0 * M_PI * th / (double)axis_divisions);

		for (int h = 0; h < height_divisions; h++)
		{
			double z = h / (double)(height_divisions - 1);
			M .SetVertex(vc++, ON_3dPoint(x, y, z));
		}
	}

	int fc = 0;
	for (int th = 0; th < axis_divisions; th++)
		for (int h = 1; h < height_divisions; h++)
			M.SetQuad(fc++,
				((th + 0) % axis_divisions) * height_divisions + (h - 1),
				((th + 1) % axis_divisions) * height_divisions + (h - 1),
				((th + 1) % axis_divisions) * height_divisions + (h + 0),
				((th + 0) % axis_divisions) * height_divisions + (h + 0)
			);
	return M;
}



// Do NOT put the definition of class CCommandMinSurfaceCommands in a header
// file. There is only ONE instance of a CCommandMinSurfaceCommands class
// and that instance is the static theMinSurfaceCommandsCommand that appears
// immediately below the class definition.

class MiniOneCommand : public CRhinoCommand
{
public:
	MiniOneCommand() = default;
	~MiniOneCommand() = default;
	UUID CommandUUID() override	{
		static const GUID MinSurfaceCommandsCommand_UUID =
		{ 0xF1C2449E, 0x7E92, 0x4DE0, { 0xA5, 0xF1, 0x5A, 0x12, 0x81, 0xE8, 0x1D, 0x62 } };
		return MinSurfaceCommandsCommand_UUID;
	}
	const wchar_t* EnglishCommandName() override { return L"MinSurfFromCircularBoundary"; }
	CRhinoCommand::result RunCommand(const CRhinoCommandContext& context) override;
};

static class MiniOneCommand theMinSurfaceCommandsCommand;

CRhinoCommand::result MiniOneCommand::RunCommand(const CRhinoCommandContext& context) {

	if (!context.IsInteractive()) {
		RhinoApp().Print(L"The minimal surface commands only work in interactive context.");
		return CRhinoCommand::failure;
	}

	ON_wString str;
	str.Format(L"MinSurface plugin version %s \n", MinSurfaceCommandsPlugIn().PlugInVersion());
	RhinoApp().Print(str);
	RhinoApp().Print(L"GPLv3 licensed, source: https://github.com/Mathias-Fuchs/MinimalSurface \n");
	RhinoApp().Print(L"Copyright Mathias Fuchs 2020 - 2021, https://mathiasfuchs.com \n");
	RhinoApp().Print(L"This is the command line interface, please make sure to try out the Grasshopper interface as well.\n");

	CRhinoGetObject go;
	int degree = 40;
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
	context.m_doc.Redraw();


	return CRhinoCommand::success;
}






class MiniTwoCommand : public CRhinoCommand
{
public:
	MiniTwoCommand() = default;
	UUID CommandUUID() override
	{
		// {874500EC-F717-4972-B992-2E606B1FE350}
		static const GUID MinSurfFromTwoBoundariesCommand_UUID =
		{ 0x874500EC, 0xF717, 0x4972, { 0xB9, 0x92, 0x2E, 0x60, 0x6B, 0x1F, 0xE3, 0x50 } };
		return MinSurfFromTwoBoundariesCommand_UUID;
	}
	const wchar_t* EnglishCommandName() override { return L"MinSurfFromTwoBoundaries"; }
	CRhinoCommand::result RunCommand(const CRhinoCommandContext& context) override;
};

static class MiniTwoCommand theMinSurfFromTwoBoundariesCommand;

CRhinoCommand::result MiniTwoCommand::RunCommand(const CRhinoCommandContext& context)
{
	if (!context.IsInteractive()) {
		RhinoApp().Print(L"The minimal surface commands only work in interactive context.");
		return CRhinoCommand::failure;
	}

	ON_wString str;
	str.Format(L"MinSurface plugin version %s \n", MinSurfaceCommandsPlugIn().PlugInVersion());
	RhinoApp().Print(str);
	RhinoApp().Print(L"GPLv3 licensed, source: https://github.com/Mathias-Fuchs/MinimalSurface \n");
	RhinoApp().Print(L"Copyright Mathias Fuchs 2020 - 2021, https://mathiasfuchs.com \n");
	RhinoApp().Print(L"This is the command line interface, please make sure to try out the Grasshopper interface as well.\n");

	CRhinoGetObject go;
	int degree = 50;
	int nrVerticesAround = 100;
	int nrVerticesVertical = 100;
	go.AddCommandOptionInteger(RHCMDOPTNAME(L"Degree"), &degree, L"Degree of output mesh as surface", 0, 10000);
	go.AddCommandOptionInteger(RHCMDOPTNAME(L"NverticesAround"), &nrVerticesAround, L"Number of vertices in the ``horizontal'' direction", 10, 100000);
	go.AddCommandOptionInteger(RHCMDOPTNAME(L"NverticesVertical"), &nrVerticesVertical, L"Number of vertices in the ``vertical'' direction", 10, 100000);

	go.SetCommandPrompt(L"Select two closed curves for the minimal surface boundary.");
	go.SetGeometryFilter(CRhinoGetObject::curve_object);
	go.SetGeometryAttributeFilter(CRhinoGetObject::closed_curve);
	while (true) {
		CRhinoGet::result res = go.GetObjects(1, 2);
		if (res == CRhinoGet::option) continue;
		if (res != CRhinoGet::object) return CRhinoCommand::cancel;
		if (go.ObjectCount() == 2) break;
	}
	const ON_Curve* targetCurve1 = go.Object(0).Curve();
	const ON_Curve* targetCurve2 = go.Object(1).Curve();

	if (!targetCurve1 || !targetCurve2) return CRhinoCommand::failure;

	//  number of boundary subdivisions for computation of the polynomials
	// one needs to fit about 2*degree real parameters, each subdivision gives one real observation
	int n = 1000 * degree;
	double t0, t1, s0, s1;
	targetCurve1->GetDomain(&t0, &t1);
	targetCurve1->GetDomain(&s0, &s1);

	std::vector<double> targetPoints1x, targetPoints1y, targetPoints1z;
	std::vector<double> targetPoints2x, targetPoints2y, targetPoints2z;

	targetPoints1x.reserve(n);
	targetPoints1y.reserve(n);
	targetPoints1z.reserve(n);

	targetPoints2x.reserve(n);
	targetPoints2y.reserve(n);
	targetPoints2z.reserve(n);

	std::vector<double> s(n);
	std::vector<double> ta1(n);
	std::vector<double> ta2(n);

	for (int i = 0; i < n; i++) s[i] = (double)i / n;

	targetCurve1->GetNormalizedArcLengthPoints(n, s.data(), ta1.data());
	targetCurve2->GetNormalizedArcLengthPoints(n, s.data(), ta2.data());

	for (int i = 0; i < n; i++) {
		auto p1 = targetCurve1->PointAt(ta1[i]);
		targetPoints1x.push_back(p1.x);
		targetPoints1y.push_back(p1.y);
		targetPoints1z.push_back(p1.z);
		auto p2 = targetCurve2->PointAt(ta2[i]);
		targetPoints2x.push_back(p2.x);
		targetPoints2y.push_back(p2.y);
		targetPoints2z.push_back(p2.z);
	}

	double R1 = 0.5;
	double R2 = 1.5;
	
	AnnularLaplaceData kx(targetPoints1x, R1, targetPoints2x, R2, degree);
	AnnularLaplaceData ky(targetPoints1y, R1, targetPoints2y, R2, degree);
	AnnularLaplaceData kz(targetPoints1z, R1, targetPoints2z, R2, degree);

	ON_Mesh disk = CreateUncappedCylinder(nrVerticesAround, nrVerticesVertical);

	int nv = disk.VertexCount();
	for (int i = 0; i < nv; i++) {
		const ON_3dPoint& p = disk.Vertex(i);
		double f = (R1 + p.z * (R2 - R1)) / sqrt(p.x * p.x + p.y * p.y);
		auto goal = ON_3dPoint(
			kx.eval(f * p.x, f * p.y),
			ky.eval(f * p.x, f * p.y),
			kz.eval(f * p.x, f * p.y)
		);

		disk.SetVertex(i, goal);
	}

	RhinoApp().ActiveDoc()->AddMeshObject(disk);
	context.m_doc.Redraw();
	
	ON_SimpleArray<ON_2dex> pairs;
	disk.GetClashingFacePairs(1, pairs);
	if (pairs.Count())
		RhinoApp().Print(L"There self-intersections. Try to rotate and/or flip the input curves. This might be easier with the MinSurface Grasshopper components.\n");

	return CRhinoCommand::success;
}
