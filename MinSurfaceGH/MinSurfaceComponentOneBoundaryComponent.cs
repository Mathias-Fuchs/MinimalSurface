using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Linq;
using System.Drawing;
using System.IO;

namespace MinSurface
{
    public class MinSurfaceComponentOneBoundaryComponent : GH_Component
    {
        public MinSurfaceComponentOneBoundaryComponent()
          : base("MiniOne", "miniOne",
              "Computes an approximate minimal surface (soap film) inside a closed curve.",
              "MinSurface", "Geometry") {}

        public override Guid ComponentGuid => new Guid("A857FF66-CD2E-4AEA-850A-555CBFE5D762");

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Input curve (required)", "inputCurve", "The first closed" +
                " input curve." +
                "The curve can either be smooth, or be a polyline, or a mixture. " +
                "The only requirement is that it be closed. " +
                "This parameter is required.", GH_ParamAccess.item);

            pManager.AddIntegerParameter("nrVerticesBoundary", "nrVerticesBoundary", "number of vertices along the boundary.", GH_ParamAccess.item, 300);
            pManager.AddIntegerParameter("Degree (optional)", "degree", "degree of the surface.", GH_ParamAccess.item, 0);
            pManager.AddAngleParameter("Rotation angle of the curve (optional)", "rotationAngle", "rotation angle of the curve", GH_ParamAccess.item, 0);

        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Output mesh", "outputMesh", "The output mesh, an approximate minimal surface.", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            System.Reflection.Assembly assembly = System.Reflection.Assembly.GetExecutingAssembly();
            System.Diagnostics.FileVersionInfo fvi = System.Diagnostics.FileVersionInfo.GetVersionInfo(assembly.Location);
            string version = fvi.FileVersion;
            Rhino.RhinoApp.WriteLine("Minimal surface component, version " + version);
            Rhino.RhinoApp.WriteLine("GPLv3 licensed, source: https://github.com/Mathias-Fuchs/MinimalSurface");
            Rhino.RhinoApp.WriteLine("Copyright Mathias Fuchs 2020 - 2021, https://mathiasfuchs.com");

            Curve tc = null;
            int nrBoundaryVertices = 0;
            int degree = 0;
            double angle = 0;

            DA.GetData(0, ref tc);
            DA.GetData(1, ref nrBoundaryVertices);
            DA.GetData(2, ref degree);
            DA.GetData(3, ref angle);

            if (tc == null || !tc.IsValid || !tc.IsClosed) { this.AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Input curve is either unvalid or not closed!"); return; }

            Curve _targetCurve = tc;

            // number of control points, tells about the complexity of the curve
            int nrCont = _targetCurve.ToNurbsCurve().Points.Count;
            int crDeg = _targetCurve.Degree;

            // if degree isn't specified, use this heuristic
            // computationally, 25 is always doable, and more than 300 doesn't usually make sense
            if (degree == 0) degree = Math.Min(Math.Max(25, nrCont * crDeg), 300);

            //  number of boundary subdivisions for computation of the polynomials
            int _n = 23 * degree;

            double[] t = _targetCurve.DivideByCount(_n, true);
            var _targetPoints = Enumerable.Range(0, _n).Select(i => _targetCurve.PointAt(t[(i + (int)(angle / Math.PI / 2.0 * _n)) % _n])).ToList();

            // now, do the actual work and compute the three complex polynomials
            // ok, let's get the coefficients
            List<double> xs1 = _targetPoints.Select(o => o.X).ToList(); // 1 ms
            LaplaceData kx = new LaplaceData(xs1, degree);

            List<double> ys1 = _targetPoints.Select(o => o.Y).ToList(); // 1 ms
            LaplaceData ky = new LaplaceData(ys1, degree);

            List<double> zs1 = _targetPoints.Select(o => o.Z).ToList(); // 1 ms
            LaplaceData kkz = new LaplaceData(zs1, degree);

            var c = new Circle(1.0);
            var kk = MeshingParameters.Default;
            kk.MinimumEdgeLength = 2 * Math.PI / (double)nrBoundaryVertices;
            kk.MaximumEdgeLength = 2 * Math.PI / (double)nrBoundaryVertices * 2;
            var MM = Mesh.CreateFromPlanarBoundary(c.ToNurbsCurve(), kk, DocumentTolerance());

            var MMM = new Mesh();
            var mvl = MM.Vertices.Select(pp =>
            {
                var p = new Point2d(pp.X, pp.Y);
                return new Point3d(kx.eval(p), ky.eval(p), kkz.eval(p));
            }
                );

            // bottleneck            
            MMM.Vertices.Capacity = MM.Vertices.Count;
            MMM.Vertices.AddVertices(mvl);
            MMM.Faces.Capacity = MM.Faces.Count;
            MMM.Faces.AddFaces(MM.Faces);

            if (MMM.Faces.GetClashingFacePairs(1).Any())
                this.AddRuntimeMessage(
                    GH_RuntimeMessageLevel.Warning,
                    "The resulting mesh has self-intersections. Modifying the rotation angle input parameter can solve this.");
            DA.SetData(0, MMM);
        }

        protected override System.Drawing.Bitmap Icon => Properties.Resources.icon03_vfe_icon.ToBitmap();

    }
}
