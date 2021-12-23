using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Linq;
using System.IO;
using System.Reflection;

namespace MinSurface
{
    public class MinSurfaceComponent : GH_Component
    {

        public MinSurfaceComponent()
          : base("MiniTwo", "minitwo",
              "Computes an approximate minimal surface (soap film) inside two given closed curves.",
              "Mathias", "Geometry")
        {
        }

        public override Guid ComponentGuid
        {
            get
            {
                return new Guid("6679fe76-1914-4cf2-a2da-a3b12cef0ff3");
            }
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Input curve (required)", "inputCurve", "The first closed" +
                " input curve." +
                "The curve can either be smooth, or be a polyline, or a mixture. " +
                "The only requirement is that it be closed. " +
                "This parameter is required.", GH_ParamAccess.item);

            pManager.AddCurveParameter("Input curve (required)", "inputCurve", "The second closed" +
    " input curve." +
    "The curve can either be smooth, or be a polyline, or a mixture. " +
    "This parameter is required.", GH_ParamAccess.item);

            pManager.AddIntegerParameter("nrVerticesVertical (optional)", "nrVerticesVertical", "number of vertices in vertical direction.", GH_ParamAccess.item, 20);
            pManager.AddIntegerParameter("nrVerticesAround (optional)", "nrVerticesAround", "number of vertices around the cylinder.", GH_ParamAccess.item, 18);
            pManager.AddIntegerParameter("Degree (optional)", "degree", "degree of the surface.", GH_ParamAccess.item, 0);

            pManager.AddAngleParameter("Rotation angle of one of the curves around itself (optional)", "roationAngle", "rotation angle of the curves around itself", GH_ParamAccess.item, 0);

            pManager.AddBooleanParameter("Flip one curve (optional)", "flip", "flip one curve? try this to prevent self-intersection.", GH_ParamAccess.item, false);

        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Output mesh", "outputMesh", "The output mesh, an approximate minimal surface.", GH_ParamAccess.item);
            //            pManager.AddNumberParameter("Non-orthogonality", "", "", GH_ParamAccess.item);
        }

        private static Mesh UncappedCylinder(int axis_divisions, int height_divisions)
        {

            Mesh M = new Mesh();
            if (axis_divisions <= 0 || height_divisions <= 0) return M;
            M.Vertices.Capacity = axis_divisions * height_divisions;
            M.Faces.Capacity = 2 * (axis_divisions * (height_divisions - 1));

            for (int th = 0; th < axis_divisions; th++)
            {
                double x = Math.Cos(2.0 * Math.PI * th / (double)axis_divisions);
                double y = Math.Sin(2.0 * Math.PI * th / (double)axis_divisions);

                for (int h = 0; h < height_divisions; h++)
                {
                    double z = (double)h / (double)(height_divisions - 1);
                    M.Vertices.Add(x, y, z);
                }
            }

            for (int th = 0; th < axis_divisions; th++)
            {
                for (int h = 1; h < height_divisions; h++)
                {
                    M.Faces.AddFace(
                     ((th + 0) % axis_divisions) * height_divisions + (h - 1),
                     ((th + 1) % axis_divisions) * height_divisions + (h - 1),
                     ((th + 0) % axis_divisions) * height_divisions + (h + 0)
                     );

                    M.Faces.AddFace(
                     ((th + 1) % axis_divisions) * height_divisions + (h - 1),
                     ((th + 1) % axis_divisions) * height_divisions + (h + 0),
                     ((th + 0) % axis_divisions) * height_divisions + (h + 0)
                     );
                }
            }
            return M;
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            System.Reflection.Assembly assembly = System.Reflection.Assembly.GetExecutingAssembly();
            System.Diagnostics.FileVersionInfo fvi = System.Diagnostics.FileVersionInfo.GetVersionInfo(assembly.Location);
            string version = fvi.FileVersion;
            Rhino.RhinoApp.WriteLine("Minimal surface component, version " + version);

// the input curve
            Curve tc = null;
            Curve tc2 = null;
            int vertical = 0;
            int around = 0;
            int degree = 0;
            double angle = 0;
            bool flip = false;

            DA.GetData(0, ref tc);
            DA.GetData(1, ref tc2);
            DA.GetData(2, ref vertical);
            DA.GetData(3, ref around);
            DA.GetData(4, ref degree);
            DA.GetData(5, ref angle);
            DA.GetData(6, ref flip);

            if (tc == null || !tc.IsValid || !tc.IsClosed) { this.AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "First input curve is either unvalid or not closed!"); return; }
            if (tc2 == null || !tc2.IsValid || !tc2.IsClosed) { this.AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Second input curve is either unvalid or not closed!"); return; }

            Curve _targetCurve = tc;
            Curve _targetCurve2 = tc2;

            // number of control points, tells about the complexity of the curve
            int nrCont = _targetCurve.ToNurbsCurve().Points.Count;
            int crDeg = _targetCurve.Degree;

            int idealDegree = nrCont * crDeg;
            int _degree = Math.Min(Math.Max(10, idealDegree), 50);

            //  number of boundary subdivisions
            int _n = 23 * _degree;
            int nrCont2 = _targetCurve2.ToNurbsCurve().Points.Count;
            int crDeg2 = _targetCurve2.Degree;

            int idealDegree2 = nrCont2 * crDeg2;
            int _degree2 = Math.Min(Math.Max(10, idealDegree2), 50);
            int _n2 = 23 * _degree2;

            int deg;
            if (degree == 0) { deg = Math.Max(_degree, _degree2); } else { deg = degree; }

            double[] t = _targetCurve.DivideByCount(_n, true);
            var _targetPoints = Enumerable.Range(0, _n).Select(i => _targetCurve.PointAt(flip ? 1 - t[i] : t[i])).ToList();
            double[] t2 = _targetCurve2.DivideByCount(_n2, true);
            var _targetPoints2 = Enumerable.Range(0, _n2).Select(i => _targetCurve2.PointAt(t2[(i + (int)((double)_n2 * (angle / 2 * Math.PI))) % _n2])).ToList();

            double R1 = 0.5;
            double R2 = 1.5;

            // now, do the actual work and compute the three complex polynomials
            // ok, let's get the coefficients
            List<double> xs1 = _targetPoints.Select(o => o.X).ToList();// 1 ms
            List<double> xs2 = _targetPoints2.Select(o => o.X).ToList();// 1 ms
                                                                        //        LaplaceData kx = new LaplaceData(xs1, deg); <-- if it's just one curve
            AnnularLaplaceData akx = new AnnularLaplaceData(xs1, R1, xs2, R2, deg);

            List<double> ys1 = _targetPoints.Select(o => o.Y).ToList();// 1 ms
            List<double> ys2 = _targetPoints2.Select(o => o.Y).ToList();// 1 ms
                                                                        //       LaplaceData ky = new LaplaceData(ys1, deg);
            AnnularLaplaceData aky = new AnnularLaplaceData(ys1, R1, ys2, R2, deg);

            List<double> zs1 = _targetPoints.Select(o => o.Z).ToList();// 1 ms
            List<double> zs2 = _targetPoints2.Select(o => o.Z).ToList();// 1 ms
                                                                        //       LaplaceData kkz = new LaplaceData(zs1, deg);
            AnnularLaplaceData akkz = new AnnularLaplaceData(zs1, R1, zs2, R2, deg);

            var MM = UncappedCylinder(around, vertical);

            // bottleneck ... so here, parallelization is worth it
            var mvl = MM.Vertices.Select(pp =>
            {
                var x = pp.X;
                var y = pp.Y;
                var z = pp.Z;
                var f = (R1 + z * (R2 - R1)) / Math.Sqrt(x * x + y * y);
                var p = new Point2d(f * x, f * y);
                return new Point3d(akx.eval(p), aky.eval(p), akkz.eval(p));
            });

            var MMM = new Mesh();
            MMM.Vertices.Capacity = MM.Vertices.Count;
            MMM.Vertices.AddVertices(mvl);

            // bottleneck            
            MMM.Faces.Capacity = MM.Faces.Count;
            MMM.Faces.AddFaces(MM.Faces);

            DA.SetData(0, MMM);

        }
        protected override System.Drawing.Bitmap Icon
            => new System.Drawing.Bitmap(Properties.Resources.icon03);
    }
}
