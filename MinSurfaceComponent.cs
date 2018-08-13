using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Linq;

namespace MinSurface
{

    public class LaplaceData
    {
        private double a0;
        private double[] an;
        private double[] bn;
        private int degree;

        // solves the equation Laplace(u) = 0 on the unit disk, 
        public LaplaceData(List<double> targets, int degree)
        {
            this.a0 = targets.Average();
            this.an = new double[degree + 1];
            this.bn = new double[degree + 1];
            this.degree = degree;

            // can obviously be parallelized
            for (int i = 1; i < degree; i++)
            {
                this.an[i] = targets.Select((double v, int ii) => 2.0 * v * Math.Cos((double)(i * ii) / (double)targets.Count * 2 * Math.PI)).Average();
                this.bn[i] = targets.Select((double v, int ii) => 2.0 * v * Math.Sin((double)(i * ii) / (double)targets.Count * 2 * Math.PI)).Average();
            }
        }

        public double eval(Point2d p)
        {
            var r = Math.Sqrt(p.X * p.X + p.Y * p.Y);
            var theta = Math.Atan2(p.Y, p.X);
            // Horner's rule could speed it up but not bottleneck
            return this.a0 + Enumerable.Range(1, this.degree).Sum(i => Math.Pow(r, i) * (this.an[i] * Math.Cos(i * theta) + this.bn[i] * Math.Sin(i * theta)));
        }
    }

    public class AnnularLaplaceData
    {
        private double a0, b0;
        private double[] an;
        private double[] bn;
        private double[] cn;
        private double[] dn;
        private int degree;

        // solves the equation Laplace(u) = 0 on an annular region R1 <= R <= R2
        public AnnularLaplaceData(List<double> targets1, double R1, List<double> targets2, double R2, int degree)
        {

            this.an = new double[degree + 1];
            this.bn = new double[degree + 1];
            this.cn = new double[degree + 1];
            this.dn = new double[degree + 1];
            this.degree = degree;

            var m1 = targets1.Average();
            var m2 = targets2.Average();
            this.b0 = (m1 - m2) / (Math.Log(R1) - Math.Log(R2));
            this.a0 = ((m1 - this.b0 * Math.Log(R1)) + (m2 - this.b0 * Math.Log(R2))) / 2.0; // these two summands should be the same.

            // can obviously be parallelized
            for (int i = 1; i < degree; i++)
            {
                var g1i = 4.0 * Math.PI * targets1.Select((double v, int ii) =>  v * Math.Cos((double)(i * ii) / (double) targets1.Count )).Average();
                var g2i = 4.0 * Math.PI * targets2.Select((double v, int ii) =>  v * Math.Cos((double)(i * ii) / (double) targets2.Count )).Average();

                var det = Math.Pow(R1 / R2, i) - Math.Pow(R2 / R1, i);
                this.an[i] = 1.0 / det * (Math.Pow(R2, -i) * g1i - Math.Pow(R1, -i) * g2i);
                this.bn[i] = 1.0 / det * (-Math.Pow(R2, i) * g1i + Math.Pow(R1, i) * g2i);

                var g1ic = 4.0 * Math.PI * targets1.Select((double v, int ii) => v * Math.Sin((double) (i * ii) / (double) targets1.Count )).Average();
                var g2ic = 4.0 * Math.PI * targets2.Select((double v, int ii) => v * Math.Sin((double) (i * ii) / (double) targets2.Count )).Average();
                this.cn[i] = 1.0 / det * (Math.Pow(R2, -i) * g1ic - Math.Pow(R1, -i) * g2ic);
                this.dn[i] = 1.0 / det * (-Math.Pow(R2, i) * g1ic + Math.Pow(R1, i) * g2ic);
            }
        }

        public double eval(Point2d p)
        {
            var r = Math.Sqrt(p.X * p.X + p.Y * p.Y);
            var theta = Math.Atan2(p.Y, p.X);
            // Horner's rule could speed it up but not bottleneck
            return
                this.a0 + this.b0 * Math.Log(r) +
                Enumerable.Range(1, this.degree).Sum(i =>
                 (Math.Pow(r, i) * this.an[i] + Math.Pow(r, -i) * this.bn[i]) * Math.Cos(i * theta) +
                 (Math.Pow(r, i) * this.cn[i] + Math.Pow(r, -i) * this.dn[i]) * Math.Sin(i * theta)
                 );
        }
    }


    public class MinSurfaceComponent : GH_Component
    {

        public MinSurfaceComponent()
          : base("Mini", "mini",
              "Computes an approximate minimal surface (soap film) inside two given closed curves.",
              "Mathias", "Geometry")
        {
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("6679fe76-1914-4cf2-a2da-a3b12cef0ff3"); }
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Input curve (required)", "inputCurve", "The first closed" +
                " input curve." +
                "The curve can either be smooth, or be a polyline, or a mixture. " +
                "The only requirement is that it be smooth. " +
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
            pManager.AddPointParameter("target points (optional)", "", "", GH_ParamAccess.list, null);
            pManager.AddPointParameter("target points (optional)", "", "", GH_ParamAccess.list, null);

        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
            => pManager.AddMeshParameter("Output mesh", "outputMesh", "The output mesh, an approximate minimal surface.", GH_ParamAccess.item);


        protected override void SolveInstance(IGH_DataAccess DA)
        {// the input curve
            Curve tc = null;
            Curve tc2 = null;
            int vertical = 0;
            int around = 0;
            int degree = 0;
            double angle = 0;
            bool flip = false;
            List<Point3d> A = new List<Point3d>();
            List<Point3d> B = new List<Point3d>();

            DA.GetData(0, ref tc);
            DA.GetData(1, ref tc2);
            DA.GetData(2, ref vertical);
            DA.GetData(3, ref around);
            DA.GetData(4, ref degree);
            DA.GetData(5, ref angle);
            DA.GetData(6, ref flip);
            DA.GetDataList(7, A);
            DA.GetDataList(8, B);
            List<Point3d> _targetPoints;
            List<Point3d> _targetPoints2;
            int deg;
            if (A.Count == 0 || B.Count == 0)
            {

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

                if (degree == 0) { deg = Math.Max(_degree, _degree2); } else { deg = degree; }

                double[] t = _targetCurve.DivideByCount(_n, true);
                _targetPoints = Enumerable.Range(0, _n).Select(i => _targetCurve.PointAt(flip ? 1 - t[i] : t[i])).ToList();
                double[] t2 = _targetCurve2.DivideByCount(_n2, true);
                _targetPoints2 = Enumerable.Range(0, _n2).Select(i => _targetCurve2.PointAt(t2[(i + (int)((double)_n2 * (angle / 2 * Math.PI))) % _n2])).ToList();

            }
            else
            {
                _targetPoints = A;
                _targetPoints2 = B;
                deg = Math.Min(Math.Max(10, Math.Max(_targetPoints.Count, _targetPoints2.Count)) , 50);
            }

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

            var MM = Mesh.CreateFromCylinder(new Cylinder(new Circle(1), 1), vertical, around);

            // bottleneck
            for (int ii = 0; ii < MM.Vertices.Count; ii++)
            {
                var x = MM.Vertices[ii].X;
                var y = MM.Vertices[ii].Y;
                var z = MM.Vertices[ii].Z;
                var f = (R1 + z * (R2 - R1)) / Math.Sqrt(x * x + y * y);
                var p = new Point2d(f * x, f * y);
                MM.Vertices.SetVertex(ii, akx.eval(p), aky.eval(p), akkz.eval(p));
            }

            DA.SetData(0, MM);
            // TODO: find the angle such that the mesh goes off perpendicularly from the ground curve
        }
    }
}
