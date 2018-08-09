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

        // solves the equation Laplace(u) = 0 on the unit disk, 
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
                var g1i = targets1.Select((double v, int ii) => 2.0 * v * Math.Cos((double)(i * ii) / (double) targets1.Count * 2.0 * Math.PI)).Average();
                var g2i = targets2.Select((double v, int ii) => 2.0 * v * Math.Cos((double)(i * ii) / (double) targets2.Count * 2.0 * Math.PI)).Average();
                var det = Math.Pow(R1 / R2, i) - Math.Pow(R2 / R1, i);
                this.an[i] = 1.0 / det * (Math.Pow(R2, -i) * g1i - Math.Pow(R1, -i) * g2i);
                this.bn[i] = 1.0 / det * (-Math.Pow(R2, i) * g1i + Math.Pow(R1, i) * g2i);

                var g1ic = targets1.Select((double v, int ii) => 2.0 * v * Math.Sin((double)(i * ii) / (double)targets1.Count * 2.0 * Math.PI)).Average();
                var g2ic = targets2.Select((double v, int ii) => 2.0 * v * Math.Sin((double)(i * ii) / (double)targets2.Count * 2.0 * Math.PI)).Average();
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
            //0

            pManager.AddCurveParameter("Input curve (required)", "inputCurve", "The closed" +
                " input curve." +
                "The curve can either be smooth, or be a polyline, or a mixture. " +
                "The only requirement is that it be smooth. " +
                "This parameter is the only required one, all others are optional.", GH_ParamAccess.item);

            pManager.AddCurveParameter("Input curve (required)", "inputCurve", "The closed" +
    " input curve." +
    "The curve can either be smooth, or be a polyline, or a mixture. " +
    "The only requirement is that it be smooth. " +
    "This parameter is the only required one, all others are optional.", GH_ParamAccess.item);

            //5
            pManager.AddIntegerParameter("Approximate number of vertices", "appNrVertices", "approximate number" +
                " of vertices of the output mesh. The actual number of vertices may be slightly higher." +
                "The default to 1000.", GH_ParamAccess.item, 1000);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
            => pManager.AddMeshParameter("Output mesh", "outputMesh", "The output mesh, an approximate minimal surface.", GH_ParamAccess.item);


        protected override void SolveInstance(IGH_DataAccess DA)
        {// the input curve
            Curve tc = null;
            Curve tc2 = null;

            // the number of vertices
            int k = 1000;

            DA.GetData(0, ref tc);
            DA.GetData(1, ref tc2);
            DA.GetData(2, ref k);



            if (tc == null || !tc.IsValid || !tc.IsClosed) { this.AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "first input curve unvalid or not closed."); return; }
            if (tc2 == null || !tc2.IsValid || !tc2.IsClosed) { this.AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "second input curve unvalid or not closed."); return; }


            Curve _targetCurve = tc;
            Curve _targetCurve2 = tc2;


            // number of control points, tells about the complexity of the curve
            int nrCont = _targetCurve.ToNurbsCurve().Points.Count;

            // the degree of the curve
            int crDeg = _targetCurve.Degree;

            int idealDegree = nrCont * crDeg;

            // Chosen degree of complex analytic polynomials
            // we cap the degree to lie within 15 and 75 ... anything beyond would me meaningless
            int _degree = Math.Min(Math.Max(15, idealDegree), 75);

            //  number of boundary subdivisions
            int _n = 23 * _degree;


            // number of control points, tells about the complexity of the curve
            int nrCont2 = _targetCurve2.ToNurbsCurve().Points.Count;

            // the degree of the curve
            int crDeg2 = _targetCurve2.Degree;

            int idealDegree2 = nrCont2 * crDeg2;

            // Chosen degree of complex analytic polynomials
            // we cap the degree to lie within 15 and 75 ... anything beyond would me meaningless
            int _degree2 = Math.Min(Math.Max(15, idealDegree2), 75);

            //  number of boundary subdivisions
            int _n2 = 23 * _degree2;





            double[] t = _targetCurve.DivideByCount(_n, true);
            var _targetPoints = new List<Point3d>();
            for (int i = 0; i < _n; i++) _targetPoints.Add(_targetCurve.PointAt(t[i]));

            double[] t2 = _targetCurve2.DivideByCount(_n2, true);
            var _targetPoints2 = new List<Point3d>();
            for (int i = 0; i < _n2; i++) _targetPoints2.Add(_targetCurve2.PointAt(t[i]));

            // now, do the actual work and compute the three complex polynomials
            // ok, let's get the coefficients
            List<double> xs = _targetPoints.Select(o => o.X).ToList();// 1 ms
            LaplaceData kx = new LaplaceData(xs, _degree);
            List<double> ys = _targetPoints.Select(o => o.Y).ToList();
            LaplaceData ky = new LaplaceData(ys, _degree);
            List<double> zs = _targetPoints.Select(o => o.Z).ToList();
            LaplaceData kkz = new LaplaceData(zs, _degree);

            var _domainMesh = Mesh.CreateFromPlanarBoundary((new Circle(1.0)).ToNurbsCurve(), new MeshingParameters(), .001);
            

            // this loop takes 170ms for 10000 vertices ... so, it is the bottleneck of the entire procedure
            // so it's worth to parallelize it
            for (int ii = 0; ii < _domainMesh.Vertices.Count; ii++)
            {
                var p = new Point2d(_domainMesh.Vertices[ii].X, _domainMesh.Vertices[ii].Y);
                double newx2 = kx.eval(p);
                double newy2 = ky.eval(p); // evalPolynomial2(ky, p);  <-- maybe it's faster because it will be inlined?
                double newz2 = kkz.eval(p); // evalPolynomial2(kkz, p);
                _domainMesh.Vertices.SetVertex(ii, newx2, newy2, newz2);
            }
            DA.SetData(0, _domainMesh);
        }
    }
}
