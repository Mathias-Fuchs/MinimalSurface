using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using MathNet.Numerics;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;
using Grasshopper.Kernel.Parameters;
// In order to load the result of this wizard, you will also need to
// add the output bin/ folder of this project to the list of loaded
// folder in Grasshopper.
// You can use the _GrasshopperDeveloperSettings Rhino command for that.

namespace MinSurface
{
    public class MinSurfaceComponent : GH_Component
    {
        int _n;
        int _degree;

        Mesh _domainMesh;
        Mesh _outputMesh;
        // required data

        Curve _targetCurve;
        Curve _targetCurve2;

        List<Point3d> _targetPoints;
        List<Point3d> _targetPoints2;

        int _nrVert;


        private struct LaplaceData
        {
            public double a0;
            public double[] an;
            public double[] bn;
        }

        public MinSurfaceComponent()
          : base("Mini", "mini",
              "Computes an approximate minimal surface (soap film) inside two given closed curves.",
              "Mathias", "Geometry")
        {
        }


        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
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

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
            => pManager.AddMeshParameter("Output mesh", "outputMesh", "The output mesh, an approximate minimal surface.", GH_ParamAccess.item);

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
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


            Init(tc, tc2, k);//22ms
            // now, do the actual work and compute the three complex polynomials
            // ok, let's get the coefficients
            List<double> xs = _targetPoints.Select(o => o.X).ToList();// 1 ms
            LaplaceData kx = solveLaplace2(xs);
            List<double> ys = _targetPoints.Select(o => o.Y).ToList();
            LaplaceData ky = solveLaplace2(ys);
            List<double> zs = _targetPoints.Select(o => o.Z).ToList();
            LaplaceData kkz = solveLaplace2(zs);

            _domainMesh = Mesh.CreateFromPlanarBoundary((new Circle(1.0)).ToNurbsCurve(), new MeshingParameters(), .001);

            // this loop takes 170ms for 10000 vertices ... so, it is the bottleneck of the entire procedure
            // so it's worth to parallelize it
            for (int ii = 0; ii < _domainMesh.Vertices.Count; ii++)
            {
                var p = new Point2d(_domainMesh.Vertices[ii].X, _domainMesh.Vertices[ii].Y);

                double newx2 = evalPolynomial2(kx, p);
                double newy2 = evalPolynomial2(ky, p);
                double newz2 = evalPolynomial2(kkz, p);

                _domainMesh.Vertices.SetVertex(ii, newx2, newy2, newz2);
            }

            DA.SetData(0, _domainMesh);

        }


        private double evalPolynomial2(LaplaceData k, Point2d p)
        {
            var r = Math.Sqrt(p.X * p.X + p.Y * p.Y);
            var theta = Math.Atan2(p.Y, p.X);
            return k.a0 + Enumerable.Range(1, _degree).Sum(i => Math.Pow(r, i) * (k.an[i] * Math.Cos(i * theta) + k.bn[i] * Math.Sin(i * theta)));
        }


        // solves the equation Laplace(u) = 0 on the unit disk, 
        private LaplaceData solveLaplace2(List<double> targets)
        {
            var k = new LaplaceData();
            k.a0 = targets.Average();
            k.an = new double[_degree + 1];
            k.bn = new double[_degree + 1];

            // can obviously be parallelized
            for (int i = 1; i < _degree; i++)
            {
                k.an[i] = targets.Select((double v, int ii) => 2.0 * v * Math.Cos((double)(i * ii) / (double)targets.Count * 2 * Math.PI)).Average();
                k.bn[i] = targets.Select((double v, int ii) => 2.0 * v * Math.Sin((double)(i * ii) / (double)targets.Count * 2 * Math.PI)).Average();
            }

            return k;
        }



        public void Init(Curve tc, Curve tc2, int k)
        {
            Control.UseBestProviders();
            
            _targetCurve = tc;
            _targetCurve2 = tc2;

            _nrVert = k;


            // number of control points, tells about the complexity of the curve
            int nrCont = Math.Max(_targetCurve.ToNurbsCurve().Points.Count, _targetCurve2.ToNurbsCurve().Points.Count);

            // the degree of the curve
            int crDeg = Math.Max(_targetCurve.Degree, _targetCurve2.Degree);

            int idealDegree = nrCont * crDeg;

            // Chosen degree of complex analytic polynomials
            // we cap the degree to lie within 15 and 75 ... anything beyond would me meaningless
            _degree = Math.Min(Math.Max(15, idealDegree), 75);

            //  number of boundary subdivisions
            _n = 23 * _degree;

            double[] t = _targetCurve.DivideByCount(_n, true);
            _targetPoints = new List<Point3d>();
            for (int i = 0; i < _n; i++)
            {
                _targetPoints.Add(_targetCurve.PointAt(t[i]));
            }

            double[] t2 = _targetCurve2.DivideByCount(_n, true);
            _targetPoints2 = new List<Point3d>();
            for (int i = 0; i < _n; i++)
            {
                _targetPoints2.Add(_targetCurve2.PointAt(t[i]));
            }

        }
    }

}
