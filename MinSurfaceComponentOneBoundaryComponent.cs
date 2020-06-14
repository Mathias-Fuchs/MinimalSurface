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
              "Computes an approximate minimal surface (soap film) inside one given closed curves.",
              "Mathias", "Geometry")
        {
        }

        public override Guid ComponentGuid
        {
            get
            {
                return new Guid("A857FF66-CD2E-4AEA-850A-555CBFE5D762");
            }
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Input curve (required)", "inputCurve", "The first closed" +
                " input curve." +
                "The curve can either be smooth, or be a polyline, or a mixture. " +
                "The only requirement is that it be closed. " +
                "This parameter is required.", GH_ParamAccess.item);

            pManager.AddIntegerParameter("nrVerticesBoundary", "nrVerticesBoundary", "number of vertices along the boundary.", GH_ParamAccess.item, 300);
            pManager.AddIntegerParameter("Degree (optional)", "degree", "degree of the surface.", GH_ParamAccess.item, 0);
            pManager.AddAngleParameter("Rotation angle of one of the curves around itself (optional)", "roationAngle", "rotation angle of the curves around itself", GH_ParamAccess.item, 0);

        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Output mesh", "outputMesh", "The output mesh, an approximate minimal surface.", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
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
            MMM.Vertices.AddVertices(mvl);
            MMM.Faces.AddFaces(MM.Faces);
            DA.SetData(0, MMM);
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                byte[] bc = System.Convert.FromBase64String("Qk14BgAAAAAAADYEAAAoAAAAGAAAABgAAAABAAgAAAAAAEICAADTIQAA0yEAAAAAAAAAAAAAAAAAAAEBAQACAgIAAwMDAAQEBAAFBQUABgYGAAcHBwAICAgACQkJAAoKCgALCwsADAwMAA0NDQAODg4ADw8PABAQEAAREREAEhISABMTEwAUFBQAFRUVABYWFgAXFxcAGBgYABkZGQAaGhoAGxsbABwcHAAdHR0AHh4eAB8fHwAgICAAISEhACIiIgAjIyMAJCQkACUlJQAmJiYAJycnACgoKAApKSkAKioqACsrKwAsLCwALS0tAC4uLgAvLy8AMDAwADExMQAyMjIAMzMzADQ0NAA1NTUANjY2ADc3NwA4ODgAOTk5ADo6OgA7OzsAPDw8AD09PQA + Pj4APz8 / AEBAQABBQUEAQkJCAENDQwBEREQARUVFAEZGRgBHR0cASEhIAElJSQBKSkoAS0tLAExMTABNTU0ATk5OAE9PTwBQUFAAUVFRAFJSUgBTU1MAVFRUAFVVVQBWVlYAV1dXAFhYWABZWVkAWlpaAFtbWwBcXFwAXV1dAF5eXgBfX18AYGBgAGFhYQBiYmIAY2NjAGRkZABlZWUAZmZmAGdnZwBoaGgAaWlpAGpqagBra2sAbGxsAG1tbQBubm4Ab29vAHBwcABxcXEAcnJyAHNzcwB0dHQAdXV1AHZ2dgB3d3cAeHh4AHl5eQB6enoAe3t7AHx8fAB9fX0Afn5 + AH9 / fwCAgIAAgYGBAIKCggCDg4MAhISEAIWFhQCGhoYAh4eHAIiIiACJiYkAioqKAIuLiwCMjIwAjY2NAI6OjgCPj48AkJCQAJGRkQCSkpIAk5OTAJSUlACVlZUAlpaWAJeXlwCYmJgAmZmZAJqamgCbm5sAnJycAJ2dnQCenp4An5 + fAKCgoAChoaEAoqKiAKOjowCkpKQApaWlAKampgCnp6cAqKioAKmpqQCqqqoAq6urAKysrACtra0Arq6uAK + vrwCwsLAAsbGxALKysgCzs7MAtLS0ALW1tQC2trYAt7e3ALi4uAC5ubkAurq6ALu7uwC8vLwAvb29AL6 + vgC / v78AwMDAAMHBwQDCwsIAw8PDAMTExADFxcUAxsbGAMfHxwDIyMgAycnJAMrKygDLy8sAzMzMAM3NzQDOzs4Az8 / PANDQ0ADR0dEA0tLSANPT0wDU1NQA1dXVANbW1gDX19cA2NjYANnZ2QDa2toA29vbANzc3ADd3d0A3t7eAN / f3wDg4OAA4eHhAOLi4gDj4 + MA5OTkAOXl5QDm5uYA5 + fnAOjo6ADp6ekA6urqAOvr6wDs7OwA7e3tAO7u7gDv7 + 8A8PDwAPHx8QDy8vIA8 / PzAPT09AD19fUA9vb2APf39wD4 + PgA + fn5APr6 + gD7 +/ sA / Pz8AP39 / QD +/ v4A////AMO+vbu6uLOurKqmpqapq7C2ur7AwsPFyLu6ubi2s66pqKWjpZ+ksbW3ur2/wMLEw7+9u7e0r7CvrKeorKmhrLW8v8HCxMO+xMK+ubWxtaaKorCusrevqbC1ubu9vLe5xr66uLW3s4N0k7e1t7u4sbOztr65tbS/yLu2uLjCl3B5ibe8ur2+tbPAvZqovLnExr66uLm6gG9pfrvDv8G9vMSthn+uvL/FxcC7ur+xdXBddcPLx8vFsJB5gKe9usTDxsC7vsS/fWBmape1r56FcHSPmbW7wMXDx8C+ubbMtHxzbGluZFphfZOIor+8xcTEx8LGrIao0dSejnxGU2p5gn2HuL7ExMTEx8LGroaJoKqTv5hMYnR5eW+kx8TGxMTFyMLGs42QioGX1qpsgY2LepbFwcTExcXFyMLEuZCPkIybwKaNoaebmcLHv7/AwsTFycLBwJeLj4yQm5SbrrCnu8e/v8DAv8DBx8K/yKKGjouDe32XrKy2xMHAwMDAwMDAw8K/yLCHjYh6aGqImq3Ewb/BwcDAwcHBxMK/wr+Qh4V1YWJ0lsLFwMDBwMDBwcHCxMK/vsmkfoNyYlZzvMnAwcHAwcHBwcLCxcK/v8XBhnZwUV6zzsDBwcHBwsLCwcLCxcK/v7/KtXxldLXQwMDBwcLCwsLDwsLDxsPAwMC/ybutyM6/wcHCwsLCwsLCwsPEx8DAwcDAwMfNxcDBwcLCwsLCwsPCw8TExMfExMTFxcXGxcXGxcbHxsbHx8fIyMfIywAA");
                System.Drawing.Bitmap bmp;
                using (var ms = new MemoryStream(bc))
                {
                    bmp = new System.Drawing.Bitmap(ms);
                }

                return bmp;
            }
        }

    }
}
