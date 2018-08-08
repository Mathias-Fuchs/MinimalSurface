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
    public enum myMeshType
    {
        grid, cobweb, custom
    }

    public class MinSurfaceComponent : GH_Component, IGH_VariableParameterComponent
    {
        int _n;
        int _degree;

        Mesh _domainMesh;
        Mesh _outputMesh;
        // required data
        List<Point3d> _targetPoints;
        int _nrVert;
        Double _weierstrassAngle;

        Curve _targetCurve;
        //Double _meanX;
        //Double _meanY;
        //Double _meanZ;

        myMeshType _mt;
        Double _rotAngle;

        List<System.Numerics.Complex> _domain;
        Boolean _parallel;



        public MinSurfaceComponent()
          : base("Mini", "mini",
              "Computes an approximate minimal surface (soap film) inside a given closed curve.",
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

            //1
            pManager.AddAngleParameter("Rotation angle", "rotationAngle", "An angle to rotate" +
                " the approximate minimal surface mesh within the resulting geometry." +
                " This parameter is sometimes helpful in optimizing the resulting mesh." +
                "The default is zero.", GH_ParamAccess.item, 0);

            //2
            pManager.AddAngleParameter("Weierstrass angle", "weierstrassAngle", "This angle " +
              "lets you try your luck by exploring the design space that Mathematics" +
              " associates with the given minimal surface for free. This parameter was baptized with the " +
              "creator of the Weierstrass parametrization of minimal surfaces in mind." +
              "The default setting is zero, resulting in the original minimal surface.", GH_ParamAccess.item, 0.0);

            //3
            pManager.AddBooleanParameter("Suppress parallel execution", "suppress parallel execution", "On some machines," +
                " parallel execution causes problems. You can use this toggle to enable parallelity nevertheless, if " +
                "you want to try to speed the solution up. " +
                "The default is true (no acceleration).",
                GH_ParamAccess.item, true);

            //4
            pManager.AddIntegerParameter("Type of mesh", "Type of mesh", "The desired type of output mesh. Defaults to grid mesh (0)." +
                "You can change this to a cobweb mesh (1). If you want to supply your own mesh, activate the parameter in the component menu, " +
                "by right-clicking the component.",
                GH_ParamAccess.item, 0);

            var param = pManager[4] as Grasshopper.Kernel.Parameters.Param_Integer;

            param.AddNamedValue("Grid mesh (default)", 0);
            param.AddNamedValue("Cobweb mesh", 1);

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
            // the number of vertices
            int k = 1000;
            // the Weierstrass angle
            Double a = 0.0;
            // the angle to rotate the mesh about
            Double rotAngle = 0.0;
            // the resulting mesh
            Mesh M = new Mesh();
            // whether to suppress parallel execution
            Boolean suppressParallel = true;

            int typeOfMesh = 0;

            myMeshType mt = myMeshType.grid;

            DA.GetData(0, ref tc);

            DA.GetData(1, ref rotAngle);

            DA.GetData(2, ref a);

            DA.GetData(3, ref suppressParallel);


            int c = this.Params.Input.Count;

            // if the custom mesh type menu point is deactivated, we have to find out which default mesh type is desired.
            if (c == 6)
            {
                DA.GetData(4, ref typeOfMesh);

                switch (typeOfMesh)
                {
                    case 0:
                        mt = myMeshType.grid;
                        break;
                    case 1:
                        mt = myMeshType.cobweb;
                        break;
                    default:
                        mt = myMeshType.grid;
                        break;
                }
                //vertices
                DA.GetData(5, ref k);
            }
            else
            {
                DA.GetData(4, ref M);
                if (M.IsValid && M.Vertices.Count > 0)
                {
                    mt = myMeshType.custom;
                }
                else
                {
                    this.AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "Something seems wrong with the supplied mesh. Falling back to default.");
                    mt = myMeshType.grid;
                    k = 1000;
                }
            }

            if (tc == null || !tc.IsValid) { this.AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No valid input curve supplied."); return; }

            Init(tc, k, a, mt, rotAngle, M, suppressParallel);//22ms
            // now, do the actual work and compute the three complex polynomials
            // ok, let's get the coefficients
            List<double> xs = _targetPoints.Select(o => o.X).ToList();// 1 ms
            MathNet.Numerics.LinearAlgebra.Vector<double> betax = solveLaplace(xs); // 8 ms

            List<double> ys = _targetPoints.Select(o => o.Y).ToList();
            MathNet.Numerics.LinearAlgebra.Vector<double> betay = solveLaplace(ys);

            List<double> zs = _targetPoints.Select(o => o.Z).ToList();
            MathNet.Numerics.LinearAlgebra.Vector<double> betaz = solveLaplace(zs);

            System.Numerics.Complex weier = new System.Numerics.Complex(Math.Cos(_weierstrassAngle), Math.Sin(_weierstrassAngle));
            // now, apply the real parts of the complex polynomials to the mesh to obtain the three transformed coordinates


            if (!_parallel)
            {
                // this loop takes 170ms for 10000 vertices ... so, it is the bottleneck of the entire procedure
                // so it's worth to parallelize it
                for (int ii = 0; ii < _domainMesh.Vertices.Count; ii++)
                {
                    System.Numerics.Complex z = new System.Numerics.Complex(_domainMesh.Vertices[ii].X, _domainMesh.Vertices[ii].Y);
                    double newx = (weier * evalPolynomial(betax, z)).Real;
                    double newy = (weier * evalPolynomial(betay, z)).Real;
                    double newz = (weier * evalPolynomial(betaz, z)).Real;
                    _outputMesh.Vertices.SetVertex(ii, newx, newy, newz);
                }
            }
            else
            {
                Parallel.For(
                    0,
                    _domainMesh.Vertices.Count,
                    ii =>
                    {
                        System.Numerics.Complex z = new System.Numerics.Complex(_domainMesh.Vertices[ii].X, _domainMesh.Vertices[ii].Y);
                        double newx = (weier * evalPolynomial(betax, z)).Real;
                        double newy = (weier * evalPolynomial(betay, z)).Real;
                        double newz = (weier * evalPolynomial(betaz, z)).Real;
                        _outputMesh.Vertices.SetVertex(ii, newx, newy, newz);
                    });
            }
            // DA.SetData(0, _domainMesh);
            DA.SetData(0, _outputMesh);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public System.Numerics.Complex evalPolynomial(MathNet.Numerics.LinearAlgebra.Vector<double> coeffs, System.Numerics.Complex z)
        {


            // TODO: optimize this using Horner's rule
            // http://web.cecs.pdx.edu/~maier/cs584/Lectures/lect07b-11-MG.pdf
            System.Numerics.Complex sum = new System.Numerics.Complex(coeffs[0], 0);
            for (int pp = 1; pp < _degree + 1; pp++)
            {
                sum += z.Power(pp) * (new System.Numerics.Complex(coeffs[pp], coeffs[pp + _degree]));
            }
            return sum;
        }


        #region dealing with the menu
        protected override void AppendAdditionalComponentMenuItems(System.Windows.Forms.ToolStripDropDown menu)
        {
            Menu_AppendItem(menu, "Supply custom mesh instead of grid or cobweb mesh.", AdvancedMenuClicked, true, CustomMeshParameterPresent);
        }

        bool CustomMeshParameterPresent = false;

        private void AdvancedMenuClicked(object sender, EventArgs e)
        {
            RecordUndoEvent("Advanced Parameters");
            if (!CustomMeshParameterPresent)
            {
                if (this.Params.Input[4].SourceCount != 0 || this.Params.Input[5].SourceCount != 0)
                {
                    this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Please disconnect the last two parameters before" +
                        "selecting this option.");
                    Params.OnParametersChanged();
                    ExpireSolution(true);
                    return;
                }
                RemoveParametersFrom(4);
                AddMyMeshParameter();
            }
            else
            {
                if (this.Params.Input[4].SourceCount != 0)
                {
                    this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Please disconnect the last parameters before" +
                        "selecting this option.");
                    Params.OnParametersChanged();
                    ExpireSolution(true);
                    return;
                }
                RemoveParametersFrom(4);
                AddMyEnumAndIntParameter();

            }
            CustomMeshParameterPresent = !CustomMeshParameterPresent;
            Params.OnParametersChanged();
            ExpireSolution(true);
        }
        private void RemoveParametersFrom(int k)
        {
            List<IGH_Param> paramsToRemove = new List<IGH_Param>();
            for (int i = k; i < Params.Input.Count; i++)
            {
                if (Params.Input[i].SourceCount == 0)
                {
                    paramsToRemove.Add(Params.Input[i]);
                }
            }

            foreach (IGH_Param p in paramsToRemove)
            {
                Params.UnregisterInputParameter(p, true);
            }
        }

        //public override void AddedToDocument(GH_Document document)
        //{
        //    PopulateParams();
        //    base.AddedToDocument(document);
        //}
        private void AddMyMeshParameter()
        {

            IGH_Param pI = new Param_Mesh();
            pI.Name = "Disk mesh in XY plane";
            pI.NickName = pI.Name;
            pI.Description = "Your own mesh to be mapped onto the minimal surface." +
                "Please provide a flat mesh whose boundary is approximately a circle" +
                " in the XY plane, and which is entirely contained in the XY plane.";
            pI.Optional = false;
            //            pI.AddVolatileData(new GH_Path(0), 0, 1);
            //            myParams.Add(pI);
            if (this.Params.IndexOfInputParam(pI.Name) < 0)
                this.Params.RegisterInputParam(pI);
        }


        private void AddMyEnumAndIntParameter()
        {
            // 4th parameter
            IGH_Param i = new Param_Integer();

            i.Name = "Type of mesh";
            i.NickName = i.Name;
            i.Description = "The desired type of output mesh. Defaults to grid mesh (0)." +
                "You can change this to a cobweb mesh (1). If you want to supply your own mesh," +
                "activate the parameter in the component menu " +
                "by right-clicking the component.";
            i.Access = GH_ParamAccess.item;
            i.Optional = true;
            //     i.AddVolatileData(new GH_Path(), 0, 0);

            var parm = i as Grasshopper.Kernel.Parameters.Param_Integer;
            parm.AddNamedValue("grid mesh (default)", 0);
            parm.AddNamedValue("cobweb mesh", 1);

            if (this.Params.IndexOfInputParam(i.Name) < 0)
                this.Params.RegisterInputParam(i);

            //5th parameter
            IGH_Param j = new Param_Integer();
            j.Name = "approximate number of vertices";
            j.NickName = j.Name;
            j.Description = "approximate number of vertices of the output mesh. Defaults to 1000.";
            j.Access = GH_ParamAccess.item;
            j.Optional = true;
            //     j.AddVolatileData(new GH_Path(), 0, 1000);

            if (this.Params.IndexOfInputParam(j.Name) < 0)
                this.Params.RegisterInputParam(j);
        }


        bool IGH_VariableParameterComponent.CanInsertParameter(GH_ParameterSide side, int index)
        {
            return false;
        }
        bool IGH_VariableParameterComponent.CanRemoveParameter(GH_ParameterSide side, int index)
        {
            return false;
        }
        IGH_Param IGH_VariableParameterComponent.CreateParameter(GH_ParameterSide side, int index)
        {
            return null;
        }
        bool IGH_VariableParameterComponent.DestroyParameter(GH_ParameterSide side, int index)
        {
            return false;
        }
        void IGH_VariableParameterComponent.VariableParameterMaintenance()
        {
        }
        #endregion


        // solves the equation Laplace(u) = 0 on the unit disk, 
        public MathNet.Numerics.LinearAlgebra.Vector<double> solveLaplace(List<double> targets)
        {

            // a real linear regression, no complex linear regression anymore.
            GH_DebugDescriptionWriter.Equals(_domain.Count, targets.Count);
            int p = 2 * _degree + 1;

            var X = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix.Create(_n, p, 0);
            var Y = MathNet.Numerics.LinearAlgebra.Double.DenseVector.Create(_n, 0);

            for (int nn = 0; nn < _n; nn++)
            {

                X[nn, 0] = 1.0;

                //  System.Numerics.Complex z = _domain[nn];
                //  System.Numerics.Complex power = new System.Numerics.Complex();
                for (int pp = 1; pp < _degree + 1; pp++)
                {
                    // previously, this was
                    // power = RealPartOf(z.MathPow(pp));
                    // but since angle = nn * 2 * Constants.Pi / _n; we replace it with
                    X[nn, pp] = Math.Cos(pp * nn * 2 * Constants.Pi / _n); // the real part of the pp-th power
                                                                           // note the minus sign!
                    X[nn, pp + _degree] = -Math.Sin(pp * nn * 2 * Constants.Pi / _n); // the negative imaginary part of the pp-th power
                    // as a crazy optimization, one can compute the sine from the cosine, up to sign
                }
                // the "response variable"
                Y[nn] = targets[nn];
            }
            // try it with something like 
            // return MathNet.Numerics.LinearRegression.MultipleRegression.NormalEquations(X, Y);
            // another crazy optimization would consist in inverting t(X) * X algebraically beforehand.
            // or making sure it is actually the unit matrix
            // also, we have basically just computed a discrete Fourier transform of the data (some fumbling with the constants is required).
            // so, these parts can potentially be sped up. However, this is not the bottleneck, so it's not that urgent.
            return X.QR().Solve(Y);
        }


        public void Init(Curve tc, int k, Double a, myMeshType mt, Double rotAngle, Mesh M, Boolean suppressParallel)
        {
            try
            {
                // what's the difference between these?
                Control.UseBestProviders();
            }
            catch { }
            _targetCurve = tc;
            _nrVert = k;
            _weierstrassAngle = a;
            _mt = mt;
            _rotAngle = rotAngle;
            _parallel = !suppressParallel;

            if (!(_targetCurve == null) && !_targetCurve.IsClosed)
            {
                throw new System.Exception("The input curve has to be closed.");
            }

            // number of control points, tells about the complexity of the curve
            int nrCont = _targetCurve.ToNurbsCurve().Points.Count;

            // the degree of the curve
            int crDeg = _targetCurve.Degree;

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

            _domain = new List<System.Numerics.Complex>();
            for (int nn = 0; nn < _n; nn++)
            {
                double angle = nn * 2 * Constants.Pi / _n;
                _domain.Add(new System.Numerics.Complex(Math.Cos(angle), Math.Sin(angle)));
            }

            switch (_mt)
            {
                // manufacture a template mesh to be mapped onto the minimal surface
                case myMeshType.grid:
                    var meshParams = MeshingParameters.Default;
                    meshParams.MaximumEdgeLength = Math.Sqrt(Math.PI / _nrVert);
                    Curve _domainCurve = new ArcCurve(new Circle(1.0));
                    _domainMesh = Mesh.CreateFromPlanarBoundary(_domainCurve, meshParams, .0001);
                    break;

                case myMeshType.cobweb:
                    int nx = Math.Max(20, (int)Math.Ceiling(Math.Sqrt(_nrVert) / 1.5));
                    int ny = Math.Max(20, (int)(_nrVert / (nx + 0.0)));
                    _domainMesh = Mesh.CreateFromPlane(new Rhino.Geometry.Plane(0, 0, 1, 0), new Interval(0, 1), new Interval(0, 2.6), nx, ny);
                    System.Numerics.Complex jj = new System.Numerics.Complex(0, 1);
                    for (int i = 0; i < _domainMesh.Vertices.Count; i++)
                    {
                        System.Numerics.Complex stripPoint = new System.Numerics.Complex(_domainMesh.Vertices[i].X, _domainMesh.Vertices[i].Y);
                        System.Numerics.Complex circlePoint = System.Numerics.Complex.Exp(2 * Math.PI * jj * stripPoint);
                        _domainMesh.Vertices.SetVertex(i, new Point3d(circlePoint.Real, circlePoint.Imaginary, 0));
                    }
                    break;

                case myMeshType.custom:

                    var orig = new Point3f(0, 0, 0);

                    // project the input mesh on the XY-plane.
                    M.Transform(Transform.PlanarProjection(new Plane(0, 0, 1, 0)));

                    // center the mesh around the origin in the XY-plane
                    var bb = M.GetBoundingBox(false);
                    M.Translate(orig - bb.Center);

                    // now, scale it so that it fits into the unit disk.
                    M.Scale(2 * Math.Sqrt(2) / bb.Diagonal.Length);

                    // project all vertices that stick out of the unit disk back onto it
                    for (int i = 0; i < M.Vertices.Count; i++)
                    {
                        var d = M.Vertices[i].DistanceTo(orig);
                        if (d > 1)
                        {
                            M.Vertices.SetVertex(i, M.Vertices[i].X / d, M.Vertices[i].Y / d, 0);
                        }
                    }

                    _domainMesh = M;
                    break;
            }

            _domainMesh.Rotate(_rotAngle, new Vector3d(0, 0, 1), new Point3d(0, 0, 0));
            _outputMesh = _domainMesh.DuplicateMesh();
        }
    }

}
