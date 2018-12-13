using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace MinSurface
{
    public class LaplaceData
    {
        private double a0;
        private double[] an;
        private double[] bn;
        private int degree;

        // solves the equation Laplace(u) = 0 on the inside of the unit disk. Input: the Dirichlet boundary values along the domain's boundary, sampled in an equispaced way
        // the domain itself or its boundary are not needed, only a list of equi-spaced 
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

        // we assume the point p lies inside the unit disk
        public double eval(Point2d p)
        {
            var r = Math.Sqrt(p.X * p.X + p.Y * p.Y);
            var theta = Math.Atan2(p.Y, p.X);
            // Horner's rule could speed it up but not bottleneck
            return this.a0 + Enumerable.Range(1, this.degree).Sum(i => Math.Pow(r, i) * (this.an[i] * Math.Cos(i * theta) + this.bn[i] * Math.Sin(i * theta)));
        }
    }
}
