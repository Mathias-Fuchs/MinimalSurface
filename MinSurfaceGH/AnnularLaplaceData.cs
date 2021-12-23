using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace MinSurface
{

    public class AnnularLaplaceData
    {
        public double R1, R2;
        private double a0, b0;
        private double[] an;
        private double[] bn;
        private double[] cn;
        private double[] dn;
        private int degree;

        // solves the equation Laplace(u) = 0 on the unit disk, 
        public AnnularLaplaceData(List<double> targets1, double R1, List<double> targets2, double R2, int degree)
        {
            this.R1 = R1;
            this.R2 = R2;
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
                var g1i = targets1.Select((double v, int ii) => 2.0 * v * Math.Cos((double)(i * ii) / (double)targets1.Count * 2.0 * Math.PI)).Average();
                var g2i = targets2.Select((double v, int ii) => 2.0 * v * Math.Cos((double)(i * ii) / (double)targets2.Count * 2.0 * Math.PI)).Average();
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

        public double ddr(Point2d p)
        {
            var r = Math.Sqrt(p.X * p.X + p.Y * p.Y);
            var theta = Math.Atan2(p.Y, p.X);

            var gg = this.b0 / r + Enumerable.Range(1, this.degree).Sum(i =>
                ((double)i * Math.Pow(r, i - 1) * this.an[i]
                - (double)i * Math.Pow(r, -i - 1) * this.bn[i]) * Math.Cos(i * theta) +
                ((double)i * Math.Pow(r, i - 1) * this.cn[i]
                - (double)i * Math.Pow(r, -i - 1) * this.dn[i]) * Math.Sin(i * theta)
      );
            return gg;
        }


        public double ddtheta(Point2d p)
        {
            var r = Math.Sqrt(p.X * p.X + p.Y * p.Y);
            var theta = Math.Atan2(p.Y, p.X);

            return Enumerable.Range(1, this.degree).Sum(i =>
                 (Math.Pow(r, i) * this.an[i] + Math.Pow(r, -i) * this.bn[i]) * (-(double)i * Math.Sin(i * theta)) +
                 (Math.Pow(r, i) * this.cn[i] + Math.Pow(r, -i) * this.dn[i]) * ((double)i * Math.Cos(i * theta))
                 );
        }

        public double drtimesdtheta(Point2d p)
        {
            var gg = ddr(p) * ddtheta(p);
            return gg;
        }


    }

}
