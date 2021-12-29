#include "LaplaceData.h"
#include <numeric>
#include <math.h>

#ifndef M_PI
#define M_PI 3.1415926535
#endif

LaplaceData::LaplaceData(const std::vector<double>& targets, int degree) : targets(targets), k(targets.size()), degree(degree) {

	this->a0 = std::reduce(targets.begin(), targets.end()) / k;

	this->an.reserve(degree); this->an.push_back(0);
	this->bn.reserve(degree); this->bn.push_back(0);

	int kk = this->k;

	// can obviously be parallelized but is not the bottleneck at all
	for (int i = 1; i < degree; i++) {
		int ii = 0;
		this->an.push_back(std::accumulate(targets.begin(), targets.end(), 0.0,
			[i, kk, &ii](double w, const double& v) { return w + 2 * v * cos((double)(i * ii++) / (double)kk * 2 * M_PI); }) / kk);

		ii = 0;
		this->bn.push_back(std::accumulate(targets.begin(), targets.end(), 0.0,
			[i, kk, &ii](double w, const double& v) { return w + 2 * v * sin((double)(i * ii++) / (double)kk * 2 * M_PI); }) / kk);
	}
}

double LaplaceData::eval(double x, double y) const {
	double r = sqrt(x * x + y * y);
	double theta = atan2(y, x);
	double sum = this->a0;
	for (int i = 1; i < this->degree; i++)
		sum += pow(r, i) * (this->an[i] * cos(i * theta) + this->bn[i] * sin(i * theta));

	return sum;
}


AnnularLaplaceData::AnnularLaplaceData(const std::vector<double>& targets1, double R1, const std::vector<double>& targets2, double R2, int degree) :
	R1(R1),
	R2(R2),
	an(std::vector<double>(degree + 1)),
	bn(std::vector<double>(degree + 1)),
	cn(std::vector<double>(degree + 1)),
	dn(std::vector<double>(degree + 1)),
	degree(degree) {

	double m1 = std::reduce(targets1.begin(), targets1.end()) / targets1.size();
	double m2 = std::reduce(targets2.begin(), targets2.end()) / targets2.size();

	b0 = (m1 - m2) / (log(R1) - log(R2));
	a0 = ((m1 - b0 * log(R1)) + (m2 - b0 * log(R2))) / 2.0; // these two summands should be the same.

	int kk1 = targets1.size();
	int kk2 = targets2.size();

	// can obviously be parallelized
	for (int i = 1; i < degree; i++) {
		int ii = 0;
		double g1i = std::accumulate(targets1.begin(), targets1.end(), 0.0,
			[i, kk1, &ii](double w, const double& v) { return w + 2 * v * cos((double)(i * ii++) / (double)kk1 * 2 * M_PI); }) / kk1;
		
		ii = 0;
		double g2i = std::accumulate(targets2.begin(), targets2.end(), 0.0,
			[i, kk2, &ii](double w, const double& v) { return w + 2 * v * cos((double)(i * ii++) / (double)kk2 * 2 * M_PI); }) / kk2;

		double det = pow(R1 / R2, i) - pow(R2 / R1, i);
		an[i] = 1.0 / det * ( pow(R2, -i) * g1i  - pow(R1, -i) * g2i );
		bn[i] = 1.0 / det * (-pow(R2,  i) * g1i  + pow(R1,  i) * g2i );

//		var g1ic = targets1.Select((double v, int ii) = > 2.0 * v * Math.Sin((double)(i * ii) / (double)targets1.Count * 2.0 * Math.PI)).Average();
//		var g2ic = targets2.Select((double v, int ii) = > 2.0 * v * Math.Sin((double)(i * ii) / (double)targets2.Count * 2.0 * Math.PI)).Average();

		ii = 0;
		double g1ic = std::accumulate(targets1.begin(), targets1.end(), 0.0,
			[i, kk1, &ii](double w, const double& v) { return w + 2 * v * sin((double)(i * ii++) / (double)kk1 * 2 * M_PI); }) / kk1;

		ii = 0;
		double g2ic = std::accumulate(targets2.begin(), targets2.end(), 0.0,
			[i, kk2, &ii](double w, const double& v) { return w + 2 * v * sin((double)(i * ii++) / (double)kk2 * 2 * M_PI); }) / kk2;


		cn[i] = 1.0 / det * ( pow(R2, -i) * g1ic - pow(R1, -i) * g2ic);
		dn[i] = 1.0 / det * (-pow(R2,  i) * g1ic + pow(R1,  i) * g2ic);
	}
}

double AnnularLaplaceData::eval(double x, double y) const {
	double r = sqrt(x * x + y * y);
	double theta = atan2(y, x);
	double sum = a0 + b0 * log(r);

	for (int i = 1; i < this->degree; i++)
		sum += 
		(pow(r, i) * an[i] + pow(r, -i) * bn[i]) * cos(i * theta) +
		(pow(r, i) * cn[i] + pow(r, -i) * dn[i]) * sin(i * theta);
	
	return sum;
}

//double AnnularLaplaceData::ddr(double x, double y) const
//{
//	double r = sqrt(x * x + y * y);
//	double theta = atan2(y, x);
//
//	var gg = this.b0 / r + Enumerable.Range(1, this.degree).Sum(i = >
//		((double)i * Math.Pow(r, i - 1) * this.an[i]
//			- (double)i * Math.Pow(r, -i - 1) * this.bn[i]) * Math.Cos(i * theta) +
//		((double)i * Math.Pow(r, i - 1) * this.cn[i]
//			- (double)i * Math.Pow(r, -i - 1) * this.dn[i]) * Math.Sin(i * theta)
//	);
//	return gg;
//}
//
//
//double AnnularLaplaceData::ddtheta(double x, double y) const
//{
//	double r = sqrt(x * x + y * y);
//	double theta = atan2(y, x);
//
//	return Enumerable.Range(1, this.degree).Sum(i = >
//		(Math.Pow(r, i) * this.an[i] + Math.Pow(r, -i) * this.bn[i]) * (-(double)i * Math.Sin(i * theta)) +
//		(Math.Pow(r, i) * this.cn[i] + Math.Pow(r, -i) * this.dn[i]) * ((double)i * Math.Cos(i * theta))
//	);
//}
//
//double AnnularLaplaceData::drtimesdtheta(Point2d p) const
//{
//	var gg = ddr(p) * ddtheta(p);
//	return gg;
//}
