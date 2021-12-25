#include "LaplaceData.h"
#include <numeric>
#include <math.h>

#ifndef M_PI
#define M_PI 3.1415926535
#endif

LaplaceData::LaplaceData(const std::vector<double>& targets, int degree) : targets(targets), k(targets.size()), degree(degree) {

	this->a0 = //std::accumulate(targets.begin(), targets.end(), 0) / k;
	std::reduce(targets.begin(), targets.end()) / k;

	this->an.reserve(degree); this->an.push_back(0);
	this->bn.reserve(degree); this->bn.push_back(0);

	int kk = this->k;
	// can obviously be parallelized but is not the bottleneck at all
	for (int i = 1; i < degree; i++)
	{
		int ii = 0;
		this->an.push_back(std::accumulate(targets.begin(), targets.end(), 0.0,
			[i, kk, &ii](double w, const double& v) { return w + 2 * v * cos((double)(i * ii++) / (double)kk * 2 * M_PI); }) / kk);

		ii = 0;
		this->bn.push_back(std::accumulate(targets.begin(), targets.end(), 0.0,
			[i, kk, &ii](double w, const double& v) { return w + 2 * v * sin((double)(i * ii++) / (double)kk * 2 * M_PI); }) / kk);

	}

}

double LaplaceData::eval(double x, double y) {
	double r = sqrt(x * x + y * y);
	double theta = atan2(y, x);
	//if (r > 1-1e-2) {
	//	int ind = (int)(this->k * (theta + 2 * M_PI) / 2 / M_PI);
	//	return this->targets[ind % this->k];
	//}
	double sum = this->a0;
	for (int i = 1; i < this->degree; i++)
		sum += pow(r, i) * (this->an[i] * cos(i * theta) + this->bn[i] * sin(i * theta));

	return sum;
}
