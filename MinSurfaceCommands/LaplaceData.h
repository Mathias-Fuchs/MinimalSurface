#pragma once
#include <vector>

class LaplaceData
{
public:

	// solves the equation Laplace(u) = 0 on the inside of the unit disk. Input: the Dirichlet boundary values along the domain's boundary, sampled in an equispaced way
// the domain itself or its boundary are not needed, only a list of equi-spaced 
	LaplaceData(const std::vector<double>& targets, int degree) ;

	// we assume the point p lies inside the unit disk
	double eval(double x, double y);


private:
	double a0;
	std::vector<double> an;
	std::vector<double> bn;
	int degree;
	int k;
	std::vector<double> targets;
};

