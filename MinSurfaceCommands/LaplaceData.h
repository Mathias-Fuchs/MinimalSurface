#pragma once
#include <vector>

class LaplaceData
{
public:
	// solves the equation Laplace(u) = 0 on the inside of the unit disk.
	LaplaceData(const std::vector<double>& targets, int degree);
	double eval(double x, double y) const;

private:
	double a0;
	int degree;
	std::vector<double> an;
	std::vector<double> bn;
	int k;
	std::vector<double> targets;
};

class AnnularLaplaceData {
public:
	double R1, R2;
private:
	double a0, b0;
	std::vector<double> an;
	std::vector<double> bn;
	std::vector<double> cn;
	std::vector<double> dn;
	int degree;
public:


	AnnularLaplaceData(const std::vector<double>& targets1, double R1, const std::vector<double>& targets2, double R2, int degree);
	double eval(double x, double y) const;
	// double ddr(double x, double y) const;
	//double ddtheta(double x, double y) const;
	//double drtimesdtheta(double x, double y) const;
};
