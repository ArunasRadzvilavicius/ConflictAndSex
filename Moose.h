// The main simulation object -- headers
#ifndef GUARD_MOOSE_H
#define GUARD_MOOSE_H
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using std::vector;
typedef vector<int> Cell;
class Moose
{
public:
	Moose();
	Moose(int Nx, int Mx, double mux, double kx, double sx, double sigmax, double Rx);
	int N;
	double R;
	double k;
	double s;
	double sigma;
	vector<Cell> Population;
	gsl_rng * rngm;
	void Selection();
	void MitoSelection();
	void Mutation();
	void Reproduction();
	vector<double> Sample(vector<double> list, int n);
	vector<Cell> Fuse(Cell cx1, Cell cx2);
	vector<Cell> ASexFuse(Cell cx1, Cell cx2);
	double AvgFitness();
	double GetFitness(Cell cx);
	double mu;
	int M;
	long GenSeed();
	double GetSexFreq();
	void PrintDistribution(int FG);
	int retval;
};
#endif
