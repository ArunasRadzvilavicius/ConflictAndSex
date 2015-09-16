// Main simulation object
#include "Moose.h"
#include <iostream>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <unistd.h>
using namespace std;

Moose::Moose(){
}

Moose::Moose(int Nx, int Mx, double mux, double kx, double sx, double sigmax, double Rx){
	N = Nx;
	M = Mx;
	R = Rx;
	mu = mux;
	sigma = sigmax;
	k = kx;
	s = sx;
	// init random nnumber generators
	rngm = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (rngm, GenSeed());
	srand(time(NULL));
	// initialize population
	Cell c(3,0);
	int m0;
	for (int i=0;i<N;i++){
		m0 = 0;
		c = Cell(3,0);
		c[0] = m0; c[1] = M; c[2] = 0;
		Population.push_back( c );
	}
	int gen = 0;
	double Fit = 0;
	// equilibrate
	for (int i=0; i<400; i++){
		Fit = AvgFitness();
		Mutation();
		MitoSelection();
		Selection();
		Reproduction();
		gen += 1;
	}
	// insert mutants
	int nmuts = 0;
	while (nmuts < N*0.05){
		int r = gsl_rng_uniform_int(rngm, N);
		if (Population[r][2]==0){
			Population[r][2] = 1;
			nmuts += 1;
		}
	}
	retval = 0;
	// run until either fixation or extinction
	while (true){
		Fit = AvgFitness();
		Mutation();
		MitoSelection();
		Selection();
		Reproduction();
		gen += 1;
		double sf = GetSexFreq();
		if (sf==0 or sf==1) {
			retval = sf;
			break;
		}
	}
}

double Moose::GetFitness(Cell cx){
	assert(cx.size()==3);
	double Fit = 0;
	Fit = 1.0 - sigma*pow(1.0*cx[0]/cx[1],s);
	return Fit;
}

void Moose::Reproduction(){
	vector<Cell> SubPop1;
	vector<Cell> SubPop2;
	random_shuffle(Population.begin(),Population.end());
	SubPop1 = vector<Cell> (Population.begin(), Population.begin() + N/2);
	SubPop2 = vector<Cell> (Population.begin() + N/2, Population.end());
	assert(SubPop1.size() == N/2); assert(SubPop2.size() == N/2);
	vector<Cell> NewPopulation;
	for (int i=0;i<N/2;i++){
		vector<Cell> cells;
		if ((SubPop1[i][2]==1 and gsl_rng_uniform(rngm)<R) or (SubPop2[i][2]==1 and gsl_rng_uniform(rngm)<R)){
			cells = Fuse(SubPop1[i],SubPop2[i]);
			cells = ASexFuse(cells[0],cells[1]);
		}
		else {
			cells = ASexFuse(SubPop1[i],SubPop2[i]);
		}
		NewPopulation.push_back(cells[0]);
		NewPopulation.push_back(cells[1]);
	}
	Population = NewPopulation;
}

void Moose::Mutation(){
	for (int i=0;i<N;i++){
		int m = Population[i][0];
		int M = Population[i][1];
		int dm = gsl_ran_binomial(rngm, mu, M-m);
		Population[i][0] += dm;
	}
}

void Moose::Selection(){
	vector<double> fits(0);
	for (int i=0;i<N;i++){
		fits.push_back(GetFitness(Population[i]));
	}
	vector <double> sample = Sample(fits, N);
	vector<Cell> newPopulation;
	for (int i=0;i<N;i++){
		newPopulation.push_back(Population[(int)sample[i]]);
	}
	Population = newPopulation;
}

void Moose::MitoSelection(){ // selection on the lower level
	for (int i=0;i<N;i++){
		int m = Population[i][0];
		int M = Population[i][1];
		int mnew = gsl_ran_binomial(rngm, m*(1+k)/(M+m*k),M);
		Population[i][0] = mnew;
	}
}

vector<double> Moose::Sample(vector<double> list, int n){
	// sample n items from the list of weights with replacement
	double total = accumulate(list.begin(),list.end(), 0.0);
	vector<double> selection;
	int j = 0;
	double w = list[0];
	while (n>0){
		double r = gsl_rng_uniform (rngm);
		double x = total*(1.0 - pow(r, 1.0/n));
		total = total - x;
		while (x > w){
			x = x-w;
			j+=1;
			w = list[j];
		}
		w-=x;
		n-=1;
		selection.push_back(j);
	}
	return selection;
}

vector<Cell> Moose::Fuse(Cell cx1, Cell cx2){
	int mn = 0;
	int Mn = 0;
	mn = cx1[0] + cx2[0];
	Mn = cx1[1] + cx2[1];
	vector<int> FL;
	FL = vector<int> (2,0);
	FL[0] = cx1[2]; FL[1] = cx2[2];
	assert(Mn == 2*M);
	int m1;
	int m2;
	m1 = gsl_ran_hypergeometric (rngm, mn, Mn-mn, M);
	m2 = mn - m1;
	Cell c1(3,0);
	Cell c2(3,0);
	random_shuffle(FL.begin(), FL.end());
	c1[0] = m1; c1[1] = M; c1[2] = FL[0];
	c2[0] = m2; c2[1] = M; c2[2] = FL[1];
	vector<Cell> cells;
	cells.push_back(c1);
	cells.push_back(c2);
	return cells;
}

vector<Cell> Moose::ASexFuse(Cell cx1, Cell cx2){
	// no fusion actually occurs
	int mn1 = 0; int mn2 = 0;
	int Mn1 = 0; int Mn2 = 0;
	mn1 = cx1[0]; mn2 = cx2[0];
	Mn1 = cx1[1]; Mn2 = cx2[1];
	assert(Mn1 == M); assert(Mn2 == M);
	int m1;
	int m2;
	m1 = gsl_ran_hypergeometric (rngm, 2*mn1, 2*Mn1-2*mn1, M);
	m2 = gsl_ran_hypergeometric (rngm, 2*mn2, 2*Mn2-2*mn2, M);
	Cell c1(3,0);
	Cell c2(3,0);
	c1[0] = m1; c1[1] = M; c1[2] = cx1[2];
	c2[0] = m2; c2[1] = M; c2[2] = cx2[2];
	vector<Cell> cells;
	cells.push_back(c1);
	cells.push_back(c2);
	return cells;
}

double Moose::GetSexFreq(){
	double s = 0;
	for (int i=0; i<N; i++){
		s += Population[i][2];
	}
	return s/N;
}

double Moose::AvgFitness(){
	double fits = 0;
	for (int i=0; i<N; i++){
		fits += GetFitness(Population[i]);
	}
	return fits/N;
}

void Moose::PrintDistribution(int FG){
	vector<int> Mutvec;
	Mutvec = vector<int> (M+1,0);
	assert((int)Mutvec.size()==M+1);
	for (int i=0; i<N; i++){
		if (Population[i][2]==FG){
			Mutvec[Population[i][0]] += 1;
		}
	}
	int mtot = accumulate(Mutvec.begin(), Mutvec.end(), 0);
	for (int i=0; i<M+1; i++){
		cout << i << " " << 1.0*Mutvec[i]/mtot << endl;
	}
}

long Moose::GenSeed(){ // generate seed based on time and process id (for parallel comp)
	long s, seed, pid;
	int sr;
	time_t seconds;
	pid = getpid();
	s = time(&seconds);
	seed = abs(((s*181)*((pid-83)*359))%104729);
	ifstream file("/dev/urandom", std::ios::binary);
	if (file.is_open()){
	char * memblock;
	int size = sizeof(int);
	memblock = new char [size];
	file.read(memblock, size);
	file.close();
	sr = *reinterpret_cast<int*>(memblock);
	delete[] memblock;
	}
	else{sr = 1 ;}
	seed += sr;
	return seed;
}
