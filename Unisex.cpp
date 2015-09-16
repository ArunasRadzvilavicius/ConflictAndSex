#include <iostream>
#include <cstdio>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <assert.h>
#include "Moose.h"
using namespace std;

int main(int argc, char **argv)
{	
	if (argc!=8) {
		cout << "unisex N Mx mux kx sx sigmax Rx" << endl;
		return 1;
	}
	int N = atoi(argv[1]);			// Population size
	int Mx = atoi(argv[2]);			// Number of edosymbionts per cell
	double mux = atof(argv[3]);		// Mutation rate cooperative->selfish
	double kx = atof(argv[4]);		// Selfish reproductive advantage
	double sx = atof(argv[5]);		// Power in the fitness function w=1-(m/M)^s
	double sigmax = atof(argv[6]);	// Strength of selection -- not used
	double Rx = atof(argv[7]);		// Fusion rate of sexual individuals
	int nfix = 0;
	int ntot = 0;
	Moose moosy;
	for (int i=0;i<200;i++){
		moosy = Moose(N,Mx,mux,kx,sx,sigmax,Rx);
		nfix += moosy.retval;
		ntot += 1;		
	}
	cout << nfix << ' ' << ntot << endl;
	return 0;
}
