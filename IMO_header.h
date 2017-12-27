// Header file for DME dihedral fitter
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <time.h> 
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <random>
#include <vector>
#include <algorithm>

using std::cout;
using std::endl;
using std::string;
using std::array;
using std::vector;
using std::ifstream;
using std::ifstream;
using std::ofstream;

#define TYPE_INT 0
#define TYPE_FLOAT 1
#define TYPE_INT_3VEC 2
#define TYPE_FLOAT_3VEC 3
#define TYPE_STRING 4

#define WORD_STRING_SIZE 64

typedef struct {
	
	char keyword[WORD_STRING_SIZE];
	int dataType;
	void* varPtr;
	char defString[WORD_STRING_SIZE];
	
} input_data_struct;

struct fitting_param_struct {
	
	vector<int> dihedralIndexMapping;
	vector<int> dihedralMappingIndex;
	vector<int> rbTypeCount; 
	vector<double> rbCoeff;
	
	vector<int> pairIndexMapping;	
	vector<int> ljTypeCount;
	vector<double> ljSigma;
	vector<double> ljEps;
	
	vector<int> surfIndexMapping;
	vector<double> uShift;
	
	void zero_params() {
		std::fill(uShift.begin(), uShift.end(), 0.0);
		std::fill(rbCoeff.begin(), rbCoeff.end(), 0.0);
		std::fill(ljEps.begin(), ljEps.end(), 0.0);
		std::fill(ljSigma.begin(), ljSigma.end(), 0.0);
	}
};

struct simplex_struct {
	
	vector<fitting_param_struct> plex;
	vector<double> plexErrors;	
	vector<int> errorRankToRow;
	
	fitting_param_struct centroid;
	fitting_param_struct reflectPoint;
	fitting_param_struct expandPoint;
	fitting_param_struct contraPoint;
	
	int numSurfs;
	int numDihedralParams;
	int numPairs;
	
	// Nelder-Mead parameters
	double nmReflect;
	double nmExpand;
	double nmContra;
	double nmShrink;
	
	void plexInitialize() {
		
		const double factorLJ = 0.01;
		const double dihedralStep = 0.5;
		
		centroid = plex[0];
		reflectPoint = plex[0];
		expandPoint = plex[0];
		contraPoint = plex[0];
		
		plexErrors.resize(plex.size(),0.0);
		errorRankToRow.resize(plex.size());
		
		bool sigEps=1;
		int pair=0;
				
		numSurfs = plex[0].uShift.size();
		numDihedralParams = plex[0].rbCoeff.size();
		numPairs = plex[0].ljSigma.size();	
		
		// Partition of simplex
		int p1 = numSurfs;
		int p2 = p1 + numDihedralParams;
		int p3 = p2 + numPairs*2;	
		
		for (int row=0; row<p1; row++) {
			plex[row+1].uShift[row] += dihedralStep;	
		}
		for (int row=p1; row<p2; row++) {
			plex[row+1].rbCoeff[row-p1] += dihedralStep;	
		}
		for (int row=p2; row<p3; row++) {
			if (sigEps) {
				plex[row+1].ljSigma[pair] *= (1.0+factorLJ);
			}
			else {
				plex[row+1].ljEps[pair++] *= (1.0-factorLJ);
			}	
			sigEps = !sigEps;
		} 
	}
	
	// Calc centroid
	void computeCentroid() {
		
		vector<double> shiftAvg(numSurfs,0.0);
		vector<double> dihedralAvg(numDihedralParams,0.0);
		vector<double> sigmaAvg(numPairs,0.0);
		vector<double> epsilonAvg(numPairs,0.0);
		
		for (int row=0; row<plex.size(); row++) {
			
			if (row != errorRankToRow[plex.size()-1]) {
				
				for (int i=0; i < numSurfs; i++) {
					shiftAvg[i] += plex[row].uShift[i];
				}
				for (int j=0; j < numDihedralParams; j++) {
					dihedralAvg[j] += plex[row].rbCoeff[j];
				}
				for (int k=0; k < numPairs; k++) {
					sigmaAvg[k] += plex[row].ljSigma[k];				
					epsilonAvg[k] += plex[row].ljEps[k];	
				}
			}			
		}
		for (int i=0; i < numSurfs; i++) {
			centroid.uShift[i] = shiftAvg[i]/(plex.size()-1);
		}
		for (int j=0; j < numDihedralParams; j++) {
			centroid.rbCoeff[j] = dihedralAvg[j]/(plex.size()-1);
		}
		for (int k=0; k < numPairs; k++) {
			centroid.ljSigma[k] = sigmaAvg[k]/(plex.size()-1);				
			centroid.ljEps[k] = epsilonAvg[k]/(plex.size()-1);	
		}
	} 
	
	// Calc reflect point
	void computeReflect() {
		//
		int worstRow = errorRankToRow[plex.size()-1];
			
		for (int i=0; i < numSurfs; i++) {
			reflectPoint.uShift[i] = centroid.uShift[i] + nmReflect*(centroid.uShift[i] - plex[worstRow].uShift[i]);
		}
		for (int j=0; j < numDihedralParams; j++) {
			reflectPoint.rbCoeff[j] = centroid.rbCoeff[j] + nmReflect*(centroid.rbCoeff[j] - plex[worstRow].rbCoeff[j]);
		}
		for (int k=0; k < numPairs; k++) {
			reflectPoint.ljSigma[k] = centroid.ljSigma[k] + nmReflect*(centroid.ljSigma[k] - plex[worstRow].ljSigma[k]);	
			reflectPoint.ljEps[k] = centroid.ljEps[k] + nmReflect*(centroid.ljEps[k] - plex[worstRow].ljEps[k]);	
		}
	}
	
	// Calc expanded point
	void computeExpand() {
		//			
		for (int i=0; i < numSurfs; i++) {
			expandPoint.uShift[i] = reflectPoint.uShift[i] + nmExpand*(reflectPoint.uShift[i] - centroid.uShift[i]);
		}
		for (int j=0; j < numDihedralParams; j++) {
			expandPoint.rbCoeff[j] = reflectPoint.rbCoeff[j] + nmExpand*(reflectPoint.rbCoeff[j] - centroid.rbCoeff[j]);
		}
		for (int k=0; k < numPairs; k++) {
			expandPoint.ljSigma[k] = reflectPoint.ljSigma[k] + nmExpand*(reflectPoint.ljSigma[k] - centroid.ljSigma[k]);	
			expandPoint.ljEps[k] = reflectPoint.ljEps[k] + nmExpand*(reflectPoint.ljEps[k] - centroid.ljEps[k] );	
		}
	}
		
	// Calc contracted point
	void computeContract() {
		//
		int worstRow = errorRankToRow[plex.size()-1];
			
		for (int i=0; i < numSurfs; i++) {
			contraPoint.uShift[i] = centroid.uShift[i] + nmContra*(plex[worstRow].uShift[i] - centroid.uShift[i]);
		}
		for (int j=0; j < numDihedralParams; j++) {
			contraPoint.rbCoeff[j] = centroid.rbCoeff[j] + nmContra*(plex[worstRow].rbCoeff[j] - centroid.rbCoeff[j]);
		}
		for (int k=0; k < numPairs; k++) {
			contraPoint.ljSigma[k] = centroid.ljSigma[k] + nmContra*(plex[worstRow].ljSigma[k] - centroid.ljSigma[k]);	
			contraPoint.ljEps[k] = centroid.ljEps[k] + nmContra*(plex[worstRow].ljEps[k] - centroid.ljEps[k]);	
		}
	}
	
	// Reduce simplex
	void reducePoint(int p) {
		//
		int bestRow = errorRankToRow[0];
			
		for (int i=0; i < numSurfs; i++) {
			plex[p].uShift[i] = plex[p].uShift[i]*nmShrink + (1.0-nmShrink)*plex[bestRow].uShift[i];
		}
		for (int j=0; j < numDihedralParams; j++) {
			plex[p].rbCoeff[j] = plex[p].rbCoeff[j]*nmShrink + (1.0-nmShrink)*plex[bestRow].rbCoeff[j];
		}
		for (int k=0; k < numPairs; k++) {
			plex[p].ljSigma[k] = plex[p].ljSigma[k]*nmShrink + (1.0-nmShrink)*plex[bestRow].ljSigma[k];
			plex[p].ljEps[k] = plex[p].ljEps[k]*nmShrink + (1.0-nmShrink)*plex[bestRow].ljEps[k];
		}
	}
	
};

struct molsize_struct {
	
	int atoms = 0;
	int bonds = 0;
	int angles = 0;
	int dihedrals = 0;
	int impDihedrals = 0;
	int pairs = 0;
};

struct constant_energy_struct {
	
	double uAngle = 0.0;
	double uDihedral = 0.0;
	double uImpDihedral = 0.0;
	double uLJ = 0.0;
	double uQQ = 0.0;
	double uLJ14 = 0.0;
	double uQQ14 = 0.0;
	double u15 = 0.0;
	double u16plus = 0.0;
	double uTotal = 0.0;
	
	void sumEnergies() {
		uTotal = uAngle+uDihedral+uImpDihedral+uLJ+uQQ;
	}
};

struct constant_struct {
	
	long double pi;
	
	bool xyzAngstroms;
	bool genXyzFileNames;
	string inputFileString;
	string parameterFile;
	string phiPairFiles;
	string numPhiPairString;
	
	int writeEnergyDist;
	
	int nrexcl;
	int gromacsCombRule;
	double scale14LJ;
	double scale14QQ;
	double sigma14factor;
	double sigma14factorCutoff;
	
	double epsQQ;
	int nRBfit;
	int numConfigs;
  
	int phiDim;
	double phi1Range[3];
	double phi2Range[3];
	double phi3Range[3];
	int numPhiSurfaces;
	int useNWsuffix;	
	int weightPhiEdges;

	int annealWrite;
	int downhillWrite;
	int atomDataSize;
	int bondDataSize;
	int angleDataSize;
	int dihedralDataSize;
	int improperDataSize;
	int pairDataSize;

	int numDihedralParams;
	int numTotalParams;
  
	bool resToZero;
	double epsK;
	double rbK;
	
	bool useBoltzIntRes;
	double kTBoltzIntegral;
	double kBoltzRes;
	
	int useAdaptiveNelderMead;
  
	double vTempInitial;
	double vTempFinal;
	double vTempFactor;
	double dihedralStep;
	double sigmaStep;
	double epsilonStep;
	double KT;
  
	double nmReflect;
	double nmExpand;
	double nmContra;
	double nmShrink;
  
	double sigmaGradFactor;
};

struct vector_struct {
	
	vector<double> phi1partition;
	vector<double> phi2partition;
	
	vector<string> xyzFileList;
	vector<string> energyFileList;
	vector<string> connectFileList;
	vector<string> phiPairFileList;
	vector<int> numPhiPairs;
	vector<string> sortedConFiles;
	
	vector<int> partitionMap;
	vector<int> integrationRule;
	vector<double> simpsonCoeffsPhi1;
	vector<double> simpsonCoeffsPhi2;
	vector<double> confIntegralsDFT;
	
	vector<vector<double> > xyzData;
	vector<vector<double> > dftData;
	vector<vector<double> > energyWeighting;

	vector<vector<constant_energy_struct> > uConst;
	
	vector<molsize_struct> molSize;
	
	vector<vector<double> > atomData;
	vector<vector<double> > bondData;
	vector<vector<double> > angleData;
	vector<vector<double> > dihedralData;
	vector<vector<double> > improperData;
	vector<vector<double> > pairData;
	
	vector<vector<double> > cosPsiData;
	vector<vector<double> > pairSepData;
	
	vector<vector<double> > energyDelta;
	
	vector<vector<int>    > bondSepMat;
	vector<vector<double> > rijMatrix;
	vector<vector<double> > sigmaMatrix;
	vector<vector<double> > epsilonMatrix;
	vector<vector<double> > qqMatrix;
	
};

// Function prototypes
int read_input_params(constant_struct &cons, vector_struct &vecs, fitting_param_struct &initialParams);
int process_input_line(string fLine, input_data_struct* inputDefaults, int inputDefaultSize, bool toPrint);

int xyz_files_read(constant_struct &cons, vector_struct &vecs, int i_s);
int connectivity_read(constant_struct &cons, vector_struct &vecs, int i_s);
int connectivity_process(constant_struct cons, vector_struct &vecs, int i_s);
int energy_read(constant_struct &cons, vector_struct &vecs, int i_s);
int constant_energy_process(constant_struct cons, vector_struct &vecs, fitting_param_struct &initialParams, int i_s);

double error_from_trial_point(constant_struct cons, vector_struct &vecs, fitting_param_struct initialParams, fitting_param_struct trialParams, bool toWrite);
int compute_boltzmann_integrals(constant_struct cons, vector_struct vecs, vector<double> energyTotal, vector<double> &confIntegralsMD, bool toPrint);
int compute_gradient_F(constant_struct cons, vector_struct vecs, fitting_param_struct initialParams, fitting_param_struct currentParams, fitting_param_struct &grad);

int simulated_annealing(constant_struct cons, vector_struct vecs, fitting_param_struct &initialParams, fitting_param_struct &currentParams, int annealIt);
int downhill_simplex(constant_struct cons, vector_struct vecs, fitting_param_struct &initialParams, fitting_param_struct &currentParams, int simplexIt);
int error_sort(vector<double> simplexErrors, vector<int> &errorRankToRow);
//int steepest_descent(constant_struct cons, vector_struct vecs, vector<double> initialParams, vector<double> &currentParams, int gradientIt);
int conjugate_gradient(constant_struct cons, vector_struct vecs, fitting_param_struct initialParams, fitting_param_struct &currentParams, int gradientIt);

int param_linear_combine(fitting_param_struct &outParams, fitting_param_struct aParams, fitting_param_struct bParams, double a, double b);
double param_scalar_product(fitting_param_struct aParams, fitting_param_struct bParams);

int print_params_console(constant_struct &cons, fitting_param_struct &printParams);
int print_grad_console(constant_struct &cons, fitting_param_struct &printParams);
int print_simplex_console(constant_struct &cons, simplex_struct &printSim);
