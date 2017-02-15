#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <time.h> 
#include <string>
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
using std::vector;
using std::ifstream;
using std::ofstream;

#include "IMO_header.h"

int read_input_params(constant_struct &cons, vector_struct &vecs)
{
	cons.inputFileString = "input_files_DME.txt";
	cons.pi = acosl(-1.0L);
		
	cons.xyzAngstroms = true;
	
	cons.nrexcl = 3;
	cons.gromacsCombRule = 2;
	cons.scale14LJ = 0.5;
	cons.scale14QQ = 0.5;
	cons.sigma14factor = 1.0; 			// Multiply sigma in 1-4 potentials by this factor
	cons.sigma14factorCutoff = 0.3;		// Only apply sigma14factor to sigma above this cutoff
	
	cons.nRBfit = 6; // Highest power term is cos^(nRBfit-1)
	
	cons.phiDim = 2;
	
	int numPhi1 = 0, numPhi2 = 0, numPhi3 = 0;
	
	if (cons.phiDim == 2)
	{
		cons.phi1min = 0.0;
		cons.phi1max = 180.0;
		cons.phi1step = 10.0;
		cons.phi2min = 0.0;
		cons.phi2max = 360.0;
		cons.phi2step = 10.0;
		
		numPhi1 = round((cons.phi1max-cons.phi1min)/cons.phi1step + 1);
		numPhi2 = round((cons.phi2max-cons.phi2min)/cons.phi2step + 1);
		cons.numConfigs = numPhi1*numPhi2;		
	}
	else if (cons.phiDim == 3)	
	{
		cons.phi1min = 0.0;
		cons.phi1max = 180.0;
		cons.phi1step = 20.0;
		cons.phi2min = 0.0;
		cons.phi2max = 340.0;
		cons.phi2step = 20.0;
		cons.phi3min = 0.0;
		cons.phi3max = 340.0;
		cons.phi3step = 20.0;	
		
		numPhi1 = round((cons.phi1max-cons.phi1min)/cons.phi1step + 1);
		numPhi2 = round((cons.phi2max-cons.phi2min)/cons.phi2step + 1);
		numPhi3 = round((cons.phi3max-cons.phi3min)/cons.phi3step + 1);
		cons.numConfigs = numPhi1*numPhi2*numPhi3;
	}
	
	vecs.simpsonCoeffsPhi1.resize(numPhi1);
	vecs.simpsonCoeffsPhi2.resize(numPhi2);
	std::fill(vecs.simpsonCoeffsPhi1.begin(),vecs.simpsonCoeffsPhi1.begin()+numPhi1,1);
	std::fill(vecs.simpsonCoeffsPhi2.begin(),vecs.simpsonCoeffsPhi2.begin()+numPhi2,1);
	
	// Partitioning of potential energy surface into conformers
	if (cons.inputFileString == "input_files_DME.txt") {

		vecs.phi1partition.push_back(0.0);
		vecs.phi1partition.push_back(120.0); 
		vecs.phi1partition.push_back(180.0);
		
		vecs.phi2partition.push_back(0.0);
		vecs.phi2partition.push_back(120.0); 
		vecs.phi2partition.push_back(240.0);
		vecs.phi2partition.push_back(360.0);
		
		// Sections are numbered in row major order, map each to a conformer
		int partitionMapDME[] = {4, 2, 5, 3, 1, 3}; // TGG, TGT, TGG', TTG, TTT, TTG'
		vecs.partitionMap.assign(partitionMapDME,partitionMapDME+6);
		
		// Integration rule for each conformer section
		vecs.integrationRule.resize((vecs.phi1partition.size()-1)*(vecs.phi2partition.size()-1));
		
		int i_section = 0;		
		for (int c1=1; c1<vecs.phi1partition.size(); c1++)
		{
			for (int c2=1; c2<vecs.phi2partition.size(); c2++)
			{
				
				int m_max = round(vecs.phi1partition[c1]/cons.phi1step);
				int m_min = round(vecs.phi1partition[c1-1]/cons.phi1step);
				int n_max = round(vecs.phi2partition[c2]/cons.phi2step);
				int n_min = round(vecs.phi2partition[c2-1]/cons.phi2step);
				
				int m = m_max-m_min;
				int n = n_max-n_min;
				
				// Determine if Simpson's rule integrable
				if ( (m%2 == 0) && (n%2 == 0)) {
					vecs.integrationRule[i_section] = 1;					
				}
				else {
					vecs.integrationRule[i_section] = 0; // Trapezoid rule
				}
				
				// Populate array of Simpsons and Trapezoid rule coefficients
				for (int i_m = m_min+1; i_m < m_max; i_m++)
				{
					if ((i_m-m_min)%2 == 0) {
						vecs.simpsonCoeffsPhi1[i_m] = 2;
					}
					else {
						vecs.simpsonCoeffsPhi1[i_m] = 4;
					}
					
				}
				for (int i_n = n_min+1; i_n < n_max; i_n++)
				{
					if ((i_n-n_min)%2 == 0) {
						vecs.simpsonCoeffsPhi2[i_n] = 2;
					}
					else {
						vecs.simpsonCoeffsPhi2[i_n] = 4;
					}			
				}
			}
		}
		for (int i=0; i<numPhi1; i++)
			cout << endl << vecs.simpsonCoeffsPhi1[i];
	}
	
	// Constraints and restraints
	cons.useBoltzIntRes = 1;
	cons.kTBoltzIntegral = 2.5;
	
	// Restraining params to default values
	cons.resToZero = 0;
	cons.dihedralK = 0.0;
	cons.epsK = 10.0; // Constraint epsK*(eps-eps_fit)^2 per configuration
	
	//cons.KT = 5.23; // 20.92/4
	//cons.KT = 15.0; // Typical max saddle point energy
	cons.KT = 20.92; // 5kcal/mol (Amber)
	//cons.KT = 83.14; // 10,000K (TraPPE)
	
	
	// Molecule specific parameters (overwrite numConfigs if needed)
	if (cons.inputFileString == "input_files_EC.txt") {
		cons.resToZero = 1;
		cons.numConfigs = 869;
		cons.dihedralK = 1.0E-6;
	}
	if (cons.inputFileString == "input_files_PC.txt") {
		cons.resToZero = 1;
		cons.numConfigs = 2976;
		cons.dihedralK = 1.0E-6;
	}
	
	cons.annealWrite = 100;	// Frequency of writing lowest error to file
	cons.downhillWrite = 10;
	
	cons.atomDataSize = 4;
	cons.bondDataSize = 4;
	cons.angleDataSize = 5;
	cons.dihedralDataSize = 11;
	cons.improperDataSize = 7;
	cons.pairDataSize = 4;
	
	cons.numDihedralFit = 0;
	
	// Simulated annealing parameters
	cons.vTempInitial = (5000 + 50*cons.KT)*double(cons.numConfigs)/1000.0;
	cons.vTempFactor = 0.9998;
	cons.dihedralStep = 5.0;
	cons.sigmaStep = 0.002;
	cons.epsilonStep = 0.001;
	
	// Nelder-Mead parameters
	cons.nmReflect = 1.0;
	cons.nmExpand = 2.0;
	cons.nmContract = 0.8;
	cons.nmShrink = 0.5;
	
	cons.sigmaGradFactor = 0.01;
	
	
	// Get filenames from inputFileString
	ifstream inputNames;
	cout << "Opening input file " << cons.inputFileString << endl;
	inputNames.open(cons.inputFileString.c_str());

	inputNames >> cons.xyzFile >> cons.energyFile >> cons.connectFile >> cons.genXyzFileNames;
	
	if (!cons.genXyzFileNames)
		inputNames >> cons.phiCoordFile;	
	
	// Read number of atoms, bonds etc. from first section of input file
	ifstream connectStream;
	connectStream.open(cons.connectFile.c_str());

	string fileHeader;
	connectStream >> fileHeader;
	assert(fileHeader == "size");
	
	// Read connectivity size
	connectStream >> cons.size[0] >> cons.size[1] >> cons.size[2] >> cons.size[3] >> cons.size[4] >> cons.size[5];
	
	// Set size of vectors
	vecs.xyzData.resize(cons.numConfigs*cons.size[0]*3);	
	vecs.energyData.resize(cons.numConfigs);
	vecs.energyWeighting.resize(cons.numConfigs);
	vecs.constantEnergy.resize(cons.numConfigs);
	
	vecs.atomData.resize(cons.size[0]*cons.atomDataSize);
	vecs.bondData.resize(cons.size[1]*cons.bondDataSize);
	vecs.angleData.resize(cons.size[2]*cons.angleDataSize);
	vecs.dihedralData.resize(cons.size[3]*cons.dihedralDataSize);
	
	if (cons.size[4] != 0)
		vecs.improperData.resize(cons.size[4]*cons.dihedralDataSize);
	if (cons.size[5] != 0)
	{
		vecs.ljPairFitData.resize(cons.size[5]*cons.pairDataSize);
		vecs.pairSepData.resize(cons.numConfigs*cons.size[5]);
	}
	
	vecs.bondSepMat.resize(cons.size[0]*cons.size[0]);
	vecs.rijMatrix.resize(cons.size[0]*cons.size[0]);
	vecs.sigmaMatrix.resize(cons.size[0]*cons.size[0]);
	vecs.epsilonMatrix.resize(cons.size[0]*cons.size[0]);
	vecs.qqMatrix.resize(cons.size[0]*cons.size[0]);
	
	vecs.phiPairSpec.resize(cons.phiDim*cons.numConfigs);
	
	return 0;
}

int main(int argc, char *argv[])
{
	int annealIt, simplexIt, downhillIt; // Number of iterations for sim. annealing, downhill simplex & conj. grad.
	
	if (argc > 3)
	{
		annealIt = atof(argv[1]);
		simplexIt = atof(argv[2]);
		downhillIt = atof(argv[3]);
	}
	else
	{
		cout << "Not enough command line arguements!\n\n";
		return 1;	
	}
	
	// Initialise data structures
	constant_struct cons; // Constants 
	vector_struct vecs;   // Vectors
		
	read_input_params(cons, vecs);
	
	srand(time(NULL));
	
	// Read xyz files
	xyz_files_read(cons, vecs);
		
	// Read connectivity data
	connectivity_read(cons, vecs);
	
	// Process connectivity data
	connectivity_process(cons, vecs);
	
	// Read energy file
	int energyRead = energy_read(cons, vecs);	
	if (energyRead == 1)
		return 1;
	
	// Deterine number of dihedrals to fit, and the starting parameters
	for (int d = 0; d < cons.size[3]; d++)
	{
		if (vecs.dihedralData[d*cons.dihedralDataSize + 4] > 0)
			cons.numDihedralFit++;
	}
	// Determine number of pairs to fit (future feature)
	
	cons.numDihedralParams = 2*(cons.nRBfit-1) + 1;	// +1 is the constant energy term
		
	if (cons.size[5] == 0)
		cons.numTotalParams = cons.numDihedralParams;
	else
		cons.numTotalParams = cons.numDihedralParams + 2;
		
	vector<double> initialParams(cons.numTotalParams);
	vector<double> currentParams(cons.numTotalParams);
	
	// Set adaptive Nelder-Mead parameters
	cons.nmReflect = 1.0;
	cons.nmExpand = 1.0 + 2.0/double(cons.numTotalParams);
	cons.nmContract = 0.9 - 0.8/double(cons.numTotalParams);
	cons.nmShrink = 1.0 - 1.0/double(cons.numTotalParams);	
	
	cons.simplexSize = cons.numTotalParams + 1;
	vecs.cosPsiData.resize(cons.numConfigs*cons.numDihedralFit);
	
	// Process constant terms in energy
	int constSuccess = constant_energy_process(cons, vecs);
	if (constSuccess == 1)
		return 1;

	// Setup vector of initial parameters
	// Loop over all dihedrals
	for (int d = 0; d < cons.size[3]; d++)
	{
		// Get dihedral type
		int dType = int(round(vecs.dihedralData[d*cons.dihedralDataSize + 4]));
		
		if (vecs.dihedralData[d*cons.dihedralDataSize + 4] != 0) // If being used in fit
		{
			// Loop over coefficients
			for (int rb=0; rb<cons.nRBfit; rb++)
			{				
				// Get coefficient
				if (!((dType > 1) && (rb == 0))) // Exclude constant terms apart from first
				{
					initialParams[(dType-1)*(cons.nRBfit-1) + rb] = vecs.dihedralData[d*cons.dihedralDataSize + 5 + rb];
					currentParams[(dType-1)*(cons.nRBfit-1) + rb] = vecs.dihedralData[d*cons.dihedralDataSize + 5 + rb];
				}
			}	
		}
	}
	if (cons.size[5] != 0)
	{
		initialParams[cons.numDihedralParams    ] = vecs.ljPairFitData[2]; // Lennard-Jones sigma and epsilon
		initialParams[cons.numDihedralParams + 1] = vecs.ljPairFitData[3];	
		currentParams[cons.numDihedralParams    ] = vecs.ljPairFitData[2];
		currentParams[cons.numDihedralParams + 1] = vecs.ljPairFitData[3];
	}
	
	
	// -------- DO OPTIMISATION ------------------------------------------------------------------------------- 
	
	simulated_annealing(cons, vecs, initialParams, currentParams, annealIt);
	
	downhill_simplex(cons, vecs, initialParams, currentParams, simplexIt);
	
	//steepest_descent(cons, vecs, initialParams, currentParams, downhillIt);
	
	conjugate_gradient(cons, vecs, initialParams, currentParams, downhillIt);
	
	// ---------------------------------------------------------------------------------------------------------
	
	// Write best energy
	error_from_trial_point(cons, initialParams, currentParams, 0, vecs, 1);
	
}

int simulated_annealing(constant_struct cons, vector_struct vecs, vector<double> initialParams, vector<double> &currentParams, int annealIt)
{
	int printFreq = 100; // Print simplex to console every printFreq steps
	
	srand(time(NULL));
	double temperature = cons.vTempInitial;
	
	ofstream annealStream;
	annealStream.open("annealing_energy.txt");
	annealStream.precision(12);
	bool toPrint, toWrite;
	
	vector<double> oldParams(cons.numTotalParams);
	vector<double> annealParams(cons.numTotalParams);
		
	for (int p=0; p<(cons.numTotalParams); p++)
	{
		oldParams[p] = currentParams[p];
		cout << oldParams[p] << " ";
	}
	cout << endl;
	
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0,1.0);
	
	for (int iT=0; iT < annealIt; iT++)
	{
		toPrint = (iT%printFreq == 0);
		toWrite = (iT%cons.annealWrite == 0);
		
		if (toPrint)
			cout << "Simulated annealing iteration " << iT << ", temp " << temperature << "\n\n";
		
		// Generate trial point
		for (int p=0; p<cons.numDihedralParams; p++)
		{
			double normalStep = distribution(generator);
			annealParams[p] = oldParams[p] + normalStep*sqrt(temperature/cons.vTempInitial)*cons.dihedralStep;
		}
		if (cons.size[5] != 0) // If fitting a LJ pair as well
		{
			// Use abs() to prevent negative LJ parameters, e.g. if step size is set too large
			annealParams[cons.numDihedralParams] = std::abs(oldParams[cons.numDihedralParams] 
				+ cons.sigmaStep*sqrt(temperature/cons.vTempInitial)*distribution(generator));
			annealParams[cons.numDihedralParams + 1] = std::abs(oldParams[cons.numDihedralParams + 1] 
				+ cons.epsilonStep*sqrt(temperature/cons.vTempInitial)*distribution(generator));	
		}
		
		// Test neighbour parameters
		long double oldError = error_from_trial_point(cons, initialParams, oldParams, 0, vecs, 0);
		annealStream << oldError << endl;
		long double trialError = error_from_trial_point(cons, initialParams, annealParams, 0, vecs, 0);
		
		// Compute acceptance probability
		double aProb = exp((oldError-trialError)/temperature);
		
		if (toPrint)
		{
			cout << "Old error   = " << oldError << endl;
			cout << "Trial error = " << trialError << endl;
			cout << "Acceptance probability = " << aProb << "\n\n";
		}	
		if (aProb > (rand()/double(RAND_MAX)))
		{
			if (toPrint)
				cout << "Accepting these params: \n";
		
			for (int col=0; col < (cons.numTotalParams); col++)
			{
				oldParams[col] = annealParams[col];
				
				if (toPrint)
					cout << std::setprecision(8) << std::setw(10) << oldParams[col] << " ";
			}
			if (toPrint)
				cout << endl;
		}
		else if (toPrint)
		{
			cout << "Rejecting these params! \n\n";
		}
		
		temperature *= cons.vTempFactor;
	}
	
	// Update final parameters
	for (int col=0; col<(cons.numTotalParams); col++)
		currentParams[col] = oldParams[col];
	
	return 0;
}

int downhill_simplex(constant_struct cons, vector_struct vecs, vector<double> initialParams, vector<double> &currentParams, int simplexIt)
{
	cout << "\nPerforming " << simplexIt << " downhill simplex iteraitons\n\n";
	long double nmEPS = 1.0E-14L;
	int printFreq = 100; // Print simplex to console every printFreq steps
	int numReflectSteps = 0;
	
	// Generate initial simplex	
	vector<double> simplex(cons.simplexSize*(cons.numTotalParams));
	vector<long double> simplexErrors(cons.simplexSize);
	vector<double> reflectParams(cons.numTotalParams);
	vector<double> expansionParams(cons.numTotalParams);
	vector<double> contractParams(cons.numTotalParams);
	vector<double> centroid(cons.numTotalParams);
	vector<int> errorRankToRow(cons.simplexSize); // Argument is error rank, output is corresponding row of simplex
	
	// Initialise error rank vector as 0,1,2...
	for (int r=0; r<cons.simplexSize; r++)
		errorRankToRow[r] = r;
	
	define_initial_simplex(cons, currentParams, simplex);
	
	bool toWrite, toPrint;
	ofstream simplexStream;
	simplexStream.precision(12);
	simplexStream.open("simplex_energy.txt");
	
	// Analyse simplex and fill simplexErrors
	for (int row=0; row<cons.simplexSize; row++)
	{
		simplexErrors[row] = error_from_trial_point(cons, initialParams, simplex, (row*cons.numTotalParams), vecs, 0);
		cout << "Row: " << row << ", error " << simplexErrors[row] << endl;
	}
	
	for (int iD=0; iD < simplexIt; iD++)
	{
		toWrite = (iD%cons.downhillWrite == 0);
		toPrint = (iD%printFreq == 0);
		
		if (toPrint)
			cout << "\nIteration: " << iD << endl << endl;

		// Sort errors, determine max, min and write best
		error_sort(cons.simplexSize, simplexErrors, errorRankToRow);
		
		int worstRow = errorRankToRow[cons.simplexSize-1];
		int maxPoint = worstRow*(cons.numTotalParams); // Start of row containing highest error point

		long double maxError = simplexErrors[ errorRankToRow[cons.simplexSize-1] ];
		long double secondMaxError = simplexErrors[ errorRankToRow[cons.simplexSize-2] ];
		long double minError = simplexErrors[ errorRankToRow[0] ];
		
		// Return if end condition met
		if (std::fabs(maxError-minError)/maxError < nmEPS)
		{
			cout << std::setprecision(12) << simplexErrors[worstRow] << endl;
			cout << std::setprecision(12) << simplexErrors[0] << endl;
			cout << std::setprecision(12) << std::fabs(maxError-minError)/maxError << endl;
			cout << "\nEnding downhill simplex routine after " << iD << " iterations\n\n";
			cout << "The percentage of reflect steps was " << 100.0*double(numReflectSteps)/double(iD) << endl;
			if (simplexIt > 0) // Update currentParams array
			{
				for (int col=0; col<(cons.numTotalParams); col++)
					currentParams[col] = simplex[ errorRankToRow[0]*cons.numTotalParams + col ];		
			}
			return 0;		
		}
		
		simplexStream << minError << endl;
		
		// Print best
		if (toPrint)
			cout << "MAX and MIN errors: " << std::setprecision(20) << maxError << " " << minError << endl;
		
		// Compute centroid of all points execpt worst
		for (int col=0; col<cons.numTotalParams; col++)
		{
			double sumCentroidCol = 0.0;
			for (int row=0; row<cons.simplexSize; row++)
			{
				// If this isn't worst row
				if (row != errorRankToRow[cons.simplexSize-1])
				{
					sumCentroidCol += simplex[row*(cons.numTotalParams) + col];
				}
			}
			centroid[col] = sumCentroidCol/(cons.simplexSize-1);
		}
		
		// Compute reflected point and output best params
		if (toPrint)
			cout << "\nComputing reflected point in opposite direction to point " << errorRankToRow[cons.simplexSize-1] << endl;
		
		for (int col=0; col<cons.numTotalParams; col++) 
		{
			reflectParams[col] = centroid[col] + cons.nmReflect*(centroid[col] - simplex[maxPoint + col]);
		}
		
		long double reflectError = error_from_trial_point(cons, initialParams, reflectParams, 0, vecs, 0);
		
		if (toPrint)
			cout << "\nTotal residual for reflected point is " << std::setprecision(14) << std::setw(16) << reflectError << endl << endl;
		
		// If neither new best nor worst (or 2nd worst), replace worst row of simplex with trial point
		if ((reflectError <= secondMaxError) && (reflectError >= minError))
		{
			if (toPrint)
			{
				cout << "Reflected point neither best nor worst.\n";
				cout << "Replacing worst point with reflected point.\n";	
			}			
			
			for (int col=0; col<cons.numTotalParams; col++)
			{
				simplex[maxPoint + col] = reflectParams[col];
			}
			simplexErrors[worstRow] = reflectError;
			numReflectSteps++;
		}
		// If best, test with expansion in direction of reflected point
		else if (reflectError < minError)
		{
			if (toPrint)
			{
				cout << "Reflected point is better than any point in simplex.\n";
				cout << "Generating new point by expansion.\n\n";
			}
			
			for (int col=0; col<cons.numTotalParams; col++)
			{
				expansionParams[col] = reflectParams[col] + cons.nmExpand*(reflectParams[col] - centroid[col]);
			}
			// Test expanded point
			long double expansionError = error_from_trial_point(cons, initialParams, expansionParams, 0, vecs, 0);
			
			if (expansionError < minError)
			{
				if (toPrint)
				{
					cout << "Expanded point even better.\n";
					cout << "Replacing worst point with expanded point.\n\n";
				}			
				for (int col=0; col < (cons.numTotalParams); col++)
				{
					simplex[maxPoint + col] = expansionParams[col];
				}
				simplexErrors[worstRow] = expansionError;
			}
			else
			{
				if (toPrint)
				{
					cout << "Expansion unsuccessful.\n";
					cout << "Replacing worst point with reflected point.\n\n";
				}
				for (int col=0; col <cons.numTotalParams; col++)
				{
					simplex[maxPoint + col] = reflectParams[col];
				}
				simplexErrors[worstRow] = reflectError;
				numReflectSteps++;
			
			}
		}
		// Else, the reflected point error must be > the second worst point
		else
		{
			if (toPrint)
				cout << "Reflection unsuccessful, performing contraction.\n";
			
			//Replace worst point by reflected point if it's an improvement
			if (reflectError < maxError)
			{
				for (int col=0; col<cons.numTotalParams; col++)
				{
					simplex[maxPoint + col] = reflectParams[col];
				}
				simplexErrors[worstRow] = reflectError;
				maxError = reflectError;
			}
			
			for (int col=0; col<cons.numTotalParams; col++)
			{
				contractParams[col] = centroid[col] + cons.nmContract*(simplex[maxPoint+col] - centroid[col]);
			}
			// Test contracted point
			long double contractError = error_from_trial_point(cons, initialParams, contractParams, 0, vecs, 0);
			
			if (contractError < maxError)
			{
				if (toPrint)
					cout << "Replacing worst point with contracted point.\n\n";
				
				for (int col=0; col<cons.numTotalParams; col++)
				{
					simplex[maxPoint + col] = contractParams[col];
				}	
				simplexErrors[worstRow] = contractError;
			}
			else // Reduction towards best point
			{
				if (toPrint)
					cout << "Contracted point is new worst point. Trying reduction towards best point.\n\n";
				
				for (int row=0; row<cons.simplexSize; row++)
				{
					for (int col=0; col<cons.numTotalParams; col++)
					{
						double bestParam = simplex[errorRankToRow[0]*(cons.numTotalParams) + col];
						double currentParam = simplex[row*(cons.numTotalParams) + col];
						simplex[row*cons.numTotalParams + col] = bestParam + cons.nmShrink*(currentParam-bestParam);
					}	
					// Update error for reduced point
					simplexErrors[row] = error_from_trial_point(cons, initialParams, simplex, (row*cons.numTotalParams), vecs, 0);
					
				}
			}	
		}
		
		// Output new simplex
		if (toPrint)
		{
			cout << "Updated simplex: " << endl;
			for (int row=0; row<cons.simplexSize; row++)
			{
				for (int col=0; col<(cons.numTotalParams); col++)
				{
					cout << std::setprecision(8) << std::setw(10) << simplex[row*(cons.numTotalParams) + col] << " ";
				}
				cout << endl;
			}
		}
	}
	
	if (simplexIt > 0) // Update currentParams array
	{
		for (int col=0; col<(cons.numTotalParams); col++)
			currentParams[col] = simplex[ errorRankToRow[0]*cons.numTotalParams + col ];		
	}
		
	cout << "\nEnding downhill simplex routine after max iterations\n\n";
	
	return 0;
}

int conjugate_gradient(constant_struct cons, vector_struct vecs, vector<double> initialParams, vector<double> &currentParams, int gradientIt)
{
	cout << "\nPerforming " << gradientIt << " conjugate gradient iteraitons\n\n";
	
	// Gradient and conjugate direction vectors
	vector<double> gradVector(cons.numTotalParams);
	vector<double> gradVectorNew(cons.numTotalParams);
	vector<double> cgVector(cons.numTotalParams);
	vector<double> cgVectorNew(cons.numTotalParams);
	
	vector<double> tempParams(cons.numTotalParams);
	
	long double currentF;
	long double previousF = error_from_trial_point(cons, initialParams, currentParams, 0, vecs, 0);
	
	ofstream cgStream;
	cgStream.precision(12);
	cgStream.open("cg_energy.txt");
	
	// Do initial steepest descent step with full line search
	compute_gradient(cons, vecs, initialParams, currentParams, gradVector);	
	
	double alpha = 0.1, tau = 0.4; // Starting step size and reduction factor
	double currentBestAlpha = alpha;
		
	const int nMax = 45; // Maximum number of iterations in line search
	long double lineSearchF[nMax];
	long double currentBestF = previousF;
		
	// Perform line search
	for (int n=0; n<nMax; n++)
	{
		for (int col=0; col<(cons.numTotalParams); col++)
			tempParams[col] = currentParams[col] - alpha*gradVector[col];
		
		// Write corresponding value of F
		lineSearchF[n] = error_from_trial_point(cons, initialParams, tempParams, 0, vecs, 0);
		
		if (lineSearchF[n] < currentBestF)
		{
			currentBestF = lineSearchF[n];
			currentBestAlpha = alpha;
		}
		
		alpha *= tau;
	}
	
	for (int col=0; col<cons.numTotalParams; col++)
	{
		cgVector[col] = gradVector[col];
		currentParams[col] -= currentBestAlpha*gradVector[col];
	}		
		
	previousF = error_from_trial_point(cons, initialParams, currentParams, 0, vecs, 0);
	
	// Do modified conjugate gradient steps -----------------------------------------------------------------------------
	for (int iG=0; iG < gradientIt; iG++)
	{
		
		// Get gradient corresponding to x_n, write to gradVectorNew
		compute_gradient(cons, vecs, initialParams, currentParams, gradVectorNew);	
		
		cout << "\nNew grad: ";
		for (int row=0; row<(cons.numTotalParams); row++)
			cout << std::setprecision(8) << std::setw(11) << gradVectorNew[row] << " ";
		cout << endl << endl;

		
		// Calculate dot-products needed for beta and theta coefficients
		double dotProd1 = 0.0, dotProd2 = 0.0, dotProd3 = 0.0;
		
		for (int row=0; row<(cons.numTotalParams); row++)
		{
			dotProd1 += gradVectorNew[row]*(gradVectorNew[row]-gradVector[row]);
			dotProd2 += gradVector[row]*gradVector[row];
			dotProd3 += gradVectorNew[row]*cgVector[row];
		}
		
		
		// Calculate coefficient beta (Polak-Ribiere)
		double betaPR = std::max(dotProd1/dotProd2, 0.0);
		//cout << "betaPR = " << betaPR << endl << endl;
		
		// Calculate coefficient theta (Zhang)
		double thetaZ = dotProd3/dotProd2;	
		
		// Update conjugate direction, advance vectors n -> n+1
		double magCgVec = 0.0;
		for (int row=0; row<(cons.numTotalParams); row++)
		{
			cgVectorNew[row] = gradVectorNew[row] + betaPR*cgVector[row] - thetaZ*(gradVectorNew[row]-gradVector[row]);
			cgVector[row] = cgVectorNew[row];
			gradVector[row] = gradVectorNew[row];
			magCgVec += cgVector[row]*cgVector[row];
		}
		
		// Perform line search in conjugate direction
		int n=0;
		alpha = 0.125;
		double deltaZ = 0.0001;
		double deltaCG = -1.0;
		while ( (n<=nMax) && (deltaCG <= 0.0) )
		{
			// Get updated params from x_i+1 = x_i - alpha*d_i+1
			for (int col=0; col<(cons.numTotalParams); col++)
				tempParams[col] = currentParams[col] - alpha*cgVectorNew[col];

			currentF = error_from_trial_point(cons, initialParams, tempParams, 0, vecs, 0);
			
			// Test Armijo-type condition 
			
			deltaCG = previousF - currentF - deltaZ*alpha*alpha*magCgVec; // should be >= 0 to continue
			
			//cout << alpha << ": previousF " << std::setw(10) << previousF << "  currentF " << std::setw(10) << currentF;
			//cout << "  delta*alpha^2*d^2" << std::setw(10) << deltaZ*alpha*alpha*magCgVec << "  CGdelta " << std::setw(10) << deltaCG << endl;
			
			alpha *= tau;
			n++;
		}
		
		// Overwrite curent parameters
		previousF = currentF;
		cgStream << currentF << endl;
		cout << "Updating parameters, F = " << currentF << endl;
		
		for (int col=0; col<(cons.numTotalParams); col++)
		{
			currentParams[col] = tempParams[col];
			cout << std::setprecision(8) << std::setw(11)  << currentParams[col] << " ";
		}	
		cout << endl;
	}		
	
	return 0;
} 

// Steepest descent function is not maintained
int steepest_descent(constant_struct cons, vector_struct vecs, vector<double> initialParams, vector<double> &currentParams, int gradientIt)
{
	vector<double> gradVector(cons.numTotalParams);
	vector<double> tempParams(cons.numTotalParams);
	long double currentF;
	long double previousF = error_from_trial_point(cons, initialParams, currentParams, 0, vecs, 0);
	
	for (int iG=0; iG < gradientIt; iG++)
	{		
		// Get gradient corresponding to current parameters
		compute_gradient(cons, vecs, initialParams, currentParams, gradVector);	
		
		cout << "\nGrad: ";
		for (int row=0; row<cons.numTotalParams; row++)
			cout << std::setw(8) << gradVector[row] << " ";
		cout << endl << endl;
		
		// Get ||g||^2
		double gradMag = 0.0;
		for (int row=0; row<cons.numTotalParams; row++)
			gradMag += gradVector[row]*gradVector[row];
		
		// Update parameters using a backtracking line search
		double alpha = 0.1, tau = 0.5, delta = -1.0;
		int n = 0, maxN = 20;
		
		while ( (n<=maxN) && (delta < 0.0) )
		{
			for (int col=0; col<cons.numTotalParams; col++)
			{
				tempParams[col] = currentParams[col] - alpha*gradVector[col];
				//cout << tempParams[col] << " ";
			}
			//cout << endl;
			
			currentF = error_from_trial_point(cons, initialParams, tempParams, 0, vecs, 0);
			
			delta = previousF - currentF; // should be >= 0 to continue	
			//cout << n << ": previousF " << std::setw(12) <<  previousF << "  currentF " << std::setw(12) << currentF;
			//cout << "  alpha " << alpha << endl;
			
			alpha *= tau;
			n++;
		}
		
		// Overwrite curent parameters
		previousF = currentF;
		cout << "Updating parameters, F = " << currentF << endl;
		
		for (int col=0; col<cons.numTotalParams; col++)
		{
			currentParams[col] = tempParams[col];
			cout << std::setw(10) << currentParams[col] << " ";
		}	
		cout << endl;
	}
	
	return 0;	
}

int xyz_files_read(constant_struct cons, vector_struct &vecs)
{
	ifstream xyzStream;
	
	// Loop over full (usually 360 degree) range at specificed step size
	if (cons.genXyzFileNames)
	{
		int fileRow=0, f=0;
		
		if (cons.phiDim == 2)
		{
			for (int d1 = cons.phi1min; d1 <= cons.phi1max; d1 = d1 + cons.phi1step)
			{
				for (int d2 = cons.phi2min; d2 <= cons.phi2max; d2 = d2 + cons.phi2step)
				{
					// Determine weighting to account for double counting of edge points
					// This assumes phi2 goes from 0->360, so the 0/360 is double counted
					if ((int(cons.phi2min) == 0) && (int(cons.phi2max) == 360))
					{
						// 0.25 for corners
						if ( ((d1 == 0) || (d1 == cons.phi1max)) && ((d2 == 0) || (d2 == cons.phi2max)) )
						{
							vecs.energyWeighting[f] = 0.25;
						}
						// 0.5 for edges
						else if ( (d1 == 0) || (d1 == cons.phi1max) || (d2 == 0) || (d2 == cons.phi2max) )
						{
							vecs.energyWeighting[f] = 0.5;
						}
						// 1.0 for central points
						else
						{
							vecs.energyWeighting[f] = 1.0;	
						}
					}
					else
					{
						vecs.energyWeighting[f] = 1.0;	
					}
					//cout << "d1, d2, weighting: " << d1 << ", " << d2 << ", " << vecs.energyWeighting[f] << endl;
					f++;
				
					bool file_open = 1;
					int ts = 0;
			
					string tsStrPrev;
			
					// This method iterates over the -###.xyz suffix to find the last step in a geometry optimisation
					// Only needed if geometries are taken from NWChem optimisation, which they probably shouldn't
					while (file_open == 1)
					{
						// Generate file name
						string d1Str = std::to_string(d1);
						string d2Str = std::to_string(d2);
						string tsStr = std::to_string(ts);
				
						// time-step code is 3-digits
						if (ts < 10)
							tsStr = "00" + tsStr;
						else if (ts < 100)
							tsStr = "0" + tsStr;
						else if (ts > 999)
							cout << "Warning: timestep has 4 digits (undefined behaviour)\n";
						
						// Generate file name		
						string fullFileName = cons.xyzFile + d1Str + "_" + d2Str + "-" + tsStr + ".xyz";

						//cout << "Opening xyz file " << fullFileName << endl;
				
						xyzStream.open(fullFileName.c_str());
				
						// Once we cannot open the file, we read data from the previous file
						if (!xyzStream.is_open())
						{
							file_open = 0;
							xyzStream.close();
							// Get previous file name
							fullFileName = cons.xyzFile + d1Str + "_" + d2Str + "-" + tsStrPrev + ".xyz";
							xyzStream.open(fullFileName.c_str());

							// Ignore first 2 lines and atom name (as per XYZ format)
							int line1;
							string line2, atom;
							xyzStream >> line1 >> line2;
					
							for (int a = 0; a < cons.size[0]; a++)
							{
								xyzStream >> atom >> vecs.xyzData[fileRow*3] >> vecs.xyzData[fileRow*3 + 1] >> vecs.xyzData[fileRow*3 + 2];
								//cout << vecs.xyzData[fileRow*3] << " " << vecs.xyzData[fileRow*3+1];
								//cout << " " << vecs.xyzData[fileRow*3+2] << endl;
								fileRow++; // Move onto next atom
						
							}
						}
						tsStrPrev = tsStr;
						xyzStream.close();
						ts++;	
					}	
				}
			}
		}
		else if (cons.phiDim == 3)
		{
			for (int d1 = cons.phi1min; d1 <= cons.phi1max; d1 = d1 + cons.phi1step)
			{
				for (int d2 = cons.phi2min; d2 <= cons.phi2max; d2 = d2 + cons.phi2step)
				{
					for (int d3 = cons.phi3min; d3 <= cons.phi3max; d3 = d3 + cons.phi3step)
					{	
						bool file_open = 1;
						int ts = 0;
						
						// Here weightings are always set to 1.0 - assumes there are no duplicate geometries
						vecs.energyWeighting[f] = 1.0;
						f++;
			
						string tsStrPrev;
			
						// This method iterates over the -###.xyz suffix to find the last step in a geometry optimisation. Only needed if geometries are taken from NWChem optimisation
						while (file_open == 1)
						{
							// Generate file name
							string d1Str = std::to_string(d1);
							string d2Str = std::to_string(d2);
							string d3Str = std::to_string(d3);
							string tsStr = std::to_string(ts);
				
							// time-step code is 3-digits
							if (ts < 10)
								tsStr = "00" + tsStr;
							else if (ts < 100)
								tsStr = "0" + tsStr;
							else if (ts > 999)
								cout << "Warning: timestep has 4 digits (undefined behaviour)\n";
						
							// Generate file name		
							string fullFileName = cons.xyzFile + d1Str + "_" + d2Str + "_" + d3Str + "-" + tsStr + ".xyz";
				
							// Open file
					
							//cout << "Opening " << fileName << endl;
				
							xyzStream.open(fullFileName);
				
							// Once we cannot open the file, we read data from the previous file
							if (!xyzStream.is_open())
							{
								file_open = 0;
								xyzStream.close();
								// Get previous file name
								fullFileName = cons.xyzFile + d1Str + "_" + d2Str + "_" + d3Str + "-" + tsStrPrev + ".xyz";
								xyzStream.open(fullFileName.c_str());
								// Read data
								cout << "Reading file " << f << ", " << fullFileName << endl;
								// Ignore first 2 lines and atom name (as per XYZ format)
								int line1;
								string line2, atom;
								xyzStream >> line1 >> line2;
					
								for (int a = 0; a < cons.size[0]; a++)
								{
									xyzStream >> atom >> vecs.xyzData[fileRow*3] >> vecs.xyzData[fileRow*3 + 1] >> vecs.xyzData[fileRow*3 + 2];
									//cout << vecs.xyzData[fileRow*3] << " " << vecs.xyzData[fileRow*3+1];
									//cout << " " << vecs.xyzData[fileRow*3+2] << endl;
									fileRow++; // Move onto next atom
						
								}
							}
							tsStrPrev = tsStr;
							xyzStream.close();
							ts++;
						
						}		
					}	
				}
			}
		}
	}
	else // Read phi coord pairs from file
	{
		assert(cons.phiDim == 2);
		
		ifstream phiPairStream, xyzStream;
		cout << "Reading phi pairs from " << cons.phiCoordFile << endl;
			
		phiPairStream.open(cons.phiCoordFile.c_str());
		if(!phiPairStream.is_open())
		{
			cout << "Error, phi pair file cannot be opened!" << endl;
			return 1;
		}
		
		int fileRow = 0;
		for (int f=0; f<cons.numConfigs; f++)
		{	
			// Set weighting
			vecs.energyWeighting[f] = 1.0;
			
			// Read phi1 and phi2 strings
			char phi1Str[16], phi2Str[16];
	
			// Read phi pair strings from phiCoordFile
			// Only supports 2D for now
			phiPairStream >> vecs.phiPairSpec[f*2] >> vecs.phiPairSpec[f*2 +1];
			
			snprintf(phi1Str, 16, "%.1f", vecs.phiPairSpec[f*2]);
			snprintf(phi2Str, 16, "%.1f", vecs.phiPairSpec[f*2 +1]);
			
			// Generate xyz filename
			// This method assumes the files end with "-000.xyz"
			string fileName = cons.xyzFile + phi1Str + "_" + phi2Str + "-000" + ".xyz";
			
			//cout << "Opening xyz file: " << fileName << endl;
			
			xyzStream.open(fileName.c_str());
			
			// Ignore first 2 lines and atom name (as per XYZ format)
			int line1;
			string line2, atom;
			xyzStream >> line1 >> line2;
		
			for (int a = 0; a < cons.size[0]; a++) 
			{
				xyzStream >> atom >> vecs.xyzData[fileRow*3] >> vecs.xyzData[fileRow*3 + 1] >> vecs.xyzData[fileRow*3 + 2];
				//cout << vecs.xyzData[fileRow*3] << " " << vecs.xyzData[fileRow*3+1];
				//cout << " " << vecs.xyzData[fileRow*3+2] << endl;
				fileRow++; // Move onto next atom			
			}
			xyzStream.close();
		}
	}
	return 0;
}

int connectivity_read(constant_struct cons, vector_struct &vecs)
{
	ifstream connectStream;	
	cout << "Opening " << cons.connectFile << endl;
	connectStream.open(cons.connectFile.c_str());
	// Skip first two lines
	char s[128];
	connectStream.getline(s,128);
	connectStream.getline(s,128);
	
	string fileSection;
	// Read atoms (with charges and van der Waals parameters)
	connectStream >> fileSection;
	if (fileSection == "atoms")
	{
		for (int a=0; a<cons.size[0]; a++)
		{
			for (int i=0; i<cons.atomDataSize; i++)
				connectStream >> vecs.atomData[a*cons.atomDataSize + i];
		}
	}
	else
	{
		cout << "Error in input file format (no atoms section)\n\n";
		return 1;
	}
	// Read bonds
	connectStream >> fileSection;
	if (fileSection == "bonds")
	{
		for (int b=0; b<cons.size[1]; b++)
		{
			for (int i=0; i<cons.bondDataSize; i++)
				connectStream >> vecs.bondData[b*cons.bondDataSize + i];
		}
	}
	else
	{
		cout << "Error in input file format (no bonds section)\n\n";
		return 1;
	}
	
	// Read angles
	connectStream >> fileSection;
	if (fileSection == "angles")
	{
		for (int c=0; c<cons.size[2]; c++)
		{
			for (int i=0; i<cons.angleDataSize; i++)
				connectStream >> vecs.angleData[c*cons.angleDataSize + i];
		}
	}
	else
	{
		cout << "Error in input file format (no angles section)\n\n";
		return 1;
	}
	
	// Read dihedrals
	connectStream >> fileSection;
	if (fileSection == "dihedrals")
	{
		for (int d=0; d<cons.size[3]; d++)
		{
			for (int i=0; i<cons.dihedralDataSize; i++)
				connectStream >> vecs.dihedralData[d*cons.dihedralDataSize + i];		
		}
	}
	else
	{
		cout << "Error in input file format (no dihedrals section)\n\n";
		return 1;
	}
	
	// Read harmonic improper dihedrals
	connectStream >> fileSection;
	if (fileSection == "impropers")
	{
		for (int id=0; id<cons.size[4]; id++)
		{
			for (int i=0; i < cons.improperDataSize; i++)
				connectStream >> vecs.improperData[id*cons.improperDataSize + i];
		}
	}
	else
	{
		cout << "Error in input file format (no impropers section)\n\n";
		return 1;
	}
	
	// Read lj pair to fit
	connectStream >> fileSection;
	if (fileSection == "lj_pair_fit")
	{
		for (int e = 0; e < cons.size[5]; e++)
		{
			for (int i = 0; i < cons.pairDataSize; i++)
			{
				connectStream >> vecs.ljPairFitData[e*cons.pairDataSize + i];
			}
		}
	}
	else
	{
		cout << "Error in input file format (no lj_pair_fit section)\n\n";
		return 1;
	}
	
	return 0;
}

int connectivity_process(constant_struct cons, vector_struct &vecs)
{
	// Read atoms separated by one bond
	for (int i=0; i < cons.size[0]; i++)
	{
		for (int j=i; j < cons.size[0]; j++)
		{
			int m1 = i*cons.size[0] + j; // Convert to 1D array coordinate
			int m2 = j*cons.size[0] + i; // Symmetric matrix
			
			int atomI = int(round(vecs.atomData[i*cons.atomDataSize])); // Get actual atom number from file
			int atomJ = int(round(vecs.atomData[j*cons.atomDataSize]));
			
			// Search bonds list
			bool bondFound = 0;
			for (int b=0; b<cons.size[1]; b++)
			{				
				if ( (vecs.bondData[b*cons.bondDataSize] == atomI) && (vecs.bondData[b*cons.bondDataSize + 1] == atomJ) )
				{
					vecs.bondSepMat[m1] = 1;
					vecs.bondSepMat[m2] = 1;
					bondFound = 1;
				}
			}
			if (bondFound == 0)
			{
				vecs.bondSepMat[m1] = 0;
				vecs.bondSepMat[m2] = 0;
			}			
		}
	}
	
	// Now loop over increasing bond separation number until array is populated
	int matUpdate = 1;
	int bondSearch = 2;
	while (matUpdate == 1)
	{
		matUpdate = 0; // This will be set to 1 if an update was made
		
		for (int i=0; i < cons.size[0]; i++)
		{
			for (int j=(i+1); j < cons.size[0]; j++) // j>i, since i=j are the same atom
			{
				int m1 = i*cons.size[0] + j; // Convert to 1D array coordinate
				int m2 = j*cons.size[0] + i; // Symmetric matrix
			
				// If number of bonds between pair i-j is not already determined,
				// look for atoms which are bondSearch-1 from atom i				
				if (vecs.bondSepMat[m1] == 0)
				{
					for (int k=0; k < cons.size[0]; k++)
					{
						int m3 = k*cons.size[0] + j;
						// If connection (k,j) is bondSearch-1, check if (k,i) is 1
						if (vecs.bondSepMat[m3] == (bondSearch - 1))
						{	
							int m4 = k*cons.size[0] + i;
							
							if(vecs.bondSepMat[m4] == 1)
							{
								vecs.bondSepMat[m1] = bondSearch;
								vecs.bondSepMat[m2] = bondSearch; 
								matUpdate = 1;
							}
						}
					}
				}
			}
		}
		bondSearch++;
	}
	
	cout << "Connectivity matrix complete:\n";
	for (int i=0; i < cons.size[0]; i++)
	{
		for (int j=0; j < cons.size[0]; j++)
		{
			cout << vecs.bondSepMat[i*cons.size[0]+j] << " ";
		}
		cout << endl;
	}
	
	// Work out sigma, epsilon and chargeQQ for each pair
	for (int i=0; i < cons.size[0]; i++)
	{
		for (int j=i; j < cons.size[0]; j++)
		{
			int m1 = i*cons.size[0]+j;
			int a1 = int(round(vecs.atomData[i*cons.atomDataSize]));
			int a2 = int(round(vecs.atomData[j*cons.atomDataSize]));
			//int m2 = j*cons.size[0]+i;
			bool toFit = 0;
			
			// Loop over pairs to work out if this is the LJ pair to be fitted
			for (int p=0; p < cons.size[5]; p++)
			{
				if ( (vecs.ljPairFitData[p*cons.pairDataSize] == a1) && (vecs.ljPairFitData[p*cons.pairDataSize + 1] == a2) )
					toFit = 1;		
			}

			assert(cons.gromacsCombRule == 2);
			double sigma = 0.5*(vecs.atomData[i*cons.atomDataSize+2]+vecs.atomData[j*cons.atomDataSize+2]);
			double epsilon = sqrt(vecs.atomData[i*cons.atomDataSize+3]*vecs.atomData[j*cons.atomDataSize+3]);
			double QQ = vecs.atomData[i*cons.atomDataSize+1]*vecs.atomData[j*cons.atomDataSize+1];
			
			// Excluded interactions
			if (vecs.bondSepMat[m1] <= cons.nrexcl)
			{
				vecs.sigmaMatrix[m1] = 0.0;
				vecs.epsilonMatrix[m1] = 0.0;
				vecs.qqMatrix[m1] = 0.0;
			}
			// 1-4 Interactions (can be removed by setting scale14LJ & scale14QQ to zero)
			if (vecs.bondSepMat[m1] == 3)
			{
				if (sigma > cons.sigma14factorCutoff) {
					vecs.sigmaMatrix[m1] = sigma*cons.sigma14factor;
				}
				else {
					vecs.sigmaMatrix[m1] = sigma;
				}
				vecs.epsilonMatrix[m1] = cons.scale14LJ*epsilon;
				vecs.qqMatrix[m1] = cons.scale14QQ*QQ;
			}
			// Full interactions
			if (vecs.bondSepMat[m1] > cons.nrexcl)
			{
				vecs.sigmaMatrix[m1] = sigma;
				vecs.epsilonMatrix[m1] = epsilon;
				vecs.qqMatrix[m1] = QQ;
			}
			if (toFit == 1)
			{
				//Remove LJ interaction for pairs being fit (usually 1-5 pairs)
				vecs.epsilonMatrix[m1] = 0.0;
			}	
		}
	}
	
	cout << "sigmaLJ matrix complete:\n";
	for (int i = 0; i < cons.size[0]; i++)
	{
		for (int j = 0; j < cons.size[0]; j++)
		{
			cout << std::setw(6) << vecs.sigmaMatrix[i*cons.size[0]+j] << " ";
		}
		cout << endl;
	}
	cout << "epsilonLJ matrix complete:\n";
	for (int i = 0; i < cons.size[0]; i++)
	{
		for (int j = 0; j < cons.size[0]; j++)
		{
			cout << std::setw(6) << vecs.epsilonMatrix[i*cons.size[0]+j] << " ";
		}
		cout << endl;
	}
	cout << "chargeQQ matrix complete:\n";
	for (int i = 0; i < cons.size[0]; i++)
	{
		for (int j = 0; j < cons.size[0]; j++)
		{
			cout << std::setw(6) << vecs.qqMatrix[i*cons.size[0]+j] << " ";
		}
		cout << endl;
	}
	return 0;
}

int constant_energy_process(constant_struct cons, vector_struct &vecs)
{
	int print_rij = 0;
	vector<double> coords(3*cons.size[0]);
	double aVec1[3], aVec2[3];
	double dVec1[3], dVec2[3], dVec3[3];
	double crossVec1[3], crossVec2[3], crossVec3[3];
	ofstream constEnStream;
	constEnStream.open("energyConst.txt");
	
	cout << "Processing constant energy terms in " << cons.numConfigs << " geometries.\n\n";
	
	// Loop over configurations
	for (int f=0; f<cons.numConfigs; f++)
	{
		double energyAngle = 0.0, energyDihedral = 0.0, energyImproper = 0.0;
		double energyLJ = 0.0, energyQQ = 0.0;
		double energy14 = 0.0, energy15 = 0.0, energyNB = 0.0;
		vecs.constantEnergy[f] = 0.0;
		
		// Read coordinates
		for (int k=0; k<cons.size[0]; k++)
		{
			coords[k*3    ] = vecs.xyzData[f*cons.size[0]*3 + k*3    ];
			coords[k*3 + 1] = vecs.xyzData[f*cons.size[0]*3 + k*3 + 1];
			coords[k*3 + 2] = vecs.xyzData[f*cons.size[0]*3 + k*3 + 2];
		}
		// Calculate bond energy (currently ignored because TraPPE uses fixed bonds)
		//
		// Calculate angle (bending) energy
		for (int c=0; c < cons.size[2]; c++)
		{
			int i = int(round(vecs.angleData[c*cons.angleDataSize    ] - 1));
			int j = int(round(vecs.angleData[c*cons.angleDataSize + 1] - 1));
			int k = int(round(vecs.angleData[c*cons.angleDataSize + 2] - 1));
			// -1 Because atoms are numbered starting from one in the input file
			
			double theta0 = vecs.angleData[c*cons.angleDataSize + 3]*cons.pi/180.0;
			double kTheta = vecs.angleData[c*cons.angleDataSize + 4];
	
			// vec j->i
			for (int m=0; m < 3; m++)
				aVec1[m] = coords[i*3 + m] - coords[j*3 + m];
			
			double modVec1 = sqrt(aVec1[0]*aVec1[0] 
				+ aVec1[1]*aVec1[1] + aVec1[2]*aVec1[2]);
			
			// vec j->k
			for (int m=0; m < 3; m++)
				aVec2[m] = coords[k*3 + m] - coords[j*3 + m];
			
			double modVec2 = sqrt(aVec2[0]*aVec2[0] + aVec2[1]*aVec2[1] + aVec2[2]*aVec2[2]);

			double dot12 = aVec1[0]*aVec2[0] + aVec1[1]*aVec2[1] + aVec1[2]*aVec2[2];

			double theta = acos(dot12/(modVec1*modVec2));
			energyAngle += 0.5*kTheta*(theta-theta0)*(theta-theta0);
			//cout << "theta, theta0: " << theta*180.0/cons.pi << ", " << theta0*180.0/cons.pi << endl;
		}
		
		// Calculate dihedral angles and constant energy contribution
		int nDfit = 0;
		for (int d=0; d < cons.size[3]; d++)
		{
			int i = int(round(vecs.dihedralData[d*cons.dihedralDataSize    ] - 1));
			int j = int(round(vecs.dihedralData[d*cons.dihedralDataSize + 1] - 1));
			int k = int(round(vecs.dihedralData[d*cons.dihedralDataSize + 2] - 1));
			int l = int(round(vecs.dihedralData[d*cons.dihedralDataSize + 3] - 1));
			//cout << "ijkl: " << i << " " << j << " " << k << " " << l << endl;

			// Calculate vectors b1 = rj-ri, b2 = rk-rj, b3 = rl-rk
			for (int m=0; m < 3; m++)
			{
				dVec1[m] = coords[j*3 + m] - coords[i*3 + m];
				dVec2[m] = coords[k*3 + m] - coords[j*3 + m];
				dVec3[m] = coords[l*3 + m] - coords[k*3 + m];
			}
				
			// Calculate normalised cross products c1 = <b1xb2>, c2 = <b2xb3>
			crossVec1[0] = dVec1[1]*dVec2[2] - dVec1[2]*dVec2[1];
			crossVec1[1] = dVec1[2]*dVec2[0] - dVec1[0]*dVec2[2];
			crossVec1[2] = dVec1[0]*dVec2[1] - dVec1[1]*dVec2[0];
			double modCrossVec1 = sqrt(crossVec1[0]*crossVec1[0] 
				+ crossVec1[1]*crossVec1[1] + crossVec1[2]*crossVec1[2]);
					
			crossVec2[0] = dVec2[1]*dVec3[2] - dVec2[2]*dVec3[1];
			crossVec2[1] = dVec2[2]*dVec3[0] - dVec2[0]*dVec3[2];
			crossVec2[2] = dVec2[0]*dVec3[1] - dVec2[1]*dVec3[0];
			double modCrossVec2 = sqrt(crossVec2[0]*crossVec2[0] 
				+ crossVec2[1]*crossVec2[1] + crossVec2[2]*crossVec2[2]);
				
			for (int m=0; m < 3; m++)
			{
				crossVec1[m] /= modCrossVec1;
				crossVec2[m] /= modCrossVec2;
			}
				
			// Basis vector: d1 = cross(c1,b2)
			crossVec3[0] = crossVec1[1]*dVec2[2] - crossVec1[2]*dVec2[1];
			crossVec3[1] = crossVec1[2]*dVec2[0] - crossVec1[0]*dVec2[2];
			crossVec3[2] = crossVec1[0]*dVec2[1] - crossVec1[1]*dVec2[0];
			double modCrossVec3 = sqrt(crossVec3[0]*crossVec3[0] 
				+ crossVec3[1]*crossVec3[1] + crossVec3[2]*crossVec3[2]);
				
			for (int m=0; m < 3; m++)
				crossVec3[m] /= modCrossVec3;
			   
			double x = crossVec1[0]*crossVec2[0] + crossVec1[1]*crossVec2[1] +
				crossVec1[2]*crossVec2[2]; //dot(c1,c2);
			double y = -(crossVec3[0]*crossVec2[0] + crossVec3[1]*crossVec2[1] +
				crossVec3[2]*crossVec2[2]); //dot(d1,c2);
        
			double phi = atan2(y,x);
			double psi = phi - cons.pi;
			
			//cout << "d, phi: " << d << ", " << phi*180.0/cons.pi << endl;
			
			// Calculate cosine psi
			double cosPsi = cos(psi);
				
			// Add cosines for dihedrals to be used in fitting procedure
			int dType = int(round(vecs.dihedralData[d*cons.dihedralDataSize + 4]));
			if (dType != 0)
			{				
				vecs.cosPsiData[f*cons.numDihedralFit + nDfit] = cosPsi;
				nDfit++;
			}
			else // Else add to constant energy
			{
				for (int rb = 0; rb < cons.nRBfit; rb++) {
					energyDihedral += vecs.dihedralData[d*cons.dihedralDataSize + 5 + rb]*pow(cosPsi,rb);
				}
			}				
		}
		
		// Calculate harmonic improper dihedral energy
		for (int id=0; id<cons.size[4]; id++)
		{
			// Calculate angle between planes ijk and jkl
			int i = int(round(vecs.improperData[id*cons.improperDataSize    ] - 1)); // -1 for 0-based arrays
			int j = int(round(vecs.improperData[id*cons.improperDataSize + 1] - 1));
			int k = int(round(vecs.improperData[id*cons.improperDataSize + 2] - 1));
			int l = int(round(vecs.improperData[id*cons.improperDataSize + 3] - 1));
			//cout << "ijkl: " << i << " " << j << " " << k << " " << l << endl;

			// Calculate vectors b1 = rj-ri, b2 = rk-rj, b3 = rl-rk
			for (int m=0; m < 3; m++)
			{
				dVec1[m] = coords[j*3 + m] - coords[i*3 + m];
				dVec2[m] = coords[k*3 + m] - coords[j*3 + m];
				dVec3[m] = coords[l*3 + m] - coords[k*3 + m];
			}
				
			// Calculate normalised cross products c1 = <b1xb2>, c2 = <b2xb3>
			crossVec1[0] = dVec1[1]*dVec2[2] - dVec1[2]*dVec2[1];
			crossVec1[1] = dVec1[2]*dVec2[0] - dVec1[0]*dVec2[2];
			crossVec1[2] = dVec1[0]*dVec2[1] - dVec1[1]*dVec2[0];
			double modCrossVec1 = sqrt(crossVec1[0]*crossVec1[0] 
				+ crossVec1[1]*crossVec1[1] + crossVec1[2]*crossVec1[2]);
					
			crossVec2[0] = dVec2[1]*dVec3[2] - dVec2[2]*dVec3[1];
			crossVec2[1] = dVec2[2]*dVec3[0] - dVec2[0]*dVec3[2];
			crossVec2[2] = dVec2[0]*dVec3[1] - dVec2[1]*dVec3[0];
			double modCrossVec2 = sqrt(crossVec2[0]*crossVec2[0] 
				+ crossVec2[1]*crossVec2[1] + crossVec2[2]*crossVec2[2]);
				
			for (int m=0; m < 3; m++)
			{
				crossVec1[m] /= modCrossVec1;
				crossVec2[m] /= modCrossVec2;
			}
				
			// Basis vector: d1 = cross(c1,b2)
			crossVec3[0] = crossVec1[1]*dVec2[2] - crossVec1[2]*dVec2[1];
			crossVec3[1] = crossVec1[2]*dVec2[0] - crossVec1[0]*dVec2[2];
			crossVec3[2] = crossVec1[0]*dVec2[1] - crossVec1[1]*dVec2[0];
			double modCrossVec3 = sqrt(crossVec3[0]*crossVec3[0] 
				+ crossVec3[1]*crossVec3[1] + crossVec3[2]*crossVec3[2]);
				
			for (int m=0; m < 3; m++)
				crossVec3[m] /= modCrossVec3;
			   
			double x = crossVec1[0]*crossVec2[0] + crossVec1[1]*crossVec2[1] 
				+ crossVec1[2]*crossVec2[2];
			double y = -(crossVec3[0]*crossVec2[0] + crossVec3[1]*crossVec2[1] 
				+ crossVec3[2]*crossVec2[2]);
        
			double phi = atan2(y,x);
			
			// Input params are in degrees and (kJ/mol)/rad^2
			double phi0 = (cons.pi/180.0)*vecs.improperData[id*cons.improperDataSize + 4];
			double kHarm = vecs.improperData[id*cons.improperDataSize + 5];
			
			// Shift phi into range phi0-180 <= phi <= phi0+180
			if (phi <= phi0-cons.pi)
				phi += 2*cons.pi;
			else if (phi >= phi0+cons.pi)
				phi -= 2*cons.pi;
			
			// Sum energy
			energyImproper += kHarm*(phi-phi0)*(phi-phi0);
			//cout << "phi, phi0, k: " << (180.0/cons.pi)*phi << " " << (180.0/cons.pi)*phi0;
			//cout << " " << kHarm << endl;
		}
		
		
		// Calculate LJ and QQ non-bonded energy
		double rIJ[3];
		
		for (int i=0; i < cons.size[0]; i++)
		{
			for (int j=(i+1); j < cons.size[0]; j++) // j>i
			{
				int m1 = i*cons.size[0]+j;
				int atomI = int(round(vecs.atomData[i*cons.atomDataSize])); // Get actual atom number from file data
				int atomJ = int(round(vecs.atomData[j*cons.atomDataSize])); 
				int atomIC = atomI - 1; // Subtract 1 for C arrays
				int atomJC = atomJ - 1;
				
				// Vector pointing from i to j
				for (int n=0; n<3; n++)
					rIJ[n] = coords[atomJC*3 + n] - coords[atomIC*3 + n];
				
				// Convert A to nm
				if (cons.xyzAngstroms)
				{
					for (int n=0; n<3; n++)
						rIJ[n] *= 0.1;
				}
				
				// Pair separation
				double r = sqrt(rIJ[0]*rIJ[0] + rIJ[1]*rIJ[1] + rIJ[2]*rIJ[2]);
				
				vecs.rijMatrix[m1] = r;
				
				// Add this pair to vecs.pairSepData if it's on the list of pairs to fit
				for (int pair=0; pair < cons.size[5]; pair++)
				{
					if ( (vecs.ljPairFitData[pair*cons.pairDataSize] == atomI) 
						&& (vecs.ljPairFitData[pair*cons.pairDataSize+1] == atomJ) )
							vecs.pairSepData[f*cons.size[5] + pair] = r;
				}
				
				double LJ6 = pow(vecs.sigmaMatrix[m1]/r,6);
				double energyPairLJ = 4.0*vecs.epsilonMatrix[m1]*(LJ6*LJ6-LJ6);
				
				// Coulomb
				// 138.9355 converts coulomb energy |q||q|/r[nm] to kJ/mol
				double energyPairQQ = 138.9355*vecs.qqMatrix[m1]/r;
				
				// Split into 1-4, 1-5 and 1-6+ parts for plotting
				if (vecs.bondSepMat[m1] == 3)
					energy14 += (energyPairQQ + energyPairLJ);
				else if (vecs.bondSepMat[m1] == 4)
					energy15 += (energyPairQQ + energyPairLJ);
				else
					energyNB += (energyPairQQ + energyPairLJ);
				
				// Add to total as well
				energyLJ += energyPairLJ;
				energyQQ += energyPairQQ;
			
			}
		}
		
		if (f == print_rij)
		{
			cout << "r_ij matrix complete for config " << print_rij << ":\n";
			for (int i = 0; i < cons.size[0]; i++)
			{
				for (int j = 0; j < cons.size[0]; j++)
				{
					cout << std::setw(10) << vecs.rijMatrix[i*cons.size[0]+j] << " ";
				}
				cout << endl;
			}
		}
		
		vecs.constantEnergy[f] = energyAngle + energyDihedral + energyImproper + energyLJ + energyQQ;
		
		//cout << "f = " << f << endl << "Angle:    " << energyAngle << endl;
		//cout << "Dihedral: " << energyDihedral << endl << "Improper: " << energyImproper << endl;
		//cout << "LJ:       " << energyLJ << endl << "QQ:       " << energyQQ << endl << endl;
		
		constEnStream << energyAngle << ", " << energyDihedral << ", " << energyImproper << ", ";
		constEnStream << energy14 << ", " << energy15 << ", " << energyNB << endl;
	}
	constEnStream.close();
	return 0;
}

int energy_read(constant_struct cons, vector_struct &vecs)
{
	ifstream energyStream;
	energyStream.open(cons.energyFile.c_str());
	
	if(!energyStream.is_open())
	{
		cout << "ERROR: Energy file was not opened!\n";
		return 1;
	}
	
	for (int f=0; f<cons.numConfigs; f++)
	{
		energyStream >> vecs.energyData[f];
		vecs.energyWeighting[f] *= exp(-vecs.energyData[f]/cons.KT); // Apply Boltzmann weighting 
		//cout << vecs.energyWeighting[f] << endl;
	}
	return 0;
}

int define_initial_simplex(constant_struct cons, vector<double> initialParams, vector<double> &simplex)
{
	int ssM1 = cons.numTotalParams; // Width of simplex, i.e. total number of params
	const double factorLJ = -0.001;
	const double dihedralStep = 0.5;
	
	// Fill in simplex, add perturbation to diagonal elements
	for (int n=0; n<cons.simplexSize; n++)
	{
		// Dihedral parameters
		for (int iD=0; iD<cons.numDihedralParams; iD++)
		{
			simplex[n*ssM1 + iD] = initialParams[iD];
			if (iD%ssM1 == n)
				simplex[n*ssM1 + iD] += dihedralStep;		
		}
		// LJ pair parameters
		if (cons.size[5] != 0)
		{
			int iD = cons.numDihedralParams;
			for (int j=0; j<2; j++)
			{
				simplex[n*ssM1 + iD] = initialParams[iD];
				if (iD%ssM1 == n)
					simplex[n*ssM1 + iD] *= (1.0 + factorLJ);
				iD++;
			}		
		}
	}	
	
	// Output initial simplex
	cout << "Starting simplex: " << endl;
	for (int row=0; row<cons.simplexSize; row++)
	{
		for (int col=0; col<cons.numTotalParams; col++)
		{
			cout << std::setw(10) << simplex[row*(cons.numTotalParams) + col] << " ";
		}
		cout << endl;
	}
	
	return 0;
}

long double error_from_trial_point(constant_struct cons, vector<double> initialParams, vector<double> &trialParams, int trialStart, vector_struct vecs, bool toWrite)
{
	ofstream energyMD;
	ofstream energyDist;
	vector<double> params(cons.numTotalParams);
	vector<double> energyTotal(cons.numConfigs);

	// Extract params from trialParams input vector
	// This could be a simplex matrix, in which case trialStart is the first element of the row of params to test
	for (int p=0; p<cons.numTotalParams; p++)
	{
		params[p] = trialParams[trialStart+p];
	}
	
	if (toWrite == 1)
	{
		energyMD.open("energyMD.txt");
		energyDist.open("energyDist.txt");
	}
	
	long double sumResid = 0.0;
	
	// Loop over configurations
	for (int f=0; f<cons.numConfigs; f++)
	{
		// Calculate dihedral energy for simplex 'd' and xyz config 'f'
		double energyDihedral = 0.0, energyPair = 0.0;
		int dd=0; 
		
		// Loop over all dihedrals
		for (int d=0; d<cons.size[3]; d++)
		{
			// Cosine(psi) for this dihedral
			double cosPsiPow = 1.0; // Start at cosPsi^0
			double coeff = 0.0;
			
			// Get diherdal type
			int dType = int(round(vecs.dihedralData[d*cons.dihedralDataSize + 4]));
			//cout << "dType = " << dType << endl;
			
			// Type 0 dihedrals are the ones not varying during the fit
			if (dType > 0)
			{	
				double cosPsi = vecs.cosPsiData[f*cons.numDihedralFit+dd];
				dd++;
				
				for (int rb=0; rb<cons.nRBfit; rb++)
				{
					// Get coefficient
					if ((dType > 1) && (rb == 0)) {
						coeff = 0.0;
					}
					else {
						coeff = params[(dType-1)*(cons.nRBfit-1) + rb];
					}				
				
					energyDihedral += cosPsiPow*coeff;
					
					// Raise cosPsiPow to the next power
					cosPsiPow *= cosPsi;
				}
			}
		}
		// Calculate improper dihedral energy
			
		// Calculate pair energy (LJ only)
		for (int pair=0; pair<cons.size[5]; pair++)
		{		
			// LJ parameters for this pair
			double sigmaPair = params[cons.numDihedralParams];
			double epsilonPair = params[cons.numDihedralParams + 1];
			double r = vecs.pairSepData[f*cons.size[5] + pair];
			double LJ6 = pow(sigmaPair/r, 6);
					
			energyPair += 4.0*epsilonPair*(LJ6*LJ6-LJ6);
		}
			
		// Total energy of configuration f for trial point
		energyTotal[f] = vecs.constantEnergy[f] + energyPair + energyDihedral;
			
		// Compare to input energy and sum weighted square of residual
		sumResid += (vecs.energyData[f]-energyTotal[f])*(vecs.energyData[f]-energyTotal[f])*vecs.energyWeighting[f];
		
		// Write to file
		if (toWrite == 1)
		{
			energyMD << energyTotal[f] << endl;
			energyDist << vecs.constantEnergy[f] << ", " << energyDihedral << ", " << energyPair << endl;	
		}
	} // End of config loop
	
	// Dihedral coefficient restraint terms
	for (int coeff=1; coeff<cons.numDihedralParams; coeff++) 
	{
		if (cons.resToZero == 1)
			sumResid += (cons.dihedralK*cons.numConfigs)*params[coeff]*params[coeff];
		else
			sumResid += (cons.dihedralK*cons.numConfigs)*(initialParams[coeff]-params[coeff])*(initialParams[coeff]-params[coeff]);
	}
	// Epsilon restraint term
	if (cons.size[5] != 0)
	{
		double epsDef = initialParams[cons.numTotalParams-1];
		sumResid += (cons.epsK*cons.numConfigs)*(epsDef - params[cons.numDihedralParams+1])*(epsDef - params[cons.numDihedralParams+1]);
	}
	
	// Boltzmann integral restraint terms
	if (cons.useBoltzIntRes == 1)
	{
		// Calculate Boltzmann integrals for each conformation
		
		
		
		
	}
	
	return sumResid;
}

int compute_gradient(constant_struct cons, vector_struct vecs, vector<double> initialParams, vector<double> currentParams, vector<double> &gradVector)
{
	// Reset grad vector
	for (int col=0; col<cons.numTotalParams; col++)
	{
		gradVector[col] = 0.0;
	}
	
	// Loop over configs
	for (int f=0; f<cons.numConfigs; f++)
	{
		// Calculate dihedral energy for simplex 'd' and xyz config 'f'
		double energyDihedral = 0.0, energyPair = 0.0, gradSigmaPair = 0.0, gradEpsPair = 0.0;
		int dd=0;
		
		// Calculate dihedral energy
		for (int d=0; d<cons.size[3]; d++)
		{
			// Cosine(psi) for this dihedral
			double cosPsiPow = 1.0; // Start at cosPsi^0
			double coeff = 0.0;
			
			// Get diherdal type
			int dType = int(round(vecs.dihedralData[d*cons.dihedralDataSize + 4]));
			
			// Type 0 dihedrals are the ones not varying during the fit
			if (dType > 0)
			{	
				double cosPsi = vecs.cosPsiData[f*cons.numDihedralFit+dd];
				dd++;
				
				for (int rb=0; rb<cons.nRBfit; rb++)
				{
					// Get coefficient
					if ((dType > 1) && (rb == 0))
						coeff = 0.0;
					else
					{
						coeff = currentParams[(dType-1)*(cons.nRBfit-1) + rb];
						//cout << "Coeff: " << f << ", " << ((dType-1)*(cons.nRBfit-1) + rb) << ", " << coeff << endl;
					}				
					energyDihedral += cosPsiPow*coeff;
					
					// Raise cosPsiPow to the next power
					cosPsiPow *= cosPsi;
				}
			}
		}
		
		// Calculate pair energy
		if (cons.size[5] != 0)
		{
			double sigmaPair = currentParams[cons.numDihedralParams];
			double epsilonPair = currentParams[cons.numDihedralParams+1];
			
			for (int pair=0; pair<cons.size[5]; pair++)
			{
				double r = vecs.pairSepData[f*cons.size[5] + pair];
				double LJ6 = pow(sigmaPair/r, 6);
			
				energyPair += 4.0*epsilonPair*(LJ6*LJ6-LJ6);
				gradSigmaPair += 4.0*epsilonPair*( LJ6*LJ6*12 - LJ6*6)/sigmaPair; 	// Derivative wrt sigma
				gradEpsPair += 4.0*(LJ6*LJ6-LJ6);									// Derivative wrt epsilon			
			}
		}
		
		// Compare total energy to DFT energy to get delta
		double energyDelta = vecs.energyData[f] - (vecs.constantEnergy[f] + energyPair + energyDihedral);
		//cout << "grad) f, energyPair = " << f << ", " << energyPair << endl;
	
		// Derivative is of form 2*(Delta)*Grad(Delta)*Weighting
		double gradWeighting = 2*vecs.energyWeighting[f];
	
		// Construct this configuration's contribution to the gradient vector
		gradVector[0] -= energyDelta*gradWeighting;
	
		// Loop over dihedrals again now that we have deltaU and weighting factors
		dd = 0;
		for (int d=0; d<cons.size[3]; d++)
		{
			int dType = int(round(vecs.dihedralData[d*cons.dihedralDataSize + 4]));			
			if (dType > 0)
			{	
				// Start at n=1 term in dihedral potential
				double cosPsi = vecs.cosPsiData[f*cons.numDihedralFit+dd];
				dd++;
				double cosPsiPow = cosPsi; // Start at cosPsi^1
			
				for (int rb=1; rb<cons.nRBfit; rb++)
				{
					gradVector[(dType-1)*(cons.nRBfit-1) + rb] -= cosPsiPow*energyDelta*gradWeighting;		
					//cout << "cosPsiPow " << cosPsiPow << endl;
					//cout << "energyDelta " << energyDelta << endl;
					//cout << "energyFactor " << energyWeighting << endl;	
									
					cosPsiPow *= cosPsi;
				}
			}
		}
		// Now gradients with respect to sigma and epsilon of pair being fit
		if (cons.size[5] != 0)
		{
			gradVector[cons.numDihedralParams    ] -= cons.sigmaGradFactor*gradSigmaPair*energyDelta*gradWeighting;
			gradVector[cons.numDihedralParams + 1] -= gradEpsPair*energyDelta*gradWeighting;	
		}
	}
	
	// Gradients due to constraint terms
	// Doesn't depend on config, so calculated once and multiplied by numConfigs
	for (int coeff=1; coeff<cons.numDihedralParams; coeff++) 
	{
		if (cons.resToZero == 1)
			gradVector[coeff] += 2*(cons.dihedralK*cons.numConfigs)*currentParams[coeff];
		else
			gradVector[coeff] += 2*(cons.dihedralK*cons.numConfigs)*(currentParams[coeff]-initialParams[coeff]);
	}	
	// Epsilon constraint term
	if (cons.size[5] != 0)
	{
		double epsDef = initialParams[cons.numTotalParams-1];
		gradVector[cons.numTotalParams-1] += 2*(cons.epsK*cons.numConfigs)*(currentParams[cons.numTotalParams-1]-epsDef);	
	}

	return 0;
}

int error_sort(int simplexSize, vector<long double> simplexErrors, vector<int> &errorRankToRow)
{	
	// Bubble sort
	bool swapPerformed = 1;
	while (swapPerformed == 1)
	{
		swapPerformed = 0;
		for (int i=0; i<(simplexSize-1); i++)
		{
			// Compare errors of point i and i+1
			if ( simplexErrors[ errorRankToRow[i] ] > simplexErrors[ errorRankToRow[i+1] ] )
			{
				int temp = errorRankToRow[i];
				errorRankToRow[i] = errorRankToRow[i+1];
				errorRankToRow[i+1] = temp;
				swapPerformed = 1;
			}
		}
	}
	
	return 0;
}

