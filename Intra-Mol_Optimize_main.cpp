
#include "IMO_header.h"

int read_input_params(constant_struct &cons, vector_struct &vecs, fitting_param_struct &initialParams)
{
	// Get filenames from inputFileString
	ifstream inputNames;
	cout << "Opening input file " << cons.inputFileString << endl;
	inputNames.open(cons.inputFileString.c_str());
	if (!inputNames.is_open()) {
		exit(EXIT_FAILURE);
	}
	
	// Get parameter file name
	inputNames >> cons.parameterFile;
	cout << "Reading parameters from " << cons.parameterFile << endl;
	
	ifstream paramStream;
	paramStream.open(cons.parameterFile.c_str());
	
	input_data_struct inputDefaults[] = {
		{"gen_xyz_file_names", TYPE_INT, &(cons.genXyzFileNames), "1"},
		{"phi_dim", TYPE_INT, &(cons.phiDim), "2"},
		{"phi1_range", TYPE_FLOAT_3VEC, &(cons.phi1Range), "0.0 180.0 10.0"},
		{"phi2_range", TYPE_FLOAT_3VEC, &(cons.phi2Range), "0.0 360.0 10.0"},
		{"phi3_range", TYPE_FLOAT_3VEC, &(cons.phi3Range), "0.0 360.0 10.0"},
		{"num_phi_surfaces", TYPE_INT, &(cons.numPhiSurfaces), "1"},
		{"kt_obj_func", TYPE_FLOAT, &(cons.KT), "20.92"},
		{"highest_order_rb_term", TYPE_INT, &(cons.nRBfit), "6"},
		{"scale_1_4_qq", TYPE_FLOAT, &(cons.scale14QQ), "0.5"},
		{"scale_1_4_lj", TYPE_FLOAT, &(cons.scale14LJ), "0.0"},
		{"LJ_comb_rule", TYPE_INT, &(cons.gromacsCombRule), "2"},
		{"eps_restraint_k", TYPE_FLOAT, &(cons.epsK), "1.0"},
		{"intra_mol_exclusion", TYPE_INT, &(cons.nrexcl), "3"},
		{"use_nwchem_suffix", TYPE_INT, &(cons.useNWsuffix), "1"},
		{"weight_edges_periodically", TYPE_INT, &(cons.weightPhiEdges), "1"},
		{"sigma_1_4_factor", TYPE_FLOAT, &(cons.sigma14factor), "1.0"},
		{"sigma_1_4_mod_cutoff", TYPE_FLOAT, &(cons.sigma14factorCutoff), "0.3"},
		{"use_adaptive_nelder_mead", TYPE_INT, &(cons.useAdaptiveNelderMead), "1"},
		{"write_energy_distribution", TYPE_INT, &(cons.writeEnergyDist), "0"}
	};
	
	int inputDefaultSize = sizeof(inputDefaults)/sizeof(inputDefaults[0]);
		
	// Set defaults
	for (int p=0; p<inputDefaultSize; p++) {
		char defaultLine[WORD_STRING_SIZE] = {0};
		sprintf(defaultLine, "%s %s", inputDefaults[p].keyword, inputDefaults[p].defString);
		process_input_line(&defaultLine[0], inputDefaults, inputDefaultSize, 0);
	}
	
	// Read parameters from file
	int nLines=0;
	string fLine;
	
    while(getline(paramStream, fLine)) {
		nLines++;
		process_input_line(fLine, inputDefaults, inputDefaultSize, 1);
		fLine[0] = '\0';
    }
	
	cons.pi = acosl(-1.0L);
		
	cons.xyzAngstroms = true;

	if (cons.sigma14factor != 1.0) { 
		cout << "Warning: Sigma-1,4 factor = " << cons.sigma14factor << ". Press any key to continue.";
		std::cin.ignore();
	}
	
	// Boltzmann conformer integral restraints (only works in 2D)
	cons.useBoltzIntRes = 0;
	cons.kTBoltzIntegral = 2.5; // Room temp
	//cons.kTBoltzIntegral = 3.3258; // 400K
	//cons.kTBoltzIntegral = 9.977292; // 1200K
	cons.kBoltzRes = 0.0;
	
	// Restraining params to default values
	cons.resToZero = 0;
	cons.rbK = 0.0;
	
	// 
	int numPhi1 = 0, numPhi2 = 0, numPhi3 = 0;
	
	if (cons.phiDim == 2 && cons.genXyzFileNames) {
		numPhi1 = round((cons.phi1Range[1]-cons.phi1Range[0])/cons.phi1Range[2] + 1);
		numPhi2 = round((cons.phi2Range[1]-cons.phi2Range[0])/cons.phi2Range[2] + 1);
		numPhi3 = 0;
		cons.numConfigs = numPhi1*numPhi2;		
	}
	else if (cons.phiDim == 3 && cons.genXyzFileNames) {
		numPhi1 = round((cons.phi1Range[1]-cons.phi1Range[0])/cons.phi1Range[2] + 1);
		numPhi2 = round((cons.phi2Range[1]-cons.phi2Range[0])/cons.phi2Range[2] + 1);
		numPhi3 = round((cons.phi3Range[1]-cons.phi3Range[0])/cons.phi3Range[2] + 1);
		cons.numConfigs = numPhi1*numPhi2*numPhi3;
		assert(cons.useBoltzIntRes==0); // Not currently supported 
	}
	else if (!cons.genXyzFileNames) {	
	}
	cout << "Num configs per surface = " << cons.numConfigs << endl;
	
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
		
		int numConformers = *std::max_element(vecs.partitionMap.begin(),vecs.partitionMap.end());
		vecs.confIntegralsDFT.resize(numConformers);
		std::fill(vecs.confIntegralsDFT.begin(),vecs.confIntegralsDFT.end(),0.0); // Ensure 0 filled
		
		// Integration rule for each conformer section
		vecs.integrationRule.resize((vecs.phi1partition.size()-1)*(vecs.phi2partition.size()-1));
		
		int i_section = 0; 
		for (int c1=1; c1<vecs.phi1partition.size(); c1++)
		{
			for (int c2=1; c2<vecs.phi2partition.size(); c2++)
			{				
				
				int m_max = round((vecs.phi1partition[c1]-cons.phi1Range[0])/cons.phi1Range[2]);
				int m_min = round((vecs.phi1partition[c1-1]-cons.phi1Range[0])/cons.phi1Range[2]);
				int n_max = round((vecs.phi2partition[c2]-cons.phi2Range[0])/cons.phi2Range[2]);
				int n_min = round((vecs.phi2partition[c2-1]-cons.phi2Range[0])/cons.phi2Range[2]);
				
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
				i_section++;
			}
		} 
	}
	
	// Molecule specific parameters (overwrite numConfigs if needed)
	if (cons.inputFileString == "input_files_EC.txt") {
		cons.resToZero = 1;
		cons.numConfigs = 869;
		cons.rbK = 1.0E-6;
	}
	if (cons.inputFileString == "input_files_PC.txt") {
		cons.resToZero = 1;
		cons.numConfigs = 2976;
		cons.rbK = 1.0E-6;
	}
	
	// Connectivity data sizes
	cons.atomDataSize = 4;
	cons.bondDataSize = 4;
	cons.angleDataSize = 5;
	cons.dihedralDataSize = 11;
	cons.improperDataSize = 7;
	cons.pairDataSize = 5;
	
	if (cons.genXyzFileNames) {
		vecs.connectFileList.resize(cons.numPhiSurfaces);
		vecs.xyzFileList.resize(cons.numPhiSurfaces);
		vecs.energyFileList.resize(cons.numPhiSurfaces);
	}
	else {
		assert(cons.numPhiSurfaces == 1);
		vecs.connectFileList.resize(1);
		vecs.xyzFileList.resize(1);
		vecs.energyFileList.resize(1);
	}
	
	initialParams.surfIndexMapping.resize(cons.numPhiSurfaces,-1);
	initialParams.dihedralIndexMapping.resize(512,-1); // Size = max number of dihedrals supported
	initialParams.pairIndexMapping.resize(512,-1);
	
	// Resize outer (i_s) loop of vector<vector> types
	vecs.xyzData.resize(cons.numPhiSurfaces);	
	vecs.dftData.resize(cons.numPhiSurfaces);	
	vecs.energyWeighting.resize(cons.numPhiSurfaces);	
	vecs.uConst.resize(cons.numPhiSurfaces);	
	
	vecs.molSize.resize(cons.numPhiSurfaces);	
	vecs.atomData.resize(cons.numPhiSurfaces);	
	vecs.bondData.resize(cons.numPhiSurfaces);	
	vecs.angleData.resize(cons.numPhiSurfaces);	
	vecs.dihedralData.resize(cons.numPhiSurfaces);	
	vecs.improperData.resize(cons.numPhiSurfaces);		
	vecs.pairData.resize(cons.numPhiSurfaces);	
	
	vecs.cosPsiData.resize(cons.numPhiSurfaces);	
	vecs.pairSepData.resize(cons.numPhiSurfaces);	
	
	vecs.energyDelta.resize(cons.numPhiSurfaces);

	vecs.bondSepMat.resize(cons.numPhiSurfaces);	
	vecs.rijMatrix.resize(cons.numPhiSurfaces);	
	vecs.sigmaMatrix.resize(cons.numPhiSurfaces);	
	vecs.epsilonMatrix.resize(cons.numPhiSurfaces);	
	vecs.qqMatrix.resize(cons.numPhiSurfaces);	
		
	for (int i_s=0; i_s<cons.numPhiSurfaces; i_s++) {
		// Connectivity file
		inputNames >> vecs.connectFileList[i_s];
		cout << "Connectivity file set to: " << vecs.connectFileList[i_s] << endl;
		connectivity_read(cons, vecs, i_s);
		
		// Process connectivity data
		connectivity_process(cons, vecs, i_s);
		
		// XYZ files (must be read before energy)
		inputNames >> vecs.xyzFileList[i_s];
		cout << "XYZ file set to: " << vecs.xyzFileList[i_s] << endl;
		xyz_files_read(cons, vecs, i_s);

		// Energy files
		inputNames >> vecs.energyFileList[i_s];
		cout << "Energy file set to: " << vecs.energyFileList[i_s] << endl;
		energy_read(cons, vecs, i_s);
		
		// Process constant terms in energy
		constant_energy_process(cons, vecs, initialParams, i_s);
	}
	
	// Hardcoded params
	//
	cons.annealWrite = 100;	// Frequency of writing lowest error to file
	cons.downhillWrite = 10;
	
	// Simulated annealing parameters
	cons.vTempInitial = (5000.0 + 50.0*cons.KT)/500.0; // Starting guess
	cons.vTempFactor = 0.9999;
	cons.dihedralStep = 2.0;
	cons.sigmaStep = 0.0002;
	cons.epsilonStep = 0.001;
	
	cons.sigmaGradFactor = 1.0;
		
	return 0;
}

int main(int argc, char *argv[])
{
	int annealIt, simplexIt, downhillIt; // Number of iterations for sim. annealing, downhill simplex & conj. grad.
	
	// Initialise data structures
	constant_struct cons; // Constants 
	vector_struct vecs;   // Vectors
	fitting_param_struct initialParams;
	fitting_param_struct currentParams;
	
	if (argc > 4) {
		annealIt = atof(argv[1]);
		simplexIt = atof(argv[2]);
		downhillIt = atof(argv[3]);
		cons.inputFileString = argv[4];
	}
	else {
		cout << "Not enough command line arguments!\n\n";
		return 1;	
	}
	
	// Read input and setup parameter struct
	read_input_params(cons, vecs, initialParams);
	currentParams = initialParams;
	
	cons.numTotalParams = currentParams.uShift.size() + currentParams.rbCoeff.size() + currentParams.ljSigma.size()*2;
	
	cout << "Starting parameters\n";
	print_params_console(cons, currentParams);
	
	srand(time(NULL));
	
	// -------- DO OPTIMISATION ------------------------------------------------------------------------------- 
	
	simulated_annealing(cons, vecs, initialParams, currentParams, annealIt);
	
	downhill_simplex(cons, vecs, initialParams, currentParams, simplexIt);
	
	conjugate_gradient(cons, vecs, initialParams, currentParams, downhillIt);
	
	// ---------------------------------------------------------------------------------------------------------
	// Write best energy
	error_from_trial_point(cons, vecs, initialParams, currentParams, 1);
	
}

int process_input_line(string fLine, input_data_struct* inputDefaults, int inputDefaultSize, bool toPrint)
{
	for (int l=0; l<inputDefaultSize; l++) {
		// Search for match
		string inputType;
		std::stringstream fLineStream(fLine);
		fLineStream >> inputType;
		
		if (!inputType.compare((inputDefaults+l)->keyword)) { // !NULL
			if (toPrint) {
				printf("%s %s\n", "Found input file line ", (inputDefaults+l)->keyword);
			}
			// Cast as appropriate data types 
			if ((inputDefaults+l)->dataType == TYPE_INT) {
				int valueInt;
				sscanf(fLine.c_str(), "%*s %d", &valueInt);
				memcpy( (inputDefaults+l)->varPtr, &valueInt, sizeof(int));
			}
			else if ((inputDefaults+l)->dataType == TYPE_FLOAT) {
				double valueFloat;
				sscanf(fLine.c_str(), "%*s %lf", &valueFloat);
				memcpy( (inputDefaults+l)->varPtr, &valueFloat, sizeof(double));
			}
			else if ((inputDefaults+l)->dataType == TYPE_INT_3VEC) {
				int value3VecInt[3];
				sscanf(fLine.c_str(), "%*s %d %d %d", value3VecInt, value3VecInt+1, value3VecInt+2);
				memcpy( (inputDefaults+l)->varPtr, value3VecInt, sizeof(double)*3);
			}
			else if ((inputDefaults+l)->dataType == TYPE_FLOAT_3VEC) {
				double value3VecFloat[3];
				sscanf(fLine.c_str(), "%*s %lf %lf %lf", value3VecFloat, value3VecFloat+1, value3VecFloat+2);
				memcpy( (inputDefaults+l)->varPtr, value3VecFloat, sizeof(double)*3);
			}
			else if ((inputDefaults+l)->dataType == TYPE_STRING) {
				char valueString[WORD_STRING_SIZE];
				sscanf(fLine.c_str(), "%*s %s", valueString);
				memcpy( (inputDefaults+l)->varPtr, valueString, sizeof(valueString) );
			}
		}
	}
	
	return 0;
}

int simulated_annealing(constant_struct cons, vector_struct vecs, fitting_param_struct &initialParams, fitting_param_struct &currentParams, int annealIt)
{
	int printFreq = 100; // Print simplex to console every printFreq steps
	
	srand(time(NULL));
	double temperature = cons.vTempInitial;
	
	ofstream annealStream;
	annealStream.open("annealing_energy.txt");
	annealStream.precision(12);
	bool toPrint, toWrite;
	
	fitting_param_struct oldParams = currentParams;
	fitting_param_struct annealParams = currentParams;
	
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0,1.0);
	
	for (int iT=0; iT < annealIt; iT++) {
		//
		toPrint = (iT%printFreq == 0);
		toWrite = (iT%cons.annealWrite == 0);
		
		if (toPrint) {
			cout << "Simulated annealing iteration " << iT << ", temp " << temperature << "\n\n";
		}
			
		// Generate trial point
		//cout << "Trial dihedral params:\n";
		for (int i_s=0; i_s < oldParams.uShift.size(); i_s++) {
			//
			double normalStep = distribution(generator);
			annealParams.uShift[i_s] = oldParams.uShift[i_s]
				+ normalStep*sqrt(temperature/cons.vTempInitial)*cons.dihedralStep;
		}
		
		for (int p=0; p < oldParams.rbCoeff.size(); p++) {
			//
			double normalStep = distribution(generator);			
			annealParams.rbCoeff[p] = oldParams.rbCoeff[p] 
				+ normalStep*sqrt(temperature/cons.vTempInitial)*cons.dihedralStep;
			//cout << annealParams.rbCoeff[p] << endl;
		}
		for (int pair=0; pair < oldParams.ljEps.size(); pair++) {
			//
			annealParams.ljSigma[pair] = std::abs(oldParams.ljSigma[pair] 
				+ cons.sigmaStep*sqrt(temperature/cons.vTempInitial)*distribution(generator));
			annealParams.ljEps[pair] = std::abs(oldParams.ljEps[pair] 
				+ cons.epsilonStep*sqrt(temperature/cons.vTempInitial)*distribution(generator));
			//cout << "Trial (sigma,epsilon) = (";
			//cout << annealParams.ljSigma[pair] << ", " << annealParams.ljEps[pair] << ")\n";
		}
		
		// Test neighbour parameters
		double oldError = error_from_trial_point(cons, vecs, initialParams, oldParams, 0);
		annealStream << oldError << endl;
		double trialError = error_from_trial_point(cons, vecs, initialParams, annealParams, 0);
		
		// Compute acceptance probability
		double aProb = exp((oldError-trialError)/temperature);
		
		if (toPrint)
		{
			cout << "Old error   = " << oldError << endl;
			cout << "Trial error = " << trialError << endl;
			cout << "Acceptance probability = " << aProb << "\n\n";
		}	
		if (aProb > (rand()/double(RAND_MAX))) {
			if (toPrint) {
				cout << "Accepting these params.\n";
				print_params_console(cons, annealParams);
			}
			oldParams = annealParams;
		}
		else if (toPrint) {
			cout << "Rejecting these params! \n\n";
		}
		temperature *= cons.vTempFactor;
	}
	
	// Update final parameters
	currentParams = oldParams;
	
	return 0;
}

int downhill_simplex(constant_struct cons, vector_struct vecs, fitting_param_struct &initialParams, fitting_param_struct &currentParams, int simplexIt)
{
	cout << "\nPerforming " << simplexIt << " downhill simplex iteraitons\n\n";
	double nmEPS = 1.0E-14L;
	int printFreq = 100; // Print simplex to console every printFreq steps
	int numReflectSteps = 0;
	
	int simpSize = cons.numTotalParams + 1;
	cout << "Simplex size: " << simpSize << endl;
	
	// Create simplex
	simplex_struct sim;
	
	// Choose parameters 
	if (cons.useAdaptiveNelderMead) { 
		// My parameters
		sim.nmReflect = 1.0;
		sim.nmExpand = 1.0 + 2.0/double(cons.numTotalParams);
		sim.nmContra = 0.9 - 0.8/double(cons.numTotalParams);
		sim.nmShrink = 1.0 - 1.0/double(cons.numTotalParams);
	}
	else {
		// Default parameters 
		sim.nmReflect = 1.0;
		sim.nmExpand = 2.0;
		sim.nmContra = 0.8;
		sim.nmShrink = 0.5;
	}

	// Initialize simplex points from current params
	sim.plex.resize(simpSize);
	for (int row=0; row<simpSize; row++) {
		sim.plex[row] = currentParams; 
	}
	sim.plexInitialize();
	
	// Print starting simplex
	cout << "--- INITIAL SIMPLEX ---\n";
	print_simplex_console(cons, sim);
	
	// Analyse simplex and fill simplexErrors
	for (int row=0; row<simpSize; row++) {
		sim.plexErrors[row] = error_from_trial_point(cons, vecs, initialParams, sim.plex[row], 0);
		sim.errorRankToRow[row] = row;
		cout << "Row: " << row << ", error " << sim.plexErrors[row] << endl;
	}
	
	for (int iD=0; iD < simplexIt; iD++) {
		//
		int toPrint = (iD%printFreq == 0);	
		if (toPrint) cout << "\nIteration: " << iD << endl << endl;

		// Sort errors, determine max, min and write best
		error_sort(sim.plexErrors, sim.errorRankToRow);		
		int bestRow = sim.errorRankToRow[0];
		int worstRow = sim.errorRankToRow[simpSize-1];
		double minError = sim.plexErrors[bestRow];
		double maxError = sim.plexErrors[worstRow];
		double secondMaxError = sim.plexErrors[sim.errorRankToRow[simpSize-2]];	
		
		// Return if end condition met
		if (std::fabs(maxError-minError)/maxError < nmEPS) {
			cout << "\nEnding downhill simplex routine after " << iD << " iterations\n\n";
			cout << "The percentage of reflect steps was " << 100.0*double(numReflectSteps)/double(iD) << endl;
			
			if (simplexIt > 0) {
				// Update currentParams array
				cout << "\n Best params from simplex:\n";			
				print_params_console(cons, sim.plex[bestRow]);	
				currentParams = sim.plex[bestRow];					
			}
			return 0;		
		}
		
		// Print best
		if (toPrint) cout << "MAX and MIN errors: " << std::setprecision(16) << maxError << " " << minError << endl;
		
		// Compute centroid of all points execpt worst		
		sim.computeCentroid();
		
		// Compute reflected point and output best params
		if (toPrint) cout << "\nComputing reflected point in opposite direction to point " << worstRow << endl;

		sim.computeReflect();
			
		double reflectError = error_from_trial_point(cons, vecs, initialParams, sim.reflectPoint, 0);
		
		// If neither new best nor worst (or 2nd worst), replace worst row of simplex with trial point
		if ((reflectError <= secondMaxError) && (reflectError >= minError)) {
			if (toPrint) {
				cout << "Reflected point neither best nor worst.\n";
				cout << "Replacing worst point with reflected point.\n";	
			}		
				
			sim.plex[worstRow] = sim.reflectPoint;	
			sim.plexErrors[worstRow] = reflectError;
			numReflectSteps++;
		}
		
		// If best, test with expansion in direction of reflected point
		else if (reflectError < minError) {
			if (toPrint) {
				cout << "Reflected point is better than any point in simplex.\n";
				cout << "Generating new point by expansion.\n\n";
			}
						
			sim.computeExpand();
			
			// Test expanded point
			double expandError = error_from_trial_point(cons, vecs, initialParams, sim.expandPoint, 0);
			
			if (expandError < minError) {
				if (toPrint) {
					cout << "Expanded point even better.\n";
					cout << "Replacing worst point with expanded point.\n\n";
				}
											
				sim.plex[worstRow] = sim.expandPoint;
				sim.plexErrors[worstRow] = expandError;
			}
			else {
				if (toPrint) {
					cout << "Expansion unsuccessful.\n";
					cout << "Replacing worst point with reflected point.\n\n";
				}			
				
				sim.plex[worstRow] = sim.reflectPoint;
				sim.plexErrors[worstRow] = reflectError;
				numReflectSteps++;			
			}
		}
		else {
			// Else, the reflected point error must be > the second worst point
			if (toPrint) {
				cout << "Reflection unsuccessful, performing contraction.\n";
			}
			//Replace worst point by reflected point if it's an improvement before contraction
			if (reflectError < maxError) {
				
				sim.plex[worstRow] = sim.reflectPoint;
				sim.plexErrors[worstRow] = reflectError;
				maxError = reflectError;
			}
			
			sim.computeContract();
						
			// Test contracted point
			double contractError = error_from_trial_point(cons, vecs, initialParams, sim.contraPoint, 0);
			
			if (contractError < maxError) {
				if (toPrint) {
					cout << "Replacing worst point with contracted point.\n\n";
				}
				
				sim.plex[worstRow] = sim.contraPoint;
				sim.plexErrors[worstRow] = contractError;
			}
			else {
				// Reduction towards best point
				if (toPrint) {
					cout << "Contracted point is new worst point. Trying reduction towards best point.\n\n";
				}
				for (int row=0; row<simpSize; row++) {
					
					if (row != bestRow) {
						sim.reducePoint(row);
						sim.plexErrors[row] = error_from_trial_point(cons, vecs, initialParams, sim.plex[row], 0);
					}					
				}
			}	
		}
		// Output new simplex
		if (toPrint) {
			print_simplex_console(cons, sim);
		}
	} // Iteration loop
	
	if (simplexIt > 0) {
		// Update currentParams array
		error_sort(sim.plexErrors, sim.errorRankToRow);
		currentParams = sim.plex[sim.errorRankToRow[0]];
	}
		
	cout << "\nEnding downhill simplex routine after max iterations\n\n";
	
	return 0;
}  

int conjugate_gradient(constant_struct cons, vector_struct vecs, fitting_param_struct initialParams, 
	fitting_param_struct &currentParams, int gradientIt)
{
	int printFreq = 100; // Print every printFreq steps
	cout << "\nPerforming " << gradientIt << " conjugate gradient iteraitons\n\n";
	
	// Initialize gradient param struct
	fitting_param_struct gradVector = currentParams;
	compute_gradient_F(cons, vecs, initialParams, currentParams, gradVector);	
	fitting_param_struct gradVectorNew = gradVector;
	fitting_param_struct gradDelta = gradVector;
	fitting_param_struct cgVector = gradVector;
	fitting_param_struct cgVectorNew = gradVector;
	
	fitting_param_struct tempParams = currentParams;
	
	// Do initial steepest descent step with line search
	double currentF;
	double previousF = error_from_trial_point(cons, vecs, initialParams, currentParams, 0);
	
	double alpha = 0.001, tau = 0.5; // Starting step size and reduction factor
	double currentBestAlpha = alpha;
		
	const int nMax = 32; // Maximum number of iterations in line search
	double lineSearchF[nMax];
	double currentBestF = previousF;
		
	// Perform initial line search
	for (int n=0; n<nMax; n++) {
		param_linear_combine(tempParams, currentParams, gradVector, 1.0, -alpha);
		
		// Write corresponding value of F
		lineSearchF[n] = error_from_trial_point(cons, vecs, initialParams, tempParams, 0);
		
		if (lineSearchF[n] < currentBestF) {
			currentBestF = lineSearchF[n];
			currentBestAlpha = alpha;
		}
		alpha *= tau;
	}
	
	cgVector = gradVector;
	param_linear_combine(currentParams, currentParams, gradVector, 1.0, -currentBestAlpha);

	previousF = error_from_trial_point(cons, vecs, initialParams, currentParams, 0);
	
	// Do modified conjugate gradient steps ----------------------------------------------------------------------------
	for (int iG=0; iG < gradientIt; iG++) {
		//
		bool toPrint = (iG%printFreq == 0);
			
		// Get gradient corresponding to x_n, write to gradVectorNew
		compute_gradient_F(cons, vecs, initialParams, currentParams, gradVectorNew);	
		
		if (toPrint) {
			cout << "\nGradient:\n";
			print_params_console(cons, gradVectorNew);
		}
		
		// Calculate dot-products needed for beta and theta coefficients
		param_linear_combine(gradDelta, gradVectorNew, gradVector, 1.0, -1.0);
		double dotProd1 = param_scalar_product(gradVectorNew, gradDelta);
		double dotProd2 = param_scalar_product(gradVector, gradVector);
		double dotProd3 = param_scalar_product(gradVectorNew, cgVector);
			
		// Calculate coefficient beta (Polak-Ribiere)
		double betaPR = std::max(dotProd1/dotProd2, 0.0);
		//cout << "betaPR = " << betaPR << endl;
		
		// Calculate coefficient theta (Zhang)
		double thetaZ = dotProd3/dotProd2;	
		//cout << "thetaZ = " << thetaZ << endl << endl;
		
		// Update conjugate direction, advance vectors n -> n+1
		param_linear_combine(cgVectorNew, gradVectorNew, cgVector, -1.0, betaPR);
		param_linear_combine(cgVectorNew, cgVectorNew, gradDelta, 1.0, -thetaZ);
		gradVector = gradVectorNew;
		cgVector = cgVectorNew;
		double magCgVec = param_scalar_product(cgVector, cgVector);
		
		// Perform line search in conjugate direction
		int n=0;
		alpha = 0.001;
		double deltaZ = 0.0001;
		double deltaCG = -1.0; // Set negative for start of while loop
		while ( (n<=nMax) && (deltaCG < 0.0) ) {
			// Get updated params from x_i+1 = x_i - alpha*d_i+1
			param_linear_combine(tempParams, currentParams, cgVectorNew, 1.0, alpha);

			currentF = error_from_trial_point(cons, vecs, initialParams, tempParams, 0);	
			// Test Armijo-type condition 
			deltaCG = previousF - currentF - deltaZ*alpha*alpha*magCgVec; // should be >= 0 to continue
			//cout << alpha << ": previousF " << std::setw(10) << previousF << "  currentF " << currentF;
			//cout << "  delta*alpha^2*d^2" << std::setw(10) << deltaZ*alpha*alpha*magCgVec;
			//cout << "  CGdelta " << std::setw(10) << deltaCG << endl;
			alpha *= tau;
			n++;
		}
		
		// Overwrite curent parameters
		previousF = currentF;
		currentParams = tempParams;
		
		if (toPrint) {
			cout << "Updating parameters, F = " << currentF << endl;
			print_params_console(cons, currentParams);
		}
	}		
	
	return 0;
}  

int connectivity_read(constant_struct &cons, vector_struct &vecs, int i_s)
{
	ifstream connectStream;	
	cout << "Opening " << vecs.connectFileList[i_s] << endl;
	
	connectStream.open(vecs.connectFileList[i_s]);
	
	// Read size of this molecule (num atoms, bonds, etc.)
	string connectLine;
	int section = -1;
	while(getline(connectStream, connectLine)) {
		if (connectLine.find("atoms") != string::npos)
			section = 0;
		else if (connectLine.find("bonds") != string::npos) {
			assert (section == 0);
			section = 1;
		}
		else if (connectLine.find("angles") != string::npos) {
			assert (section == 1);
			section = 2;
		}
		else if (connectLine.find("dihedrals") != string::npos) {
			assert (section == 2);
			section = 3;
		}
		else if (connectLine.find("impropers") != string::npos) {
			assert (section == 3);
			section = 4;
		}
		else if (connectLine.find("lj_pair_fit") != string::npos) {
			assert (section == 4);
			section = 5;
		}
		else if (!connectLine.empty() && section >= 0) {
			
			switch(section) {
				case 0 : vecs.molSize[i_s].atoms++; break;
				case 1 : vecs.molSize[i_s].bonds++; break;
				case 2 : vecs.molSize[i_s].angles++; break;
				case 3 : vecs.molSize[i_s].dihedrals++; break;
				case 4 : vecs.molSize[i_s].impDihedrals++; break;
				case 5 : vecs.molSize[i_s].pairs++; break;
			}
		}
	}
	cout << vecs.molSize[i_s].atoms << " atoms found.\n";
	cout << vecs.molSize[i_s].bonds << " bonds found.\n";
	cout << vecs.molSize[i_s].angles << " angles found.\n";
	cout << vecs.molSize[i_s].dihedrals << " dihedrals found.\n";
	cout << vecs.molSize[i_s].impDihedrals << " improper dihedrals found.\n";
	cout << vecs.molSize[i_s].pairs << " special pairs found.\n";
	
	// Resize inner loop of vector<vector> types for this molecule and energy surface
	vecs.xyzData[i_s].resize(cons.numConfigs*vecs.molSize[i_s].atoms*3);	
	vecs.dftData[i_s].resize(cons.numConfigs);
	vecs.energyWeighting[i_s].resize(cons.numConfigs);
	vecs.uConst[i_s].resize(cons.numConfigs);
	
	vecs.atomData[i_s].resize(vecs.molSize[i_s].atoms*cons.atomDataSize);
	vecs.bondData[i_s].resize(vecs.molSize[i_s].bonds*cons.bondDataSize);
	vecs.angleData[i_s].resize(vecs.molSize[i_s].angles*cons.angleDataSize);
	vecs.dihedralData[i_s].resize(vecs.molSize[i_s].dihedrals*cons.dihedralDataSize);
	
	if (vecs.molSize[i_s].impDihedrals != 0){
		vecs.improperData[i_s].resize(vecs.molSize[i_s].impDihedrals*cons.dihedralDataSize);
	}
	if (vecs.molSize[i_s].pairs != 0) {
		vecs.pairData[i_s].resize(vecs.molSize[i_s].pairs*cons.pairDataSize);
	}
	
	vecs.energyDelta[i_s].resize(cons.numConfigs);		
	vecs.bondSepMat[i_s].resize(vecs.molSize[i_s].atoms*vecs.molSize[i_s].atoms);
	vecs.rijMatrix[i_s].resize(vecs.molSize[i_s].atoms*vecs.molSize[i_s].atoms);
	vecs.sigmaMatrix[i_s].resize(vecs.molSize[i_s].atoms*vecs.molSize[i_s].atoms);
	vecs.epsilonMatrix[i_s].resize(vecs.molSize[i_s].atoms*vecs.molSize[i_s].atoms);
	vecs.qqMatrix[i_s].resize(vecs.molSize[i_s].atoms*vecs.molSize[i_s].atoms);
	
	string fileSection;
	connectStream.close();
	connectStream.open(vecs.connectFileList[i_s]);
	
	// Read atoms (with charges and van der Waals parameters)
	connectStream >> fileSection;
	if (fileSection == "atoms") {
		for (int a=0; a<vecs.molSize[i_s].atoms; a++) {
			for (int i=0; i<cons.atomDataSize; i++) {
				connectStream >> vecs.atomData[i_s][a*cons.atomDataSize + i];
			}
		}
	}
	else {
		cout << "Error in input file format (no atoms section)\n\n";
		return 1;
	}
	
	// Read bonds
	connectStream >> fileSection;
	if (fileSection == "bonds") {
		for (int b=0; b<vecs.molSize[i_s].bonds; b++) {
			for (int i=0; i<cons.bondDataSize; i++) {
				connectStream >> vecs.bondData[i_s][b*cons.bondDataSize + i];
			}
		}
	}
	else {
		cout << "Error in input file format (no bonds section)\n\n";
		return 1;
	}
	
	// Read angles
	connectStream >> fileSection;
	if (fileSection == "angles") {
		for (int c=0; c < vecs.molSize[i_s].angles; c++) {
			for (int i=0; i < cons.angleDataSize; i++) {
				connectStream >> vecs.angleData[i_s][c*cons.angleDataSize + i];
			}
		}
	}
	else {
		cout << "Error in input file format (no angles section)\n\n";
		return 1;
	}
	
	// Read dihedrals
	connectStream >> fileSection;
	if (fileSection == "dihedrals") {
		for (int d=0; d < vecs.molSize[i_s].dihedrals; d++) {
			for (int i=0; i < cons.dihedralDataSize; i++) { 
				connectStream >> vecs.dihedralData[i_s][d*cons.dihedralDataSize + i];		
			}
		}
	}
	else {
		cout << "Error in input file format (no dihedrals section)\n\n";
		return 1;
	}
	
	// Read harmonic improper dihedrals
	connectStream >> fileSection;
	if (fileSection == "impropers") {
		for (int id=0; id < vecs.molSize[i_s].impDihedrals; id++) {
			for (int i=0; i < cons.improperDataSize; i++) { 
				connectStream >> vecs.improperData[i_s][id*cons.improperDataSize + i];
			}
		}
	}
	else {
		cout << "Error in input file format (no impropers section)\n\n";
		return 1;
	}
	
	// Read lj pair to fit
	connectStream >> fileSection;
	if (fileSection == "lj_pair_fit") {
		for (int e = 0; e < vecs.molSize[i_s].pairs; e++) {
			for (int i = 0; i < cons.pairDataSize; i++) {
				connectStream >> vecs.pairData[i_s][e*cons.pairDataSize + i];
			}
		}
	}
	else {
		cout << "Error in input file format (no lj_pair_fit section)\n\n";
		return 1;
	}
	
	return 0;
}

int connectivity_process(constant_struct cons, vector_struct &vecs, int i_s)
{
	// Read atoms separated by one bond
	for (int i=0; i < vecs.molSize[i_s].atoms; i++)
	{
		for (int j=i; j < vecs.molSize[i_s].atoms; j++)
		{
			int m1 = i*vecs.molSize[i_s].atoms + j; // Convert to 1D array coordinate
			int m2 = j*vecs.molSize[i_s].atoms + i; // Symmetric matrix
			
			int atomI = int(round(vecs.atomData[i_s][i*cons.atomDataSize])); // Get actual atom number from file
			int atomJ = int(round(vecs.atomData[i_s][j*cons.atomDataSize]));
			
			// Search bonds list
			bool bondFound = 0;
			for (int b=0; b<vecs.molSize[i_s].bonds; b++) {				
				if ( (vecs.bondData[i_s][b*cons.bondDataSize] == atomI) && (vecs.bondData[i_s][b*cons.bondDataSize + 1] == atomJ) ) {
					vecs.bondSepMat[i_s][m1] = 1;
					vecs.bondSepMat[i_s][m2] = 1;
					bondFound = 1;
				}
			}
			if (bondFound == 0) {
				vecs.bondSepMat[i_s][m1] = 0;
				vecs.bondSepMat[i_s][m2] = 0;
			}			
		}
	}
	
	// Now loop over increasing bond separation number until array is populated
	int matUpdate = 1;
	int bondSearch = 2;
	while (matUpdate == 1) {
		matUpdate = 0; // This will be set to 1 if an update was made
		
		for (int i=0; i < vecs.molSize[i_s].atoms; i++) {
			for (int j=(i+1); j < vecs.molSize[i_s].atoms; j++) { // j>i, since i=j are the same atom
				int m1 = i*vecs.molSize[i_s].atoms + j; // Convert to 1D array coordinate
				int m2 = j*vecs.molSize[i_s].atoms + i; // Symmetric matrix
			
				// If number of bonds between pair i-j is not already determined,
				// look for atoms which are bondSearch-1 from atom i				
				if (vecs.bondSepMat[i_s][m1] == 0) {
					
					for (int k=0; k < vecs.molSize[i_s].atoms; k++) {
						
						int m3 = k*vecs.molSize[i_s].atoms + j;
						// If connection (k,j) is bondSearch-1, check if (k,i) is 1
						if (vecs.bondSepMat[i_s][m3] == (bondSearch - 1)) {
								
							int m4 = k*vecs.molSize[i_s].atoms + i;
							
							if(vecs.bondSepMat[i_s][m4] == 1) {
								vecs.bondSepMat[i_s][m1] = bondSearch;
								vecs.bondSepMat[i_s][m2] = bondSearch; 
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
	for (int i=0; i < vecs.molSize[i_s].atoms; i++) {
		for (int j=0; j < vecs.molSize[i_s].atoms; j++) {
			cout << vecs.bondSepMat[i_s][i*vecs.molSize[i_s].atoms+j] << " ";
		}
		cout << endl;
	}
	
	// Work out sigma, epsilon and chargeQQ for each pair
	for (int i=0; i < vecs.molSize[i_s].atoms; i++) {
		for (int j=i; j < vecs.molSize[i_s].atoms; j++) {
			
			int m1 = i*vecs.molSize[i_s].atoms+j;
			int a1 = int(round(vecs.atomData[i_s][i*cons.atomDataSize]));
			int a2 = int(round(vecs.atomData[i_s][j*cons.atomDataSize]));
			//int m2 = j*vecs.molSize[i_s].atoms+i;
			bool toFit = 0;
			
			// Loop over pairs to work out if this is the LJ pair to be fitted
			for (int p=0; p < vecs.molSize[i_s].pairs; p++) {
				if ( (vecs.pairData[i_s][p*cons.pairDataSize] == a1) && (vecs.pairData[i_s][p*cons.pairDataSize + 1] == a2) )
					toFit = 1;		
			}
			
			assert(cons.gromacsCombRule == 2);
			
			double sigma_i = vecs.atomData[i_s][i*cons.atomDataSize+2];
			double sigma_j = vecs.atomData[i_s][j*cons.atomDataSize+2];
			
			// Apply sigma sclaing factor for 1-4 interactions	
			if (vecs.bondSepMat[i_s][m1] == 3) {
				sigma_i = sigma_i > cons.sigma14factorCutoff ? 
					sigma_i*cons.sigma14factor : sigma_i;
				sigma_j = sigma_j > cons.sigma14factorCutoff ? 
					sigma_j*cons.sigma14factor : sigma_j;
			}	
			
			double sigma = 0.5*(sigma_i+ sigma_j);
			double epsilon = sqrt(vecs.atomData[i_s][i*cons.atomDataSize+3]*vecs.atomData[i_s][j*cons.atomDataSize+3]);
			double QQ = vecs.atomData[i_s][i*cons.atomDataSize+1]*vecs.atomData[i_s][j*cons.atomDataSize+1];
			
			// Excluded interactions
			if (vecs.bondSepMat[i_s][m1] <= cons.nrexcl) {
				vecs.sigmaMatrix[i_s][m1] = 0.0;
				vecs.epsilonMatrix[i_s][m1] = 0.0;
				vecs.qqMatrix[i_s][m1] = 0.0;
			}
			// 1-4 Interactions with scaling factor
			if (vecs.bondSepMat[i_s][m1] == 3) {
				vecs.sigmaMatrix[i_s][m1] = sigma;
				vecs.epsilonMatrix[i_s][m1] = cons.scale14LJ*epsilon;
				vecs.qqMatrix[i_s][m1] = cons.scale14QQ*QQ;
			}
			// Full interactions
			if (vecs.bondSepMat[i_s][m1] > cons.nrexcl) {
				vecs.sigmaMatrix[i_s][m1] = sigma;
				vecs.epsilonMatrix[i_s][m1] = epsilon;
				vecs.qqMatrix[i_s][m1] = QQ;
			}
			if (toFit == 1) {
				//Remove LJ interaction for pairs being fit (usually 1-5 pairs)
				vecs.epsilonMatrix[i_s][m1] = 0.0;
			}	
		}
	}
	
	cout << "sigmaLJ matrix complete:\n";
	for (int i = 0; i < vecs.molSize[i_s].atoms; i++) {
		for (int j = 0; j < vecs.molSize[i_s].atoms; j++) {
			cout << std::setw(6) << vecs.sigmaMatrix[i_s][i*vecs.molSize[i_s].atoms+j] << " ";
		}
		cout << endl;
	}
	cout << "epsilonLJ matrix complete:\n";
	for (int i = 0; i < vecs.molSize[i_s].atoms; i++) {
		for (int j = 0; j < vecs.molSize[i_s].atoms; j++) {
			cout << std::setw(6) << vecs.epsilonMatrix[i_s][i*vecs.molSize[i_s].atoms+j] << " ";
		}
		cout << endl;
	}
	cout << "chargeQQ matrix complete:\n";
	for (int i = 0; i < vecs.molSize[i_s].atoms; i++) {
		for (int j = 0; j < vecs.molSize[i_s].atoms; j++) {
			cout << std::setw(6) << vecs.qqMatrix[i_s][i*vecs.molSize[i_s].atoms+j] << " ";
		}
		cout << endl;
	}
	return 0;
}

int xyz_files_read(constant_struct &cons, vector_struct &vecs, int i_s)
{
	ifstream xyzStream;
	
	// Loop over full (usually 360 degree) range at specificed step size
	if (cons.genXyzFileNames) {
		int i_f=0;
		double phi3Min, phi3Max, phi3Step;
		
		if (cons.phiDim == 2){ 			
			// Set phi3 range so inner loop executed once
			cons.phi3Range[0] = 0.0; 
			cons.phi3Range[1] = 360.0;
			cons.phi3Range[2] = 10.0;
			phi3Min = 180.0; // Could be anything 0<x<360
			phi3Max = 180.0;
			phi3Step = 10.0; // anything >0 will work
		}
		else if (cons.phiDim == 3) {
			phi3Min = cons.phi3Range[0];
			phi3Max = cons.phi3Range[1];
			phi3Step = cons.phi3Range[2];
		}
			
		for (int d1=cons.phi1Range[0]; d1<=cons.phi1Range[1]; d1 = d1+cons.phi1Range[2])
		{
			for (int d2=cons.phi2Range[0]; d2<=cons.phi2Range[1]; d2 = d2+cons.phi2Range[2])
			{
				for (int d3=phi3Min; d3<=phi3Max; d3 = d3+phi3Step)
				{
					// Determine weighting to account for double counting of edge points
					// This assumes phi2 goes from 0->360, so the 0/360 is double counted
					if (cons.weightPhiEdges) {
					
						int numBoundaries = ((d1 == cons.phi1Range[0]) || (d1 == cons.phi1Range[1])) + 
							((d2 == cons.phi2Range[0]) || (d2 == cons.phi2Range[1])) + 
								((d3 == cons.phi3Range[0]) || (d3 == cons.phi3Range[1]));
					
						// 0.125 for vertices (3D)
						if (numBoundaries==3) {
							vecs.energyWeighting[i_s][i_f] = 0.125;
						}
						// 0.25 for edges (3D) corners (2D)
						else if (numBoundaries==2) {
							vecs.energyWeighting[i_s][i_f] = 0.25;
						}
						// 0.5 for faces (3D) or edges (2D)
						else if (numBoundaries==1) {
							vecs.energyWeighting[i_s][i_f] = 0.5;
						}
						// 1.0 for central points
						else {
							vecs.energyWeighting[i_s][i_f] = 1.0;	
						}
					}
					else {
						vecs.energyWeighting[i_s][i_f] = 1.0;	
					}
					//cout << "d1, d2, weighting: " << d1 << ", " << d2 << ", " << vecs.energyWeighting[i_s][i_f] << endl;
				
					bool file_open = 1;
					int ts = 0;
			
					string tsStrPrev, fullFileName;
					
					string d1Str = std::to_string(d1);
					string d2Str = std::to_string(d2);
					string d3Str = std::to_string(d3);
			
					if (cons.useNWsuffix) {
						// This method iterates over the -###.xyz suffix to find the last step in a geometry optimisation
						// Only needed if geometries are taken from NWChem optimisation, which they probably shouldn't
						while (file_open == 1) {
							// Generate file name
							string tsStr = std::to_string(ts);
							
							// time-step code is 3-digits
							if (ts < 10)
								tsStr = "00" + tsStr;
							else if (ts < 100)
								tsStr = "0" + tsStr;
							else if (ts > 999)
								cout << "Warning: timestep has 4 digits (undefined behaviour)\n";
							
							// Generate file name		
							if (cons.phiDim==3) {
								fullFileName = vecs.xyzFileList[i_s] + d1Str + "_" + d2Str + "_" + d3Str + "-" + tsStr + ".xyz";
							}
							else {
								fullFileName = vecs.xyzFileList[i_s] + d1Str + "_" + d2Str + "-" + tsStr + ".xyz";
							}
							xyzStream.open(fullFileName.c_str());
				
							// Once we cannot open the file, we read data from the previous file
							if (!xyzStream.is_open()) {
								file_open = 0;
								// Get previous file name
								if (cons.phiDim==3) {
									fullFileName = vecs.xyzFileList[i_s] + d1Str + "_" + d2Str + "_" + d3Str + "-" + tsStrPrev + ".xyz";
								}
								else {
									fullFileName = vecs.xyzFileList[i_s] + d1Str + "_" + d2Str + "-" + tsStrPrev + ".xyz";
								}
								
							}
							tsStrPrev = tsStr;
							xyzStream.close();
							ts++;
						}	
					}					
					else {
						if (cons.phiDim==3) {
							fullFileName = vecs.xyzFileList[i_s] + d1Str + "_" + d2Str + "_" + d3Str + ".xyz";
						}
						else {
							fullFileName = vecs.xyzFileList[i_s] + d1Str + "_" + d2Str + ".xyz";
						}
					}	
					
					//cout << "Opening xyz file " << fullFileName << endl;
					xyzStream.open(fullFileName.c_str());

					// Read but ignore first two lines
					string numAtoms, atom, xyzHeader;
					std::getline(xyzStream, numAtoms);
					std::getline(xyzStream, xyzHeader);
					
					for (int a=0; a<vecs.molSize[i_s].atoms; a++) {
						int k = (i_f*vecs.molSize[i_s].atoms+a)*3;
						xyzStream >> atom >> vecs.xyzData[i_s][k] >> vecs.xyzData[i_s][k+1] >> vecs.xyzData[i_s][k+2];
						//cout << "K = " << k << endl;
						//cout << vecs.xyzData[i_s][k] << " " << vecs.xyzData[i_s][k+1] << " " << vecs.xyzData[i_s][k+2] << endl;	
					}
					i_f++;
				}		
			}
		}
	}
	else {
		assert(cons.phiDim == 2);
		
		ifstream phiPairStream, xyzStream;
		cout << "Reading phi pairs from " << vecs.xyzFileList[0] << endl;
			
		phiPairStream.open(vecs.xyzFileList[0].c_str());
		if( !phiPairStream.is_open() ) {
			cout << "Error, phi pair file cannot be opened!" << endl;
			return 1;
		}

		for (int i_f=0; i_f<cons.numConfigs; i_f++) {	
			// Set weighting
			vecs.energyWeighting[i_s][i_f] = 1.0;
			
			// Read phi1 and phi2 strings
			char phi1Str[16], phi2Str[16];
	
			// Read phi pair strings from phiCoordFile
			// Only supports 2D for now
			double phi1Spec, phi2Spec;
			phiPairStream >> phi1Spec >> phi2Spec;
			
			snprintf(phi1Str, 16, "%.1f", phi1Spec);
			snprintf(phi2Str, 16, "%.1f", phi2Spec);
			
			// Generate xyz filename
			// This method assumes the files end with "-000.xyz"
			string fileName = vecs.xyzFileList[0] + phi1Str + "_" + phi2Str + "-000" + ".xyz";
			
			//cout << "Opening xyz file: " << fileName << endl;
			
			xyzStream.open(fileName.c_str());
			
			// Ignore first 2 lines and atom name (as per XYZ format)
			string numAtoms, atom, xyzHeader;
			std::getline(xyzStream, numAtoms);
			std::getline(xyzStream, xyzHeader);
		
			for (int a = 0; a < vecs.molSize[i_s].atoms; a++) {
				int k = (i_f*vecs.molSize[i_s].atoms+a)*3;
				xyzStream >> atom >> vecs.xyzData[i_s][k] >> vecs.xyzData[i_s][k+1] >> vecs.xyzData[i_s][k+2];
				//cout << vecs.xyzData[i_s][k] << " " << vecs.xyzData[i_s][k+1] << " " << vecs.xyzData[i_s][k+2] << endl;			
			}
			xyzStream.close();
		}
	}
	return 0;
}

int energy_read(constant_struct &cons, vector_struct &vecs, int i_s)
{
	ifstream energyStream;
	energyStream.open(vecs.energyFileList[i_s]);
	
	if(!energyStream.is_open()) {
		cout << "ERROR: Energy file was not opened!\n";
		exit(EXIT_FAILURE);
		return 1;
	}
	
	for (int i_u=0; i_u<cons.numConfigs; i_u++) {
		energyStream >> vecs.dftData[i_s][i_u];
		// Apply Boltzmann weighting 
		vecs.energyWeighting[i_s][i_u] *= exp(-vecs.dftData[i_s][i_u]/cons.KT); 
		//cout << "Energy, weighting: " << i_u << " " << vecs.dftData[i_s][i_u] << ", ";
		//cout << vecs.energyWeighting[i_s][i_u] << endl;
	}
	
	/* Integrate DFT energy over conformers
	if (cons.useBoltzIntRes == 1) {
		assert(cons.numPhiSurfaces == 1);
		
		int numSections = (vecs.phi1partition.size()-1)*(vecs.phi2partition.size()-1);	
		vector<double> sectionIntegrals(numSections, 0.0);
	
		// Calculate Boltzmann distribution integrals for each section
		int i_section = 0; 
		assert(vecs.integrationRule[i_section] == 1); // Until other integral methods supported
			
		for (int c1=1; c1<vecs.phi1partition.size(); c1++) {
			for (int c2=1; c2<vecs.phi2partition.size(); c2++) {
								
				int m_max = round((vecs.phi1partition[c1]-cons.phi1Range[0])/cons.phi1Range[2]);
				int m_min = round((vecs.phi1partition[c1-1]-cons.phi1Range[0])/cons.phi1Range[2]);
				int n_max = round((vecs.phi2partition[c2]-cons.phi2Range[0])/cons.phi2Range[2]);
				int n_min = round((vecs.phi2partition[c2-1]-cons.phi2Range[0])/cons.phi2Range[2]);
			
				for (int i_m = m_min; i_m <= m_max; i_m++) {
					for (int i_n = n_min; i_n <= n_max; i_n++) {						
						// Convert (i_m,i_n) to 1D index f, phi1-major order
					
						int numPhi2 = round((cons.phi2Range[1]-cons.phi2Range[0])/cons.phi2Range[2] + 1);
						int f = i_m*numPhi2 + i_n;
					
						if (vecs.integrationRule[i_section] == 1) {
							double intCoeff = vecs.simpsonCoeffsPhi1[i_m]*vecs.simpsonCoeffsPhi2[i_n];					
							sectionIntegrals[i_section] += intCoeff*exp(-vecs.dftData[f]/cons.kTBoltzIntegral);
						}
					}	
				}
				// Rescale sum according to integral formula
				if (vecs.integrationRule[i_section] == 1) {
					sectionIntegrals[i_section] *= (cons.phi1Range[2]*cons.phi2Range[2])/9.0;
				}
				// Add integral for this section to relevant conformer
				vecs.confIntegralsDFT[vecs.partitionMap[i_section]-1] += sectionIntegrals[i_section];
			
				i_section++;
			}		
		}
		
		cout << "\n DFT conformer integrals:\n";
		for (int i_conf = 0; i_conf < vecs.confIntegralsDFT.size(); i_conf++) {
			cout << vecs.confIntegralsDFT[i_conf] << endl;
		}
	} */
	
	return 0;
}

int constant_energy_process(constant_struct cons, vector_struct &vecs, fitting_param_struct &initialParams, int i_s)
{
	cout << "Processing constant energy terms in " << cons.numConfigs;
	cout << " geometries of surface " << i_s << ".\n\n";

	double aVec1[3], aVec2[3];
	double dVec1[3], dVec2[3], dVec3[3];
	double crossVec1[3], crossVec2[3], crossVec3[3];
	
	// Reinitialize type counters
	int totalUniqueDihedrals = initialParams.rbTypeCount.size();
	int totalUniquePairs = initialParams.ljSigma.size();
	int totalUniqueMols = initialParams.uShift.size();
		
	// Map molecule connectivity file to molecule ID

	// Search file list for this molecule
	bool molNew = 1;
	int molID;
	for (int mol=0; mol<vecs.sortedConFiles.size(); mol++) {
		if (vecs.connectFileList[i_s] == vecs.sortedConFiles[mol]) {
			molNew = 0;
			molID = mol;
		}
	}
	if (molNew) {
		cout << "Connect file " << vecs.connectFileList[i_s] << " is new molecule type. ID=" << totalUniqueMols << endl;
		vecs.sortedConFiles.push_back(vecs.connectFileList[i_s]);
		initialParams.surfIndexMapping[i_s] = totalUniqueMols++;
		initialParams.uShift.push_back(0.0);
	}
	else {
		initialParams.surfIndexMapping[i_s] = molID;
		cout << "Connect file " << vecs.connectFileList[i_s] << " is existing molecule type. ID=" << molID << endl;
	}
		
	// Read coordinates of all configs in i_s
	vector<double> coords(3*vecs.molSize[i_s].atoms*cons.numConfigs);
	
	for (int i_f=0; i_f < cons.numConfigs; i_f++) {
		for (int atom=0; atom<vecs.molSize[i_s].atoms; atom++) {
			int offSet = i_f*vecs.molSize[i_s].atoms*3 + atom*3;
			coords[offSet    ] = vecs.xyzData[i_s][offSet    ];
			coords[offSet + 1] = vecs.xyzData[i_s][offSet + 1];
			coords[offSet + 2] = vecs.xyzData[i_s][offSet + 2];
		}
	}
	// Bond energy currently ignored because TraPPE uses fixed bonds
	
	// Calculate angle (bending) energy
	for (int c=0; c < vecs.molSize[i_s].angles; c++) {
		
		int i = int(round(vecs.angleData[i_s][c*cons.angleDataSize    ] - 1)); 
		int j = int(round(vecs.angleData[i_s][c*cons.angleDataSize + 1] - 1));
		int k = int(round(vecs.angleData[i_s][c*cons.angleDataSize + 2] - 1));
		// -1 Because atoms are numbered starting from one in the input file
			
		double theta0 = vecs.angleData[i_s][c*cons.angleDataSize + 3]*cons.pi/180.0;
		double kTheta = vecs.angleData[i_s][c*cons.angleDataSize + 4];
	
		for (int i_f=0; i_f < cons.numConfigs; i_f++) {
			int offSet = 3*vecs.molSize[i_s].atoms*i_f;
			
			// vec j->i
			for (int xyz=0; xyz < 3; xyz++) {
				aVec1[xyz] = coords[offSet + i*3 + xyz] - coords[offSet + j*3 + xyz];
				aVec2[xyz] = coords[offSet + k*3 + xyz] - coords[offSet + j*3 + xyz];
			}
			double modVec1 = sqrt(aVec1[0]*aVec1[0] + aVec1[1]*aVec1[1] + aVec1[2]*aVec1[2]);
			double modVec2 = sqrt(aVec2[0]*aVec2[0] + aVec2[1]*aVec2[1] + aVec2[2]*aVec2[2]);

			double dot12 = aVec1[0]*aVec2[0] + aVec1[1]*aVec2[1] + aVec1[2]*aVec2[2];

			double theta = acos(dot12/(modVec1*modVec2));
			
			vecs.uConst[i_s][i_f].uAngle += 0.5*kTheta*(theta-theta0)*(theta-theta0);
			//cout << "theta, theta0: " << theta*180.0/cons.pi << ", " << theta0*180.0/cons.pi << endl;	
			//cout << "i_f, k_theta, U_a: " << i_f << ", " << kTheta << ", " << vecs.uConst[i_s][i_f].uAngle << endl;
		}
	}
		
	// Calculate dihedral angles and constant energy contribution
	for (int d=0; d < vecs.molSize[i_s].dihedrals; d++) {
		// Get dihedral type index
		int dType = int(round(vecs.dihedralData[i_s][d*cons.dihedralDataSize + 4]));
			
		// If this dihedral's parameters are being varied during fit
		if (dType != 0) {			
			// If new type (mapping array should be initialized to -1)
			if (initialParams.dihedralIndexMapping[dType] < 0) {	
				// Index dType maps to this dihedral number in fitting
				initialParams.dihedralIndexMapping[dType] = totalUniqueDihedrals++;	
				initialParams.rbTypeCount.push_back(1);	// First of this type			
				cout << "Type " << dType << " is new dihedral type " << totalUniqueDihedrals << endl;		
				// Add coeffs for new type
				for (int cosExp=1; cosExp < cons.nRBfit; cosExp++) {
					double newCoeff = vecs.dihedralData[i_s][d*cons.dihedralDataSize + 5 + cosExp];
					initialParams.rbCoeff.push_back(newCoeff);
				}	
			}
			else {
				// If not a new type, simply add one to count and average coeffs
				size_t count = ++initialParams.rbTypeCount[initialParams.dihedralIndexMapping[dType]];	
				
				for (int cosExp=1; cosExp < cons.nRBfit; cosExp++) {
					
					int offSet = initialParams.dihedralIndexMapping[dType]*(cons.nRBfit-1);
					double newCoeff = vecs.dihedralData[i_s][d*cons.dihedralDataSize + 5 + cosExp];
					
					// Take average of this and all the others of this type
					initialParams.rbCoeff[offSet + (cosExp-1)] = 
						(newCoeff/count) + (count-1)*initialParams.rbCoeff[offSet + (cosExp-1)]/count;
				}
			}
		}
		
		int i = int(round(vecs.dihedralData[i_s][d*cons.dihedralDataSize    ] - 1));
		int j = int(round(vecs.dihedralData[i_s][d*cons.dihedralDataSize + 1] - 1));
		int k = int(round(vecs.dihedralData[i_s][d*cons.dihedralDataSize + 2] - 1));
		int l = int(round(vecs.dihedralData[i_s][d*cons.dihedralDataSize + 3] - 1));
		//cout << "ijkl: " << i << " " << j << " " << k << " " << l << endl;

		for (int i_f=0; i_f < cons.numConfigs; i_f++) {
			//
			int offSet = 3*vecs.molSize[i_s].atoms*i_f;		
			// Calculate vectors b1 = rj-ri, b2 = rk-rj, b3 = rl-rk
			for (int xyz=0; xyz < 3; xyz++) {
				dVec1[xyz] = coords[offSet + j*3 + xyz] - coords[offSet + i*3 + xyz];
				dVec2[xyz] = coords[offSet + k*3 + xyz] - coords[offSet + j*3 + xyz];
				dVec3[xyz] = coords[offSet + l*3 + xyz] - coords[offSet + k*3 + xyz];
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
				
			for (int m=0; m < 3; m++) {
				crossVec1[m] /= modCrossVec1;
				crossVec2[m] /= modCrossVec2;
			}
				
			// Basis vector: d1 = cross(c1,b2)
			crossVec3[0] = crossVec1[1]*dVec2[2] - crossVec1[2]*dVec2[1];
			crossVec3[1] = crossVec1[2]*dVec2[0] - crossVec1[0]*dVec2[2];
			crossVec3[2] = crossVec1[0]*dVec2[1] - crossVec1[1]*dVec2[0];
			double modCrossVec3 = sqrt(crossVec3[0]*crossVec3[0] 
				+ crossVec3[1]*crossVec3[1] + crossVec3[2]*crossVec3[2]);
				
			for (int m=0; m < 3; m++) {
				crossVec3[m] /= modCrossVec3;
			}
			   
			double x = crossVec1[0]*crossVec2[0] + crossVec1[1]*crossVec2[1] +
				crossVec1[2]*crossVec2[2]; //dot(c1,c2);
			double y = -(crossVec3[0]*crossVec2[0] + crossVec3[1]*crossVec2[1] +
				crossVec3[2]*crossVec2[2]); //dot(d1,c2);
        
			double phi = atan2(y,x);
			double psi = phi - cons.pi;
			double cosPsi = cos(psi);
			//cout << "i_f, d, phi: " << i_f << ", " << d << ", " << phi*180.0/cons.pi << endl;
				
			// Write powers of cosPsi or store in constant energy 
			for (int cosExp=0; cosExp < cons.nRBfit; cosExp++) {
				
				double cosPsiPow = pow(cosPsi,cosExp);
				
				if (dType <= 0) {
					double constCoeff = vecs.dihedralData[i_s][d*cons.dihedralDataSize + 5 + cosExp];
					vecs.uConst[i_s][i_f].uDihedral += constCoeff*cosPsiPow;
				}
				else if (cosExp > 0) {
					vecs.cosPsiData[i_s].push_back(cosPsiPow);	
				}
			}	
		}	
	}
	// Calculate harmonic improper dihedral energy
	for (int id=0; id<vecs.molSize[i_s].impDihedrals; id++) {
		// Calculate angle between planes ijk and jkl
		int i = int(round(vecs.improperData[i_s][id*cons.improperDataSize    ] - 1)); // -1 for 0-based arrays
		int j = int(round(vecs.improperData[i_s][id*cons.improperDataSize + 1] - 1));
		int k = int(round(vecs.improperData[i_s][id*cons.improperDataSize + 2] - 1));
		int l = int(round(vecs.improperData[i_s][id*cons.improperDataSize + 3] - 1));
		//cout << "ijkl: " << i << " " << j << " " << k << " " << l << endl;

		for (int i_f=0; i_f < cons.numConfigs; i_f++) {
			
			int offSet = 3*vecs.molSize[i_s].atoms*i_f;
		
			// Calculate vectors b1 = rj-ri, b2 = rk-rj, b3 = rl-rk
			for (int xyz=0; xyz < 3; xyz++) {
				dVec1[xyz] = coords[offSet + j*3 + xyz] - coords[offSet + i*3 + xyz];
				dVec2[xyz] = coords[offSet + k*3 + xyz] - coords[offSet + j*3 + xyz];
				dVec3[xyz] = coords[offSet + l*3 + xyz] - coords[offSet + k*3 + xyz];
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
				
			for (int m=0; m < 3; m++) {
				crossVec1[m] /= modCrossVec1;
				crossVec2[m] /= modCrossVec2;
			}
				
			// Basis vector: d1 = cross(c1,b2)
			crossVec3[0] = crossVec1[1]*dVec2[2] - crossVec1[2]*dVec2[1];
			crossVec3[1] = crossVec1[2]*dVec2[0] - crossVec1[0]*dVec2[2];
			crossVec3[2] = crossVec1[0]*dVec2[1] - crossVec1[1]*dVec2[0];
			double modCrossVec3 = sqrt(crossVec3[0]*crossVec3[0] 
				+ crossVec3[1]*crossVec3[1] + crossVec3[2]*crossVec3[2]);
				
			for (int xyz=0; xyz < 3; xyz++) {
				crossVec3[xyz] /= modCrossVec3;
			}
			   
			double x = crossVec1[0]*crossVec2[0] + crossVec1[1]*crossVec2[1] + crossVec1[2]*crossVec2[2];
			double y = -(crossVec3[0]*crossVec2[0] + crossVec3[1]*crossVec2[1] + crossVec3[2]*crossVec2[2]);
        
			double phi = atan2(y,x);
			
			// Input params are in degrees and (kJ/mol)/rad^2
			double phi0 = (cons.pi/180.0)*vecs.improperData[i_s][id*cons.improperDataSize + 4];
			double kHarm = vecs.improperData[i_s][id*cons.improperDataSize + 5];
			
			// Shift phi into range phi0-180 <= phi <= phi0+180
			if (phi <= phi0-cons.pi) {
				phi += 2*cons.pi;
			}
			else if (phi >= phi0+cons.pi) {
				phi -= 2*cons.pi;
			}
			// Sum energy
			vecs.uConst[i_s][i_f].uImpDihedral += kHarm*(phi-phi0)*(phi-phi0);
		}
	}
		
	// Calculate interatomic distances, LJ and QQ non-bonded energy
	double rIJ[3];
		
	for (int i=0; i < vecs.molSize[i_s].atoms; i++) {
		for (int j=(i+1); j < vecs.molSize[i_s].atoms; j++) {
			
			int m1 = i*vecs.molSize[i_s].atoms+j;
			int atomI = int(round(vecs.atomData[i_s][i*cons.atomDataSize]))-1; // Get actual atom number from file data
			int atomJ = int(round(vecs.atomData[i_s][j*cons.atomDataSize]))-1; 		
			
			double ljSigma = vecs.sigmaMatrix[i_s][m1];
			double ljEps = vecs.epsilonMatrix[i_s][m1];
			
			// Loop over list of special pairs
			int fitPair = 0; 
			for (int pair=0; pair < vecs.molSize[i_s].pairs; pair++) {
				// If pair (i,j) is on the list of pairs to fit
				if ((vecs.pairData[i_s][pair*cons.pairDataSize] == (atomI+1)) 
				&& (vecs.pairData[i_s][pair*cons.pairDataSize+1] == (atomJ+1))) {
					// Get pair type index
					int pairType = int(round(vecs.pairData[i_s][pair*cons.pairDataSize + 2]));
			
					double newSigma = vecs.pairData[i_s][pair*cons.pairDataSize + 3];
					double newEpsilon = vecs.pairData[i_s][pair*cons.pairDataSize + 4];
					// If this pair's parameters are being varied during fit
					if (pairType != 0) {	
						
						fitPair = 1;								
						// If new type (mapping array should be initialized to -1)
						if (initialParams.pairIndexMapping[pairType] < 0) {	
							// Index pairType maps to this pair number in fitting
							initialParams.pairIndexMapping[pairType] = totalUniquePairs++;	
							initialParams.ljTypeCount.push_back(1.0);	// First of this type			
							
							// Add sigma and eps
							initialParams.ljSigma.push_back(newSigma);
							initialParams.ljEps.push_back(newEpsilon);
							cout << "Type " << pairType << " is new pair type " << totalUniquePairs << endl;		
						}
						else {
							// If not a new type, simply add one to count and average sigma/eps
							int mapping = initialParams.pairIndexMapping[pairType];
							size_t count = ++initialParams.ljTypeCount[initialParams.pairIndexMapping[pairType]];	
							// Take average of this and all the others of this type
							initialParams.ljSigma[mapping] = 
								(newSigma/count) + (count-1)*initialParams.ljSigma[mapping]/count;
							initialParams.ljEps[mapping] = 
								(newEpsilon/count) + (count-1)*initialParams.ljEps[mapping]/count;
						}
					}
					else {
						// Pair has special value of sigma and epsilon, but is not being fit
						ljSigma = newSigma;
						ljEps = newEpsilon;
						cout << "WARNING: pair (" << i << ", " << j << ") LJ parameters overwritten/n"; 
					}
				}
			}
			
			// If pair (i,j) wasn't on the list of pairs to fit, add to constant energy	
			for (int i_f=0; i_f < cons.numConfigs; i_f++) {
		
				int offSet = 3*vecs.molSize[i_s].atoms*i_f;
			
				// Vector pointing from i to j
				for (int xyz=0; xyz<3; xyz++) {
					rIJ[xyz] = coords[offSet + atomJ*3 + xyz] - coords[offSet + atomI*3 + xyz];
				}
			
				// Convert A to nm
				if (cons.xyzAngstroms) {
					for (int n=0; n<3; n++) {
						rIJ[n] *= 0.1;
					}
				}
			
				// Pair separation
				double r_ij = sqrt(rIJ[0]*rIJ[0] + rIJ[1]*rIJ[1] + rIJ[2]*rIJ[2]);
			
				vecs.rijMatrix[i_s][m1] = r_ij; // Just for printing
			
				// Coulomb energy,
				// 138.9355 converts q[e]q[e]/r[nm] to coulomb energy in kJ/mol
				double energyPairQQ = 138.93546*vecs.qqMatrix[i_s][m1]/r_ij;
				vecs.uConst[i_s][i_f].uQQ += energyPairQQ;
				//cout << i << ", " << j << ":\nr,qq    " << r_ij << ", " << vecs.qqMatrix[i_s][m1] << endl;
			
				if (fitPair == 1) {
					vecs.pairSepData[i_s].push_back(r_ij);					
				}
				else {
					double LJ6 = pow(ljSigma/r_ij,6);
					double energyPairLJ = 4.0*ljEps*(LJ6*LJ6-LJ6);
					//cout << "r,sigma " << r_ij << ", " << ljSigma << endl;
					
					// Split into 1-4, 1-5 and 1-6+ parts for plotting
					if (vecs.bondSepMat[i_s][m1] == 3) {
						vecs.uConst[i_s][i_f].uLJ14 += energyPairLJ;
						vecs.uConst[i_s][i_f].uQQ14 += energyPairQQ;
					} 
					else if (vecs.bondSepMat[i_s][m1] == 4) {
						vecs.uConst[i_s][i_f].u15 += (energyPairQQ + energyPairLJ);
					}
					else {
						vecs.uConst[i_s][i_f].u16plus += (energyPairQQ + energyPairLJ);
					}	
					// Add to total LJ 
					vecs.uConst[i_s][i_f].uLJ += energyPairLJ;
				}		
			}
		}
	} 

	cout << "r_ij matrix complete for final config of surface " << i_s << " :\n";
	for (int i = 0; i < vecs.molSize[i_s].atoms; i++) {
		for (int j = 0; j < vecs.molSize[i_s].atoms; j++) {
			cout << std::setw(10) << vecs.rijMatrix[i_s][i*vecs.molSize[i_s].atoms+j] << " ";
		}
		cout << endl;
	}
	
	ofstream outDist;
	if (cons.writeEnergyDist) {
		outDist.open("const_energy_dist.txt");
	}
	
	// Sum total constant energies
	for (int i_f=0; i_f < cons.numConfigs; i_f++) {
		vecs.uConst[i_s][i_f].sumEnergies();
		//cout << "Total constant energy (" << i_s << ", " << i_f << ") = " << vecs.uConst[i_s][i_f].uTotal << endl;
		if (cons.writeEnergyDist) {
			outDist << vecs.uConst[i_s][i_f].uAngle << " " << vecs.uConst[i_s][i_f].uImpDihedral << " ";
			outDist << vecs.uConst[i_s][i_f].uLJ14 << " " << vecs.uConst[i_s][i_f].uQQ14 << " ";
			outDist << vecs.uConst[i_s][i_f].uLJ - vecs.uConst[i_s][i_f].uLJ14 << " ";
			outDist << vecs.uConst[i_s][i_f].uQQ - vecs.uConst[i_s][i_f].uQQ14 << " " << endl;
		}
	}
	
	

	return 0;
}

double error_from_trial_point(constant_struct cons, vector_struct &vecs, fitting_param_struct initialParams, fitting_param_struct trialParams, bool toWrite)
{
	ofstream energyMD;
	ofstream energyDist;

	vector<vector<double> > energyTotal;
	vector<vector<double> > energyDihedral;
	vector<vector<double> > energyPair;
	
	energyTotal.resize(cons.numPhiSurfaces);
	energyDihedral.resize(cons.numPhiSurfaces);
	energyPair.resize(cons.numPhiSurfaces);
	
	if (toWrite == 1) {
		energyMD.open("energyMD.txt");
		energyDist.open("energyDist.txt");
	}
	//print_params_console(cons, trialParams);
	
	// Total sum of all contributions to objective function (F) per surface
	vector<double> sumF(cons.numPhiSurfaces,0.0);
	double totalSumF = 0.0;
	
	// Loop over energy surfaces
	for (int i_s=0; i_s<cons.numPhiSurfaces; i_s++) {
		//
		energyTotal[i_s].resize(cons.numConfigs,0.0);
		energyDihedral[i_s].resize(cons.numConfigs,0.0);
		energyPair[i_s].resize(cons.numConfigs,0.0);
		int i_psiPow = 0, i_pairSep = 0;
				
		// Loop over all dihedrals
		for (int d=0; d < vecs.molSize[i_s].dihedrals; d++) {
			// Get diherdal type
			int dType = int(round(vecs.dihedralData[i_s][d*cons.dihedralDataSize + 4]));
			// Map dihedral type to dihedral ID
			int rbID = trialParams.dihedralIndexMapping[dType];			
			//cout << "dType, dID " << dType << " " << rbID << endl;
			if (dType > 0) {
				//
				for (int i_f=0; i_f<cons.numConfigs; i_f++) {
					// This loop order is less efficient, but makes it easier to access cosPsiData
					for (int cosExp=1; cosExp<cons.nRBfit; cosExp++) {
						double coeff = trialParams.rbCoeff[rbID*(cons.nRBfit-1)+(cosExp-1)];
						double cosPsiPow = vecs.cosPsiData[i_s][i_psiPow++];
						//cout << "d, i_f, coeff, cosPsiPow " << d << ", " << i_f <<;
						//cout << ", " << coeff << ", " << cosPsiPow << endl; 	
						energyDihedral[i_s][i_f] += cosPsiPow*coeff;
					}
				}	
				for (int cosExp=1; cosExp<cons.nRBfit; cosExp++) {	
					double coeff = trialParams.rbCoeff[rbID*(cons.nRBfit-1)+(cosExp-1)];
					// Dihedral coeff restraint term added directly to sumF
					if (cons.resToZero == 1) {
						sumF[i_s] += cons.rbK*coeff*coeff;
					}
					else {
						double coeffDelta = coeff-initialParams.rbCoeff[rbID*(cons.nRBfit-1)+(cosExp-1)];
						sumF[i_s] += cons.rbK*coeffDelta*coeffDelta*cons.numConfigs;
					}
				}				
			}
		}
		// Calculate improper dihedral energy
			
		// Calculate pair energy (LJ only)
		for (int pair=0; pair < vecs.molSize[i_s].pairs; pair++) {		
			// Get pair type
			int pairType = int(round(vecs.pairData[i_s][pair*cons.pairDataSize + 2]));
				
			if (pairType > 0) {
				// Map pair type to pair ID
				int pairID = trialParams.pairIndexMapping[pairType];
				
				for (int i_f=0; i_f<cons.numConfigs; i_f++) {	
					double r_ij = vecs.pairSepData[i_s][i_pairSep++];
					double sigma = trialParams.ljSigma[pairID];
					double epsilon = trialParams.ljEps[pairID];
					double LJ6 = pow(sigma/r_ij, 6);
					
					energyPair[i_s][i_f] += 4.0*epsilon*(LJ6*LJ6-LJ6);
					//cout << "pair, sigma, epsilon, r_ij, energyPair: ";
					//cout << pair << ", " << sigma << ", " << epsilon <<;
					//cout << ", " << r_ij << ", " << energyPair[i_s][i_f] << endl;
				} 
				// LJ restraint term added directly to F
				double epsDelta = trialParams.ljEps[pairID] - initialParams.ljEps[pairID];
				sumF[i_s] += cons.epsK*epsDelta*epsDelta*double(cons.numConfigs)/double(trialParams.ljTypeCount[pairID]);
				//if (toWrite) cout << "epsDelta, epsK " << epsDelta << ", " << cons.epsK << endl;
			}	
		}	
		
		int molID = trialParams.surfIndexMapping[i_s];
		
		// Now total energy is known, sum F over configs
		for (int i_f=0; i_f<cons.numConfigs; i_f++) {	
			//
			energyTotal[i_s][i_f] = vecs.uConst[i_s][i_f].uTotal + trialParams.uShift[molID] 
				+ energyPair[i_s][i_f] + energyDihedral[i_s][i_f];
						
			// F += weighting*(delta*U)^2 
			vecs.energyDelta[i_s][i_f] = energyTotal[i_s][i_f]-vecs.dftData[i_s][i_f];
			sumF[i_s] += vecs.energyDelta[i_s][i_f]*vecs.energyDelta[i_s][i_f]*vecs.energyWeighting[i_s][i_f];
		
			// Write to file
			if (toWrite == 1) {
				//cout << "Total: " << energyTotal[i_s][i_f] << "\nDFT:   " << vecs.dftData[i_s][i_f] << endl;
				energyMD << energyTotal[i_s][i_f] << endl;
			}
		}
		// Normalize by number of configs on surface i_s
		sumF[i_s] /= double(cons.numConfigs);
		
		totalSumF += sumF[i_s];
	} // i_s loop
	
	/* Calculate Boltzmann distribution integrals -----------------------------------------
	if (cons.useBoltzIntRes == 1)
	{
		int numConformers = *std::max_element(vecs.partitionMap.begin(),vecs.partitionMap.end());
			
		vector<double> confIntegralsMD(numConformers, 0.0);
		compute_boltzmann_integrals(cons, vecs, energyTotal, confIntegralsMD, toWrite);
	
		// Calculate Boltzmann integral restraint terms
		assert(confIntegralsMD.size() == vecs.confIntegralsDFT.size());
		
		double sumIntegralsDFT = std::accumulate(vecs.confIntegralsDFT.begin(), vecs.confIntegralsDFT.end(), 0.0);
		double sumIntegralsMD = std::accumulate(confIntegralsMD.begin(), confIntegralsMD.end(), 0.0);
			
		double sumRMS = 0.0;
		
		if (toWrite == 1) {
			cout << "FracDFT    FracMD" << endl;
		}
		
		for (int i_conf = 0; i_conf < confIntegralsMD.size(); i_conf++)
		{
			double fracDFT = vecs.confIntegralsDFT[i_conf]/sumIntegralsDFT;
			double fracMD = confIntegralsMD[i_conf]/sumIntegralsMD;		
			
			sumResid += cons.kBoltzRes*(fracDFT-fracMD)*(fracDFT-fracMD);
			
			if (toWrite == 1) {
				sumRMS += (fracDFT-fracMD)*(fracDFT-fracMD);
				cout << std::setprecision(6) << fracDFT << " " << fracMD << " " << endl;			
			}
		}
		sumRMS = sqrt(sumRMS/confIntegralsMD.size());
		if (toWrite == 1) {
		cout << "RMS = " << sumRMS << endl;
		}
	} */
	
	return totalSumF;
}

int compute_gradient_F(constant_struct cons, vector_struct vecs, fitting_param_struct initialParams, 
fitting_param_struct currentParams, fitting_param_struct &grad)
{
	grad.zero_params();
	
	// Run error_from_trial_point to fill vecs.energyDelta (U_DFT - U_MD)
	error_from_trial_point(cons, vecs, initialParams, currentParams, 0);
	
	// Loop over energy surfaces
	for (int i_s=0; i_s<cons.numPhiSurfaces; i_s++) {
		//
		int i_psiPow = 0, i_pairSep = 0;
				
		// Loop over all dihedrals
		for (int d=0; d < vecs.molSize[i_s].dihedrals; d++) {
			// Get diherdal type
			int dType = int(round(vecs.dihedralData[i_s][d*cons.dihedralDataSize + 4]));
			// Map dihedral type to dihedral ID
			int rbID = currentParams.dihedralIndexMapping[dType];			
			//cout << "dType, dID " << dType << " " << rbID << endl;
			if (dType > 0) {
				//
				for (int i_f=0; i_f < cons.numConfigs; i_f++) {
					// Weighted delta 
					double wd = vecs.energyWeighting[i_s][i_f]*vecs.energyDelta[i_s][i_f];
					for (int cosExp=1; cosExp < cons.nRBfit; cosExp++) {					
						grad.rbCoeff[rbID*(cons.nRBfit-1)+(cosExp-1)] += 2*wd*vecs.cosPsiData[i_s][i_psiPow++];
					}
				}	
				// Dihedral coeff restraint term
				for (int cosExp=1; cosExp < cons.nRBfit; cosExp++) {
					//
					if (cons.resToZero == 1) {
						double rbDelta = currentParams.rbCoeff[rbID];
						grad.rbCoeff[rbID*(cons.nRBfit-1)+(cosExp-1)] += 
							2*cons.rbK*rbDelta*double(cons.numConfigs)/double(currentParams.rbTypeCount[rbID]);	
					}
					else {
						double rbDelta = currentParams.rbCoeff[rbID] - initialParams.rbCoeff[rbID];
						grad.rbCoeff[rbID*(cons.nRBfit-1)+(cosExp-1)] += 
							2*cons.rbK*rbDelta*double(cons.numConfigs)/double(currentParams.rbTypeCount[rbID]);	
					}
				}				
			}
		}
		// Calculate improper dihedral energy
			
		// Calculate pair energy (LJ only)
		for (int pair=0; pair < vecs.molSize[i_s].pairs; pair++) {		
			// Get pair type
			int pairType = int(round(vecs.pairData[i_s][pair*cons.pairDataSize + 2]));
				
			if (pairType > 0) {
				// Map pair type to pair ID
				int ljID = currentParams.pairIndexMapping[pairType];
				
				for (int i_f=0; i_f<cons.numConfigs; i_f++) {	
					double r_ij = vecs.pairSepData[i_s][i_pairSep++];
					double sigma = currentParams.ljSigma[ljID];
					double epsilon = currentParams.ljEps[ljID];
					double LJ6 = pow(sigma/r_ij, 6);
					
					// Weighted delta
					double wd = vecs.energyWeighting[i_s][i_f]*vecs.energyDelta[i_s][i_f];					
					grad.ljSigma[ljID] += 2*wd*4.0*epsilon*(LJ6*LJ6*12 - LJ6*6)/sigma;
					grad.ljEps[ljID] += 2*wd*4.0*(LJ6*LJ6-LJ6);
				} 
				// LJ restraint term added directly to F
				double epsDelta = currentParams.ljEps[ljID]-initialParams.ljEps[ljID];
				grad.ljEps[ljID] += 2*cons.epsK*epsDelta*double(cons.numConfigs)/double(currentParams.ljTypeCount[ljID]);				
			}	
		}	
		
		int molID = currentParams.surfIndexMapping[i_s];
		for (int i_f=0; i_f<cons.numConfigs; i_f++) {				
			grad.uShift[molID] += 2*vecs.energyWeighting[i_s][i_f]*vecs.energyDelta[i_s][i_f];
		}
	} // i_s loop
	
	// Apply scaling factor for sigma gradient (which has different units)
	for (int k=0; k < grad.ljSigma.size(); k++) {
		grad.ljSigma[k] *= cons.sigmaGradFactor;
	}
	
	//cout << "Current:\n";
	//print_params_console(cons, currentParams);
	//cout << "Grad:\n";
	//print_params_console(cons, grad);

	return 0;
} 

int param_linear_combine(fitting_param_struct &outParams, fitting_param_struct aParams, fitting_param_struct bParams,  
double a, double b)
{
	for (int i=0; i < outParams.uShift.size(); i++) {
		outParams.uShift[i] = a*aParams.uShift[i] + b*bParams.uShift[i];
	}
	for (int j=0; j < outParams.rbCoeff.size(); j++) {
		outParams.rbCoeff[j] = a*aParams.rbCoeff[j] + b*bParams.rbCoeff[j];
	}
	for (int k=0; k < outParams.ljEps.size(); k++) {
		outParams.ljSigma[k] = a*aParams.ljSigma[k] + b*bParams.ljSigma[k];
		outParams.ljEps[k] = a*aParams.ljEps[k] + b*bParams.ljEps[k];
	}
	return 0;
}

double param_scalar_product(fitting_param_struct aParams, fitting_param_struct bParams)
{
	double product = 0.0;
	
	for (int i=0; i < aParams.uShift.size(); i++) {
		product += aParams.uShift[i]*bParams.uShift[i];
	}
	for (int j=0; j < aParams.rbCoeff.size(); j++) {
		product += aParams.rbCoeff[j]*bParams.rbCoeff[j];
	}
	for (int k=0; k < aParams.ljEps.size(); k++) {
		product += aParams.ljSigma[k]*bParams.ljSigma[k];
		product += aParams.ljEps[k]*bParams.ljEps[k];
	}
	
	return product;
}

int error_sort(vector<double> simplexErrors, vector<int> &errorRankToRow)
{
	// Bubble sort (should be good for nearly sorted lists)
	int simplexSize = simplexErrors.size();
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

int compute_boltzmann_integrals(constant_struct cons, vector_struct vecs, vector<double> energyTotal, vector<double> &confIntegralsMD, bool toPrint)
{
	int numSections = (vecs.phi1partition.size()-1)*(vecs.phi2partition.size()-1);
	
	vector<double> sectionIntegrals(numSections, 0.0);
	
	// Calculate Boltzmann distribution integrals for each section
	int i_section = 0; 
	assert(vecs.integrationRule[i_section] == 1); // Until other integral methods supported
	
	for (int c1=1; c1<vecs.phi1partition.size(); c1++)
	{
		for (int c2=1; c2<vecs.phi2partition.size(); c2++)
		{				
			int m_max = round((vecs.phi1partition[c1]-cons.phi1Range[0])/cons.phi1Range[2]);
			int m_min = round((vecs.phi1partition[c1-1]-cons.phi1Range[0])/cons.phi1Range[2]);
			int n_max = round((vecs.phi2partition[c2]-cons.phi2Range[0])/cons.phi2Range[2]);
			int n_min = round((vecs.phi2partition[c2-1]-cons.phi2Range[0])/cons.phi2Range[2]);
			
			for (int i_m = m_min; i_m <= m_max; i_m++)
			{
				for (int i_n = n_min; i_n <= n_max; i_n++)
				{						
					// Convert (i_m,i_n) to 1D index f, phi1-major order					
					int numPhi2 = round((cons.phi2Range[1]-cons.phi2Range[0])/cons.phi2Range[2] + 1);
					int f = i_m*numPhi2 + i_n;
					
					if (vecs.integrationRule[i_section] == 1) // Simpsons rule
					{
						double intCoeff = vecs.simpsonCoeffsPhi1[i_m]*vecs.simpsonCoeffsPhi2[i_n];					
						sectionIntegrals[i_section] += intCoeff*exp(-energyTotal[f]/cons.kTBoltzIntegral);
					}
				}	
			}
			
			// Rescale sum according to integral formula
			if (vecs.integrationRule[i_section] == 1) // Simpsons rule
			{
				sectionIntegrals[i_section] *= (cons.phi1Range[2]*cons.phi2Range[2])/9.0;
			}
			// Map integral for this section to relevant conformer
			confIntegralsMD[vecs.partitionMap[i_section]-1] += sectionIntegrals[i_section];
			
			i_section++;
		}		
	}
	
	if (toPrint == 1)
	{
		cout << "Complete section and conformer integrals:\n\n";
		for (int i_sec = 0; i_sec < sectionIntegrals.size(); i_sec++) {
			cout << i_sec << ": " << sectionIntegrals[i_sec] << endl;
		}
		cout << endl;
		for (int i_conf = 0; i_conf < confIntegralsMD.size(); i_conf++) {
			cout << i_conf << ": " << confIntegralsMD[i_conf] << endl;
		}
	}
	return 0;
}

int print_params_console(constant_struct &cons, fitting_param_struct &printParams) 
{
	const int rbPerRow = 3;
	
	for (int i=0; i < printParams.uShift.size(); i++) {
		printf("%7.3f ", printParams.uShift[i]);
	}
	cout << endl;
	for (int j=0; j < printParams.rbCoeff.size(); j++) {
		printf("%7.3f ", printParams.rbCoeff[j]);
		if ((j+1)%((cons.nRBfit-1)*rbPerRow)==0) cout << endl; 
	}
	for (int k=0; k < printParams.ljEps.size(); k++) {
		printf("%7.4f ", printParams.ljSigma[k]);
		printf("%7.4f ", printParams.ljEps[k]);
	}
	cout << endl;
	return 0;
}

int print_simplex_console(constant_struct &cons, simplex_struct &printSim)
{
	for (int row=0; row<printSim.plex.size(); row++) {
		print_params_console(cons, printSim.plex[row]);	
	}	
	return 0;
}
