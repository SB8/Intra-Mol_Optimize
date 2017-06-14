
int define_initial_simplex(constant_struct cons, vector<double> initialParams, vector<double> &simplex)
{
	int ssM1 = cons.numTotalParams; // Width of simplex, i.e. total number of params
	const double factorLJ = 0.01;
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
		if (vecs.molSize[i_s].pairs != 0)
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

int downhill_simplex(constant_struct cons, vector_struct vecs, fitting_param_struct &initialParams, fitting_param_struct &currentParams, int simplexIt)
{
	cout << "\nPerforming " << simplexIt << " downhill simplex iteraitons\n\n";
	double nmEPS = 1.0E-14L;
	int printFreq = 100; // Print simplex to console every printFreq steps
	int numReflectSteps = 0;
	
	simplex_struct sStruc;
	
	// Generate initial simplex	
	vector<double> simplex(cons.simplexSize*(cons.numTotalParams));
	vector<double> simplexErrors(cons.simplexSize);
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
	for (int row=0; row<cons.simplexSize; row++) {
		simplexErrors[row] = error_from_trial_point(cons, vecs, initialParams, simplex[row], 0);
		cout << "Row: " << row << ", error " << simplexErrors[row] << endl;
	}
	
	for (int iD=0; iD < simplexIt; iD++) {
		//
		toWrite = (iD%cons.downhillWrite == 0);
		toPrint = (iD%printFreq == 0);
		
		if (toPrint) {
			cout << "\nIteration: " << iD << endl << endl;
		}

		// Sort errors, determine max, min and write best
		error_sort(cons.simplexSize, simplexErrors, errorRankToRow);
		
		int worstRow = errorRankToRow[cons.simplexSize-1];
		int maxPoint = worstRow*(cons.numTotalParams); // Start of row containing highest error point

		double maxError = simplexErrors[ errorRankToRow[cons.simplexSize-1] ];
		double secondMaxError = simplexErrors[ errorRankToRow[cons.simplexSize-2] ];
		double minError = simplexErrors[ errorRankToRow[0] ];
		
		// Return if end condition met
		if (std::fabs(maxError-minError)/maxError < nmEPS) {
			cout << std::setprecision(12) << simplexErrors[worstRow] << endl;
			cout << std::setprecision(12) << simplexErrors[0] << endl;
			cout << std::setprecision(12) << std::fabs(maxError-minError)/maxError << endl;
			cout << "\nEnding downhill simplex routine after " << iD << " iterations\n\n";
			cout << "The percentage of reflect steps was " << 100.0*double(numReflectSteps)/double(iD) << endl;
			
			if (simplexIt > 0) {
				// Update currentParams array
				cout << "\n Best params from simplex:\n";
				for (int col=0; col<(cons.numTotalParams); col++)
				{
					currentParams[col] = simplex[ errorRankToRow[0]*cons.numTotalParams + col ];		
					cout << currentParams[col] << " ";
				}
					
			}
			return 0;		
		}
		
		simplexStream << minError << endl;
		
		// Print best
		if (toPrint)
			cout << "MAX and MIN errors: " << std::setprecision(20) << maxError << " " << minError << endl;
		
		// Compute centroid of all points execpt worst
		for (int col=0; col<cons.numTotalParams; col++) {
			//
			double sumCentroidCol = 0.0;
			for (int row=0; row<cons.simplexSize; row++) {
				// If this isn't worst row
				if (row != errorRankToRow[cons.simplexSize-1]) {
					sumCentroidCol += simplex[row*(cons.numTotalParams) + col];
				}
			}
			centroid[col] = sumCentroidCol/(cons.simplexSize-1);
		}
		
		// Compute reflected point and output best params
		if (toPrint) {
			cout << "\nComputing reflected point in opposite direction to point " << errorRankToRow[cons.simplexSize-1] << endl;
		}
		
		for (int col=0; col<cons.numTotalParams; col++) {
			reflectParams[col] = centroid[col] + cons.nmReflect*(centroid[col] - simplex[maxPoint + col]);
		}
		
		double reflectError = error_from_trial_point(cons, vecs, initialParams, reflectParams, 0);
		
		if (toPrint) {
			cout << "\nTotal residual for reflected point is " << std::setprecision(14) << std::setw(16) << reflectError << endl << endl;
		}	
		
		// If neither new best nor worst (or 2nd worst), replace worst row of simplex with trial point
		if ((reflectError <= secondMaxError) && (reflectError >= minError)) {
			if (toPrint) {
				cout << "Reflected point neither best nor worst.\n";
				cout << "Replacing worst point with reflected point.\n";	
			}			
			for (int col=0; col<cons.numTotalParams; col++) {
				simplex[maxPoint + col] = reflectParams[col];
			}
			simplexErrors[worstRow] = reflectError;
			numReflectSteps++;
		}
		// If best, test with expansion in direction of reflected point
		else if (reflectError < minError) {
			if (toPrint) {
				cout << "Reflected point is better than any point in simplex.\n";
				cout << "Generating new point by expansion.\n\n";
			}
			
			for (int col=0; col<cons.numTotalParams; col++) {
				expansionParams[col] = reflectParams[col] + cons.nmExpand*(reflectParams[col] - centroid[col]);
			}
			// Test expanded point
			double expansionError = error_from_trial_point(cons, vecs, initialParams, expansionParams, 0);
			
			if (expansionError < minError) {
				if (toPrint) {
					cout << "Expanded point even better.\n";
					cout << "Replacing worst point with expanded point.\n\n";
				}			
				for (int col=0; col < (cons.numTotalParams); col++) {
					simplex[maxPoint + col] = expansionParams[col];
				}
				simplexErrors[worstRow] = expansionError;
			}
			else {
				if (toPrint) {
					cout << "Expansion unsuccessful.\n";
					cout << "Replacing worst point with reflected point.\n\n";
				}
				for (int col=0; col <cons.numTotalParams; col++) {
					simplex[maxPoint + col] = reflectParams[col];
				}
				simplexErrors[worstRow] = reflectError;
				numReflectSteps++;
			
			}
		}
		// Else, the reflected point error must be > the second worst point
		else {
			if (toPrint) {
				cout << "Reflection unsuccessful, performing contraction.\n";
			}
			//Replace worst point by reflected point if it's an improvement
			if (reflectError < maxError) {
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
			double contractError = error_from_trial_point(cons, vecs, initialParams, contractParams, 0);
			
			if (contractError < maxError) {
				if (toPrint) {
					cout << "Replacing worst point with contracted point.\n\n";
				}
				
				for (int col=0; col<cons.numTotalParams; col++)
				{
					simplex[maxPoint + col] = contractParams[col];
				}	
				simplexErrors[worstRow] = contractError;
			}
			else {
				// Reduction towards best point
				if (toPrint) {
					cout << "Contracted point is new worst point. Trying reduction towards best point.\n\n";
				}
				for (int row=0; row<cons.simplexSize; row++) {
					
					for (int col=0; col<cons.numTotalParams; col++)
					{
						double bestParam = simplex[errorRankToRow[0]*(cons.numTotalParams) + col];
						double currentParam = simplex[row*(cons.numTotalParams) + col];
						simplex[row*cons.numTotalParams + col] = bestParam + cons.nmShrink*(currentParam-bestParam);
					}	
					// Update error for reduced point
					simplexErrors[row] = error_from_trial_point(cons, vecs, initialParams, simplex[row], vecs, 0);
					
				}
			}	
		}
		
		// Output new simplex
		if (toPrint) {
			cout << "Updated simplex: " << endl;
			for (int row=0; row<cons.simplexSize; row++) {
				
			}
		}
	}
	
	if (simplexIt > 0) {
		// Update currentParams array
		
		for (int col=0; col<(cons.numTotalParams); col++) {
			currentParams[col] = simplex[errorRankToRow[0]*cons.numTotalParams + col];
		}				
	}
		
	cout << "\nEnding downhill simplex routine after max iterations\n\n";
	
	return 0;

}  