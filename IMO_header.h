// Header file for DME dihedral fitter

struct constant_struct
{
	long double pi;
	
	bool xyzAngstroms;
	bool genXyzFileNames;
	string inputFileString;	
	string xyzFile;
	string energyFile;
	string connectFile;
	string phiCoordFile;
	
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
	double phi1min;
	double phi1max;
	double phi1step;
	double phi2min;
	double phi2max;
	double phi2step;
	double phi3min;
	double phi3max;
	double phi3step;
  
	int annealWrite;
	int downhillWrite;
	int atomDataSize;
	int bondDataSize;
	int angleDataSize;
	int dihedralDataSize;
	int improperDataSize;
	int pairDataSize;
	int size[6];
	int numDihedralFit;
	int numDihedralParams;
	int numTotalParams;
	int simplexSize;
  
	bool resToZero;
	double epsK;
	double dihedralK;
	
	bool useBoltzIntRes;
	double kTBoltzIntegral;
	double kBoltzRes;
  
	double vTempInitial;
	double vTempFinal;
	double vTempFactor;
	double dihedralStep;
	double sigmaStep;
	double epsilonStep;
	double KT;
  
	double nmReflect;
	double nmExpand;
	double nmContract;
	double nmShrink;
  
	double sigmaGradFactor;
};

struct vector_struct
{
	
	vector<double> phi1partition;
	vector<double> phi2partition;
	
	vector<int> partitionMap;
	vector<int> integrationRule;
	vector<double> simpsonCoeffsPhi1;
	vector<double> simpsonCoeffsPhi2;
	vector<double> confIntegralsDFT;
	
	vector<double> xyzData;
	vector<double> energyData;
	vector<double> energyWeighting;
	vector<double> constantEnergy;
	
	vector<double> atomData;
	vector<double> bondData;
	vector<double> angleData;
	vector<double> dihedralData;
	vector<double> improperData;
	vector<double> ljPairFitData;
	
	vector<double> pairSepData;
	vector<double> phiPairSpec;
	vector<double> cosPsiData;
	vector<int> bondSepMat;
	vector<double> rijMatrix;
	vector<double> sigmaMatrix;
	vector<double> epsilonMatrix;
	vector<double> qqMatrix;
	
};

// Function prototypes
int read_input_params(constant_struct &cons, vector_struct &vecs);
int xyz_files_read(constant_struct cons, vector_struct &vecs);
int connectivity_read(constant_struct cons, vector_struct &vecs);
int connectivity_process(constant_struct cons, vector_struct &vecs);
int energy_read(constant_struct cons, vector_struct &vecs);
int constant_energy_process(constant_struct cons, vector_struct &vecs);
int define_initial_simplex(constant_struct cons, vector<double> initialParams, vector<double> &simplex);
long double error_from_trial_point(constant_struct cons, vector<double> initialParams, vector<double> &trialParams, int trialStart, vector_struct vecs, bool toWrite);
int compute_boltzmann_integrals(constant_struct cons, vector_struct vecs, vector<double> energyTotal, vector<double> &confIntegralsMD, bool toPrint);
int compute_gradient_F(constant_struct cons, vector_struct vecs, vector<double> initialParams, vector<double> currentParams, vector<double> &gradVector);
int error_sort(int simplexSize, vector<long double> simplexErrors, vector<int> &errorRankToRow);

int simulated_annealing(constant_struct cons, vector_struct vecs, vector<double> initialParams, vector<double> &currentParams, int annealIt);
int downhill_simplex(constant_struct cons, vector_struct vecs, vector<double> initialParams, vector<double> &currentParams, int simplexIt);
int steepest_descent(constant_struct cons, vector_struct vecs, vector<double> initialParams, vector<double> &currentParams, int gradientIt);
int conjugate_gradient(constant_struct cons, vector_struct vecs, vector<double> initialParams, vector<double> &currentParams, int gradientIt);

