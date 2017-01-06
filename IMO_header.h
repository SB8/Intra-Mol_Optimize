// Header file for DME dihedral fitter

#define xyzFormat "Angstroms" 
#define inputFileString "input_files_DME.txt"

const double PI_ = 3.14159265359;

struct constant_struct
{
	int nrexcl;
	int gromacsCombRule;
	double scale14LJ;
	double scale14QQ;
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
  
	double vTempInitial;
	double vTempFinal;
	double vTempFactor;
	double dihedralStep;
	double sigmaStep;
	double epsilonStep;
	double KT;
  
	double simplexAlpha;
	double simplexGamma;
	double simplexBeta;
	double simplexSigma;
  
	double sigmaGradFactor;
};

struct vector_struct
{
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
int read_input_params(constant_struct &cons, vector_struct &vecs, bool genFileNames);
int xyz_files_read(constant_struct &cons, vector_struct &vecs, string xyzFile, bool genFileNames, string phiCoordFile);
int connectivity_read(constant_struct &cons, vector_struct &vecs, string connectFile);
int connectivity_process(constant_struct &cons, vector_struct &vecs);
int energy_read(constant_struct &cons, vector_struct &vecs, string energyFile);
int constant_energy_process(constant_struct &cons, vector_struct &vecs);
int define_initial_simplex(constant_struct &cons, vector<double> initialParams, vector<double> &simplex);
double error_from_trial_point(constant_struct cons, vector<double> initialParams, vector<double> &trialParams, int trialStart, vector_struct vecs, bool toWrite);
int compute_gradient(constant_struct cons, vector_struct vecs, vector<double> initialParams, vector<double> currentParams, vector<double> &gradVector);
int error_sort(int simplexSize, vector<double> simplexErrors, vector<int> &errorRankToRow);

int simulated_annealing(constant_struct cons, vector_struct vecs, vector<double> initialParams, vector<double> &currentParams, int annealIt);
int downhill_simplex(constant_struct cons, vector_struct vecs, vector<double> initialParams, vector<double> &currentParams, int simplexIt);
int steepest_descent(constant_struct cons, vector_struct vecs, vector<double> initialParams, vector<double> &currentParams, int gradientIt);
int conjugate_gradient(constant_struct cons, vector_struct vecs, vector<double> initialParams, vector<double> &currentParams, int gradientIt);

