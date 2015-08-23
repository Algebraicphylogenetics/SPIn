#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

// Defining Variables ////////////////////////////////////////////////////////
typedef long long int LONG_INTEGER;
typedef int INTEGER;
//declaring the functions:
//void  average(std::vector<std::vector<int> > &counts, std::vector<double> &results);

// Program options ///////////////////////////////////////////////////////////
enum ProgramOptions {
	CREATE_FILE, LOAD_FILE, CALC_RULES_VALUES, SHOW_RESULTS, OPERATE, EXIT
};

// Global Variables //////////////////////////////////////////////////////////
#define MAX_CHAIN_LENGTH 100000
int MAX_NUM_SPECIES;
int SHIFT_CORRIMENT;
static int g_mask;
static int g_chainLength = -1;
static int g_numSpecies = -1;
static char* g_chainCharsSpeciesADN = 0; // Matrix dim: g_chainLength*g_numSpecies

static INTEGER g_valsSaved; // The chains that are not repeated...
static LONG_INTEGER *g_chainValsADN = 0; // Matrix dim: g_chainLength*NUM_ADN_OPTIONS
static INTEGER *g_numCountsData = 0; // Matrix dim: g_chainLength*NUM_ADN_OPTIONS
static std::string *g_nameSpecies = 0; // Vector dim: g_numSpecies
static double *g_average = 0; // Vector dim: g_chainLenght
double loglikSBD = 0.0;

// Info for the bases of the ADN /////////////////////////////////////////////
#define NUM_ADN_BASES 4
#define NUM_ADN_OPTIONS 4
enum ADNBases {
	ADENINA, CITOSINA, TIMINA, GUANINA, NONE
};
static std::string changeADNBases[1][6] = { { "(at)(cg)", "(ac)(gt)",
		"(ag)(ct)", "(ag)", "(ac)", "(at)" } };


//*********************************************************AIC/BIC for the SBD and ATR
double computeAICSBDATR( std::vector<LONG_INTEGER>  &counts, std::vector<double>  &mle,  double K, double n, double &bic) {

	double aic = 0.0;
	
	double loglik = 0.0;
	for (int i = 0; i < counts.size(); i++) {
	    if (mle[i]>0)
		loglik += counts[i] * log(mle[i]);
		}
	K = K - 1;
	aic = -2 * loglik + 2 * K + (2 * K * (K + 1)) / (n - K - 1);	
	bic = -2 * loglik + log(n) * K ;	

	return aic;
	
	}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERAL function

void get_average(std::vector<std::vector<int> > &counts,
		std::vector<double> &results, std::vector<double> &orbitSum) {
	if (results.size() != 0) {
		results.clear();
	}
	if (orbitSum.size() != 0) {
		orbitSum.clear();
	}

	for (int i = 0; i < counts.size(); i++) {
		double accum = 0.0;
		for (int j = 0; j < counts[i].size(); j++) {
			accum += counts[i][j];
		}
		orbitSum.push_back(accum);
		accum /= counts[i].size();
		results.push_back(accum);
	}
}

double computeAIC(std::vector<std::vector<int> > &counts, double K, double n, double &bic) {
	double aic = 0.0;
	std::vector<double> averageVec; // vector of averages
	std::vector<double> cumSumVec; // vector of cumm sums

	get_average(counts, averageVec, cumSumVec);

	double loglik = 0.0;

	for (int i = 0; i < counts.size(); i++) {
		loglik += cumSumVec[i] * log((double)averageVec[i]/(double)g_chainLength);
	}
	K = K - 1;

	aic = -2 * loglik + 2 * K + (2 * K * (K + 1)) / (n - K - 1);	
	bic = -2 * loglik +  log(n) * K ;	

	return aic;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Return the max number of species the system can work with
int MaxNumSpecies() {
	int numBits = sizeof(LONG_INTEGER) * 8;
	int val = -1;

	if (NUM_ADN_BASES + 1 <= 0) {
		std::cout << "Num of ADN Bases cannot be negative or zero."
				<< std::endl;
	} else if (NUM_ADN_BASES + 1 <= 2) {
		val = numBits;
	} else if (NUM_ADN_BASES + 1 <= 4) {
		val = numBits / 2;
	} else if (NUM_ADN_BASES + 1 <= 8) {
		val = numBits / 3;
	} else if (NUM_ADN_BASES + 1 <= 16) {
		val = numBits / 4;
	} else {
		std::cout
				<< "Num of ADN Bases superior to 16. Contact the program administrator."
				<< std::endl;
	}

	return val;
}

// Get SHIFT corriment info for bit manipulation
int GetShiftCorriment() {
	int numBits = sizeof(LONG_INTEGER) * 8;
	int val = -1;

	if (NUM_ADN_BASES + 1 <= 0) {
		std::cout << "Num of ADN Bases cannot be negative or zero."
				<< std::endl;
	} else if (NUM_ADN_BASES + 1 <= 2) {
		val = 1;
	} else if (NUM_ADN_BASES + 1 <= 4) {
		val = 2;
	} else if (NUM_ADN_BASES + 1 <= 8) {
		val = 3;
	} else if (NUM_ADN_BASES + 1 <= 16) {
		val = 4;
	} else {
		std::cout
				<< "Num of ADN Bases superior to 16. Contact the program administrator."
				<< std::endl;
	}

	return val;
}

void SetGlobalMask() {
	if (NUM_ADN_BASES + 1 <= 0) {
		std::cout << "Num of ADN Bases cannot be negative or zero."
				<< std::endl;
	} else if (NUM_ADN_BASES + 1 <= 2) {
		g_mask = 1;
	} else if (NUM_ADN_BASES + 1 <= 4) {
		g_mask = 3;
	} else if (NUM_ADN_BASES + 1 <= 8) {
		g_mask = 7;
	} else if (NUM_ADN_BASES + 1 <= 16) {
		g_mask = 15;
	} else {
		std::cout
				<< "Num of ADN Bases superior to 16. Contact the program administrator."
				<< std::endl;
	}
}

inline char ADNBaseValToChar(int adnVal) {
	if (adnVal == 0)
		return 'a';
	else if (adnVal == 1)
		return 'c';
	else if (adnVal == 2)
		return 'g';
	else if (adnVal == 3)
		return 't';

	// Error
	std::cout << "Error converting ADNBase in ADNBaseValToChar." << std::endl;
	return 0;
}

inline int ADNBaseCharToVal(int adnChar) {
	if (adnChar == 'a')
		return 0;
	else if (adnChar == 'c')
		return 1;
	else if (adnChar == 'g')
		return 2;
	else if (adnChar == 't')
		return 3;
	// Error
	std::cout << "Error converting ADNBase in ADNBaseCharToVal." << std::endl;
	return -1;
}

LONG_INTEGER change_adn_bases_simple(LONG_INTEGER value, int col, int fil) {
	std::string change = changeADNBases[fil][col];
	LONG_INTEGER val1 = ADNBaseCharToVal(change[1]);
	LONG_INTEGER val2 = ADNBaseCharToVal(change[2]);

	LONG_INTEGER new_value = 0;
	for (int i = 0; i < g_numSpecies; i++) {
		LONG_INTEGER val = value >> i * SHIFT_CORRIMENT & g_mask;
		if (val == val1) {
			new_value |= (val2 << i * SHIFT_CORRIMENT);
		} else if (val == val2) {
			new_value |= (val1 << i * SHIFT_CORRIMENT);
		} else {
			new_value |= (val << i * SHIFT_CORRIMENT);
		}
	}

	return new_value;
}

LONG_INTEGER change_adn_bases(LONG_INTEGER value, int col, int fil) {
	std::string change = changeADNBases[fil][col];
	LONG_INTEGER val1 = ADNBaseCharToVal(change[1]);
	LONG_INTEGER val2 = ADNBaseCharToVal(change[2]);
	LONG_INTEGER val3 = ADNBaseCharToVal(change[5]);
	LONG_INTEGER val4 = ADNBaseCharToVal(change[6]);

	LONG_INTEGER new_value = 0;
	for (int i = 0; i < g_numSpecies; i++) {
		LONG_INTEGER val = value >> i * SHIFT_CORRIMENT & g_mask;
		if (val == val1) {
			new_value |= (val2 << i * SHIFT_CORRIMENT);
		} else if (val == val2) {
			new_value |= (val1 << i * SHIFT_CORRIMENT);
		} else if (val == val3) {
			new_value |= (val4 << i * SHIFT_CORRIMENT);
		} else if (val == val4) {
			new_value |= (val3 << i * SHIFT_CORRIMENT);
		}
	}

	return new_value;
}

std::string transform_adn_chain_val_to_string(LONG_INTEGER val) {
	std::string adnChain;
	for (int i = g_numSpecies - 1; i >= 0; i--) {
		LONG_INTEGER value = val >> i * SHIFT_CORRIMENT & g_mask;
		adnChain.push_back(ADNBaseValToChar(value));
	}

	return adnChain;
}
// Free System Memory
void FreeMemory() {
	if (g_average) {
		delete[] g_average;
		g_average = 0;
	}
	if (g_chainCharsSpeciesADN) {
		delete[] g_chainCharsSpeciesADN;
		g_chainCharsSpeciesADN = 0;
	}
	if (g_chainValsADN) {
		delete[] g_chainValsADN;
		g_chainValsADN = 0;
	}
	if (g_numCountsData) {
		delete[] g_numCountsData;
		g_numCountsData = 0;
	}
	if (g_nameSpecies) {
		delete[] g_nameSpecies;
		g_nameSpecies = 0;
	}
}

// Reserve Memory
void ReserveMemory(int numSpecies, int chainLength) {
	g_chainCharsSpeciesADN = new char[numSpecies * chainLength];
	g_chainValsADN = new LONG_INTEGER[chainLength * NUM_ADN_OPTIONS];
	g_numCountsData = new INTEGER[chainLength * NUM_ADN_OPTIONS];
	g_nameSpecies = new std::string[numSpecies];
}

// Crate a random file
void create_file(std::string &fileName, int numSpecies, int chainLength) {
	std::ofstream myfile;

	// Assert the demanding is possible
	if (chainLength > MAX_CHAIN_LENGTH) {
		std::cout << "ADN chain too long." << std::endl;
		return;
	}
	if (numSpecies > MAX_NUM_SPECIES) {
		std::cout << "Number of species cannot be superior to "
				<< MAX_NUM_SPECIES << "." << std::endl;
		return;
	}

	// Assert Name is correct
	int pos;
	if ((pos = fileName.find(".fa")) != -1)
		fileName.replace(pos, 3, ".fa");
	if (fileName.find(".fa") == -1)
		fileName.append(".fa");

	// Save random Information
	myfile.open(fileName.c_str());
	myfile << numSpecies << " " << chainLength << std::endl;
	for (int i = 0; i < numSpecies; i++) {
		myfile << "Taxon" << i << std::endl;
		for (int j = 0; j < chainLength; j++) {
			myfile << ADNBaseValToChar(rand() % NUM_ADN_BASES) << std::endl;
		}
	}
	myfile.close();
}

// Load a random file
void load_file(std::string &fileName) {
	std::ifstream myfile;

	// Assert the name is correct
	int pos;
	if ((pos = fileName.find(".fa")) != -1)
		fileName.replace(pos, 3, ".fa");
	if (fileName.find(".fa") == -1)
		fileName.append(".fa");

	myfile.open(fileName.c_str());
	int numSpecies;
	int chainLength;
	 
	//if (data.size() > 0) {
	myfile >> numSpecies;
	myfile >> chainLength;
	std::cout << "For species " << numSpecies << " len " << chainLength << " myfile " << myfile << std::endl;
 	if (numSpecies > 0 && chainLength > 0) {
		// Free Memory
		FreeMemory();

		// Reserve Memory
		g_numSpecies = numSpecies;
		g_chainLength = chainLength;
		ReserveMemory(numSpecies, chainLength);
		// Load Info
		std::string data;
		bool read = true;
		INTEGER specie = -1;
		INTEGER valSavedForSpecie = 0;
		while (read) {
			myfile >> data;
			if (data.size() > 0) {
				if (data[0] == '>') {
					specie++;
					valSavedForSpecie = 0;
					g_nameSpecies[specie] = data;
				} else if (myfile.eof()) {
					read = false;
				} else {
					for (unsigned int i = 0; i < data.size(); i++) {
						if (valSavedForSpecie < g_chainLength) {
							g_chainCharsSpeciesADN[specie * g_chainLength
									+ valSavedForSpecie] = tolower(data[i]);
							valSavedForSpecie++;
							if (valSavedForSpecie == 1000)
								int s = 0;
						} else {
							std::cout << "For species "
									<< g_nameSpecies[specie]
									<< " the bases truncated for " << data[i]
									<< "." << std::endl;
						}
					}
				}
			}
		}
		
		/*for(int i = 0; i < g_numSpecies; i++)
		 {
		 myfile >> g_nameSpecies[i];
		 for(int j = 0; j < g_chainLength; j++)
		 {
		 myfile >> g_chainCharsSpeciesADN[i*g_chainLength + j];
		 }
		 }*/
		// Set Values to Zero
		memset(g_numCountsData, 0, NUM_ADN_OPTIONS * g_chainLength
				* sizeof(INTEGER));
		memset(g_chainValsADN, 0, NUM_ADN_OPTIONS * g_chainLength
				* sizeof(LONG_INTEGER));
		g_valsSaved = 0;
	}

	myfile.close();
}

// Save Information
void SaveValue(LONG_INTEGER value, INTEGER &valsSaved) {
	bool valFound = false;
	int valPosCol = -1;
	int valPosFil = -1;
	for (int i = 0; i < valsSaved && !valFound; i++) {
		for (int j = 0; j < NUM_ADN_OPTIONS && !valFound; j++) {
			if (g_chainValsADN[j * g_chainLength + i] == value) {
				valPosFil = j;
				valPosCol = i;
				valFound = true;
			}
		}
	}
	if (valFound) {
		g_numCountsData[valPosFil * g_chainLength + valPosCol] += 1;
#ifdef SHOW_ALG_INFO
		std::cout << "-Val Found----------------------------" << std::endl;
		std::cout << transform_adn_chain_val_to_string(g_chainValsADN[valPosFil*g_chainLength + valPosCol]) << std::endl;
#endif
	} else {
#ifdef SHOW_ALG_INFO
		std::cout << "---------------------------------------" << std::endl;
#endif
		g_numCountsData[valsSaved] += 1;
		g_chainValsADN[valsSaved] = value;
#ifdef SHOW_ALG_INFO
		std::cout << transform_adn_chain_val_to_string(value) << std::endl;
#endif
		for (int j = 1; j < NUM_ADN_OPTIONS; j++) {
			g_chainValsADN[j * g_chainLength + valsSaved] = change_adn_bases(
					value, j - 1, 0);
#ifdef SHOW_ALG_INFO
			std::cout << transform_adn_chain_val_to_string(change_adn_bases(value, j-1, 0)) << std::endl;
#endif
		}
		valsSaved++;
	} 
}

// Calculate the system probabilities
void calc_rules_and_values() {
	if (!g_chainCharsSpeciesADN) {
		std::cout << "Error: No file Loaded." << std::endl;
		return;
	}

	INTEGER valsSaved = 0;
	for (int i = 0; i < g_chainLength; i++) {
		LONG_INTEGER valChunk = 0;
		for (int j = 0; j < g_numSpecies; j++) {
			int val = ADNBaseCharToVal(g_chainCharsSpeciesADN[j * g_chainLength
					+ i]);
			valChunk = (valChunk << SHIFT_CORRIMENT) | val;
		}
		// Save the val
		SaveValue(valChunk, valsSaved);
	}
	g_valsSaved = valsSaved;
}

void show_results() {
	if (!g_chainCharsSpeciesADN) {
		std::cout << "Error: No file Loaded." << std::endl;
		return;
	}
	if (g_valsSaved == 0) {
		std::cout << "No computation done." << std::endl;
	}

	int option = -1;
	std::cout << " - (0) Show To Screen." << std::endl;
	std::cout << " - (1) Save To File." << std::endl;
	std::cin >> option;

	if (option == 0) {
		std::cout << "Total Options Found: " << g_valsSaved << std::endl;
		LONG_INTEGER expected = 1;
		for (int i = 1; i < g_numSpecies; i++)
			expected *= 4;
		std::cout << "Total Expected: " << expected << std::endl;
		for (int i = 0; i < g_valsSaved; i++) {
			std::cout
					<< "------------------------------------------------------------------"
					<< std::endl;
			for (int j = 0; j < NUM_ADN_OPTIONS; j++) {
				std::cout << "P(" << transform_adn_chain_val_to_string(
						g_chainValsADN[j * g_chainLength + i]) << ") = "
						<< g_numCountsData[j * g_chainLength + i] << std::endl;
			}
		}
	} else if (option == 1) {
		std::string fileName;
		std::cout << "Name of the file to save: ";
		std::cin >> fileName;

		int pos;
		if ((pos = fileName.find(".txt")) != -1)
			fileName.replace(pos, 4, ".txt");
		if (fileName.find(".txt") == -1)
			fileName.append(".txt");

		// Save
		std::ofstream myfile;
		myfile.open(fileName.c_str());
		myfile << "Total Options Found: " << g_valsSaved << std::endl;
		INTEGER expected = 1;
		for (int i = 1; i < g_numSpecies; i++)
			expected *= 4;
		myfile << "Total Expected: " << expected << std::endl;
		for (int i = 0; i < g_valsSaved; i++) {
			myfile
					<< "------------------------------------------------------------------"
					<< std::endl;
			for (int j = 0; j < NUM_ADN_OPTIONS; j++) {
				myfile << "P(" << transform_adn_chain_val_to_string(
						g_chainValsADN[j * g_chainLength + i]) << ") = "
						<< g_numCountsData[j * g_chainLength + i] << std::endl;
			}
		}
		myfile.close();
	}
}

void average() {
	if (g_average != 0) {
		delete[] g_average;
		g_average = 0;
	}
	g_average = new double[g_valsSaved];

	for (int i = 0; i < g_valsSaved; i++) {
		double accum = 0.0;
		for (int j = 0; j < NUM_ADN_OPTIONS; j++) {
			accum += g_numCountsData[j * g_chainLength + i];
		}
		accum /= NUM_ADN_OPTIONS;
		g_average[i] = accum;
	}

	for (int i = 0; i < g_valsSaved; i++) {
		std::cout << "Average Orbital " << i << ": " << g_average[i]
				<< std::endl;
	}
}

//Models
class ModelSMM;
class ModelK81;
class ModelK80;
class ModelJC;
class ModelSBD;
class ModelATR;

class ModelGMM {
public:

	ModelGMM() :
		_calculated(false) {
		if (!g_chainCharsSpeciesADN) {
			std::cout << "No File Loaded.\n" << std::endl;
			return;
		}
		_numSpecies = g_numSpecies;
		_adnChainLength = g_chainLength;
		CalcRulesAndValues();
	}

	void ShowResults() {
		if (!g_chainCharsSpeciesADN) {
			std::cout << "Error: No file Loaded." << std::endl;
			return;
		}
		if (!_calculated) {
			std::cout << "No computation done." << std::endl;
		}

		int option = 1;
	//	std::cout << " - (0) Show To Screen." << std::endl;
	//	std::cout << " - (1) Save To File." << std::endl;
		//std::cin >> option;

		if (option == 0) {
			std::cout << "Observed number of orbits under the GMM: " << _orbitals.size()
					<< std::endl;
			LONG_INTEGER expected = 1;
			for (int i = 0; i < _numSpecies; i++)
				expected *= 4;
			std::cout << "Expected number of orbits under the GMM: " << expected << std::endl;
			for (int i = 0; i < _orbitals.size(); i++) {
				std::cout
						<< "------------------------------------------------------------------"
						<< std::endl;
				std::cout << "P(" << transform_adn_chain_val_to_string(
						_orbitals[i]) << ") = " << _counts[i] << std::endl;

			}
		} else if (option == 1) {
			std::string fileName;
			std::cout << "Name of the file to save: ";
			std::cin >> fileName;

			int pos;
			if ((pos = fileName.find(".dat")) != -1)
				fileName.replace(pos, 4, ".txt");
			if (fileName.find(".txt") == -1)
				fileName.append(".txt");

			// Save
			std::ofstream myfile;
			myfile.open(fileName.c_str());
			myfile << "Observed number of orbits under the GMM: " << _orbitals.size() << std::endl;


			INTEGER expected = pow(4, _numSpecies);

			myfile << "Expected number of orbits under the GMM: " << expected << std::endl;
			for (int i = 0; i < _orbitals.size(); i++) {
				myfile
						<< "------------------------------------------------------------------"
						<< std::endl;
				myfile << "P(" << transform_adn_chain_val_to_string(
						_orbitals[i]) << ") = " << _counts[i] << std::endl;
			}
			myfile.close();
		}
	}
	bool IsCalculated() {
		return _calculated;
	}

	INTEGER _numSpecies;
	INTEGER _adnChainLength;

	std::vector<LONG_INTEGER> _orbitals;
	std::vector<INTEGER> _counts;

	bool _calculated;

	// Private methods
	void CalcRulesAndValues() {
		for (int i = 0; i < g_chainLength; i++) {
			LONG_INTEGER valOrbital = 0;
			for (int j = 0; j < g_numSpecies; j++) {
				int val = ADNBaseCharToVal(g_chainCharsSpeciesADN[j
						* g_chainLength + i]);
				valOrbital = (valOrbital << SHIFT_CORRIMENT) | val;
			}
			// Save the val
			SaveOrbital(valOrbital);
		}
		_calculated = true;
	}

	void SaveOrbital(LONG_INTEGER &newOrbital) {
		bool orbFound = false;
		int orbPosCol = -1;
		for (int i = 0; i < _orbitals.size() && !orbFound; i++) {
			if (_orbitals[i] == newOrbital) {
				orbPosCol = i;
				orbFound = true;
			}
		}
		if (orbFound) {
			_counts[orbPosCol] += 1;
		} else {
			_counts.push_back(1);
			_orbitals.push_back(newOrbital);
		}
	}
};

static ModelGMM *g_modelGMM = 0;

class ModelSMM {
public:

	ModelSMM(ModelGMM *pmodelGMM) :
		_calculated(false) {
		if (!pmodelGMM || !pmodelGMM->IsCalculated()) {
			std::cout << "Error: Model GMM not created or calculated.\n"
					<< std::endl;
			return;
		}
		_numSpecies = g_numSpecies;
		_adnChainLength = g_chainLength;
		CalcRulesAndValues(pmodelGMM);
	}

	double GetTheoreticalOrbits() {
		double expected = pow(4, _numSpecies) / 2.0;
		return expected;
	}

	void ShowResults() {
		if (!g_chainCharsSpeciesADN) {
			std::cout << "Error: No file Loaded." << std::endl;
			return;
		}
		if (!_calculated) {
			std::cout << "No computation done." << std::endl;
		}

		int option = 1;
//		std::cout << " - (0) Show To Screen." << std::endl;
	//	std::cout << " - (1) Save To File." << std::endl;
		//std::cin >> option;

		if (option == 0) {
			std::cout << "Observed number of orbits under the SMM: " << _orbitals.size()
					<< std::endl;
			LONG_INTEGER expected = 1;
		
			expected = GetTheoreticalOrbits();
			std::cout << "Expected number of orbits under the SMM: " << expected << std::endl;
			for (int i = 0; i < _orbitals.size(); i++) {
				std::cout
						<< "------------------------------------------------------------------"
						<< std::endl;
				for (int j = 0; j < _orbitals[i].size(); j++)
					std::cout << "P(" << transform_adn_chain_val_to_string(
							_orbitals[i][j]) << ") = " << _counts[i][j]
							<< std::endl;

			}
		} else if (option == 1) {
			std::string fileName;
			std::cout << "Name of the file to save: ";
			std::cin >> fileName;

			int pos;
			if ((pos = fileName.find(".dat")) != -1)
				fileName.replace(pos, 4, ".txt");
			if (fileName.find(".txt") == -1)
				fileName.append(".txt");

			// Save
			std::ofstream myfile;
			myfile.open(fileName.c_str());
			myfile << "Observed number of orbits under the SMM: " << _orbitals.size() << std::endl;
			INTEGER expected = 1;
			for (int i = 0; i < _numSpecies; i++)
				expected *= 4;
			expected /= 2.0;
			myfile << "Expected number of orbits under the SMM: " << expected << std::endl;
			for (int i = 0; i < _orbitals.size(); i++) {
				myfile
						<< "------------------------------------------------------------------"
						<< std::endl;
				for (int j = 0; j < _orbitals[i].size(); j++)
					myfile << "P(" << transform_adn_chain_val_to_string(
							_orbitals[i][j]) << ") = " << _counts[i][j]
							<< std::endl;
			}
			myfile.close();
		}
	}

	bool IsCalculated() {
		return _calculated;
	}

	INTEGER _numSpecies;
	INTEGER _adnChainLength;

	std::vector<std::vector<LONG_INTEGER> > _orbitals;
	std::vector<std::vector<INTEGER> > _counts;

	bool _calculated;

	// Private methods
	void CalcRulesAndValues(ModelGMM *pmodelGMM) {
		LONG_INTEGER valOrbital;
		INTEGER *foundBefore = new INTEGER[pmodelGMM->_orbitals.size()];
		memset(foundBefore, 0, sizeof(INTEGER) * pmodelGMM->_orbitals.size());
		for (int i = 0; i < pmodelGMM->_orbitals.size(); i++) {
			if (!foundBefore[i]) {
				int pos = -1;
				SearchTwinOrbital(pmodelGMM, i, pos);
				if (pos != -1)
					foundBefore[pos] = 1;
				SaveOrbital(pmodelGMM, i, pos);
			}
		}
		delete foundBefore;
		_calculated = true;
	}

	void SearchTwinOrbital(ModelGMM *pmodelGMM, int posMain, int &posTwin) {
		posTwin = -1;
		LONG_INTEGER mainOrbVal = pmodelGMM->_orbitals[posMain];
		LONG_INTEGER twinOrbVal = change_adn_bases(mainOrbVal, 0, 0);

		bool bFound = false;
		for (int i = posMain + 1; i < pmodelGMM->_orbitals.size() && !bFound; i++) {
			if (pmodelGMM->_orbitals[i] == twinOrbVal) {
				posTwin = i;
				bFound = true;
			}
		}
	}

	void SaveOrbital(ModelGMM *pmodelGMM, int posMain, int posTwin) {
		LONG_INTEGER mainOrbVal = pmodelGMM->_orbitals[posMain];
		LONG_INTEGER twinOrbVal = change_adn_bases(mainOrbVal, 0, 0);

		std::vector<LONG_INTEGER> newOrbital;
		newOrbital.push_back(mainOrbVal);
		newOrbital.push_back(twinOrbVal);
		_orbitals.push_back(newOrbital);

		std::vector<INTEGER> newCounts;
		newCounts.push_back(pmodelGMM->_counts[posMain]);
		if (posTwin != -1)
			newCounts.push_back(pmodelGMM->_counts[posTwin]);
		else
			newCounts.push_back(0);
		_counts.push_back(newCounts);
	}
};

static ModelSMM *g_modelSMM = 0;

//*****************************************************

class ModelSBD {
public:

	ModelSBD(ModelGMM *pmodelGMM) :
		_calculated(false) {
		if (!pmodelGMM || !pmodelGMM->IsCalculated()) {
			std::cout << "Error: Model GMM not created or calculated.\n"
					<< std::endl;
			return;
		}
		memset(_eqSBD, 0, 3*sizeof(INTEGER*));
		_numSpecies = g_numSpecies;
		_adnChainLength = g_chainLength;

		CalcRulesAndValues(pmodelGMM);

	}

    INTEGER *_eqSBD[3]; // for all 6 combinations of two leaves, A has 4 choices, C has 2 two choices, G has 1 choice  
	INTEGER _numClasses;
	
	INTEGER _numSpecies;
	INTEGER _adnChainLength;
	double mleSBD;


	std::vector<LONG_INTEGER> _orbitals;
	std::vector<double> _counts;
    std::vector<LONG_INTEGER> _patterns;
	
	void CreatePairs()
	{
		FreePairs();
		_numClasses = g_numSpecies-1;
		for(int i = 0; i < 3; i++)
    	{	
         		_eqSBD[i] = new INTEGER[_numClasses];  

         		for(int k = 0; k < _numClasses; k++)
         			_eqSBD[i][k] = 0; 
    	}	
	}
	
	void FreePairs()
	{
		for(int i = 0; i < 3; i++)
    	{	       		
       			if(_eqSBD[i])
       			{
       				delete []_eqSBD[i];
       			}   
    	}
	memset(_eqSBD, 0, 3*sizeof(INTEGER*));
	}	
	~ModelSBD()
	{
		FreePairs();
	}


	double GetTheoreticalOrbits() { // dimension of the model
		double expected = pow(4, g_numSpecies)-3*(g_numSpecies-1);
		return expected;
	}

	void ShowResults() {
		if (!g_chainCharsSpeciesADN) {
			std::cout << "Error: No file Loaded." << std::endl;
			return;
		}
		if (!_calculated) {
			std::cout << "No computation done." << std::endl;
		}

		int option = 0;
			std::cout << "Number of unique observations (GMM orbits) " << _orbitals.size()
					<< std::endl;
			LONG_INTEGER expected = 1;

			expected=GetTheoreticalOrbits();

			std::cout << "Dimension of the SBD: " << expected  << std::endl;
				}


	bool IsCalculated() {
		return _calculated;
	}

	

	bool _calculated;

	// Private methods
	void CalcRulesAndValues(ModelGMM *pmodelGMM) {
			SaveOrbital(pmodelGMM);
			_calculated = true;
		}

	void GetEquations(ModelGMM *pmodelGMM) {
					for (int i = 0; i < pmodelGMM->_orbitals.size(); i++) {

    				for (int j = 0; j < _numSpecies; j++) {
    					
				    	LONG_INTEGER value = pmodelGMM->_orbitals[i] >> j * SHIFT_CORRIMENT & g_mask; // j is the species
						
						if (value<3)
						{
							 if(j==0){
								_eqSBD[value][_numClasses-1] -= (pmodelGMM->_counts[i]);

													}
							else if(j==(_numClasses)){
								_eqSBD[value][0] += (pmodelGMM->_counts[i]);  
								  }

							else {
							 	_eqSBD[value][_numClasses-j-1] -= (pmodelGMM->_counts[i]);   
								_eqSBD[value][_numClasses-j] += (pmodelGMM->_counts[i]);  
								  }				    
							}
    					}}
						
						}

    void SaveOrbital(ModelGMM *pmodelGMM) {
			CreatePairs();

			GetEquations(pmodelGMM);
		    int norm= 6 * pow(4, _numSpecies-2);
						for (int i = 0; i < pmodelGMM->_orbitals.size(); i++) 
						{
						mleSBD=0.0;
         					for (int j = 0; j < _numSpecies; j++) 
							{
	    			    	LONG_INTEGER value = (pmodelGMM->_orbitals[i]) >> j * SHIFT_CORRIMENT & g_mask; // j is the species
			                 	 
							if(value<3)
							{	 
								 if(j==0){
								 
		                         mleSBD = - _eqSBD[value][_numClasses-1];
								 }		
								 else if(j==(_numClasses))
								 {
								 mleSBD += _eqSBD[value][0];
								 }
								 else
								 {
								  mleSBD -= _eqSBD[value][j-1];
								  mleSBD += _eqSBD[value][j];
								 }
		    				}
							}
              
							mleSBD = ((double)pmodelGMM->_counts[i] - mleSBD/(double)norm)/(double)g_chainLength; 
						_orbitals.push_back(pmodelGMM->_counts[i]);
    				    _counts.push_back(mleSBD);
						_patterns.push_back(pmodelGMM->_orbitals[i]);
    				  }
					
					}
					
    
};

static ModelSBD *g_modelSBD = 0;


//*******************************************************
class ModelATR {
public:

	ModelATR(ModelGMM *pmodelGMM) :
		_calculated(false) {
		if (!pmodelGMM || !pmodelGMM->IsCalculated()) {
			std::cout << "Error: Model GMM not created or calculated.\n"
					<< std::endl;
			return;
		}
		memset(_eqATR, 0, 16*sizeof(INTEGER*));
		_numSpecies = g_numSpecies;
		_adnChainLength = g_chainLength;
		CalcRulesAndValues(pmodelGMM);

	}
~ModelATR()
	{
		FreePairs();
	}

	double GetTheoreticalOrbits() { // dimension of the model
		double expected= pow(4, _numSpecies) - 3.0*(_numSpecies)*(_numSpecies-1) ;
		return expected;
	}

    INTEGER _numSpecies;
	INTEGER _adnChainLength;
    INTEGER *_eqATR[4][4]; // for all 6 combinations of two leaves, A has 4 choices, C has 2 two choices, G has 1 choice  
	INTEGER _numClasses;
	
	void CreatePairs()
	{
		FreePairs();
		_numClasses = g_numSpecies*(g_numSpecies-1)/2;
		for(int i = 0; i < 4; i++)
    	{	
    		for(int j = 0; j < 4; j++)
       		{
         		_eqATR[i][j] = new INTEGER[_numClasses];  
         		for(int k = 0; k < _numClasses; k++)
         			_eqATR[i][j][k] = 0; 
       		}
    	}	
	}
	
	void FreePairs()
	{
		for(int i = 0; i < 4; i++)
    	{	
    		for(int j = 0; j < 4; j++)
       		{
       			if(_eqATR[i][j])
       			{
       				delete []_eqATR[i][j];
       			}   
       		}
    	}
	memset(_eqATR, 0, 16*sizeof(INTEGER*));
	}	

	void ShowResults() {
		
		if (!g_chainCharsSpeciesADN) {
			std::cout << "Error: No file Loaded." << std::endl;
			return;
		}
		if (!_calculated) {
			std::cout << "No computation done." << std::endl;
		}

		int option = 0;
		LONG_INTEGER expected = 1;

	    expected=GetTheoreticalOrbits();
	    std::cout << "Expected dimension of the ATR: " << expected  <<  std::endl;
		
					
		}

	bool IsCalculated() {
		return _calculated;
	}

	std::vector<LONG_INTEGER> _orbitals;
	std::vector<double> _counts;
    std::vector<INTEGER> _patterns;
	
	int index;
    double mleATR;
	bool _calculated;

	// Private methods
	void CalcRulesAndValues(ModelGMM *pmodelGMM) {
			//INTEGER ntVal=0; // letter a
			SaveOrbital(pmodelGMM);
			_calculated = true;
		}

        void GetEquations(ModelGMM *pmodelGMM) {
    	         CreatePairs();
         			
         			index = _numClasses-1; 
         			
    				for (int j = 0; j < _numSpecies; j++) 
					{
					for (int k = j+1; k < _numSpecies; k++) 
					{
						for (int i = 0; i < pmodelGMM->_orbitals.size(); i++) { //across the data

					    INTEGER value1 = pmodelGMM->_orbitals[i] >> j * SHIFT_CORRIMENT & g_mask; // j is the species
						INTEGER value2 = pmodelGMM->_orbitals[i] >> k * SHIFT_CORRIMENT & g_mask; // k 
						bool found = false;
                        
						if(value1!=value2){		
									
						if(value1 > value2) {
						  _eqATR[value2][value1][index] += (pmodelGMM->_counts[i]);
						}
						else 
						
						  	_eqATR[value1][value2][index] -= (pmodelGMM->_counts[i]);
							 }
				        }
    					index -= 1;
    				
					}
					}
					
		    }	
			void SaveOrbital(ModelGMM *pmodelGMM) {
					    GetEquations(pmodelGMM);
						int norm= 2 * pow(4, _numSpecies-2);

						for (int i = 0; i < pmodelGMM->_orbitals.size(); i++) 
						{ 
						mleATR=0.0;
						index = _numClasses-1; 

						for (int j = 0; j < _numSpecies; j++) 
						{
							for (int k = j+1; k < _numSpecies; k++) 
							{
							INTEGER value1 = pmodelGMM->_orbitals[i] >> j * SHIFT_CORRIMENT & g_mask; // j is the species
							INTEGER value2 = pmodelGMM->_orbitals[i] >> k * SHIFT_CORRIMENT & g_mask; // k	
				
							if(value1!=value2)
							{		

								if (value1>value2)
								{
								mleATR += _eqATR[value2][value1][index];
															//						std::cout <<  " mleATR1 " << mleATR << std::endl;

								}
								else
								{
						
								mleATR -= _eqATR[value1][value2][index];
															//						std::cout <<  " mleATR2 " << mleATR << std::endl;

								}
							}
							index -= 1;	
							
							}
						}
						
						mleATR = ((double)pmodelGMM->_counts[i] - mleATR/(double)norm)/(double)g_chainLength; 
						 						
						_orbitals.push_back(pmodelGMM->_counts[i]);
    				    _counts.push_back(mleATR);
						_patterns.push_back(pmodelGMM->_orbitals[i]);
				        }
						}
						


};

static ModelATR *g_modelATR = 0;
//*******************************************************

class ModelK81 {
public:

	ModelK81(ModelSMM *pmodelSMM) :
		_calculated(false) {
		if (!pmodelSMM || !pmodelSMM->IsCalculated()) {
			std::cout << "Error: Model SMM not created or calculated.\n"
					<< std::endl;
			return;
		}
		_numSpecies = g_numSpecies;
		_adnChainLength = g_chainLength;
		CalcRulesAndValues(pmodelSMM);
	}

	LONG_INTEGER GetTheoreticalOrbits() {
		LONG_INTEGER expected = 1;
		expected = pow(4, _numSpecies - 1);
		return expected;
	}
	void ShowResults() {
		if (!g_chainCharsSpeciesADN) {
			std::cout << "Error: No file Loaded." << std::endl;
			return;
		}
		if (!_calculated) {
			std::cout << "No computation done." << std::endl;
		}

		int option = -1;
		std::cout << " - (0) Show To Screen." << std::endl;
		std::cout << " - (1) Save To File." << std::endl;
		std::cin >> option;

		if (option == 0) {
			std::cout << "Total Options Found: " << _orbitals.size()
					<< std::endl;
			LONG_INTEGER expected = 1;
			for (int i = 1; i < _numSpecies; i++)
				expected *= 4;

			std::cout << "Total Expected: " << expected << std::endl;
			for (int i = 0; i < _orbitals.size(); i++) {
				std::cout
						<< "------------------------------------------------------------------"
						<< std::endl;
				for (int j = 0; j < _orbitals[i].size(); j++) {
					if (j == _baseOrb[i])
						std::cout << "P(" << transform_adn_chain_val_to_string(
								_orbitals[i][j]) << ") ** = " << _counts[i][j]
								<< std::endl;
					else
						std::cout << "P(" << transform_adn_chain_val_to_string(
								_orbitals[i][j]) << ")    = " << _counts[i][j]
								<< std::endl;
				}
			}
		} else if (option == 1) {
			std::string fileName;
			std::cout << "Name of the file to save: ";
			std::cin >> fileName;

			int pos;
			if ((pos = fileName.find(".dat")) != -1)
				fileName.replace(pos, 4, ".txt");
			if (fileName.find(".txt") == -1)
				fileName.append(".txt");

			// Save
			std::ofstream myfile;
			myfile.open(fileName.c_str());
			myfile << "Total Options Found: " << _orbitals.size() << std::endl;
			INTEGER expected = 1;
			for (int i = 1; i < _numSpecies; i++)
				expected *= 4;

			myfile << "Total Expected: " << expected << std::endl;
			for (int i = 0; i < _orbitals.size(); i++) {
				myfile
						<< "------------------------------------------------------------------"
						<< std::endl;
				for (int j = 0; j < _orbitals[i].size(); j++) {
					if (j == _baseOrb[i])
						myfile << "P(" << transform_adn_chain_val_to_string(
								_orbitals[i][j]) << ") ** = " << _counts[i][j]
								<< std::endl;
					else
						myfile << "P(" << transform_adn_chain_val_to_string(
								_orbitals[i][j]) << ")    = " << _counts[i][j]
								<< std::endl;
				}
			}
			myfile.close();
		}
	}

	bool IsCalculated() {
		return _calculated;
	}

	INTEGER _numSpecies;
	INTEGER _adnChainLength;

	std::vector<std::vector<LONG_INTEGER> > _orbitals;
	std::vector<std::vector<INTEGER> > _counts;
	std::vector<INTEGER> _baseOrb;

	bool _calculated;

	// Private methods
	void CalcRulesAndValues(ModelSMM *pmodelSMM) {
		LONG_INTEGER valOrbital;
		INTEGER *foundBefore = new INTEGER[pmodelSMM->_orbitals.size()];
		memset(foundBefore, 0, sizeof(INTEGER) * pmodelSMM->_orbitals.size());
		for (int i = 0; i < pmodelSMM->_orbitals.size(); i++) {
			if (!foundBefore[i]) {
				int pos = -1;
				SearchTwinOrbital(pmodelSMM, i, pos);
				if (pos != -1)
					foundBefore[pos] = 1;
				SaveOrbital(pmodelSMM, i, pos);
			}
		}
		delete foundBefore;
		_calculated = true;
	}

	void SearchTwinOrbital(ModelSMM *pmodelSMM, int posMain, int &posTwin) {
		posTwin = -1;
		LONG_INTEGER mainOrbVal = pmodelSMM->_orbitals[posMain][0];
		LONG_INTEGER twinOrbVal = change_adn_bases(mainOrbVal, 1, 0);

		bool bFound = false;
		for (int i = posMain + 1; i < pmodelSMM->_orbitals.size() && !bFound; i++) {
			if (pmodelSMM->_orbitals[i][0] == twinOrbVal) {
				posTwin = i;
				bFound = true;
			} else if (pmodelSMM->_orbitals[i][1] == twinOrbVal) {
				posTwin = i;
				bFound = true;
			}
		}
	}

	void SaveOrbital(ModelSMM *pmodelSMM, int posMain, int posTwin) {
		LONG_INTEGER mainOrbVal = pmodelSMM->_orbitals[posMain][0];

		std::vector<LONG_INTEGER> newOrbital;
		newOrbital.push_back(pmodelSMM->_orbitals[posMain][0]);
		newOrbital.push_back(pmodelSMM->_orbitals[posMain][1]);
		if (posTwin != -1) {
			newOrbital.push_back(pmodelSMM->_orbitals[posTwin][0]);
			newOrbital.push_back(pmodelSMM->_orbitals[posTwin][1]);
		} else {
			newOrbital.push_back(change_adn_bases(mainOrbVal, 1, 0));
			newOrbital.push_back(change_adn_bases(mainOrbVal, 2, 0));
		}
		_orbitals.push_back(newOrbital);

		std::vector<INTEGER> newCounts;
		newCounts.push_back(pmodelSMM->_counts[posMain][0]);
		newCounts.push_back(pmodelSMM->_counts[posMain][1]);
		if (posTwin != -1) {
			newCounts.push_back(pmodelSMM->_counts[posTwin][0]);
			newCounts.push_back(pmodelSMM->_counts[posTwin][1]);
		} else {
			newCounts.push_back(0);
			newCounts.push_back(0);
		}
		_counts.push_back(newCounts);
		_baseOrb.push_back(0);
		SortLastOrbital();
	}

	void SortLastOrbital() {
		INTEGER lastOrb = _orbitals.size() - 1;
		for (int i = 0; i < _orbitals[lastOrb].size(); i++) {
			for (int j = i + 1; j < _orbitals[lastOrb].size(); j++) {
				if (_orbitals[lastOrb][i] > _orbitals[lastOrb][j]) {
					LONG_INTEGER tempValOrb = _orbitals[lastOrb][i];
					_orbitals[lastOrb][i] = _orbitals[lastOrb][j];
					_orbitals[lastOrb][j] = tempValOrb;

					INTEGER tempCounts = _counts[lastOrb][i];
					_counts[lastOrb][i] = _counts[lastOrb][j];
					_counts[lastOrb][j] = tempCounts;

					if (i == _baseOrb[lastOrb])
						_baseOrb[lastOrb] = j;
					else if (j == _baseOrb[lastOrb])
						_baseOrb[lastOrb] = i;
				}
			}
		}
		int s = 0;
	}
};

static ModelK81 *g_modelK81 = 0;

class ModelK80 {
public:

	ModelK80(ModelK81 *pmodelK81) :
		_calculated(false) {
		if (!pmodelK81 || !pmodelK81->IsCalculated()) {
			std::cout << "Error: Model SMM not created or calculated.\n"
					<< std::endl;
			return;
		}
		_numSpecies = g_numSpecies;
		_adnChainLength = g_chainLength;
		CalcRulesAndValues(pmodelK81);
	}

	LONG_INTEGER GetTheoreticalOrbits() {
		LONG_INTEGER expected = pow(2.0, 2 * _numSpecies - 3) + pow(2.0,
				_numSpecies - 2);

		return expected;
	}

	void ShowResults() {
		if (!g_chainCharsSpeciesADN) {
			std::cout << "Error: No file Loaded." << std::endl;
			return;
		}
		if (!_calculated) {
			std::cout << "No computation done." << std::endl;
		}

		int option = -1;
		std::cout << " - (0) Show To Screen." << std::endl;
		std::cout << " - (1) Save To File." << std::endl;
		std::cin >> option;

		if (option == 0) {
			std::cout << "Total Options Found: " << _orbitals.size()
					<< std::endl;

			double expected = pow(2.0, 2 * _numSpecies - 3) + pow(2.0,
					_numSpecies - 2);
			std::cout << "Total Expected: " << expected << std::endl;
			for (int i = 0; i < _orbitals.size(); i++) {
				std::cout
						<< "------------------------------------------------------------------"
						<< std::endl;
				for (int j = 0; j < _orbitals[i].size(); j++) {
					std::cout << "P(" << transform_adn_chain_val_to_string(
							_orbitals[i][j]) << ")    = " << _counts[i][j]
							<< std::endl;
				}
			}
		} else if (option == 1) {
			std::string fileName;
			std::cout << "Name of the file to save: ";
			std::cin >> fileName;

			int pos;
			if ((pos = fileName.find(".dat")) != -1)
				fileName.replace(pos, 4, ".txt");
			if (fileName.find(".txt") == -1)
				fileName.append(".txt");

			// Save
			std::ofstream myfile;
			myfile.open(fileName.c_str());
			myfile << "Total Options Found: " << _orbitals.size() << std::endl;
			double expected = pow(2.0, 2 * _numSpecies - 3) + pow(2.0,
					_numSpecies - 2);

			myfile << "Total Expected: " << expected << std::endl;
			for (int i = 0; i < _orbitals.size(); i++) {
				myfile
						<< "------------------------------------------------------------------"
						<< std::endl;
				for (int j = 0; j < _orbitals[i].size(); j++) {
					myfile << "P(" << transform_adn_chain_val_to_string(
							_orbitals[i][j]) << ")    = " << _counts[i][j]
							<< std::endl;
				}
			}
			myfile.close();
		}
	}

	bool IsCalculated() {
		return _calculated;
	}

	INTEGER _numSpecies;
	INTEGER _adnChainLength;

	std::vector<std::vector<LONG_INTEGER> > _orbitals;
	std::vector<std::vector<INTEGER> > _counts;

	bool _calculated;

	// Private methods
	void CalcRulesAndValues(ModelK81 *pmodelK81) {
		INTEGER *foundBefore = new INTEGER[pmodelK81->_orbitals.size()];
		memset(foundBefore, 0, sizeof(INTEGER) * pmodelK81->_orbitals.size());
		for (int i = 0; i < pmodelK81->_orbitals.size(); i++) {
			if (!foundBefore[i]) {
				int pos = -1;
				SearchTwinOrbital(pmodelK81, i, pos);
				if (pos != -1 && pos != -2)
					foundBefore[pos] = 1;
				SaveOrbital(pmodelK81, i, pos);
			}
		}
		delete foundBefore;
		_calculated = true;
	}

	void SearchTwinOrbital(ModelK81 *pmodelK81, int posMain, int &posTwin) {
		posTwin = -1;
		LONG_INTEGER mainOrbVal = pmodelK81->_orbitals[posMain][0];
		LONG_INTEGER twinOrbVal = change_adn_bases_simple(mainOrbVal, 3, 0);

		if (twinOrbVal == pmodelK81->_orbitals[posMain][2]) {
			posTwin = -2;
			return;
		}

		bool bFound = false;
		for (int i = posMain + 1; i < pmodelK81->_orbitals.size() && !bFound; i++) {
			if (pmodelK81->_orbitals[i][2] == twinOrbVal) {
				posTwin = i;
				bFound = true;
			}
		}
	}

	void SaveOrbital(ModelK81 *pmodelK81, int posMain, int posTwin) {
		LONG_INTEGER mainOrbVal = pmodelK81->_orbitals[posMain][0];

		std::vector<LONG_INTEGER> newOrbital;
		std::vector<INTEGER> newCounts;
		for (int i = 0; i < pmodelK81->_orbitals[posMain].size(); i++) {
			newOrbital.push_back(pmodelK81->_orbitals[posMain][i]);
			newCounts.push_back(pmodelK81->_counts[posMain][i]);
		}

		if (posTwin != -1 && posTwin != -2) {
			for (int i = 0; i < pmodelK81->_orbitals[posTwin].size(); i++) {
				newOrbital.push_back(pmodelK81->_orbitals[posTwin][i]);
				newCounts.push_back(pmodelK81->_counts[posTwin][i]);
			}
		}

		if (posTwin == -1) {
			LONG_INTEGER twinOrbVal = change_adn_bases_simple(
					pmodelK81->_orbitals[posMain][0], 3, 0);
			LONG_INTEGER fourthOrbVal = change_adn_bases(twinOrbVal, 1, 0);
			LONG_INTEGER secondOrbVal = change_adn_bases(twinOrbVal, 0, 0);
			LONG_INTEGER firstOrbVal = change_adn_bases(twinOrbVal, 2, 0);

			newOrbital.push_back(firstOrbVal);
			newOrbital.push_back(secondOrbVal);
			newOrbital.push_back(twinOrbVal);
			newOrbital.push_back(fourthOrbVal);

			newCounts.push_back(0);
			newCounts.push_back(0);
			newCounts.push_back(0);
			newCounts.push_back(0);
		}

		_orbitals.push_back(newOrbital);
		_counts.push_back(newCounts);
	}
};

static ModelK80 *g_modelK80 = 0;

class ModelJC {
public:

	ModelJC(ModelK80 *pmodelK80) :
		_calculated(false) {
		if (!pmodelK80 || !pmodelK80->IsCalculated()) {
			std::cout << "Error: Model K80 not created or calculated.\n"
					<< std::endl;
			return;
		}
		_numSpecies = g_numSpecies;
		_adnChainLength = g_chainLength;
		CalcRulesAndValues(pmodelK80);
	}

	double GetTheoreticalOrbits() {
		double expected = (1.0 / 3.0) * pow(2.0, 2 * _numSpecies - 3) + pow(
				2.0, _numSpecies - 2) + 1.0 / 3.0;

		return expected;
	}

	void ShowResults() {
		if (!g_chainCharsSpeciesADN) {
			std::cout << "Error: No file Loaded." << std::endl;
			return;
		}
		if (!_calculated) {
			std::cout << "No computation done." << std::endl;
		}

		int option = -1;
		std::cout << " - (0) Show To Screen." << std::endl;
		std::cout << " - (1) Save To File." << std::endl;
		std::cin >> option;

		if (option == 0) {
			std::cout << "Total Options Found: " << _orbitals.size()
					<< std::endl;

			double expected = (1.0 / 3.0) * pow(2.0, 2 * _numSpecies - 3)
					+ pow(2.0, _numSpecies - 2) + 1.0 / 3.0;
			std::cout << "Total Expected: " << expected << std::endl;
			for (int i = 0; i < _orbitals.size(); i++) {
				std::cout
						<< "------------------------------------------------------------------"
						<< std::endl;
				for (int j = 0; j < _orbitals[i].size(); j++) {
					std::cout << "P(" << transform_adn_chain_val_to_string(
							_orbitals[i][j]) << ")    = " << _counts[i][j]
							<< std::endl;
				}
			}
		} else if (option == 1) {
			std::string fileName;
			std::cout << "Name of the file to save: ";
			std::cin >> fileName;

			int pos;
			if ((pos = fileName.find(".dat")) != -1)
				fileName.replace(pos, 4, ".txt");
			if (fileName.find(".txt") == -1)
				fileName.append(".txt");

			// Save
			std::ofstream myfile;
			myfile.open(fileName.c_str());
			myfile << "Total Options Found: " << _orbitals.size() << std::endl;
			double expected = (1.0 / 3.0) * pow(2.0, 2 * _numSpecies - 3)
					+ pow(2.0, _numSpecies - 2) + 1.0 / 3.0;

			myfile << "Total Expected: " << expected << std::endl;
			for (int i = 0; i < _orbitals.size(); i++) {
				myfile
						<< "------------------------------------------------------------------"
						<< std::endl;
				for (int j = 0; j < _orbitals[i].size(); j++) {
					myfile << "P(" << transform_adn_chain_val_to_string(
							_orbitals[i][j]) << ")    = " << _counts[i][j]
							<< std::endl;
				}
			}
			myfile.close();
		}
	}

	bool IsCalculated() {
		return _calculated;
	}

	INTEGER _numSpecies;
	INTEGER _adnChainLength;

	std::vector<std::vector<LONG_INTEGER> > _orbitals;
	std::vector<std::vector<INTEGER> > _counts;

	bool _calculated;

	// Private methods
	void CalcRulesAndValues(ModelK80 *pmodelK80) {
		INTEGER *foundBefore = new INTEGER[pmodelK80->_orbitals.size()];
		memset(foundBefore, 0, sizeof(INTEGER) * pmodelK80->_orbitals.size());
		for (int i = 0; i < pmodelK80->_orbitals.size(); i++) {
			if (!foundBefore[i]) {
				int pos1 = -1;
				int pos2 = -1;
				SearchTwinOrbital(pmodelK80, i, pos1, pos2);
				if (pos1 != -1 && pos1 != -2) {
					foundBefore[pos1] = 1;
				}
				if (pos2 != -1 && pos2 != -2) {
					foundBefore[pos2] = 1;
				}
				SaveOrbital(pmodelK80, i, pos1, pos2);
			}
		}
		delete foundBefore;
		_calculated = true;
	}

	void SearchTwinOrbital(ModelK80 *pmodelK80, int posMain, int &posTwin1,
			int &posTwin2) {
		posTwin1 = -1;
		posTwin2 = -1;
		LONG_INTEGER mainOrbVal = pmodelK80->_orbitals[posMain][0];
		LONG_INTEGER twinOrbVal = change_adn_bases_simple(mainOrbVal, 4, 0);
		LONG_INTEGER twinOrbVal2 = change_adn_bases_simple(mainOrbVal, 5, 0);
		if (pmodelK80->_orbitals[posMain].size() == 4) {
			if (twinOrbVal == pmodelK80->_orbitals[posMain][1]) {
				// We are inside the same orbital nothing must be done
				posTwin1 = -1;
				posTwin2 = -1;
				return;
			}
			posTwin1 = -2;
			posTwin2 = -2;
			for (int i = posMain + 1; i < pmodelK80->_orbitals.size(); i++) {
				if (pmodelK80->_orbitals[i].size() == 4) {
					if (pmodelK80->_orbitals[i][1] == twinOrbVal) {
						posTwin1 = i;
						if (posTwin2 >= 0)
							return;
					}
					if (pmodelK80->_orbitals[i][3] == twinOrbVal2) {
						posTwin2 = i;
						if (posTwin1 >= 0)
							return;
					}
				} else if (pmodelK80->_orbitals[i].size() == 8) {
					if (pmodelK80->_orbitals[i][1] == twinOrbVal
							|| pmodelK80->_orbitals[i][5] == twinOrbVal) {
						posTwin1 = i;
						posTwin2 = -1;
						return;
					}
				} else {
					std::cout << "Orbital in K80 model of size: "
							<< pmodelK80->_orbitals[i].size() << std::endl;
					return;
				}
			}
		} else if (pmodelK80->_orbitals[posMain].size() == 8) {
			posTwin1 = -2;
			posTwin2 = -2;
			for (int i = posMain + 1; i < pmodelK80->_orbitals.size(); i++) {
				if (pmodelK80->_orbitals[i].size() == 4) {
					if (pmodelK80->_orbitals[i][1] == twinOrbVal) {
						posTwin1 = i;
						posTwin2 = -1;
						return;
					}
					if (pmodelK80->_orbitals[i][3] == twinOrbVal2) {
						posTwin1 = -1;
						posTwin2 = i;
					}
				}
				if (pmodelK80->_orbitals[i].size() == 8) {
					LONG_INTEGER twinOrbVal2 = change_adn_bases_simple(
							mainOrbVal, 5, 0);
					if (pmodelK80->_orbitals[i][1] == twinOrbVal
							|| pmodelK80->_orbitals[i][5] == twinOrbVal) {
						posTwin1 = i;
						if (posTwin2 >= 0)
							return;
					}
					if (pmodelK80->_orbitals[i][3] == twinOrbVal2
							|| pmodelK80->_orbitals[i][7] == twinOrbVal2) {
						posTwin2 = i;
						if (posTwin1 >= 0)
							return;
					}
				}
			}
		}
	}

	void SaveOrbital(ModelK80 *pmodelK80, int posMain, int posTwin1,
			int posTwin2) {
		LONG_INTEGER mainOrbVal = pmodelK80->_orbitals[posMain][0];
		std::vector<LONG_INTEGER> newOrbital;
		std::vector<INTEGER> newCounts;

		// Copy the previous orbital
		for (int i = 0; i < pmodelK80->_orbitals[posMain].size(); i++) {
			newOrbital.push_back(pmodelK80->_orbitals[posMain][i]);
			newCounts.push_back(pmodelK80->_counts[posMain][i]);
		}

		if (pmodelK80->_orbitals[posMain].size() == 4) {
			if (posTwin1 >= 0 && posTwin2 >= 0) {
				// Copy both 4 Orbitals
				for (int i = 0; i < pmodelK80->_orbitals[posTwin1].size(); i++) {
					newOrbital.push_back(pmodelK80->_orbitals[posTwin1][i]);
					newCounts.push_back(pmodelK80->_counts[posTwin1][i]);
				}
				for (int i = 0; i < pmodelK80->_orbitals[posTwin2].size(); i++) {
					newOrbital.push_back(pmodelK80->_orbitals[posTwin2][i]);
					newCounts.push_back(pmodelK80->_counts[posTwin2][i]);
				}
			} else if (posTwin1 >= 0 && posTwin2 == -1) {
				// Copy 8 Orbital
				for (int i = 0; i < pmodelK80->_orbitals[posTwin1].size(); i++) {
					newOrbital.push_back(pmodelK80->_orbitals[posTwin1][i]);
					newCounts.push_back(pmodelK80->_counts[posTwin1][i]);
				}
			} else if (posTwin1 >= 0 && posTwin2 == -2) {
				// Copy 4 orbital and create the other
				for (int i = 0; i < pmodelK80->_orbitals[posTwin1].size(); i++) {
					newOrbital.push_back(pmodelK80->_orbitals[posTwin1][i]);
					newCounts.push_back(pmodelK80->_counts[posTwin1][i]);
				}
				Create4Orbital(mainOrbVal, newOrbital, 5);
				for (int i = 0; i < 4; i++)
					newCounts.push_back(0);
			} else if (posTwin1 == -2 && posTwin2 >= 0) {
				// Copy 4 orbital and create the other
				for (int i = 0; i < pmodelK80->_orbitals[posTwin2].size(); i++) {
					newOrbital.push_back(pmodelK80->_orbitals[posTwin2][i]);
					newCounts.push_back(pmodelK80->_counts[posTwin2][i]);
				}
				Create4Orbital(mainOrbVal, newOrbital, 4);
				for (int i = 0; i < 4; i++)
					newCounts.push_back(0);
			} else if (posTwin1 == -2 && posTwin2 == -2) {
				AddOrbitalsTo4Orbital(posTwin1, posTwin2, mainOrbVal,
						newOrbital, newCounts);
			}
		} else if (pmodelK80->_orbitals[posMain].size() == 8) {
			if (posTwin1 >= 0 && posTwin2 == -1) {
				// Copy 4 orbital
				for (int i = 0; i < pmodelK80->_orbitals[posTwin1].size(); i++) {
					newOrbital.push_back(pmodelK80->_orbitals[posTwin1][i]);
					newCounts.push_back(pmodelK80->_counts[posTwin1][i]);
				}
			} else if (posTwin1 == -1 && posTwin2 >= 0) {
				// Copy 4 orbital
				for (int i = 0; i < pmodelK80->_orbitals[posTwin2].size(); i++) {
					newOrbital.push_back(pmodelK80->_orbitals[posTwin2][i]);
					newCounts.push_back(pmodelK80->_counts[posTwin2][i]);
				}
			} else if (posTwin1 >= 0 && posTwin2 >= 0) {
				// Copy 8 orbital
				for (int i = 0; i < pmodelK80->_orbitals[posTwin1].size(); i++) {
					newOrbital.push_back(pmodelK80->_orbitals[posTwin1][i]);
					newCounts.push_back(pmodelK80->_counts[posTwin1][i]);
				}
				// Copy 8 orbital
				for (int i = 0; i < pmodelK80->_orbitals[posTwin2].size(); i++) {
					newOrbital.push_back(pmodelK80->_orbitals[posTwin2][i]);
					newCounts.push_back(pmodelK80->_counts[posTwin2][i]);
				}
			} else if (posTwin1 >= 0 && posTwin2 == -2) {
				// Copy 8 orbital
				for (int i = 0; i < pmodelK80->_orbitals[posTwin1].size(); i++) {
					newOrbital.push_back(pmodelK80->_orbitals[posTwin1][i]);
					newCounts.push_back(pmodelK80->_counts[posTwin1][i]);
				}
				// Create 8 Orbital
				Create8Orbital(mainOrbVal, newOrbital, 5);
				for (int i = 0; i < 8; i++) {
					newCounts.push_back(0);
				}
			} else if (posTwin1 == -2 && posTwin2 >= 0) {
				// Copy 8 orbital
				for (int i = 0; i < pmodelK80->_orbitals[posTwin2].size(); i++) {
					newOrbital.push_back(pmodelK80->_orbitals[posTwin2][i]);
					newCounts.push_back(pmodelK80->_counts[posTwin2][i]);
				}
				// Create 8 Orbital
				Create8Orbital(mainOrbVal, newOrbital, 4);
				for (int i = 0; i < 8; i++) {
					newCounts.push_back(0);
				}
			} else if (posTwin1 == -2 && posTwin2 == -2) {
				AddOrbitalsTo8Orbital(posTwin1, posTwin2, mainOrbVal,
						newOrbital, newCounts);
			}
		}

		_orbitals.push_back(newOrbital);
		_counts.push_back(newCounts);
	}

	void AddOrbitalsTo8Orbital(INTEGER posTwin1, INTEGER posTwin2,
			LONG_INTEGER mainOrbVal, std::vector<LONG_INTEGER> &newOrbital,
			std::vector<INTEGER> &newCounts) {
		LONG_INTEGER acGenOrbVal2 = change_adn_bases_simple(mainOrbVal, 4, 0);
		LONG_INTEGER orbVal3 = change_adn_bases(acGenOrbVal2, 0, 0);
		LONG_INTEGER orbVal1 = change_adn_bases(acGenOrbVal2, 1, 0);
		LONG_INTEGER orbVal4 = change_adn_bases(acGenOrbVal2, 2, 0);
		if (acGenOrbVal2 == newOrbital[1] || acGenOrbVal2 == newOrbital[5]) {
			LONG_INTEGER atGenOrbVal4 = change_adn_bases_simple(mainOrbVal, 5,
					0);
			orbVal1 = change_adn_bases(atGenOrbVal4, 0, 0);
			orbVal3 = change_adn_bases(atGenOrbVal4, 1, 0);
			LONG_INTEGER orbVal2 = change_adn_bases(atGenOrbVal4, 2, 0);

			// This is a 4 orbital
			newOrbital.push_back(orbVal1);
			newOrbital.push_back(orbVal2);
			newOrbital.push_back(orbVal3);
			newOrbital.push_back(atGenOrbVal4);

			for (int i = 0; i < 4; i++)
				newCounts.push_back(0);

			return;
		}

		LONG_INTEGER agGenOrbVal3 = change_adn_bases_simple(acGenOrbVal2, 3, 0);

		if (agGenOrbVal3 == orbVal4 || agGenOrbVal3 == acGenOrbVal2) {
			// This is a 4 orbital
			newOrbital.push_back(orbVal1);
			newOrbital.push_back(acGenOrbVal2);
			newOrbital.push_back(orbVal3);
			newOrbital.push_back(orbVal4);

			for (int i = 0; i < 4; i++)
				newCounts.push_back(0);
		} else {
			Create8Orbital(mainOrbVal, newOrbital, 4);
			for (int i = 0; i < 8; i++)
				newCounts.push_back(0);
			Create8Orbital(mainOrbVal, newOrbital, 5);
			for (int i = 0; i < 8; i++)
				newCounts.push_back(0);
		}
	}

	void AddOrbitalsTo4Orbital(INTEGER posTwin1, INTEGER posTwin2,
			LONG_INTEGER mainOrbVal, std::vector<LONG_INTEGER> &newOrbital,
			std::vector<INTEGER> &newCounts) {
		LONG_INTEGER acGenOrbVal2 = change_adn_bases_simple(mainOrbVal, 4, 0);
		LONG_INTEGER orbVal3 = change_adn_bases(acGenOrbVal2, 0, 0);
		LONG_INTEGER orbVal1 = change_adn_bases(acGenOrbVal2, 1, 0);
		LONG_INTEGER orbVal4 = change_adn_bases(acGenOrbVal2, 2, 0);

		LONG_INTEGER agGenOrbVal3 = change_adn_bases_simple(acGenOrbVal2, 3, 0);
		if (agGenOrbVal3 != orbVal4 && agGenOrbVal3 != acGenOrbVal2) {
			LONG_INTEGER orbVal6 = change_adn_bases(agGenOrbVal3, 0, 0);
			LONG_INTEGER orbVal8 = change_adn_bases(agGenOrbVal3, 1, 0);
			LONG_INTEGER orbVal5 = change_adn_bases(agGenOrbVal3, 2, 0);

			newOrbital.push_back(orbVal1);
			newOrbital.push_back(acGenOrbVal2);
			newOrbital.push_back(orbVal3);
			newOrbital.push_back(orbVal4);

			newOrbital.push_back(orbVal5);
			newOrbital.push_back(orbVal6);
			newOrbital.push_back(agGenOrbVal3);
			newOrbital.push_back(orbVal8);

			for (int i = 0; i < 8; i++)
				newCounts.push_back(0);
		} else {
			LONG_INTEGER atGenOrbVal4 = change_adn_bases_simple(mainOrbVal, 5,
					0);
			LONG_INTEGER orbVal5 = change_adn_bases(atGenOrbVal4, 0, 0);
			LONG_INTEGER orbVal7 = change_adn_bases(atGenOrbVal4, 1, 0);
			LONG_INTEGER orbVal6 = change_adn_bases(atGenOrbVal4, 2, 0);

			newOrbital.push_back(orbVal1);
			newOrbital.push_back(acGenOrbVal2);
			newOrbital.push_back(orbVal3);
			newOrbital.push_back(orbVal4);

			newOrbital.push_back(orbVal5);
			newOrbital.push_back(orbVal6);
			newOrbital.push_back(orbVal7);
			newOrbital.push_back(atGenOrbVal4);

			for (int i = 0; i < 8; i++)
				newCounts.push_back(0);
		}
	}

	void Create4Orbital(LONG_INTEGER mainOrbVal,
			std::vector<LONG_INTEGER> &newOrbital, int rule) {
		LONG_INTEGER orbVal1 = change_adn_bases_simple(mainOrbVal, rule, 0);
		LONG_INTEGER orbVal2 = change_adn_bases(orbVal1, 0, 0);
		LONG_INTEGER orbVal3 = change_adn_bases(orbVal1, 1, 0);
		LONG_INTEGER orbVal4 = change_adn_bases(orbVal1, 2, 0);

		newOrbital.push_back(orbVal1);
		newOrbital.push_back(orbVal2);
		newOrbital.push_back(orbVal3);
		newOrbital.push_back(orbVal4);

		// Order the orbital
		for (int i = newOrbital.size() - 4; i < newOrbital.size(); i++) {
			for (int j = i + 1; j < newOrbital.size(); j++) {
				if (newOrbital[i] > newOrbital[j]) {
					LONG_INTEGER temp = newOrbital[i];
					newOrbital[i] = newOrbital[j];
					newOrbital[j] = temp;
				}
			}
		}
	}

	void Create8Orbital(LONG_INTEGER mainOrbVal,
			std::vector<LONG_INTEGER> &newOrbital, int rule) {
		LONG_INTEGER orbVal1 = change_adn_bases_simple(mainOrbVal, rule, 0);
		LONG_INTEGER orbVal2 = change_adn_bases(orbVal1, 0, 0);
		LONG_INTEGER orbVal3 = change_adn_bases(orbVal1, 1, 0);
		LONG_INTEGER orbVal4 = change_adn_bases(orbVal1, 2, 0);

		newOrbital.push_back(orbVal1);
		newOrbital.push_back(orbVal2);
		newOrbital.push_back(orbVal3);
		newOrbital.push_back(orbVal4);

		// Order the values introduced
		for (int i = newOrbital.size() - 4; i < newOrbital.size(); i++) {
			for (int j = i + 1; j < newOrbital.size(); j++) {
				if (newOrbital[i] > newOrbital[j]) {
					LONG_INTEGER temp = newOrbital[i];
					newOrbital[i] = newOrbital[j];
					newOrbital[j] = temp;
				}
			}
		}

		LONG_INTEGER orbVal7 = change_adn_bases_simple(orbVal1, 3, 0);
		LONG_INTEGER orbVal5 = change_adn_bases(orbVal7, 1, 0);
		LONG_INTEGER orbVal6 = change_adn_bases(orbVal7, 0, 0);
		LONG_INTEGER orbVal8 = change_adn_bases(orbVal7, 2, 0);

		newOrbital.push_back(orbVal5);
		newOrbital.push_back(orbVal6);
		newOrbital.push_back(orbVal7);
		newOrbital.push_back(orbVal8);

		// Order the values introduced
		for (int i = newOrbital.size() - 4; i < newOrbital.size(); i++) {
			for (int j = i + 1; j < newOrbital.size(); j++) {
				if (newOrbital[i] > newOrbital[j]) {
					LONG_INTEGER temp = newOrbital[i];
					newOrbital[i] = newOrbital[j];
					newOrbital[j] = temp;
				}
			}
		}
	}
};

static ModelJC *g_modelJC = 0;
//************************************************************************
int main(int argc, char* argv[]) {
	// Get initial constant data
	MAX_NUM_SPECIES = MaxNumSpecies();
	SHIFT_CORRIMENT = GetShiftCorriment();
	SetGlobalMask();

	if (argc >= 2) {
		std::string fileName = argv[1];
		std::cout << "Loading file " << fileName << ".." << std::endl;
		load_file(fileName);
		std::cout << "Calculating rules under the General Markov Model... " << std::endl;
		if (g_modelGMM) {
			delete g_modelGMM;
			g_modelGMM = 0;
		}
		g_modelGMM = new ModelGMM();

		//std::cout << "Calculating rules under the Stable Base Distribution model... " << std::endl;
		//if (g_modelSBD) {
		//delete g_modelSBD;
		//g_modelSBD = 0;
		//}
		//g_modelSBD = new ModelSBD(g_modelGMM);

		//std::cout << "Calculating rules under the Algebraic Time Reversible model... " << std::endl;
		//if (g_modelATR) {
		//			delete g_modelATR;
		//			g_modelATR = 0;
		//		}
		//g_modelATR = new ModelATR(g_modelGMM);


	std::cout << "Calculating rules under the Strand Symmetric model... " << std::endl;
		if (g_modelSMM) {
			delete g_modelSMM;
			g_modelSMM = 0;
		}
		g_modelSMM = new ModelSMM(g_modelGMM);

		std::cout << "Calculating rules under the Kimura 3-parameter model... " << std::endl;
		if (g_modelK81) {
			delete g_modelK81;
			g_modelK81 = 0;
		}
		g_modelK81 = new ModelK81(g_modelSMM);

		std::cout << "Calculating rules under the Kimura 2-parameter model... " << std::endl;
		if (g_modelK80) {
			delete g_modelK80;
			g_modelK80 = 0;
		}
		g_modelK80 = new ModelK80(g_modelK81);

		std::cout << "Calculating rules under the Jukes-Cantor model... " << std::endl;
		if (g_modelJC) {
			delete g_modelJC;
			g_modelJC = 0;
		}
		g_modelJC = new ModelJC(g_modelK80);
        
		double bicSBD;
		double bicATR;
		double bicSMM;
		double bicK81;
		double bicK80;
		double bicJC;
		
		//double aicSBD = computeAICSBDATR( g_modelSBD->_orbitals, g_modelSBD->_counts,
				//		g_modelSBD->GetTheoreticalOrbits(), g_chainLength, bicSBD);
		//double aicATR = computeAICSBDATR( g_modelATR->_orbitals, g_modelATR->_counts,
				//g_modelATR->GetTheoreticalOrbits(), g_chainLength, bicATR);
		double aicSMM = computeAIC(g_modelSMM->_counts,
				g_modelSMM->GetTheoreticalOrbits(), g_chainLength,  bicSMM);
		double aicK81 = computeAIC(g_modelK81->_counts,
				g_modelK81->GetTheoreticalOrbits(), g_chainLength,  bicK81);
		double aicK80 = computeAIC(g_modelK80->_counts,
				g_modelK80->GetTheoreticalOrbits(), g_chainLength,  bicK80);
		double aicJC = computeAIC(g_modelJC->_counts,
				g_modelJC->GetTheoreticalOrbits(), g_chainLength,  bicJC);

		//std::cout << "AIC score SBD: " << aicSBD << std::endl;
		//std::cout << "AIC score ATR: " << aicATR << std::endl;
		std::cout << "AIC score SMM: " << aicSMM << std::endl;
		std::cout << "AIC score K81: " << aicK81 << std::endl;
		std::cout << "AIC score K80: " << aicK80 << std::endl;
		std::cout << "AIC score JC: " << aicJC << std::endl;
		int pos = -1;
		pos = fileName.find(".fa");
		fileName.replace(pos, 6, ".ResAIC");
		//} else {
//			fileName.append("ResAIC");
		//}
		fileName.append(".txt");

		// Save
		std::ofstream myfile;
		myfile.open(fileName.c_str());
		//myfile << aicSBD << std::endl;
		//myfile << aicATR << std::endl;
		myfile << aicSMM << std::endl;
		myfile << aicK81 << std::endl;
		myfile << aicK80 << std::endl;
		myfile << aicJC << std::endl;
		myfile.close();

        //pos = -1;
		//if ((pos = fileName.find(".fa")) != -1) {
		pos = fileName.find("AIC");
	    fileName.replace(pos, 3, "BIC");
	//	} else {
	//		fileName.append("ResBIC");
	//	}
		

        myfile.open(fileName.c_str());
		//myfile << bicSBD << std::endl;
		//myfile << bicATR << std::endl;
		myfile << bicSMM << std::endl;
		myfile << bicK81 << std::endl;
		myfile << bicK80 << std::endl;
		myfile << bicJC << std::endl;
		myfile.close();
 

	} else {
		int optionChosen = EXIT;

		do {
			std::cout << "Choose an option." << std::endl;
			std::cout << " - (0) Create a File." << std::endl;
			std::cout << " - (1) Load a File." << std::endl;
			std::cout << " - (2) Calc Rules and Values." << std::endl;
			std::cout << " - (3) Show Results." << std::endl;
			std::cout << " - (4) Operate." << std::endl;
			std::cout << " - (5) Exit." << std::endl;
			std::cin >> optionChosen;
			switch (optionChosen) {
			case CREATE_FILE: {
				int numSpecies;
				int chainLength;
				std::string fileName;
				std::cout << "Creating File..." << std::endl;
				std::cout << "Number of species (max " << MAX_NUM_SPECIES
						<< "): ";
				std::cin >> numSpecies;
				std::cout << "Chains Length: ";
				std::cin >> chainLength;
				std::cout << "Name of the file: ";
				std::cin >> fileName;
				create_file(fileName, numSpecies, chainLength);
				break;
			}
			case LOAD_FILE: {
				std::string fileName;
				std::cout << "Loading File.." << std::endl;
				std::cout << "Name of the file: ";
				std::cin >> fileName;
				load_file(fileName);
				break;
			}
			case CALC_RULES_VALUES: {
				if (g_modelGMM) {
					delete g_modelGMM;
					g_modelGMM = 0;
				}
				g_modelGMM = new ModelGMM();
				std::cout << "Model GMM specified." << std::endl;

				if (g_modelSBD) {
					delete g_modelSBD;
					g_modelSBD = 0;
					}
			    g_modelSBD = new ModelSBD(g_modelGMM);
				std::cout << "Model SBD specified." << std::endl;

				if (g_modelATR) {
					delete g_modelATR;
					g_modelATR = 0;
					}
				g_modelATR = new ModelATR(g_modelGMM);
				std::cout << "Model ATR specified." << std::endl;

				if (g_modelSMM) {
					delete g_modelSMM;
					g_modelSMM = 0;
				}
				g_modelSMM = new ModelSMM(g_modelGMM);
				std::cout << "Model SMM specified." << std::endl;

				if (g_modelK81) {
					delete g_modelK81;
					g_modelK81 = 0;
				}
				g_modelK81 = new ModelK81(g_modelSMM);
				std::cout << "Model K81 specified." << std::endl;

				if (g_modelK80) {
					delete g_modelK80;
					g_modelK80 = 0;
				}
				g_modelK80 = new ModelK80(g_modelK81);
				std::cout << "Model K80 specified." << std::endl;
				if (g_modelJC) {
					delete g_modelJC;
					g_modelJC = 0;
				}
				g_modelJC = new ModelJC(g_modelK80);
				std::cout << "Model JC69 specified." << std::endl;
				break;
			}
			case SHOW_RESULTS: {
				if (g_modelGMM)
					g_modelGMM->ShowResults();
				if (g_modelSBD)
					g_modelSBD->ShowResults();
				if (g_modelATR)
					g_modelATR->ShowResults();
				if (g_modelSMM)
					g_modelSMM->ShowResults();
				if (g_modelK81)
					g_modelK81->ShowResults();
				if (g_modelK80)
					g_modelK80->ShowResults();
				if (g_modelJC)
					g_modelJC->ShowResults();
				break;
			}
			case OPERATE: {
			
			  
		double bicSBD;
		double bicATR;
		double bicSMM;
		double bicK81;
		double bicK80;
		double bicJC;

				//computeAIC(std::vector<std::vector<int> > &counts, double K, double n)
				double aicSBD = computeAICSBDATR( g_modelSBD->_orbitals, g_modelSBD->_counts,
						g_modelSBD->GetTheoreticalOrbits(), g_chainLength,  bicSBD);
				double aicATR = computeAICSBDATR( g_modelATR->_orbitals, g_modelATR->_counts,
						g_modelATR->GetTheoreticalOrbits(), g_chainLength,  bicATR);
				double aicSMM = computeAIC(g_modelSMM->_counts,
						g_modelSMM->GetTheoreticalOrbits(), g_chainLength,  bicSMM);
				double aicK81 = computeAIC(g_modelK81->_counts,
						g_modelK81->GetTheoreticalOrbits(), g_chainLength,  bicK81);
				double aicK80 = computeAIC(g_modelK80->_counts,
						g_modelK80->GetTheoreticalOrbits(), g_chainLength,  bicK80);
				double aicJC = computeAIC(g_modelJC->_counts,
						g_modelJC->GetTheoreticalOrbits(), g_chainLength,  bicJC);
				std::cout << "AIC score under the SBD model: " << aicSBD << std::endl;
				std::cout << "AIC score under the ATR model: " << aicATR << std::endl;
				std::cout << "AIC score under the SMM model: " << aicSMM << std::endl;
				std::cout << "AIC score under the K81 model: " << aicK81 << std::endl;
				std::cout << "AIC score under the K80 model: " << aicK80 << std::endl;
				std::cout << "AIC score under the JC69 model: " << aicJC << std::endl;
				break;
			}
			case EXIT: {
				break;
			}
			default: {
				std::cout << "Incorrect selection, try again." << std::endl;
				break;
			}
			}
		} while (optionChosen != EXIT);
	}

	// Free Memory
	FreeMemory();

	return -1;
}
