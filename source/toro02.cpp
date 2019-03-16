//#include <iostream>
//#include <fstream>
//#include <vector>
//#include <string>
//#include <sstream>
//#include <iomanip>
//#include <map>
//#include <cstdlib>
//#include <iterator>
//#include <algorithm>
//#include <cmath>
//
//#ifndef WIDTH	
//#define WIDTH 16
//#endif // !WIDTH
//#ifndef PGA
//#define PGA 100.0
//#endif // !PGA
//
//#ifndef PRECISION
//#define PRECISION 8
//#endif // !PRECISION
//#ifndef PRECISION2
//#define PRECISION2 3
//#endif // !PRECISION2
//#ifndef e
//#define e 2.7182818284590
//#endif // !e
//
//
//
//using namespace std;
//
//double closest(std::vector<double> const& vec, double value) {
//	auto const it = std::lower_bound(vec.begin(), vec.end(), value);
//	if (it == vec.end()) { return -1; }
//	return *it;
//}
//
//int main() {
//
//	typedef std::vector<vector<double>> Matrix;
//	typedef std::vector<vector<vector<double>>> Matrix2;
//	typedef std::vector<double> Vector;
//	double R;
//	double MINF; //Lower limit og magnitude given in the tabl
//	double MSUP; //Upper limit of magnitude given in the table
//	int NMAG; //Number of magnitudes for which intensity is given
//	double DMAG;
//	double RINF; //Lower limit of distance given in the table
//	double RSUP; //Upper limit of distance given in the table
//	int NRAD; //Number of distances for which intensity is given
//	int TYPEMD; //An integer indicating the type of distance used by the attenuation table
//	double DLRAD;
//	double AMAX = 0; //Type of stadistic distribution
//	double sigmadepth, sigmamodeling;
//	double STRESSI;//Inicial stress
//	char OPTION1, OPTION2;
//
//	//--- ATTENUATION VARIABLES ----
//	Vector periodsreq;
//	Vector frequenciesreq;
//	double aux1, aux2;
//	double A, B, C, D, E, F;
//	//-------------------------------
//
//	ifstream periods;
//	periods.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/toro02/periods.txt");
//
//	//If we can read/write great
//	while (periods.good())
//	{
//		periods >> aux1;
//		periodsreq.push_back(aux1);
//	}
//
//#if 0
//	for (size_t i = 0; i < periodsreq.size(); i++)
//	{
//		cout << periodsreq[i] << endl;
//	}
//
//#endif // 0
//
//	for (size_t i = 0; i < periodsreq.size(); i++) // Changing periods to frequencies
//	{
//		aux2 = 1 / periodsreq.at(i);
//		frequenciesreq.push_back(aux2);
//	}
//#if 0
//	cout << periodsreq.size() << " " << frequenciesreq.size() << endl;
//	for (size_t i = 0; i < frequenciesreq.size(); i++)
//	{
//		cout << frequenciesreq[i] << endl;
//	}
//#endif // 1
//
//
//	// ---------------------------------- ASKING FOR DATA ------------------------------------------------
//
//	//Enter the minimum magnitude, maximum and the number of intermediate quantities
//	MINF = 5.0;
//	MSUP = 8.0;
//	NMAG = 16;
//	//Enter the minimum distance, maximum and the number of intermediate distances
//	RINF = 1.0;
//	RSUP = 500.0;
//	NRAD = 100;
//
//	cout << "CHOOSE THE TYPE OF COEFFICIENTS: " << endl;
//	cout << "-Midcontinent, equations using Moment Magnitude (A)" << endl;
//	cout << "-Gulfcontinent, equations using Moment Magnitude (B)" << endl;
//	cin >> OPTION1;
//	cout << endl;
//	cout << "CHOOSE THE TYPE OF APPROACH: " << endl;
//	cout << "-Empirical approach - Rjb (A)" << endl;
//	cout << "-Modelling approach - Rrup (B)" << endl;
//	cin >> OPTION2;
//
//	// -------------------------------------- FINISH ------------------------------------------------------
//
//
//
//	//-------------------------------- GENERATING COEFICENTS -----------------------------------------------
//	Matrix results(periodsreq.size());//Vertical size of Matrix
//	Matrix aleatoryMvalues(periodsreq.size());
//	Matrix aleatoryDvalues(periodsreq.size());
//	Matrix table2(8);//T97, Table 2 required
//	Matrix table3(8);//T97, Table 3
//	Matrix table4(8);//T97, Table 4
//	Vector coficentesreq(15);
//	Vector aleatoriesreq_M(3);
//	Vector aleatoriesreq_D(2);
//	double aleatoryuncertainty;
//	double epistemicuncertainty;
//	double standarddesviation;
//	int index1 = 0, index2 = 0;//Positions to interpolate
//
//	ifstream toro2002;
//	ifstream aleatorymagnitude;
//	ifstream aleatorydistance;
//	ofstream coeficentsT2; // Archive
//	Vector f2(8);
//	Vector f3(8);
//	Vector f4(8);
//	Vector coefT2(16);
//	Vector coefT3(4);
//	Vector coefT4(3);
//	double nearfrequency, nearfrequency2;
//	double required, required2;
//
//#if 1
//
//	switch (OPTION1)
//	{
//	case 'A':
//		toro2002.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/toro02/table2mcmm.txt");
//		if (toro2002.good())
//		{
//			for (size_t i = 0; i < 8; i++) {
//				for (size_t j = 0; j < 16; j++)	toro2002 >> coefT2.at(j);
//				table2.at(i) = coefT2;
//			}
//		}
//		break;
//	case 'B':
//		toro2002.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/toro02/table2gulfmm.txt");
//		if (toro2002.good())
//		{
//			for (size_t i = 0; i < 8; i++) {
//				for (size_t j = 0; j < 16; j++)	toro2002 >> coefT2.at(j);
//				table2.at(i) = coefT2;
//			}
//		}
//		break;
//
//	}
//
//
//
//	for (size_t i = 0; i < 8; i++)
//		f2.at(i) = table2.at(i).at(0);
//
//#if 0
//	cout << "Table 2 (Toro et al., 2002)" << endl;
//	for (size_t i = 0; i < 8; i++) {
//		for (size_t j = 0; j < 16; j++)	cout << setw(WIDTH) << table2.at(i).at(j);
//		cout << endl;
//	}
//	cout << endl;
//
//	for (size_t i = 0; i < frequenciesreq.size(); i++)
//		cout << frequenciesreq.at(i) << endl;
//
//
//#endif // 0
//
//	for (size_t i = 0; i < frequenciesreq.size(); i++)
//	{
//		required = frequenciesreq[i];
//		if (required >= PGA) {
//			periodsreq[i] = 0.0;
//			for (size_t j = 0; j < 15; j++)	coficentesreq.at(j) = table2.at(7).at(j + 1);
//			results.at(i) = coficentesreq;
//		}
//		else {
//			nearfrequency = closest(f2, required);
//			std::vector<double>::iterator it = std::find(f2.begin(), f2.end(), nearfrequency);
//			index1 = std::distance(f2.begin(), it);
//
//			if ((nearfrequency <= required))
//			{
//				if (index1 == 7) { index2 = index1 - 1; }
//				else { index2 = index1 + 1; }
//			}
//
//			if ((nearfrequency >= required))
//			{
//				if (index1 == 0) { index2 = index1 + 1; }
//				else { index2 = index1 - 1; }
//			}
//
//#if 0
//			cout << index1 << " " << index2 << endl;
//#endif // 0
//
//
//			for (size_t j = 0; j < 15; j++)
//			{
//				coficentesreq.at(j) = (((frequenciesreq[i] - table2.at(index1).at(0))*table2.at(index2).at(j + 1)) + ((table2.at(index2).at(0) - frequenciesreq[i])*table2.at(index1).at(j + 1))) / (table2.at(index2).at(0) - table2.at(index1).at(0));
//			}
//			results.at(i) = coficentesreq;
//		}
//	}
//
//#if 1
//
//	coeficentsT2.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/toro02/resultsT2.dat");
//	coeficentsT2 << "#Coefficents for differents frequencies (table2mcmm - Toro et al., 2002)" << endl;
//	coeficentsT2 << setw(WIDTH) << "c1" << setw(WIDTH) << "c2" << setw(WIDTH) << "c4" << setw(WIDTH) << "c5" << setw(WIDTH) << "c6" << setw(WIDTH) << "c7" << setw(WIDTH) << "c8" << setw(WIDTH) << "c10" << setw(WIDTH) << "P. sigma" << setw(WIDTH) << "T. sigma" << endl;
//
//	for (size_t i = 0; i < periodsreq.size(); i++)
//	{
//
//		for (size_t j = 0; j < 15; j++)
//			coeficentsT2 << setw(WIDTH) << setprecision(6) << results.at(i).at(j);
//		coeficentsT2 << endl;
//	}
//
//
//#endif // 0
//
//	cout << "NICE" << endl;
//#endif // 0
//
//#if 1
//	aleatorymagnitude.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/toro02/table3.txt");
//	if (aleatorymagnitude.good())
//	{
//		for (size_t i = 0; i < 8; i++) {
//			for (size_t j = 0; j < 4; j++)	aleatorymagnitude >> coefT3.at(j);
//			table3.at(i) = coefT3;
//		}
//	}
//
//	aleatorydistance.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/toro02/table4.txt");
//	if (aleatorydistance.good())
//	{
//		for (size_t i = 0; i < 8; i++) {
//			for (size_t j = 0; j < 3; j++)	aleatorydistance >> coefT4.at(j);
//			table4.at(i) = coefT4;
//		}
//	}
//
//#if 1
//	cout << "Table 4 (Toro et al., 2002)" << endl;
//	for (size_t i = 0; i < 8; i++) {
//		for (size_t j = 0; j < 3; j++)	cout << setw(WIDTH) << table4.at(i).at(j);
//		cout << endl;
//	}
//	cout << endl;
//
//#endif // 0
//
//
//	for (size_t i = 0; i < 8; i++)
//	{
//		f3.at(i) = table3.at(i).at(0);
//	}
//	for (size_t i = 0; i < 8; i++)
//	{
//		f4.at(i) = table4.at(i).at(0);
//	}
//
//
//	for (size_t i = 0; i < frequenciesreq.size(); i++)
//	{
//		required = frequenciesreq[i];
//		if (required >= PGA) {
//			for (size_t j = 0; j < 3; j++)	aleatoriesreq_M.at(j) = table3.at(7).at(j + 1);
//			aleatoryMvalues.at(i) = aleatoriesreq_M;
//		}
//		else {
//			nearfrequency = closest(f3, required);
//			std::vector<double>::iterator it2 = std::find(f3.begin(), f3.end(), nearfrequency);
//			index1 = std::distance(f3.begin(), it2);
//
//			if ((nearfrequency <= required))
//			{
//				if (index1 == 7) { index2 = index1 - 1; }
//				else { index2 = index1 + 1; }
//			}
//
//			if ((nearfrequency >= required))
//			{
//				if (index1 == 0) { index2 = index1 + 1; }
//				else { index2 = index1 - 1; }
//			}
//
//#if 0
//			cout << index1 << " " << index2 << endl;
//#endif // 0
//
//
//			for (size_t j = 0; j < 3; j++)
//			{
//				aleatoriesreq_M.at(j) = (((frequenciesreq[i] - table3.at(index1).at(0))*table3.at(index2).at(j + 1)) + ((table3.at(index2).at(0) - frequenciesreq[i])*table3.at(index1).at(j + 1))) / (table3.at(index2).at(0) - table3.at(index1).at(0));
//			}
//			aleatoryMvalues.at(i) = aleatoriesreq_M;
//		}
//	}
//
//
//
//	for (size_t i = 0; i < frequenciesreq.size(); i++)
//	{
//		required2 = frequenciesreq[i];
//		if (required2 >= PGA) {
//			for (size_t j = 0; j < 2; j++)	aleatoriesreq_D.at(j) = table4.at(7).at(j + 1);
//			aleatoryDvalues.at(i) = aleatoriesreq_D;
//		}
//		else {
//			nearfrequency2 = closest(f4, required2);
//			std::vector<double>::iterator it4 = std::find(f4.begin(), f4.end(), nearfrequency2);
//			index1 = std::distance(f4.begin(), it4);
//
//			if ((nearfrequency2 <= required2))
//			{
//				if (index1 == 7) { index2 = index1 - 1; }
//				else { index2 = index1 + 1; }
//			}
//
//			if ((nearfrequency2 >= required2))
//			{
//				if (index1 == 0) { index2 = index1 + 1; }
//				else { index2 = index1 - 1; }
//			}
//
//#if 0
//			cout << index1 << " " << index2 << endl;
//#endif // 0
//
//
//			for (size_t j = 0; j < 2; j++)
//			{
//				aleatoriesreq_D.at(j) = (((frequenciesreq[i] - table4.at(index1).at(0))*table4.at(index2).at(j + 1)) + ((table4.at(index2).at(0) - frequenciesreq[i])*table4.at(index1).at(j + 1))) / (table4.at(index2).at(0) - table4.at(index1).at(0));
//			}
//			aleatoryDvalues.at(i) = aleatoriesreq_D;
//		}
//	}
//
//
//
//#if 0
//
//	coeficentsT2.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/toro02/resultsT2.dat");
//	coeficentsT2 << "#Coefficents for differents frequencies (table2mcmm - Toro et al., 2002)" << endl;
//	coeficentsT2 << setw(WIDTH) << "c1" << setw(WIDTH) << "c2" << setw(WIDTH) << "c4" << setw(WIDTH) << "c5" << setw(WIDTH) << "c6" << setw(WIDTH) << "c7" << setw(WIDTH) << "c8" << setw(WIDTH) << "c10" << setw(WIDTH) << "P. sigma" << setw(WIDTH) << "T. sigma" << endl;
//
//	for (size_t i = 0; i < periodsreq.size(); i++)
//	{
//
//		for (size_t j = 0; j < 15; j++)
//			coeficentsT2 << setw(WIDTH) << setprecision(6) << results.at(i).at(j);
//		coeficentsT2 << endl;
//	}
//
//
//#endif // 0
//
//
//#endif // 1
//
//	//--------------------------------------- FINISH -------------------------------------------------------
//
//	//-------------------------- CHECKING ----------------------------------------
//#if 0
//	for (size_t i = 0; i < frequenciesreq.size(); i++)
//	{
//		for (size_t j = 0; j < 15; j++)
//			cout << setw(WIDTH) << setprecision(6) << results.at(i).at(j);
//		cout << endl;
//	}
//
//#endif // 0
//	//----------------------- FINISH CHECKING -------------------------------------
//
//	//---------------------- GENERATING VALUES -------------------------------
//
//	Vector aceleraciones(NRAD);
//	Vector distances(NRAD); //Rjb
//	Vector magnitudes(NMAG);
//	DMAG = (MSUP - MINF) / (NMAG - 1);
//	DLRAD = (log10(RSUP) - log10(RINF)) / (NRAD - 1);
//
//
//	for (size_t i = 0; i < NRAD; i++)
//	{
//		distances[i] = pow(10, log10(RINF) + i*DLRAD);
//	}
//	distances[NRAD - 1] = RSUP; //Changing last value
//	for (size_t i = 0; i < NMAG; i++)
//	{
//		magnitudes[i] = MINF + i*DMAG;
//	}
//	magnitudes[NMAG - 1] = MSUP;//Changing last value
//
//								//---------------------------- FINISH -----------------------------------
//								//---------------------------- OUTPUT -----------------------------------
//								//Value Type of distance
//								//	1 (or blank) Focal
//								//	2 Epicentral
//								//	3 Joyner and Boore
//								//	4 Closest to rupture area(Rrup)
//	TYPEMD = 3;
//#if 1
//	cout << "MAGNITUDES:" << endl;
//	for (size_t i = 0; i < NMAG; i++)
//	{
//		cout << magnitudes[i] << endl;
//	}
//	cout << endl;
//#endif // 0
//#if 1
//	cout << "DISTANCES:" << endl;
//	for (size_t i = 0; i < NRAD; i++)
//	{
//		cout << setw(WIDTH) << distances[i];
//	}
//	cout << endl;
//#endif // 0
//
//	//OUTPUT 1
//	ofstream attenueationtableT02;
//
//	attenueationtableT02.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/toro02/toro02.atn");
//
//	attenueationtableT02 << setprecision(PRECISION2);
//	//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Description" << setw(WIDTH) << ": Sample attenuation file constructed for illustration purposes (2008)" << endl;	
//	//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Units" << setw(WIDTH) << ": cm/sec/sec" << endl;
//	//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Distribution" << setw(WIDTH) << ": 2" << endl;
//	//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Dimension" << setw(WIDTH) << ": Aceleration" << endl;
//	attenueationtableT02 << setw(WIDTH) << MINF << setw(WIDTH) << MSUP << setw(WIDTH) << NMAG << endl;
//	attenueationtableT02 << setw(WIDTH) << RINF << setw(WIDTH) << RSUP << setw(WIDTH) << NRAD << setw(WIDTH) << TYPEMD << endl;
//
//	attenueationtableT02 << setprecision(PRECISION);
//
//	for (size_t i = 0; i < periodsreq.size(); i++)//Loop over periods
//	{
//		attenueationtableT02 << setw(WIDTH) << periodsreq.at(i) << setw(WIDTH) << "0.7" << setw(WIDTH) << AMAX << endl;
//		for (size_t j = 0; j < NMAG; j++)//Loop over magnitudes
//		{
//
//			if (magnitudes[j]<5.5)
//			{
//				sigmamodeling = aleatoryMvalues.at(i).at(1) - ((aleatoryMvalues.at(i).at(1) - aleatoryMvalues.at(i).at(0))*(5.5 - magnitudes[j])) / 0.5;
//			}
//			else
//			{
//				sigmamodeling = aleatoryMvalues.at(i).at(2) - ((aleatoryMvalues.at(i).at(2) - aleatoryMvalues.at(i).at(1))*(8 - magnitudes[j])) / 2.5;
//			}
//
//
//			if (periodsreq.at(i) <= 1.0)
//			{
//				epistemicuncertainty = 0.36 + 0.07*(magnitudes[j] - 6);
//			}
//			if (periodsreq.at(i) == 2.0)
//			{
//				epistemicuncertainty = 0.34 + 0.06*(magnitudes[j] - 6);
//
//			}
//			if ((periodsreq.at(i) > 1.0) && (periodsreq.at(i) < 2.0))
//			{
//				epistemicuncertainty = 0.34 + 0.06*(magnitudes[j] - 6) - (((0.34 + 0.06*(magnitudes[j] - 6)) - (0.36 + 0.07*(magnitudes[j] - 6)))*(2 - periodsreq.at(i))) / 1;
//			}
//
//			for (size_t k = 0; k < NRAD; k++)//Loop over coeficents
//			{
//				switch (OPTION2) {
//				case 'A':
//					R = pow(pow(distances[k], 2) + (pow(results.at(i).at(14), 2))*(pow(pow(e, (-1.25 + 0.227*magnitudes[j])), 2)), 0.5);
//					break;
//				case 'B':
//					R = distances[k] + 0.089*pow(e, 0.6*magnitudes[j]);
//					break;
//				}
//				if (distances[k] <= 5.0)
//				{
//					sigmadepth = aleatoryDvalues.at(i).at(0);
//
//				}
//				if ((distances[k] > 5) && (distances[k] < 20))
//				{
//					sigmadepth = (15 * (aleatoryDvalues.at(i).at(1)) - (aleatoryDvalues.at(i).at(1) - aleatoryDvalues.at(i).at(0))*(20 - distances[k])) / 15;
//				}
//				if (distances[k] >= 20)
//				{
//					sigmadepth = aleatoryDvalues.at(i).at(1);
//				}
//
//				R = pow(pow(distances[k], 2) + pow(results.at(i).at(14), 2), 0.5);
//				A = results.at(i).at(0);
//				B = results.at(i).at(1)*(magnitudes[j] - 6.0);
//				C = results.at(i).at(10)*pow((magnitudes[j] - 6.0), 2);
//				D = results.at(i).at(11)*log(R);
//				E = (results.at(i).at(12) - results.at(i).at(11))*max(log(R / 100.0), 0.0);
//				F = results.at(i).at(13)*R;
//				aleatoryuncertainty = pow(pow(sigmadepth, 2) + pow(sigmamodeling, 2), 0.5);
//				standarddesviation = pow(pow(aleatoryuncertainty, 2) + pow(epistemicuncertainty, 2), 0.5);
//				aceleraciones[k] = pow(e, A + B + C - D - E - F - standarddesviation);
//
//				attenueationtableT02 << setw(WIDTH) << aceleraciones[k];//Saving values
//				
//			}
//			attenueationtableT02 << endl;
//		}
//	}
//
//	system("pause");
//	return 0;
//}