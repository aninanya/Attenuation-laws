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
//	double sigma;
//	double AMAX = 0; //Type of stadistic distribution
//	double V30 = 800;//Shear wave velocity
//	double STRESS=120;//Stress asked
//	char OPTION;
//					
////--- ATTENUATION VARIABLES ----
//	Vector periodsreq;
//	Vector frequenciesreq;
//	double aux1, aux2;
//	double A, B, C, D, E, F;
//	//-------------------------------
//
//	ifstream periods;
//	periods.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/darragh15/periods.txt");
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
//	//::std::cout << "Enter the minimum magnitude, maximum and the number of intermediate quantities: " << endl;
//	//cin >> MINF >> MSUP >> NMAG;
//	MINF = 4.5;
//	MSUP = 8.5;
//	NMAG = 16;
//	//cout << "Enter the minimum distance, maximum and the number of intermediate distances: " << endl;
//	//cin >> RINF >> RSUP >> NRAD;
//	RINF = 1.0;
//	RSUP = 1000.0;
//	NRAD = 10000;
//
//	cout << "CHOOSE THE TYPE OF REGRESSION COEFFICIENTS: "<<endl;
//	cout << "-Single corner model with variable medium stress parameter as a function of Moment Magnitude (A)" << endl;
//	cout << "-Single corner model with constant medium stress parameter as a function of Moment Magnitude (B)" << endl;
//	cout << "-Double corner model with variable medium stress parameter as a function of Moment Magnitude (C)" << endl;
//	cout << "-Double corner model with constant medium stress parameter as a function of Moment Magnitude (D)" << endl;
//	cin >> OPTION;
//
//	// -------------------------------------- FINISH ------------------------------------------------------
//
//
//
//	//-------------------------------- GENERATING COEFICENTS -----------------------------------------------
//	Matrix results(periodsreq.size());//Vertical size of Matrix
//	Matrix table3(24);//DARRAGH15, Table 3.5
//	Vector coficentesreq(10);
//	int index1, index2;//Positions to interpolate
//
//	ifstream darragh2015;
//	ofstream coefT3; // Archive
//	Vector tfrequencies(24);//Vector of table 3.5a's frequencies
//	Vector coef(11);
//	double nearfrequency;
//	double required;
//
//
//#if 1
//	switch (OPTION)
//	{
//	case 'A':
//		darragh2015.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/darragh15/table3-5a.txt");
//		if (darragh2015.good())
//		{
//			for (size_t i = 0; i < 24; i++) {
//				for (size_t j = 0; j < 11; j++)	darragh2015 >> coef.at(j);
//				table3.at(i) = coef;
//			}
//		}
//		
//		for (size_t i = 0; i < 24; i++)
//			tfrequencies.at(i) = table3.at(i).at(0);
//		break;
//	
//	case 'B':
//		darragh2015.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/darragh15/table3-5b.txt");
//		if (darragh2015.good())
//		{
//			for (size_t i = 0; i < 24; i++) {
//				for (size_t j = 0; j < 11; j++)	darragh2015 >> coef.at(j);
//				table3.at(i) = coef;
//			}
//		}
//
//		for (size_t i = 0; i < 24; i++)
//			tfrequencies.at(i) = table3.at(i).at(0);
//		break;
//	case 'C':
//		darragh2015.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/darragh15/table3-5c.txt");
//		if (darragh2015.good())
//		{
//			for (size_t i = 0; i < 24; i++) {
//				for (size_t j = 0; j < 11; j++)	darragh2015 >> coef.at(j);
//				table3.at(i) = coef;
//			}
//		}
//
//		for (size_t i = 0; i < 24; i++)
//			tfrequencies.at(i) = table3.at(i).at(0);
//		break;
//	case 'D':
//		darragh2015.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/darragh15/table3-5d.txt");
//		if (darragh2015.good())
//		{
//			for (size_t i = 0; i < 24; i++) {
//				for (size_t j = 0; j < 11; j++)	darragh2015 >> coef.at(j);
//				table3.at(i) = coef;
//			}
//		}
//
//		for (size_t i = 0; i < 24; i++)
//			tfrequencies.at(i) = table3.at(i).at(0);
//		break;
//	}
//
//#if 0
//	cout << "Table 3.5 (Darragh et al., 2015)" << endl;
//	for (size_t i = 0; i < 24; i++) {
//		for (size_t j = 0; j < 11; j++)	cout << setw(WIDTH) << table3.at(i).at(j);
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
//		if (required > PGA) {
//			for (size_t j = 0; j < 10; j++)	coficentesreq.at(j) = table3.at(23).at(j + 1);
//			results.at(i) = coficentesreq;
//			periodsreq[i] = 0.0;
//		}
//		else {
//			nearfrequency = closest(tfrequencies, required);
//			std::vector<double>::iterator it = std::find(tfrequencies.begin(), tfrequencies.end(), nearfrequency);
//			index1 = std::distance(tfrequencies.begin(), it);
//
//			if ((nearfrequency <= required))
//			{
//				if (index1 == 23) { index2 = index1 - 1; }
//				else { index2 = index1 + 1; }
//			}
//
//			if ((nearfrequency >= required))
//			{
//				if (index1 == 0) { index2 = index1 + 1; }
//				else { index2 = index1 - 1; }
//			}
//
//
//
//#if 0
//			cout << index1 << " " << index2 << endl;
//#endif // 0
//
//
//			for (size_t j = 0; j < 10; j++)
//			{
//				coficentesreq.at(j) = (((frequenciesreq[i] - table3.at(index1).at(0))*table3.at(index2).at(j + 1)) + ((table3.at(index2).at(0) - frequenciesreq[i])*table3.at(index1).at(j + 1))) / (table3.at(index2).at(0) - table3.at(index1).at(0));
//			}
//			results.at(i) = coficentesreq;
//		}
//	}
//
//#if 1
//
//	coefT3.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/darragh15/resultsT3.dat");
//	coefT3 << "#Coefficents for differents frequencies (table 3.5 - Darragh et al., 2015)" << endl;
//	coefT3 << setw(WIDTH) << "c1" << setw(WIDTH) << "c2" << setw(WIDTH) << "c4" << setw(WIDTH) << "c5" << setw(WIDTH) << "c6" << setw(WIDTH) << "c7" << setw(WIDTH) << "c8" << setw(WIDTH) << "c10" << setw(WIDTH) << "P. sigma" << setw(WIDTH) << "T. sigma" << endl;
//
//	for (size_t i = 0; i < periodsreq.size(); i++)
//	{
//
//		for (size_t j = 0; j < 10; j++)
//			coefT3 << setw(WIDTH) << setprecision(6) << results.at(i).at(j);
//		coefT3 << endl;
//	}
//
//
//#endif // 0
//
//	cout << "NICE" << endl;
//#endif // 0
//
//	
//	//--------------------------------------- FINISH -------------------------------------------------------
//
//	//-------------------------- CHECKING ----------------------------------------
//#if 1
//	for (size_t i = 0; i < frequenciesreq.size(); i++)
//	{
//		for (size_t j = 0; j < 10; j++)
//			cout << setw(WIDTH) << setprecision(6) << results.at(i).at(j);
//		cout << endl;
//	}
//#endif // 0
//	//----------------------- FINISH CHECKING -------------------------------------
//
//	//---------------------- GENERATING VALUES -------------------------------
//
//	Vector aceleraciones(NRAD);
//	Vector distances(NRAD);
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
//	//---------------------------- FINISH -----------------------------------
//
//	//---------------------------- OUTPUT -----------------------------------
//	//Value Type of distance
//	//	1 (or blank) Focal
//	//	2 Epicentral
//	//	3 Joyner and Boore
//	//	4 Closest to rupture area(Rrup)
//	TYPEMD = 3;
//#if 0
//	cout << "MAGNITUDES:" << endl;
//	for (size_t i = 0; i < NMAG; i++)
//	{
//		cout << magnitudes[i] << endl;
//	}
//	cout << endl;
//#endif // 0
//#if 0
//	cout << "DISTANCES:" << endl;
//	for (size_t i = 0; i < NRAD; i++)
//	{
//		cout << distances[i] << endl;
//	}
//	cout << endl;
//#endif // 0
//
//	//OUTPUT 1
//	ofstream attenueationtableD15;
//
//	switch (OPTION)
//	{
//	case 'A':
//		attenueationtableD15.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/darragh15/darragh15_SC-VS-S.atn");
//		break;
//	case 'B':
//		attenueationtableD15.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/darragh15/darragh15_SC-CS-S.atn");
//		break;
//	case 'C':
//		attenueationtableD15.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/darragh15/darragh15_DC-VS-S.atn");
//		break;
//	case 'D':
//		attenueationtableD15.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/darragh15/darragh15_DC-CS-S.atn");
//		break;
//	}
//
//	attenueationtableD15 << setprecision(PRECISION2);
//	//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Description" << setw(WIDTH) << ": Sample attenuation file constructed for illustration purposes (2008)" << endl;	
//	//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Units" << setw(WIDTH) << ": cm/sec/sec" << endl;
//	//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Distribution" << setw(WIDTH) << ": 2" << endl;
//	//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Dimension" << setw(WIDTH) << ": Aceleration" << endl;
//	attenueationtableD15 << setw(WIDTH) << MINF << setw(WIDTH) << MSUP << setw(WIDTH) << NMAG << endl;
//	attenueationtableD15 << setw(WIDTH) << RINF << setw(WIDTH) << RSUP << setw(WIDTH) << NRAD << setw(WIDTH) << TYPEMD << endl;
//
//	attenueationtableD15 << setprecision(PRECISION);
//
//	for (size_t i = 0; i < periodsreq.size(); i++)//Loop over periods
//	{
//		sigma = results.at(i).at(9);
//		attenueationtableD15 << setw(WIDTH) << periodsreq.at(i) << setw(WIDTH) << sigma << setw(WIDTH) << AMAX << endl;
//		for (size_t j = 0; j < NMAG; j++)//Loop over magnitudes
//		{
//
//
//			for (size_t k = 0; k < NRAD; k++)//Loop over coeficents
//			{
//				R = distances[k];
//				A = results.at(i).at(0);
//				B = results.at(i).at(1)*magnitudes[j];
//				C = results.at(i).at(4) + results.at(i).at(5)*magnitudes[j];
//				D = log(R + pow(e, results.at(i).at(2)));
//				E = (results.at(i).at(7))*pow((magnitudes[j] - 6), 2.0);
//				F = results.at(i).at(6)*R;
//				aceleraciones[k] = pow(e, A + B + C*D + E);
//
//				attenueationtableD15 << setw(WIDTH) << aceleraciones[k];//Saving values
//			}
//			attenueationtableD15 << endl;
//		}
//	}
//	system("pause");
//	return 0;
//}