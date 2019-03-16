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
//	double STRESSI;//Inicial stress
//	char OPTION;
//					
////--- ATTENUATION VARIABLES ----
//	Vector periodsreq;
//	Vector frequenciesreq;
//	double aux1, aux2;
//	double A, B, C, D, E;
//	//-------------------------------
//
//	ifstream periods;
//	periods.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/silva02/periods.txt");
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
//	RSUP = 400.0;
//	NRAD = 10000;
//
//	cout << "CHOOSE THE TYPE OF REGRESSION COEFFICIENTS: "<<endl;
//	cout << "-Single corner with variable stress drop (A)" << endl;
//	cout << "-Single corner with constant stress drop (B)" << endl;
//	cout << "-Single corner model with constant stress drop and saturation (C)" << endl;
//	cout << "-Double corner model (D)" << endl;
//	cout << "-Double corner model with saturation (E)" << endl;
//	cout << "-Soil BC (F)" << endl;
//	
//	cin >> OPTION;
//
//	// -------------------------------------- FINISH ------------------------------------------------------
//
//
//
//	//-------------------------------- GENERATING COEFICENTS -----------------------------------------------
//	Matrix results(periodsreq.size());//Vertical size of Matrix
//	Matrix table3(27);//SILVA02, Table 3
//	Matrix table4(27);//SILVA02, Table 4
//	Matrix table5(27);//SILVA02, Table 5
//	Matrix table6(27);//SILVA02, Table 6
//	Matrix table7(27);//SILVA02, Table 7
//	Matrix tableBC(8);//SILVA02, Table BC
//	Vector coficentesreq(10);
//	int index1, index2;//Positions to interpolate
//
//	ifstream silva2002;
//	ofstream coefT3; // Archive
//	ofstream coefT4; // Archive
//	ofstream coefT5; // Archive
//	ofstream coefT6; // Archive
//	ofstream coefT7; // Archive
//	ofstream coefTBC; // Archive
//	Vector f3(27);//Vector of table 3's frequencies
//	Vector f4(27);//Vector of table 4's frequencies
//	Vector f5(27);//Vector of table 5's frequencies
//	Vector f6(27);//Vector of table 6's frequencies
//	Vector f7(27);//Vector of table 7's frequencies
//	Vector fBC(8);//Vector of table BC's frequencies
//	Vector coef(11);
//	double nearfrequency;
//	double required;
//
//	switch (OPTION)
//	{
//	case 'A':
//		silva2002.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/silva02/table3.txt");
//		if (silva2002.good())
//		{
//			for (size_t i = 0; i < 27; i++) {
//				for (size_t j = 0; j < 11; j++)	silva2002 >> coef.at(j);
//				table3.at(i) = coef;
//			}
//		}
//
//
//		for (size_t i = 0; i < 27; i++)
//			f3.at(i) = table3.at(i).at(0);
//
//#if 0
//		cout << "Table 3 (Silva et al., 2002)" << endl;
//		for (size_t i = 0; i < 27; i++) {
//			for (size_t j = 0; j < 11; j++)	cout << setw(WIDTH) << table3.at(i).at(j);
//			cout << endl;
//		}
//		cout << endl;
//
//		for (size_t i = 0; i < frequenciesreq.size(); i++)
//			cout << frequenciesreq.at(i) << endl;
//
//
//#endif // 0
//
//		for (size_t i = 0; i < frequenciesreq.size(); i++)
//		{
//			required = frequenciesreq[i];
//			if (required > PGA) {
//				for (size_t j = 0; j < 10; j++)	coficentesreq.at(j) = table3.at(26).at(j + 1);
//				results.at(i) = coficentesreq;
//				periodsreq[i] = 0.0;
//			}
//			else {
//				nearfrequency = closest(f3, required);
//				std::vector<double>::iterator it = std::find(f3.begin(), f3.end(), nearfrequency);
//				index1 = std::distance(f3.begin(), it);
//
//				if ((nearfrequency <= required))
//				{
//					if (index1 == 26) { index2 = index1 - 1; }
//					else { index2 = index1 + 1; }
//				}
//
//				if ((nearfrequency >= required))
//				{
//					if (index1 == 0) { index2 = index1 + 1; }
//					else { index2 = index1 - 1; }
//				}
//
//
//
//#if 0
//				cout << index1 << " " << index2 << endl;
//#endif // 0
//
//
//				for (size_t j = 0; j < 10; j++)
//				{
//					coficentesreq.at(j) = (((frequenciesreq[i] - table3.at(index1).at(0))*table3.at(index2).at(j + 1)) + ((table3.at(index2).at(0) - frequenciesreq[i])*table3.at(index1).at(j + 1))) / (table3.at(index2).at(0) - table3.at(index1).at(0));
//				}
//				results.at(i) = coficentesreq;
//			}
//		}
//
//#if 1
//
//		coefT3.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/silva02/resultsT3.dat");
//		coefT3 << "#Coefficents for differents frequencies (Table3 - Silva et al., 2002)" << endl;
//		coefT3 << setw(WIDTH) << "c1" << setw(WIDTH) << "c2" << setw(WIDTH) << "c4" << setw(WIDTH) << "c5" << setw(WIDTH) << "c6" << setw(WIDTH) << "c7" << setw(WIDTH) << "c8" << setw(WIDTH) << "c10" << setw(WIDTH)<<"P. sigma" << setw(WIDTH)<<"T. sigma" << endl;
//
//		for (size_t i = 0; i < periodsreq.size(); i++)
//		{
//
//			for (size_t j = 0; j < 10; j++)
//				coefT3 << setw(WIDTH) << setprecision(6) << results.at(i).at(j);
//			coefT3 << endl;
//		}
//
//
//#endif // 0
//
//		cout << "NICE" << endl;
//		break;
//	case 'B':
//		silva2002.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/silva02/table4.txt");
//		if (silva2002.good())
//		{
//			for (size_t i = 0; i < 27; i++) {
//				for (size_t j = 0; j < 11; j++)	silva2002 >> coef.at(j);
//				table4.at(i) = coef;
//			}
//		}
//
//
//		for (size_t i = 0; i < 27; i++)
//			f4.at(i) = table4.at(i).at(0);
//
//#if 0
//		cout << "Table 4 (Silva et al., 2002)" << endl;
//		for (size_t i = 0; i < 27; i++) {
//			for (size_t j = 0; j < 11; j++)	cout << setw(WIDTH) << table4.at(i).at(j);
//			cout << endl;
//		}
//		cout << endl;
//
//		for (size_t i = 0; i < frequenciesreq.size(); i++)
//			cout << frequenciesreq.at(i) << endl;
//
//
//#endif // 0
//
//		for (size_t i = 0; i < frequenciesreq.size(); i++)
//		{
//			required = frequenciesreq[i];
//			if (required > PGA) {
//				for (size_t j = 0; j < 10; j++)	coficentesreq.at(j) = table4.at(26).at(j + 1);
//				results.at(i) = coficentesreq;
//				periodsreq[i] = 0.0;
//			}
//			else {
//				nearfrequency = closest(f4, required);
//				std::vector<double>::iterator it = std::find(f4.begin(), f4.end(), nearfrequency);
//				index1 = std::distance(f4.begin(), it);
//
//				if ((nearfrequency <= required))
//				{
//					if (index1 == 26) { index2 = index1 - 1; }
//					else { index2 = index1 + 1; }
//				}
//
//				if ((nearfrequency >= required))
//				{
//					if (index1 == 0) { index2 = index1 + 1; }
//					else { index2 = index1 - 1; }
//				}
//
//
//
//#if 0
//				cout << index1 << " " << index2 << endl;
//#endif // 0
//
//
//				for (size_t j = 0; j < 10; j++)
//				{
//					coficentesreq.at(j) = (((frequenciesreq[i] - table4.at(index1).at(0))*table4.at(index2).at(j + 1)) + ((table4.at(index2).at(0) - frequenciesreq[i])*table4.at(index1).at(j + 1))) / (table4.at(index2).at(0) - table4.at(index1).at(0));
//				}
//				results.at(i) = coficentesreq;
//
//			}
//		}
//
//#if 1
//
//		coefT4.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/silva02/resultsT4.dat");
//		coefT4 << "#Coefficents for differents frequencies (Table4 - Silva et al., 2002)" << endl;
//		coefT4 << setw(WIDTH) << "c1" << setw(WIDTH) << "c2" << setw(WIDTH) << "c4" << setw(WIDTH) << "c5" << setw(WIDTH) << "c6" << setw(WIDTH) << "c7" << setw(WIDTH) << "c8" << setw(WIDTH) << "c10" << setw(WIDTH) << "P. sigma" << setw(WIDTH) << "T. sigma" << endl;
//
//		for (size_t i = 0; i < periodsreq.size(); i++)
//		{
//
//			for (size_t j = 0; j < 10; j++)
//				coefT4 << setw(WIDTH) << setprecision(6) << results.at(i).at(j);
//			coefT4 << endl;
//		}
//
//
//#endif // 0
//
//		cout << "NICE" << endl;
//		break;
//	case'C':
//		silva2002.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/silva02/table5.txt");
//		if (silva2002.good())
//		{
//			for (size_t i = 0; i < 27; i++) {
//				for (size_t j = 0; j < 11; j++)	silva2002 >> coef.at(j);
//				table5.at(i) = coef;
//			}
//		}
//
//
//		for (size_t i = 0; i < 27; i++)
//			f5.at(i) = table5.at(i).at(0);
//
//#if 0
//		cout << "Table 5 (Silva et al., 2002)" << endl;
//		for (size_t i = 0; i < 27; i++) {
//			for (size_t j = 0; j < 11; j++)	cout << setw(WIDTH) << table5.at(i).at(j);
//			cout << endl;
//		}
//		cout << endl;
//
//		for (size_t i = 0; i < frequenciesreq.size(); i++)
//			cout << frequenciesreq.at(i) << endl;
//
//
//#endif // 0
//
//		for (size_t i = 0; i < frequenciesreq.size(); i++)
//		{
//			required = frequenciesreq[i];
//			if (required > PGA) {
//				for (size_t j = 0; j < 10; j++)	coficentesreq.at(j) = table5.at(26).at(j + 1);
//				results.at(i) = coficentesreq;
//				periodsreq[i] = 0.0;
//			}
//			else {
//				nearfrequency = closest(f5, required);
//				std::vector<double>::iterator it = std::find(f5.begin(), f5.end(), nearfrequency);
//				index1 = std::distance(f5.begin(), it);
//
//				if ((nearfrequency <= required))
//				{
//					if (index1 == 26) { index2 = index1 - 1; }
//					else { index2 = index1 + 1; }
//				}
//
//				if ((nearfrequency >= required))
//				{
//					if (index1 == 0) { index2 = index1 + 1; }
//					else { index2 = index1 - 1; }
//				}
//
//
//
//#if 0
//				cout << index1 << " " << index2 << endl;
//#endif // 0
//
//
//				for (size_t j = 0; j < 10; j++)
//				{
//					coficentesreq.at(j) = (((frequenciesreq[i] - table5.at(index1).at(0))*table5.at(index2).at(j + 1)) + ((table5.at(index2).at(0) - frequenciesreq[i])*table5.at(index1).at(j + 1))) / (table5.at(index2).at(0) - table5.at(index1).at(0));
//				}
//				results.at(i) = coficentesreq;
//
//			}
//		}
//
//#if 1
//
//		coefT5.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/silva02/resultsT5.dat");
//		coefT5 << "#Coefficents for differents frequencies (Table5 - Silva et al., 2002)" << endl;
//		coefT5 << setw(WIDTH) << "c1" << setw(WIDTH) << "c2" << setw(WIDTH) << "c4" << setw(WIDTH) << "c5" << setw(WIDTH) << "c6" << setw(WIDTH) << "c7" << setw(WIDTH) << "c8" << setw(WIDTH) << "c10" << setw(WIDTH) << "P. sigma" << setw(WIDTH) << "T. sigma" << endl;
//
//		for (size_t i = 0; i < periodsreq.size(); i++)
//		{
//
//			for (size_t j = 0; j < 10; j++)
//				coefT5 << setw(WIDTH) << setprecision(6) << results.at(i).at(j);
//			coefT5 << endl;
//		}
//
//
//#endif // 0
//
//		cout << "NICE" << endl;
//			break;
//	case 'D':
//		silva2002.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/silva02/table6.txt");
//
//		if (silva2002.good())
//		{
//			for (size_t i = 0; i < 27; i++) {
//				for (size_t j = 0; j < 11; j++)	silva2002 >> coef.at(j);
//				table6.at(i) = coef;
//			}
//		}
//
//
//		for (size_t i = 0; i < 27; i++)
//			f6.at(i) = table6.at(i).at(0);
//
//#if 0
//		cout << "Table 6 (Silva et al., 2002)" << endl;
//		for (size_t i = 0; i < 27; i++) {
//			for (size_t j = 0; j < 11; j++)	cout << setw(WIDTH) << table6.at(i).at(j);
//			cout << endl;
//		}
//		cout << endl;
//
//		for (size_t i = 0; i < frequenciesreq.size(); i++)
//			cout << frequenciesreq.at(i) << endl;
//
//
//#endif // 0
//
//		for (size_t i = 0; i < frequenciesreq.size(); i++)
//		{
//			required = frequenciesreq[i];
//			if (required > PGA) {
//				for (size_t j = 0; j < 10; j++)	coficentesreq.at(j) = table6.at(26).at(j + 1);
//				results.at(i) = coficentesreq;
//				periodsreq[i] = 0.0;
//			}
//			else {
//				nearfrequency = closest(f6, required);
//				std::vector<double>::iterator it = std::find(f6.begin(), f6.end(), nearfrequency);
//				index1 = std::distance(f6.begin(), it);
//
//				if ((nearfrequency <= required))
//				{
//					if (index1 == 26) { index2 = index1 - 1; }
//					else { index2 = index1 + 1; }
//				}
//
//				if ((nearfrequency >= required))
//				{
//					if (index1 == 0) { index2 = index1 + 1; }
//					else { index2 = index1 - 1; }
//				}
//
//
//
//#if 0
//				cout << index1 << " " << index2 << endl;
//#endif // 0
//
//
//				for (size_t j = 0; j < 10; j++)
//				{
//					coficentesreq.at(j) = (((frequenciesreq[i] - table6.at(index1).at(0))*table6.at(index2).at(j + 1)) + ((table6.at(index2).at(0) - frequenciesreq[i])*table6.at(index1).at(j + 1))) / (table6.at(index2).at(0) - table6.at(index1).at(0));
//				}
//				results.at(i) = coficentesreq;
//			}
//		}
//
//#if 1
//
//		coefT6.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/silva02/resultsT6.dat");
//		coefT6 << "#Coefficents for differents frequencies (Table6 - Silva et al., 2002)" << endl;
//		coefT6 << setw(WIDTH) << "c1" << setw(WIDTH) << "c2" << setw(WIDTH) << "c4" << setw(WIDTH) << "c5" << setw(WIDTH) << "c6" << setw(WIDTH) << "c7" << setw(WIDTH) << "c8" << setw(WIDTH) << "c10" << setw(WIDTH) << "P. sigma" << setw(WIDTH) << "T. sigma" << endl;
//
//		for (size_t i = 0; i < periodsreq.size(); i++)
//		{
//
//			for (size_t j = 0; j < 10; j++)
//				coefT6 << setw(WIDTH) << setprecision(6) << results.at(i).at(j);
//			coefT6 << endl;
//		}
//
//
//#endif // 0
//
//		cout << "NICE" << endl;
//
//		break;
//
//	case 'E':
//		silva2002.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/silva02/table7.txt");
//
//		if (silva2002.good())
//		{
//			for (size_t i = 0; i < 27; i++) {
//				for (size_t j = 0; j < 11; j++)	silva2002 >> coef.at(j);
//				table7.at(i) = coef;
//			}
//		}
//
//
//		for (size_t i = 0; i < 27; i++)
//			f7.at(i) = table7.at(i).at(0);
//
//#if 0
//		cout << "Table 7 (Silva et al., 2002)" << endl;
//		for (size_t i = 0; i < 27; i++) {
//			for (size_t j = 0; j < 11; j++)	cout << setw(WIDTH) << table7.at(i).at(j);
//			cout << endl;
//		}
//		cout << endl;
//
//		for (size_t i = 0; i < frequenciesreq.size(); i++)
//			cout << frequenciesreq.at(i) << endl;
//
//
//#endif // 0
//
//		for (size_t i = 0; i < frequenciesreq.size(); i++)
//		{
//			required = frequenciesreq[i];
//			if (required > PGA) {
//				for (size_t j = 0; j < 10; j++)	coficentesreq.at(j) = table7.at(26).at(j + 1);
//				results.at(i) = coficentesreq;
//				periodsreq[i] = 0.0;
//			}
//			else {
//				nearfrequency = closest(f7, required);
//				std::vector<double>::iterator it = std::find(f7.begin(), f7.end(), nearfrequency);
//				index1 = std::distance(f7.begin(), it);
//
//				if ((nearfrequency <= required))
//				{
//					if (index1 == 26) { index2 = index1 - 1; }
//					else { index2 = index1 + 1; }
//				}
//
//				if ((nearfrequency >= required))
//				{
//					if (index1 == 0) { index2 = index1 + 1; }
//					else { index2 = index1 - 1; }
//				}
//
//
//
//#if 0
//				cout << index1 << " " << index2 << endl;
//#endif // 0
//
//
//				for (size_t j = 0; j < 10; j++)
//				{
//					coficentesreq.at(j) = (((frequenciesreq[i] - table7.at(index1).at(0))*table7.at(index2).at(j + 1)) + ((table7.at(index2).at(0) - frequenciesreq[i])*table7.at(index1).at(j + 1))) / (table7.at(index2).at(0) - table7.at(index1).at(0));
//				}
//				results.at(i) = coficentesreq;
//			}
//		}
//
//#if 1
//
//		coefT7.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/silva02/resultsT7.dat");
//		coefT7 << "#Coefficents for differents frequencies (Table7 - Silva et al., 2002)" << endl;
//		coefT7 << setw(WIDTH) << "c1" << setw(WIDTH) << "c2" << setw(WIDTH) << "c4" << setw(WIDTH) << "c5" << setw(WIDTH) << "c6" << setw(WIDTH) << "c7" << setw(WIDTH) << "c8" << setw(WIDTH) << "c10" << setw(WIDTH) << "P. sigma" << setw(WIDTH) << "T. sigma" << endl;
//
//		for (size_t i = 0; i < periodsreq.size(); i++)
//		{
//
//			for (size_t j = 0; j < 10; j++)
//				coefT7 << setw(WIDTH) << setprecision(6) << results.at(i).at(j);
//			coefT7 << endl;
//		}
//
//
//#endif // 0
//
//		cout << "NICE" << endl;
//
//		break;
//	case 'F':
//		silva2002.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/silva02/tableBC.txt");
//
//		if (silva2002.good())
//		{
//			for (size_t i = 0; i < 8; i++) {
//				for (size_t j = 0; j < 8; j++)	silva2002 >> coef.at(j);
//				tableBC.at(i) = coef;
//			}
//		}
//
//
//		for (size_t i = 0; i < 8; i++)
//			fBC.at(i) = tableBC.at(i).at(0);
//
//#if 0
//		cout << "Table BC (Silva et al., 2002)" << endl;
//		for (size_t i = 0; i < 8; i++) {
//			for (size_t j = 0; j < 8; j++)	cout << setw(WIDTH) << tableBC.at(i).at(j);
//			cout << endl;
//		}
//		cout << endl;
//
//		for (size_t i = 0; i < frequenciesreq.size(); i++)
//			cout << frequenciesreq.at(i) << endl;
//
//
//#endif // 0
//
//		for (size_t i = 0; i < frequenciesreq.size(); i++)
//		{
//			required = frequenciesreq[i];
//			if (required > PGA) {
//				for (size_t j = 0; j < 7; j++)	coficentesreq.at(j) = tableBC.at(7).at(j + 1);
//				results.at(i) = coficentesreq;
//				periodsreq[i] = 0.0;
//			}
//			else {
//				nearfrequency = closest(fBC, required);
//				std::vector<double>::iterator it = std::find(fBC.begin(), fBC.end(), nearfrequency);
//					index1 = std::distance(fBC.begin(), it);
//
//				if ((nearfrequency <= required))
//				{
//					if (index1 == 7) { index2 = index1 - 1; }
//					else { index2 = index1 + 1; }
//				}
//
//				if ((nearfrequency >= required))
//				{
//					if (index1 == 0) { index2 = index1 + 1; }
//					else { index2 = index1 - 1; }
//				}
//
//
//
//#if 0
//				cout << index1 << " " << index2 << endl;
//#endif // 0
//
//
//				for (size_t j = 0; j < 7; j++)
//				{
//					coficentesreq.at(j) = (((frequenciesreq[i] - tableBC.at(index1).at(0))*tableBC.at(index2).at(j + 1)) + ((tableBC.at(index2).at(0) - frequenciesreq[i])*tableBC.at(index1).at(j + 1))) / (tableBC.at(index2).at(0) - tableBC.at(index1).at(0));
//				}
//				results.at(i) = coficentesreq;
//			}
//		}
//		
//
//#if 1
//
//		coefTBC.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/silva02/resultsTBC.dat");
//		coefTBC << "#Coefficents for differents frequencies (TableBC - Silva et al., 2002)" << endl;
//		coefTBC << setw(WIDTH) << "c1" << setw(WIDTH) << "c2" << setw(WIDTH) << "c4" << setw(WIDTH) << "c5" << setw(WIDTH) << "c6" << setw(WIDTH) << "c7" << setw(WIDTH) << "c8" << setw(WIDTH) << "c10" << setw(WIDTH) << "P. sigma" << setw(WIDTH) << "T. sigma" << endl;
//
//		for (size_t i = 0; i < periodsreq.size(); i++)
//		{
//
//			for (size_t j = 0; j < 7; j++)
//				coefTBC << setw(WIDTH) << setprecision(6) << results.at(i).at(j);
//			coefTBC << endl;
//		}
//
//
//#endif // 0
//
//		cout << "NICE" << endl;
//
//		break;
//	default:
//		cout << "Choose another one" << endl;
//		break;
//
//	}
//
//	
//	//--------------------------------------- FINISH -------------------------------------------------------
//
//	//-------------------------- CHECKING ----------------------------------------
//#if 1
//	switch (OPTION)
//	{
//	case 'F':
//		for (size_t i = 0; i < frequenciesreq.size(); i++)
//		{
//			for (size_t j = 0; j < 7; j++)
//				cout << setw(WIDTH) << setprecision(6) << results.at(i).at(j);
//			cout << endl;
//		}
//			break;
//
//	default:
//		for (size_t i = 0; i < frequenciesreq.size(); i++)
//		{
//			for (size_t j = 0; j < 10; j++)
//				cout << setw(WIDTH) << setprecision(6) << results.at(i).at(j);
//			cout << endl;
//		}
//		break;
//	}
//	
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
//	ofstream attenueationtableS02;
//
//	switch (OPTION)
//	{
//	case 'F':
//		attenueationtableS02.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/silva02/silva02bc.atn");
//
//		attenueationtableS02 << setprecision(PRECISION2);
//		//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Description" << setw(WIDTH) << ": Sample attenuation file constructed for illustration purposes (2008)" << endl;	
//		//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Units" << setw(WIDTH) << ": cm/sec/sec" << endl;
//		//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Distribution" << setw(WIDTH) << ": 2" << endl;
//		//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Dimension" << setw(WIDTH) << ": Aceleration" << endl;
//		attenueationtableS02 << setw(WIDTH) << MINF << setw(WIDTH) << MSUP << setw(WIDTH) << NMAG << endl;
//		attenueationtableS02 << setw(WIDTH) << RINF << setw(WIDTH) << RSUP << setw(WIDTH) << NRAD << setw(WIDTH) << TYPEMD << endl;
//
//		attenueationtableS02 << setprecision(PRECISION);
//
//		for (size_t i = 0; i < periodsreq.size(); i++)//Loop over periods
//		{
//			sigma = results.at(i).at(9);
//			attenueationtableS02 << setw(WIDTH) << periodsreq.at(i) << setw(WIDTH) << sigma << setw(WIDTH) << AMAX << endl;
//			for (size_t j = 0; j < NMAG; j++)//Loop over magnitudes
//			{
//
//
//				for (size_t k = 0; k < NRAD; k++)//Loop over coeficents
//				{
//					R = distances[k];
//					A = results.at(i).at(0);
//					B = results.at(i).at(1)*magnitudes[j];
//					C = results.at(i).at(3) + results.at(i).at(4)*magnitudes[j];
//					D = log(R + pow(e, results.at(i).at(2)));
//					E = (results.at(i).at(5))*pow((magnitudes[j] - 6), 2.0);
//					aceleraciones[k] = pow(e, A + B + C*D + E);
//
//					attenueationtableS02 << setw(WIDTH) << aceleraciones[k];//Saving values
//				}
//				attenueationtableS02 << endl;
//			}
//		}
//
//		break;
//	default:
//		attenueationtableS02.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/silva02/silva02hr.atn");
//
//		attenueationtableS02 << setprecision(PRECISION2);
//		//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Description" << setw(WIDTH) << ": Sample attenuation file constructed for illustration purposes (2008)" << endl;	
//		//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Units" << setw(WIDTH) << ": cm/sec/sec" << endl;
//		//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Distribution" << setw(WIDTH) << ": 2" << endl;
//		//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Dimension" << setw(WIDTH) << ": Aceleration" << endl;
//		attenueationtableS02 << setw(WIDTH) << MINF << setw(WIDTH) << MSUP << setw(WIDTH) << NMAG << endl;
//		attenueationtableS02 << setw(WIDTH) << RINF << setw(WIDTH) << RSUP << setw(WIDTH) << NRAD << setw(WIDTH) << TYPEMD << endl;
//
//		attenueationtableS02 << setprecision(PRECISION);
//
//		for (size_t i = 0; i < periodsreq.size(); i++)//Loop over periods
//		{
//			sigma = results.at(i).at(9);
//			attenueationtableS02 << setw(WIDTH) << periodsreq.at(i) << setw(WIDTH) << sigma << setw(WIDTH) << AMAX << endl;
//			for (size_t j = 0; j < NMAG; j++)//Loop over magnitudes
//			{
//
//
//				for (size_t k = 0; k < NRAD; k++)//Loop over coeficents
//				{
//					R = distances[k];
//					A = results.at(i).at(0);
//					B = results.at(i).at(1)*magnitudes[j];
//					C = results.at(i).at(4) + results.at(i).at(5)*magnitudes[j];
//					D = log(R + pow(e, results.at(i).at(2)));
//					E = (results.at(i).at(7))*pow((magnitudes[j] - 6), 2.0);
//					aceleraciones[k] = pow(e, A + B + C*D + E);
//
//					attenueationtableS02 << setw(WIDTH) << aceleraciones[k];//Saving values
//				}
//				attenueationtableS02 << endl;
//			}
//		}
//		break;
//	}
//
//	
//	system("pause");
//	return 0;
//}