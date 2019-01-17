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
//
//#ifndef WIDTH	
//#define WIDTH 13
//#endif // !WIDTH
//#ifndef PRECISION
//#define PRECISION 6
//#endif // !PRECISION
//
//
//using namespace std;
//double closest(std::vector<double> const& vec, double value);
//
//int main() {
//	typedef std::vector<vector<double>> Matrix;
//	typedef std::vector<double> Vector;
//
//	//TABLE 6
//
//	Matrix resultadosT6(25);//vertical size of Matrix
//	Matrix table6(26);//vertical size of Matrix
//	Vector f6(26);
//	// Open our file tabla6.txt
//	ifstream inFile1;
//	inFile1.open("table6.txt");
//	
//	vector <double> coef1(12);
//	vector <double> required6(25); 
//	vector <double> period6(25);
//	vector <double> coefreq6(12);
//
//	// If we can read/write great
//	if (inFile1.good())
//	{
//		for (size_t i = 0; i < 26; i++) {
//			for (size_t j = 0; j < 12; j++)	inFile1 >> coef1.at(j);
//			table6.at(i) = coef1;
//		}
//	}
//
//	for (size_t i = 0; i < 26; i++) 
//		f6.at(i) = table6.at(i).at(0);
//
//	
//		
//		// Plot table, only if is necessary 1 = show, 2 = hide
//#if 0
//	cout << "Tabla 6:" << endl;
//	for (size_t i = 0; i < 26; i++) {
//		for (size_t j = 0; j < 12; j++)	cout << setw(WIDTH) << table6.at(i).at(j) << " ";
//		cout << endl;
//	}
//
//#endif // 0
//	ifstream inFreq6; //h
//	inFreq6.open("freq6.txt"); //h
//
//	if (inFreq6.good())
//	{
//		for (size_t i = 0; i < 25; i++)
//			inFreq6 >> required6.at(i); 
//	}
//#if 0
//	cout << endl;
//	for (size_t i = 0; i < 25; i++)
//		cout << required6.at(i) << endl;
//#endif // 1
//
//	for (size_t i = 0; i < 25; i++)
//		period6.at(i) = 1/(required6.at(i));
//	
//	for (size_t k = 0; k < 25; k++)
//	{
//		double frequency1 = 0;
//
//		double required1 = required6[k]; 
//		frequency1 = closest(f6, required1);
//		//cout << frequency1 << endl;
//
//		std::vector<double>::iterator it1 = std::find(f6.begin(), f6.end(), frequency1);
//		int index1t6 = std::distance(f6.begin(), it1);
//		int index2t6;
//
//		if (frequency1 < required1) { index2t6 = index1t6 + 1; }
//		else { index2t6 = index1t6 - 1; }
//
//
//#if 0
//		cout << index1t6 << " " << index2t6 << endl;
//#endif // 0
//		coefreq6[0] = required6[k];
//		coefreq6[1] = period6[k];
//
//		for (size_t j = 0; j < 10; j++)
//		{
//			coefreq6.at(j+2) = (((required6.at(k) - table6.at(index1t6).at(0))*table6.at(index2t6).at(j + 2)) +
//				((table6.at(index2t6).at(0) - required6.at(k))*table6.at(index1t6).at(j + 2))) / (table6.at(index2t6).at(0) - table6.at(index1t6).at(0));
//		}
//		resultadosT6.at(k) = coefreq6;
//	}
//
//#if 1
//	ofstream coefT6;
//	coefT6.open("resultsT6.dat");
//	coefT6 << "#Coefficents for differents frequencies (Table6 - Atkinson and Boore, 2006)" << endl;
//	coefT6 << setw(WIDTH) << "f(Hz)" << setw(WIDTH) <<"T(sec)"<< setw(WIDTH) <<"c1" << setw(WIDTH) <<"c2" << setw(WIDTH) <<"c3" << setw(WIDTH) <<"c4" << setw(WIDTH) <<"c5" << setw(WIDTH) <<"c6" << setw(WIDTH) <<"c7" << setw(WIDTH) <<"c8" << setw(WIDTH) <<"c9" << setw(WIDTH) <<"c10"<< endl;
//	
//	for (size_t i = 0; i < 25; i++)
//	{
//	
//		for (size_t j = 0; j < 12; j++)	
//			coefT6 << setw(WIDTH) << setprecision(4) << resultadosT6.at(i).at(j);
//		coefT6 << endl;
//	}
//	
//
//#endif // 0
//
//	//TABLE 7
//
//	Matrix resultadosT7(25);//vertical size of Matrix
//	Matrix table7(26);//vertical size of Matrix
//	Vector f7(26);
//	// Open our file tabla7.txt
//	ifstream inFile2;
//	inFile2.open("table7.txt");
//
//	vector <double> coef2(4);
//	vector <double> required7(25);
//	vector <double> period7(25);
//	vector <double> coefreq7(5);
//
//	// If we can read/write great
//	if (inFile2.good())
//	{
//		for (size_t i = 0; i < 26; i++) {
//			for (size_t j = 0; j < 4; j++)	inFile2 >> coef2.at(j);
//			table7.at(i) = coef2;
//		}
//	}
//
//	for (size_t i = 0; i < 26; i++)
//		f7.at(i) = table7.at(i).at(0);
//
//
//
//	// Plot table, only if is necessary 1 = show, 2 = hide
//#if 0
//	cout << "Tabla 7:" << endl;
//	for (size_t i = 0; i < 26; i++) {
//		for (size_t j = 0; j < 4; j++)	cout << setw(WIDTH) << table7.at(i).at(j) << " ";
//		cout << endl;
//	}
//
//#endif // 0
//	ifstream inFreq7; //h
//	inFreq7.open("freq7.txt"); //h
//
//	if (inFreq7.good())
//	{
//		for (size_t i = 0; i < 25; i++)
//			inFreq7 >> required7.at(i);
//	}
//#if 0
//	cout << endl;
//	for (size_t i = 0; i < 25; i++)
//		cout << required7.at(i) << endl;
//#endif // 1
//
//	for (size_t i = 0; i < 25; i++)
//		period7.at(i) = 1 / (required7.at(i));
//
//	for (size_t k = 0; k < 25; k++)
//	{
//		double frequency2 = 0;
//
//		double required2 = required7[k];
//		frequency2 = closest(f7, required2);
//		//cout << frequency2 << endl;
//
//		std::vector<double>::iterator it2 = std::find(f7.begin(), f7.end(), frequency2);
//		int index1t7 = std::distance(f7.begin(), it2);
//		int index2t7;
//
//		if (frequency2 < required2) { index2t7 = index1t7 + 1; }
//		else { index2t7 = index1t7 - 1; }
//
//
//#if 0
//		cout << index1t7 << " " << index2t7 << endl;
//#endif // 0
//		coefreq7[0] = required7[k];
//		coefreq7[1] = period7[k];
//
//		for (size_t j = 0; j < 3; j++)
//		{
//			coefreq7.at(j + 2) = (((required7.at(k) - table7.at(index1t7).at(0))*table7.at(index2t7).at(j + 1)) +
//				((table7.at(index2t7).at(0) - required7.at(k))*table7.at(index1t7).at(j + 1))) / (table7.at(index2t7).at(0) - table7.at(index1t7).at(0));
//		}
//		resultadosT7.at(k) = coefreq7;
//	}
//
//#if 1
//	ofstream coefT7;
//	coefT7.open("resultsT7.dat");
//	coefT7 << "#Coefficents for differents frequencies (Table7 - Atkinson and Boore, 2006)" << endl;
//	coefT7 << setw(WIDTH) << "f(Hz)" << setw(WIDTH) << "T(sec)" << setw(WIDTH) << "delta" << setw(WIDTH) << "Ml" << setw(WIDTH) << "Mh"<< endl;
//
//	for (size_t i = 0; i < 25; i++)
//	{
//			for (size_t j = 0; j < 5; j++)
//			coefT7 << setw(WIDTH) << setprecision(4) << resultadosT7.at(i).at(j);
//		coefT7 << endl;
//	}
//
//
//#endif // 0
//	
//	//TABLE 8
//
//	Matrix resultadosT8(25);//vertical size of Matrix
//	Matrix table8(26);//vertical size of Matrix
//	Vector f8(26);
//	// Open our file tabla8.txt
//	ifstream inFile3;
//	inFile3.open("table8.txt");
//
//	vector <double> coef3(4);
//	vector <double> required8(25);
//	vector <double> period8(25);
//	vector <double> coefreq8(5);
//
//	// If we can read/write great
//	if (inFile3.good())
//	{
//		for (size_t i = 0; i < 26; i++) {
//			for (size_t j = 0; j < 4; j++)	inFile3 >> coef3.at(j);
//			table8.at(i) = coef3;
//		}
//	}
//
//	for (size_t i = 0; i < 26; i++)
//		f8.at(i) = table8.at(i).at(0);
//
//
//
//	// Plot table, only if is necessary 1 = show, 2 = hide
//#if 0
//	cout << "Tabla 8:" << endl;
//	for (size_t i = 0; i < 26; i++) {
//		for (size_t j = 0; j < 4; j++)	cout << setw(WIDTH) << table8.at(i).at(j) << " ";
//		cout << endl;
//	}
//
//#endif // 0
//	ifstream inFreq8; //h
//	inFreq8.open("freq8.txt"); //h
//
//	if (inFreq8.good())
//	{
//		for (size_t i = 0; i < 25; i++)
//			inFreq8 >> required8.at(i);
//	}
//#if 0
//	cout << endl;
//	for (size_t i = 0; i < 25; i++)
//		cout << required8.at(i) << endl;
//#endif // 1
//
//	for (size_t i = 0; i < 25; i++)
//		period8.at(i) = 1 / (required8.at(i));
//
//	for (size_t k = 0; k < 25; k++)
//	{
//		double frequency3 = 0;
//
//		double required3 = required8[k];
//		frequency3 = closest(f8, required3);
//		//cout << frequency3 << endl;
//
//		std::vector<double>::iterator it3 = std::find(f8.begin(), f8.end(), frequency3);
//		int index1t8 = std::distance(f8.begin(), it3);
//		int index2t8;
//
//		if (frequency3 < required3) { index2t8 = index1t8 + 1; }
//		else { index2t8 = index1t8 - 1; }
//
//
//#if 0
//		cout << index1t8 << " " << index2t8 << endl;
//#endif // 0
//		coefreq8[0] = required8[k];
//		coefreq8[1] = period8[k];
//
//		for (size_t j = 0; j < 3; j++)
//		{
//			coefreq8.at(j + 2) = (((required8.at(k) - table8.at(index1t8).at(0))*table8.at(index2t8).at(j + 1)) +
//				((table8.at(index2t8).at(0) - required8.at(k))*table8.at(index1t8).at(j + 1))) / (table8.at(index2t8).at(0) - table8.at(index1t8).at(0));
//		}
//		resultadosT8.at(k) = coefreq8;
//	}
//
//#if 1
//	ofstream coefT8;
//	coefT8.open("resultsT8.dat");
//	coefT8 << "#Coefficents for differents frequencies (Table8 - Atkinson and Boore, 2006)" << endl;
//	coefT8 << setw(WIDTH) << "f(Hz)" << setw(WIDTH) << "T(sec)" << setw(WIDTH) << "bln" << setw(WIDTH) << "b1" << setw(WIDTH) << "b2" << endl;
//
//	for (size_t i = 0; i < 25; i++)
//	{
//		for (size_t j = 0; j < 5; j++)
//			coefT8 << setw(WIDTH) << setprecision(4) << resultadosT8.at(i).at(j);
//		coefT8 << endl;
//	}
//
//
//#endif // 0
//
//	//TABLE 8
//
//	Matrix resultadosT9(25);//vertical size of Matrix
//	Matrix table9(26);//vertical size of Matrix
//	Vector f9(26);
//	// Open our file tabla9.txt
//	ifstream inFile4;
//	inFile4.open("table9.txt");
//
//	vector <double> coef4(12);
//	vector <double> required9(25);
//	vector <double> period9(25);
//	vector <double> coefreq9(12);
//
//	// If we can read/write great
//	if (inFile4.good())
//	{
//		for (size_t i = 0; i < 26; i++) {
//			for (size_t j = 0; j < 12; j++)	inFile4 >> coef4.at(j);
//			table9.at(i) = coef4;
//		}
//	}
//
//	for (size_t i = 0; i < 26; i++)
//		f9.at(i) = table9.at(i).at(0);
//
//
//
//	// Plot table, only if is necessary 1 = show, 2 = hide
//#if 0
//	cout << "Tabla 9:" << endl;
//	for (size_t i = 0; i < 26; i++) {
//		for (size_t j = 0; j < 12; j++)	cout << setw(WIDTH) << table9.at(i).at(j) << " ";
//		cout << endl;
//	}
//
//#endif // 0
//	ifstream inFreq9; //h
//	inFreq9.open("freq9.txt"); //h
//
//	if (inFreq9.good())
//	{
//		for (size_t i = 0; i < 25; i++)
//			inFreq9 >> required9.at(i);
//	}
//#if 0
//	cout << endl;
//	for (size_t i = 0; i < 25; i++)
//		cout << required9.at(i) << endl;
//#endif // 1
//
//	for (size_t i = 0; i < 25; i++)
//		period9.at(i) = 1 / (required9.at(i));
//
//	for (size_t k = 0; k < 25; k++)
//	{
//		double frequency4 = 0;
//
//		double required4 = required9[k];
//		frequency4 = closest(f9, required4);
//		//cout << frequency4 << endl;
//
//		std::vector<double>::iterator it4 = std::find(f9.begin(), f9.end(), frequency4);
//		int index1t9 = std::distance(f9.begin(), it4);
//		int index2t9;
//
//		if (frequency4 < required4) { index2t9 = index1t9 + 1; }
//		else { index2t9 = index1t9 - 1; }
//
//
//#if 0
//		cout << index1t9 << " " << index2t9 << endl;
//#endif // 0
//		coefreq9[0] = required9[k];
//		coefreq9[1] = period9[k];
//
//		for (size_t j = 0; j < 10; j++)
//		{
//			coefreq9.at(j + 2) = (((required9.at(k) - table9.at(index1t9).at(0))*table9.at(index2t9).at(j + 2)) +
//				((table9.at(index2t9).at(0) - required9.at(k))*table9.at(index1t9).at(j + 2))) / (table9.at(index2t9).at(0) - table9.at(index1t9).at(0));
//		}
//		resultadosT9.at(k) = coefreq9;
//	}
//
//#if 1
//	ofstream coefT9;
//	coefT9.open("resultsT9.dat");
//	coefT9 << "#Coefficents for differents frequencies (Table9 - Atkinson and Boore, 2006)" << endl;
//	coefT9 << setw(WIDTH) << "f(Hz)" << setw(WIDTH) << "T(sec)" << setw(WIDTH) << "c1" << setw(WIDTH) << "c2" << setw(WIDTH) << "c3" << setw(WIDTH) << "c4" << setw(WIDTH) << "c5" << setw(WIDTH) << "c6" << setw(WIDTH) << "c7" << setw(WIDTH) << "c8" << setw(WIDTH) << "c9" << setw(WIDTH) << "c10" << endl;
//
//	for (size_t i = 0; i < 25; i++)
//	{
//
//		for (size_t j = 0; j < 12; j++)
//			coefT9 << setw(WIDTH) << setprecision(4) << resultadosT9.at(i).at(j);
//		coefT9 << endl;
//	}
//
//
//#endif // 0
//
//
//	system("pause");
//	return 0;
//}
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//double closest(std::vector<double> const& vec, double value) {
//	auto const it = std::lower_bound(vec.begin(), vec.end(), value);
//	if (it == vec.end()) { return -1; }
//
//	return *it;
//}