#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <map>
#include <cstdlib>
#include <iterator>
#include <algorithm>
#include <cmath>

#ifndef WIDTH	
#define WIDTH 13
#endif // !WIDTH
#ifndef PRECISION
#define PRECISION 6
#endif // !PRECISION


using namespace std;

double closest(std::vector<double> const& vec, double value) {
	auto const it = std::lower_bound(vec.begin(), vec.end(), value);
	if (it == vec.end()) { return -1; }
	return *it;
}

int main() {

	typedef std::vector<vector<double>> Matrix;
	typedef std::vector<double> Vector;
	double MINF; //Lower limit og magnitude given in the tabl
	double MSUP; //Upper limit of magnitude given in the table
	double NMAG; //Number of magnitudes for which intensity is given
	double DMAG;
	double RINF; //Lower limit of distance given in the table
	double RSUP; //Upper limit of distance given in the table
	double NRAD; //Number of distances for which intensity is given
	double TYPE; //An integer indicating the type of distance used by the attenuation table
	double DLRAD;
	double NT; //Number of periods
	double S;
	double sigma = 0.7;
	double AMAX = 0; //Type of stadistic distribution

	int index1, index2;

	//--- ATTENUATION VARIABLES ----

	double M, f0, f1, f2, R0, Rcd;

	S = 0.0;
	Rcd = 4.0;
	R0 = 3.0;

	Vector periodsreq;
	Vector frequenciesreq;
	double aux1, aux2;

	//-------------------------------
	
	ifstream periods;
	periods.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/p6.txt");

	//If we can read/write great
	while (periods.good())
	{
		periods >> aux1;
		periodsreq.push_back(aux1);
	}

#if 0
	for (size_t i = 0; i < periodsreq.size(); i++)
	{
		cout << periodsreq[i] << endl;
	}

#endif // 0


	for (size_t i = 0; i < periodsreq.size(); i++) // Changing periods to frequencies
	{
		aux2 = 1 / periodsreq.at(i);
		frequenciesreq.push_back(aux2);
	}


	// ---------------------------------- ASKING FOR DATA ------------------------------------------------
	::std::cout << "Enter the minimum magnitude, maximum and the number of intermediate quantities: " << endl;
	/*cin >> MINF >> MSUP >> NMAG;*/
	MINF = 4.0;
	MSUP = 8.5;
	NMAG = 12.0;
	//cout << "Enter the minimum distance, maximum and the number of intermediate distances: " << endl;
	//cin >> RINF >> RSUP >> NRAD;
	RINF = 70.0;
	RSUP = 140.0;
	NRAD = 20.0;

	// -------------------------------------- FINISH ------------------------------------------------------


	//-------------------------------- GENRATING COEFICENTS -----------------------------------------------
	if (S == 0.0) {
			Matrix results(periodsreq.size());//vertical size of Matrix
			Matrix table6(26);//vertical size of Matrix
			
			// Open our file tabla6.txt
			ifstream inFile1;
			inFile1.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/table6.txt");


			vector <double> f6(26);// Vector of table 6's frequencies
			vector <double> coef(12);// Vector of table 6's coeficents 
			vector <double> coefreq(10);
			// If we can read/write great
			if (inFile1.good())
			{
				for (size_t i = 0; i < 26; i++) {
					for (size_t j = 0; j < 12; j++)	inFile1 >> coef.at(j);
					table6.at(i) = coef;
				}
			}
		
			for (size_t i = 0; i < 26; i++) 
				f6.at(i) = table6.at(i).at(0);
		
			
				
				// Plot table, only if is necessary 1 = show, 2 = hide
		#if 0
			cout << "Tabla 6:" << endl;
			for (size_t i = 0; i < 26; i++) {
				for (size_t j = 0; j < 12; j++)	cout << setw(WIDTH) << table6.at(i).at(j) << " ";
				cout << endl;
			}		


			cout << endl;
			for (size_t i = 0; i < frequenciesreq.size(); i++)
				cout << frequenciesreq.at(i) << endl;
		#endif // 1
					
			for (size_t k = 0; k < frequenciesreq.size(); k++)
			{
				double nearfrequency = 0;
		
				double required = frequenciesreq[k];
				nearfrequency = closest(f6, required);
				//cout << frequency << endl;
		
				std::vector<double>::iterator it = std::find(f6.begin(), f6.end(), nearfrequency);
				index1 = std::distance(f6.begin(), it);
	
				if (nearfrequency < required) { index2 = index1 + 1; }
				else { index2 = index1 - 1; }
						
		
		#if 0
				cout << index1 << " " << index2 << endl;
		#endif // 0

				for (size_t j = 0; j < 10; j++)
				{
					coefreq.at(j) = (((frequenciesreq.at(k) - table6.at(index1).at(0))*table6.at(index2).at(j + 2)) +
						((table6.at(index2).at(0) - frequenciesreq.at(k))*table6.at(index1).at(j + 2))) / (table6.at(index2).at(0) - table6.at(index1).at(0));
				}
				results.at(k) = coefreq;
			}
		
		#if 0
			ofstream coefT6; // Archive
			coefT6.open("resultsT6.dat");
			coefT6 << "#Coefficents for differents frequencies (Table6 - Atkinson and Boore, 2006)" << endl;
			coefT6 << setw(WIDTH) << "f(Hz)" << setw(WIDTH) <<"T(sec)"<< setw(WIDTH) <<"c1" << setw(WIDTH) <<"c2" << setw(WIDTH) <<"c3" << setw(WIDTH) <<"c4" << setw(WIDTH) <<"c5" << setw(WIDTH) <<"c6" << setw(WIDTH) <<"c7" << setw(WIDTH) <<"c8" << setw(WIDTH) <<"c9" << setw(WIDTH) <<"c10"<< endl;
			
			for (size_t i = 0; i < periodsreq.size(); i++)
			{
			
				for (size_t j = 0; j < 12; j++)	
					coefT6 << setw(WIDTH) << setprecision(4) << results.at(i).at(j);
				coefT6 << endl;
			}
			
		
		#endif // 0


	}
	else {
		Matrix results(periodsreq.size());//vertical size of Matrix
		Matrix table9(26);//vertical size of Matrix

						  // Open our file tabla9.txt
		ifstream inFile2;
		inFile2.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/table9.txt");


		vector <double> f9(26);// Vector of table 9's frequencies
		vector <double> coef(12);// Vector of table 9's coeficents 
		vector <double> coefreq(10);// Vector of results
		// If we can read/write great
		if (inFile2.good())
		{
			for (size_t i = 0; i < 26; i++) {
				for (size_t j = 0; j < 12; j++)	inFile2 >> coef.at(j);
				table9.at(i) = coef;
			}
		}

		for (size_t i = 0; i < 26; i++)
			f9.at(i) = table9.at(i).at(0);



		// Plot table, only if is necessary 1 = show, 2 = hide
#if 0
		cout << "Tabla 9:" << endl;
		for (size_t i = 0; i < 26; i++) {
			for (size_t j = 0; j < 12; j++)	cout << setw(WIDTH) << table9.at(i).at(j) << " ";
			cout << endl;
		}


		cout << endl;
		for (size_t i = 0; i < frequenciesreq.size(); i++)
			cout << frequenciesreq.at(i) << endl;
#endif // 1

		for (size_t k = 0; k < frequenciesreq.size(); k++)
		{
			double nearfrequency = 0;

			double required = frequenciesreq[k];
			nearfrequency = closest(f9, required);
			//cout << frequency << endl;

			std::vector<double>::iterator it = std::find(f9.begin(), f9.end(), nearfrequency);
			int index1 = std::distance(f9.begin(), it);
			int index2;

			if (nearfrequency < required) { index2 = index1 + 1; }
			else { index2 = index1 - 1; }


#if 0
			cout << index1t6 << " " << index2t6 << endl;
#endif // 0

			for (size_t j = 0; j < 10; j++)
			{
				coefreq.at(j) = (((frequenciesreq.at(k) - table9.at(index1).at(0))*table9.at(index2).at(j + 2)) +
					((table9.at(index2).at(0) - frequenciesreq.at(k))*table9.at(index1).at(j + 2))) / (table9.at(index2).at(0) - table9.at(index1).at(0));
			}
			results.at(k) = coefreq;
		}

#if 0
		ofstream coefT9; // Archive
		coefT9.open("resultsT9.dat");
		coefT9 << "#Coefficents for differents frequencies (Table9 - Atkinson and Boore, 2006)" << endl;
		coefT9 << setw(WIDTH) << "f(Hz)" << setw(WIDTH) << "T(sec)" << setw(WIDTH) << "c1" << setw(WIDTH) << "c2" << setw(WIDTH) << "c3" << setw(WIDTH) << "c4" << setw(WIDTH) << "c5" << setw(WIDTH) << "c6" << setw(WIDTH) << "c7" << setw(WIDTH) << "c8" << setw(WIDTH) << "c9" << setw(WIDTH) << "c10" << endl;

		for (size_t i = 0; i < periodsreq.size(); i++)
		{

			for (size_t j = 0; j < 12; j++)
				coefT9 << setw(WIDTH) << setprecision(4) << results.at(i).at(j);
			coefT9 << endl;
		}

#endif // 0

	}
	
	//---------------------------- FINISH ----------------------------------


	//---------------------- GENERATING VALUES -------------------------------
	Vector distances(NRAD);
	Vector magnitudes(NMAG);
	DMAG = (MSUP - MINF) / (NMAG - 1);
	DLRAD = (log10(RSUP) - log10(RINF)) / (NRAD - 1);
	for (size_t i = 0; i < NRAD; i++)
	{
		distances[i] = pow(10.0, log(RINF) + i*DLRAD);
	}

	for (size_t i = 0; i < NMAG; i++)
	{
		magnitudes[i] = MINF + i*DMAG;
	}


	//---------------------------- FINISH -----------------------------------

	//---------------------------- OUTPUT -----------------------------------

	::std::cout << setw(WIDTH) << MINF << setw(WIDTH) << MSUP << setw(WIDTH) << NMAG << endl;
	::std::cout << setw(WIDTH) << RINF << setw(WIDTH) << RSUP << setw(WIDTH) << NRAD << endl;
	ofstream prueba;
	prueba.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/prueba.txt");
	prueba << "escgfxdgxdfges"<<endl;
	cout << "Hola" << endl;

	for (size_t i = 0; i < periodsreq.size(); i++)
	{
		f0 = max(log(R0 / Rcd), 0.0);
		f1 = min(log(Rcd), log(RINF));
		f2 = max(log(Rcd / RSUP), 0.0);


		::std::cout << setw(WIDTH) << periodsreq.at(i) << setw(WIDTH) << sigma << setw(WIDTH) << AMAX << endl;
		
	}
	
	system("pause");
	return 0;
}