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
#define WIDTH 16
#endif // !WIDTH
#ifndef PGA
#define PGA 100.0
#endif // !PGA

#ifndef PRECISION
#define PRECISION 8
#endif // !PRECISION
#ifndef PRECISION2
#define PRECISION2 3
#endif // !PRECISION2
#ifndef e
#define e 2.7182818284590
#endif // !e



using namespace std;

double closest(std::vector<double> const& vec, double value) {
	auto const it = std::lower_bound(vec.begin(), vec.end(), value);
	if (it == vec.end()) { return -1; }
	return *it;
}

int main() {

	typedef std::vector<vector<double>> Matrix;
	typedef std::vector<vector<vector<double>>> Matrix2;
	typedef std::vector<double> Vector;
	double MINF; //Lower limit og magnitude given in the tabl
	double MSUP; //Upper limit of magnitude given in the table
	int NMAG; //Number of magnitudes for which intensity is given
	double DMAG;
	double RINF; //Lower limit of distance given in the table
	double RSUP; //Upper limit of distance given in the table
	int NRAD; //Number of distances for which intensity is given
	int TYPEMD; //An integer indicating the type of distance used by the attenuation table
	double DLRAD;
	double AMAX = 0; //Type of stadistic distribution
	double STRESS;//Inicial stress
	double VS30; //Average shear-wave velocity in the upper 30 m of the geological profile, in m/s
	double BDEPTH=1;//Basin depth under the site in km

	//--- ATTENUATION VARIABLES ----
	Vector periodsreq;
	Vector frequenciesreq;
	double aux1, aux2;
	char OPTION;
	double lnG1, lnG2, lnG3, lnG4, lnG5;
	//-------------------------------

	ifstream periods;
	periods.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/data/graizer15/periods.txt");

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
#if 0
	cout << periodsreq.size() << " " << frequenciesreq.size() << endl;
	for (size_t i = 0; i < frequenciesreq.size(); i++)
	{
		cout << frequenciesreq[i] << endl;
	}
#endif // 1


	// ---------------------------------- ASKING FOR DATA ------------------------------------------------

	//Enter the minimum magnitude, maximum and the number of intermediate quantities
	MINF = 4.0;
	MSUP = 8.5;
	NMAG = 12;
	//Enter the minimum distance, maximum and the number of intermediate distances
	RINF = 1;
	RSUP = 1000.0;
	NRAD = 100;

	cout << "ENTER THE SHEAR WAVE VELOCITY IN m/s:" << endl;
	cin >> VS30;
	cout << "CHOOSE THE TYPE OF REGION:" << endl;
	cout << "-CEUS (A)" << endl;
	cout << "-Gulf Coast Region (B)" << endl;
	cin >> OPTION;

	// -------------------------------------- FINISH ------------------------------------------------------



	//-------------------------------- GENERATING COEFICENTS -----------------------------------------------
	Matrix resultsT4(periodsreq.size());//Vertical size of Matrix
	Matrix table4(96);//T97, Table 4
	Vector coefT4(4);
	Vector coficentesreq(3);
	double c1 = 0.4, c2 = -6.25, c3 = 0.55, c4 = 2.237, c5 = -7.542, c6 = -0.125, c7 = 1.19, c8 = -6.13, c9 = 0.6, c11 = 3.9, c12 = -0.3445, D2=0.7,F = 2.232, sigma = 0.848, R2, r2, Q0;
	double SAnorm, Lin_amp,I,  mu, varsigma, Tsp, S,kVs30,fVs30;
	double m1 = -0.002, m2 = -0.12, m3 = 3, a1 = 0.0347, a2 = -0.5542, a3 = 3.694, Dsp = 0.75, t1 = 0.0008, t2 = 0.16, t3 = -0.4875, s1 = 0, s2 = 0.077, s3 = 0.3251;
	int index1 = 0, index2 = 0;//Positions to interpolate

	ifstream graizer2015;
	ofstream coeficentsT4; // Archive
	Vector f4(96);
	double nearfrequency, nearfrequency2;
	double required;
	switch (OPTION)
	{
	case 'A':
		Q0 = 650.0;
		break;
	case 'B':
		Q0 = 300.0;
		break;

	}

	//----------------------- FINISH CHECKING -------------------------------------

	//---------------------- GENERATING VALUES -------------------------------

	Vector acelerations(NRAD);
	Vector distances(NRAD); //Rjb
	Vector magnitudes(NMAG);
	DMAG = (MSUP - MINF) / (NMAG - 1);
	DLRAD = (log10(RSUP) - log10(RINF)) / (NRAD - 1);


	for (size_t i = 0; i < NRAD; i++)
	{
		distances[i] = pow(10, log10(RINF) + i*DLRAD);
	}
	distances[NRAD - 1] = RSUP; //Changing last value
	for (size_t i = 0; i < NMAG; i++)
	{
		magnitudes[i] = MINF + i*DMAG;
	}
	magnitudes[NMAG - 1] = MSUP;//Changing last value

								//---------------------------- FINISH -----------------------------------
								//---------------------------- OUTPUT -----------------------------------
								//Value Type of distance
								//	1 (or blank) Focal
								//	2 Epicentral
								//	3 Joyner and Boore
								//	4 Closest to rupture area(Rrup)
	TYPEMD = 4;
#if 1
	cout << "MAGNITUDES:" << endl;
	for (size_t i = 0; i < NMAG; i++)
	{
		cout << magnitudes[i] << endl;
	}
	cout << endl;
#endif // 0
#if 1
	cout << "DISTANCES:" << endl;
	for (size_t i = 0; i < NRAD; i++)
	{
		cout << setw(WIDTH) << distances[i];
	}
	cout << endl;
#endif // 0	

	//OUTPUT 1
	ofstream attenueationtableG15;

	attenueationtableG15.open("C:/Users/Hugo Ninnanya/Documents/GibHub/aninanya/BOORE/BOORE/results/graizer15/graizer15.atn");

	attenueationtableG15 << setprecision(PRECISION2);
	//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Description" << setw(WIDTH) << ": Sample attenuation file constructed for illustration purposes (2008)" << endl;	
	//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Units" << setw(WIDTH) << ": cm/sec/sec" << endl;
	//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Distribution" << setw(WIDTH) << ": 2" << endl;
	//ab06 << setw(WIDTH) << "#" << setw(WIDTH) << ": Dimension" << setw(WIDTH) << ": Aceleration" << endl;
	attenueationtableG15 << setw(WIDTH) << MINF << setw(WIDTH) << MSUP << setw(WIDTH) << NMAG << endl;
	attenueationtableG15 << setw(WIDTH) << RINF << setw(WIDTH) << RSUP << setw(WIDTH) << NRAD << setw(WIDTH) << TYPEMD << endl;

	attenueationtableG15 << setprecision(PRECISION);

	for (size_t i = 0; i < periodsreq.size(); i++)//Loop over periods
	{
		attenueationtableG15 << setw(WIDTH) << periodsreq.at(i) << setw(WIDTH) << "0.1" << setw(WIDTH) << AMAX << endl;
		for (size_t j = 0; j < NMAG; j++)//Loop over magnitudes
		{
			R2 = c4*magnitudes[j] + c5;
			for (size_t k = 0; k < NRAD; k++)//Loop over coeficents
			{
				r2 = distances[k] / R2;
				lnG1 = log((c1*atan(magnitudes[j] + c2) + c3)*F);
				lnG2 = -0.5*log(pow(1-r2,2)+4*D2*D2*r2);
				lnG3 = -((c11+c12*magnitudes[j]) / Q0)*distances[k];
				lnG4 = 1 - 0.5*log(VS30 / 2800);
				lnG5 = 0;

				if (frequenciesreq.at(i)>PGA)
				{
					acelerations[k] = pow(e, lnG1 + lnG2 + lnG3 + lnG4 + lnG5)/pow(e,-0.000257*distances[k]+0.27);
				}
				else {
					mu = m1*distances[k] + m2*magnitudes[j] + m3;
					I = 1.4;
					S = s1*distances[k] - (s2*magnitudes[j] + s3);
					Tsp = t1*distances[k] + t2*magnitudes[j] + t3;
					varsigma = a1*magnitudes[j] * magnitudes[j] + a2*magnitudes[j] + a3;
					kVs30 = -0.5*log(VS30 / 2800);
					fVs30 = VS30 / (120 - 1.6);
					SAnorm = I*pow(e, -0.5*pow((log(periodsreq.at(i)) + mu) / S, 2)) + pow(pow(1 - pow(periodsreq.at(i) / Tsp, varsigma), 2) + 4 * Dsp*Dsp*pow(periodsreq.at(i) / Tsp, varsigma), -0.5);
					Lin_amp = 1 + kVs30 / pow(pow(1 - (fVs30 / frequenciesreq.at(i)), 2) + 1.96*(fVs30 / frequenciesreq.at(i)), 0.5);

					acelerations[k] = pow(e, lnG1 + lnG2 + lnG3 + lnG4 + lnG5)*SAnorm*Lin_amp / (pow(e, -0.000257*distances[k] + 0.27));
				}
				attenueationtableG15 << setw(WIDTH) << acelerations[k];//Saving values

			}
			attenueationtableG15 << endl;
		}
	}




	system("pause");
	return 0;
}
