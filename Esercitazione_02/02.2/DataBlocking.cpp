#include "DataBlocking.h"
#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

//Somma degli elementi di un vettore tra due indici del vettore stesso
double Sum(vector<double> vector, int start, int end){
	double sum = 0;
	for(int i=start; i<end; i++) sum += vector.at(i);
	return sum;
}

//Media di un vettore di dati
double Mean(vector<double> vector) {
	double sum = 0;
	for(int i=0; i<vector.size(); i++) sum += vector.at(i);
	return sum/vector.size();
}

//Varianza di un vettore di dati
double Variance(vector<double> vector) {
	double sum=0;
	for(unsigned int i=0; i<vector.size(); i++){
		sum+=pow(Mean(vector)-vector.at(i), 2.);
	}
	return sum/(vector.size()-1) ;
 }

//Deviazione standard (avendo gia vettore delle medie e vettore dei quadrati delle medie)
double Error(vector<double> Mean, vector<double> Mean2, int n,int i){
	if (n==0){
    return 0;
	}
	else{
  	return sqrt((Mean2.at(n) - pow(Mean.at(n),2))/i);
	}
}

//Data blocking adattato al Random Walk 
void Analysis(int NumberOfBlocks, vector<double> VectorOfData, const char*OutputFile){

	ofstream results;
	results.open(OutputFile);

	int NumberOfData = VectorOfData.size();
	int NumberOfDataPerBlocks = int(NumberOfData/NumberOfBlocks);
	int NumerOfSubIntervals = NumberOfDataPerBlocks/NumberOfBlocks;	

	for(int i=0; i<NumberOfBlocks; i++){

		vector<double> vmean(NumerOfSubIntervals,0);
		vector<double> vmean2(NumerOfSubIntervals,0);

		for(int j=0; j<NumerOfSubIntervals;j++){

			double temp = sqrt( Sum(VectorOfData, i*NumberOfDataPerBlocks + j*NumerOfSubIntervals, i*NumberOfDataPerBlocks + (j+1)*NumerOfSubIntervals)/NumerOfSubIntervals);
			vmean.at(j) = temp;				//Media del sottoblocco
			vmean2.at(j) = pow(temp,2);		//Media al quadrato del sotto blocco

		}

		double mean_cum = 0.;
		double mean2_cum = 0.;
		double err_cum = 0.;

		mean_cum = Mean(vmean);
		mean2_cum = Mean(vmean2);
		err_cum = sqrt( ( mean2_cum - pow(mean_cum,2) ) / (NumerOfSubIntervals - 1)  ) ;

		if(results.is_open()) results << i << " " << setprecision(12) << mean_cum << " " << setprecision(12) << err_cum << endl;
		else cerr << "PROBLEM: Unable to open result.txt" << endl;
	}

	results.close();

}

void ChiQuadro(int NumberOfIntervals, int TotalNumberOfData, vector<double> VectorOfData, const char* OutputFile){

	ofstream results;
  results.open(OutputFile);

	//Numero di dati sui cui calcolo il chi quadro singolo
	int NumberOfDataPerInterval = int(TotalNumberOfData/NumberOfIntervals);

	vector<double> chi2_i(NumberOfIntervals,0); //Vettore contenente i chi quadri i-esimi
	vector<double> intervals(NumberOfIntervals+1,0);		//Vettore degli intervalli

	//Vettore che contiene l'estremo destro di ogni intervallo in cui è diviso [0,1)
	for(int i=0; i<NumberOfIntervals+1; i++){
		intervals.at(i) = (1.*i)/NumberOfIntervals;
	}

	for(int i=0; i<NumberOfIntervals; i++){
		//Vettore che conterrà gli osservati per ogni intervallo
		vector<double> observed(NumberOfIntervals,0);

		//Conteggio dei dati nei singoli intervalli
    for(int j=i*NumberOfDataPerInterval; j<(i+1)*NumberOfDataPerInterval; j++){
      for(int a=0; a<NumberOfIntervals+1; a++){
        if( (VectorOfData.at(j) >= intervals.at(a) ) && ( VectorOfData.at(j) < intervals.at(a+1) ) ){
          observed.at(a) ++ ;
        }
      }
    }

		//Calclo del chi quadro per l'i-esimo blocco
    double sum = 0;
    double expected = ((1.*NumberOfDataPerInterval)/(1.*NumberOfIntervals));		//Conteggi attesi nel signolo blocco
    for(int k=0; k<NumberOfIntervals; k++){
      sum +=  pow( observed.at(k) - expected ,2) / expected;
    }
    chi2_i.at(i) = sum;
  }

	for(int i=1; i < chi2_i.size(); i++){
  	if(results.is_open()) results << i << " " << chi2_i.at(i) << endl;
  	else cerr << "PROBLEM: Unable to open data.out" << endl;
  }
  results.close();
}
