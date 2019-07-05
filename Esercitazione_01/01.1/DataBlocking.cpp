#include "DataBlocking.h"
#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>	//Per settare la il numero di cifre significative

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

//Data blocking applicato ad un vettore di dati con un certo numero di blocchi scelto da main
void Analysis(int NumberOfBlocks, vector<double> VectorOfData, const char*OutputFile){

	ofstream results;
	results.open(OutputFile);

	unsigned int NumberOfData = VectorOfData.size();
	int NumberOfDataPerBlocks = int(NumberOfData/NumberOfBlocks);
	vector<double> vmean(NumberOfBlocks,0);   //Vettore che mi contiene le medie dei singoli blocchi
	vector<double> vmean2(NumberOfBlocks,0);   //Vettore che contiene i quadrati delle medie dei singoli blocchi
	vector<double> mean_cum(NumberOfBlocks,0);  //Vettore che mi contiene la media cumulativa dei blocchi
	vector<double> mean2_cum(NumberOfBlocks,0);   //Vettore che contiene i quadrati delle medie cumulative dei blocchi
	vector<double> err_cum(NumberOfBlocks,0);      //Vettore che mi contiene gli errori cumulativi

	double sum;
	int k = 0;
	for(int i=0; i<NumberOfBlocks; i++){
		sum = 0;
		for(int j=0; j<NumberOfDataPerBlocks; j++){
			k = j+i*NumberOfDataPerBlocks;
			sum += VectorOfData.at(k);
		}
		vmean.at(i) = sum/NumberOfDataPerBlocks;
		vmean2.at(i) = pow(vmean.at(i),2);
	}

	for(int i=0; i<NumberOfBlocks; i++){
		for(int j=0; j<i+1; j++){
			mean_cum.at(i) += vmean.at(j);
			mean2_cum.at(i) += vmean2.at(j);
		}
		mean_cum.at(i) /= (i+1);
		mean2_cum.at(i) /= (i+1);
		err_cum.at(i) = Error(mean_cum, mean2_cum, i, i+1);
	}

	for(int i=0;i<NumberOfBlocks;i++){
		if (results.is_open()){
			results <<i*NumberOfDataPerBlocks<< " "<< setprecision(12) <<mean_cum[i] << " " <<setprecision(12) << err_cum[i]<< endl;
	 	}
		else cerr << "PROBLEM: Unable to open data.out" << endl;
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
  	if(results.is_open()) results << i << " " << setprecision(12) << chi2_i.at(i) << endl;
  	else cerr << "PROBLEM: Unable to open data.out" << endl;
  }
  results.close();
}
