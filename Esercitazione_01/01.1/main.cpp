#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "DataBlocking.h"

using namespace std;

int main (){

  Random rnd;
  //inzializzo il seme del generatore di numeri casuali
  SetRandomGenerator(&rnd);

  unsigned int NumberOfData = 100000;
  unsigned int NumberOfBlocks = 1000;

  vector<double> Data1(NumberOfData,0);
  vector<double> Data2(NumberOfData,0);

  //Genero i dati
    //Primo esercizio
  for(int j=0; j<NumberOfData; j++){
    Data1.at(j) = rnd.Rannyu();
  }
    //Secondo esercizio
  for(int i=0; i<NumberOfData; i++){
    Data2.at(i) = pow((rnd.Rannyu() -0.5),2);
  }

  //Primo punto dell'esercizio
  Analysis(NumberOfBlocks, Data1, "../Results/01.1/Results01.1.1.out");

  //Secondo punto dell'esercizio
  Analysis(NumberOfBlocks, Data2, "../Results/01.1/Results01.1.2.out");

  //Terzo punto dell'esercizio
  int NumberOfIntervals = 100;  //numero di intervalli in cui divido l'intervallo [0,1)
  int TotalNumberOfData = 1000000;  //numero totale di numeri casuali da generare

  vector<double> chi2(NumberOfIntervals);
  vector<double> Data3(TotalNumberOfData);

  for(int i=0; i<TotalNumberOfData; i++){
	  Data3.at(i) = rnd.Rannyu();
  }

  ChiQuadro(NumberOfIntervals, TotalNumberOfData, Data3,"../Results/01.1/Results01.1.3.out");

  rnd.SaveSeed();
  return 0;
}
