#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "DataBlocking.h"


using namespace std;

int main (int argc, char *argv[]){

  Random rnd;
  //inzializzo il seme del generatore di numeri casuali
  SetRandomGenerator(&rnd);

  ofstream results;
  results.open("../Results/01.2/Results01.2.out");

  int NumberOfThrows = 10000;
  int NumberOfConfigurations[4] = {1,2,10,100};
  //Cosi posso utilizzare la classe vector
  vector<int> VectorOfConfiguration (NumberOfConfigurations,NumberOfConfigurations + sizeof(NumberOfConfigurations) / sizeof(int));

  vector<double> UniformMean(NumberOfThrows);
  vector<double> ExponentialMean(NumberOfThrows);
  vector<double> LorentianMean(NumberOfThrows);

  //Ciclo sul numero di configuarazioni
  for(int d=0; d<VectorOfConfiguration.size(); d++){

    //Vettori di appoggio
    vector<double> std(VectorOfConfiguration.at(d));
    vector<double> exp(VectorOfConfiguration.at(d));
    vector<double> lor(VectorOfConfiguration.at(d));

    //Ciclo sul numero di lanci da effettuare
    for(int i=0;i<NumberOfThrows;i++){
    	for(int j=0; j<VectorOfConfiguration.at(d); j++){
        std.at(j) = rnd.Rannyu();
        exp.at(j) = rnd.Exponential(1); // Lambda = 1
        lor.at(j) = rnd.Lorentzian(1,0); // Gamma = 1, Mu = 0
      }
      UniformMean.at(i) = Mean(std);
      ExponentialMean.at(i) = Mean(exp);
      LorentianMean.at(i) = Mean(lor);

      if (results.is_open()) results << setprecision(12) << UniformMean.at(i) << "  " << setprecision(12) << ExponentialMean.at(i) << "  " << setprecision(12) << LorentianMean.at(i) << endl;
      else cerr << "PROBLEM: Unable to open results2.txt" << endl;
    }
  }

  results.close();
  rnd.SaveSeed();
  return 0;
}
