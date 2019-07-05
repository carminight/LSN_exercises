#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "random.h"
#include "Option.h"
#include "DataBlocking.h"

using namespace std;

int main (int argc, char *argv[]){

  double InitialPrice = 100.;
  double StrikePrice = 100.;
  double RiskFreeRate = 0.1;
  double Volatility = 0.25;
  double TimeToMaturity = 1;
  double InitialTime = 0.;

  int NumberOfSimulation = 1000000;
  int NumberOfBlocks = 1000;
  int NumberOfSteps = 100;

  Option asset(InitialPrice,  StrikePrice,  RiskFreeRate, Volatility);

  //Vettori che contengono i prezzi ottenuti dalle simulazioni continue e discrete
  vector<double> putDirectly(NumberOfSimulation,0);
  vector<double> callDirectly(NumberOfSimulation,0);
  vector<double> putDiscretized(NumberOfSimulation,0);
  vector<double> callDiscretized(NumberOfSimulation,0);

  for(int i = 0; i< NumberOfSimulation; i++){
    putDirectly.at(i) = asset.PutPrice(InitialTime,TimeToMaturity,1);
    callDirectly.at(i) = asset.CallPrice(InitialTime,TimeToMaturity,1);
    putDiscretized.at(i) = asset.PutPrice(InitialTime,TimeToMaturity,NumberOfSteps);
    callDiscretized.at(i) = asset.CallPrice(InitialTime,TimeToMaturity,NumberOfSteps);
  }

  Analysis(NumberOfBlocks, putDirectly, "../Results/ResultsPutDirectly.out");
  Analysis(NumberOfBlocks, callDirectly, "../Results/ResultsCallDirectly.out");
  Analysis(NumberOfBlocks, putDiscretized, "../Results/ResultsPutDiscretized.out");
  Analysis(NumberOfBlocks, callDiscretized, "../Results/ResultsCallDiscretized.out");

  return 0;
}
