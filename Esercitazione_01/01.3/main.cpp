#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "random.h"
#include "Buffon.h"
#include "DataBlocking.h"

using namespace std;

int main (int argc, char *argv[]){

	int NumberOfThrows = 20000;
	int NumberOfBlocks = 100;

	double _DistanceBetweenLines = 6.;
	double LengthOfNeedle = 4.;

	Buffon Test(LengthOfNeedle,_DistanceBetweenLines);
	//VEttore che conterrà i risultati dell'esperimento
	vector<double> pi_vector(NumberOfThrows,0);

  	for(int j=0; j<NumberOfThrows; j++){
    	pi_vector.at(j) = Test.Experiment(NumberOfThrows);
	}

	//Data blocking
  	Analysis(NumberOfBlocks, pi_vector, "../Results/01.3/Results01.3.out");

  	return 0;
}
