#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "random.h"
#include "RW.h"
#include "DataBlocking.h"

using namespace std;

int main (int argc, char *argv[]){

  int NumberOfIterations = 1000000;   //Numero di Random Walk totali
  int NumberOfSteps = 100;            //Numero di passi del Random Walk e coincide con il numero di blocchi che utilizzo per il data blocking

  //RANDOM WALK DISCRETO
  RW randomwalk1(1);                  //Random walk di passo reticolare 1

  int NumberOfRWPerStep = int(NumberOfIterations/NumberOfSteps);    //Numero di random walk che voglio simulare per ogni passo
  vector<double> squaredDistances1(NumberOfIterations,0);          //Vettore delle distanze al quarato

  for(int i=0; i<NumberOfSteps; i++){
		for(int j=0; j<NumberOfRWPerStep; j++){
      int position = j +i*NumberOfRWPerStep;
			squaredDistances1[position] = (double) pow(randomwalk1.DiffusionLattice(i) , 2);     
		}
	}
  Analysis(NumberOfSteps, squaredDistances1, "../Results/02.2/Results02.2.1" );

  //RANDOM WALK CONTINUO
  RW randomwalk2(1);   //Random walk di passo 1

  vector<double> squaredDistances2(NumberOfIterations,0);

  for(int i=0; i<NumberOfSteps; i++){
		for(int j=0; j<NumberOfRWPerStep; j++){
      int position = j +i*NumberOfRWPerStep;
			squaredDistances2[position] = (double) pow(randomwalk2.IsotropicDiffusion(i) , 2);     
		}
	}
  Analysis(NumberOfSteps, squaredDistances2, "../Results/02.2/Results02.2.2" );

  return 0;
}
