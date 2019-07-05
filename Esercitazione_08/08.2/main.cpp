#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "DataBlocking.h"
#include "funcvirtual.h"
#include "Metropolis.h"

using namespace std;

int main (int argc, char *argv[]){

  int dimension = 1;
  double metropolisStep = 2.6;
  int numberOfSteps = 10000;
  int numberOfIntegrals = 10000;
  int numberOfBlocks = 100;
  int numberOfBins = 200;
  int numberOfEquilibrationStep = 100;
  
  WaveFunction* PsiTrial = new WaveFunction(); 
  //Impongo i parametri ottimali ottenuti da Variational MonteCarlo
  PsiTrial->SetMi(0.795);
  PsiTrial->SetSigma(0.61);
  Function* H = new Hamiltonian(PsiTrial);

	vector<double> InitialPositions(dimension,0);
	InitialPositions.at(0) = 0;
  Metropolis walker(PsiTrial, InitialPositions, metropolisStep, numberOfSteps, dimension, numberOfIntegrals, numberOfBins);
  vector<double> results_integration(numberOfIntegrals,0);

  walker.Equilibrate(numberOfEquilibrationStep, "uniform");
  for(int i=0; i<numberOfIntegrals; i++){

    walker.RunUniform();
    results_integration.at(i) = walker.Integrate(H);
    if(i%10000){
      walker.AcceptanceRate();
    }
    walker.Print("../Results/08.2/psi.out");
    walker.Histogram(i);
    walker.Restart();
  }

  walker.AnalysisHistogram(numberOfBlocks,"../Results/08.2/configurations.out");                //DataBlocking applicato all'instogramma
  Analysis(numberOfBlocks, results_integration,"../Results/08.2/hamiltonian_results.out");      //DataBlocking degli integrali

  delete PsiTrial;
  delete H;
  return 0;

}
