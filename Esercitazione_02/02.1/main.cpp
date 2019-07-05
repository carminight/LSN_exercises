#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "random.h"
#include "DataBlocking.h"
#include "Integrate.h"
#include "funcvirtual.h"

using namespace std;

int main (int argc, char *argv[]){

  int NumberOfIntegrals = 10000;
  int NumberOfBlocks = 100;

  //Primo punto
  //Funzione da integrare
  Function* coseno = new Coseno();
  Integrate Integrale1(0,1, coseno); 

  vector<double> vectorOfIntegrals1(NumberOfIntegrals,0);
  for(int i=0; i<NumberOfIntegrals; i++){

    vectorOfIntegrals1.at(i) = Integrale1.Media(NumberOfIntegrals);

  }

  Analysis(NumberOfBlocks, vectorOfIntegrals1, "../Results/02.1/Results02.1.1.out");

  //Secondo punto
  //Funzione da integrare
  Function* g = new gfunc(); 
  Integrate Integrale2(0,1, g); 

  vector<double> vectorOfIntegrals2(NumberOfIntegrals,0);
  for(int i=0; i<NumberOfIntegrals; i++){

    vectorOfIntegrals2.at(i) = Integrale2.MediaRetta(NumberOfIntegrals);

  }

  Analysis(NumberOfBlocks, vectorOfIntegrals2, "../Results/02.1/Results02.1.2.out");

  delete coseno;
  delete g;

  return 0;
}
