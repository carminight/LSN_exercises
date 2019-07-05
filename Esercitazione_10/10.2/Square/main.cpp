#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include "random.h"
#include "DataBlocking.h"
#include "TravelingSalesman.h"
#include "mpi.h"

using namespace std;

int main (int argc,char* argv[]){

  MPI::Init();
  int size = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();
  double tstart = MPI_Wtime();

  int numberOfCities = 30; 
  int numberOfParents = 1; 
  int numberOfGenerations = 250; 
  int numberOfMutations = 3; 
  int side = 1.;

  double InitialBeta = 0.;
  double FinalBeta = 80.;
  double increase = 0.0005;

  //Vettori necessari per raccogliere i dati provenienti dai nodi
  int path_irecv[size][numberOfCities];       //matrice che raccoglie il migliore cammino di ciascun nodo
  vector<double> distance_irecv(size,0);      //vettore che raccoglie la minor distanza di ciascun nodo

  //Con questa inizializzazione ogni nodo possiede un seed diverso e quindi diverse posizioni delle città
  TravelingSalesman walker(numberOfCities, numberOfParents, numberOfGenerations, rank+10);
  walker.Square(side);
  vector<vector<double>> positionsOfCities = walker.GetPositions();

  for(double beta = InitialBeta; beta<FinalBeta; beta += increase){
	  walker.RunMetropolis(beta, numberOfMutations);
  }

  //Da ciascun nodo prendo la migliore distanza e il miglior cammino
  double bestDistance = walker.MinimumDistance();
  int *bestPath = new int[numberOfCities];
  bestPath = walker.BestWalk();

  //Tutti i risultati vengono passati ad un singolo nodo (in questo caso il no. 0)
  MPI_Gather(bestPath, numberOfCities, MPI_INT, path_irecv[rank], numberOfCities, MPI_INT, 0, MPI::COMM_WORLD);
  MPI_Gather(&bestDistance, 1, MPI_DOUBLE, &distance_irecv[rank], 1, MPI_DOUBLE, 0, MPI::COMM_WORLD);

  //Il nodo no. 0 completa l'analisi

  if(rank==0){
	  int indexMin = min_element(distance_irecv.begin(),distance_irecv.end()) - distance_irecv.begin(); 

	  ofstream results;
    results.open("../../Results/10.2/distanceAndpath_square.out");
		results << distance_irecv[indexMin] << " " << distance_irecv[indexMin] << endl;  
  
    for (int i=0; i <numberOfCities;i++){
      results << positionsOfCities[path_irecv[indexMin][i]][1] << " " << positionsOfCities[path_irecv[indexMin][i]][2] << endl;
    }
    //Hometown
    results << positionsOfCities[path_irecv[indexMin][0]][1] << " " << positionsOfCities[path_irecv[indexMin][0]][2] << endl;


    results.close();
  }
  double tend = MPI_Wtime();
  double dt = tend - tstart;
  MPI::Finalize();

  cout << endl << "Tempo di esecuzione: " << dt << endl;
  
  delete[] bestPath;
  return 0;

}
