#include "Buffon.h"
#include "random.h"
#include <math.h>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

//Costruttore di default
Buffon::Buffon(){

	_DistanceBetweenLines = 0.;
	_LengthOfNeedle = 0.;

}

//Costruttore
Buffon::Buffon(double LengthOfNeedle , double DistanceBetweenLines){

	Start();

	_DistanceBetweenLines = DistanceBetweenLines;
	_LengthOfNeedle = LengthOfNeedle;

}
//Distruttore
Buffon:: ~Buffon(){
	_rnd.SaveSeed();
}

void Buffon::Start(){

	int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
     Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.in");
  string property;
  if (input.is_open()){
     while ( !input.eof() ){
        input >> property;
        if( property == "RANDOMSEED" ){
           input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
           _rnd.SetRandom(seed,p1,p2);
        }
     }
     input.close();
  } else cerr << "PROBLEM: Unable to open seed.in" << endl;
}

double Buffon::Sin_Angle(){
	double x;
	double y;

	do{
		x = _rnd.Rannyu(-1,1);
		y = _rnd.Rannyu(0,1);

	}
	//Controllo se il numero generato è all'interno della semicirconferenza
	while((pow(x,2)+pow(y,2)) >= 1);

	return y/sqrt(pow(x,2)+pow(y,2));
}

double Buffon::Experiment(int NumberOfThrows){

	double center;
	double sin;
	int NumberOfHit = 0;
	double y1,y2;

	for(int i=0; i< NumberOfThrows; i++){
		//Genero il centro dell'ago
		center = _rnd.Rannyu(-_DistanceBetweenLines * 0.5, _DistanceBetweenLines * 0.5);
		//Genero l'angolo rispetto alla normale
		sin = Sin_Angle();

		//Ricavo gli estremi dell'ago
		y1 = center + sin * _LengthOfNeedle * 0.5;
		y2 = center - sin * _LengthOfNeedle * 0.5;

		//Controllo se l'ago ha intersecato una riga
		if((y1>=0 && y2<=0) || (y1<=0 && y2>=0)){
			NumberOfHit++;
		}
	}
	//Estimatore del pigreco
	return (2. * _LengthOfNeedle * NumberOfThrows) / (NumberOfHit * _DistanceBetweenLines);
}
