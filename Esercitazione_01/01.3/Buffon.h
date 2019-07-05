#ifndef _Buffon_h_
#define _Buffon_h_

#include "random.h"

class Buffon{

	private:

		Random _rnd;
		double _DistanceBetweenLines;
		double _LengthOfNeedle;

	public:

		//Costruttore di default
		Buffon();
		//Costruttore
		Buffon(double LengthOfNeedle , double DistanceBetweenLines);
		//Distruttore
		~Buffon();
		//Inizializza il seme del generatore di numeri casuali
		void Start();
		//Metodo per effetture l'esperimento di Buffon
		double Experiment(int NumberOfThrows);
		//Metodo per stimare il seno di un angolo
		double Sin_Angle();

};

#endif
