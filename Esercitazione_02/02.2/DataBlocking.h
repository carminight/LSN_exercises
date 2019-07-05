#ifndef _DataBlocking_h_
#define _DataBlocking_h_

#include "random.h"
#include <vector>
using namespace std;

//Somma degli elementi di un vettore tra due indici del vettore stesso
double Sum(vector<double> vector, int start, int end);
//Media di un vettore di dati
double Mean(vector<double> vector);
//Varianza di un vettore di dati
double Variance(vector<double> vector);
//Deviazione standard
double Error(vector<double> Mean, vector<double> Mean2, int n, int i);
//Data blocking applicato ad un vettore di dati con un certo numero di blocchi scelto da main
void Analysis(int NumberOfBlocks, vector<double> VectorOfData, const char*OutputFile);
//Test del chi quadro
void ChiQuadro(int NumberOfData, int NumberOfBlocks, vector<double> VectorOfData, const char* OutputFile);
#endif
