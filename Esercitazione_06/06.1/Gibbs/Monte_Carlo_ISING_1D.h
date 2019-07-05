#ifndef __Ising_
#define __Ising_

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props,iu,ic,im,ix,ig;
double nbins;
double walker[m_props];
double Kb = 8.62e-5;

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_u,stima_c,stima_m,stima_x,stima_g,stima_u2;
double err_u,err_c,err_m,err_x,err_g;
double energy, energy2;

//configuration
const int m_spin=50;
double s[m_spin];
int mode;                   //Se è 0 la configurazione iniziale è quello a T infinito, se è 1 parto da una vecchia configurazione

// thermodynamical state
int nspin;
double beta,temp,J,h;

// simulation
int nstep, nblk, metro;

//functions
void Input(void);
void Reset(int);
void Equilibration(int);
void Accumulate(void);
void Averages(int);
void Move(int);
void ConfFinal(void);
void OldConfig(void);
void Measure(void);
int Pbc(int);
double Error(double,double,int);
void Metropolis(int o);
void Gibbs(int o);
void Energy(int);
void Susceptibility(int);
void SpecificHeath(int);
void Magnetization(int);
void Magnetization_h(int);

#endif
