#ifndef __MolDyn_NVE__
#define __MolDyn_NVE__

#include "random.h"

//Parameters, observables
const int m_props=3000;
int n_props, iv, iw, it,ik,itot;
double vtail,ptail;
double walker[m_props];

//Random generator
Random rnd;
int seed[4];

//Numero di blocchi
int nblk;

//Configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part];
double xold[m_part],yold[m_part],zold[m_part];
double xold2[m_part],yold2[m_part],zold2[m_part];
double vx[m_part],vy[m_part],vz[m_part];
double t; //variabile per il calcolo della velocità quadratica media del sistema

//Averages
double acc,att;
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_pot,stima_pres,err_pot,err_press,stima_temp,temperature,err_temp,stima_ekin,err_ekin,stima_etot,err_etot;

//Thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut, mode;

//Simulation
int nstep, iprint;
double delta;

//Functions
void Input();
void Move(void);
void Restart(void);
void ConfFinal(void);
void OldConfig1(void);
void OldOldConfig(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
void Thermalization(void);
void Square_Velocity(void);
void Velocity(void);
void Kinetic_Energy(void);
void Total_Energy(void);
void Temperature(void);
void Print_Istant(void);

//Blocking
void Reset(int);
void Accumulate(void);
void Averages(int);
double Error(double,double,int);

#endif // __MolDyn_NVE__
