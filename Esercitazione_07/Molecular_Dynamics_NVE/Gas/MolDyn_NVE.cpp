#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <ostream>
#include <iomanip>
#include <vector>
#include "MolDyn_NVE.h"
#include "random.h"
#define _USE_MATH_DEFINES

using namespace std;

int main(){

  Input();   //Inizialization

  cout << "Temperature before the thermalization: " << temp << endl;

  int thermalizationSteps = 30000;
  ofstream Term;
  Term.open("../../Results/Molecular_Dynamics_NVE/Gas/output_restart.out",ios::app);

  //Termalizzazione del sistema partendo da una temperature più bassa e lasciando evolvere il sistema
  cout << "Thermalization of the system" << endl;
  for(int i = 0; i<thermalizationSteps;i++){
    Move();
     
    if(i%20 == 0){
	    Thermalization();
	    Term << stima_temp << endl;
    }
  }
  Term.close();
	cout << "Temperature after thermalization: " << stima_temp << endl;
	cout << endl << endl;
	cout << "Starting simulation" << endl;

  int nconf = 1;
  //Evoluzione del sistema all'equilibrio
  for(int iblk=1; iblk <= nblk; ++iblk){ 
	  cout <<endl<< "Block number " << iblk << endl;
	  Reset(iblk);   //Reset block averages
   	for(int istep=1; istep <= nstep; ++istep){
   		Move();           //Move particles with Verlet algorithm
      if(istep%10 == 0){
        Measure();   //Properties measurement
        //Valori istantanei delle grandezze
        //Print_Istant(); 
			  Accumulate(); //Update block averages
        //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
        nconf += 1;
     	}
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal();         //Write final configuration to restart
  return 0;
}

//Prepare all stuff for the simulation
void Input(){
  
  ifstream ReadInput,ReadConf;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  //Read seed for random numbers
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  input.close();

  //Lettura dei valori di input
  ReadInput.open("input.gas"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> nblk;

  cout << "Number of blocks for data blocking = " << nblk << endl;

  vtail = (8.0*M_PI*rho)/(9.0*pow(rcut,9)) - (8.0*M_PI*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*M_PI*rho)/(9.0*pow(rcut,9)) - (16.0*M_PI*rho)/(3.0*pow(rcut,3));
  cout << "Tail correction for the potential energy = " << vtail << endl;
  cout << "Tail correction for the virial           = " << ptail << endl;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

  //Prepare arrays for measurements
  iv = 0; //Potential energy
  iw = 1; //Virial
  it = 2; //Temperatue
  ik = 3; //Kinetic Energy
  itot = 4; //Total Energy
  n_props = 5; //Number of observables

  //measurement of g(r)
  igofr = 5;
  numberOfBins = 100;
  n_props = n_props + numberOfBins;
  binSize = (box/2.0)/(double)numberOfBins;

  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();
  

  //Prepare initial velocities
  cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
  double sumv[3] = {0.0, 0.0, 0.0};
  for (int i=0; i<npart; ++i){
    vx[i] = rnd.Rannyu() - 0.5;
    vy[i] = rnd.Rannyu() - 0.5;
    vz[i] = rnd.Rannyu() - 0.5;

    sumv[0] += vx[i];
    sumv[1] += vy[i];
    sumv[2] += vz[i];
  }
  for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
  double sumv2 = 0.0, fs;
  for (int i=0; i<npart; ++i){
    vx[i] = vx[i] - sumv[0];
    vy[i] = vy[i] - sumv[1];
    vz[i] = vz[i] - sumv[2];

    sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
  }
  sumv2 /= (double)npart;

  fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
  for (int i=0; i<npart; ++i){
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;

    xold[i] = x[i] - vx[i] * delta;
    yold[i] = y[i] - vy[i] * delta;
    zold[i] = z[i] - vz[i] * delta;


  }
  return;
}

//Metodo che termalizza il sistema a partire da una temperature più alta o bassa a seconda della fase che si considera lasciandolo evolvere
void Thermalization(){    
  double sumv[3] = {0.0, 0.0, 0.0};

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    vx[i] = Pbc(x[i] - xold2[i])/(2.0 * delta);
    vy[i] = Pbc(y[i] - yold2[i])/(2.0 * delta);
    vz[i] = Pbc(z[i] - zold2[i])/(2.0 * delta);
    sumv[0] += vx[i];
    sumv[1] += vy[i];
    sumv[2] += vz[i];
  }

  for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];

      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;

    Square_Velocity();
    temp = (2.0 / 3.0) * t/(double)npart;
    stima_temp = temp;

    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
    for (int i=0; i<npart; ++i){
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;
    }
  return;
}

void Velocity(void){

  for(int i=0; i<npart; ++i){
    vx[i] = Pbc(x[i] - xold2[i])/(2.0 * delta);    //r(t+dt)-r(t-dt)
    vy[i] = Pbc(y[i] - yold2[i])/(2.0 * delta);
    vz[i] = Pbc(z[i] - zold2[i])/(2.0 * delta);
  }
  return;
}

void Square_Velocity(void) { 
  t=0.;
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  return;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    //r(t-dt) 
    xold2[i] = xold[i];
    yold2[i] = yold[i];
    zold2[i] = zold[i];

    //r(t+dt)
    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    //r(t)
    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;  
    y[i] = ynew;
    z[i] = znew;
  }
  
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }

  return f;
}

void Measure(){

  double v = 0.0, w = 0.0;
  double vij, wij;
  double dx, dy, dz, dr;

  //reset the hystogram of g(r)
  for (int k=igofr; k<igofr+numberOfBins; ++k) walker[k]=0.0;

  //cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

      // distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      walker[igofr+int(dr/binSize)]+=2;

      if(dr < rcut){
        vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
        wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

        // contribution to energy and virial
        v += vij;
        w += wij;
      }
    }
  }

  walker[iv] = 4.0 * v;
  walker[iw] = 48.0 * w / 3.0;
  Velocity();
  Square_Velocity();
  Temperature();
  Kinetic_Energy();
}

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();

  return;
}

void Temperature(void){
  stima_temp = (2.0 / 3.0) * t/(double)npart; 
  walker[it] = stima_temp;
  return;
}

void Kinetic_Energy(void){
  walker[ik] = t/(double)npart;
}

void Print_Istant(){  //Stampa i valori istantanei di pressione, energia potenziale, temperature, energia totale, energia cinetica
  ofstream Epot, Pres, Ekin, Temp, Etot;
  Epot.open("../Results/Gas/output_epot.dat",ios::app);
  Pres.open("../Results/Gas/output_pres.dat",ios::app);
  Temp.open("../Results/Gas/output_temp.dat",ios::app);
  Etot.open("../Results/Gas/output_etot.dat",ios::app);
  Ekin.open("../Results/Gas/output_ekin.dat",ios::app);

  double potential = walker[iv]/(double)npart + vtail;
  double pres =  rho * temp + (walker[iw]+ ptail * (double)npart) / vol;

  //Istant Potential energy per particle
  Epot << potential << endl;
  //Istant Pressure
  Pres << pres << endl;
  //istant Temperature
  Temp << walker[it] << endl;
  //Istant kinetic energy
  Ekin << walker[ik] <<endl;
  //Istant total energy
  Etot << walker[ik]+potential <<endl;

  Epot.close();
  Pres.close();
  Ekin.close();
  Etot.close();
  Temp.close();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
  rnd.SaveSeed();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
  return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk){
  return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void Reset(int iblk){ //Reset block averages

  if(iblk == 1){
    for(int i=0; i<n_props; ++i){
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }

  for(int i=0; i<n_props; ++i){
    blk_av[i] = 0;
  }
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}

void Accumulate(void){ //Update block averages

  for(int i=0; i<n_props; ++i){
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;
}

void Averages(int iblk){ //Print results for current block

  double r, deltaV;
  ofstream Epot,Pres,Temp,Ekin,Etot, Gofr, Gave;
  Epot.open("../../Results/Molecular_Dynamics_NVE/Gas/ave_epot.out",ios::app);
  Pres.open("../../Results/Molecular_Dynamics_NVE/Gas/ave_pres.out",ios::app);
  Temp.open("../../Results/Molecular_Dynamics_NVE/Gas/ave_temp.out",ios::app);
  Ekin.open("../../Results/Molecular_Dynamics_NVE/Gas/ave_ekin.out",ios::app);
  Etot.open("../../Results/Molecular_Dynamics_NVE/Gas/ave_etot.out",ios::app);
  Gofr.open("../../Results/Molecular_Dynamics_NVE/Gas/output.gofr.0",ios::app);
  Gave.open("../../Results/Molecular_Dynamics_NVE/Gas/output.gave.0",ios::app);

  //Potential energy
  stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; 
  glob_av[iv] += stima_pot;
  glob_av2[iv] += stima_pot*stima_pot;
  err_pot = Error(glob_av[iv],glob_av2[iv],iblk);

  //Pressure
  stima_pres = rho * stima_temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; 
  glob_av[iw] += stima_pres;
  glob_av2[iw] += stima_pres*stima_pres;
  err_press = Error(glob_av[iw],glob_av2[iw],iblk);

  //Temperatura
  temperature = blk_av[it]/blk_norm; 
  glob_av[it] += temperature;
  glob_av2[it] += temperature*temperature;
  err_temp = Error(glob_av[it],glob_av2[it],iblk);


  stima_ekin = blk_av[ik]/blk_norm; //Energia Cinetica
  glob_av[ik] += stima_ekin;
  glob_av2[ik] += stima_ekin*stima_ekin;
  err_ekin = Error(glob_av[ik],glob_av2[ik],iblk);

  stima_etot = stima_ekin+stima_pot; //Energia totale
  glob_av[itot] += stima_etot;
  glob_av2[itot] += stima_etot*stima_etot;
  err_etot = Error(glob_av[itot],glob_av2[itot],iblk);

  //Potential energy per particle
  Epot << "  " << iblk <<  "  " << setprecision(12) << glob_av[iv]/(double)iblk << "  " << setprecision(12) << err_pot << endl;
  //Pressure
  Pres << "  " << iblk <<  "  " << setprecision(12) << glob_av[iw]/(double)iblk << "  " << setprecision(12) << err_press << endl;
  //Average Temperature
  Temp << "  " << iblk <<  "  " << setprecision(12) << glob_av[it]/(double)iblk << "  " << setprecision(12) << err_temp << endl;
  //Average kinetic energy
  Ekin << "  " << iblk <<  "  " << setprecision(12) << glob_av[ik]/(double)iblk << "  " << setprecision(12) << err_ekin << endl;
  //Average total energy
  Etot << "  " << iblk <<  "  " << setprecision(12) << glob_av[itot]/(double)iblk << "  " << setprecision(12) << err_etot << endl;

  //g(r)   Radial function distribution
  for (int k=igofr; k<igofr+numberOfBins; ++k) {
    r = (k-3)*binSize; 
    deltaV = (4./3.)*M_PI*(pow(r+binSize,3)-pow(r,3)); 
    stima_gdir = (1./(rho*deltaV*(double)npart))*(blk_av[k]/blk_norm);
    glob_av[k] += stima_gdir;
    glob_av2[k] += stima_gdir*stima_gdir;
    err_gdir=Error(glob_av[k],glob_av2[k],iblk);
    Gofr   <<  "  " << setprecision(12) << r << "  " << setprecision(12) << stima_gdir << endl;  
    if(iblk == nblk){  
      Gave << "  " << setprecision(12) << r <<  "  " << setprecision(12) << glob_av[k]/(double)iblk << "  " << setprecision(12) << err_gdir<< endl; //stampo l'ultimo blocco di g(r) al variare del raggio con rispettivo errore
    }
  }

  cout <<endl<< "----------------------------" << endl << endl;

  Epot.close();
  Pres.close();
  Temp.close();
  Ekin.close();
  Etot.close();
  Gofr.close();
  Gave.close();
}
