#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(){

  Input();                        //Inizialization
  int equilibrationSteps = 1000;
//  int oldconfig = 10;           //Indice per indicare quando voglio stampare una vecchia configurazione. Da scommentare qualora lo si desiderasse
  cout << endl << "Equilibration of the system..." << endl;
  Equilibration(equilibrationSteps);

	do{
		beta = 1./temp;  
		for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
			accepted = 0;
			Reset(iblk);   //Reset block averages
			for(int istep=1; istep <= nstep; ++istep){
				Move(metro);
				Measure();
				Accumulate(); //Update block averages
			}
			if(iblk == nblk){
			  cout << "Temperature " << temp << endl;
			  if(metro ==1) cout << "Acceptance rate " << accepted/attempted << endl << endl;
      }
			Averages(iblk);   //Print results for current block
			if(iblk == nblk) cout << "----------------------------" << endl << endl;
//    if(iblk == oldconfiguration) OldConfig();     //Da scommentare qualora lo si volesse
		}
		temp += 0.05; 
	}
	while(temp <= 2.1);

  cout << endl << endl << "Switching on the external field..." << endl;

  temp = 0.5; //Inizializzazione della temperature per le misure a campo acceso
	h = 0.02;   //Accensione del campo esterno
	beta = 1./temp;  
	cout << endl << "External field : " << h << endl;
  cout << endl << "Equilibration of the system..." << endl;
  Equilibration(equilibrationSteps);

	do{
		beta = 1./temp;   //aggiornamento del parametro beta
		for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
			Reset(iblk);   //Reset block averages
			for(int istep=1; istep <= nstep; ++istep){
				Move(metro);
				Measure();
				Accumulate(); //Update block averages
			}
			if(iblk == nblk){
			cout << "Temperature " << temp << endl;
		  if(metro ==1) cout << "Acceptance rate " << accepted/attempted << endl << endl;
      }
			Magnetization_h(iblk); 
			if(iblk == nblk) cout << "----------------------------" << endl << endl;
		}
		temp += 0.05; 
	}
  while(temp <= 2.1);
  
  ConfFinal(); //Write final configuration

    return 0;
}


void Input(void){

  ifstream ReadInput;

  cout << endl << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  input.close();

//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;

  ReadInput >> metro; //If==1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> mode;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //SpecificHeat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility

  n_props = 4; //Number of observables

  //Scelgo se partire da una configurazione random oppure una precedente
  if(mode == 1) {  
		cout << "Read old configuration from file old.0 " << endl << endl;
		ifstream ReadConf;
		ReadConf.open("old.0");

		for (int i=0; i<nspin; ++i){
			ReadConf >> s[i] ;
		}
		ReadConf.close();

	}

	else if(mode == 0){ //Random configuration                                   
	cout << "Initial configuration with T = infinity" << endl;
	for (int i=0; i<nspin; ++i){
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
}

//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro){
  int o;

  for(int i=0; i<nspin; ++i){
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1){ //Metropolis
      Metropolis(o);
    }
    else{ //Gibbs sampling
      Gibbs(o);
    }
  }
}

//Equilibrazione del sistema
void Equilibration(int numberOfequilibrationSteps){
  for(int istep=1; istep <= numberOfequilibrationSteps; ++istep){ 
		Move(metro);
	}

}

void Measure(){

  double u = 0.0, s_tot=0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i){
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);   
    s_tot += s[i];
  }

  walker[iu] = u;
  walker[ic] = pow(u,2);
  walker[im] = s_tot;              //Magnetizzazione
  walker[ix] = pow(s_tot,2);       //Suscettività
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
	Energy(iblk);
	SpecificHeath(iblk);
	Susceptibility(iblk);
  Magnetization(iblk);
}

void Energy(int iblk) {
	ofstream Ene, ResultsU;

	Ene.open("../../Results/Metropolis/output.energy.0",ios::app);
  stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
  glob_av[iu]  += stima_u;
  glob_av2[iu] += stima_u*stima_u;
  err_u=Error(glob_av[iu],glob_av2[iu],iblk);
	energy = glob_av[iu]/(double)iblk;
  Ene << "  " << iblk <<  "  " << setprecision(12) << stima_u << "  " << setprecision(12) << energy << "  " << setprecision(12) << err_u << endl;
  Ene.close();

	if(iblk == nblk){
   	ResultsU.open("../../Results/Metropolis/energy.result.0",ios::app);  
    ResultsU << "  " << temp <<  "  " << setprecision(12) << glob_av[iu]/(double)iblk << "  " << setprecision(12) << err_u << endl;
		ResultsU.close();
	}
}

void SpecificHeath(int iblk) {
	ofstream SpecificHeat, ResultsC;

	SpecificHeat.open("../../Results/Metropolis/output.specific_heat.0",ios::app);
  
  stima_u2 = blk_av[ic]/blk_norm/(double)nspin; 
  stima_c = pow(beta,2) * (stima_u2*(double)nspin - pow(stima_u,2) * (double)nspin*nspin)/(double)nspin; //Colore specifico
  glob_av[ic]  += stima_c;
  glob_av2[ic] += stima_c*stima_c;
  err_c=Error(glob_av[ic],glob_av2[ic],iblk);
  SpecificHeat << "  " << iblk <<  "  " << setprecision(12) << stima_c << "  " << setprecision(12) << glob_av[ic]/(double)iblk << "  " << setprecision(12) << err_c << endl;
  SpecificHeat.close();

  if(iblk == nblk){
    ResultsC.open("../../Results/Metropolis/specific_heat.result.0",ios::app);  
    ResultsC << "  " << temp <<  "  " << setprecision(12) << glob_av[ic]/(double)iblk << "  " << setprecision(12) << err_c << endl;
    ResultsC.close();
	}
}

void Susceptibility(int iblk){
	ofstream Susceptibility, ResultsX;

	Susceptibility.open("../../Results/Metropolis/output.susceptibility.0",ios::app);
  stima_x = beta*blk_av[ix]/blk_norm/(double)nspin; //Suscettività
  glob_av[ix]  += stima_x;
  glob_av2[ix] += stima_x*stima_x;
  err_x=Error(glob_av[ix],glob_av2[ix],iblk);
  Susceptibility << "  " << iblk << setprecision(12) <<  "  " << stima_x << "  " << setprecision(12) << glob_av[ix]/(double)iblk << "  " << setprecision(12) << err_x << endl;
  Susceptibility.close();
  if(iblk == nblk){
    ResultsX.open("../../Results/Metropolis/susceptibility.result.0",ios::app); 
    ResultsX << "  " << temp <<  "  " << setprecision(12) << glob_av[ix]/(double)iblk << "  " << setprecision(12) << err_x << endl;
    ResultsX.close();
  }
}

void Magnetization(int iblk){
	ofstream Magnetization, ResultsM;

	Magnetization.open("../../Results/Metropolis/output.magnetization.0",ios::app);
	stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetizzazione
	glob_av[im]  += stima_m;
	glob_av2[im] += stima_m*stima_m;
	err_m=Error(glob_av[im],glob_av2[im],iblk);
	Magnetization << "  " << iblk <<  "  " << setprecision(12) << stima_m << "  " << setprecision(12) << glob_av[im]/(double)iblk << "  " << setprecision(12) << err_m << endl;
	Magnetization.close();

	if(iblk == nblk){
   	ResultsM.open("../../Results/Metropolis/magnetization.result.0",ios::app);  
		ResultsM << "  " << temp <<  "  " << setprecision(12) << glob_av[im]/(double)iblk << "  " << setprecision(12) << err_m << endl;
		ResultsM.close();
	}
}

void Magnetization_h(int iblk){
	ofstream Magnetization_h, ResultsM_h;

	Magnetization_h.open("../../Results/Metropolis/output.magnetization_h.0",ios::app);
	stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetizzazione
	glob_av[im]  += stima_m;
	glob_av2[im] += stima_m*stima_m;
	err_m=Error(glob_av[im],glob_av2[im],iblk);
	Magnetization_h << "  " << iblk <<  "  " << setprecision(12) << stima_m << "  " << setprecision(12) << glob_av[im]/(double)iblk << "  " << setprecision(12) << err_m << endl;
	Magnetization_h.close();

	if(iblk == nblk){
   	ResultsM_h.open("../../Results/Metropolis/magnetization_h.result.0",ios::app);  
		ResultsM_h << "  " << temp <<  "  " << setprecision(12) << glob_av[im]/(double)iblk << "  " << setprecision(12) << err_m << endl;
		ResultsM_h.close();
	}


}

void ConfFinal(void){
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i){
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

//Metodo per stampare la configurazione attuale in un file
void OldCongig(void){
  ofstream WriteConf;

  cout << "Print the configuration to file old.0 " << endl << endl;
  WriteConf.open("old.0");
  for (int i=0; i<nspin; ++i){
    WriteConf << s[i] << endl;
  }
  WriteConf.close();
}

int Pbc(int i){  //Algorithm for periodic boundary conditions

  if(i >= nspin) i = i - nspin;
  else if(i < 0) i = i + nspin;
  return i;
}

double Error(double sum, double sum2, int iblk){
  return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void Metropolis(int o){

  double deltaE = 2.* J * s[o] * ( s[Pbc(o-1)] + s[Pbc(o+1)] ) + h * 2. * s[o];
  if(deltaE < 0) s[o] = s[o]*(-1), accepted++;
  
  else{
    double a = rnd.Rannyu();
    double boltzmannweight = exp(-beta * deltaE);
		if( a <= boltzmannweight){
			s[o] = s[o]*(-1);
      accepted++;
    }
  }
  attempted++;
}

void Gibbs(int o){
  
  //Prendo in considerazione solo il caso -1 -> +1 poichè il contrario è (1 - p)
  double deltaE = 2. * J * ( s[Pbc(o-1)] + s[Pbc(o+1)]) + h * 2. ;
  double p_up = 1. / (1. + exp(- beta * deltaE));
  double a = rnd.Rannyu();

  if( a <= p_up) s[o] = 1;

  else s[o] = -1;

}
