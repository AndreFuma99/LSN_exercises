#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "random.h"
#include <vector>

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
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
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
   
// L'obiettivo è generare un 3D RW su un reticolo cubico

   Walker wlk;
   int a = 1;   // costante reticolare 
   int M = 10000;   // M rappresenta il numero di "lanci"
   int N = 100;   // N rappresenta il numero di blocchi in cui divido M
   int L = M/N;   // L rappresenta il numero di lanci in ogni blocco
   std::vector<std::vector<double>> ave(N, std::vector<double>(N));   // ave è un vettore di vettori che contiene, per ogni blocco, i valori medi di distanza dall'origine per ogni passo
   std::vector<std::vector<double>> av2(N, std::vector<double>(N));   // av2 è un vettore di vettori che contiene, per ogni blocco, i valori quadratici medi di distanza dall'origine per ogni passo
   std::vector<double> Mean(N);   // Vettore ottenuto facendo la media sui blocchi mantenendo fisso il passo
   std::vector<double> Err(N);   // Vettore ottenuto facendo la deviazione standard rispetto ai blocchi mantenendo fisso il passo
   
   for(int i=0; i<N; i++){   // per ogni blocco
      std::vector<double> walks(N);    
      for(int j=0; j<L; j++){   // per ogni lancio (all'interno di un blocco)
         wlk.setx(0);
         wlk.sety(0);
         wlk.setz(0);
         for(int l=0; l<100; l++){   // per ogni passo
            double k = rnd.Rannyu();
      	    double t = rnd.Rannyu();
            int r = 0;
            if (t>0.5)   r++;   // questa condizione mi dice se il walker si sposterà in direzione positiva o negativa
            int rw = 2*a*(r-0.5);
         
      	    if (k<=1./3){
      	       wlk.setx(wlk.getx()+rw);
      	    } else if (k>=2./3){
      	       wlk.setz(wlk.getz()+rw);
      	    } else {
      	       wlk.sety(wlk.gety()+rw);
      	    }
      	    walks[l] += sqrt(wlk.getx()*wlk.getx()+wlk.gety()*wlk.gety()+wlk.getz()*wlk.getz());   // Aggiorno, a passo fissato, la distanza dall'origine
      	 }
      }
      for(int l=0; l<100; l++){   // per ogni passo, a blocco fissato
         ave[i][l] = walks[l]/L;
         av2[i][l] = ave[i][l]*ave[i][l];
      }
   }   
   
   for(int i=0; i<100; i++){   // per ogni passo
      double sum = 0.;
      double sum2 = 0.;
      for(int j=0; j<N; j++){   // per ogni blocco
         sum += ave[j][i];
         sum2 += av2[j][i];
      }
      Mean[i] = sum/N;
      Err[i] = sum2/N - Mean[i]*Mean[i];      
   }
   
   ofstream output_1;
   output_1.open("Es02_2_01.csv");
   if (output_1.is_open()){
      output_1 << "Average_Value" << "," << "Statistical_Uncertainty" << endl;
      for(int i=0; i<100; i++){
         output_1 << Mean[i] << "," << Err[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es02_2_01.csv" << endl;
   output_1.close();
   


// L'obiettivo è generare un 3D RW nel continuo

   Walker newlk;

   for(int i=0; i<N; i++){   // per ogni blocco
      for (int j=0; j<100; j++){   // per ogni passo
         ave[i][j] = 0.;
         av2[i][j] = 0.;
      }
      Mean[i] = 0.;
      Err[i] = 0.;
   }
   
   for(int i=0; i<N; i++){   // per ogni blocco
      std::vector<double> walks(N);    
      for(int j=0; j<L; j++){   // per ogni lancio (all'interno di un blocco)
         newlk.setx(0);
         newlk.sety(0);
         newlk.setz(0);
         for(int l=0; l<100; l++){   // per ogni passo
            double t = rnd.Rannyu(0., M_PI);   // sto generando l'angolo theta
      	    double f = rnd.Rannyu(0., 2*M_PI);   // sto generando l'angolo phi
      	    newlk.setx(newlk.getx()+a*sin(t)*cos(f));
      	    newlk.sety(newlk.gety()+a*sin(t)*sin(f));
      	    newlk.setz(newlk.getz()+a*cos(t));
      	                
      	    walks[l] += sqrt(newlk.getx()*newlk.getx()+newlk.gety()*newlk.gety()+newlk.getz()*newlk.getz());   // Aggiorno, a passo fissato, la distanza dall'origine
      	 }
      }
      for(int l=0; l<100; l++){   // per ogni passo, a blocco fissato
         ave[i][l] = walks[l]/L;
         av2[i][l] = ave[i][l]*ave[i][l];
      }
   }
   
   for(int i=0; i<100; i++){   // per ogni passo
      double sum = 0.;
      double sum2 = 0.;
      for(int j=0; j<N; j++){   // per ogni blocco
         sum += ave[j][i];
         sum2 += av2[j][i];
      }
      Mean[i] = sum/N;
      Err[i] = sum2/N - Mean[i]*Mean[i];      
   }
   
   ofstream output_2;
   output_2.open("Es02_2_02.csv");
   if (output_2.is_open()){
      output_2 << "Average_Value" << "," << "Statistical_Uncertainty" << endl;
      for(int i=0; i<100; i++){
         output_2 << Mean[i] << "," << Err[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es02_2_02.csv" << endl;
   output_2.close();

   rnd.SaveSeed();
   return 0;
}


