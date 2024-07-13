#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>

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
   
// Il mio obiettivo è quello di verificare il teorema del limite centrale

   double L = 1.;   // tasso di decadimento
   double G = 1.;   // larghezza della lorentziana
   int M = 10000;   // numero di trials
   int N = 100;   // numero di lanci
   std::vector<double> std_dice(M);   // dado standard
   std::vector<double> exp_dice(M);   // dado esponenziale
   std::vector<double> lor_dice(M);   // dado lorentziano

   for(int i=0; i<M; i++){
   
      double sum_s = 0.;
      double sum_e = 0.;
      double sum_l = 0.;
      
      for(int j=0; j<N; j++){
         double k = rnd.Rannyu();
         sum_s += k;
         sum_e += -log(1-k)/L;
         sum_l += G*tan(M_PI*(k-0.5));
      }
      
      std_dice[i] = sum_s/N;
      exp_dice[i] = sum_e/N;
      lor_dice[i] = sum_l/N;
      
   }
   
   ofstream outputS;
   outputS.open("Es01_2_std.csv");
   if (outputS.is_open()){
      outputS << "Standard_Sums" << endl;
      for(int i=0; i<M; i++){
         outputS << std_dice[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es01_2_std.csv" << endl;
   outputS.close();
   
   ofstream outputE;
   outputE.open("Es01_2_exp.csv");
   if (outputE.is_open()){
      outputE << "Exponential_Sums" << endl;
      for(int i=0; i<M; i++){
         outputE << exp_dice[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es01_2_exp.csv" << endl;
   outputE.close();
   
   ofstream outputL;
   outputL.open("Es01_2_lor.csv");
   if (outputL.is_open()){
      outputL << "Lorentzian_Sums" << endl;
      for(int i=0; i<M; i++){
         outputL << lor_dice[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open Es01_2_lor.csv" << endl;
   outputL.close();

   rnd.SaveSeed();
   return 0;
}


