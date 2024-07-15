#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "random.h"
#include "city.h"
#include <vector>
#include <algorithm>
#include <numeric>

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
   


   // Genero le città
   
   vector<City> cities;
   cities.reserve(34);
   
   for (int i=0; i<34; i++){
      double theta = rnd.Rannyu(0,2*M_PI);
      double x = cos(theta);
      double y = sin(theta);
      cities.emplace_back(x, y, i+1);
   }
   
   // Creo la popolazione iniziale -> genero randomicamente 900 percorsi, ossia 900 vector che contengono gli indici delle città
   
   int dim = 900;
   vector<vector<int>> population(dim, vector<int>(35));
   
   // Generazione di una sequenza di indici delle città da 1 a 34
   vector<int> cityIndices(35);
   iota(cityIndices.begin(), cityIndices.end(), 1); // Riempie cityIndices con {1, 2, ..., 34}
   
   // Definizione dell'indice comune per la prima posizione di ogni vettore
   int commonIndex = 1;
   
   for (int i = 0; i < dim; ++i) {
      // Copiare gli indici delle città
      vector<int> indices = cityIndices;

      // Mescolare gli indici per ottenere una sequenza casuale
      double g = rnd.Rannyu();
      Shuffle(indices, g);

      // Impostare il primo e l'ultimo elemento al valore comune
      indices[0] = commonIndex;
      indices[34] = commonIndex;

      // Assegnare la sequenza al vettore i-esimo di population
      population[i] = indices;
    }
   
   double p_m1 = 0.08;   // probabilità della prima mutazione
   double p_m2 = 0.06;   // probabilità della seconda mutazione
   double p_m3 = 0.07;   // probabilità della terza mutazione
   double p_m4 = 0.09;   // probabilità della quarta mutazione
   double p_c = 0.8;   // probabilità del crossover
   vector<double> LOSS_BEST(120);   // vettore che contiene la loss function del miglior cammino per ogni generazione
   vector<double> LOSS_HALF(120);   // vettore che contiene la loss function, mediata sulla migliore metà della popolazione, per ogni generazione
   
   for (int i=0; i<120; i++){   // Faccio 120 generazioni perchè ho visto che dopo circa 95 L si stabilizza
      vector<double> Loss(dim);
      for (int j=0; j<dim; j++){   // calcolo la loss function per ogni cammino
         vector<double> X(35);
         vector<double> Y(35);
         for (int k=0; k<35; k++){
            X[k] = cities[population[j][k]].getX();
            Y[k] = cities[population[j][k]].getY();
         }
         Loss[j]=loss(X, Y);
      }
      
      sortPopulationByLoss(Loss, population);
      LOSS_BEST[i] = Loss[0];
      double sum = 0;
      for (int q=0; q<(dim/2); q++)   sum += Loss[q];
      LOSS_HALF[i] = sum/(dim/2.);
      
      vector<vector<int>> new_population(dim, vector<int>(35));
      
      for (int j=0; j<(dim/2); j++){   // Applico la selezione 900/2 volte per ottenere 900 figli 
         int m = selection(rnd.Rannyu(), dim);
         int n = selection(rnd.Rannyu(), dim);
         while(m==n){
            m = selection(rnd.Rannyu(), dim);
            n = selection(rnd.Rannyu(), dim);
         }
         
         new_population[2*j] = population[m];
         new_population[2*j + 1] = population[n];
         
         // Ad ogni coppia di figli applico crossover e mutazione
         double C = rnd.Rannyu();
         if (C<p_c){
            double cut_point = rnd.Rannyu(1, 34);   // genero la posizione in cui tagliare i cammini dei genitori
            int CUT_POINT = (int) cut_point;
            auto [child1, child2] = crossover(new_population[2*j], new_population[2*j + 1], CUT_POINT);
            new_population[2*j] = child1;
            new_population[2*j + 1] = child2;
         }
         
         double M1 = rnd.Rannyu();
         if (M1<p_m1){
            double a = rnd.Rannyu(1, 33);
            int A = (int) a;
            double b = rnd.Rannyu(1, 33);
            int B = (int) b;
            new_population[2*j] = mutation_1(new_population[2*j], A, B);
            new_population[2*j + 1] = mutation_1(new_population[2*j + 1], A, B);
         }
         
         double M2 = rnd.Rannyu();
         if (M2<p_m2){
            double start = rnd.Rannyu(1, 10);
            int START = (int) start;
            double m = rnd.Rannyu(1, 20);
            int M = (int) m;
            double n = rnd.Rannyu(1, 5);
            int N = (int) n;
            new_population[2*j] = mutation_2(new_population[2*j], START, M, N);
            new_population[2*j + 1] = mutation_2(new_population[2*j + 1], START, M, N);
         }
         
         double M3 = rnd.Rannyu();
         if (M3<p_m3){
            double start1 = rnd.Rannyu(1, 10);
            int START1 = (int) start1;
            double start2 = rnd.Rannyu(21, 25);
            int START2 = (int) start2;
            double m = rnd.Rannyu(1, 10);
            int M = (int) m;
            new_population[2*j] = mutation_3(new_population[2*j], START1, START2, M);
            new_population[2*j + 1] = mutation_3(new_population[2*j + 1], START1, START2, M);
         }
         
         double M4 = rnd.Rannyu();
         if (M4<p_m4){
            double start = rnd.Rannyu(1, 25);
            int START = (int) start;
            double m = rnd.Rannyu(1, 10);
            int M = (int) m;
            new_population[2*j] = mutation_4(new_population[2*j], START, M);
            new_population[2*j + 1] = mutation_4(new_population[2*j + 1], START, M);
         }
      }
      
      for (const auto& vec : new_population) {   // controllo che ogni individuo della nuova popolazione rispetti i vincoli
         if (!check(vec)) cout << "Il vettore di interi non soddisfa tutte le condizioni." << endl;
      }
      
      population = new_population;
      
   }
   vector<double> Loss(dim);
   for (int j=0; j<dim; j++){   // calcolo la loss function per l'ultima popolazione
      vector<double> X(35);
      vector<double> Y(35);
      for (int k=0; k<35; k++){
         X[k] = cities[population[j][k]].getX();
         Y[k] = cities[population[j][k]].getY();
      }
      Loss[j]=loss(X, Y);
   }  
   
   sortPopulationByLoss(Loss, population);   // ordino l'ultima popolazione

   // Stampa della struttura population
   for (const auto& vec : population) {
        for (int index : vec) {
            cout << index << " ";
        }
        cout << (check(vec) ? " - Valid" : " - Invalid") << endl;
    }

   ofstream output;
   output.open("Es09_01.csv");   // file di output contenente le coordinate del best path
   if (output.is_open()){
      output << "X" << "," << "Y" << endl;
      for(int i=0; i<35; i++){
         for(int j=0; j<33; j++){
            if(population[0][i]==cities[j].getIndex()){
               output << cities[j+1].getX() << "," << cities[j+1].getY() << endl;
            }
         }
      }
   } else cerr << "PROBLEM: Unable to open Es09_01.csv" << endl;
   output.close();
  
   ofstream Output;
   Output.open("Es09_Loss_Best.csv");   // file di output contenente la loss del best path per ogni generazione
   if (Output.is_open()){
      Output << "Generazione" << "," << "BestLoss" << endl;
      for(int i=0; i<120; i++)   Output << i+1 << "," << LOSS_BEST[i] << endl;
   } else cerr << "PROBLEM: Unable to open Es09_Loss_Best.csv" << endl;
   Output.close();

   ofstream OutPut;
   OutPut.open("Es09_Loss_Half.csv");   // file di output contenente la loss, mediata sulla migliore metà della popolazione, per ogni generazione
   if (OutPut.is_open()){
      OutPut << "Generazione" << "," << "HalfLoss" << endl;
      for(int i=0; i<120; i++)   OutPut << i+1 << "," << LOSS_HALF[i] << endl;
   } else cerr << "PROBLEM: Unable to open Es09_Loss_Half.csv" << endl;
   OutPut.close();

   rnd.SaveSeed();
   return 0;
}


