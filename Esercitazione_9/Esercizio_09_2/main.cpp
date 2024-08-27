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
   int n_city = 34;
   cities.reserve(n_city);
   
   for (int i=0; i<n_city; i++){
      double x = rnd.Rannyu(-1,1);
      double y = rnd.Rannyu(-1,1);
      cities.emplace_back(x, y, i+1);
   }
   
   // Creo la popolazione iniziale -> genero randomicamente dim percorsi, ossia dim vector che contengono gli indici delle città
   
   int dim = 900;
   vector<vector<int>> population(dim, vector<int>(n_city));
   
   // Genero una sequenza di indici delle città da 1 a 34
   vector<int> cityIndices(n_city);
   iota(cityIndices.begin(), cityIndices.end(), 1); // Riempie cityIndices con {1, 2, ..., 34}
   
   // Definizione dell'indice comune per la prima posizione di ogni vettore
   int commonIndex = 1;
   
   for (int i = 0; i < dim; ++i) {
      // Copio gli indici delle città
      vector<int> indices = cityIndices;

      // Mescolo gli indici per ottenere una sequenza casuale
      double g = rnd.Rannyu();
      Shuffle(indices, g);

      // Imposto il primo elemento al valore comune
      indices[0] = commonIndex;

      // Assegno la sequenza al vettore i-esimo di population
      population[i] = indices;
   }
   
   double p_m1 = 0.08;   // probabilità della permutazione di una coppia
   double p_m3 = 0.07;   // probabilità della permutazione di due gruppi
   double p_m4 = 0.09;   // probabilità dell'inversione dell'ordine
   double p_c = 0.8;   // probabilità del crossover
   
   int gen = 120;   // numero di generazioni
   vector<double> LOSS_BEST(gen);   // vettore che contiene la loss function del miglior cammino per ogni generazione
   vector<double> LOSS_HALF(gen);   // vettore che contiene la loss function, mediata sulla migliore metà della popolazione, per ogni generazione
   
   for (int i=0; i<gen; i++){   // Faccio 120 generazioni perchè ho visto che dopo circa 85 la Loss del miglior cammino si stabilizza
   
      vector<double> Loss(dim);
      for (int j=0; j<dim; j++){   // calcolo la loss function per ogni cammino
         vector<double> X(n_city);
         vector<double> Y(n_city);
         for (int k=0; k<n_city; k++){
            X[k] = cities[population[j][k]].getX();
            Y[k] = cities[population[j][k]].getY();
         }
         Loss[j]=loss(X, Y, n_city);
      }
      
      sortPopulationByLoss(Loss, population);
      
      LOSS_BEST[i] = Loss[0];
      double sum = 0;
      for (int q=0; q<(dim/2); q++)   sum += Loss[q];
      LOSS_HALF[i] = sum/(dim/2.);
      
      vector<vector<int>> new_population(dim, vector<int>(n_city));
      
      for (int j=0; j<(dim/2); j++){   // Applico la selezione dim/2 volte per ottenere dim figli 
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
            double cut_point = rnd.Rannyu(1, n_city-1);   // genero la posizione in cui tagliare i cammini dei genitori
            int CUT_POINT = (int) cut_point;
            auto [child1, child2] = crossover(new_population[2*j], new_population[2*j + 1], CUT_POINT);
            new_population[2*j] = child1;
            new_population[2*j + 1] = child2;
         }
         
         double M1 = rnd.Rannyu();
         if (M1<p_m1){
            double a = rnd.Rannyu(1, n_city);
            int A = (int) a;
            double b = rnd.Rannyu(1, n_city);
            int B = (int) b;
            new_population[2*j] = mutation_1(new_population[2*j], A, B);
            new_population[2*j + 1] = mutation_1(new_population[2*j + 1], A, B);
         }
         
         double M3 = rnd.Rannyu();
         if (M3<p_m3){
            double m = rnd.Rannyu(2, (n_city/2)-1);
            int M = (int) m;
            double a = rnd.Rannyu(1, n_city);
            int A = (int) a;
            double b = rnd.Rannyu(1, n_city);
            int B = (int) b;
            new_population[2*j] = mutation_3(new_population[2*j], M, A, B, n_city);
            new_population[2*j + 1] = mutation_3(new_population[2*j + 1], M, A, B, n_city);
         }
         
         double M4 = rnd.Rannyu();
         if (M4<p_m4){
            double m = rnd.Rannyu(1, n_city);
            int M = (int) m;
            double a = rnd.Rannyu(1, n_city);
            int A = (int) a;
            new_population[2*j] = mutation_4(new_population[2*j], M, A, n_city);
            new_population[2*j + 1] = mutation_4(new_population[2*j + 1], M, A, n_city);
         }
      }
      
      for (const auto& vec : new_population) {   // controllo che ogni individuo della nuova popolazione rispetti i vincoli
         if (!check(vec)) cout << "Il vettore di interi non soddisfa tutte le condizioni." << endl;
      }
      
      population = new_population;
      
   }
   
   vector<double> Loss(dim);
   for (int j=0; j<dim; j++){   // calcolo la loss function per l'ultima popolazione
      vector<double> X(n_city);
      vector<double> Y(n_city);
      for (int k=0; k<n_city; k++){
         X[k] = cities[population[j][k]].getX();
         Y[k] = cities[population[j][k]].getY();
      }
      Loss[j]=loss(X, Y, n_city);
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
   output.open("Es09_02.csv");   // file di output contenente le coordinate del best path
   if (output.is_open()){
      output << "X" << "," << "Y" << endl;
      for(int i=0; i<n_city; i++){
         for(int j=0; j<n_city-1; j++){
            if(population[0][i]==cities[j].getIndex()){
               output << cities[j+1].getX() << "," << cities[j+1].getY() << endl;
            }
         }
      }
   } else cerr << "PROBLEM: Unable to open Es09_02.csv" << endl;
   output.close();
  
   ofstream Output;
   Output.open("Es09_Loss_Best.csv");   // file di output contenente la loss del best path per ogni generazione
   if (Output.is_open()){
      Output << "Generazione" << "," << "BestLoss" << endl;
      for(int i=0; i<gen; i++)   Output << i+1 << "," << LOSS_BEST[i] << endl;
   } else cerr << "PROBLEM: Unable to open Es09_Loss_Best.csv" << endl;
   Output.close();

   ofstream OutPut;
   OutPut.open("Es09_Loss_Half.csv");   // file di output contenente la loss, mediata sulla migliore metà della popolazione, per ogni generazione
   if (OutPut.is_open()){
      OutPut << "Generazione" << "," << "HalfLoss" << endl;
      for(int i=0; i<gen; i++)   OutPut << i+1 << "," << LOSS_HALF[i] << endl;
   } else cerr << "PROBLEM: Unable to open Es09_Loss_Half.csv" << endl;
   OutPut.close();

   rnd.SaveSeed();
   return 0;
}


