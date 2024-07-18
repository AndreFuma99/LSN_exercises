#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "random.h"
#include "city.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include <mpi.h>

using namespace std;
 
int main (int argc, char *argv[]){

   MPI_Init(&argc, &argv);   // per inizializzare l'ambiente MPI

   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);   // determina il rank del processo corrente all'interno del comunicatore specificato
   MPI_Comm_size(MPI_COMM_WORLD, &size);   // determina il numero totale di processi all'interno del comunicatore specificato
   MPI_Status stat1, stat2;

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      for(int i=0; i<=rank; i++) {   // importo primes diversi per ogni rank
         Primes >> p1 >> p2 ;
      }
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


   
   vector<City> cities;
   int n_city = 110;
   cities.reserve(n_city);

   double x, y;
   ifstream Caprov("cap_prov_ita.dat");
   if (Caprov.is_open()){
      for (int i=0; i<n_city; i++){
         Caprov >> x >> y;
         cities.emplace_back(x, y, i+1);
      }
   } else cerr << "PROBLEM: Unable to open cap_prov_ita.dat" << endl;
   Caprov.close();

   // Creo la popolazione iniziale
   
   int dim = 1000;
   vector<vector<int>> population(dim, vector<int>(n_city));
   
   // Genero una sequenza di indici delle città da 1 a 110
   vector<int> cityIndices(n_city);
   iota(cityIndices.begin(), cityIndices.end(), 1); // Riempie cityIndices con {1, 2, ..., 110}
   
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
   
   double p_m1 = 0.1;   // probabilità della permutazione di una coppia
   double p_m3 = 0.1;   // probabilità della permutazione di due gruppi
   double p_m4 = 0.1;   // probabilità dell'inversione dell'ordine
   double p_c = 0.9;   // probabilità del crossover
   
   int N_migr = 3; // Intervallo di generazioni tra le migrazioni
   
   vector <int> imesg(dim);   // vettore d'appoggio per la migrazione
   vector <int> imesg2(dim);  // vettore d'appoggio per la migrazione
   int itag=1; int itag2=2;
   vector <int> which_swap = {0, 1, 2, 3};   // vettore contenente l'indice di ogni rank
   
   int gen = 750;   // numero di generazioni
   vector<double> LOSS_BEST(gen);   // vettore che contiene la loss function del miglior cammino per ogni generazione
   vector<double> LOSS_HALF(gen);   // vettore che contiene la loss function, mediata sulla migliore metà della popolazione, per ogni generazione
   
   for (int i=0; i<gen; i++){
   
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
      
      for (int j=0; j<(dim/2); j++){
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
      
      // Migrazione ogni N_migr generazioni
      if (i % N_migr == 0){
      
         cout << "Generation: " << i << endl;
         
         random_shuffle(which_swap.begin(), which_swap.end());   // per mischiare gli indici del vettore which_swap
    
         for (int j = 0; j < n_city; j++) {
            imesg[j] = population[0][j]; // Prendo il miglior individuo della popolazione locale
            imesg2[j] = population[0][j]; // Prendo il miglior individuo della popolazione locale
         }

         // Eseguo lo scambio tra i processi
         if (rank == which_swap[1]) {
         
               MPI_Send(&imesg[0], n_city, MPI_INTEGER, which_swap[0], itag, MPI_COMM_WORLD);   // ws[1] invia imesg a ws[0]
               MPI_Recv(&imesg2[0], n_city, MPI_INTEGER, which_swap[0], itag2, MPI_COMM_WORLD, &stat2);   // ws[1] riceve imesg2 da ws[0]
               
         } else if (rank == which_swap[0]) {
            
               MPI_Send(&imesg2[0], n_city, MPI_INTEGER, which_swap[1], itag2, MPI_COMM_WORLD);   // ws[0] invia imesg2 a ws[1]
               MPI_Recv(&imesg[0], n_city, MPI_INTEGER, which_swap[1], itag, MPI_COMM_WORLD, &stat1);   // ws[0] riceve imesg da ws[1]
            
         }
         
         // Sostituisco gli individui locali con quelli ricevuti
         if (rank == which_swap[1]) {
            for (int j = 0; j < n_city; j++)   population[0][j] = imesg2[j];
         } else if (rank == which_swap[0]) {
            for (int j = 0; j < n_city; j++)   population[0][j] = imesg[j];
         }
         
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

      }
      
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

   // Stampo i risultati localmente per ogni processo
    string filename1 = "Es10_rank_" + std::to_string(rank) + ".csv";
    ofstream output(filename1);
    if (output.is_open()) {
        output << "X,Y" << endl;
        for (int i = 0; i < n_city; ++i) {
            for (int j = 0; j < n_city-1; ++j) {
                    if (population[0][i] == cities[j].getIndex()) {
                        output << cities[j + 1].getX() << "," << cities[j + 1].getY() << endl;
                    }
            }
        }
    } else {
        cerr << "PROBLEM: Unable to open " << filename1 << endl;
    }
    output.close();

    string filename2 = "Es10_Loss_Best_rank_" + std::to_string(rank) + ".csv";
    ofstream Output(filename2);
    if (Output.is_open()) {
        Output << "Generazione,BestLoss" << endl;
        for (int i = 0; i < gen; ++i) {
            Output << i + 1 << "," << LOSS_BEST[i] << endl;
        }
    } else {
        cerr << "PROBLEM: Unable to open " << filename2 << endl;
    }
    Output.close();

    string filename3 = "Es10_Loss_Half_rank_" + std::to_string(rank) + ".csv";
    ofstream OutPut(filename3);
    if (OutPut.is_open()) {
        OutPut << "Generazione,HalfLoss" << endl;
        for (int i = 0; i < gen; ++i) {
            OutPut << i + 1 << "," << LOSS_HALF[i] << endl;
        }
    } else {
        cerr << "PROBLEM: Unable to open " << filename3 << endl;
    }
    OutPut.close();


    rnd.SaveSeed();
   
    MPI_Finalize();   // per finalizzare l'ambiente MPI
   
    return 0;
}


