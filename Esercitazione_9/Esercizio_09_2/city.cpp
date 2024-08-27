#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "city.h"
#include "random.h"
#include <unordered_set>
#include <algorithm>

using namespace std;

City::City(double x, double y, int index) : x(x), y(y), index(index) {}
// Default constructor, does not perform any action

City::~City(){}
// Default destructor, does not perform any action

double City::getX() const { return x; }
double City::getY() const { return y; }
int City::getIndex() const { return index; }

// Metodo per stampare le informazioni della città
void City::print() const {
    cout << "City Index: " << index << ", X: " << x << ", Y: " << y << endl;
}

// Funzione che mescola gli indici delle città per ottenere una sequenza casuale
void Shuffle(std::vector<int>& vec, double g) {
    for (size_t i = vec.size() - 1; i > 0; --i) {
        int j = 1 + static_cast<int>((i) * g);
        swap(vec[i], vec[j]);
    }
}

// Check function che verifica se ogni individuo rispetta i vincoli
bool check(const std::vector<int>& vec) {
    // Verifica che ci siano almeno due elementi nel vettore
    if (vec.size() < 2) {
        return false;
    }

    // Verifica che il primo elemento sia uguale a 1
    if (vec.front() != 1) {
        return false;
    }

    // Verifica che ogni numero intero si presenti una e una sola volta
    std::unordered_set<int> uniqueNumbers;
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i == 0 || i == vec.size() - 1) {
            continue; // Salta il primo e l'ultimo elemento
        }
        // Inserisce l'elemento nel set e verifica se è già presente
        if (!uniqueNumbers.insert(vec[i]).second) {
            return false; // Numero duplicato trovato
        }
    }

    return true; // Tutte le condizioni sono soddisfatte
}

// Loss function
double loss(std::vector<double> x, std::vector<double> y, int n_city){
  double sum = 0;
  for (int i=0; i<(n_city-1); i++)   sum += (x[i]-x[i+1])*(x[i]-x[i+1])+(y[i]-y[i+1])*(y[i]-y[i+1]);
  sum += (x[n_city-1]-x[0])*(x[n_city-1]-x[0])+(y[n_city-1]-y[0])*(y[n_city-1]-y[0]);   // valuto anche il cammino dall'ultima città alla prima
  return sum;
}

// Periodic Boundary Conditions -> per fare agire le mutazioni sull'intero individuo
int pbc(int a, int n_city){
	if(a<n_city){
	return a;
	}
	else
	 return a%n_city + 1;
}

// Funzione che ordina la popolazione in base alla loss function crescente
void sortPopulationByLoss(std::vector<double>& Loss, std::vector<std::vector<int>>& population) {
    size_t n = Loss.size();

    // Crea un vettore di indici
    std::vector<size_t> indices(n);
    for (size_t l = 0; l < n; ++l) {
        indices[l] = l;
    }

    // Ordina gli indici in base ai valori di Loss
    std::sort(indices.begin(), indices.end(), [&Loss](size_t l1, size_t l2) {
        return Loss[l1] < Loss[l2];
    });

    // Crea copie ordinate di Loss e population
    std::vector<double> sorted_Loss(n);
    std::vector<std::vector<int>> sorted_population(n);
    for (size_t l = 0; l < n; ++l) {
        sorted_Loss[l] = Loss[indices[l]];
        sorted_population[l] = population[indices[l]];
    }

    // Sostituisci Loss e population con le versioni ordinate
    Loss = std::move(sorted_Loss);
    population = std::move(sorted_population);
}

// Funzione di scambio
void exchange(std::vector<int>&v, int i, int j) {

	int dim = v.size();
	if(i < dim and j < dim) {
	
		int appo;
		appo = v[i];
		v[i] = v[j];
		v[j] = appo;
		
	}
}

// Operatore di selezione
int selection(double r, int dim){
  double p = 3.;
  float index = dim*pow(r,p);
  int intero = (int) index;
  return intero;
}

// Mutazione 1 -> Permutazione di una coppia di città
std::vector<int> mutation_1(std::vector<int> path, int A, int B){
  std::vector<int> mutated_path = path;
  
  int temp = path[A];
  mutated_path[A] = path[B];
  mutated_path[B] = temp;
  
  return mutated_path;
}

// Mutazione 3 -> Permutazione di m città contigue (a partire da a) con altre m città contigue (a partire da b)
std::vector<int> mutation_3(std::vector<int> path, int m, int a, int b, int n_city) {
    std::vector<int> mutated_path = path; 
		
	for(int i=0; i<m; i++)   exchange(mutated_path, pbc(a+i, n_city), pbc(b+i, n_city));
    
    return mutated_path;
}

// Mutazione 4 -> Inversione dell'ordine in cui appaiono nel cammino di m città, a partire da a
std::vector<int> mutation_4(std::vector<int> path, int m, int a, int n_city) {
    std::vector<int> mutated_path = path;

    for(int i=0; i<m/2; i++)   exchange(mutated_path, pbc(a+i, n_city), pbc(a+m-1-i, n_city));
    
    return mutated_path;
}

// Operatore di crossover
std::pair<std::vector<int>, std::vector<int>> crossover(const std::vector<int> mother, const std::vector<int> father, int cut_point) {
    int parentSize = mother.size(); // stessa lunghezza di mother e father
    
    // Creare i nuovi percorsi con la prima parte conservata
    std::vector<int> child1(mother.begin(), mother.begin() + cut_point);
    std::vector<int> child2(father.begin(), father.begin() + cut_point);
    
    // Funzione per completare il percorso con le città mancanti
    auto complete_path = [](std::vector<int>& child, const std::vector<int>& partner) {
        for (int city : partner) {
            if (std::find(child.begin(), child.end(), city) == child.end()) {
                child.push_back(city);
            }
        }
    };
    
    // Completare i percorsi
    complete_path(child1, father);
    complete_path(child2, mother);
    
    // Aggiungere l'ultimo elemento dei genitori ai figli
    child1.push_back(mother.back());
    child2.push_back(father.back());
    
    // Se la lunghezza dei figli è maggiore di quella dei genitori, rimuovi un elemento dalla fine
    if (child1.size() > parentSize) {
        child1.pop_back();
    }
    if (child2.size() > parentSize) {
        child2.pop_back();
    }
    
    return {child1, child2};
}


