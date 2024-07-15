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
    for (size_t i = vec.size() - 2; i > 0; --i) {
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

    // Verifica che il primo e l'ultimo elemento siano uguali
    if (vec.front() != vec.back()) {
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
double loss(std::vector<double> x, std::vector<double> y){
  double sum = 0;
  for (int i=0; i<34; i++)   sum += (x[i]-x[i+1])*(x[i]-x[i+1])+(y[i]-y[i+1])*(y[i]-y[i+1]);
  return sum;
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

// Mutazione 2 -> Shift di n posizioni per m città contigue, a partire da start
std::vector<int> mutation_2(std::vector<int> path, int start, int m, int n){
  std::vector<int> mutated_path = path;

  // Verifica che m sia appropriato e start + m - 1 non ecceda i limiti
    if (m > 0 && m < 35 && start >= 1 && start + m - 1 < 35) {
      // Calcola la posizione finale dopo lo spostamento
      int shiftStart = start;
      int shiftEnd = start + m;

      // Normalizza n per evitare operazioni inutili
      n = n % m;
      if (n == 0) return mutated_path; // Nessuno shift necessario se n è multiplo di m

      // Per eseguire lo shift
      rotate(mutated_path.begin() + shiftStart, mutated_path.begin() + shiftEnd - n, mutated_path.begin() + shiftEnd);
    } else {
      cerr << "Mutazione 2: Parametri non validi: m deve essere < 34, start e m devono essere validi per il vettore di dimensione 35." << endl;
    }
    
    return mutated_path;
}

// Mutazione 3 -> Permutazione di m città contigue (a partire da start1) con altre m città contigue (a partire da start2)
std::vector<int> mutation_3(std::vector<int> path, int start1, int start2, int m) {
    std::vector<int> mutated_path = path;
    int N = path.size(); // Determina la dimensione del vettore

    // Verifica che i parametri siano validi
    if (m > 0 && m < N/2 &&
        start1 > 0 && start1 + m <= N && 
        start2 > 0 && start2 + m <= N && 
        (start1 + m <= start2 || start2 + m <= start1)) { // Assicura che i segmenti non si sovrappongano
        // Scambia gli elementi dei due segmenti
        for (int i = 0; i < m; ++i) {
            std::swap(mutated_path[start1 + i], mutated_path[start2 + i]);
        }
    } else {
        std::cerr << "Mutazione: Parametri non validi: m deve essere < N/2, start1 e start2 devono essere validi e non sovrapporsi." << std::endl;
    }
    
    return mutated_path;
}

// Mutazione 4 -> Inversione dell'ordine in cui appaiono nel cammino di m città, a partire da start
std::vector<int> mutation_4(std::vector<int> path, int start, int m) {
    std::vector<int> mutated_path = path;
    int N = path.size();

    // Verifica che i parametri siano validi
    if (m > 0 && m < N && start > 0 && start + m <= N) {
        // Per invertire il segmento
        std::reverse(mutated_path.begin() + start, mutated_path.begin() + start + m);
    } else {
        std::cerr << "Mutazione 4: Parametri non validi: m deve essere < N, start deve essere valido per il vettore di dimensione " << N << "." << std::endl;
    }
    
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


