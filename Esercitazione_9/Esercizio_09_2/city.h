#include <vector>

#ifndef __City__
#define __City__

class City {

private:
  double x;
  double y;
  int index; // Indice della città

public:
  // Default constructor
  City(double x, double y, int index);
  // Destructor
  ~City();
  
  double getX() const;
  double getY() const;
  int getIndex() const;
  
  // Metodo per stampare le informazioni della città
  void print() const;

};

void Shuffle(std::vector<int>& vec, double g);

bool check(const std::vector<int>& vec);

double loss(std::vector<double>, std::vector<double>, int);

int pbc(int, int);

void sortPopulationByLoss(std::vector<double>&, std::vector<std::vector<int>>&);

void exchange(std::vector<int>&, int, int);

int selection(double, int);

std::vector<int> mutation_1(std::vector<int>, int, int);

std::vector<int> mutation_3(std::vector<int>, int, int, int, int);

std::vector<int> mutation_4(std::vector<int>, int, int, int);

std::pair<std::vector<int>, std::vector<int>> crossover(const std::vector<int> mother, const std::vector<int> father, int cut_point);

#endif // __City__


