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
  //void setX(double newX);
  //void setY(double newY);
  //void setIndex(int newIndex);
  
  // Metodo per stampare le informazioni della città
  void print() const;

  // Function for statistical uncertainty estimation
  //double error(std::vector<double>, std::vector<double>, int);
};

void Shuffle(std::vector<int>& vec, double g);

bool check(const std::vector<int>& vec);   // Function to verify bonds

double loss(std::vector<double>, std::vector<double>);

void sortPopulationByLoss(std::vector<double>&, std::vector<std::vector<int>>&);

int selection(double, int);

std::vector<int> mutation_1(std::vector<int>, int, int);

std::vector<int> mutation_2(std::vector<int>, int, int, int);

std::vector<int> mutation_3(std::vector<int>, int, int, int);

std::vector<int> mutation_4(std::vector<int>, int, int);  

std::pair<std::vector<int>, std::vector<int>> crossover(const std::vector<int> mother, const std::vector<int> father, int cut_point);

#endif // __City__

