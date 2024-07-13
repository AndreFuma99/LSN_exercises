#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  int nconf = 1;
  System SYS;
  SYS.initialize();
  SYS.initialize_properties();
  SYS.block_reset(0);
  
  // Fase di equilibrazione
  for(int i=0; i < 1000; i++){
    SYS.step();
    SYS.measure();
  }
  
  // Faccio partire la simulazione
  for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
    SYS.step();
    SYS.measure();
    
    SYS.averages(i+1);
    SYS.block_reset(i+1);
  }
  SYS.finalize();

  return 0;
}


