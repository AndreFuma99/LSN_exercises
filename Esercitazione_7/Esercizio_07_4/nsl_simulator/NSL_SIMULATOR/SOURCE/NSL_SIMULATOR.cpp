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
  for(int i=0; i < 17000; i++){
    SYS.step();
    SYS.measure();
  }

  // Faccio partire la simulazione
  for(int i = 0; i < SYS.get_nbl(); i++){ //loop over blocks
    SYS.block_reset(i+1);
    for(int j = 0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
      SYS.step();
      SYS.measure();
        if(j % 10 == 0){
       // SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
          nconf++;
        }
    }
    SYS.averages(i+1);
    //SYS.block_reset(i+1);
  }
  SYS.finalize();

  return 0;
}


