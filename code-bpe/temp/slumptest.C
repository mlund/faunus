#include <iostream>
#include "classes/slump.h"

using namespace std;

int main () {

  float test1=0.5;
  int test=0;
  slump2 sl;
  for (int i=0; i<1.e6; i++) {
    if (sl.random_one()<test1)
    test++;
  }
  cout <<"Trial fraction = "<< (test1*100) <<"%"<<endl;
  cout <<"The generated fraction = " << test/(1.e6)*100<<"%"<<endl;
  return 0;
} 
