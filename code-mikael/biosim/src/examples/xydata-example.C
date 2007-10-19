#include "../xydata.h"
#include<iostream>

using namespace std;

int main() {
  xydata xy(-5,5,0.2);
  xy.add( -0.3, 100.23 );
  xy.list();
};
