#ifndef GROUNDPIN_H
#define GROUNDPIN_H

#include <Arduino.h>

class GroundPin {
  
public:
  
  GroundPin(int pin_num);
    
protected:
  
  int _pin_num;
  
};

#endif