#include "GroundPin.h"

GroundPin::GroundPin(int pin_num){
 
  _pin_num=pin_num;
  
  pinMode(_pin_num, OUTPUT);
  digitalWrite(_pin_num, LOW);
  
}