#ifndef VOLTAGEPIN_H
#define VOLTAGEPIN_H

#include <Arduino.h>
#include "Parser.h"

class VoltagePin {
  
public:
  
  VoltagePin(int pin_num);
    
  virtual bool getState();
  virtual bool getStatus();
  bool getInversion();
    
  void toggle();
  void setState(bool state);
  void setInversion(bool invert=1);
  
  void parseState(char *arg);
    
protected:
  
  int _pin_num;
  bool _state;
  bool _invert;
  
};

#endif