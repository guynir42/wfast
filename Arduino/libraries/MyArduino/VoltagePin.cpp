#include "VoltagePin.h"

VoltagePin::VoltagePin(int pin_num){
 
  _pin_num=pin_num;
  
  pinMode(_pin_num, OUTPUT);
  digitalWrite(_pin_num, HIGH);
  
  _state=1;
  _invert=0;
  
}

bool VoltagePin::getState(){
  
  return _state;
  
}

bool VoltagePin::getStatus(){
  
  return getState();
  
}

bool VoltagePin::getInversion(){
  
  return _invert;
  
}

void VoltagePin::toggle(){
  
  setState(!getState());
  
}

void VoltagePin::setState(bool state){
  
  _state=state;
  if(_invert) digitalWrite(_pin_num, !state);
  else digitalWrite(_pin_num, state);
  
}

void VoltagePin::setInversion(bool invert){
 
  _invert=invert;
  _state=!_state; // physical state doesn't change. only the way we label it is reversed! 
  
}

void VoltagePin::parseState(char *arg){
 
  Parser::cleanString(arg);
  
  if(Parser::partialMatch(arg, "on") || Parser::partialMatch(arg, "1") || Parser::partialMatch(arg, "high")) setState(1);
  if(Parser::partialMatch(arg, "off") || Parser::partialMatch(arg, "0") || Parser::partialMatch(arg, "low")) setState(0);
  
}