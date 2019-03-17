#include <InputPinArray.h>

InputPinArray::InputPinArray(){
  
  initialize();
  
}

InputPinArray::~InputPinArray(){
  
}

void InputPinArray::initialize(){
  
  _num_pins=0;
  update();
  
}

char *InputPinArray::printout(){
 
  snprintf(_out_buffer, Parser::MAXN, "IN-PINS: ");
  
  for(int i=0;i<_num_pins;i++){
   
    snprintf(_out_buffer, Parser::MAXN, "%s%d:%d", _out_buffer, i, getState(i));
    
  }
  
  return _out_buffer;
  
}

// getters
InputPin &InputPinArray::getPin(int num){
 
  if(num<0 || num>=_num_pins){
   
    snprintf(_out_buffer, Parser::MAXN, "cannot access pin %d (max: %d)", num, _num_pins);
    return _pins[0];
  }
  
  return _pins[num];
    
}

InputPin &InputPinArray::operator[](int num){
  
  return getPin(num);
  
}

bool InputPinArray::getState(int num){
  
  return getPin(num).getPinState();
  
}

bool InputPinArray::getStateAll(){
 
  if(_num_pins<1) return 0;
  
  bool result=1;
  
  for(int i=0;i<_num_pins;i++){
   
    if(!_pins[i].getPinState()) result=0;
    
  }
  
  return result; 
  
}

bool InputPinArray::getStateAny(){
  
  if(_num_pins<1) return 0;
  
  bool result=0;
  
  for(int i=0;i<_num_pins;i++){
   
    if(_pins[i].getPinState()) result=1;
    
  }
  
  return result; 
    
}

int InputPinArray::getNumPins() const { 
  
  return _num_pins;
  
}

char *InputPinArray::getOutBuffer(){
  
  return _out_buffer;
  
}

// setters
void InputPinArray::update(){
 
  for(int i=0;i<_num_pins;i++){
   
    _pins[i].update();
    
  }
  
}

void InputPinArray::setMode(int mode){
 
  for(int i=0;i<_num_pins;i++){
   
    _pins[i].setMode(mode);
    
  }
  
}