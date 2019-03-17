#include <OutputPinArray.h>

OutputPinArray::OutputPinArray(){
  
  initialize();
  
}

OutputPinArray::~OutputPinArray(){
  
}

void OutputPinArray::initialize(){
 
  _num_pins=0;
  update();
  
}

char *OutputPinArray::printout(char *buf){
  
  return buf;
  
}

OutputPin &OutputPinArray::getPin(int num){
 
  num=num%_num_pins;
  
  return _pins[num];
  
}

OutputPin &OutputPinArray::operator[](int num){
 
  return getPin(num);
  
}

int OutputPinArray::getNumPins() const {
  
  return _num_pins;
  
}

void OutputPinArray::update(){
 
  for(int i=0;i<_num_pins;i++){
   
    _pins[i].update();
    
  }
  
}

void OutputPinArray::setMode(int mode){
   
  for(int i=0;i<_num_pins;i++){
   
    _pins[i].setMode(mode);
    
  }
  
}

void OutputPinArray::setTimerMode(int mode){
  
  for(int i=0;i<_num_pins;i++){
   
    _pins[i].timer.setMode(mode);
    
  }
  
}

void OutputPinArray::setInterval(float time){
  
  for(int i=0;i<_num_pins;i++){
   
    _pins[i].timer.setInterval(time);
    
  }
  
}

void OutputPinArray::setPulseLength(float time){
  
  for(int i=0;i<_num_pins;i++){
   
    _pins[i].timer.setPulseLength(time);
    
  }
  
}

void OutputPinArray::setPinNumbers(int array[], int size){
 
  _num_pins=size;
  
  for(int i=0;i<_num_pins;i++){
   
    _pins[i].setPinNum(array[i]);
    _pins[i].initialize();
    _pins[i].update();
    
  }
  
}
