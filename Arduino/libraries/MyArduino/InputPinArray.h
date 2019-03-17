#ifndef INPUTPINARRAY_H
#define INPUTPINARRAY_H

#include <Arduino.h>
#include <InputPin.h>
#include <Parser.h>

class InputPinArray {
  
public:
  
  InputPinArray();
  ~InputPinArray();
  
  void initialize();
  char *printout();
  
  // getters
  InputPin &getPin(int num);
  InputPin &operator[](int num);
  bool getState(int num);
  bool getStateAll();
  bool getStateAny();
  int getNumPins() const;
  char *getOutBuffer();
  
  // setters
  void update();
  void setMode(int mode);
  
  #ifndef PIN_ARR_SIZE
    static const int MAX_NUMBER=10;
  #else
    static const int MAX_NUMBER=PIN_ARR_SIZE;
  #endif
    
protected:
    
  char _out_buffer[Parser::MAXN];  
  
  int _num_pins;
  InputPin _pins[MAX_NUMBER];
  
  
    
};

#endif