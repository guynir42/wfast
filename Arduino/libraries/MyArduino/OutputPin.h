#ifndef OUTPUTPIN_H
#define OUTPUTPIN_H

#include <Arduino.h>
#include <Timer.h>
#include <Parser.h>

class OutputPin {
  
public:
  
  OutputPin(int pin=0);
  ~OutputPin();
  
  void initialize();
  char *printout();
  
  Timer timer;
  
  // getters
  int getPinNum() const;
  int getMode() const;
  char *getModeString() const;
  int getStatus();
  int getState();
  char *getOutBuffer();
  bool getInversionState() const;
  
  // setters
  void setPinNum(int number);
  void setMode(int mode);
  void setInverse();
  void setInverse(bool state);
  void update();
  
  // serial commands
  void parse(char *arg);
  
  enum Mode {OFF=0, ON, WATCH}; // removed UNWATCH because there's invert now...
  
  static bool debug_bit;
  
protected:
  
  int _pin_num;
  bool _inverse;
  // bool _state;
  
  int _current_mode;
  char _out_buf[Parser::MAXN];
  
};

#endif