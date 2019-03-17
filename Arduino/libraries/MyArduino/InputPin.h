#ifndef INPUTPIN_H
#define INPUTPIN_H

#include <Arduino.h>
#include "Parser.h"

class InputPin {

public:
  
  InputPin(int pin_number=0);
  ~InputPin();
  
  void initialize();
  char *printout();
  
  // getters
  int getPinNumber();
  bool getState();
  bool getPinState();
  int getMode();
  char *getModeString();
  char *getOutBuffer();
  int getValue();
  int getPushNumber();
  int getPushDelay();
  bool hasChanged();
  
  // setters
  void setMode(int mode);
  void setPushNumber(int num);
  void setPinState(bool state=1);
  void setThreshold(int thresh);
  void update(); // in PUSH mode, you need to either getPinState or update frequently
  
  enum Inputs{DIGITAL=0, PUSH, SWITCH, ANALOG, ANTI_ANALOG, INTEGRATOR};
  
  static const int STRN=20;
  
  static bool debug_bit;
  
protected:
  
  int _pin_number;
  
  char _out_buffer[STRN];
  
  int _current_mode;
  
  bool _last_input;
  int _push_num;
  static const int _push_delay=1; // millisecs
  long _last_push_time; // millisecs
  bool _has_changed;
  
  int _threshold; // for ANALOG and ANTI_ANALOG. when above (or below for anti) this value, getState() returns true
  int _integrator_length;
  int _integrator_sum;
  int _integrator_index;
  int _integrator_result;
  
  
};

#endif