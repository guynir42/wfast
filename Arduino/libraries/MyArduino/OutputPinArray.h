#ifndef OUTPUTPINARRAY_H
#define OUTPUTPINARRAY_H

#include <OutputPin.h>


class OutputPinArray {
  
public:

  OutputPinArray();
  ~OutputPinArray();
  
  void initialize();
  char *printout(char *buf);
  
  // getters
  OutputPin &getPin(int num);
  OutputPin &operator[](int num);
  int getNumPins() const;  
  
  // setters
  void update();
  void setMode(int mode);
  void setTimerMode(int mode);
  void setInterval(float time);
  void setPulseLength(float time);
  void setPinNumbers(int array[], int size);
  
#ifndef PIN_ARR_SIZE
  static const int MAX_NUMBER=10;
#else
  static const int MAX_NUMBER=PIN_ARR_SIZE;
#endif
  
protected:
  
  int _num_pins;
  OutputPin _pins[MAX_NUMBER];
  
};

#endif