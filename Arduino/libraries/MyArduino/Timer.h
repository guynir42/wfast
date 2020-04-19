#ifndef TIMER_H
#define TIMER_H

#include <Arduino.h>
#include <EEPROM.h>
#include <Parser.h>

class Timer{
  
public:
  
  Timer();
  ~Timer();
  
  void initialize();
  void time2hms(char *buf, float time, int units=-1);
  char * printout();
  
  // getters
  float getRunningTime();
  float getTimeInInterval();
  float getRemainingInInterval();
  float getTimeInPulse();
  float getRemainingInPulse();
  int getState();
  float getInterval() const;
  float getPulseLength() const;
  int getMode() const;
  char *getModeString() const;
  char *getOutBuffer();
  bool getFirstPulse() const;
  void getTimeUnitStr(char *buf);
  char getTimeUnitsChar(); 
  
  // setters
  void resetTime();

  void setMode(int mode);
  void setUnit(int unit);
  void setInterval(float interval);
  void setPulseLength(float pulse);
    
  void setupCountdown(float interval=-1);
  void setupExpiration(float interval=-1);
  void setupAlternator(float interval=-1);
  void setupPulsed(float interval=-1, float pulse=-1);
  void setupSingle(float interval=-1, float pulse=-1);
  void setFirstPulse(bool first=1);
  
  void setRunningTime(float time);
  void beginInterval();
  void beginPulse();
  
  void parse(char *arg);
  void parseTiming(char *arg);
  
  void loadEEPROM();
  void saveEEPROM();
  void attachEEPROM(int position=1, int mode=0, int unit=0, float interval=5, float pulse_length=1, bool pulse_first=0);
  
  static bool debug_bit;
  static int mem_size;
  
  enum mode_type {COUNT=0, EXPIRE, ALTER, PULSED, SINGLE, OFF, ON};
  enum time_unit {SEC=0, MIN, HOUR};
  
protected:
  
  unsigned long long int _start_time_milli;
  
  int _current_mode;
  int _current_unit;
  
  float _interval; // in seconds!!
  float _pulse_length; // in seconds!!
  
  bool _first_pulse; // first do the pulse, then the interval!
  
  char _out_buf[Parser::MAXN];
  
  int _eeprom_pos_start;
  int _eeprom_pos_end;
  
};

#endif