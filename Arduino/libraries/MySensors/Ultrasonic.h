#ifndef ULTRASONIC_H
#define ULTRASONIC_H

#include "Arduino.h"

class Ultrasonic {
	
public:

	Ultrasonic(int trig_pin, int echo_pin);
	
	float measure(); // send a pulse, wait a few microseconds, then read the response from sensor... 
	
	float time2cm(float duration_microseconds);
	float getDuration(); // in microseconds
	float getDistance(); // in cm
	
	
private:
	int _trig_pin;
	int _echo_pin;
	
	int _delay_microseconds=12;
	int _short_delay_microseconds=2;
	
	float _duration;
	float _cm; 
	
	
};

#endif