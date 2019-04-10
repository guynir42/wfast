#include "Ultrasonic.h"

Ultrasonic::Ultrasonic(int trig_pin, int echo_pin){
	
	_trig_pin=trig_pin;
	_echo_pin=echo_pin;
	
	measure();
	
}

float Ultrasonic::measure(){
	
  // The sensor is triggered by a HIGH pulse of 10 or more microseconds.
  // Give a short LOW pulse beforehand to ensure a clean HIGH pulse:
  pinMode(_trig_pin, OUTPUT);
  digitalWrite(_trig_pin, LOW);
  delayMicroseconds(_short_delay_microseconds);
  digitalWrite(_trig_pin, HIGH);
  delayMicroseconds(_delay_microseconds);
  digitalWrite(_trig_pin, LOW);
 
  // Read the signal from the sensor: a HIGH pulse whose
  // duration is the time (in microseconds) from the sending
  // of the ping to the reception of its echo off of an object.
  pinMode(_echo_pin, INPUT);
  _duration=pulseIn(_echo_pin, HIGH);
  _cm=time2cm(_duration);
  return _cm;
	
}

float Ultrasonic::time2cm(float duration_microseconds){
	
	// speed of sound 340 m/s = 29 cm/microsec 
	// also divide by 2 for return trip of echo
	return duration_microseconds/2/29; 
	
}

float Ultrasonic::getDuration(){
	
	return _duration;
	
}

float Ultrasonic::getDistance(){
	
	return _cm;
	
}