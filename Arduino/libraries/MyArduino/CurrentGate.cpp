#include "CurrentGate.h"

CurrentGate::CurrentGate(short int pin1, short int pin2){
	
	blockCurrent();
	_state=0;
	
}

bool CurrentGate::getState() const {
	
	return _state;
	
}

bool CurrentGate::getStatus() const {
	
	return getState();
	
}

bool CurrentGate::getInverse() const {
	
	return _inverse;
	
}

void CurrentGate::setState(bool state){
	
	_state=state;
	
	if(_state && !_inverse) allowCurrent();
	else blockCurrent();
	
}

void CurrentGate::setInverse(bool inverse){
	
	_inverse = inverse;
	
}

void CurrentGate::allowCurrent(){
	
	digitalWrite(_pin1, LOW);
	pinMode(_pin1, OUTPUT);
	digitalWrite(_pin2, LOW);
	pinMode(_pin2, OUTPUT);
	
}

void CurrentGate::blockCurrent(){
	
	digitalWrite(_pin1, LOW);
	pinMode(_pin1, INPUT);
	digitalWrite(_pin2, LOW);
	pinMode(_pin2, INPUT);
	
}

