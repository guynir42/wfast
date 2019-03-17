#include "InputPin.h"

bool InputPin::debug_bit=1;

InputPin::InputPin(int pin_number){
 
  _pin_number=pin_number;
  
  initialize();
  
}

InputPin::~InputPin(){
  
}

void InputPin::initialize(){
 
  setMode(DIGITAL);
  pinMode(_pin_number, INPUT);
  digitalWrite(_pin_number, LOW); // cancel the pull-up
  
  _push_num=0;
  _last_input=0;
  _last_push_time=millis();
  
  _has_changed=1;
  
  for(int i=0;i<STRN;i++) _out_buffer[i]=0;
  
  _threshold=700;
  _integrator_length=1000;
  _integrator_sum=0;
  _integrator_index=0;
  _integrator_result=0;
  
}

char *InputPin::printout(){
 
  snprintf(_out_buffer, STRN, "mode: %s val: %d", getModeString(), getValue());
  
  return _out_buffer;
  
}

int InputPin::getPinNumber(){
 
  return _pin_number;
  
}

bool InputPin::getState(){
  
  return getPinState();
  
}

bool InputPin::getPinState(){

  switch(_current_mode){
   
    case DIGITAL: return digitalRead(_pin_number);
    
    case PUSH: update(); return _push_num%2; 
    
    case SWITCH: return !digitalRead(_pin_number);
    
    case ANALOG: return getValue()>=_threshold;
      
    case ANTI_ANALOG: return getValue()>=_threshold;
      
    default: return 0;
  
  }
  
}

int InputPin::getMode(){
  
  return _current_mode;
  
}

char *InputPin::getModeString(){
  
  switch(_current_mode){
    
    case DIGITAL: return "DIGITAL";    
    case PUSH: return "PUSH";
    case SWITCH: return "SWITCH";
    default: return "unknonw";
  }
  
}

char *InputPin::getOutBuffer(){
  
  return _out_buffer;
  
}

int InputPin::getValue(){
 
  switch(_current_mode){
   
    case ANALOG: return analogRead(_pin_number);
    case ANTI_ANALOG: return 1023-analogRead(_pin_number);
    case INTEGRATOR: 
      
      _integrator_sum+=analogRead(_pin_number);
      _integrator_index++;
      if(_integrator_index>=_integrator_length){ 
	
	_integrator_result=_integrator_sum/_integrator_length;
	_integrator_sum=0;
	_integrator_index=0;
	
      }
      
      return _integrator_result;
    
    default : return 0;
  
  }  
  
}

bool ::InputPin::hasChanged(){
 
  if(_has_changed){
    
    _has_changed=0;
    return 1;
  }
  else return 0;
  
}

void InputPin::setMode(int mode){
 
  _current_mode=mode;
  
  switch(_current_mode){
   
    case DIGITAL: pinMode(_pin_number, INPUT); return;
    case PUSH: pinMode(_pin_number, INPUT); digitalWrite(_pin_number, HIGH); return;
    case SWITCH: pinMode(_pin_number, INPUT); digitalWrite(_pin_number, HIGH); return;
    case ANALOG: case ANTI_ANALOG: case INTEGRATOR: pinMode(_pin_number, INPUT); return;
    
  }
  
}

void InputPin::setPushNumber(int num){
  
  _push_num=num;
  
}

void InputPin::setPinState(bool state){
 
  // this works for push button only!
  _push_num=2*(_push_num%2)+state;
  
}

void InputPin::setThreshold(int thresh){
  
  _threshold=thresh;
  
}

void InputPin::update(){
         
  // wait until enough time passed. 
  if(_last_push_time+_push_delay<millis() && _last_input!=digitalRead(_pin_number)){
    
    _last_input=digitalRead(_pin_number); // update last value
    _last_push_time=millis(); // update last time
    if(!_last_input) _push_num++; // count only releases of ON input
    
  }
  
}

