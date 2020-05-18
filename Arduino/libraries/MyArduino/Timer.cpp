#include "Timer.h"

bool Timer::debug_bit=1;
int Timer::mem_size=25;

Timer::Timer(){
  
  initialize();
  
}

Timer::~Timer(){
  
}

void Timer::initialize(){
  
  _current_mode=0; // default behavior
  _current_unit=0; // default time unit (seconds)
  _interval=5;
  _pulse_length=1;
  
  for(int i=0;i<Parser::MAXN;i++) _out_buf[i]=0;
  
  _first_pulse=0;
  
  _eeprom_pos_start=0;
  _eeprom_pos_end=0;
  
  resetTime();
  
}

void Timer::time2hms(char *buf, float time, int units){
	
	if(units<0) units=_current_unit;
	
	bool minus=0;
	
	if(time<0) minus=1;
	
	int time_sec;
	switch(units){
		case SEC: time_sec=abs(time); break;
		case MIN: time_sec=abs(time)*60; break;
		case HOUR: time_sec=abs(time)*3600; break;
	}
	
	int hours=time_sec/3600;
	int minutes= (time_sec-hours*3600)/60;
	int seconds= time_sec%60;
	
	snprintf(buf, Parser::STRN, "");
	if(minus) snprintf(buf, Parser::STRN, "-");
	if(hours<=0 && minutes<=0 && seconds<=0) snprintf(buf, Parser::STRN, "0");
	if(hours>0) snprintf(buf, Parser::STRN, "%s%dh", buf, hours);
	if(minutes>0 || (hours>0 && seconds>0) ) snprintf(buf, Parser::STRN, "%s%dm", buf, minutes);
	if(seconds>0) snprintf(buf, Parser::STRN, "%s%ds", buf, seconds);
	
}

char *Timer::printout(){
  
  char interval[Parser::STRN];
  // Parser::ftoa(interval, getInterval(), 2);
  time2hms(interval, getInterval());
  
  char pulse[Parser::STRN];
  // Parser::ftoa(pulse, getPulseLength(), 2);
  time2hms(pulse, getPulseLength());
  
  //char unit_str[10];
  //getTimeUnitStr(unit_str);
  
  char remaining[Parser::STRN];
    
  switch(_current_mode){
  
	case COUNT:
		if(getState()){
			snprintf(_out_buf, Parser::MAXN, "%s (ON) | int= %s (passed)", getModeString(), interval);
		}
		else{
			time2hms(remaining, getRemainingInInterval());
			snprintf(_out_buf, Parser::MAXN, "%s (OFF) | int= %s of %s", getModeString(), remaining, interval);
		}
		break;
	
	case EXPIRE:

		if(getState()){
			time2hms(remaining, getRemainingInInterval());
			snprintf(_out_buf, Parser::MAXN, "%s (ON) | int= %s of %s", getModeString(), remaining, interval);
			
		}
		else{
			snprintf(_out_buf, Parser::MAXN, "%s (OFF) | int= %s (passed)", getModeString(), interval);
		}
		break;
		
	case ALTER:
		time2hms(remaining, getTimeInInterval());
		if(getState()){
			snprintf(_out_buf, Parser::MAXN, "%s (ON) | int= %s of %s", getModeString(), remaining, interval);
		}
			
		else{
			snprintf(_out_buf, Parser::MAXN, "%s (OFF) | int= %s of %s", getModeString(), remaining, interval);
		}
		break;
	
	case PULSED:	
	
		if(getState()){
		
			time2hms(remaining, getTimeInPulse());
			snprintf(_out_buf, Parser::MAXN, "%s (ON) | int= %s | pul: %s of %s ", getModeString(), interval, remaining, pulse);
	  
		}
		else{
		  
			time2hms(remaining, getTimeInInterval());
			snprintf(_out_buf, Parser::MAXN, "%s (OFF) | int= %s of %s | pul: %s ", getModeString(), remaining, interval, pulse);
		  
		}
	  
		break;

	case SINGLE:
	  
		if(getState()){
		  
			time2hms(remaining, getRemainingInPulse());
			snprintf(_out_buf, Parser::MAXN, "%s (ON) | int: %s | pul: %s of %s ", getModeString(), interval, remaining, pulse);
	  
		}
		else if(getFirstPulse()==0 && getRunningTime()<getInterval()){ // if we are before the one and only pulse
		  
			time2hms(remaining, getRemainingInInterval());
			snprintf(_out_buf, Parser::MAXN, "%s (OFF) | int: %s of %s | pul: %s ", getModeString(), remaining, interval, pulse);
		  
		}
		else{
			
			snprintf(_out_buf, Parser::MAXN, "%s (OFF) | int: %s | pul: %s (passed)", getModeString(), interval, pulse);
			
		}
	  
		break;
	
	case ON:
	
		snprintf(_out_buf, Parser::MAXN, "timer is always ON!");
		break;
		
	case OFF:
	
		snprintf(_out_buf, Parser::MAXN, "timer is always OFF!");
		break;
	  
  } // switch
  
  return _out_buf;
  
}

// getters
float Timer::getRunningTime(){
 
	long int milli_seconds_mod = (millis()-_start_time_milli); 
	
	if (_current_mode==ALTER) milli_seconds_mod = milli_seconds_mod%((long int) (_interval*2*1000) ); // keep this number inside the range of ON and OFF intervals (convert _interval to ms)
	if (_current_mode==PULSED) milli_seconds_mod = milli_seconds_mod%((long int)((_interval+_pulse_length)*1000)); // keep this number inside the range of OFF interval and ON pulse (convert _interval/_pulse_length to ms)
	
	float time=((float) milli_seconds_mod)/1000.0; // convert to seconds!
	if(_current_unit==Timer::MIN) time=time/60;
	else if(_current_unit==Timer::HOUR) time=time/3600;
	if(time<0) resetTime(); // to handle when the Arduino built in timer resets itself
	return time;
  
}

float Timer::getTimeInInterval(){
	
	switch(_current_mode){
		
		case COUNT:
		case EXPIRE:
			if(getRunningTime()>_interval) return _interval;
			else return getRunningTime();
			
		case ALTER:
			if(getRunningTime()<=_interval) return getRunningTime();
			if(getRunningTime()<=2*_interval) return getRunningTime() - _interval;
			
		case PULSED: 
		case SINGLE: 
			if(_first_pulse==0){
				if(getRunningTime()<=_interval) return getRunningTime(); // still in the interval
				else return _interval; // interval has passed...
			}
			else{
				if(getRunningTime()>=_pulse_length && getRunningTime()<_interval+_pulse_length) return getRunningTime() - _pulse_length; // passed the initial pulse, now inside the interval...
				else if (getRunningTime()<_pulse_length) return 0;
				else return _pulse_length + _interval;
			}
		
		default:
			return 0;	

	}// switch 
	
				
	
}

float Timer::getRemainingInInterval(){
	
	switch(_current_mode){
				
		case COUNT:
		case EXPIRE:
			if(getRunningTime()>_interval) return 0;
			else return _interval-getRunningTime();
			
		case ALTER:
			if(getRunningTime()<_interval) return _interval-getRunningTime();
			if(getRunningTime()<2*_interval) return 2*_interval-getRunningTime();
			
		case PULSED: 
		case SINGLE: 
			if(_first_pulse==0){
				if(getRunningTime()<=_interval) return _interval-getRunningTime(); // still in the interval
				else return 0; // interval has passed...
			}
			else{
				if(getRunningTime()>=_pulse_length && getRunningTime()<_interval+_pulse_length) return _interval+_pulse_length-getRunningTime(); // passed the initial pulse, now inside the interval...
				else return 0;
			}
		
		default:
			return 0;	

	}// switch 
	
}

float Timer::getTimeInPulse(){
		
		switch(_current_mode){
				
		case COUNT:
		case EXPIRE:
			return 0;
			
		case ALTER:
			return 0;
			
		case PULSED: 
		case SINGLE: 
			if(_first_pulse==0){
				if(getRunningTime()>_interval && getRunningTime()<=_interval+_pulse_length) return getRunningTime() - _interval; // passed the interval, now in the pulse
				else if(getRunningTime()<_interval) return 0;
				else return _pulse_length;
			}
			else{
				if(getRunningTime()>=0 && getRunningTime()<_pulse_length) return getRunningTime(); // still inside the initial pulse
				else return _pulse_length;
			}
		
		default:
			return 0;	

	}// switch 
	
	
}

float Timer::getRemainingInPulse(){
		
		switch(_current_mode){
				
		case COUNT:
		case EXPIRE:
			return 0;
			
		case ALTER:
			return 0;
			
		case PULSED: 
		case SINGLE: 
			if(_first_pulse==0){
				if(getRunningTime()>_interval && getRunningTime()<=_interval+_pulse_length) return _interval+_pulse_length-getRunningTime(); // passed the interval, now in the pulse
				else return 0;
			}
			else{
				if(getRunningTime()>=0 && getRunningTime()<_pulse_length) return _pulse_length-getRunningTime(); // still inside the initial pulse
				else return 0;
			}
		
		default:
			return 0;	

	}// switch 
	
	
}

int Timer::getState(){
  
  switch(_current_mode){
   
    case COUNT:
      
      if(_interval<=getRunningTime()) return 1;
      else return 0;
         
    case EXPIRE:
      
      if(_interval>=getRunningTime()) return 1;
      else return 0;
      
    case ALTER:
      
      if(((int)(getRunningTime()/_interval))%2==0) return _first_pulse;
      if(getRunningTime()>=_interval*2){ return _first_pulse; }
      else return !_first_pulse;
      
    case PULSED:      
    case SINGLE:
	
      if(_first_pulse==0){
		if(getRunningTime()>=_interval && getRunningTime()<=_interval+_pulse_length) return 1;
		else if(getRunningTime()>_interval+_pulse_length){ return 0; }
		else return 0;
      }
      else{
		if(getRunningTime()>=0 && getRunningTime()<=_pulse_length) return 1;
		else if(getRunningTime()>_pulse_length+_interval){ return 0; }
		else return 0;
      }
	
    case OFF:
      return 0;
      
    case ON:
      return 1;
      
    default:
      
      return 0;
            
  }
  
}

float Timer::getInterval() const {
  
  return _interval;
  
}

float Timer::getPulseLength() const { 
  
  return _pulse_length;
  
}

int Timer::getMode() const {
  
  return _current_mode;
  
}

char *Timer::getModeString() const {
  
  switch(_current_mode){
    
    case COUNT:
      return "COUNT";
    case EXPIRE:
      return "EXPIRE";
    case ALTER:
      return "ALTER";
    case PULSED:
      return "PULSED";
    case SINGLE:
      return "SINGLE";
    case OFF:
      return "OFF";
    case ON:
      return "ON";
    default:
      return "";
    
  }
  
}

char *Timer::getOutBuffer(){
  
  return _out_buf;
  
}

bool Timer::getFirstPulse() const {
 
  return _first_pulse;
  
}

void Timer::getTimeUnitStr(char *buf){
	
	if(_current_unit==Timer::SEC) snprintf(buf, 2, "s");
	else if(_current_unit==Timer::MIN) snprintf(buf, 2, "m");
	else if(_current_unit==Timer::HOUR) snprintf(buf, 2, "h");
		
	
}

char Timer::getTimeUnitsChar(){
	
	if(_current_unit==Timer::SEC) return 's';
	else if(_current_unit==Timer::MIN) return 'm';
	else if(_current_unit==Timer::HOUR) return 'h';
	
}

// setters
void Timer::resetTime(){
  
  // Serial.println("reset time");
  
  _start_time_milli=millis();
  
}

void Timer::setMode(int mode){
 
  resetTime();
  _current_mode=mode;
  
	if(_eeprom_pos_start>0 && _eeprom_pos_end>0){ 
	
		double mode_double=mode;
		EEPROM.put(_eeprom_pos_start+0, mode_double);
		
	}
	
}

void Timer::setUnit(int unit){
	
	resetTime();
	_current_unit=unit;
		
	if(_eeprom_pos_start>0 && _eeprom_pos_end>0){ 
	
		double unit_double=unit;
		EEPROM.put(_eeprom_pos_start+4, unit_double);
		
	}
	
}

void Timer::setInterval(float interval){
  
  _interval=interval;
  
  if(_eeprom_pos_start>0 && _eeprom_pos_end>0){

	double interval_double=interval;
	EEPROM.put(_eeprom_pos_start+8, interval_double);
	
  }
  
}

void Timer::setPulseLength(float pulse){
 
  _pulse_length=pulse;
  
  if(_eeprom_pos_start>0 && _eeprom_pos_end>0){ 
	
	double pulse_double=pulse;
	EEPROM.put(_eeprom_pos_start+16, pulse_double);
	
  }
  
}

void Timer::setupCountdown(float interval){
  
  if(interval<=0) interval=_interval;
  
  resetTime();
  setInterval(interval);
  setMode(COUNT);
  
}

void Timer::setupExpiration(float interval){
    
  if(interval<=0) interval=_interval;
  
  resetTime();
  setInterval(interval);
  setMode(EXPIRE);
}

void Timer::setupAlternator(float interval){
  
  if(interval<=0) interval=_interval;
  
  resetTime();
  setInterval(interval);
  setMode(ALTER);
  
}

void Timer::setupPulsed(float interval, float pulse){
  
  if(interval<=0) interval=_interval;
  if(pulse<=0) pulse=_pulse_length;
  
  resetTime();
  setInterval(interval);
  setPulseLength(pulse);
  setMode(PULSED);
  
}

void Timer::setupSingle(float interval, float pulse){
     
  if(interval<=0) interval=_interval;
  if(pulse<=0) pulse=_pulse_length;
  
  resetTime();
  setInterval(interval);
  setPulseLength(pulse);
  setMode(SINGLE);
  
}

void Timer::setFirstPulse(bool first){
  
  _first_pulse=first;
  
    if(_eeprom_pos_start>0 && _eeprom_pos_end>0){ 
	
		double first_double=first;
		EEPROM.put(_eeprom_pos_start+24, first_double);
	
	}

  
}

void Timer::setRunningTime(float time){
 
  float time_ms = time*1000;
  if(_current_unit==Timer::MIN) time_ms=time_ms*60;
  else if(_current_unit==Timer::HOUR) time_ms=time_ms*3600;
  
  _start_time_milli = millis()-time_ms;
  
}

void Timer::beginInterval(){
 
  switch(_current_mode){
  
    case COUNT: setRunningTime(getInterval()); break;
    
    case ALTER: setRunningTime(_first_pulse*_interval);
      
    case PULSED: 
      if(_first_pulse) setRunningTime(getPulseLength());
      else resetTime();
      break;
    
    // resetTime(); break;
    
    default: resetTime(); break;
    
  }
  
}

void Timer::beginPulse(){
 
  switch(_current_mode){
  
    case COUNT: resetTime(); break;
    
    case ALTER: setRunningTime((!_first_pulse)*_interval);
      
    case PULSED: 
      if(!_first_pulse) setRunningTime(getInterval());
      else resetTime();
      break;
    
    default: setRunningTime(getInterval()); break;
    
  }
  
}

void Timer::parse(char *arg){
     
  char str1[Parser::STRN];
  char str2[Parser::MAXN];
  
  Parser::splitStr(arg, str1, str2);
  Parser::lowerString(str1);
  Parser::cleanString(str1);
  
  if(debug_bit){ Serial.print("Timer::parse "); Serial.print(str1); Serial.print(" | "); Serial.println(str2); }  

  if(Parser::partialMatch(str1, "countdown")){ 
    
    setMode(COUNT); 
    parseTiming(str2);
    char buf[12];
    Parser::ftoa(buf, getInterval(), 2);
    snprintf(_out_buf, Parser::MAXN, "timer COUNT: %s s", buf); 
    
  }
  else if(Parser::partialMatch(str1, "expire") || Parser::partialMatch(str1, "expiration")){ 
    
    setMode(EXPIRE); 
    parseTiming(str2);
    char buf[12];
    Parser::ftoa(buf, getInterval(), 2);
    snprintf(_out_buf, Parser::MAXN, "timer EXPIRE: %s s", buf); 
    
  }
  else if(Parser::partialMatch(str1, "alternate") || Parser::partialMatch(str1, "alternator")){ 
    
    setMode(ALTER);
    parseTiming(str2);
    char buf[12];
    Parser::ftoa(buf, getInterval(), 2);
    snprintf(_out_buf, Parser::MAXN, "timer ALTER: %s s", buf); 
    
  }
  else if(Parser::partialMatch(str1, "pulsed")){ 
    
    setMode(PULSED); 
    parseTiming(str2);
    char buf[12];
    Parser::ftoa(buf, getInterval(), 2);
    char buf2[12];
    Parser::ftoa(buf2, getPulseLength(), 2);
    snprintf(_out_buf, Parser::MAXN, "timer PULSED: %s/%s s", buf, buf2); 
    
  }
  
  else if(Parser::partialMatch(str1, "single")){ 
    
    setMode(SINGLE); 
    parseTiming(str2);
    char buf[12];
    Parser::ftoa(buf, getInterval(), 2);
    char buf2[12];
    Parser::ftoa(buf2, getPulseLength(), 2);
    snprintf(_out_buf, Parser::MAXN, "timer SINGLE: %s/%s s", buf, buf2); 
    
  }
  
  else if(Parser::partialMatch(str1, "off")){ 
    
    setMode(OFF); 
    snprintf(_out_buf, Parser::MAXN, "timer is OFF"); 
    
  }
  
  else if(Parser::partialMatch(str1, "on")){ 
    
    setMode(ON); 
    snprintf(_out_buf, Parser::MAXN, "timer ON"); 
    
  }
  
  else{
    snprintf(_out_buf, Parser::STRN, "bad command: %s", str1);
  }
  
  if(debug_bit) Parser::println(_out_buf);
  
}

void Timer::parseTiming(char *arg){
         
  char str1[Parser::STRN];
  char str2[Parser::MAXN];
  
  Parser::splitStr(arg, str1, str2);
  Parser::cleanString(str1);
  Parser::cleanString(str2);
  
  if(strlen(str1)>0 && atof(str1)>0){ 
    
    setInterval(atof(str1));
    
    if(strlen(str2)>0 && atof(str2)>0){
     
      setPulseLength(atof(str2));
      
    }
    
  }
  
}

void Timer::loadEEPROM(){
	
	// double mode, unit, interval, pulse_length, pulse_first;
	int mode, unit;
	float interval, pulse_length;
	bool pulse_first;
	
	EEPROM.get(_eeprom_pos_start+0, mode);
	EEPROM.get(_eeprom_pos_start+4, unit);
	EEPROM.get(_eeprom_pos_start+8, interval);
	EEPROM.get(_eeprom_pos_start+16, pulse_length);
	EEPROM.get(_eeprom_pos_start+24, pulse_first);
	
	if(mode>=0 && mode<7) _current_mode=mode;
	if(unit>=0 && unit<3) _current_unit=unit;
	if(interval>0 && interval<1e8) _interval=interval;
	if(pulse_length>0 && pulse_length<1e8) _pulse_length=pulse_length;
	if(pulse_first>=0 || pulse_first<2) _first_pulse=pulse_first;
	
}

void Timer::saveEEPROM(){
	
	EEPROM.put(_eeprom_pos_start+0, _current_mode);
	EEPROM.put(_eeprom_pos_start+4, _current_unit);
	EEPROM.put(_eeprom_pos_start+8, _interval);
	EEPROM.put(_eeprom_pos_start+16, _pulse_length);
	EEPROM.put(_eeprom_pos_start+24, _first_pulse);
	
	double mem_check=0;
	EEPROM.get(0, mem_check);
	if(mem_check<_eeprom_pos_end || mem_check>1000 || ((int) mem_check)!=mem_check ){// update mem_check (first EEPROM position) if needed
		
		mem_check=_eeprom_pos_end;
		EEPROM.put(0, mem_check);
		
	}
	
}

void Timer::attachEEPROM(int position, int mode, int unit, float interval, float pulse_length, bool pulse_first){
	
	_eeprom_pos_start=position;
	_eeprom_pos_end=_eeprom_pos_start+mem_size;
	
	int mem_check=0;
	EEPROM.get(0, mem_check);
	
	if(((int) mem_check)==mem_check && mem_check<1000 && mem_check>=_eeprom_pos_end){ // if mem_check value looks ok, load data from memory
		// Serial.print("loading from EEPROM... memory starts at: ");
		// Serial.println(_eeprom_pos_start);
		loadEEPROM();
	}
	else{ // assume memory is not formatted, use defaults or given parameters
		
		_current_mode=mode;
		_current_unit=unit;
		_interval=interval;
		_pulse_length=pulse_length;
		_first_pulse=pulse_first;
		saveEEPROM();
	}
		
}