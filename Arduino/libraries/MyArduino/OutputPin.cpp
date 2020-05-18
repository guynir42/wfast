#include <OutputPin.h>

bool OutputPin::debug_bit=1;

OutputPin::OutputPin(int pin){
  
  setPinNum(pin);
  
  initialize();
    
}

OutputPin::~OutputPin(){
  
}

void OutputPin::initialize(){
  
  _current_mode=0;
  _inverse=0;
  //_state=0;
   digitalWrite(getPinNum(), LOW); // make sure it comes alive at LOW  
   pinMode(getPinNum(), OUTPUT);
  
  for(int i=0;i<Parser::MAXN;i++) _out_buf[i]=0;
  
  timer.resetTime();
  
}

char *OutputPin::printout(){
  
  update();
  
  snprintf(_out_buf, Parser::MAXN, "pin= %d | mode: %s | state= %d", getPinNum(), getModeString(), getStatus());
  
  return _out_buf;
  
}

// getters
int OutputPin::getPinNum() const {
  
  return _pin_num;
  
}

int OutputPin::getMode() const {
  
  return _current_mode;
  
}

char *OutputPin::getModeString() const {
  
  switch(getMode()){
   
    case OFF   : return "OFF";
    case ON    : return "ON";
    case WATCH : return "WATCH";
    
    default    : return "";
    
  }
  
}

int OutputPin::getStatus(){
  
  update();
  
  return digitalRead(getPinNum()); // This is unreliable!
  
}

int OutputPin::getState(){
  
  return getStatus();
  
}

char *OutputPin::getOutBuffer(){
  
  return _out_buf;
  
}

bool OutputPin::getInversionState() const {
  
  return _inverse;
  
}

// setters
void OutputPin::setPinNum(int number){
  
  _pin_num=number;
  
}

void OutputPin::setMode(int mode){
  
  // add verification!!
  _current_mode=mode;
  update();
  
}

void OutputPin::setInverse(){
  
  _inverse=!_inverse;
  
}

void OutputPin::setInverse(bool state){
  
  _inverse=state;
  
}

void OutputPin::update(){
      
  bool state;
  
  if(getMode()==ON) state=1;
  if(getMode()==OFF) state=0;
  
  if(getMode()==WATCH){
   
    if(timer.getState()) state=1;
    else state=0;
    
  }
  
  state=(state!=getInversionState()); // XOR for booleans! 
  
  if(state) digitalWrite(getPinNum(), HIGH);
  else digitalWrite(getPinNum(), LOW);
  
}

void OutputPin::parse(char *arg){
 
  char str1[Parser::STRN];
  char str2[Parser::MAXN];
 
  Parser::splitStr(arg, str1, str2);
  Parser::lowerString(str1);
  Parser::cleanString(str1);

  if(debug_bit){ Serial.print("OutputPin::parse "); Serial.print(str1); Serial.print(" | "); Serial.println(str2); Serial.flush();}  
  
  if(Parser::partialMatch(str1, "on") || Parser::partialMatch(str1, "yes") || Parser::partialMatch(str1, "1")){ setMode(ON); snprintf(_out_buf, Parser::STRN, "pin is ON"); }
  else if(Parser::partialMatch(str1, "off") || Parser::partialMatch(str1, "no") || Parser::partialMatch(str1, "0")){ setMode(OFF); snprintf(_out_buf, Parser::STRN, "pin is OFF"); }
  else if(Parser::partialMatch(str1, "watch") || Parser::partialMatch(str1, "timer") || Parser::partialMatch(str1, "2")){ setMode(WATCH); snprintf(_out_buf, Parser::STRN, "pin on WATCH"); }
  
  else if(Parser::partialMatch(str1, "interval")){ 
    timer.setInterval(atof(str2)); 
    char buf[Parser::STRN];
    Parser::ftoa(buf, timer.getInterval(), 3);
    snprintf(_out_buf, Parser::MAXN, "interval: %s", buf);    
  }
  else if(Parser::partialMatch(str1, "pulse_length")){
    timer.setPulseLength(atof(str2)); 
    char buf[Parser::STRN];
    Parser::ftoa(buf, timer.getPulseLength(), 3);
    snprintf(_out_buf, Parser::MAXN, "pulse: %s", buf);    
  }
  else if(Parser::partialMatch(str1, "mode")){
    timer.parse(str2);
    snprintf(_out_buf, Parser::MAXN, "com->timer: %s", str2);
  }
  else if(Parser::partialMatch(str1, "invert") || Parser::partialMatch(str1, "inverse")){
    if(Parser::parseBool(str2)==0) _inverse=0;
    else if(Parser::parseBool(str2)==1) _inverse=1;
    else _inverse=!_inverse;
  }
  
  else if(Parser::partialMatch(str1, "printout")) printout(); 
  
  else{
    snprintf(_out_buf, Parser::MAXN, "bad command: %s", str1);
  }
    
  Parser::println(_out_buf);
  
}










