#include "Parser.h"

bool Parser::debug_bit=1;

Parser::Parser(){
  
  initialize();
  
    
}

Parser::Parser(char *name){
  
  initialize();
  
  setThisSensorName(name);
  
}

Parser::Parser(char *name, char *version){
  
  initialize();
  
  setThisSensorName(name);
  setVersion(version);
  
}

Parser::~Parser(){
  
  clearCommands();
    
}

void Parser::initialize(){
  
  _num_commands=0;
  _is_done_reading=0;
  _input_pos=0;
    
  setThisSensorName("");
  setVersion("");

}

void Parser::clearCommands(){
  
  _num_commands=0;
  
}

void Parser::printout(){
  
  for(int i=0;i<_num_commands;i++){
    
    print("command ");
    print(i);
    println(_commands[i].keyword);
    
    
//     char buf[50];
//     snprintf(buf, 50, "command %d: %s", i, _commands[i].keyword);
//     println(buf);
    
  }// for i
  
}

void Parser::print(int number){
 
  Serial.print(number);
  Serial.flush();
  
}

void Parser::print(char *str){
 
  Serial.print(str);
  Serial.flush();
  
}

void Parser::print(const char *str){
   
  Serial.print(str);
  Serial.flush();
  
}

void Parser::println(int number){
  
  Serial.println(number);
  Serial.flush();
  
}

void Parser::println(char *str){
 
  Serial.println(str);
  Serial.flush();
  
}

void Parser::println(const char *str){
 
  Serial.println(str);
  Serial.flush();
  
}

void Parser::addCommand(char *keyword, FunctionPointer fp){
  
  snprintf(_commands[_num_commands].keyword, 16, "%s", keyword);
  _commands[_num_commands].fp=fp;
  
  _num_commands++;
  
  char buf[30];
  snprintf(buf, 30, "added keyword= %s", keyword);
  if(debug_bit) println(buf);
  snprintf(buf, 30, "_num_commands= %d", _num_commands);
  if(debug_bit) println(buf);
  
}

void Parser::readingInput(){
     
  while(Serial.available()){
    
    // get the new byte:
    _input_char = (char)Serial.read();

    if (_input_char=='\n' || _input_char=='\r' || _input_char==line_ending || _input_pos>=MAXN-1) {
      _is_done_reading = true;
      _input_string[_input_pos]=0;
      _input_pos=0;
      break;
    }
    // add it to the inputString:
    _input_string[_input_pos]=_input_char;
    _input_pos++;

  }
  
}

bool Parser::isStringComplete() const {
   
  return _is_done_reading;
  
}

void Parser::checkInput(){
 
  if(isStringComplete()){

    if(debug_bit){
      print("input: ");
      println(_input_string);
    }
    
    parse(_input_string);
    _is_done_reading=0;
    
  }
  
}

void Parser::getVersion() const {
 
  println(_version);
  
}

void Parser::getThisSensorName() const {
 
  print("Ard-Name: ");
  println(sensor_name);
  
}

void Parser::setVersion(char *version){
 
//   cleanString(version);
  
  snprintf(_version, MAXN, "%s", version);
  
}

void Parser::setThisSensorName(char *name){
 
  cleanString(name);
  
  snprintf(sensor_name, STRN, "%s", name);
//   println(sensor_name);
  
}

void Parser::parse(char *text){
  
  if(debug_bit){ 
    print("text: ");
    println(text);
  }
  
  char str1[STRN];
  char str2[MAXN];
  
  splitStr(text, str1, str2); // 2nd and 3rd arguments are output arguments...
  
  if(debug_bit){  
    print("str1: "); println(str1);
    print("str2: "); println(str2);
  }
  
  cleanString(str1);
  
  for(int i=0;i<_num_commands;i++){ // run over list of commands
    
    if(strcmp(_commands[i].keyword, str1)==0){
      
      FunctionPointer func=_commands[i].fp;
      
      if(debug_bit){
	print("func: ");
	println(_commands[i].keyword);
      }
      
      func(str2);
      return; // don't do anything else after any command is found 
      
    }
    
  }// for i
  
  // if no commands are listed, try the default commands...
  if(strcmp(str1, "version")==0){
   
    if(debug_bit) println("func: version (default)");
    getVersion();
    
  }
  
  if(strcmp(str1, "getname")==0){
   
    if(debug_bit) println("func: getname (default)");
    getThisSensorName();
    
  }
  
  if(strcmp(str1, "setname")==0){
   
    if(debug_bit) println("func: setname (default)");
    setThisSensorName(str2);
    
  }
  
  
  
}

void Parser::splitStr(char *input_str, char *out_first_word){
      
  strcpy(out_first_word, "");
  
  if(strlen(input_str)==0) return;  // treat empty inputs
  
  int i=0;
  for(;i<strlen(input_str);i++){
   
    if(i>=STRN-1 || i>=MAXN-1) break;
    
    int check=0;
    if(input_str[i]==line_ending || input_str[i]=='\n' || input_str[i]=='\r' || input_str[i]==0) check=1;
    for(int j=0;j<strlen(separators);j++) if(input_str[i]==separators[j]) check=1;
    
    if(check) break;
    
    out_first_word[i]=input_str[i];
    
  }
  
  out_first_word[i]=0; // add terminating character
  
}

void Parser::splitStr(char *input_str, char *out_first_word, char *out_second_word){
  
  strcpy(out_first_word, ""); 
  strcpy(out_second_word, "");
  
  if(strlen(input_str)==0) return;  // treat empty inputs
  
  int i=0;
  for(;i<strlen(input_str);i++){
   
    if(i>=STRN-1 || i>=MAXN-1) break;
    
    int check=0;
    if(input_str[i]==line_ending || input_str[i]=='\n' || input_str[i]=='\r' || input_str[i]==0) check=1;
    for(int j=0;j<strlen(separators);j++) if(input_str[i]==separators[j]) check=1;
    
    if(check) break;
    
    out_first_word[i]=input_str[i];
    
  }
  
  out_first_word[i]=0; // add terminating character
  
  i++; // skip the separator
  
  int j=0;
  
  for(;i<strlen(input_str);i++){
  
    if(j>=STRN-1 || i>=MAXN-1) break;
    
    if(input_str[i]==line_ending || input_str[i]=='\n' || input_str[i]=='\r' || input_str[i]==0) break;
    
    out_second_word[j]=input_str[i];
    j++;
        
  }

  out_second_word[j]=0;
  
}

int Parser::partialMatch(char *str, char *comp_str){
 
  cleanString(str);
  
  if(strlen(str)==0) return 0;
  
  for(int i=0;i<strlen(str); i++){
   
    if(i>=strlen(comp_str)) return 1; // if you finished reading the comp_str that's a win
    
    if(str[i]!=comp_str[i]) return 0;// if any of the few letters in str mismatch, you lose
    
  }
  
  return 1; // if you finished all letters in str, you win
  
}

void Parser::itoa(char *buf, int number){
  
  snprintf(buf, STRN, "%d", number);
   
}

void Parser::ftoa(char *buf, float number, int precision){

  long p[] = {0,10,100,1000,10000,100000,1000000,10000000,100000000};
 
  long integers = (long)number;
  snprintf(buf, STRN, "%ld", integers);
  
  long decimal = abs((long)((number-integers)*p[precision]));
  
  char buf2[STRN]={0};
  
  switch(precision){
    case 1: snprintf(buf2, STRN, ".%01ld", decimal); break;
    case 2: snprintf(buf2, STRN, ".%02ld", decimal); break;
    case 3: snprintf(buf2, STRN, ".%03ld", decimal); break;
    case 4: snprintf(buf2, STRN, ".%04ld", decimal); break;
    case 5: snprintf(buf2, STRN, ".%05ld", decimal); break;
    case 6: snprintf(buf2, STRN, ".%06ld", decimal); break;
    case 7: snprintf(buf2, STRN, ".%07ld", decimal); break;
    case 8: snprintf(buf2, STRN, ".%08ld", decimal); break;
    case 9: snprintf(buf2, STRN, ".%09ld", decimal); break;
  }
  
  snprintf(buf, STRN, "%s%s", buf, buf2);
    
  
}

int Parser::parseBool(char *str){
  
  cleanString(str);
  lowerString(str);
  
  if(strcmp(str, "no")==0 || strcmp(str, "off")==0 || str[0]=='0') return 0;
  if(strcmp(str, "yes")==0 || strcmp(str, "on")==0 || str[0]=='1') return 1;
  
  return -1; // if no command is recognized!
    
}

void Parser::lowerString(char *str){
 
  for(int i=0;i<strlen(str);i++){
   
    str[i]=tolower(str[i]);
    
  }
  
}

void Parser::cleanString(char *str){
 
  if(strlen(str)==0) return;
  
  char new_str[STRN+1];
  int N=strlen(str);
  for(;N>0;N--) if(str[N-1]!=' '&&str[N-1]!='\n'&&str[N-1]!='\r') break; // find end of string without spaces
  
  int i=0;
  for(;i<strlen(str);i++) if(str[i]!=' '&&str[i]!='\n'&&str[i]!='\r') break;// find start of string without spaces
   
  int j=0;
  for(;i<N;i++){ 
    new_str[j]=str[i];
    j++;
  }
  
  new_str[j]=0; // add end of string
  
  snprintf(str, STRN, "%s", new_str);
  
}

char Parser::separators[]=",|";
char Parser::line_ending=';';

