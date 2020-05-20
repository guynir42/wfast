#include "Seeed_BME280.h"
#include <Wire.h>

#include "Parser.h"
#include "OutputPin.h"
#include "GroundPin.h";
#include "VoltagePin.h";
#include "InputPin.h";

const int buf_size=10; // buffer the last few measurements

BME280 bme280; // the sensor SDA and SCL should be connected to the ATMega2560 pins 20 and 21
GroundPin bmeGND(19);
VoltagePin bmeVCC(18); 

InputPin button(A7);
VoltagePin buttonVCC(A6);
GroundPin buttonGND(A5); 

OutputPin relay1(12); 
OutputPin relay2(11); 
OutputPin relay3(10); 
OutputPin relay4(9); 
OutputPin relay5(8); 
OutputPin relay6(7); 
OutputPin relay7(6); 
OutputPin relay8(5); 

// remote control
VoltagePin remoteVCC(53); // the end of that row is use as ground for the remote
InputPin remoteD0(51); 
InputPin remoteD1(49); 
InputPin remoteD2(47); 
InputPin remoteD3(45); 
InputPin remoteALL(43); 


int index=0; // running index is looped back when reaching buf_size

Parser parser("assist", "DomeAssistant v1.00");

// buffers for keeping latest measurements
float pressure[buf_size]={0};
float temperature[buf_size]={0};
float humidity[buf_size]={0};

void setup() {
  Serial.begin(9600);
    
  if (!bme280.init()) {
    Serial.println("Device error!");
  }

  Parser::debug_bit=0;
  OutputPin::debug_bit=0;
  Timer::debug_bit=0;

  parser.addCommand("measure", &measure);
  parser.addCommand("status", &statusReport);
  parser.addCommand("button", &press_button); 
  parser.addCommand("relay1", &relay1func);
  parser.addCommand("relay2", &relay2func);
  parser.addCommand("relay3", &relay3func);
  parser.addCommand("relay4", &relay4func);
  parser.addCommand("relay5", &relay5func);
  parser.addCommand("relay6", &relay6func);
  parser.addCommand("relay7", &relay7func);
  parser.addCommand("relay8", &relay8func);
  parser.addCommand("help", &help);

  relay1.timer.setMode(Timer::OFF);
  relay1.timer.setInterval(300); 
  relay2.timer.setMode(Timer::OFF);
  relay3.timer.setMode(Timer::OFF);
  relay4.timer.setMode(Timer::OFF);
  relay5.timer.setMode(Timer::OFF);
  relay6.timer.setMode(Timer::OFF);
  relay7.timer.setMode(Timer::OFF);
  relay8.timer.setMode(Timer::OFF);

  
  
}

void serialEvent() {

  parser.readingInput();

}

void loop() {

    parser.checkInput(); // check if any command is issued...
    
    pressure[index]=bme280.getPressure()/100; // translate pressure from Pascal to mbar
    temperature[index]=bme280.getTemperature();
    humidity[index]=bme280.getHumidity();

    if(button.getState()){ // button is only noticed when lights are off! 
      press_button("");
    }
    
    relay1.update();
    relay2.update();
    relay3.update();
    relay4.update();
    relay5.update();
    relay6.update();
    relay7.update();
    relay8.update();

    index++;
    if(index>=buf_size) index=0;
    
//    delay(1000);

//    measure(""); 

    if(remoteALL.getState()){
      press_button(""); 
    }
    
    
}

void measure(char *arg){

  Serial.print("P= ");
  Serial.print(average_pressure()); 
  Serial.print(" | ");
  
  Serial.print("T= ");
  Serial.print(average_temperature()); 
  Serial.print(" | ");
  
  Serial.print("H= ");
  Serial.print(average_humidity()); 
  Serial.println();

}

void statusReport(char *arg) {

  Serial.println("This is dome-assistant status report!"); 
  Serial.println("Relays: "); 
  Serial.print(relay1.printout()); if(relay1.getMode()==OutputPin::WATCH) { Serial.print(" | "); Serial.print(relay1.timer.printout()); } Serial.println(); 
  Serial.print(relay2.printout()); if(relay2.getMode()==OutputPin::WATCH) { Serial.print(" | "); Serial.print(relay2.timer.printout()); } Serial.println();  
  Serial.print(relay3.printout()); if(relay3.getMode()==OutputPin::WATCH) { Serial.print(" | "); Serial.print(relay3.timer.printout()); } Serial.println();  
  Serial.print(relay4.printout()); if(relay4.getMode()==OutputPin::WATCH) { Serial.print(" | "); Serial.print(relay4.timer.printout()); } Serial.println();  
  Serial.print(relay5.printout()); if(relay5.getMode()==OutputPin::WATCH) { Serial.print(" | "); Serial.print(relay5.timer.printout()); } Serial.println();  
  Serial.print(relay6.printout()); if(relay6.getMode()==OutputPin::WATCH) { Serial.print(" | "); Serial.print(relay6.timer.printout()); } Serial.println();  
  Serial.print(relay7.printout()); if(relay7.getMode()==OutputPin::WATCH) { Serial.print(" | "); Serial.print(relay7.timer.printout()); } Serial.println();  
  Serial.print(relay8.printout()); if(relay8.getMode()==OutputPin::WATCH) { Serial.print(" | "); Serial.print(relay8.timer.printout()); } Serial.println();  

  measure(arg); 
  
}

void press_button(char *arg){

    //if(relay1.getState()==0){ // only do something if light is currently off...
      relay1.setMode(OutputPin::WATCH); 
      relay1.timer.setupExpiration(); // default is 5 minutes 
      if(relay1.timer.getInterval()<0) relay1.timer.setInterval(300); 
      Serial.print("Turning on the lights for "); 
      Serial.print(relay1.timer.getInterval()/60);
      Serial.println(" minutes!");      
    //}

}

void help(char *arg){

  parser.getVersion();

  Serial.println("This is the help section for DomeAssistant"); 
  Serial.println("To control individual relays, use e.g., 'relay1, on', or 'relay1, off'"); 
  Serial.println("To enable timer use e.g., 'relay1, watch'...");
  Serial.println("...then use 'relay1, mode, <mode>' to choose the timer mode. ");
  Serial.println("Timer modes are: 'countdown', 'expire', 'alternate', 'pulsed', 'single', 'on' and 'off'"); 
  Serial.println("");
  Serial.println("To set the interval or pulse length use, e.g., 'relay1, interval, 30'"); 
  Serial.println("You can also add interval (and optionally the pulse length) when setting timer mode, "); 
  Serial.println("e.g., 'relay1, timer, alternate, 10' or 'relay2, timer, pulsed, 300, 2'"); 
  Serial.println(""); 
  Serial.println("You can also use 'button' to turn on the lights for some interval"); 
  Serial.println("(defined in relay1, defaults to 5 minutes)"); 
  Serial.println("");
  Serial.println("Use 'measure' to get the pressure/temperature/humidity from the sensor."); 
  Serial.println("");
  Serial.println("");
  
  // ... fill this! 
  
}

void relay1func(char *arg){  relay1.parse(arg); }
void relay2func(char *arg){  relay2.parse(arg); }
void relay3func(char *arg){  relay3.parse(arg); }
void relay4func(char *arg){  relay4.parse(arg); }
void relay5func(char *arg){  relay5.parse(arg); }
void relay6func(char *arg){  relay6.parse(arg); }
void relay7func(char *arg){  relay7.parse(arg); }
void relay8func(char *arg){  relay8.parse(arg); }

float average_pressure(){

  float sum=0; 
  
  for(int i=0;i<buf_size;i++){

    sum+=pressure[i];
  
  }

  return sum/buf_size;
  
}

float average_temperature(){
  
  float sum=0; 
  
  for(int i=0;i<buf_size;i++){

    sum+=temperature[i];
  
  }

  return sum/buf_size;
  
}

float average_humidity(){
  

  float sum=0; 
  
  for(int i=0;i<buf_size;i++){

    sum+=humidity[i];
  
  }

  return sum/buf_size;
  
}




