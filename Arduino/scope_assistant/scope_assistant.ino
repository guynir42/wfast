/***************************************
 * This device needs to be hooked up to the front of the telescope plate. 
 * It has 3 functions:
 * 1) Measure acceleration to give the orientation of the telescope (and measure strong vibrations). 
 * 2) MEasure the distance to nearest object in front of telescope using ultrasonic sensor. 
 * 3) Turn on/off some devices on top of the telescope (e.g., heaters, lens cap). 
 * 
 * Usage / control: 
 * 
 */


#include "Parser.h"
#include "VoltagePin.h"
#include "GroundPin.h"
#include "OutputPin.h"
#include "Timer.h"

#include "Ultrasonic.h"
#include <SparkFun_ADXL345.h>         // SparkFun ADXL345 Library

Parser parser("scope", "ScopeAssistant v1.01");
// v1.01 -- now ouput timer interval as 5th output of "measure"


// define the accelerometer 
// pinout MEGA2560: SCL->SCL and SDA->SDA and VCC to 3.3v and GND->GND)
// pinout NANO328:  SCL->A5 and SDA->A4 and VCC to 3.3v and GND->GND)
VoltagePin accVCC(19);
GroundPin accGND(18);
ADXL345 adxl = ADXL345();             // USE FOR I2C COMMUNICATION
int x, y, z; // acceleration values

// define the ultrasonic sensor and its voltage/ground pins
VoltagePin ultraVCC(43);
GroundPin ultraGND(37);
Ultrasonic ultra(41,39);

Timer timer; // timer if off by default, can activate it to get periodic status reports via serial

// relay setup (must connect VCC pin to board 5v pin and GND to GND pin)
OutputPin relay1(10); // this relay will control the heater
OutputPin relay2(11); // maybe open lid?
OutputPin relay3(12); // maybe close lid?
OutputPin relay4(13); // what else?


void setup() {
  
  Serial.begin(9600);                 // Start the serial terminal

  Parser::debug_bit=0;
  Timer::debug_bit=0;
  parser.addCommand("measure", &measure);
  parser.addCommand("status", &statusReport);
  parser.addCommand("timer", &setTimer);
  parser.addCommand("duration", &setTimer);
  parser.addCommand("interval", &setTimer);
  parser.addCommand("period", &setTimer);
  parser.addCommand("heater", &heater);
  parser.addCommand("help", &help);

  timer.setInterval(1);
  timer.setMode(Timer::OFF);

  adxl.powerOn();                     // Power on the ADXL345

  adxl.setRangeSetting(2);           // Give the range settings
  // Accepted values are 2g, 4g, 8g or 16g
  // Higher Values = Wider Measurement Range
  // Lower Values = Greater Sensitivity

  relay1.timer.setMode(Timer::OFF);
  relay2.timer.setMode(Timer::OFF);
  relay3.timer.setMode(Timer::OFF);
  relay4.timer.setMode(Timer::OFF);

  relay1.setMode(OutputPin::WATCH);
  relay2.setMode(OutputPin::WATCH);
  relay3.setMode(OutputPin::WATCH);
  relay4.setMode(OutputPin::WATCH);

  relay1.setInverse();
  relay2.setInverse();
  relay3.setInverse();
  relay4.setInverse(0); // make sure the LED is off by default

}

void serialEvent() {

  parser.readingInput();

}

/****************** MAIN CODE ******************/
/*     Accelerometer Readings and Interrupt    */
void loop() {

  adxl.readAccel(&x, &y, &z);         // Read the accelerometer values and store them in variables declared above x,y,z  
  ultra.measure();         // Read the ultrasonic sensor 
  
  parser.checkInput(); // check if any command is issued...

  if(timer.getMode()==Timer::EXPIRE && timer.getState()==0){ // if timer triggers, send status report and reset it
    measure("");
    timer.setupExpiration();
  }

  relay1.update();
  relay2.update();
  relay3.update();
  relay4.update();

}


void measure(char *arg){

  Serial.print(x);
  Serial.print(", ");
  Serial.print(y);
  Serial.print(", ");
  Serial.print(z);
  Serial.print(", ");
  Serial.print(ultra.getDistance());
  Serial.print(", ");
  if(timer.getMode()==Timer::OFF) Serial.println("0");
  else Serial.println(timer.getInterval());
  
}

void statusReport(char *arg) {

  Serial.print("Acceleration: ");
  Serial.print(x);
  Serial.print(", ");
  Serial.print(y);
  Serial.print(", ");
  Serial.println(z);
  
  Serial.print("Distance: ");
  Serial.println(ultra.getDistance());

  Serial.print("Relay1 (heater): ");
  Serial.println(relay1.timer.printout());

  Serial.print("Feedback Timer: ");
  Serial.println(timer.printout());
  
}

void setTimer(char *arg){
  
  char str1[Parser::STRN];
  char str2[Parser::STRN];

  Parser::splitStr(arg, str1, str2);
  Parser::cleanString(str1);
  if(strlen(str1)>0 && atof(str1)>0){ // legal timer interval
    timer.parseTiming(str1);
    timer.setupExpiration();

  }
  else if(strlen(str1)==0){ // only output the interval
    Serial.print("Timer interval is: ");
    Serial.println(timer.getInterval());
    if(timer.getMode()==Timer::OFF){
      Serial.println("Timer is off");
    }
  }
  else{ // turn off interval
    timer.setMode(Timer::OFF);
  }
  
  
}

void heater(char *arg){

  relay1.timer.parse(arg);
  Serial.println("updating heater timer...");
  Serial.println(relay1.timer.printout());
  
}

void help(char *arg){

  parser.getVersion();
  Serial.println(" Use 'status' to get updated a_x, a_y, a_z, dist measurements. ");
  Serial.println(" Use 'timer' or 'interval' or 'duration' to set timer");
  Serial.println(" When timer is set to non-zero T, will ouput status automatically every T");
  Serial.println(" Set timer to zero to turn it off");
  Serial.println(" Use 'timer' (or other commands) without argument to get the timer interval");
  
}


