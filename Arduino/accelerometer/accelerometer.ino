/*  *********************************************
    SparkFun_ADXL345_Example
    Triple Axis Accelerometer Breakout - ADXL345
    Hook Up Guide Example

    Utilizing Sparkfun's ADXL345 Library
    Bildr ADXL345 source file modified to support
    both I2C and SPI Communication

    E.Robert @ SparkFun Electronics
    Created: Jul 13, 2016
    Updated: Sep 06, 2016

    Development Environment Specifics:
    Arduino 1.6.11

    Hardware Specifications:
    SparkFun ADXL345
    Arduino Uno

    https://learn.sparkfun.com/tutorials/adxl345-hookup-guide/all#assembly

    
 *  *********************************************/

#include <SparkFun_ADXL345.h>         // SparkFun ADXL345 Library
#include "Parser.h"
#include "VoltagePin.h"
#include "GroundPin.h"
#include "Timer.h"

Parser parser("accelerometer", "accelerometer v1.02");
// v1.02 added option to add timer that sends out status report at set interval

int x, y, z; // acceleration values

/*********** COMMUNICATION SELECTION ***********/
/*    Comment Out The One You Are Not Using    */
//ADXL345 adxl = ADXL345(10);           // USE FOR SPI COMMUNICATION, ADXL345(CS_PIN);
ADXL345 adxl = ADXL345();             // USE FOR I2C COMMUNICATION
VoltagePin adxlVCC(2); // power for the CS pin 
GroundPin adxlGND(3); // ground for the SDO pin
// additional pins on the ADXL chip:
// VCC and GND (in addition to those marked above)
// SDA is connected to A4
// SCL is connected to A5

Timer timer;

void setup() {
  Serial.begin(9600);                 // Start the serial terminal
//  Serial.println("SparkFun ADXL345 Accelerometer Hook Up Guide Example");
//  Serial.println();

  Parser::debug_bit = 0;
  parser.addCommand("status", &statusReport);
  parser.addCommand("timer", &setTimer);
  parser.addCommand("duration", &setTimer);
  parser.addCommand("interval", &setTimer);
  parser.addCommand("period", &setTimer);
  parser.addCommand("help", &help);

  parser.setThisSensorName("accel");

  adxl.powerOn();                     // Power on the ADXL345

  adxl.setRangeSetting(2);           // Give the range settings
  // Accepted values are 2g, 4g, 8g or 16g
  // Higher Values = Wider Measurement Range
  // Lower Values = Greater Sensitivity

  // adxl.setSpiBit(0);                  // Configure the device to be in 4 wire SPI mode when set to '0' or 3 wire SPI mode when set to 1
  // Default: Set to 1
  // SPI pins on the ATMega328: 11, 12 and 13 as reference in SPI Library

  timer.setInterval(1);
  timer.setMode(Timer::OFF);

}

void serialEvent() {

  parser.readingInput();

}

/****************** MAIN CODE ******************/
/*     Accelerometer Readings and Interrupt    */
void loop() {

  // Accelerometer Readings
  adxl.readAccel(&x, &y, &z);         // Read the accelerometer values and store them in variables declared above x,y,z  
  
  parser.checkInput(); // check if any command is issued...

  if(timer.getMode()==Timer::EXPIRE && timer.getState()==0){
    statusReport();
    timer.setupExpiration();
  }

}

void statusReport() {
  statusReport("");
}

void statusReport(char *arg) {

  Serial.print(x);
  Serial.print(", ");
  Serial.print(y);
  Serial.print(", ");
  Serial.println(z);
  // Serial.println(timer.getInterval()); // for debugging only

}

void setTimer(char *arg){
  
  char str1[Parser::STRN];
  char str2[Parser::STRN];

  Parser::splitStr(arg, str1, str2);
  Parser::cleanString(str1);
  if(strlen(str1)>0 && atof(str1)>0){ // legal timer interval
    timer.parseTiming(str1);
    timer.setupExpiration();

    // debugging outputs:
    Serial.print("Interval set to ");
    Serial.print(timer.getInterval());
    Serial.println(" seconds.");
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
    
    // debugging output:
    Serial.println("Timer is off");
  }
  
  
}

void help(char *arg){

  parser.getVersion();
  Serial.println(" Use 'status' to get updated x,y,z position");
  Serial.println(" Use 'timer' or 'interval' or 'duration' to set timer");
  Serial.println(" When timer is set to non-zero T, will ouput status automatically every T");
  Serial.println(" Set timer to zero to turn it off");
  Serial.println(" Use 'timer' (or other commands) without argument to get the timer interval");
  
}


