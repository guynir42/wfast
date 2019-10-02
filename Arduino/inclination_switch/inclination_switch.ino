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

Parser parser("inclination", "InclinationSwitch v1.00");
// v1.00 -- hard limit switch by measuring inclination 


// define the accelerometer 
// pinout MEGA2560: SCL->SCL and SDA->SDA and VCC to 3.3v and GND->GND)
// pinout NANO328:  SCL->A5 and SDA->A4 and VCC to 3.3v and GND->GND)
VoltagePin accVCC(19);
GroundPin accGND(18);
ADXL345 adxl = ADXL345();             // USE FOR I2C COMMUNICATION
int x, y, z; // acceleration values

Timer timer; // timer if off by default, can activate it to get periodic status reports via serial

// relay setup (must connect VCC pin to board 5v pin and GND to GND pin)
VoltagePin relay1(A1); // kill the telescope when under horizon
VoltagePin relay2(A2); // add later
VoltagePin relay3(A3); // add later
VoltagePin relay4(A4); // add later

char state[5] = "up";

int index=0;
float projections[10]={0};

void setup() {
  
  Serial.begin(9600);                 // Start the serial terminal

  Parser::debug_bit=0;
  Timer::debug_bit=0;
//  parser.addCommand("measure", &measure);
  parser.addCommand("status", &statusReport);
  parser.addCommand("timer", &setTimer);
  parser.addCommand("duration", &setTimer);
  parser.addCommand("interval", &setTimer);
  parser.addCommand("period", &setTimer);
//  parser.addCommand("help", &help);

  timer.setInterval(0.01);
  timer.setMode(Timer::EXPIRE);

  adxl.powerOn();                     // Power on the ADXL345

  adxl.setRangeSetting(2);           // Give the range settings
  // Accepted values are 2g, 4g, 8g or 16g
  // Higher Values = Wider Measurement Range
  // Lower Values = Greater Sensitivity

  // relay1.timer.setMode(Timer::ON);
  // relay2.timer.setMode(Timer::OFF);
  // relay3.timer.setMode(Timer::OFF);
  // relay4.timer.setMode(Timer::OFF);

  // relay1.setMode(OutputPin::WATCH);
  // relay2.setMode(OutputPin::WATCH);
  // relay3.setMode(OutputPin::WATCH);
  // relay4.setMode(OutputPin::WATCH);

  relay1.setInversion();
  relay2.setInversion();
  relay3.setInversion();
  relay4.setInversion();

}

void serialEvent() {

  parser.readingInput();

}

/****************** MAIN CODE ******************/
/*     Accelerometer Readings and Interrupt    */
void loop() {

  adxl.readAccel(&x, &y, &z);         // Read the accelerometer values and store them in variables declared above x,y,z  
  
  parser.checkInput(); // check if any command is issued...

  if(timer.getMode()==Timer::EXPIRE && timer.getState()==0){ // if timer triggers, send status report and reset it

    float X,Y,Z;
    X=x;
    Y=y;
    Z=z;

    float norm = sqrt(X*X+Y*Y+Z*Z); // normalize the acceleration vector
    
    projections[index]=Z/norm; // projection of normalized vector on down direction
    index++;
    if(index>=10) index=0; // keep track of the last 10 measurements

    // calculate the mean projection
    float average=0;
    for(int i=0; i<10;i++) if(!isnan(projections[i])) average+=projections[i];
    average/=10;
    
    if(average<0){
      relay1.setState(0); // vector is pointing down, kill the switch
      snprintf(state,5,"down");
    }
    else{
      relay1.setState(1); // everything is ok, switch can turn on
      snprintf(state,5,"up");
    }
    
    timer.setupExpiration(); // reset the timer

    statusReport("");
    
  }

}


//void measure(char *arg){
//
//  Serial.print(x);
//  Serial.print(", ");
//  Serial.print(y);
//  Serial.print(", ");
//  Serial.print(z);
//  Serial.print(", ");
//  Serial.print(state);
//  Serial.print(", t= ");
//  if(timer.getMode()==Timer::OFF) Serial.println("0");
//  else Serial.println(timer.getInterval());
//  
//}

void statusReport(char *arg) {

  Serial.print("Acc: ");
  Serial.print(x);
  Serial.print(", ");
  Serial.print(y);
  Serial.print(", ");
  Serial.print(z);
  Serial.print(", projections: ");
  for(int i=0;i<10;i++){
    Serial.print(projections[i]);
    Serial.print(", ");
  }
  
  Serial.print(state);
  
  Serial.print(", Timer: ");
  Serial.println(timer.getInterval()); 
//  Serial.println(timer.printout());
  
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

//void help(char *arg){
//
//  parser.getVersion();
//  Serial.println(" Use 'measure' to get updated a_x, a_y, a_z. ");
//  Serial.println(" Use 'timer' or 'interval' or 'duration' to set timer.");
//  Serial.println(" Must set a non-zero timer value for limit switch to work!");
//  Serial.println(" Set timer to zero to turn it off (limit disabled).");
//  Serial.println(" Use 'timer' (or other commands) without argument to get the timer interval");
//  
//}


