/* This device uses one accelerometers (ADXL345) to make sure the telescope 
 *  doesn't go to a bad direction, in both RA and Dec. 
 *  Please make sure to set transmit=1 and make sure the raw x,y,z match the two projections
 *  Uses code from: https://howtomechatronics.com/tutorials/arduino/how-to-track-orientation-with-arduino-and-adxl345-accelerometer/
 *  2 accelerometer code from: https://forum.arduino.cc/index.php?topic=610345.0 (I couldn't get this to work!)
 *  

*/

#include "Parser.h"
#include "VoltagePin.h"
#include "GroundPin.h"
#include "OutputPin.h"
#include "Timer.h"

#define BUF_SIZE 5

#include<Wire.h>

// pinout MEGA2560: SCL->SCL and SDA->SDA and VCC to 3.3v and GND->GND)
// pinout NANO328:  SCL->A5 (blue) and SDA->A4 (yellow) and VCC to 3.3v and GND->GND)

Parser parser("inclination", "DoubleSwitch v1.00");

bool transmit=0; // do we want a constant serial output with the current data?
float pause=10; // delay between measurements, in milliseconds
unsigned long int anti_flip_delay=10; 

// thresholds for inclination (value of z-component of normalized acceleration vector)
float thresh1=0.05; 
float thresh2=-0.1;

// defining the accelerometers
// #define accel_module1 (0x1d)         // SDO-> Vcc
#define accel_module1 (0x53)         // SDO-> GND

float x, y, z; // acceleration values

// first vector in space to make sure front end is above horizon (must put X direction towards front end!)
float X1=1.00; 
float Y1=0.0;
float Z1=0.0; 

// second vector to make sure counterweight is below horizon 
float X2=0;
float Y2=0;
float Z2=1;

char state1[5] = "up"; // if projection of x,y,z on X1,Y1,Z1 is above threshold
char state2[5] = "up"; // if projection of x,y,z on X2,Y2,Z2 is above threshold

int index=0;
float projections1[BUF_SIZE]={0};
float projections2[BUF_SIZE]={0};

float average1=0;
float average2=0;

VoltagePin led(13); // indication light! 

// relay setup (must connect VCC pin to board 5v pin and GND to GND pin)
VoltagePin relay(A1); // kill the telescope when under horizon OR when counterweight above horizon

void setup() {

  Serial.begin(9600); // Initiate serial communication for printing the results on the Serial monitor

  Parser::debug_bit=0;
  parser.addCommand("transmit", &setTransmission); 
  parser.addCommand("pause", &setPause); 
  parser.addCommand("status", &statusReport);

  relay.setInversion(); // invert the switch so that on means the relay is on

  Wire.begin(); // Initiate the Wire library
  
  // Set ADXL345 in measuring mode (sensor 1)
  Wire.beginTransmission(accel_module1); // Start communicating with the device 
  Wire.write(0x2D); // Access/ talk to POWER_CTL Register - 0x2D
  // Enable measurement
  Wire.write(8); // (8dec -> 0000 1000 binary) Bit D3 High for measuring enable 
  Wire.endTransmission();

}


void serialEvent() {

  parser.readingInput();

}

void loop() {

  parser.checkInput(); // check if any command is issued...

  // === Read acceleromter data === //
  readAcceleration(accel_module1, x, y, z);

  float norm=sqrt(x*x+y*y+z*z); // normalize the current measurements
  float norm1 = sqrt(X1*X1+Y1*Y1+Z1*Z1); // normalize the first vector
  float norm2 = sqrt(X2*X2+Y2*Y2+Z2*Z2); // normalize the second vector

  projections1[index]=(x*X1+y*Y1+z*Z1)/norm/norm1; // dot product of current measurement on the first vector
  projections2[index]=(x*X2+y*Y2+z*Z2)/norm/norm2; // dot product of current measurement on the second vector
  
  index++; // looping index
  if(index>=BUF_SIZE) index=0; // keep track of the last 10 measurements (choose BUF_SIZE)

  average1=0;
  average2=0;
  
  for(int i=0; i<BUF_SIZE;i++){ 
    if(!isnan(projections1[i])) average1+=projections1[i];
    if(!isnan(projections2[i])) average2+=projections2[i];
  }
  average1/=BUF_SIZE;
  average2/=BUF_SIZE;

  if(average1<thresh1) snprintf(state1,5,"down");
  else snprintf(state1,5,"up");
  
  if(average2<thresh2) snprintf(state2,5,"down");
  else snprintf(state2,5,"up");

  if(average1<thresh1 || average2<thresh2){ // one of the accelerometers is pointing down
  
    if(relay.getState()) // only when turning off from the on state
      anti_flip_delay*=2; // each time we kill the current the delay grows larger, so this doesn't endlessly flip the power to the mount (power-off the arduino to reset)

    relay.setState(0); // kill the switch
    led.setState(1); // indicate this on the built-in LED
    delay(anti_flip_delay); // this delay is to prevent flickering: when the acceleration is borderline for turning off the current it might start fliping the switch on/off very fast
    
  }
  else{
    relay.setState(1); // everything is ok, switch can turn on
    led.setState(0); // turn off the LED
  }

  delay(pause); 
  
  if(transmit) statusReport("");

}


void readAcceleration(int address, float &x, float &y, float &z){

  Wire.beginTransmission(address);
  Wire.write(0x32); // Start with register 0x32 (ACCEL_XOUT_H)
  Wire.endTransmission(false);
  Wire.requestFrom(address, 6, true); // Read 6 registers total, each axis value is stored in 2 registers
  x = ( Wire.read()| Wire.read() << 8); // X-axis value
  x = x/256; //For a range of +-2g, we need to divide the raw values by 256, according to the datasheet
  y = ( Wire.read()| Wire.read() << 8); // Y-axis value
  y = y/256;
  z = ( Wire.read()| Wire.read() << 8); // Z-axis value
  z = z/256;
  
}


void statusReport(char *arg){
  
  Serial.print("Raw values: ");
  Serial.print(x);
  Serial.print(", ");
  Serial.print(y);
  Serial.print(", ");
  Serial.print(z);
  
  Serial.print(" | vector1: ");
  Serial.print(X1);
  Serial.print(", ");
  Serial.print(Y1);
  Serial.print(", ");
  Serial.print(Z1);
  
  Serial.print(" | average1= ");
  Serial.print(average1); 
  
  Serial.print(" | state1= ");
  Serial.print(state1);
  
  Serial.print(" | vector2: ");
  Serial.print(X2);
  Serial.print(", ");
  Serial.print(Y2);
  Serial.print(", ");
  Serial.print(Z2);

  Serial.print(" | average2= ");
  Serial.print(average2); 
   
  Serial.print(" | state2= ");
  Serial.println(state2); 

}

void setTransmission(char *arg){

  transmit=parser.parseBool(arg); 
  
}

void setPause(char *arg){

  pause=atoi(arg); 
  
}


