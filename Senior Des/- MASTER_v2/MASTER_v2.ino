/*
  -------------------------------------------------------
                        MASTER
  -------------------------------------------------------

  Functions
  ---------------------------------- Weather Sensing
  SparkFun Weather Shield
        - Wind Speed | ADA1733
        - Humidity | Si7021
        - Temperature | Si7021
        - Light Level | ALS-PT19
        - Atmospheric Pressure | MPL3115A2
        - GPS | GP635T
  ---------------------------------- Data Reading/Logging
  SparkFun microSD Shield

  ----------------------------------

  NOTES
  - upload first, then connect battery
*/
#include <SPI.h>                              // SPI for SD card               - Pre-installed in Library Manager
#include <SD.h>                               // SD card                       - Pre-installed in Library Manager
#include <Wire.h>                             // I2C needed for sensors        - Pre-installed in Library Manager
#include <SoftwareSerial.h>                   // Serial communication for GPS  - Pre-installed in Library Manager
//#include <DualTB9051FTGMotorShield.h>         // Motor Shield                  - Search "TB9051" and install from Library Manager
#include <TinyGPS++.h>                        // GPS                           - Included in .libraries 
#include "SparkFunMPL3115A2.h"                // Pressure sensor               - Search "SparkFun MPL3115" and install from Library Manager
#include "SparkFun_Si7021_Breakout_Library.h" // Humidity sensor               - Search "SparkFun Si7021" and install from Library Manager

MPL3115A2 myPressure;               // Create an instance of the pressure sensor
Weather myHumidity;                 // Create an instance of the humidity sensor
TinyGPSPlus gps;                    // Creates instance of gps
//DualTB9051FTGMotorShield md; // Create instance of trap motor

static const int RXPin = 5, TXPin = 4; //GPS is attached to pin 4(TX from GPS) and pin 5(RX into GPS)
SoftwareSerial ss(RXPin, TXPin);

// Hardware Pins
//------------------------------------------------------------
const byte STAT_BLUE = 7;
const byte STAT_GREEN = 8;
const byte GPS_PWRCTL = 6; // Pulling pin puts gps to lower power mode
const byte chipSelect = 8;

const byte buttonPin = 9;
const byte LIGHT = A1;
const byte WIND = A2;
const byte REFERENCE_3V3 = A3;


//  Global Variables
//------------------------------------------------------------
long lastSecond;                    //The millis counter to see when a second rolls by
float resolution = 5.0 / 1024.0;

// Weather Variables
//-------------------------------------------------------------
float humidity;
float temp_h;
float pressure;
byte years;
byte hours;
byte months;
byte days;
byte minutes;
byte seconds;

//-------------------------------------- ANEMOMETER SETUP
float wind_speed;                 // Wind speed variable
float windVolts;                  // Wind voltage
const float voltageMin = 0.45;    // Lowest voltage output

void setup()
{
  Serial.begin(9600); // Begin serial communication

 // md.init();
  
  ss.begin(9600);
  SD.begin();
  //md.enableDrivers();

  pinMode(buttonPin, INPUT);
  pinMode(chipSelect, OUTPUT); //
  pinMode(STAT_BLUE, OUTPUT); // Status LED Blue
  pinMode(STAT_GREEN, OUTPUT); // Status LED Green

  digitalWrite(STAT_BLUE, HIGH);

  File input = SD.open("input.csv");
  if (input) {
    // read from the file until there's nothing else in it:
    byte times[10];
    int i = 0;
    while (input.available()) {
      times[i] = input.read();
      Serial.write(times[i]);
      i = i + 1;
    }
    Serial.write(1);
    // close the file:
    input.close();
  } else {
    // if the file didn't open, print an error:
    Serial.println("FOCK");
  }

  if (SD.exists("weather.csv")) {
    SD.remove("weather.csv");
  }

  // Grab Date
  months = gps.date.month();
  days = gps.date.day();
  years = gps.date.year();

  String date1 = String(months);
  String date2 = String(days);
  String date3 = String(years);

  String date_string = (date1 + "-" + date2 + "-" + date3);

  // Print CSV Header
  File dataFile = SD.open("weather.csv", FILE_WRITE);
  dataFile.println(date_string);
  dataFile.print("Pressure [Pa]");
  dataFile.print(",");
  dataFile.print("Temperature [F]");
  dataFile.print(",");
  dataFile.print("Humidity [%]");
  dataFile.print(",");
  dataFile.print("Wind Speed [m/s]");
  dataFile.print(",");
  dataFile.print("Time,");
  dataFile.println();
  dataFile.close();


  // Sensor Configuration
  //-----------------------------------------------------------------------------------
  // Configure pressure sensor
  myPressure.begin();              // Get sensor online
  myPressure.setModeBarometer();   // Measure pressure in Pascals from 20 to 110 kPa
  myPressure.setOversampleRate(7); // Set Oversample to the recommended 128
  myPressure.enableEventFlags();  // Enable all three pressure and temp event flags

  // Configure humidity sensor
  myHumidity.begin();

  // Get humidity and temperature
  humidity = myHumidity.getRH();
  temp_h = myHumidity.readTemp();

  // Get Pressure
  pressure = myPressure.readPressure();
}

void loop()
{
  
//      delay(1); // wait for drivers to be enabled so fault pins are no longer low
//  //md.setM2Speed(400);
//  int buttonVoltage = digitalRead(buttonPin);//*(5.0/1023.0);
//  Serial.println(buttonVoltage);
//  md.setM2Speed(400);
//  while (buttonVoltage == 1) {
//    delay(265);        //delay to allow switch to open again
//    md.setM2Speed(0);
//    delay(2000);
//    break;
//    delay(5);
//}
  // Collect Data
  //----------------------------------------------------------------------------------
  // Grab weather values
  pressure = myPressure.readPressure(); // Pressure
  temp_h = myHumidity.readTempF(); // Temperature from Humidity sensor
  humidity = myHumidity.getRH();
  wind_speed = get_wind_speed();

  // Grab Time
  hours = gps.time.hour();
  minutes = gps.time.minute();
  seconds = gps.time.second();

  // Print the values to the serial monitor
//  Serial.print(pressure);
//  Serial.print(" Pa");
//  Serial.print("\t");
//  Serial.print(temp_h);
//  Serial.print(" F");
//  Serial.print("\t");
//  Serial.print(humidity);
//  Serial.print("%");
//  Serial.print("\t");
//  Serial.print(wind_speed);
//  Serial.print("m/s");
//  Serial.print("\t");
//  Serial.print(hours);
//  Serial.print(":");
//  Serial.print(minutes);
//  Serial.print(":");
//  Serial.print(seconds);
//  Serial.println();

  // Write to SD card
  //---------------------------------------------------------------------------------
  File dataFile = SD.open("weather.csv", FILE_WRITE);
  if (dataFile) {
    dataFile.print(String(pressure, 2));
    dataFile.print(",");
    dataFile.print(String(temp_h, 2));
    dataFile.print(",");
    dataFile.print(String(humidity, 2));
    dataFile.print(",");
    dataFile.print(String(wind_speed, 2));
    dataFile.print(",");
    dataFile.print(String(hours));
    dataFile.print(":");
    dataFile.print(String(minutes));
    dataFile.print(":");
    dataFile.print(String(seconds));
    dataFile.print(",");
    dataFile.println();
    dataFile.close();
  }
  else
  {
    Serial.println("error opening csv");
  }
}

//Returns the voltage of the light sensor based on the 3.3V rail
//This allows us to ignore what VCC might be (an Arduino plugged into USB has VCC of 4.5 to 5.2V)
float get_light_level()
{
  float operatingVoltage = analogRead(REFERENCE_3V3);
  float lightSensor = analogRead(LIGHT);

  operatingVoltage = 3.3 / operatingVoltage; //The reference voltage is 3.3V

  lightSensor = operatingVoltage * lightSensor;

  return (lightSensor);
}

//Returns the voltage of the raw pin based on the 3.3V rail
//This allows us to ignore what VCC might be (an Arduino plugged into USB has VCC of 4.5 to 5.2V)
//Battery level is connected to the RAW pin on Arduino and is fed through two 5% resistors:
//3.9K on the high side (R1), and 1K on the low side (R2)
float get_wind_speed()
{
  windVolts = analogRead(WIND) * resolution;

  if (windVolts > voltageMin) {
    wind_speed = 1.2821 * exp(0.328 * windVolts) ;
  }
  else {
    wind_speed = 0.00;
  }

  delay(100);

  return (wind_speed);
}


    

   
 
