/*
  created 24 Nov 2010
  modified 9 Apr 2012
  by Tom Igoe

  modified 18 Sep 2014
  by Bobby Chan @ SparkFun Electronics Inc.

  SD Card Datalogger

  This example is based off an example code from Arduino's site
  http://arduino.cc/en/Tutorial/Datalogger and it shows how to
  log data from three analog sensors with a timestamp based on when
  the Arduino began running the current program to an SD card using
  the SD library https://github.com/greiman/SdFat by William
  Greiman. This example code also includes an output to the
  Serial Monitor for debugging.

  The circuit:
   analog sensors on analog pins 0, 1, and 2
   SD card attached to SPI bus as follows:
 ** MOSI - pin 11
 ** MISO - pin 12
 ** CLK - pin 13
 ** CS - pin 4

  This example code is in the public domain.
*/

#include <SD.h>

// On the Ethernet Shield, CS is pin 4. Note that even if it's not
// used as the CS pin, the hardware CS pin (10 on most Arduino boards,
// 53 on the Mega) must be left as an output or the SD library
// functions will not work.

// Chip Select pin is tied to pin 8 on the SparkFun SD Card Shield
const int chipSelect = 8;

void setup()
{
  // Open serial communications and wait for port to open:
  Serial.begin(9600);
  Serial.print("Initializing SD card...");

  // make sure that the default chip select pin is set to
  // output, even if you don't use it:
  pinMode(chipSelect, OUTPUT);

  // see if the card is present and can be initialized:
  if (!SD.begin(chipSelect)) {
    Serial.println("Card failed, or not present");
    // don't do anything more:
    return;
  }
  Serial.println("card initialized.");

  File input = SD.open("input.csv");
  if (input) {
    // read from the file until there's nothing else in it:
    int i = 0;          // counter
    int nonsense = 100;
    int temp_int;
    int CR = 10;
    int LINE = 13;
    char temp_array[35];
    byte temp;

    while (input.available()) {
      temp = input.read();
      temp_int = int(temp);
      // Serial.print("Temp: ");
      //Serial.println(temp_int);
      if (temp_int < nonsense & temp_int != CR & temp_int != LINE) {
        temp_array[i] = temp_int;
        //Serial.println(temp_array[i]);
        i = i + 1;
      }
    }

    String time1;
    String time2;
    String times[15];
    int j = 0;
    int k = 0;

//    for (int j = 0; j < i; j++){
//      if ((j+1)% 2 != 0){
//        time1 = String(temp_array[j]);
//      }
//      else{
//        time2 = String(temp_array[j]);
//      }
//    }

    while (j < i){
      while (k < 2){
        time1 = temp_array[j];
        j = j + 1;
        k = k + 1;
        time2 = temp_array[j];
        k = k + 1;
        times[j] = time1 + time2;
        Serial.println(times[j]);
      }
      k = 0;
      j = j + 1;
    }
    
    Serial.print("length array: ");
    Serial.println(i);
    Serial.println("----------------------------");

    char a = 65;
    char b = 66;
    Serial.print(String(a) + String(b));
    // close the file:
    input.close();
  } else {
    // if the file didn't open, print an error:
    Serial.println("FOCK");
  }
}

void loop()
{
  //  byte a = 1;
  //  for (int i = 0; i <= 10; i++){
  //    Serial.println(times[i]);
  //  }
}
