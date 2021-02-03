// includes
#include <DualTB9051FTGMotorShield.h>

DualTB9051FTGMotorShield md;

// declares
const int buttonPin = 9;
int buttonState = 0;

void setup()
{
  Serial.begin(9600);

  // initialization
  pinMode(buttonPin, INPUT);
  md.init();

  // Uncomment to flip a motor's direction:
  //md.flipM1(true);
  //md.flipM2(true);
  //md.setM2Speed(200);
}

void loop()
{
  md.enableDrivers();
  delay(1); // wait for drivers to be enabled so fault pins are no longer low
  //md.setM2Speed(400);
  int buttonVoltage = digitalRead(buttonPin);//*(5.0/1023.0);
  Serial.println(buttonVoltage);
  md.setM2Speed(400);
  while (buttonVoltage == 1) {
    Serial.println("Button OPEN ");
//    md.setM2Speed(0);
//    delay(5000);        //delay for interval wait (taking measurements)
//    md.setM2Speed(400);
    delay(265);        //delay to allow switch to open again
    md.setM2Speed(0);
    delay(2000);
    break;
    delay(5);
  }

//  Serial.print("M1 current: ");
//  Serial.println(buttonVoltage);

  md.disableDrivers();
}
