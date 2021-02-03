int buttonPin = 9;
int buttonState = 0;

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  pinMode(buttonPin, INPUT);
}

void loop() {
  // put your main code here, to run repeatedly:
  int buttonVoltage = digitalRead(buttonPin); //* (5.0/1023.0);
//  delay(1000);
  Serial.println(buttonVoltage);
}
