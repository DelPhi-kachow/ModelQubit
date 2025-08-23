/*
modelQubit.ino
Benjamin Krawciw
8/17/2025

This program runs on an arduino microcontroller connected to a circuit which models a qubit. The arduino sets a measurement axis via voltages 
s (sine) and c(cosine) and observes voltages x and y, which will "collapse" to either align or antialign with c and s. The arduino records these
(and other relevant) voltages over time and prints them to the serial monitor, which can be recorded to a csv file. 
*/


// Set pin numbers
const int cPin = A0;
const int sPin = A1;
const int thetPin = A2;
const int phiPin = A3;
const int xPin = A4;
const int yPin = A5;

// Other parameters

// Measurement tracking
float c = 0;
float s = 0;
float theta = 0;
float phi = 0;
float x = 0;
float y = 0;
float t = 0;
bool on = false;

// Set measurement axis
void setAxis_abs(float angle) {
  //In absolute terms, with input as an angle in radians
  c = cos(angle);
  s = sin(angle);
}

void setAxis_rel(float angle) {
  //Set c and s relative to x and y, by angle in radians
  int cPrime = cos(angle);
  int sPrime = sin(angle);

  //angle addition
  c = (x * cPrime) - (y * sPrime);
  s = (y * cPrime) + (x * sPrime);
}

// Switch measurement on/off
void measureOn() {
  // Turn on the measurement coupling
  on = true;
  analogWrite(cPin, int(255 * c));
  analogWrite(sPin, int(255 * s));
}

void measureOff() {
  //turn off the measurement coupling
  on = false;
  analogWrite(cPin, 0);
  analogWrite(sPin, 0);
}

// Record state to variables
void readState() {
  // timer
  t = millis() / 1000.0;
  // analog reads
  theta = analogRead(thetPin) / 1023;
  phi = analogRead(phiPin) / 1023;
  x = analogRead(xPin) / 1023;
  y = analogRead(yPin) / 1023;
}

// Print state to serial
void tellState() {
  // Record state in array
  float state[8] = {t, c, s, on, theta, phi, x, y};

  // Print each element of the array
  for (int i = 0; i < 8; i++) {
    Serial.print(state[i]);
    Serial.print(", ");
  }

  // Take a measurement
  void measure() {
    //Read state
    readState();

    //Report to serial
    tellState();
  }

  //End line
  Serial.println();
}



// Runs once when device turned on
void setup() {
  // Setup pins

  // Begin serial
  Serial.begin(9600);

  // Print metadata
  Serial.println("Voltages are recorded as fractions of 5V")
  Serial.println("t (s), c, s, on/off, theta, phi, x, y");

}

// Runs continuously after setup
void loop() {
  // Cycle through c/s settings

    // At beginning of measurement cycle, change c/s

    // Turn measurement on

    // Turn measurement off

  // All the while, record data

}
