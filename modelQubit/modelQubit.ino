/*
modelQubit.ino
Benjamin Krawciw
8/17/2025

This program runs on an arduino microcontroller connected to a circuit which models a qubit. The arduino sets a measurement axis via voltages 
s (sine) and c(cosine) and observes voltages x and y, which will "collapse" to either align or antialign with c and s. The arduino records these
(and other relevant) voltages over time and prints them to the serial monitor, which can be recorded to a csv file. 
*/


// Set pin numbers
const int cPin = 6;
const int sPin = 5;
const int thetPin = A1;
const int phiPin = A0;
const int xPin = A2;
const int yPin = A3;

// Other parameters
const int dt = 100; //The number of milliseconds between consecutive readings

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

  // The scaling is a bit weird because we need to output negative voltages too
  int cOut = int( 255.0 * ((c + 1.0) / 2.0));
  int sOut = int( 255.0 * ((s + 1.0) / 2.0));
  analogWrite(cPin, cOut);
  analogWrite(sPin, sOut);
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
  // analog reads as fraction of max voltage
  float thetaRaw = analogRead(thetPin) / 1023.0;
  float phiRaw = analogRead(phiPin) / 1023.0;
  float xRaw = analogRead(xPin) / 1023.0;
  float yRaw = analogRead(yPin) / 1023.0;

  //Convert raw readings to appropriate coordinates
  //A lot of the complexity here is because there is a level shifting circuit to deal with negative voltages
  theta = 2.0 * thetaRaw - 1.0;
  phi = phiRaw;
  x = 2.0 * xRaw - 1.0;
  y = 2.0 * yRaw - 1.0;
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
  //End line
  Serial.println();
}

// Take a measurement
void measure() {
  //Read state
  readState();

  //Report to serial
  tellState();
}

// Take a suite of measurements at regular intervals
void measurements(int N) {
  for (int i = 0; i < N; i++) {
    measure();
    delay(dt);
  }
}

// Runs once when device turned on
void setup() {
  // Setup pins

  // Begin serial
  Serial.begin(9600);

  // Print metadata
  Serial.println("Voltages are recorded as fractions of 5V");
  Serial.println("t (s), c, s, on/off, theta, phi, x, y");

}

// Runs continuously after setup
void loop() {
  // Cycle through c/s settings

    // At beginning of measurement cycle, change c/s
    setAxis_rel(PI / 2.0);

    // Turn measurement on
    measureOn();
    measurements(20); //Record

    // Turn measurement off
    measureOff(); 
    measurements(20); //Record

}
