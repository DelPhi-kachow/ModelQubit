/*
modelQubit.ino
Benjamin Krawciw
8/17/2025

This program runs on an arduino microcontroller connected to a circuit which models a qubit. The arduino sets a measurement axis via voltages 
s (sine) and c(cosine) and observes voltages x and y, which will "collapse" to either align or antialign with c and s. The arduino also sets 
the phi voltage (same drive circuit as c/s). It records these (and other relevant) voltages over time and prints them to the serial monitor, 
which can be recorded to a csv file. 
*/


// Set pin numbers
const int cPin = 6;
const int sPin = 5;
const int phiPin = 3; // PWM output (same drive circuit as c/s)
const int thetPin = A1;
const int xPin = A2;
const int yPin = A3;

// Other parameters
const int dt = 50; //The number of milliseconds between consecutive readings
const int N = 10; //Number of phi samples (and trials per angle). Sets precision of measured probabilities.
const int m = 10; //Number of angles to measure probabilities on.

// Phi sweep: PWM 0 → -5 V, PWM 255 → +5 V. We only use the old triangle range.
const float V_RAIL = 5.0;     // Full-scale op-amp voltage at the PWM extremes
const float V_PHI_LO = 0.0;   // Intended sweep low (V)
const float V_PHI_HI = 1.848; // Intended sweep high (V)
const float V_CLIP = 0.277;   // Clip both ends to avoid unphysical levels
// Stored/commanded phi is in units of V_RAIL (same as c/s), so 0 V → 0, 5 V → 1
const float phiMin = (V_PHI_LO + V_CLIP) / V_RAIL;
const float phiMax = (V_PHI_HI - V_CLIP) / V_RAIL;

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
  float cPrime = cos(angle);
  float sPrime = sin(angle);

  //angle addition
  c = (x * cPrime) - (y * sPrime);
  s = (y * cPrime) + (x * sPrime);
}

// Set phi (normalized to ±V_RAIL, same scaling as c/s: -1 → -5 V, +1 → +5 V)
void setPhi(float phiIn) {
  phi = phiIn;
  int phiOut = int(255.0 * ((phi + 1.0) / 2.0));
  analogWrite(phiPin, phiOut);
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
  analogWrite(cPin, 128);
  analogWrite(sPin, 128);
}

// Set (x, y) indirectly by measuring along (xIn, yIn), waiting to equilibrate, then releasing.
// Phi is parked slightly negative during the collapse so the state can settle, then restored.
void setXY(float xIn, float yIn) {
  float phiPrev = phi;
  setPhi(-V_CLIP / V_RAIL); // slightly below 0 V

  c = xIn;
  s = yIn;
  measureOn();
  delay(1000); // Fixed settle time; could later wait until (x, y) is close enough
  c = 0.0;
  s = 0.0;
  measureOff();

  setPhi(phiPrev);
}

// Record state to variables
void readState() {
  // timer
  t = millis() / 1000.0;
  // analog reads as fraction of max voltage
  // phi is commanded by setPhi (no longer read back)
  float thetaRaw = analogRead(thetPin) / 1023.0;
  float xRaw = analogRead(xPin) / 1023.0;
  float yRaw = analogRead(yPin) / 1023.0;

  //Convert raw readings to appropriate coordinates
  //A lot of the complexity here is because there is a level shifting circuit to deal with negative voltages
  theta = 2.0 * thetaRaw - 1.0;
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

// Perform an individual trial for a parameter set
void trial() {
  // Turn measurement on
  measureOn();
  measurements(20); //Record

  // Turn measurement off
  measureOff(); 
  measurements(20); //Record
}

// Perform trials for each angle
void angleSuite() {
  //Iterate through relative angles
  for (int i = 0; i < m; i++){
    float thet = i * ((PI / 2.0) / (m - 1.0));
    setAxis_rel(thet);
    trial();
  }
}

// Perform the full suite: N phi samples × m relative angles
void fullSuite(){
  for (int i = 0; i < N; i++) {
    float phiVal = phiMin + i * (phiMax - phiMin) / (N - 1.0);
    setPhi(phiVal);
    angleSuite();
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
  fullSuite();
}
