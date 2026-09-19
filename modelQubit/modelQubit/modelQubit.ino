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
const int settleMs = 50; //Wait after switching c/s so the 100k*100nF (10 ms) PWM filter settles (~5 tau)
const int printDigits = 4; //Decimal places in serial output (Serial.print default is 2 -> 0.01 steps)
const long baudRate = 115200; //At 9600 baud one ~60-char line takes ~60 ms to send; set the logger to match
const int N = 10; //Number of phi samples (and trials per angle). Sets precision of measured probabilities.
const int m = 10; //Number of angles to measure probabilities on.

// Phi sweep: PWM 0 → -5 V, PWM 255 → +5 V.
// Endpoints from circuit analysis (not a fit to data):
//   theta circuit: V_theta = (50/47)(1 + u) - V_phi,  u = x c + y s  (fractions of 5 V)
//   The radial circuit settles at r^2 = 14/19 (not 1), so with a unit axis
//   u ranges over [-R_EQ, +R_EQ], R_EQ = sqrt(14/19) ≈ 0.858.
//   theta changes sign at V_phi = (50/47)(1 + u), so the sweep must span
//   V_phi in [(50/47)(1 - R_EQ), (50/47)(1 + R_EQ)] = [0.151, 1.977] V,
//   i.e. phi in [0.03013, 0.39540] of 5 V.
// The THETA CALIBRATION section of mqDataAnalysis.py cross-checks these against
// the measured theta transfer. (PWM resolution is ~39 mV per step.)
const float V_RAIL = 5.0;      // Full-scale op-amp voltage at the PWM extremes (c/s drive)
const float R_EQ = sqrt(14.0 / 19.0);  // Equilibrium radius set by the 1-r^2 circuit
const float THETA_GAIN = 10.0 / 47.0;  // (50/47 V) / V_RAIL: theta-circuit scale in units of 5 V

// Phi drive: 2:1 divider into a unity-gain Sallen-Key low pass (22k, 220nF/470nF).
// V_phi = PHI_GAIN * duty * V_PWM_HIGH, so 0 V .. 2.5 V, one PWM step = 9.8 mV.
// Measure the pin's high level and update V_PWM_HIGH; it is the board's 5 V rail.
const float PHI_GAIN = 0.5;       // Divider gain ahead of the filter
const float V_PWM_HIGH = 5.0;     // Measured PWM high level (V)
const int PHI_SETTLE_MS = 250;    // Filter settles to <0.1 LSB in ~90 ms; 250 ms is margin
// Stored/commanded phi is in units of V_RAIL (same as c/s), so 0 V → 0, 5 V → 1
const float phiMin = THETA_GAIN * (1.0 - R_EQ);  // 0.03013 (0.151 V)
const float phiMax = THETA_GAIN * (1.0 + R_EQ);  // 0.39540 (1.977 V)

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

  //Normalize the current state so (c, s) is a unit vector (5 V axis).
  //Logged x, y sit at |r| ~ 0.87 (14/19 bias), which would otherwise shrink the axis.
  float r = sqrt(x * x + y * y);
  if (r < 0.1) {
    //State too close to the origin to define a direction; fall back to absolute angle
    setAxis_abs(angle);
    return;
  }
  float xn = x / r;
  float yn = y / r;

  //angle addition
  c = (xn * cPrime) - (yn * sPrime);
  s = (yn * cPrime) + (xn * sPrime);
}

// Set phi, stored in units of V_RAIL like c and s, but phi now has its own
// drive circuit: V_phi = PHI_GAIN * (code / 255) * V_PWM_HIGH, so phi >= 0 only.
// The stored value is the quantized one actually delivered, so the log records
// what the circuit did rather than what was requested.
void setPhi(float phiIn) {
  float vWanted = phiIn * V_RAIL;
  int code = int(round(255.0 * vWanted / (PHI_GAIN * V_PWM_HIGH)));
  if (code < 0) code = 0;
  if (code > 255) code = 255;
  analogWrite(phiPin, code);
  phi = (code / 255.0) * PHI_GAIN * V_PWM_HIGH / V_RAIL;
  delay(PHI_SETTLE_MS); // Let the Sallen-Key filter settle before anything is recorded
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
// Phi is parked at 0 V during the collapse, which keeps theta > 0 for every (x, y),
// so the state settles on +(c, s). (The new phi circuit cannot go negative.)
void setXY(float xIn, float yIn) {
  float phiPrev = phi;
  setPhi(0.0);

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

  //Convert raw readings to appropriate coordinates (all stored as V / 5V)
  //A lot of the complexity here is because there is a level shifting circuit to deal with negative voltages
  //theta interface: V_pin = 2V + 0.8*V_theta  ->  V_theta/5V = (thetaRaw - 0.4) / 0.8
  theta = 1.25 * thetaRaw - 0.5;
  //x, y interface: V_pin = 2.5V + 0.5*V_x  ->  V_x/5V = 2*xRaw - 1
  x = 2.0 * xRaw - 1.0;
  y = 2.0 * yRaw - 1.0;
}

// Print state to serial
void tellState() {
  // Record state in array
  float state[8] = {t, c, s, on, theta, phi, x, y};

  // Print each element of the array
  for (int i = 0; i < 8; i++) {
    Serial.print(state[i], printDigits);
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
  delay(settleMs); //Let c/s reach their set values so the first theta reading is valid
  measurements(20); //Record

  // Turn measurement off
  measureOff();
  delay(settleMs);
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
  Serial.begin(baudRate);

  // Print metadata ("fw vN" tells mqDataAnalysis.py which corrections apply)
  // v2: theta = V_theta/5V, 4-decimal, settled reads; v3: unit-length axis, circuit-derived phi sweep
  // v4: dedicated phi drive (gain 0.5 + Sallen-Key), phi logged as the quantized value delivered
  Serial.println("Voltages are recorded as fractions of 5V (fw v4: theta = V_theta/5V, settled reads, unit axis, phi gain 0.5)");
  Serial.println("t (s), c, s, on/off, theta, phi, x, y");

}

// Runs continuously after setup
void loop() {
  fullSuite();
}
