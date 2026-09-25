/*
circuitDebug.ino
Benjamin Krawciw

Interactive serial console for debugging the modelQubit circuit.
Monitors the same voltages as the trial suite, and accepts commands to
set x, y, c, s, and phi from the Serial Monitor.

Values are fractions of 5 V (same convention as modelQubit.ino):
  -1 → -5 V,  0 → 0 V,  +1 → +5 V
For c and s that is the full ±5 V range. Phi has its own drive circuit
(2:1 divider into a unity-gain Sallen-Key filter), so phi covers 0 .. +0.5
only (0 .. 2.5 V) in steps of 9.8 mV; values outside that range are clipped.

Commands (one per line; case-insensitive names):
  help              Print this help
  status            Print one reading now
  stream on|off     Continuously print readings (default: on)
  c <val>           Set c and drive the c pin
  s <val>           Set s and drive the s pin
  cs <c> <s>        Set both c and s
  phi <val>         Set phi and drive the phi pin (0 .. 0.5, i.e. 0 .. 2.5 V)
  x <val>           Collapse toward (val, last y); see note below
  y <val>           Collapse toward (last x, val)
  xy <x> <y>        Collapse toward (x, y)
  on                Drive current c/s onto the pins
  off               Park c/s at midscale (~0 V)

Note: x and y have no DAC. Setting them measures along that axis for 1 s
(with phi briefly parked slightly negative), then releases measurement.
*/


// Pin numbers (match modelQubit.ino)
const int cPin = 9;  // Timer1 PWM at 31 kHz (see setupCsPwm)
const int sPin = 10;
const int phiPin = 3;
const int thetPin = A1;
const int xPin = A2;
const int yPin = A3;

const int dt = 50; // ms between streamed readings
// A streamed line is ~73 chars: 76 ms at 9600 baud (longer than dt), 6.3 ms at 115200.
// Must match modelQubit.ino and the Serial Monitor setting.
const long baudRate = 115200;
const float V_RAIL = 5.0;

// Phi drive (must match modelQubit.ino): V_phi = PHI_GAIN * duty * V_PWM_HIGH
const float PHI_GAIN = 0.5;       // Divider gain ahead of the Sallen-Key filter
const float V_PWM_HIGH = 5.0;     // Measured PWM high level (V)
const int PHI_SETTLE_MS = 250;    // Filter settles to <0.1 LSB in ~90 ms; 250 ms is margin

float c = 0;
float s = 0;
float theta = 0;
float phi = 0;
float x = 0;
float y = 0;
float t = 0;
bool on = false;

// Last commanded (x, y) targets for single-axis updates
float xCmd = 0;
float yCmd = 0;

bool streaming = true;
unsigned long lastStreamMs = 0;

String lineBuf;

int normToPwm(float v) {
  return int(255.0 * ((v + 1.0) / 2.0));
}

// Phi has its own drive circuit: V_phi = PHI_GAIN * (code/255) * V_PWM_HIGH.
// Stores the quantized value actually delivered, and waits for the filter.
void setPhi(float phiIn) {
  float vWanted = phiIn * V_RAIL;
  int code = int(round(255.0 * vWanted / (PHI_GAIN * V_PWM_HIGH)));
  if (code < 0) code = 0;
  if (code > 255) code = 255;
  analogWrite(phiPin, code);
  phi = (code / 255.0) * PHI_GAIN * V_PWM_HIGH / V_RAIL;
  delay(PHI_SETTLE_MS);
}

void measureOn() {
  on = true;
  analogWrite(cPin, normToPwm(c));
  analogWrite(sPin, normToPwm(s));
}

void measureOff() {
  on = false;
  analogWrite(cPin, 128);
  analogWrite(sPin, 128);
}

// Collapse toward (xIn, yIn), then release. Same pattern as the trial suite.
void setXY(float xIn, float yIn) {
  xCmd = xIn;
  yCmd = yIn;

  float phiPrev = phi;
  setPhi(0.0); // 0 V keeps theta > 0 for every (x, y); phi can no longer go negative

  c = xIn;
  s = yIn;
  measureOn();
  delay(1000);
  c = 0.0;
  s = 0.0;
  measureOff();

  setPhi(phiPrev);
}

void readState() {
  t = millis() / 1000.0;
  analogRead(thetPin); //Dummy read to address internal capacitor timing
  float thetaRaw = analogRead(thetPin) / 1023.0;
  float xRaw = analogRead(xPin) / 1023.0;
  float yRaw = analogRead(yPin) / 1023.0;
  //theta interface: V_pin = 2V + 0.8*V_theta  ->  V_theta/5V = (thetaRaw - 0.4) / 0.8
  theta = 1.25 * thetaRaw - 0.5;
  x = 2.0 * xRaw - 1.0;
  y = 2.0 * yRaw - 1.0;
}

void tellState() {
  float state[8] = {t, c, s, on, theta, phi, x, y};
  for (int i = 0; i < 8; i++) {
    Serial.print(state[i], 4); // 4 decimals: one phi step is 0.00196 of 5 V
    Serial.print(", ");
  }
  Serial.println();
}

void printHelp() {
  Serial.println(F("circuitDebug — voltages as fractions of 5V"));
  Serial.println(F("t (s), c, s, on/off, theta, phi, x, y"));
  Serial.println(F("Commands: help, status, stream on|off,"));
  Serial.println(F("  c <v>, s <v>, cs <c> <s>, phi <v>  (phi: 0..0.5),"));
  Serial.println(F("  x <v>, y <v>, xy <x> <y>, on, off"));
}

String nextToken(String &src) {
  src.trim();
  int sp = src.indexOf(' ');
  if (sp < 0) {
    String tok = src;
    src = "";
    return tok;
  }
  String tok = src.substring(0, sp);
  src = src.substring(sp + 1);
  src.trim();
  return tok;
}

bool parseFloatToken(String &src, float &out) {
  String tok = nextToken(src);
  if (tok.length() == 0) return false;
  out = tok.toFloat();
  return true;
}

void handleCommand(String cmdLine) {
  cmdLine.trim();
  if (cmdLine.length() == 0) return;

  String rest = cmdLine;
  String cmd = nextToken(rest);
  cmd.toLowerCase();

  if (cmd == "help" || cmd == "?") {
    printHelp();
    return;
  }

  if (cmd == "status") {
    readState();
    tellState();
    return;
  }

  if (cmd == "stream") {
    String arg = nextToken(rest);
    arg.toLowerCase();
    if (arg == "on") {
      streaming = true;
      Serial.println(F("stream on"));
    } else if (arg == "off") {
      streaming = false;
      Serial.println(F("stream off"));
    } else {
      Serial.println(F("usage: stream on|off"));
    }
    return;
  }

  if (cmd == "on") {
    measureOn();
    Serial.println(F("measurement on"));
    return;
  }

  if (cmd == "off") {
    measureOff();
    Serial.println(F("measurement off"));
    return;
  }

  if (cmd == "c") {
    float v;
    if (!parseFloatToken(rest, v)) {
      Serial.println(F("usage: c <val>"));
      return;
    }
    c = v;
    measureOn();
    Serial.print(F("c="));
    Serial.println(c);
    return;
  }

  if (cmd == "s") {
    float v;
    if (!parseFloatToken(rest, v)) {
      Serial.println(F("usage: s <val>"));
      return;
    }
    s = v;
    measureOn();
    Serial.print(F("s="));
    Serial.println(s);
    return;
  }

  if (cmd == "cs") {
    float cv, sv;
    if (!parseFloatToken(rest, cv) || !parseFloatToken(rest, sv)) {
      Serial.println(F("usage: cs <c> <s>"));
      return;
    }
    c = cv;
    s = sv;
    measureOn();
    Serial.print(F("c="));
    Serial.print(c);
    Serial.print(F(" s="));
    Serial.println(s);
    return;
  }

  if (cmd == "phi") {
    float v;
    if (!parseFloatToken(rest, v)) {
      Serial.println(F("usage: phi <val>"));
      return;
    }
    if (v < 0.0 || v > PHI_GAIN) {
      Serial.print(F("phi clipped to 0.."));
      Serial.println(PHI_GAIN, 4);
    }
    setPhi(v);
    Serial.print(F("phi="));
    Serial.print(phi, 4);
    Serial.print(F(" ("));
    Serial.print(phi * V_RAIL, 4);
    Serial.println(F(" V)"));
    return;
  }

  if (cmd == "x") {
    float v;
    if (!parseFloatToken(rest, v)) {
      Serial.println(F("usage: x <val>"));
      return;
    }
    Serial.println(F("collapsing via setXY..."));
    setXY(v, yCmd);
    Serial.print(F("xCmd="));
    Serial.print(xCmd);
    Serial.print(F(" yCmd="));
    Serial.println(yCmd);
    return;
  }

  if (cmd == "y") {
    float v;
    if (!parseFloatToken(rest, v)) {
      Serial.println(F("usage: y <val>"));
      return;
    }
    Serial.println(F("collapsing via setXY..."));
    setXY(xCmd, v);
    Serial.print(F("xCmd="));
    Serial.print(xCmd);
    Serial.print(F(" yCmd="));
    Serial.println(yCmd);
    return;
  }

  if (cmd == "xy") {
    float xv, yv;
    if (!parseFloatToken(rest, xv) || !parseFloatToken(rest, yv)) {
      Serial.println(F("usage: xy <x> <y>"));
      return;
    }
    Serial.println(F("collapsing via setXY..."));
    setXY(xv, yv);
    Serial.print(F("xCmd="));
    Serial.print(xCmd);
    Serial.print(F(" yCmd="));
    Serial.println(yCmd);
    return;
  }

  Serial.print(F("unknown command: "));
  Serial.println(cmd);
  Serial.println(F("type help"));
}

void pollSerial() {
  while (Serial.available() > 0) {
    char ch = (char)Serial.read();
    if (ch == '\r') continue;
    if (ch == '\n') {
      handleCommand(lineBuf);
      lineBuf = "";
    } else {
      lineBuf += ch;
      // Guard against runaway input
      if (lineBuf.length() > 80) lineBuf = "";
    }
  }
}

// c/s PWM on Timer1 (pins 9/10). Arduino's default Timer1 setup is 8-bit
// phase-correct PWM with prescaler 64 (490 Hz). Dropping the prescaler to 1 gives
// 16 MHz / (2 * 255) = 31.37 kHz, which the c/s Sallen-Key filter (f_-3dB ~ 113 Hz)
// attenuates by ~98 dB. analogWrite() keeps working (same 0-255 duty scale), and
// millis()/delay() are untouched (they run on Timer0). Pin 3 (phi) stays on Timer2
// at 490 Hz, which is what the phi filter was designed for. Don't use the Servo
// library: it also takes over Timer1.
#if !defined(__AVR_ATmega328P__) && !defined(__AVR_ATmega168__)
#error "c/s PWM setup assumes an ATmega328P (Uno/Nano), where pins 9/10 are Timer1"
#endif
void setupCsPwm() {
  TCCR1B = (TCCR1B & 0b11111000) | 0b001; // CS12:0 = 001 -> clk/1
}

void setup() {

  // To avoid floating analog pins, it is imperative to set them up as inputs
  pinMode(thetPin, INPUT);
  pinMode(xPin, INPUT);
  pinMode(yPin, INPUT);

  setupCsPwm();
  Serial.begin(baudRate);
  measureOff();
  setPhi(0.0);

  printHelp();
  Serial.println(F("Voltages are recorded as fractions of 5V"));
  Serial.println(F("t (s), c, s, on/off, theta, phi, x, y"));
}

void loop() {
  pollSerial();

  if (streaming) {
    unsigned long now = millis();
    if (now - lastStreamMs >= (unsigned long)dt) {
      lastStreamMs = now;
      readState();
      tellState();
    }
  }
}
