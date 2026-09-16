/*
circuitDebug.ino
Benjamin Krawciw

Interactive serial console for debugging the modelQubit circuit.
Monitors the same voltages as the trial suite, and accepts commands to
set x, y, c, s, and phi from the Serial Monitor.

Values are fractions of 5 V (same convention as modelQubit.ino):
  -1 → -5 V,  0 → 0 V,  +1 → +5 V

Commands (one per line; case-insensitive names):
  help              Print this help
  status            Print one reading now
  stream on|off     Continuously print readings (default: on)
  c <val>           Set c and drive the c pin
  s <val>           Set s and drive the s pin
  cs <c> <s>        Set both c and s
  phi <val>         Set phi and drive the phi pin
  x <val>           Collapse toward (val, last y); see note below
  y <val>           Collapse toward (last x, val)
  xy <x> <y>        Collapse toward (x, y)
  on                Drive current c/s onto the pins
  off               Park c/s at midscale (~0 V)

Note: x and y have no DAC. Setting them measures along that axis for 1 s
(with phi briefly parked slightly negative), then releases measurement.
*/


// Pin numbers (match modelQubit.ino)
const int cPin = 6;
const int sPin = 5;
const int phiPin = 3;
const int thetPin = A1;
const int xPin = A2;
const int yPin = A3;

const int dt = 50; // ms between streamed readings
const float V_RAIL = 5.0;
const float V_CLIP = 0.05;

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

void setPhi(float phiIn) {
  phi = phiIn;
  analogWrite(phiPin, normToPwm(phi));
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
  setPhi(-V_CLIP / V_RAIL);

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
  theta = 2.0 * thetaRaw - 1.0;
  x = 2.0 * xRaw - 1.0;
  y = 2.0 * yRaw - 1.0;
}

void tellState() {
  float state[8] = {t, c, s, on, theta, phi, x, y};
  for (int i = 0; i < 8; i++) {
    Serial.print(state[i]);
    Serial.print(", ");
  }
  Serial.println();
}

void printHelp() {
  Serial.println(F("circuitDebug — voltages as fractions of 5V"));
  Serial.println(F("t (s), c, s, on/off, theta, phi, x, y"));
  Serial.println(F("Commands: help, status, stream on|off,"));
  Serial.println(F("  c <v>, s <v>, cs <c> <s>, phi <v>,"));
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
    setPhi(v);
    Serial.print(F("phi="));
    Serial.println(phi);
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

void setup() {

  // To avoid floating analog pins, it is imperative to set them up as inputs
  pinMode(thetPin, INPUT);
  pinMode(xPin, INPUT);
  pinMode(yPin, INPUT);

  Serial.begin(9600);
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
