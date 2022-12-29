//Proyecto: Arduino para control de raster
//S. Viñals - CMAM
//8-Jun 2021
//
// PIN10: Señal de trigger para generador de funciones

#define pinTrigger 10

/////
//Definición variables
/////
int nPulsos = 0;  
int pulses_dutyCicle = 0;  
int pulses_widthOn = 0; //us
int pulses_widthOff = 0;
bool pulsado = false; 
//
char unChar;
String state;
char lastAlpha;
int NumberInt;
int parameters[4];

/////
void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  pinMode(pinTrigger, OUTPUT);
 
}

/////
void loop() {

  if(Serial.available()){
    readNumberInt();

    nPulsos = parameters[0];
    pulses_dutyCicle = parameters[1];
    pulses_widthOn = parameters[2];
    
       //Cálculo delay entre pulsos a partir de dutyCicle
     //pulses_widthOff=int((100-pulses_dutyCicle)*pulses_widthOn/pulses_dutyCicle);
    
    pulses_widthOff=pulses_dutyCicle;
    Serial.println(pulses_widthOff);
     pulsado=true;
  }

  while(pulsado){
    for(int i=0; i<nPulsos; i++) {
      digitalWrite(pinTrigger,HIGH);
      delayMicroseconds(pulses_widthOn);
      digitalWrite(pinTrigger,LOW);
      delayMicroseconds(pulses_widthOff);
      }
    Serial.println("Trigger finished");
    pulsado=false;
    }
}

/////
void readNumberInt(){
  delayMicroseconds(300);
  state = "";
  int  i=0;
  while (Serial.available())
  { //Check if there is an available byte to read
    delay(100); //Delay added to make thing stable
    char c = Serial.read(); //Conduct a serial read
    if(c==','){
      if(state.length() >0){
        parameters[i] = state.toInt();
        state = "";
        i++;
      }
     }    
     else{ state += c; }
    }
 }
