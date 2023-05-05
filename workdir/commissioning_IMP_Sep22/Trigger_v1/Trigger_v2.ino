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
float pulses_dutyCicle = 0;  
float f_pulses_widthOn = 0; //us
float f_pulses_widthOff = 0;
int pulses_widthOff = 0;
int pulses_widthOn = 0;
bool pulsado = false; 
bool isUs_off = true;
bool isUs_on = true;

//
char unChar;
String state;
char lastAlpha;
int NumberInt;
float parameters[4];

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
    f_pulses_widthOff = parameters[1];
    f_pulses_widthOn = parameters[2];
    
       //Cálculo delay entre pulsos a partir de dutyCicle
     //pulses_widthOff=int((100-pulses_dutyCicle)*pulses_widthOn/pulses_dutyCicle);
    
    Serial.println(nPulsos);
    Serial.println(f_pulses_widthOff);
    Serial.println(f_pulses_widthOn);
    pulsado = true;
    if(f_pulses_widthOff<10000) {
        isUs_off=true;}
    else{
        isUs_off=false;
        f_pulses_widthOff = f_pulses_widthOff/1e3;
      }
     pulses_widthOff = int(f_pulses_widthOff);

         if(f_pulses_widthOn<10000) {
        isUs_on=true;}
    else{
        isUs_on=false;
        f_pulses_widthOn = f_pulses_widthOn/1e3;
      }
     pulses_widthOn = int(f_pulses_widthOn);

     
    Serial.println(pulses_widthOff);
    Serial.println(pulses_widthOn);
  }
 

  while(pulsado){
    for(int i=0; i<nPulsos; i++) {
      digitalWrite(pinTrigger,HIGH);
      if (isUs_on)
        delayMicroseconds(pulses_widthOn);
      else
        delay(pulses_widthOn);
        
      digitalWrite(pinTrigger,LOW);
      if (isUs_off)
        delayMicroseconds(pulses_widthOff);
      else
        delay(pulses_widthOff);
      }
    pulsado=false;
    
    Serial.println("Pulsed finished");
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
        parameters[i] = state.toFloat();
        state = "";
        i++;
      }
     }    
     else{ state += c; }
    }
 }
