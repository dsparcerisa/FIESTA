instrreset;

%% Initialize the serial port in MACOS
% arduinoPort= '/dev/tty.usbserial-AR0KKSET';
%% Initialize the serial port in Windows
arduinoPort= 'COM5';

COMtrigger = serial(arduinoPort,'BaudRate',9600);

%% Opens the port
fopen(COMtrigger);
pause(2) % necessary to work

%% Set the parameters of the irradiadiation
display("Working with the Arduino as an external trigger");
aux_npul=input('Set # of pulses: ')%, 's');
%npul = str2num(aux_npul);
aux_dc=input('Set Width Off (us): ')%;, 'i');
aux_width=input('Set width pulses (us): ')%, 'i');


    