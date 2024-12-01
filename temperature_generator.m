%function y=temperature_generator(t,A1)
%A0=10+273.15;
%A1=20;
%w=2*pi/(365*86400);
%y=A1.*sin(w.*t-pi/2);


function y=temperature_generator(t,A1)
%A0=10+273.15;
A1=1;                       % Amplitude of the temp. change - 1 degree + or - average, e.g., A1=5 is 5 degrees either side of the average.
w=2*pi/(365*86400);         % dys/ second - this seems to give it an annual temperature change - i.e., unimodal temperature curve over 400 dys.
y=A1.*sin(w.*t-pi/2);       % element wise multiplication of the A1 value. t is in the mcbiology script 