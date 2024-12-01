function y=nh3_mineralization(t,A1)
%A0=10+273.15;
%A1=20;
w=2*pi/(365*86400);
y=A1.*sin(w.*t-pi/2);

id=(y<0);
y=y.*(1-id);
end

%function y=nh3_mineralization(t,A1)
%A0=10+273.15;
%A1=1;                     % Amplitude of the temp. change - 1 degree + or - average
%w=1*pi/(365*86400);         %dys/ second
%y=A1.*sin(w.*t-pi/2);