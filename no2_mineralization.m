function y=no2_mineralization(t,A1)
%A0=10+273.15;
%A1=2;
w=2*pi/(365*86400);
y=A1.*sin(w.*t-pi/2);

id=(y<0);
y=y.*(1-id);
end