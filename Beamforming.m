%BeamForming parameter
Antenna = 8; 
thetaa = 60;
index = (0 : Antenna-1)';

angle_Range=(0 : 0.01 : 180);
%angle_Range= (1 : 1 : 360);
setss = angle_Range/ 180 * pi;

Hset = exp(1j * pi * index * cos(setss));
x = exp(1j * pi * index * cos(thetaa / 180 * pi));
BeamForming_ULA =  Hset' * x;

polarplot(setss, abs(BeamForming_ULA))
title('Beamforming Antenna:8 Steering Angle:60')

