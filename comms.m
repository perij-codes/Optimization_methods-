%% Communications Architecture
% Sruti Vutukury 5/20/20
%% Link Budget
%%% solve for eta, modulation, rolloff factor, Gt (transmitter gain)
Eb_N0 = 5.4; %dB;
Eb_N0min = 2.4; %dB
syms Gt
Rb = 4.3E6; %bps
Pt = 33; %transmitter power in W
Gr = 10^(31.4/10); %ground receiver gain; converting dB to Watts
La = 10^(-1.5/10); %atmospheric losses; converting dB to Watts
Ll = 10^(-166.8/10); %line losses; converting dB to Watts
c = 3E8; %speed of light
lambda = (c/2290E6); %wavelength in m; S-band frequency
d = 1350E3; %m; mean altitude of orbit
Kb = 1.38E-23; %Boltzmann constant J/K )
Tsys = 273+50; %K operating temperature of the antenna

eqn = ((Eb_N0) - (Pt*Gt*Gr*La*Ll*lambda^2)/(((4*pi*d)^2)*Kb*Tsys*Rb) == 0);
Gt = double(solve(eqn, Gt)) %gain in Watts

%% BlackBird Song Downlink
clc; clear all
f = 40000;
bdepth = 16;
channels = 1;
length = (2*60)+19; %lenght of song in seconds
Rb = f*bdepth*channels
comp = 4.4;
size_uncomp = (Rb*length)
size_comp = (Rb*length)/comp

%% Number of Ground Stations
%% method 1
clc; clear all; close all;
beamwidth = 5.2*(pi/180); %in rad
h = 1350E3; %LEO orbit altitude in m
re = 6378000; %earth radius in m

r_ground = h*tan(beamwidth/2);
cov_ind  = pi*r_ground^2; %coverage of each ground station
cov_needed = 4*pi*(re+h)^2; %for full coverage at 1350 km
num_stations = cov_needed/cov_ind

%% method 2
clc; clear all; close all;
beamwidth = 5.2*(pi/180); %in rad
h = 1350E3; %LEO orbit altitude in m
re = 6378000; %earth radius in m

area_each = beamwidth^2*h^2
cov_needed = 4*pi*(re+h)^2; %for full coverage at 1350 km
num_stations = cov_needed/area_each

%% method 3
clc; clear all; close all;
beamwidth = 5.2*(pi/180); %in rad
re = 6378000; %earth radius in m
num_stations = 10000; 

syms h r_ground cov_ind cov_needed real positive
eqns = num_stations - (4*pi*(re+h)^2)/(pi*(h*tan(beamwidth/2))^2) == 0
h = double(solve(eqns,h))/10^3 %altitude in km

