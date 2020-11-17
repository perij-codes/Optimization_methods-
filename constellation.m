%% Number of Satellites
% Sruti Vutukury 5/20/20
R_earth= 6378137.0;
h_leo = ((1350-26)*1E3); %altitude of orbit minus max altitude of BB
eps = 10; %deg
alpha = asind((R_earth/(R_earth+h_leo))*cos(eps));
beta = 90 - alpha - eps;

S_coverage = 2*pi*R_earth^2*(1-cos(beta));
S_earth = 4*pi*R_earth^2;
coverage = S_coverage/S_earth;

num_sats = (100/coverage)
with_spares = num_sats*1.1