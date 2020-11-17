%%%%% prop system sizing

%ion propulsion
Force = 0.5; %N
delv1 = 1.2686993E11; %m/s for entry
delv2 = 2.18408832493985E9; %m/s for flyby phase
m = 11100; %kg
dist1 = 1.496E11; %1 AU
dist2 = 8.228E13; %550 AU

del_time1 = (m*delv1)/Force;
del_time2 = (m*delv2)/Force;

yearsec = 60*60*24*365;
fprintf('time from earth to entry in years: \n')
disp(del_time1/yearsec);
fprintf('time to get to 550 AU in years: \n')
disp(del_time2/yearsec);

%nuclear pulse propulsion
Force = 1E12; %N
delv1 = 1.269E11; %m/s for entry
delv2 = 2.184E9; %m/s for flyby phase
dist1 = 1.496E11; %1 AU
dist2 = 8.228E13; %550 AU
time_limit = 3.154E9;
yearsec = 60*60*24*365;
m = 11100; %kg

syms del_time1 del_time2
eqns = [ del_time1 - (m*delv1/Force) == 0, del_time2 - (m*delv2/Force) == 0;];
A = solve(eqns, [del_time1 del_time2]);
fprintf('burn time from earth to entry in sec: \n')
disp(double(A.del_time1)) %in sec
fprintf('burn time to get to 550 AU in sec: \n')
disp(double(A.del_time2)) %in sec


function days = earthSunTime()
    %time to travel from parking orbit to solar fly-by entry including burn 1
    m = 11100; %mass of spacecraft in kg
    vi = 1.726E4; %vehicle ejection velocity
    vf = 2.953218208E4; %vehicle velocity after burn 1
    delv = 1.22719680E4; %delv of burn 1
    burn_time = 1.4086E3; %burn time of burn 1
    thrust = (m*delv)/burn_time; %thrust from impulse
    dtot = 1.496e+11; %distance from earth to sun at max perihelion; from database
    d1 = (vf^2 - vi^2)/(2*(thrust/m)); %distance covered during burn 1
    d2 = dtot - d1;
    drift_time = d2/vf; %drift time after burn 1

    total_time = burn_time + drift_time;
    days = total_time/(60*60*24); %total earth-sun travel time in days

end