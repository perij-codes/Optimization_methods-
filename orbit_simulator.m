%% Orbit Simulator
%% Sruti Vutukury, 5/20/20
%%% initialization
tic
initialize_const()
global const

%%%initialize parent satellite
startr_eci = [-6.782157568250546e+06;8.912514036068704e+05;7.900984236126156e+05]; %in m?
startv_eci = [8.610796068677932e+02;-2.452337305886469e+02;7.551208733440835e+03]; %in km?
quat_body_ecef = [0.452238172076652;0.593352106329713;0.080808724751932;-0.660971908355286];
quat_body_eci = [-0.569316809235482;0.482756659159648;-0.655272469379955;-0.115941233526476];
quat_ecef_eci = [9.172194587818468e-04;3.884327781096487e-05;0.998612074831165;0.052660053181324];
quat_eci_body = [0.569316809235482;-0.482756659159648;0.655272469379955;-0.115941233526476];
quat_eci_ecef = [-9.172194587818468e-04;-3.884327781096487e-05;-0.998612074831165;0.052660053181324];

[a0,e0,E0,I0,omega0,Omega0] = rv2kepler(startr_eci,startv_eci,const.mu); 

%% 1 satellite
const.n = 1; %number of satellites in constellation

intial_positions = zeros(3,length(const.n));
intial_velocities = zeros(3,length(const.n));
intial_positions(:,1) = startr_eci;
intial_velocities(:,1) = startv_eci;

%duration of the mission you want to simulate
num_states = 1.5*60*60; %1.5 hours i.e approximately 1 orbit %length(data) for 1 day
%num_states = 10;
starttime = 17193600; %start time in GPS time in seconds, duration of mission to look at in seconds

%%% propogate orbit (operational phase; in operational phase, does
spherical_earth()
hold on
[t_array, states] = orbit_propagator(intial_positions,intial_velocities,starttime,num_states);

toc
%% constellation
%initialize other satellites in constellation

const.n = 186; %number of satellites in constellation
smaxes = zeros(const.n,1); %semi major axes
eccs = zeros(const.n,1); %eccentricities
anomalies = zeros(const.n,1); %eccentric anomallies
inclins = zeros(const.n,1); %inclinations
omegas = zeros(const.n,1); %argument of periapses
Omegas = zeros(const.n,1); %longitude of the ascending node aka RAAN

smaxes(1) = a0; eccs(1) = e0; anomalies(1) = E0; 
inclins(1) = I0; omegas(1) = omega0; Omegas(1) = Omega0;

%%%deviations between orbits
for s = 1:const.n
    smaxes(s) = a0;
    eccs(s) = e0;
    anomalies(s) = E0;
    inclins(s) = 55*(pi/180); %55 degrees)
    omegas(s) = omega0;
    RAAN =  Omega0 + (60*(pi/180))*s;
    while RAAN > 2*pi
        RAAN  = RAAN - 2*pi;
    end
    Omegas(s) = RAAN;
end

intial_positions = zeros(3,length(const.n));
intial_velocities = zeros(3,length(const.n));
intial_positions(:,1) = startr_eci;
intial_velocities(:,1) = startv_eci;

for i = 2:const.n
    [pos, vel] = kepler2rv(smaxes(i), eccs(i), inclins(i),...
        Omegas(i), omegas(i), anomalies(i), const.mu);
    intial_positions(1:3,i) = pos;
    intial_velocities(1:3,i) = vel;
end

%duration of the mission you want to simulate
%num_states = 1.5*60*60; %1.5 hours i.e approximately 1 orbit %length(data) for 1 day
num_states = 10;
starttime = 17193600; %start time in GPS time in seconds, duration of mission to look at in seconds

%%% propogate orbit (operational phase; in operational phase, does
spherical_earth()
hold on
[t_array, states] = orbit_propagator(intial_positions,intial_velocities,starttime,num_states);

figure(1)
%plotting three satellites only
times =linspace(1,length(states),1);
for i = times
    plot3(states(i,1),states(i,2),states(i,3),'k',...
        states(i,4), states(i,5), states(i,6),'b',...
        states(i,7), states(i,8), states(i,9),'r');
    hold on
end
toc

function [t_array, states] = orbit_propagator(r,v,start_time,duration)
    % INPUT r, v in ECI coordinates in m and m/s
    % start time of mission; gps_time in continuous seconds past 01-Jan-2000 11:59:47 UTC
    % duration is mission duration you want to look at; in seconds
    % OUTPUTS r, v in ECI coordinates in m and m/s at time = current_time + delta_time

    global const

    %package initial state to call ODE45
    state0 = zeros(6,const.n);
    state0(1:3,:)= r;
    state0(4:6,:)= v;
      
    %create tarray
    opts = odeset('RelTol', 1E-12, 'AbsTol', 1E-4, 'OutputFcn','odephas3');
    %opts = odeset('RelTol', 1E-12, 'AbsTol', 1E-4);
    tspan = [0,duration];
    %[t_array, states] = ode113(@state_dot, tspan, state0, opts, start_time);
    [t_array, states] = ode45(@state_dot, tspan, state0, opts, start_time);
    
    orb_elemsf = zeros(length(t_array)); 
    
    for count = 1:length(t_array)
        r  = states(count,1:3)'; %new position 
        v = states(count,4:6)'; %new velocity 
        %quatf = states(count,7:10); %new quaternions
        
        %calculate orbital elements
        [a,e,E,I,omega,Omega] = rv2kepler(r,v,mu);
        orb_elemsf(count,1) = a; orb_elemsf(count,2) = e; orb_elemsf(count,3) = E;
        orb_elemsf(count,4) = I; orb_elemsf(count,5) = omega; orb_elemsf(count,6) = Omega;
        orb_elemsf(count,7) = T; orb_elemsf(count,8) = tp;
    end       
end

    function statef = state_dot(t, state0, start_time)
        %y = [x; y; z; xdot; ydot; zdot]
        fprintf('%f \n',t)

        global const

        % unpack position and velocity in ECI
        r = state0(1:const.n*3);  v = state0(const.n*3+1:const.n*6);
        
        F_envdrag = zeros(3,const.n);
        F_solrad = zeros(3,const.n);
        %acc_tbmoon = zeros(3,24);
        %acc_tbsun = zeros(3,24);
        acc_Js = zeros(3,const.n);
        
        for j = 1:const.n
            
            F_envdrag(1:3,j) = env_atmospheric_drag(r(j*3-2:j*3),v(j*3-2:j*3));
            %F_solrad = env_solradpressure(v_body_unit); %returns acceleration in ECI
            F_solrad(1:3,j) = env_solradpressure();

            %secular perturbations from Moon (third body in circular orbit)
            %returns acceleration in ECI
            now = start_time + t; %current time in seconds
            datetime = utl_time2datetime(now,const.INITGPS_WN); 
            [rp_earth_moon,~] = planetEphemeris(juliandate(datetime),'Moon','Earth','421');
            rp_earth_moon = 1E3*rp_earth_moon'; %positional vector from Moon to Earth; used for 3rd body perturb calcs
            acc_tbmoon(1:3,j) = -const.mu_moon*(((r(j*3-2,j*3)+rp_earth_moon)/norm(r(j*3-2,j*3)+rp_earth_moon)^3) - (rp_earth_moon/norm(rp_earth_moon)^3)); 

            %secular perturbations from Sun (third body in circular orbit)
            %returns acceleration in ECI
            [rp_earth,~] = planetEphemeris(juliandate(datetime),'Sun','Earth','421');
            rp_earth = 1E3*rp_earth; %positional vector from Sun to Earth; used for 3rd body perturb and solar radiation pressure calcs
            acc_tbsun(1:3,j) = -const.mu_sun*(((r(j*3-2,j*3)+rp_earth')/norm(r(j*3-2,j*3)+rp_earth')^3) - (rp_earth'/norm(rp_earth')^3)); 

             % gravity + perturbations due to J-coefficients; returns acceleration in ECEF
            [quat_ecef_eci,~]=env_earth_attitude(now);
            quat_eci_ecef= utl_quat_conj(quat_ecef_eci);
            pos_ecef=utl_rotateframe(quat_ecef_eci,r(j*3-2:j*3));

            [g_ecef(1),g_ecef(2),g_ecef(3)]  = gravitysphericalharmonic(pos_ecef', 'EGM96',2);
            %g_ecef=env_gravity(now,pos_ecef);

            %convert back to to ECI
            %acc_Js(1:3,j) =utl_rotateframe(quat_eci_ecef,g_ecef);
            acc_Js(1:3,j) = [8.29068446288587;-1.08941247031678;-0.968588264273039];
        end

        statef = zeros(length(state0),1);
        statef(1:const.n*3)=state0(const.n*3+1:const.n*6);

        %acceration in m/s^2, ECI coords
        %statef(73:144) = ((F_envdrag+F_solrad)./(const.MASS)) + acc_tbmoon + acc_tbsun + acc_Js;
        %statef(const.n*3+1:const.n*6) = ((F_envdrag+F_solrad)./(const.MASS)) + acc_Js;
        statef(const.n*3+1:const.n*6) = acc_Js;
    end 


function F_envdrag = env_atmospheric_drag(r,v)
    %input: mission time in sec
    % r position in ECEF frame in m
    % v velocity in ECEF frame in m/s

    % outputs: Fenv_drag in direction of velocity ECEF frame
    % assumes a simple, fully static exponentially decaying model (Vallado)
    %%convert ECEF coordinates to latitude, longitude, altitude (LLA) geodetic coordinates
    lla = ecef2lla(r'); 

    h = abs(lla(3)); %altitude of Position (m); geocentric altitude = geodetic altitude
    
        % in meters
    hmin = 1E3*[0;25;30;35;40;45;50;55;60;65;...
        70;75;80;85;90;95;100;110;120;130;...
        140;150;160;180;200;250;300;350;400;...
        450;500;600;700;800;900;1000];
    % in meters
    hmax = 1E3*[25;30;35;40;45;50;55;60;65;70;75;...
        80;85;90;95;100;110;120;130;140;150;160;...
        180;200;250;300;350;400;450;500;600;...
        700;800;900;1000;inf];
    % in meters
    h0 = 1E3*[0;25;30;35;40;45;50;55;60;...
        65;70;75;80;85;90;95;100;110;...
        120;130;140;150;160;180;200;250;300;...
        350;400;450;500;600;700;800;900;1000];
    % kg/m^3
    rho0 = [1.225; 3.899E-2; 1.774E-2;...
        8.279E-3; 3.972E-3; 1.995E-3;...
        1.057E-3; 5.821E-4; 3.206E-4; ...
        1.718E-4; 8.770E-5; 4.178E-5;...
        1.905E-5; 8.337E-6; 3.396E-6; ...
        1.343E-6; 5.297E-7; 9.661E-8; ...
        2.438E-8; 8.484E-9; 3.845E-9; ...
        2.070E-9; 1.224E-9; 5.464E-10;...
        2.789E-10; 7.248E-11; 2.418E-11;...
        9.158E-12; 3.725E-12; 1.585E-12;...
        6.967E-13; 1.454E-13; 3.614E-14;...
        1.170E-14; 5.245E-15; 3.019E-15];
    % in meters
    H = 1E3*[8.44; 6.49; 6.75; 7.07; 7.47; 7.83;...
        7.95; 7.73; 7.29; 6.81; 6.33; 6.00;...
        5.70; 5.41; 5.38; 5.74; 6.15; 8.06;...
        11.6; 16.1; 20.6; 24.6; 26.3; 33.2;...
        38.5; 46.9; 52.5; 56.4; 59.4; 62.2;...
        65.8; 79.0; 109.0; 164.0; 225.0; 268.0];   
    idx = 0;
    for i = 1:length(hmin)
        if abs(h) >= hmin(i) && abs(h) < hmax(i)
            idx = i;
            break
        end
    end    
    rho = rho0(idx)*exp(-(h-h0(idx))/H(idx));
    
    %[rho,~] = get_rho(h);

    %v_body = utl_rotateframe(quat_body_ecef,v);
    %v_body_unit = v_body/norm(v_body);

    Cd = 1.15; %conservative drag coeff
    %A = utl_area(v_body_unit);
    A = 10;

    F_envdrag = -0.5*rho*Cd*A*(v*v')*(v./norm(v));
    F_envdrag = reshape(F_envdrag,[3,1]);

end

function [rho,H] = get_rho(h)
    
    % in meters
    hmin = 1E3*[0;25;30;35;40;45;50;55;60;65;...
        70;75;80;85;90;95;100;110;120;130;...
        140;150;160;180;200;250;300;350;400;...
        450;500;600;700;800;900;1000];
    % in meters
    hmax = 1E3*[25;30;35;40;45;50;55;60;65;70;75;...
        80;85;90;95;100;110;120;130;140;150;160;...
        180;200;250;300;350;400;450;500;600;...
        700;800;900;1000;inf];
    % in meters
    h0 = 1E3*[0;25;30;35;40;45;50;55;60;...
        65;70;75;80;85;90;95;100;110;...
        120;130;140;150;160;180;200;250;300;...
        350;400;450;500;600;700;800;900;1000];
    % kg/m^3
    rho0 = [1.225; 3.899E-2; 1.774E-2;...
        8.279E-3; 3.972E-3; 1.995E-3;...
        1.057E-3; 5.821E-4; 3.206E-4; ...
        1.718E-4; 8.770E-5; 4.178E-5;...
        1.905E-5; 8.337E-6; 3.396E-6; ...
        1.343E-6; 5.297E-7; 9.661E-8; ...
        2.438E-8; 8.484E-9; 3.845E-9; ...
        2.070E-9; 1.224E-9; 5.464E-10;...
        2.789E-10; 7.248E-11; 2.418E-11;...
        9.158E-12; 3.725E-12; 1.585E-12;...
        6.967E-13; 1.454E-13; 3.614E-14;...
        1.170E-14; 5.245E-15; 3.019E-15];
    % in meters
    H = 1E3*[8.44; 6.49; 6.75; 7.07; 7.47; 7.83;...
        7.95; 7.73; 7.29; 6.81; 6.33; 6.00;...
        5.70; 5.41; 5.38; 5.74; 6.15; 8.06;...
        11.6; 16.1; 20.6; 24.6; 26.3; 33.2;...
        38.5; 46.9; 52.5; 56.4; 59.4; 62.2;...
        65.8; 79.0; 109.0; 164.0; 225.0; 268.0];   
    idx = 0;
    for i = 1:length(hmin)
        if abs(h) >= hmin(i) && abs(h) < hmax(i)
            idx = i;
            break
        end
    end    
    rho = rho0(idx)*exp(-(h-h0(idx))/H(idx));
end

function F_solrad = env_solradpressure()
    %  Computes the acceleration due to solar radiation pressure assuming the spacecraft surface normal to the Sun direction
    %  assumes a cylindrical shadow model
    % Output: solar radiation pressure force in ECI

    Fs = 1367; %solar constant
    c = 3E8; %speed of light
    %A = utl_area(v_body_unit);
    A = 10;
    q = 1; %conservative
    i = 0; %conservative

    F_solrad = (Fs/c)*A*(1+q)*cos(i);
end

function initialize_const()
    global const

    %Time
    const.INITGPS_WN= 2045;% positive int
    % initial gps week number, epoch for time.
    const.INIT_DYEAR= decyear(utl_time2datetime(0.0,const.INITGPS_WN));
    % Earth's gravitational constant (m^3/s^2)
    const.mu = 3986004.415e8;%3.986e14;% positive scalar
    % Moon's gravitational constant (m^3/s^2)
    const.mu_moon = 4.9048695E12; % positive scalar
    %largest planar area of satellite in m^2
    const.satArea = 0.1*sqrt(2)*0.3;
    % Sun's gravitational constant (m^3/s^2)
    const.mu_sun = 1.32712440018E20; % positive scalar
    %Equatorial Radius of Earth (m)*/
    const.R_EARTH= 6378137.0;
    % Simulation timestep            (ns)
    const.dt = int64(1e8);% positive int64
    % Earth's eccentricity.
    const.e_earth = 0.0167086;
    % Earth orbital period (s)
    const.period_earth = 365.256363004*24*60*60;
    % Astronomical unit [m]; DE430
    const.AU = 149597870700.000000;
    % Speed of light  [m/s]; DE430
    const.c_light   = 299792458.000000000;
    % Solar radiation pressure at 1 AU
    const.P_Sol = 1367/const.c_light; % [N/m^2] (1367 W/m^2); IERS 96
    %Solar radiation pressure coefficient
    const.Cr = 1; %dimensionless
    perihelion_date = datetime(2019,1,3,5,20,0,'TimeZone','UTCLeapSeconds');
    [rp_earth_moon,vp_earth_moon] = planetEphemeris(juliandate(perihelion_date),'Moon','Earth');
    const.rp_earth_moon = 1E3*rp_earth_moon'; %positional vector from Moon to Earth; used for 3rd body perturb calcs
    [rp_earth,vp_earth] = planetEphemeris(juliandate(perihelion_date),'Sun','Earth');
    const.rp_earth = 1E3*rp_earth; %positional vector from Sun to Earth; used for 3rd body perturb and solar radiation pressure calcs
    rp_earth = rp_earth';
    vp_earth = vp_earth';
    h_earth = cross(rp_earth,vp_earth);
    % Quat between earth's perifocal and eci frame.
    T0=utl_time2datetime(0,const.INITGPS_WN);%pan epoch
    T5=T0+years(5);%5 years after pan epoch
    dcm_ECEF0_ECI=dcmeci2ecef('IAU-2000/2006',[year(T0),month(T0),day(T0),hour(T0),minute(T0),second(T0)]);
    dcm_ECEF5_ECI=dcmeci2ecef('IAU-2000/2006',[year(T5),month(T5),day(T5),hour(T5),minute(T5),second(T5)]);
    polarprecessionaxis= -cross((dcm_ECEF0_ECI*dcm_ECEF5_ECI'*[0;0;1;]),[0;0;1;]);
    const.MASS = 1000;
    const.PRECESSION_RATE= polarprecessionaxis/seconds(T5-T0);% 3 vector
    % earth's axis precession rate (rad/s)
    const.quat_ecef0_eci= utl_quaternion2array(quaternion(dcm_ECEF0_ECI,'rotmat','frame'));
    %ecef0 is ecef frame at time 0 inertialy stuck.
    const.earth_rate_ecef=[sin(0.2/3600*pi/180);sin(-0.30/3600*pi/180);1;]*7.2921158553E-5;% 3 vector
    % earth's inertial rotation rate in ecef frame (rad/s)
end


function A = utl_area(v_body_unit)
    %%% utl_area calculates the total area of satellite exposed that is normal
    %%% to inputted velocity vector
    %%% input: unit velocty column vector in body frame, output: area in m^2
    A = dot(abs(v_body_unit),[0.03,0.03,0.01]);
end

function u = utl_rotateframe(q,v)
    %UTL_ROTATEFRAME quaternion frame rotation
    %  U = ROTATEFRAME(Q,V) rotates the frame of reference for the three vector V
    %  using quaternion Q stored as an array with 4th component real. 

    u= v+cross(2*q(1:3),cross(q(1:3),v)-q(4)*v);
end

function [r_ECI,v_ECI] = env_ECEFtoECI(time,r_ECEF,v_ECEF)
%Returns a position and velocity in ECI from ECEF, units s, m, m/s
%  The shape of the vectors should be (3,1)
%   time(double): seconds since const.INITGPS_WN
%   Modified 12/12/2019 by Sruti
[quat_ecef_eci,rate_ecef_eci_ecef]=env_earth_attitude(time);
quat_eci_ecef= utl_quat_conj(quat_ecef_eci);
r_ECI=utl_rotateframe(quat_eci_ecef,r_ECEF);
v_ECEF= v_ECEF+ cross(rate_ecef_eci_ecef,r_ECEF);
v_ECI=utl_rotateframe(quat_eci_ecef,v_ECEF);
end


function qout = utl_quat_conj(q)
%utl_quat_conj Congugate q. This reverses the direction of rotation.
% Quaternions have the forth component as scalar. 
qout= [-q(1:3);q(4)];
end

function [quat_ecef_eci,rate_ecef]= env_earth_attitude(time)
%earth_attitude Return the quaternion and anguar rate to rotate vectors from eci to ecef
%   time(double): time in seconds since const.INITGPS_WN
%#codegen
global const
rate_ecef=const.earth_rate_ecef;
quat_ecef0p_ecef0= [const.PRECESSION_RATE*time; 1.0;];
quat_ecef0p_ecef0= quat_ecef0p_ecef0/norm(quat_ecef0p_ecef0);
earth_axis= const.earth_rate_ecef/norm(const.earth_rate_ecef);
theta= norm(const.earth_rate_ecef)*time;% earth rotation angle
quat_ecef_ecef0p= [earth_axis*sin(theta/2);cos(theta/2);];
quat_ecef0_eci=const.quat_ecef0_eci;
quat_ecef_eci=utl_quat_cross_mult(quat_ecef_ecef0p,utl_quat_cross_mult(quat_ecef0p_ecef0,quat_ecef0_eci));
end

function out_datetime = utl_time2datetime(time,init_GPS_week_number)
%time2datetime converts time to a datetime with timezone 'UTCLeapSeconds'
%   time(double): time since init_GPS_week_number in seconds
%   init_GPS_week_number(int): initial GPS week number.
pan_epoch=seconds(init_GPS_week_number*7*24*3600)+datetime(1980,1,6,'TimeZone','UTCLeapSeconds');
out_datetime = seconds(time)+pan_epoch;
end

function [a,e,true_anom,I,omega,Omega] = rv2kepler(r,v,mu) 

%Convert orbital position and velocity to orbital elements
%Inputs:
%   r - 3x1 array: orbital position vector in ECI
%   v - 3x1 array: orbital velocity vector in ECI
%   mu - scalar: gravitational parameter
%
%Outputs:
%   a - scalar: semi-major axis
%   e - scalar: eccentricity
%   E - scalar: eccentric anomaly (in radians)
%   nu - scalar: true anomaly
%   I - scalar: inclination (in radians)
%   omega - scalar: argument of periapsis (in radians)
%   Omega - scalar: longitude of the ascending node (in radians)
%   T - scalar: orbital period
%   tp - scalar: time of periapse passage

%ensure that the inputs are column vectors
r = r(:);
v = v(:);
r_mag  = norm(r);
v_mag  = norm(v);

%calculate a using Energy
Energy = (v_mag^2/2)-(mu/r_mag);
a = -mu/(2*Energy)   ; %semi-major axis

h = cross(r,v);
h_mag = norm(h);

n = cross([0;0;1],h);
n_mag = norm(n);

e_vec = cross((v./mu),h)-(r./r_mag);
e = norm(e_vec); %eccentricity

cosE = (1/e)*(1-(r_mag/a)) ; %cosine of the eccentric anomaly (from radius eq)
true_anom = acos((a*cosE-a*e)/r_mag);
%sinE = (sin(true_anom)*sqrt(1-e^2))/(1+e*cos(true_anom)); %sine of eccentric anomaly (from velocity eq.)
%E = mod(atan2(sinE,cosE),2*pi); %ecentric anomaly

I = acos(dot((h./h_mag),[0;0;1]))  ; %inclination
I = mod(I,2*pi); %inclination must be in proper range

omega = acos(dot(n,e_vec)/(n_mag*e))  ; %argument of periapsis
omega = 2*pi - mod(omega,2*pi); %omega 

Omega = acos(dot([1;0;0],(n./n_mag)))    ; %longitude of ascending node
Omega = 2*pi - mod(Omega,2*pi);

T = (2*pi/(sqrt(mu)))*a^(3/2); %orbital period

%M = (E - e*sinE);
mean_motion = (2*pi)/T;
%tp = -M/mean_motion; %time of periapsis passage

% %we'd like the time of periapsis passage to strictly positive
% if tp < 0
%     tp = T+tp;
% end

end


function [r, v] = kepler2rv(a, e, i, O, o, nu, mu)
    %convert Keplerian Orbital Elements to state vectors in ECI                                                 %
    % Inputs:                                                                 %
    %--------                                                                  
    % a semi-major axis
    %e Eccentricity Magnitude (unitless)
    %i nclination (radians)
    %O Right Ascention of the ascending node (radians)
    %
    %o argument of periapse (radians)
    %nu true anomaly (radians)
    %truLon True Longitude (radians)
    %
    %argLat argument of latitude (radians)
    %lonPer Longitude of Periapse (radians)
    %mu Gravitational Constant

    % Outputs:
    %r_ECI  position vecotr ECI
    %v_ECI  velocity vector ECI   
    %------------------------------------------------------------------       %
    % Programed by Darin Koblick  03-04-2012                                  %
    %------------------------------------------------------------------       %
    p = a*(1-e^2); %Semilatus Rectum km
    %Find rPQW and vPQW
    rPQW = cat(2,p.*cos(nu)./(1 +e.*cos(nu)),p.*sin(nu)./(1+e.*cos(nu)),zeros(size(nu)));
    vPQW = cat(2,-sqrt(mu./p).*sin(nu),sqrt(mu./p).*(e+cos(nu)),zeros(size(nu)));
    %Create Transformation Matrix
    PQW2IJK = NaN(3,3,size(p,1));
    cO = cos(O); sO = sin(O); co = cos(o); so = sin(o); ci = cos(i); si = sin(i);
    PQW2IJK(1,1,:) = cO.*co-sO.*so.*ci; PQW2IJK(1,2,:) = -cO.*so-sO.*co.*ci; PQW2IJK(1,3,:) = sO.*si;
    PQW2IJK(2,1,:) = sO.*co+cO.*so.*ci; PQW2IJK(2,2,:) = -sO.*so+cO.*co.*ci; PQW2IJK(2,3,:) = -cO.*si;
    PQW2IJK(3,1,:) = so.*si;            PQW2IJK(3,2,:) = co.*si;             PQW2IJK(3,3,:) = ci;
    %Transform rPQW and vPQW to rECI and vECI
    r = multiDimMatrixMultiply(PQW2IJK,rPQW)';  
    v = multiDimMatrixMultiply(PQW2IJK,vPQW)';
end

function c = multiDimMatrixMultiply(a, b)
    c = NaN(size(b));
    c(:,1) = sum(bsxfun(@times,a(1,:,:),permute(b,[3 2 1])),2);
    c(:,2) = sum(bsxfun(@times,a(2,:,:),permute(b,[3 2 1])),2);
    c(:,3) = sum(bsxfun(@times,a(3,:,:),permute(b,[3 2 1])),2);
end

function out = utl_quaternion2array(quat)
    %UTL_QUATERNION2ARRAY converts a quaternion to a column vector length 4
    %   the 4th component is real
    assert(isa(quat,'quaternion'),'input must a quaternion')
    q=compact(quat);
    out= [q(2);q(3);q(4);q(1)];
end