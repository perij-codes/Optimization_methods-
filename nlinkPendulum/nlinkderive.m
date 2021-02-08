%%%derive lagrange rhs 
%%%sruti vutukury

syms ihat jhat khat t real
ihat = [1 0 0]; jhat=[0 1 0]; khat=cross(ihat,jhat); %standard unit vectors

%%% INITIALIZE PARAMETERS
p.N = 2; %number of links
p.g = 9.81;
p.d = sym('d',[1 p.N],'real');
p.L = sym('L',[1 p.N],'real');
p.m = sym('m',[1 p.N],'real');
p.Ig = sym('Ig',[1 p.N],'real');

%initial symbolic variables
ths = sym('th',[1 p.N],'real');
thdots = sym('thdot',[1 p.N],'real');
thddots = sym('thddot',[1 p.N],'real');
er = sym([],'real');
et = sym([],'real');
rOG = sym([],'real');
rOE = sym([],'real');
rEG = sym([],'real');
vOG = sym([],'real');
vEG = sym([],'real');
vOE = sym([],'real');
EKs = sym([],'real');
EPs = sym([],'real');

for i = 1:p.N
%     ths(i) = sym(sprintf('th%d', i));
%     thdots(i) = sym(sprintf('th%ddot', i));
%     thddots(i) = sym(sprintf('th%dddot', i));

    %polar coordinates for each link
    er(i,:) = [cos(ths(i)) sin(ths(i)) 0];
    et(i,:) = cross(khat,er(i,:));

end
% ths = ths';
% thdots = thdots';
% thddots = thddots';
% rOG = rOG';
% rOE = rOE';
% rEG = rEG';
% vOG = vOG';
% vEG = vEG';
% vOE = vOE';

for i = 1:p.N
    %position vectors for each link
    if i == 1
        rOG(i,:) = p.d(i)*er(i,:); %MAIN
        rOE(i,:) = p.L(i)*er(i,:); %i.e. from E1 to O
    else
        %end of prev link position + rel position to link CG
        rEG(i,:) = p.d(i)*er(i,:); %relative position
        rOG(i,:) = rOE(i-1,:) + rEG(i,:); %MAIN
        %end of prev link position + rel position to end of current link
        rOE(i,:) = rOE(i-1,:) + p.L(i)*er(i,:);
    end

    %velocity vectors for each link
    if i == 1
        vOG(i,:) = cross(thdots(i)*khat, rOG(i,:)); %MAIN
        vOE(i,:) = cross(thdots(i)*khat, rOE(i,:));
%         vOG(i,:) = thdots(i)*p.d(i)*et(1,:); %MAIN
%         vOE(i,:) = thdots(i)*p.L(i)*et(1,:);
    else
        vEG(i,:) = cross(thdots(i)*khat, rEG(i,:));
        vOE(i,:) = cross(thdots(i)*khat, rOE(i,:));
        vOG(i,:) = vEG(i,:) + vOE(i-1,:); %MAIN

%         vEG(i,:) = thdots(i)*p.d(i)*et(1,:);
%         vOE(i,:) = cross(thdots(i)*khat, rOE(i));
%         vOG(i,:) = vEG(i,:) + vOE(i,:); %MAIN
    end

    %potential energy per link
    EPs(i) = p.m(i)*p.g*rOG(i,:)*jhat'; %or *ihat.'???
    %kinetic energy per link: trans + rotational
    EKs(i) = p.m(i)*dot(vOG(i,:),vOG(i,:)) + p.Ig(i)*thdots(i)^2;
end

%total energies
EP = abs(sum(EPs));
EK = abs(0.5*sum(EKs));

%lagrangian
L = EK-EP; %THE PROBLEM IS HERE. SHOULD BE TH1 NOT THDOT1

dLdthdots = jacobian(L, thdots);  % fs = [jacobian(L, th1dot),jacobian(L, th2dot)t,...]

%Lagrange equations
%The first three lines are calculating d/dt   using the chain rule
eqns =   jacobian(dLdthdots,ths   )*thdots'  ...
       + jacobian(dLdthdots,thdots)*thddots' ...
       + jacobian(dLdthdots,t)               ...
       - jacobian(L,ths)';

%get these eqns into a RHS file;
[M,b] = equationsToMatrix(eqns,thddots); 
M_L = simplify(M);
b_L = simplify(b);

matlabFunction(M_L,'file','tryM');
matlabFunction(b_L,'file','tryR');

disp('done')