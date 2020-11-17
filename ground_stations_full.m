%%%%plotting ground stations and field of view
%%%%Sruti Vutukury, 5/20/2020

clc; clear all; close all
%%%earth mesh
r = 6378000;
[X,Y,Z] = sphere;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;
surf(X2+5,Y2-5,Z2,'FaceAlpha',0)
axis equal
hold on

%%%ground stations
N = 61;
[X,Y,Z,N_new] = mySphere(N);
X = X.*6378010;
Y = Y.*6378010;
Z = Z.*6378010;
scatter3(X,Y,Z,'filled','k')

function [X,Y,Z,N_new] = mySphere(N)
    % Generate Node xyz positions
    % Used 2004 paper by Markus Deserno, Max-Planck-Institut:
    % "How to generate equidistributed points on the surface of a sphere"
    % Enforces constant intervales d_theta ~ d_phi
    % Assumes unit radius
    % Does not replace MATLAB "sphere" function
    % Create Sphere 3D Geometry Centered at (x,y,z) = (0,0,0)
    %
    % N: target number of nodes
    % N_new: final number of nodes
    % X,Y,Z: column vectors of length N_new containing node coordinates
    r_unit = 1;
    Area = 4*pi*r_unit^2/N;
    Distance = sqrt(Area);
    M_theta = round(pi/Distance);
    d_theta = pi/M_theta;
    d_phi = Area/d_theta;
    N_new = 0;
    for m = 0:M_theta-1

        Theta = pi*(m+0.5)/M_theta;
        M_phi = round(2*pi*sin(Theta)/d_phi); % not exact

        for n = 0:M_phi-1        
            Phi = 2*pi*n/M_phi;    

            N_new = N_new + 1;

            X(N_new) = sin(Theta)*cos(Phi);
            Y(N_new) = sin(Theta)*sin(Phi);
            Z(N_new) = cos(Theta);

        end
    end
end