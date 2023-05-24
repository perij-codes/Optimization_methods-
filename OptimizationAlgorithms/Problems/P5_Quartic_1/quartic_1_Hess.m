% IOE 511/MATH 562, University of Michigan
% Code written by: Albert S. Berahas
% Modified by: Janani Peri
% Problem Number: 4
% Problem Name: quad_1000_1000
% Problem Description: A randomly generated convex quadratic function; the 
%                      random seed is set so that the results are 
%                      reproducable. Dimension n = 1000; Condition number
%                      kappa = 1000

% function that computes the Hessian of the quartic_1 function
function [H] = quartic_1_Hess(x)

% Matrix Q
Q = [5 1 0 0.5;
     1 4 0.5 0;
     0 0.5 3 0;
     0.5 0 0 2];
 
% Set sigma value
sigma = 1e-4;

% compute hessian value
H = eye(length(Q)) + sigma*2*(Q*Q)*(x*x');
%H = eye(length(Q)) + sigma*(((x'*Q)*x)*Q + 2*Q*x*(Q*x)');
end