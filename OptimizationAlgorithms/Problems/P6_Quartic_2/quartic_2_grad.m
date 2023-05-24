% IOE 511/MATH 562, University of Michigan

% Code written by: Janani Peri

% Problem Number: 6
% Problem Name: quartic_2
% Problem Description: A quartic function. Dimension n = 4

% function that computes the gradient value of the quartic_2 function
function [g] = quartic_2_grad(x)

% Matrix Q
Q = [5 1 0 0.5;
     1 4 0.5 0;
     0 0.5 3 0;
     0.5 0 0 2];
 
% Set sigma value
sigma = 1e4;

% compute gradient value
g = x + sigma*(Q*x)*x'*(Q*x);

end