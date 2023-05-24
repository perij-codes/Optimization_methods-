% IOE 511/MATH 562, University of Michigan

% Code written by: Janani Peri
% Problem Number: 5
% Problem Name: quartic_1
% Problem Description: A quartic function. Dimension n = 4

% function that computes the gradient of the quartic_1 function
function [g] = quartic_1_grad(x)

% Matrix Q
Q = [5 1 0 0.5;
     1 4 0.5 0;
     0 0.5 3 0;
     0.5 0 0 2];
 
% Set sigma value
sigma = 1e-4;

% compute gradient value
%g = x + sigma*(x'*Q*x)*(Q*x);
g = x + sigma*(Q*x)*x'*(Q*x);

end