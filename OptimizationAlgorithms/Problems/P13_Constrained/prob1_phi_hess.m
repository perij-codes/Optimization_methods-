% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [h] = prob1_phi_hess(x,nu,problem)
    %evaluates phi hessian value of problem 1
    x1 = x(1);
    x2 = x(2);
    h = nu*[(6*x1^2 + 2*x2^2 -4) (4*x1*x2);
        (4*x2*x1) (2*x1^2 +6*x2^2 -4)];
end