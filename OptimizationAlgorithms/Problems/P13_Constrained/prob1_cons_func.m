% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [f] = prob1_cons_func(x)
    %evaluates function value of the constraint of problem 1
    f = (x(1)^2 + x(2)^2 -2)';
end