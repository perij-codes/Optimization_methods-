% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [g] = prob1_cons_grad(x)
    %evaluates gradient of constraint of problem 1
    g = [2*x(1);2*x(2)];
end