% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [g] = exponential_10_grad(x)
    g = 4*(x-1).^3;
    g(1) = 2*exp(x(1))/(exp(x(1))+1)^2 - 0.1*exp(-x(1));
end