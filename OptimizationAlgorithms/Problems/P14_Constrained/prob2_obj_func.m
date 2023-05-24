% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [f] = prob2_obj_func(x)
    %evaluates function value of problem 2
    f = exp(prod(x))-0.5*(x(1)^3 + x(2)^3 + 1)^2;
end