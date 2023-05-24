% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [f] = prob2_cons_func(x)
    %evaluates constraints of problem 2
    f = [sum(x.^2)-10, (x(2)*x(3)-5*x(4)*x(5)), (x(1)^3 + x(2)^3 + 1)];
end