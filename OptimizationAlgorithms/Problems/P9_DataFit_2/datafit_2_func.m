% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [f] = datafit_2_func(x)
    y = [1.5, 2.25, 2.625]';
    w = x(1);
    z = x(2);
    f = 0;
    for i = 1:3
        s = (y(i) - w*(1-z^i))^2;
        f = f+s;
    end
end