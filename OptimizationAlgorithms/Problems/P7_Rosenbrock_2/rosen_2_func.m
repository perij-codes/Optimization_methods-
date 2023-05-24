%written by svutukury

function [f] = rosen_2_func(x)
    f = (1 - x(1))^2 + 100*(x(2) - x(1)^2)^2;
end