%written by svutukury

function [f] = rosen_100_func(x)
    %returns function value of the rosenbrock_100 problem
    %x = 100 x 1
    f = sum((1 - x(1:99)).^2 + 100*(x(2:100)-x(1:99).^2).^2);
end