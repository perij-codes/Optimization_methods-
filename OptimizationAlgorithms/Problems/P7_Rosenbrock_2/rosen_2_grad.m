%written by svutukury

function [g] = rosen_2_grad(x)
    g = [(-2 + 2*x(1) - 400*x(2)*x(1) + 400*x(1)^3); (200*x(2) - 200*x(1)^2)];
end