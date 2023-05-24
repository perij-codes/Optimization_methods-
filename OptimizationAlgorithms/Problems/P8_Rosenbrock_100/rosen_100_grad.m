%written by svutukury

function [g] = rosen_100_grad(x)
    %returns gradient value of the rosenbrock_100 problem
    %x = 100 x 1
    g = zeros(length(x),1);
    g(1) = -2*(1-x(1)) -400*x(1)*(x(2)-x(1)^2);  
    g(2:99) = 200*(x(2:99)-(x(1:98)).^2)-2*(1-x(2:99)) - 400*x(2:99).*(x(3:100)-x(2:99).^2);
end