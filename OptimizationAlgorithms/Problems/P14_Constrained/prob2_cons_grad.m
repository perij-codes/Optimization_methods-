% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [g] = prob2_cons_grad(x)
    %evaluates gradient of constraints of problem 2
    gradc1 = 2*x;
    gradc2 = [0; x(3); x(2); -5*x(5); -5*x(4)];
    gradc3 = [3*x(1)^2; 3*x(2)^2; 0; 0; 0];
    g = [gradc1 gradc2 gradc3];
end