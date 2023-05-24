% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [h] = prob2_cons_hess(x)
    %evaluates hessian  of constraints of problem 2
    [g] = prob2_cons_grad(x);
    h1 = gradient(g(:,1));
    h2 = gradient(g(:,2));
    h3 = gradient(g(:,3));

    h = [h1, h2, h3];
end