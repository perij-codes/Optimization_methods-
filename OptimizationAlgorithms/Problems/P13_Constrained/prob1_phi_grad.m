% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [g] = prob1_phi_grad(x,nu,problem)
    %evaluates phi gradient_x of problem 1

    g = problem.compute_obj_grad(x) + ...
        nu*problem.compute_cons_func(x)*problem.compute_cons_grad(x);
    
    %hardcoded
%     x1 = x(1);
%     x2 = x(2);
%     g = [(nu*(2*x1^3 + 2*x1*x2^2 -4*x1) +1); ...
%          (nu*(2*x2*x1^2 + 2*x2^3 -4*x2) +1)];

end