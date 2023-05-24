% Code written by: Sruti Vutukury

function [g] = lagrange_gradx(x,lamb,problem)
    %evaluates gradient of the lagrangian wrt x
    %used specifically for constrained optimization algorithms: QPM, SQP

    g = problem.compute_obj_grad(x);
    grad_cons = problem.compute_cons_grad(x);
    for i = 1:length(lamb) %number of constraints
        g = g + lamb(i)*grad_cons(:,i);
    end
end