% Code written by: Sruti Vutukury

function [L] = lagrange(x,lamb,problem)
    %evaluates lagrangian, returns scalar
    %used specifically for constrained optimization algorithms: QPM, SQP
    L = problem.compute_obj_func(x) + lamb'*problem.compute_func_cons(x);
end