% Code written by: Sruti Vutukury

function [h] = lagrange_hessxx(x,lamb,problem)
    %evaluates gradient of the lagrangian wrt x
    %used specifically for constrained optimization algorithms: QPM, SQP
    
    h = problem.compute_obj_hess(x);
    hess_cons = problem.compute_cons_hess(x);

    if problem.m > 1
        for i = 1:length(lamb) %number of constraints
            h = h + lamb(i)*hess_cons(:,i);
        end
    else
        h = h + lamb*hess_cons;
    end
end