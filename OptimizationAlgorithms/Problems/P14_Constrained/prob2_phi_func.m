% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [f] = prob2_phi_func(x,nu,problem)
    %evaluates phi of problem 2
    f = problem.compute_obj_func(x) + (nu/2)*norm(problem.compute_cons_func(x),2)^2;
end