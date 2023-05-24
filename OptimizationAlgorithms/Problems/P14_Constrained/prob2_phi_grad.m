% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [g] = prob2_phi_grad(x,nu,problem)
    %problem 2: computes grad_x of phi
    
    %g = problem.compute_grad(x) + nu*problem.compute_func_cons(x);
    cons = problem.compute_cons_func(x); %gradient of the constraints 1*3
    grad_cons = problem.compute_cons_grad(x); %gradient of the constraints 5*3
    g = problem.compute_obj_grad(x) + ...
        nu*cons(1)*grad_cons(:,1) +....
        nu*cons(2)*grad_cons(:,2) +....
        nu*cons(3)*grad_cons(:,3);
end