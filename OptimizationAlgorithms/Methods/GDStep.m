% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [x_new,phi_new,phi_grad_new] = GDStep(x,phi,phi_grad,nu,problem,method)
    %{
    implementation of a Gradient Descent step
    uses either backtracking line search or constant step size
    modified to include slack variables
    %}

    d = -phi_grad; % search direction
    
    % determine step size
    switch method.step_type
        case 'Constant'
            alpha = method.alpha;
            x_new = x + alpha*d;
            phi_new = problem.compute_phi_func(x_new,nu,problem);
            phi_grad_new = problem.compute_phi_grad(x_new,nu,problem);
        case 'Backtracking'
            alpha = method.alpha;
            while problem.compute_phi_func(x + alpha*d,nu,problem) > phi + method.c1*alpha*phi_grad'*d
                alpha = method.tau*alpha;
            end
            x_new = x + alpha*d;
            phi_new = problem.compute_phi_func(x_new,nu,problem);
            phi_grad_new = problem.compute_phi_grad(x_new,nu,problem);
    end
end

