% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [x_new,phi_new,phi_grad_new,H_new] = NewtonModStep(x,phi,phi_grad,H,nu,problem,method)
    %{
    Modified Newton's Method
    uses backtracking line search
    %}
    
    %%subroutine to calculate eta: Hessian Approximation
    Hii = diag(H); 
    if min(Hii) > 0
        eta = 0;
    else
        eta = -min(Hii) + method.beta;
    end
    check = false;
    while check == false
        try chol(H+eta*eye(problem.n),"upper");
            %what should happen here?
            L = chol(H+eta*eye(problem.n),"upper");
            break
        catch
            eta = max(2*eta, method.beta);
        end
    end

    %%step update
    d = -(eta*eye(problem.n) + H)\phi_grad; %descent direction

    %%use bls to calculate alpha
    alpha = method.alphabar;
    while problem.compute_phi_func(x + alpha*d,nu,problem) > phi + method.c1*alpha*phi_grad'*d
        alpha = method.tau*alpha;
    end

    x_new = x + alpha*d;
    phi_new = problem.compute_phi_func(x_new,nu,problem);
    phi_grad_new = problem.compute_phi_grad(x_new,nu,problem);
    H_new = problem.compute_phi_hess(x_new,nu,problem);
end