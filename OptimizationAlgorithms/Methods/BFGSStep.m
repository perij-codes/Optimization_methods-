% Code written by: Sruti Vutukury

function [x_new,phi_new,phi_grad_new,H_new] = BFGSStep(x,phi,phi_grad,H,nu,problem,method)
    %{
    implementation of a BFGS step, quasi-newton method
    uses either backtracking line search or weak wolfe linesearch
    inputs: problem, method, options
    %}
    
    d = -H*phi_grad;
    %calculate alpha to satisfy armijo-wolfe
    alpha = method.alphabar;
    while problem.compute_phi_func(x+alpha*d,nu,problem) > (phi + method.c1*alpha*phi_grad'*d)
        alpha = method.tau*alpha;
    end
    x_new = x + alpha*d; %alpha > 0 for Armijo-Wolfe conditions
    s = x_new - x;
    phi_grad_new = problem.compute_phi_grad(x_new,nu,problem);
    y = phi_grad_new-phi_grad;
    phi_new = problem.compute_phi_func(x_new,nu,problem);

    %compute BFGS Hessian approximation
    if s'*y <= method.bfgseps*norm(s')*norm(y')
        %skip update if s'k is not suff. positive
        H_new = H;
    else
        rho = 1/(s'*y);
        V = eye(length(x)) - rho*y*s';
        H_new = V'*H*V + rho*(s*s');
    end
end
