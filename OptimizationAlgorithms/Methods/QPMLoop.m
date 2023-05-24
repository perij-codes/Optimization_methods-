% Code written by: Sruti Vutukury
% Algorithm theory discussed in "Numerical Optimization", J.Nocedal

function [xfinal,f,norm_c,fstore] = QPMLoop(problem,method,options)
    %{
    quadratic penality method for constrained optimization
    inputs: problem, method, options
    outputs: final interate x, final function value f, other outputs
    %}

    % initialization 
    xbar = problem.x0bar; %initial point
    nu = method.nu0; %initial penalty parameter
    chi = 1/nu; %initial slack variable

    lamb = -problem.compute_cons_grad(xbar)\problem.compute_obj_grad(xbar); %lagrangian; scalar output
    norm_c = norm(problem.compute_cons_func(xbar),"inf"); %norm of the constraint

    %approximate solution xk to the unconstrained optimization subproblem
    k = 0;
    term1 = max(norm(lagrange_gradx(problem.x0bar,lamb,problem),"inf"),1); %termination condition 1, based on KKT conditions
    term2 = max(norm(problem.compute_cons_func(problem.x0bar),"inf"),1); %termination condition 2, based on KKT conditions
    fstore = [];
    f = problem.compute_obj_func(xbar);

    %keep searching until you get a point that satisfies BOTH KKT conditions
    %if one of the KKT conditions is not met, continue
    while (k < 1000) && ((norm(lagrange_gradx(xbar,lamb,problem),"inf") > method.eps*term1)...
            || (norm(problem.compute_cons_func(xbar),"inf") >  method.eps*term2))

        phi = problem.compute_phi_func(xbar,nu,problem);
        grad_phi0 = problem.compute_phi_grad(xbar,nu,problem);
        grad_phi = grad_phi0;
        j = 0;
        x = xbar;

        %solve subproblem
        H = eye(problem.n); %change and see for problem2
        nu_kp1 = nu*method.gamma;
        chi_kp1 = 1/nu_kp1;

        while (j < 200) && norm(grad_phi,"inf") > min(method.bfgseps*max(norm(grad_phi0,"inf"),1),chi_kp1)
        %while (j < 100) && (norm(grad_phi,"inf") > chi_kp1)
            if method.subname == "BFGS"
                [x_new,phi_new,grad_phi_new,H_new] = BFGSStep(x,phi,grad_phi,H,nu,problem,method);
                H = H_new;
            elseif method.subname == "GradientDescent"
                [x_new,phi_new,grad_phi_new] = GDStep(x,phi,grad_phi,nu,problem,method);
            elseif method.subname == "NewtonModStep"
                [x_new,phi_new,grad_phi_new,H_new] = NewtonModStep(x,phi,grad_phi,H,nu,problem,method);
                H = H_new;
            end
            
            j = j + 1;
            x = x_new;
            phi = phi_new;
            grad_phi = grad_phi_new;
        end
            
        %update iterate
        xbar = x;
        nu = nu_kp1;
        chi = chi_kp1;
        k = k + 1;

        nu_kp1 = nu*method.gamma;
        chi_kp1 = 1/nu_kp1;

        f = problem.compute_obj_func(xbar); %objective function value
        fstore = [fstore f];
        lamb = -problem.compute_cons_grad(xbar)\problem.compute_obj_grad(xbar);
        norm_c = norm(problem.compute_cons_func(xbar),"inf"); %2 or inf norm?

    end %end outer loop QPM

    xfinal = xbar;
end