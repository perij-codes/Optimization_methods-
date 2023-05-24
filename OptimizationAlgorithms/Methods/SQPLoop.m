% Code written by: Sruti Vutukury
% Algorithm theory discussed in "Numerical Optimization", J.Nocedal

function [xfinal,f,norm_c,fstore] = SQPLoop(problem,method,options)
    %{
    sequential quadratic programming algorithm
    inputs: problem, method, options
    outputs: final interate x, final function value f, other outputs
    %}

    %initialization
    xk = problem.x0bar;
    lamb = -problem.compute_cons_grad(xk)\problem.compute_obj_grad(xk); %lagrangian; scalar
    k = 0;
    fstore = [];
    termination1 = max(norm(lagrange_gradx(problem.x0bar,lamb,problem),"inf"),1); %termination condition 1, based on KKT conditions
    termination2 = max(norm(problem.compute_cons_func(problem.x0bar),"inf"),1); %termination condition 2, based on KKT conditions
    norm_c = norm(problem.compute_cons_func(xk),"inf"); %norm of the constraints

    while (k < 1000) && ((norm(lagrange_gradx(xk,lamb,problem),"inf") > method.eps*termination1)...
            || (norm(problem.compute_cons_func(xk),"inf") >  method.eps*termination2))
        
        Hk = eye(problem.n);
        %Hk = lagrange_hessxx(xk,lamb,problem);
        term2 = problem.compute_cons_grad(xk);
        term3 = term2';
        term4 = zeros(problem.m);

        % solve the subproblem
        A = [Hk term2; term3 term4];
        b = -[problem.compute_obj_grad(xk); problem.compute_cons_func(xk)'];
        direction = A\b; %(n+m)x(n+m)

        %update iterate
        xk = xk + direction(1:problem.n);
        lamb = direction(problem.n+1:end);

        f = problem.compute_obj_func(xk);
        fstore = [fstore f];
        norm_c = norm(problem.compute_cons_func(xk),"inf"); %2 or inf norm?
        k = k+1;

    end %end SQP loop with KKT term conditions

    xfinal = xk;
end