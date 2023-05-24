% Code written by: Sruti Vutukury
% Code modified by: Janani Peri
% Algorithm theory discussed in "Numerical Optimization", J.Nocedal

function [x,f,outputs] = BFGSLoop(problem,method,options)
    %{
    implementation of BFGS, quasi-newton method
    uses either backtracking line search or weak wolfe linesearch
    inputs: problem, method, options
    outputs: final interate x, final function value f, other outputs
    %}

    fc = 0; %counts how many function evaluations
    gc = 0; %counts how many gradient evaluations
    hc = 0; %counts how many hessian evaluations

    x = problem.x0;
    f = problem.compute_f(x); fc = fc+1;
    g = problem.compute_g(x); gc = gc+1;
    H = eye(problem.n);
    g0 = g;
    xhold = [];
    fhold = [];
    
    tic
    tstart = cputime; %start timer to measure performance
    k = 0;
    while (k < options.max_iterations) && (norm(g,"inf") > 1e-6*max(norm(g0,"inf"),1)) && (toc < options.time_limit)
        d = -H*g;

        switch method.step_type
            case "Backtracking"
                %calculate alpha to satisfy armijo-wolfe
                alpha = method.alphabar;
                while problem.compute_f(x-alpha*g) > (f - method.c1*alpha*norm(g)^2)
                    fc = fc+1;
                    alpha = method.tau*alpha;
                end
    
            case "Weak_Wolfe_Linesearch"
                %calculate alpha to satisfy Weak Wolfe line search conditions:
                alpha = method.alpha;
                while 1 && (toc < options.time_limit)
                    
                    if problem.compute_f(x + alpha*d) <= f + method.c1*alpha*g'*d
                        fc = fc+1;
                        if (problem.compute_g(x + alpha*d))'*d >= method.c2*g'*d
                             gc = gc+1;
                            break
                        end    
                    end
                    if problem.compute_f(x + alpha*d) <= f + method.c1*alpha*g'*d
                       method.alpha_low = alpha;
                       fc = fc+1;
    
                    else 
                       method.alpha_high = alpha;
                    end
                    alpha = method.c*method.alpha_low + (1-method.c)*method.alpha_high;
                end
        end

        x_new = x + alpha*d; %alpha > 0 for Armijo-Wolfe conditions
        s = x_new - x;
        g_new = problem.compute_g(x_new); gc = gc+1;
        y = g_new-g;
        f_new = problem.compute_f(x_new); fc = fc+1;
    
        %compute BFGS Hessian approximation
        rho = 1/(s'*y);
        V = eye(length(x)) - rho*y*s';
        H_new = V'*H*V + rho*(s*s');
    
        % iteration updates
        x = x_new; f = f_new; g = g_new;
        xhold = [xhold x];
        fhold = [fhold f];
        k = k + 1;
    end

    %final outputs after termination
    outputs.xhold = xhold;
    outputs.fhold = fhold;
    outputs.k = k; %iterations
    outputs.time = cputime - tstart; %time to find a solution
    outputs.evals = 0; %total number of gradient and hessian evaluations
    outputs.fc = fc; %number of function evaluations
    outputs.gc = gc; %number of gradient evaluations
    outputs.hc = hc; %number of hessian evaluations

end
