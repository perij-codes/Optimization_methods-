% Code written by: Sruti Vutukury
% Code modified by: Janani Peri
% Algorithm theory discussed in "Numerical Optimization", J.Nocedal

function [x,f,outputs] = GDLoop(problem,method,options)
    %{
    Gradient Descent
    uses either constant, backtracking line search, or weak wolfe linesearch
    
    inputs: problem, method, options
    outputs: final interate x, final function value f, other outputs
    %}

    fc = 0; %counts how many function evaluations
    gc = 0; %counts how many gradient evaluations
    hc = 0; %counts how many hessian evaluations

    x = problem.x0;
    f = problem.compute_f(x); fc = fc+1;
    g = problem.compute_g(x); gc = gc+1;
    g0 = g;
    xhold = [];
    fhold = [];
    
    tic
    tstart = cputime; %start timer to measure performance
    k = 0;
    while (k < options.max_iterations) && (norm(g,"inf") > 1e-6*max(norm(g0,"inf"),1)) && (toc < options.time_limit)
        d = -g; % search direction
        % determine step size
        switch method.step_type
            case 'Constant'
                alpha = method.alpha;
            
            case 'Backtracking'
                alpha = method.alpha;
                while problem.compute_f(x + alpha*d) > f + method.c1*alpha*g'*d
                    fc = fc+1;
                    alpha = method.tau*alpha;
                end
    
            case 'Weak_Wolfe_Linesearch'
                alpha = method.alpha;
%                 x_new = x + alpha*d;
%                 f_new = problem.compute_f(x_new); fc = fc+1;
%                 g_new = problem.compute_g(x_new); gc = gc+1;
                while 1 && (toc < options.time_limit)
                    %disp('inside while')
                    if (problem.compute_f(x+alpha*d) <= f + method.c1*alpha*g'*d) 
                       fc = fc+1;
                       if (problem.compute_g(x+alpha*d)'*d >= method.c2*g'*d)
                           gc = gc + 1;
                           break
                       end
                    end    
                    if problem.compute_f(x+alpha*d) <= f + method.c1*alpha*g'*d
                       fc = fc + 1;
                       method.alpha_low = alpha;
                    else 
                       method.alpha_high = alpha;
                    end
                    alpha = method.c*method.alpha_low + (1-method.c)*method.alpha_high;
                end

        end
        x_new = x + alpha*d;
        f_new = problem.compute_f(x_new); fc = fc+1;
        g_new = problem.compute_g(x_new); gc = gc+1;
        
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

