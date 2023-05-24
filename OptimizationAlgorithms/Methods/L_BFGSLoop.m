% Code written by: Sruti Vutukury
% Code modified by: Janani Peri
% Algorithm theory discussed in "Numerical Optimization", J.Nocedal

function [x,f,outputs] = L_BFGSLoop(problem,method,options)
    %{
    L-BFGS
    uses either constant or backtracking line search
    
    inputs: problem, method, options
    outputs: final interate x, final function value f, other outputs
    %}

    % compute initial function/gradient/Hessian
    fc = 0; %counts how many function evaluations
    gc = 0; %counts how many gradient evaluations
    hc = 0; %counts how many hessian evaluations
    x = problem.x0;
    f = problem.compute_f(x); fc = fc + 1;
    g = problem.compute_g(x); gc = gc + 1;
    
    % alpha_low = method.alpha_low;
    % alpha_high = method.alpha_high;
    g0 = g;
    s = [];
    y = [];
    xhold = [];
    fhold = [];

    % set initial iteration counter
    tic
    tstart = cputime; %start timer to measure performance
    k = 0;
    while (k < options.max_iterations) && (norm(g,"inf") > 1e-6*max(norm(g0,"inf"),1)) && (toc < options.time_limit)
        H = eye(problem.n);
        
        %%L-BFGS two-loop recursion to calculate dk
        q = g;
        [~,col] = size(s);
    
        counter = 0;
        alpha_hold = [];
        rho_hold = [];
        if ~isempty(s) %k > method.m
            for i = 1:col %k-1:-1:k-method.m
                rho = 1/(s(:,col-counter)'*y(:,col-counter));
                alpha = rho*s(:,col-counter)'*q;
                q = q - alpha*y(:,col-counter);
                counter = counter+1;
                alpha_hold = [alpha alpha_hold];
                rho_hold = [rho rho_hold];
            end
        end
        r = (H*q);
        if ~isempty(s) %k > method.m
            for i = 1:col %k-method.m:1:k-1
                beta = rho_hold(i)*y(:,i)'*r;
                r = r + s(:,i)*(alpha_hold(i)-beta);
            end
        end
    
        d = -r;

        switch method.step_type
            case "Backtracking"
                %%use bls to calculate alpha
                alpha = method.alphabar;
                while problem.compute_f(x + alpha*d) > f + method.c1*alpha*g'*d
                    fc = fc + 1;
                    alpha = method.tau*alpha;
                end

            case "Constant"
                alpha = method.alphabar;
        end

        %iteration updates
        
        x_new = x + alpha*d;
        f_new = problem.compute_f(x_new); fc = fc + 1;
        g_new = problem.compute_g(x_new); gc = gc + 1;
        
        s_new = x_new-x;
        y_new = g_new-g;
        if s_new'*y_new <= method.eps*norm(s_new)*norm(y_new)
            %skip update if s'k is not suff. positive
        else
            s = [s s_new];
            y = [y y_new];
        end
        if k+1 > method.m && length(s) ~= 0
            s(:,1) = [];
            y(:,1) = [];
        end

        x = x_new; f = f_new; g = g_new;

        xhold = [xhold x_new];
        fhold = [fhold f_new];

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
