% Code written by: Sruti Vutukury
% Code modified by: Janani Peri
% Algorithm theory discussed in "Numerical Optimization", J.Nocedal

function [x,f,outputs] = NewtonModLoop(problem,method,options)
    %{
    Modified Newton's Method
    uses either backtracking line search or weak wolfe line search
    
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
    g0 = g;
    H = problem.compute_h(x); hc = hc + 1;
    xhold = [];
    fhold = [];

    % set initial iteration counter
    tic
    tstart = cputime; %start timer to measure performance
    k = 0;
    while (k < options.max_iterations) && (norm(g,"inf") > 1e-6*max(norm(g0,"inf"),1)) && (toc < options.time_limit)
        switch method.step_type
            case 'Backtracking'
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
                d = -(eta*eye(problem.n) + H)\g; %descent direction
    
                %%use bls to calculate alpha
                alpha = method.alphabar;
                while problem.compute_f(x + alpha*d) > f + method.c1*alpha*g'*d
                    fc = fc + 1;
                    alpha = method.tau*alpha;
                end
    
                x_new = x + alpha*d;
                f_new = problem.compute_f(x_new); fc = fc + 1;
                g_new = problem.compute_g(x_new); gc = gc + 1;
                H_new = problem.compute_h(x_new); hc = hc + 1;
    
            case 'Weak_Wolfe_Linesearch'
                d = -inv(H)*g; %check if we need this
    
                alpha = method.alpha;
                while 1 && (toc < options.time_limit)
%                     x_new = x + alpha*d;
%                     f_new = problem.compute_f(x_new); fc = fc + 1;
%                     g_new = problem.compute_g(x_new); gc = gc + 1;
%                     H_new = problem.compute_h(x_new); hc = hc + 1;
%     
                    %disp('inside while')
                    if problem.compute_f(x + alpha*d) <= f + method.c1*alpha*g'*d
                        %disp('inside first if')
                        fc = fc+1;
                        if (problem.compute_g(x + alpha*d))'*d >= method.c2*g'*d
                             %disp('breaking critireon')
                             gc = gc+1;
                            break
                        end    
                    end
                    if problem.compute_f(x + alpha*d) <= f + method.c1*alpha*g'*d
                       method.alpha_low = alpha;
                       fc = fc + 1;
                    else 
                       method.alpha_high = alpha;
                    end
                    alpha = method.c*method.alpha_low + (1-method.c)*method.alpha_high;
                end
                x_new = x + alpha*d;
                f_new = problem.compute_f(x_new); fc = fc + 1;
                g_new = problem.compute_g(x_new); gc = gc + 1;
                H_new = problem.compute_h(x_new); hc = hc + 1;
        end %end switch case

        % iteration updates
        x = x_new; f = f_new; g = g_new; H = H_new;
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