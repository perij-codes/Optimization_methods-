% Code written by: Sruti Vutukury
% Algorithm theory discussed in "Numerical Optimization", J.Nocedal

function [x,f,outputs] = TRNewtonCGLoop(problem,method,options)
    %{
    trust region Newton with CG subproblem solver, uses the true Hessian
    Inputs: problem, method, options (structs)
    Outputs: final iterate (x), final function value (f)
    %}

    % compute initial function/gradient/Hessian
    fc = 0; %counts how many function evaluations
    gc = 0; %counts how many gradient evaluations
    hc = 0; %counts how many hessian evaluations
    x = problem.x0;
    f = problem.compute_f(x); fc = fc + 1;
    g = problem.compute_g(x); gc = gc + 1;
    g0 = g;
    delta = method.delta0; %initial step length, must be between (0,max step length)
    xhold = [];
    fhold = [];
    Bk = problem.compute_h(x); hc = hc + 1; %true Hessian
    hc = hc + 1;
    mk0 = f; %model value at start

    % set initial iteration counter
    tic
    tstart = cputime; %start timer to measure performance
    k = 0;
    while (k < options.max_iterations) && (norm(g,"inf") > 1e-6*max(norm(g0,"inf"),1)) && (toc < options.time_limit)
        
        %%%%% get step pk via solving TR subproblem with CG, Algo 7.2
        pk = CG_subproblem_solver(g,Bk,delta,problem,options);
        mkpk = f + g'*pk + 0.5*pk'*Bk*pk;

        %%%%% TR iteration update, Algo 4.1, Nocedal
        f_new = problem.compute_f(x+pk); fc = fc + 1;
        rhok = (f - f_new)/(mk0-mkpk); %solve for rho from Eqn 4.4
        
        %if-else logic from TR algo slides, simple
        if rhok > method.c1
            x_new = x + pk;
            if (rhok > method.c2)
                delta = 2.25*delta;
            end
        else
            x_new = x;
            delta = 0.45*delta;
        end
        
        f_new = problem.compute_f(x+pk); fc = fc + 1;
        mk0 = f_new;
        g_new = problem.compute_g(x+pk); gc = gc + 1;
        Bk = problem.compute_h(x+pk); hc = hc + 1;
        x = x_new;
        f = f_new;
        g = g_new;
        xhold = [xhold x];
        fhold = [fhold f];
        k = k + 1;
    end

    %final outputs after termination
    outputs.xhold = xhold;
    outputs.fhold = fhold;
    outputs.k = k; %iterations
    outputs.time = tstart-cputime; %time to find a solution
    outputs.evals = 0; %total number of gradient and hessian evaluations
    outputs.fc = fc; %number of function evaluations
    outputs.gc = gc; %number of gradient evaluations
    outputs.hc = hc; %number of hessian evaluations
end

function pk = CG_subproblem_solver(g,H,delta,problem,options)
    %{
    solve the TR subproblem approximately --> pk, Algorithm 7.2, Numerical Optimization, Nocedal
    inputs: gradient g, Hessian H
    outputs: step length pk
    %}

    zj = zeros(problem.n,1); %iterate generated by CG solver
    rj = g; %parameter 
    dj = -g; %search direction
    Bk = H; %Hessian
    epsk = 1e-4; %norm(g,2); from lecture %termination tolerance for CG; aka options.term_tol_CG 
    found = false; %did not find a pk, continue j loop

    if norm(rj,2) < epsk
        pk = zeros(problem.n,1);
        found = true;
    else
        j= 0;
        while j < options.max_iterations && ~found %if a pk is found, leave CG solver
            if (dj'*Bk)*dj <= 0
                % stops the method if its current search direction
                % dj is a direction of nonpositive curvature along Bk
                % squared pk = zj + tau*dj and took positive root to find tau
                num = -2*zj'*dj + sqrt((2*zj'*dj)^2-4*(dj'*dj)*(zj'*zj-delta^2));
                den = 2*(dj'*dj);
                tau = num/den;
                pk = zj + tau*dj;
                found = true;
                break
            end

            alphaj = (rj'*rj)/((dj'*Bk)*dj);
            zjp1 = zj + alphaj*dj;
            %zjp1 = sum(alphaj*dj);
            
            if norm(zjp1,2) >= delta
                % stops the method if zj+1 violates the trust-region bound
                num = -2*zj'*dj + sqrt((2*zj'*dj)^2-4*(dj'*dj)*(zj'*zj-delta^2));
                den = 2*(dj'*dj);
                tau = num/den;
                pk = zj + tau*dj;
                found = true;
                break
            end

            rjp1 = rj + alphaj*Bk*dj;
            if norm(rjp1,2) <= epsk
                pk = zjp1;
                found = true;
                break
            end

            zj = zjp1;
            betajp1 = (rjp1'*rjp1)/(rj'*rj);
            rj = rjp1;
            djp1 = -rjp1 + betajp1*dj;
            dj = djp1;
            j  = j+1;
        end %end j iters 
    end
    if found == false
       check2 = 1;
    end
end
