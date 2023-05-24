% driver script for all constrained and unconstrained optimization problems and algorithms
% written by Sruti Vutukury
% algorithms from "Numerical Optimization" by J.Nocedal

close all; clear all; clc
addpath(pwd,"Problems/")
addpath(pwd,"Methods/")
T = table(); %stores outputs for analysis

%% unconstrained optimization problems and algorithms
% both line search and trust region algorithms

for p = 1:12 %iterate problem
    % set problem (minimal requirement: name of problem)
    problem.name = strcat("problem",string(p));
    disp(problem.name)
    for m = 1:11 %iterate method
        % set method (minimal requirement: name of method
        if m == 1
            method.name = 'GradientDescent';
            method.step_type = 'Backtracking';
            method.alpha = 1; %alpha0
            method.tau = 0.5;
            method.c1 = 1e-4;
        elseif m == 2
            method.name = 'GradientDescent';
            method.step_type = 'Weak_Wolfe_Linesearch';
            method.alpha = 1; %constant step size
            method.tau = 1e-3;
            method.c1 = 1e-4;
            method.c2 = 0.9;
            method.alpha_low =0;
            method.alpha_high =1000;
            method.c = 0.5;
        elseif m == 3
            method.name = 'NewtonMod';
            method.step_type = 'Backtracking';
            method.beta = 1e-6;
            method.alphabar = 1;
            method.tau = 0.5;
            method.c1 = 1e-4;
        elseif m == 4
            method.name = 'NewtonMod';
            method.step_type = 'Weak_Wolfe_Linesearch';
            method.beta = 1e-6;
            method.alphabar = 1;
            method.alpha = 1;
            method.tau = 1e-3;
            method.c1 = 1e-4;
            method.c2 = 0.9;
            method.alpha_low =0;
            method.alpha_high =1000;
            method.c = 0.5;
        elseif m == 5
            method.name = "TRNewtonCG";
            method.step_type = "";
            method.delta0 = 1; %1, overall bound on step lengths; if you increase NM should converge in one step
            method.c1 = 0.25; %from lecture notes
            method.c2 = 0.75; %piazza
        elseif m == 6
            method.name = "TRSR1CG";
            method.step_type = "";
            method.eps = 0.1; %convergence tolerance
            method.delta0 = 0.9; %overall bound on step lengths
            method.eta = 1e-2; % Algo 6.2
            method.r = 0.5; % Algo 6.2 % is this supposed to change?
        elseif m == 7
            method.name = "BFGS";
            method.step_type = "Backtracking";
            method.alphabar = 1;
            method.tau = 0.5;
            method.c1 = 1e-4;
        elseif m == 8
            method.name = "BFGS";
            method.step_type = 'Weak_Wolfe_Linesearch';
            method.beta = 1e-6;
            method.alphabar = 1;
            method.alpha = 1;
            method.tau = 1e-3;
            method.c1 = 1e-4;
            method.c2 = 0.9;
            method.alpha_low =0;
            method.alpha_high =1000;
            method.c = 0.5;
        elseif m == 9
            method.name = "DFP";
            method.step_type = "Backtracking";
            method.alphabar = 1;
            method.tau = 0.5;
            method.c1 = 1e-4;
        elseif m == 10
            method.name = "DFP";
            method.step_type = 'Weak_Wolfe_Linesearch';
            method.beta = 1e-6;
            method.alphabar = 1;
            method.alpha = 1;
            method.tau = 1e-3;
            method.c1 = 1e-4;
            method.c2 = 0.9;
            method.alpha_low =0;
            method.alpha_high =1000;
            method.c = 0.5;
        elseif m == 11
            method.name = "L-BFGS";
            method.step_type= "Backtracking";
            method.m = 2; %memory
            method.alphabar = 1; %constant step size
            method.tau = 0.5;
            method.c1 = 1e-4;
            method.eps = 1e-6;
        end
        disp(strcat(method.name,", ",method.step_type))
        
        % set options
        options.term_tol = 1e-6; %optimality tolerance
        options.max_iterations = 1e3; %iteration limit
        options.time_limit = 300; %terminate method after 5 minutes, useful for problems 3,4
        
        % run method and return x^* and f^*
        [x,f,outputs] = optSolver(problem,method,options);

        runname = strcat("problem",string(p),"_method",string(m),".mat");
        fullFileName = fullfile(pwd,"/Outputs", runname);
        save(fullFileName,"method","problem","options","x","f","outputs")

        T(p,m) = {f}; %export to excel
        
    end
end
filename = 'run_outputs.xlsx';
writetable(T,filename,'Sheet',1,'Range','D1')

%% constrained optimization
for p = 13:14 %iterate problem
    % set problem (minimal requirement: name of problem)
    if p == 13
        problem.name = 'problem13';
        problem.x0bar = [2;2];
        problem.xstar = [-1;-1];
        problem.n = 2; %number of decision vars
        problem.m = 1; %number of constraints
    elseif p == 14
        problem.name = 'problem14';
        problem.x0bar = [-1.8; 1.7; 1.9; -0.8; -0.8];
        problem.xstar = [-1.71; 1.59; 1.82; -0.763; -0.763];
        problem.n = 5; %number of decision vars
        problem.m = 3; %number of constraints
    end
    disp(problem.name)
    for m = 12:13 %iterate method
        if m == 12
            method.name = 'QPM';
            method.gamma = 10;
            method.nu0 = 1;
            method.eps = 1e-5; %for termination conditions
            for s = 1:2 %iterate subproblem solver
                %can expand to include other unconstraned subproblem solver
                if s == 1
                    method.subname = 'BFGS';
                    method.step_type = 'Backtracking';
                    method.bfgseps = 1e-6;
                    method.alphabar = 1;
                    method.tau = 0.5;
                    method.c1 = 1e-4;
                elseif s == 2
                    method.subname = 'GradientDescent';
                    method.step_type = 'Backtracking';
                    method.bfgseps = 1e-6; %subproblem solver termination
                    method.alpha = 1; %constant step size
                    method.tau = 0.5;
                    method.c1 = 1e-4;
                end
            end
        elseif m == 13
            method.name = 'SQP';
            method.step_type = 'Constant';
            method.alpha = 1e-3; %constant step size
            method.eps = 1e-5; %for termination conditions
        end        
        disp(strcat(method.name,", ",method.subname,", ", method.step_type))

        % set options
        options.term_tol = 1e-6;
        options.max_iterations = 1e3;
        
        % run method and return x^* and f^*
        [x,f,norm_c] = optSolver(problem,method,options);

        runname = strcat("problem",string(p),"_method",string(m),".mat");
        fullFileName = fullfile(pwd,"/Outputs", runname);
        save(fullFileName,"method","problem","options","x","f","outputs")

        T(p,m) = {f}; %export to excel
        
    end
end
filename = 'run_outputs.xlsx';
writetable(T,filename,'Sheet',2,'Range','D1')