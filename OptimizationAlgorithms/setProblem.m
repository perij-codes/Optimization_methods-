% Code written by: Albert S. Berahas
% Modified by Sruti Vutukury

function [problem] = setProblem(problem)
    % Function that specifies the problem. Specifically, a way to compute: 
    %    (1) function values; (2) gradients; and, (3) Hessians (if necessary).
    %
    %           Input: problem (struct), required (problem.name)
    %           Output: problem (struct)
    %
    % Error(s): 
    %       (1) if problem name not specified;
    %       (2) function handles (for function, gradient, Hessian) not found
    %

    % check is problem name available
    if ~isfield(problem,'name')
        error('Problem name not defined!!!')
    end
    
    % set function handles according the the selected problem
    switch problem.name

        case 'problem1'     
            problem.compute_f = @quad_10_10_func;
            problem.compute_g = @quad_10_10_grad;
            problem.compute_h = @quad_10_10_Hess;
            rng(0)
            problem.x0=20*rand(10,1)-10;
            problem.n = length(problem.x0);
               
        case 'problem2'     
            problem.compute_f = @quad_10_1000_func;
            problem.compute_g = @quad_10_1000_grad;
            problem.compute_h = @quad_10_1000_Hess;
            rng(0)
            problem.x0=20*rand(10,1)-10;
            problem.n = length(problem.x0);   
            
        case 'problem3'     
            problem.compute_f = @quad_1000_10_func;
            problem.compute_g = @quad_1000_10_grad;
            problem.compute_h = @quad_1000_10_Hess;
            rng(0)
            problem.x0=20*rand(1000,1)-10;
            problem.n = length(problem.x0);
           
        case 'problem4'     
            problem.compute_f = @quad_1000_1000_func;
            problem.compute_g = @quad_1000_1000_grad;
            problem.compute_h = @quad_1000_1000_Hess; 
            rng(0)
            problem.x0=20*rand(1000,1)-10;
            problem.n = length(problem.x0);
            
        case 'problem5'     
            problem.compute_f = @quartic_1_func;
            problem.compute_g = @quartic_1_grad;
            problem.compute_h = @quartic_1_Hess;
            problem.x0=[cos(70) sin(70) cos(70) sin(70)]' ;
            problem.n = length(problem.x0);
                 
        case 'problem6'     
            problem.compute_f = @quartic_2_func;
            problem.compute_g = @quartic_2_grad;
            problem.compute_h = @quartic_2_Hess;
            problem.x0=[cos(70) sin(70) cos(70) sin(70)]' ;
            problem.n = length(problem.x0);
           
        case 'problem7'     
            problem.compute_f = @rosen_2_func;
            problem.compute_g = @rosen_2_grad;
            problem.compute_h = @rosen_2_hess;
            problem.x0 = [-1.2; 1];
            problem.n = 2;
            problem.n = length(problem.x0);
            
        case "problem8"
            problem.compute_f = @rosen_100_func;
            problem.compute_g = @rosen_100_grad;
            problem.compute_h = @rosen_100_hess;
            problem.n = 100;
            problem.x0 = ones(problem.n,1);
            problem.x0(1) = -1.2;
    
        case "problem9"
            problem.compute_f = @datafit_2_func;
            problem.compute_g = @datafit_2_grad;
            problem.compute_h = @datafit_2_hess;
            problem.y = [1.5; 2.25; 2.625];
            problem.n = 2;
            problem.x0 = [1; 1];

        case "problem10"
            problem.compute_f = @exponential_10;
            problem.compute_g = @exponential_10_grad;
            problem.compute_h = @exponential_10_hess;
            problem.n = 10;
            problem.x0 = zeros(problem.n,1);
            problem.x0(1) = 1;
    
         case "problem11"
            problem.compute_f = @exponential_1000;
            problem.compute_g = @exponential_1000_grad;
            problem.compute_h = @exponential_1000_hess;
            problem.n = 100;
            problem.x0 = zeros(problem.n,1);
            problem.x0(1) = 1;

        case "problem12"
            problem.compute_f = @genhumps_5_func;
            problem.compute_g = @genhumps_5_grad;
            problem.compute_h = @genhumps_5_Hess;
            problem.n = 5;
            problem.x0 = 506.2*ones(problem.n,1);
            problem.x0(1) = problem.x0(1)*-1;

        case 'problem13'
            problem.compute_phi_func = @prob1_phi_func; %phi, used for subproblem
            problem.compute_phi_grad = @prob1_phi_grad;  
            problem.compute_phi_hess = @prob1_phi_hess;
    
            problem.compute_obj_func = @prob1_obj_func; %just f
            problem.compute_obj_grad = @prob1_obj_grad;
            problem.compute_obj_hess = @prob1_obj_hess;
    
            problem.compute_cons_func = @prob1_cons_func; %constraints
            problem.compute_cons_grad = @prob1_cons_grad;
            problem.compute_cons_hess = @prob1_cons_hess;

        case 'problem14'
            problem.compute_phi_func = @prob2_phi_func; %phi, used for subproblem
            problem.compute_phi_grad = @prob2_phi_grad;  
            problem.compute_phi_hess = @prob2_phi_hess;
    
            problem.compute_obj_func = @prob2_obj_func; %just f
            problem.compute_obj_grad = @prob2_obj_grad;
            problem.compute_obj_hess = @prob2_obj_hess;
    
            problem.compute_cons_func = @prob2_cons_func; %constraints
            problem.compute_cons_grad = @prob2_cons_grad;
            problem.compute_cons_hess = @prob2_cons_hess;

        otherwise
            error('Problem not defined!!!')
end