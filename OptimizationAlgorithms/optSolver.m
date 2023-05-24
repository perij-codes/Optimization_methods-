% Code written by: Albert S. Berahas
% Code modified by: Sruti Vutukury

function [x,f,outputs] = optSolver(problem,method,options)
    addpath(pwd,'Project_Problems_MATLAB')
    addpath(pwd,'Methods')
    % Inputs: problem, method, options (structs)
    % Outputs: final iterate (x), final function value (f)

    % set problem, method and options
    [problem] = setProblem(problem);
    [method] = setMethod(method);
    [options] = setOptions(options);
    
    if method.name == "GradientDescent" %methods 1,2
        [x,f,outputs] = GDLoop(problem,method,options);
    elseif method.name == "NewtonMod" %methods 3,4
        [x,f,outputs] = NewtonModLoop(problem,method,options);
    elseif method.name == "TRNewtonCG" %method 5
        [x,f,outputs] = TRNewtonCGLoop(problem,method,options);
    elseif method.name == "TRSR1CG" %method 6
        [x,f,outputs] = TRSR1CGLoop(problem,method,options);
    elseif method.name == "BFGS" %method 7,8
        [x,f,outputs] = BFGSLoop(problem,method,options);
    elseif method.name == "DFP" %method 9,10
        [x,f,outputs] = DFPLoop(problem,method,options);
    elseif method.name == "L-BFGS" %method 11
        [x,f,outputs] = L_BFGSLoop(problem,method,options);
    elseif method.name == "QPM" %method 12
        [x,f,outputs] = QPMLoop(problem,method,options);
    elseif method.name == "SQP" %method 13
        [x,f,outputs] = SQPLoop(problem,method,options);
    end
end