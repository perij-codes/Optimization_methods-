% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [H] = exponential_10_hess(x)
    term1 = (2*exp(x(1))*(-exp(x(1))+1))/(exp(x(1))+1)^3 +0.1*exp(-x(1));
    %termrest = 12*(x-1).^2;
    termrest = 12*(x(2:end)-1).^2;
    v = [term1; termrest];
    H = diag(v,0);
end