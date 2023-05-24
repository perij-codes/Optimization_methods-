% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [f] = exponential_10(x)
    n = length(x);
    term3 = 0;
    for i = 2:n
        term3 = term3 + (x(i)-1)^4;
    end
    f = (exp(x(1))-1)/(exp(x(1))+1) + 0.1*exp(-x(1)) + term3; 
end