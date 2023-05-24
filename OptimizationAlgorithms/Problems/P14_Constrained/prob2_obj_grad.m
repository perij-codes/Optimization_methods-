% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [g] = prob2_obj_grad(x)
    %evaluates gradient of problem 2
%     s = exp(prod(x));
%     g = zeros(5,1);
%     g(1) = s*prod(x(2:end)) - (x(1)^3 + x(2)^3 + 1)*(3*x(1)^2);
%     g(2) = s*(x(1)*x(3)*x(4)*x(5)) - (x(1)^3 + x(2)^3 + 1)*(3*x(2)^2);
%     g(3) = s*(x(1)*x(2)*x(4)*x(5));
%     g(4) = s*(x(1)*x(2)*x(3)*x(5));
%     g(5) = s*(x(1)*x(2)*x(3)*x(4));

      x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4); x5 = x(5);
      g1 = x2*x3*x4*x5*exp(x1*x2*x3*x4*x5) - 3*x1^2*(x1^3 + x2^3 + 1);
      g2 = x1*x3*x4*x5*exp(x1*x2*x3*x4*x5) - 3*x2^2*(x1^3 + x2^3 + 1);
      g3 = x1*x2*x4*x5*exp(x1*x2*x3*x4*x5);
      g4 = x1*x2*x3*x5*exp(x1*x2*x3*x4*x5);
      g5 = x1*x2*x3*x4*exp(x1*x2*x3*x4*x5);
      g = [g1;g2;g3;g4;g5];
end