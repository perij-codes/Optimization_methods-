% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [H] = prob2_obj_hess(x)
    %evaluates hessian of objective function problem 2
    H = eye(5);
%     s = exp(prod(x));
%     H(1,1) = prod(x(2:5))^2*s - (x(1)^3 + x(2)^3 + 1)*(6*x(1)) -(3*x(1)^2)^2;
%     H(2,1) = s*(prod(x(3:5))) + prod(x(2:5))*s*prod(x(1),x(3:5)) - 9*x(1)^2*x(2)^2;
%     H(3,1) = s*(x(2)*x(4)*x(5)) + prod(x(2:5))*s*(x(1)*x(2)*x(4)*x(5));
%     H(4,1) = s*(x(2)*x(3)*x(5)) + prod(x(2:5))*s*(x(1)*x(2)*x(3)*x(5));
%     H(5,1) = s*(x(2)*x(3)*x(4)) + prod(x(2:5))*s*(x(1)*x(2)*x(3)*x(4));
% 
%     H(1,2) = s*prod(x(1),x(3:5),x(2:5)) - 9*x(2)^2*x(1)^2;
%     H(2,2) = s*prod(x(1),x(3:5))^2 -9*x(2)^4 - 6*x(2)*(x(1)^3+x(2)^3+1);
%     H(3,2) = s*x(1)*prod(x(4:5))+prod(x(1),x(3:5),x(1:2),x(4:5),s);
%     H(4,2) = s*x(1)*x(3)*x(5)+prod(x(1),x(3:5),x(1:3),x(5),s);
%     H(5,2) = s*x(1)*prod(x(3:4))+prod(x(1),x(3:5),x(1:4),s);

%     if issymmetric(H)
%         disp("hessian is symmetric")
%     else
%         disp("hessian not symmetric")
%     end

    % these equations are those derived via symbolic matlab on the driver script
    x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4); x5 = x(5);
    H(1,1) = x2^2*x3^2*x4^2*x5^2*exp(x1*x2*x3*x4*x5) - 6*x1*(x1^3 + x2^3 + 1) - 9*x1^4;
    H(1,2) = x3*x4*x5*exp(x1*x2*x3*x4*x5) - 9*x1^2*x2^2 + x1*x2*x3^2*x4^2*x5^2*exp(x1*x2*x3*x4*x5);
    %H(1,3) = x3*x4*x5*exp(x1*x2*x3*x4*x5) - 9*x1^2*x2^2 + x1*x2*x3^2*x4^2*x5^2*exp(x1*x2*x3*x4*x5);
    H(1,3) = x2*x4*x5*exp(x1*x2*x3*x4*x5) + x1*x2^2*x3*x4^2*x5^2*exp(x1*x2*x3*x4*x5);
    H(1,4) = x2*x3*x5*exp(x1*x2*x3*x4*x5) + x1*x2^2*x3^2*x4*x5^2*exp(x1*x2*x3*x4*x5);
    H(1,5) = x2*x3*x4*exp(x1*x2*x3*x4*x5) + x1*x2^2*x3^2*x4^2*x5*exp(x1*x2*x3*x4*x5);
    
    H(2,2) = x1^2*x3^2*x4^2*x5^2*exp(x1*x2*x3*x4*x5) - 6*x2*(x1^3 + x2^3 + 1) - 9*x2^4;
    H(2,3) = x1*x4*x5*exp(x1*x2*x3*x4*x5) + x1^2*x2*x3*x4^2*x5^2*exp(x1*x2*x3*x4*x5);
    H(2,4) = x1*x3*x5*exp(x1*x2*x3*x4*x5) + x1^2*x2*x3^2*x4*x5^2*exp(x1*x2*x3*x4*x5);
    H(2,5) = x1*x3*x4*exp(x1*x2*x3*x4*x5) + x1^2*x2*x3^2*x4^2*x5*exp(x1*x2*x3*x4*x5);

    H(3,3) = x1^2*x2^2*x4^2*x5^2*exp(x1*x2*x3*x4*x5);
    H(3,4) = x1*x2*x5*exp(x1*x2*x3*x4*x5) + x1^2*x2^2*x3*x4*x5^2*exp(x1*x2*x3*x4*x5);
    H(3,5) = x1*x2*x4*exp(x1*x2*x3*x4*x5) + x1^2*x2^2*x3*x4^2*x5*exp(x1*x2*x3*x4*x5);
    
    H(4,4) = x1^2*x2^2*x3^2*x5^2*exp(x1*x2*x3*x4*x5);
    H(4,5) = x1*x2*x3*exp(x1*x2*x3*x4*x5) + x1^2*x2^2*x3^2*x4*x5*exp(x1*x2*x3*x4*x5);

    H(5,5) = x1^2*x2^2*x3^2*x4^2*exp(x1*x2*x3*x4*x5);

    H = H'+triu(H',1)'; %flips the upper triangle to the lower triangle to make the hessian symmetric
end