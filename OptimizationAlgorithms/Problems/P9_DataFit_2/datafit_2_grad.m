% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [g] = datafit_2_grad(x)
    y = [1.5, 2.25, 2.625];
    w = x(1);
    z = x(2);

    term1 = 0;
    term2 = 0;
    for i = 1:3
        s1 = -2*y(i)*(1-z^i)+(1-z^i)^2*2*w;
        term1 = term1+s1;

        s2 = -2*y(i)*w*(-i*z^(i-1)) + 2*w^2*(1-z^i)*(-i*z^(i-1));
        term2 = term2 + s2;
    end

%     term1 = -2*y(1)*(1-x(2)) + (1-x(2))^2*2*x(1)+ ...
%             -2*y(2)*(1-x(2)^2) + (1-x(2)^2)^2*2*x(1)+...
%             -2*y(3)*(1-x(2)^3) + (1-x(2)^3)^2*2*x(1);
%     
%     term2 = 2*y(1)*x(1) - 2*x(1)^2*(1-x(2)) + ...
%             2*y(2)*x(1)*(-2*x(2)) - 2*x(1)^2*(1-x(2)^2)*(-2*x(2)) + ...
%             2*y(3)*x(1)*(-3*x(2)^2) - 2*x(1)^2*(1-x(2)^3)*(-3*x(2)^2);
    
    g = [term1; term2];

end