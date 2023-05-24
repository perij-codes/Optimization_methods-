% IOE 511/MATH 562, University of Michigan
% Code written by: Sruti Vutukury

function [H] = datafit_2_hess(x)
    y = [1.5, 2.25, 2.625];
    w = x(1);
    z = x(2);

    term1 = 0;
    term2 = 0;
    term3 = 0;
    term4 = 0;

    for i = 1:3
        s1 = 2*(1+z^i)^2;
        term1 = term1+s1;

        s2 = -2*y(i)*i*z^(i-1)-4*i*w*(1-z^i)*(z^(i-1));
        term2 = term2 + s2;

        s3 = 2*y(i)*(i*z^(i-1))-4*i*w*(1-z^i)*(z^(i-1));
        term3 = term3 + s3;

        s4 = 2*y(i)*w*(i*(i-1)*z^(i-2)) - 2*w^2*(1-z^i)*(i*(i-1)*z^(i-2)) + ...
            (i*z^(i-1))*(2*w^2)*(i*z^(i-1));
        term4 = term4 + s4;

    end

%     term1 = 2*(1+x(2))^2 + 2*(1+x(2)^2)^2 + 2*(1+x(2)^3)^2;
%     term2 = 2*y(1)-4*x(1)*(1-x(2)) + 4*y(2)*x(2)-8*x(1)*(1-x(2)^2)*x(2) + ...
%             6*y(3)*(x(2)^2)-12*x(1)*(1-x(2)^3)*(x(2)^2);
%     term3 = 2*y(1)-4*x(1)*(1-x(2)) + 4*y(2)*x(2)-8*x(1)*(1-x(2)^2)*x(2) + ...
%             6*y(3)*(x(2)^2)-12*x(1)*(1-x(2)^3)*(x(2)^2);
%     term4i1 = 2*(x(2)^2);
%     term4i2 = 4*y(2)*x(1)-4*(x(2)^2)*(1-x(1)^2)+(8*x(2)^2*(x(1)^2));
%     term4i3 = 12*y(3)*x(1)*x(2)-12*(x(1)^2)*(1-x(2)^3)*x(2)+(18*x(2)^4*(x(1)^2));
%     term4 = term4i3 + term4i2 + term4i1;
    
    H = [term1 term2; term3 term4];

end