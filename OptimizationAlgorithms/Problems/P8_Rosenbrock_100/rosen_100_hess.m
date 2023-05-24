%written by svutukury

function [H] = rosen_100_hess(x)
    %returns gradient value of the rosenbrock_100 problem
    %x = 100 x 100

    %unvectorized! this works
%     H = eye(length(x));
%     for i = 2:99
%         ii = 202-400*(x(i+1)-3*x(i)^2);
%         ip1i = -400*x(i);
%         im1i = -400*x(i-1);
%         H(i,i) = ii;
%         H(i-1,i) = im1i;
%         H(i+1,i) = ip1i;
%     end
%     H(1,1) = 2 - 400*x(2)+1200*x(1)^2;
%     H(2,1) = -400*x(1);

    %vectorized, this works!
    n = length(x);
    H = eye(n);
    ii = 202 - 400*(x(3:n).^2 - x(2:n-1)).';
    ip1i = -400*x(2:n-1).';
    im1i = -400*x(3:n).';
    %H(2:n-1,2:n-1) = diag(ii) + diag(ip1i,-1) + diag(im1i,1);
    H(1,1) = 2 - 400*x(2) + 1200*x(1)^2;
    H(2:n,2:n) = diag([H(1,1) ii]) + diag(ip1i,-1) + diag(im1i,1);
    H(2,1) = -400*x(1);
    H(1,2) = H(2,1);
end