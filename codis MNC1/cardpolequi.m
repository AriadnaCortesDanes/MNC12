%% Code 3 : Lagrange cardinal polynomials (equispaced nodes)
% Input:    1. a & b: interval [a,b]
%           2. n: with x_0 = a and x_n = b (n+1 interp nodes)
%           3. m: number of points in dense grid
%
% Output:   1. P: matrix of card. polynomials.
%           2. xn: equispaced interpolation nodes
%           3. z: dense grid
function [P,xn,z] = cardpolequi(a, b, n, m)
    xn = a + [0:n]'*(b-a)/n;
    z = linspace(a, b, m)';
    
    % Matrix containing cardinal polynomials
    P = ones(m, n+1);
    
    for jj = 0:n
        knj = [0:jj-1 jj+1:n];
        
        for kk = knj
            P(:,jj+1) = P(:,jj+1).*(z-xn(kk+1))/(xn(jj+1)-xn(kk+1)); 
        end
    end
end
% If fn is a column vector containing f(x) at the nodes xn, then P*fn 
% provides the interpolant P_n f(z) on the dense grid z