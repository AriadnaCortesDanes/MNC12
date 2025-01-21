% Code 4: Chebyshev Barycentric Interpolation
% Chebyshev nodes of second kind (practical abscissas)
function [ff,xn,z] = barifun(a,b,n,m,f)
% Input
    % a, b: limits inferior i superior del interval    
    % n = nombre de punts equiespaiats (x)
    % m = nombre de punts de la malla (z)
    % f = funcio a interpolar
% Output
    % ff: interpolant evaluat a la malla z
    % xn: nodes de chevychev
    % z: malla
    xn = ( a + (b-a)*(1+cos([0:n]*pi/n))/2)';
    z = linspace(a,b,m)';
    fn = f(xn);
    n = length(xn); ff = 0*z; num = ff ; den = ff ;
    s = (-1).^[0:n-1]';
    sn = s.*[fn(1)/2; fn(2:end-1); fn(end)/2];
    sd = [s(1)/2; s(2:end-1); s(end)/2];
    for ii = 1:length(num);
        num = ((z(ii) - xn).^(-1))'*sn;
        den = ((z(ii) - xn).^(-1))'*sd;
        ff(ii) = num/den;
    end
end