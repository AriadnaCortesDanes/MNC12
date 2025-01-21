% Code 10: tanh-rule for 2nd kind improper integrals (-1, +1)
% Input:    n (3 abcisses)
%           a-b (integration domain)
%           c (tanh scaling factor)
% Output:   I_n(f)
function In = qtanh(n, a, b, c, fun)
    h = c/sqrt(n);
    u = [-n:n]*h/2;
    
    x = (b-a)*.5*(tanh(u) + 1) + a;
    f = feval(fun, x);
    
    In = (b-a)*(.25*h)*sum(f./(cosh(u).^2));
end