%Code 9: Cot-Map for 1st kind improper integrals [0,+inf)
% Input: n: #abscissas; L: scaling factor
% Output: I_n(f)
function In = qcot(n,L,fun)
    m = [1:n];
    t = m'*pi/(n+1);
    x = L*(cot(t/2)).^ 2;
    w = 0*x.';
    f = feval(fun,x);
    for ii=1:n
        c1 = (4*L)*sin(t(ii))/((n+1)*(1-cos(t(ii)))^2);
        w(ii) = c1*sin(m*t(ii))*((1-cos(m'*pi))./m');
    end
    In = w*f;
end