% Code 16: QR - Solving Ax=b (requires bs.m)
% Input:    1) UR matrix provided by Code 15, i.e. UR = myqr(A)
%           2) b right-hand side
% Output:   x solution minimizing ||Ax-b||
function x = qrsolve(UR, b)
    [n, m] = size(UR);

    if n > m
        M = m;
    else
        M = m-1;
    end

    for jj = 1:M
        u = [1; UR(jj+1:n, jj)];
        gam = 2/norm(u)^2;

        b(jj:n) = b(jj:n) - gam*u*(u'*b(jj:n));
    end

    x = bs(triu(UR(1:m,1:m)), b(1:m));
end