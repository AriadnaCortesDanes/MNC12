function [x, GSB] = gaussSeidel(A, b, tol)
    D = diag(diag(A));
    L = tril(A,-1);
    U = triu(A,1);
    aux = inv(L+D);
    GSB = -aux*U;
    GSu = aux*b;
    rhoGS = abs(eigs(GSB,1));
    
    if rhoGS > 1
        fprintf('Gauss Seidel no converge');
    end
    nb = norm(b,"inf");
    x = zeros(size(b));
    x_prev = zeros(size(b));
    test = 1;

    while test > tol
        x = GSB*x_prev+GSu;
        r = norm(b-A*x,"inf");
        test = r/nb;
        x_prev = x;
    end

    fprintf("GS (n=%d) r-inf: %d", length(b), r(end));
end