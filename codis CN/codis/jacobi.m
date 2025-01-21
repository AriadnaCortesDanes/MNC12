function [x, JB] = jacobi(A, b, tol)
    D = diag(diag(A));
    aux = diag(1./diag(A));
    JB = -aux*(A-D);
    Jc = aux*b;
    rhoJ = abs(eigs(JB,1));
    if rhoJ > 1
        fprintf('Jacobi no converge');
    end
    y = [1 2 -1 1].';
    k = 1;
    nb = norm(b,"inf");
    x = zeros(size(b));
    r(k) = norm(b-A*x,"inf");
    test = 1;
    while test > tol
        x = JB*x+Jc;
        r(k+1) = norm(b-A*x,"inf");
        test = r(k+1)/nb;
        k = k+1;
    end

    fprintf("JB (n=%d) r-inf: %d", length(b), r(end));
end