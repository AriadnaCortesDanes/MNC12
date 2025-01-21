%% Practical 4 - Cortés & García

clear all;
close all;
format long;

%% (a) Kirchoff

% In Practical 1 (5.1) we built the matrix A as:
n = 30;
r = 3; v = 1;

A = zeros(2*n-1,2*n-1);
A(1,1:2) = [2*r r];
A(end, end-1:end) = [-r 3*r];

for ii = 2:2*n-2
    if (mod(ii,2) == 0)
        A(ii, ii-1:ii+1) = [1 -1 -1];
    else
        A(ii, ii-1:ii+1) = [-r 2*r r];
    end
end

% Also, b was built as:
b = zeros(2*n-1,1);
b(1) = v;

% Ax coincides in both cases
x = rand(2*n-1, 1);
fprintf("The norm of A*x - Afun_pract14(x) is %d \n", norm(A*x - Afun_pract14(x)))

%%
% We can clearly see that, as expected, the given function applied to an x
% vector returns the same result as the matrix we built in Practical 1.

%% (b) GMRES algorithm

tolVec = 10.^(-1*[2:8]);
for m_max = 10:5:30
    errVec = [];
    kVec = [];
    for tol = tolVec
        [x, k] = mygmres(@Afun_pract14, b, tol, m_max);

        errVec = [errVec ; norm(Afun_pract14(x) - b)];
        kVec = [kVec; k];
    end

    figure(1)
    semilogx(tolVec, kVec, "-o","LineWidth", 1);
    hold on;
    
    figure(2)
    loglog(tolVec, errVec, "-o", "LineWidth", 1);
    hold on;
end

figure(1)
title("Iterations needed as a function of asked tolerance", "Interpreter","latex")
xlabel("$\epsilon$", "Interpreter","latex")
ylabel("$m(\epsilon)$", "Interpreter","latex")
legend("$m = 10$", "$m = 15$", "$m = 20$","$m = 25$","$m = 30$", "Interpreter","latex")
drawnow;

figure(2)
title("Error as a function of asked tolerance", "Interpreter","latex")
xlabel("$\epsilon$", "Interpreter","latex")
ylabel("$||Ax-b||$", "Interpreter","latex")
legend("$m = 10$", "$m = 15$", "$m = 20$","$m = 25$","$m = 30$", "Interpreter","latex", "Location","southeast")
drawnow;


hold off;

%%
% In Figure 1 we can see that when both the error of the approximation x
% and the tolerance limit set are in a logarithmic scale, each 
% value of m from 10 to 30 produces different resuts. m in this figure is
% the dimension of the Krylov space we have reached for every tolerance
% value. Hence, it was expected to see that for larger tolerances, we'd
% need a lower value of m (given the way in which we are constructing the
% Krylov Space). This is because, in the end, the tolerance is a threshold
% to whether we consider two vectors orthogonals or not.
%
% In Figure 2, on the other hand, we can see how, when we ask for a smaller
% tolerance, the error of the approximation with respect to the actual 
% result also reduces. We see, though that there is a moment for each value
% of m (maximum dimension of the Krylov Space) where the error cannot be
% reduced anymore and we can see it remains constant even when we ask for a
% smaller tolerance.

%% Auxiliar codes

%%

% function Ax = Afun_pract14(x)
%     R = 3; V = 1; n = 30;
%     Ax = [];
% 
%     Ax = [Ax; 2*R*x(1) + R*x(2)];
%     Ax = [Ax; x(1)-x(2)-x(3)];
% 
%     for j = 3:n
%         Ax = [Ax; 2*R*x(2*j-3) + R*x(2*j-2) - R*x(2*j-4)];
%         Ax = [Ax; x(2*j-3) - x(2*j-2) - x(2*j-1)];
%     end
%
%     Ax = [Ax; 3*R*x(2*n-1) - R*x(2*n-2)];
% end

%%

% Code 18: GMRES (Arnoldi iteration - Reorthogonalization)
% function [x, k] = mygmres(Afun, b, tol, dimkryl)
% Requires: mytqr.m and qrsolve.m (Codes 15 & 16)
% Input:    Afun (x --> Ax function) | b (r.h.s.)
%           tol (iteration stops if |A*x_m-b| < tol)
%           dimkryl: max. dimension of krylov space.
% Output:   x: approximation found
%           k: no. of iterations used
%     hkp1 = 1;
%     resd = 1;
%     k = 1;
% 
%     Q(:, 1) = b/norm(b);                        % q_1 = b/|b|
%     H = [];
% 
%     while resd > tol && hkp1 > eps && k < dimkryl
%         Q(:, k+1) = feval(Afun, Q(:, k));       % q_{k+1} = Aq_k
%         h = Q(:, 1:k)'*Q(:, k+1);               % h_{ij} - coeffs.
%         Q(:, k+1) = Q(:, k+1) - Q(:, 1:k)*h;    % See (5.225)
%         S = Q(:, 1:k)'*Q(:, k+1);
%         Q(:, k+1) = Q(:, k+1) - Q(:, 1:k)*S;
%         h = h + S;                              % Reorth.
%         hkp1 = norm(Q(:, k+1));
%         H(1:k+1, k) = [h; hkp1];
%         Q(:, k+1) = Q(:, k+1)/hkp1;
% 
%         UR = myqr(H(1:k+1, 1:k));
%         e1b = [norm(b);  zeros(k, 1)];
% 
%         y = qrsolve(UR, e1b);
%         x = Q(:, 1:k)*y;
% 
%         k = k+1;
%         resd = norm(feval(Afun, x)-b);
%     end
% end

%%

% Code 15: QR - Householder factorization
% Input:    A (n x m matrix with n >= m and rank(A) = m)
% Output:   UR matrix containing (see text):
%           1) Householder vectors u_j (low triang. part)
%           2) R (upper triangular matrix)
% function UR = myqr(A)
%     [n, m] = size(A);
% 
%     if n > m
%         M = m;
%     else
%         M = m-1;
%     end
% 
%     for jj = 1:M
%         B = A(jj:n,jj:m);
%         x = B(:, 1);
% 
%         tau = sign(x(1))*norm(x);
%         gam = 1 + x(1)/tau;
% 
%         u = [1 ; x(2:end)/(tau + x(1))];
%         vauxT = gam*u'*B;
% 
%         A(jj:n,jj:m) = B - u*vauxT;
%         A(jj+1:n, jj) = u(2:end);
%     end
% 
%     UR = A;
% end

%%

% Code 16: QR - Solving Ax=b (requires bs.m)
% Input:    1) UR matrix provided by Code 15, i.e. UR = myqr(A)
%           2) b right-hand side
% Output:   x solution minimizing ||Ax-b||
% function x = qrsolve(UR, b)
%     [n, m] = size(UR);
% 
%     if n > m
%         M = m;
%     else
%         M = m-1;
%     end
% 
%     for jj = 1:M
%         u = [1; UR(jj+1:n, jj)];
%         gam = 2/norm(u)^2;
% 
%         b(jj:n) = b(jj:n) - gam*u*(u'*b(jj:n));
%     end
% 
%     x = bs(triu(UR(1:m,1:m)), b(1:m));
% end

%%

% Code 12: Backward Substitution for Upper Triangular Systems
% Input:    U: Upp. Triangular non-singular square matrix
%           b: column right-hand side
% Output:   x: solution of Ux=b
% function x = bs(U, b)
%     x = 0*b;
%     n = length(b);
%     x(n) = b(n)/U(n,n);
% 
%     for ii = n-1:-1:1
%         x(ii) = (b(ii)-U(ii, ii+1:n)*x(ii+1:n))/U(ii,ii);
%     end
% end