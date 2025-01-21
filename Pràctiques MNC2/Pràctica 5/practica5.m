%% Practical 5 - Cortés & García

clear all;
close all;
format long;

%%

p = @(v, T) 8*T./(3.*v - 1) - 3./(v.^2); % Van der Waals gas equation

% We'll use vk = [vl_k, vg_k] for iteration k
% Initialize starting vector v0 for T = 0.99
v0 = [.8; 1.2];

vdw_fun = @(v, T) [log((3*v(2)-1)/(3*v(1)-1)) + (9/(4*T))*(1/v(2)-1/v(1)) - 1/(3*v(2)-1) + 1/(3*v(1)-1);
    8*T/3*(1/(3*v(2)-1) - 1/(3*v(1)-1)) - 1/(v(2)^2) + 1/(v(1)^2)];

% Initialize vectors to keep track of the solutions of vl and vg and their
% respective solutions pl and pg.
vVec = [];
pVec = [];

% Iterate for {0.99, 0.98, ..., 0.85}
for T = .99:-.01:.85
    new_vdw = @(v) vdw_fun(v, T);

    [Vk, resd, it] = newtonn(v0, eps, 50, new_vdw);

    vk = Vk(:, end);
    v0 = vk;

    pl = p(vk(1), T);
    % We only compute pl because of the condition that pl = pg

    vVec = [vVec, vk(1), vk(2)];
    pVec = [pVec, pl, pl];

end


v0 = [.8; 1.2]; % Reinitialize v0

% Iterate for {0.9900, 0.9903, ..., 0.9999}
% We want to get closer to T=T_critical but we can't really get there
for T = .99:.0003:.9999
    new_vdw = @(v) vdw_fun(v, T);

    [Vk, resd, it] = newtonn(v0, eps, 50, new_vdw);

    vk = Vk(:, end);
    v0 = vk;

    pl = p(vk(1), T);
    % We only compute pl because of the condition that pl = pg

    vVec = [vVec, vk(1), vk(2)];
    pVec = [pVec, pl, pl];
end

% Sort vVec and pVec to display them properly
[vVec, I] = sort(vVec);
pVec = pVec(I);

vol = linspace(0, 5, 500); % Define a mesh for v

% Plot several isotherm curves
figure(1)
for T = linspace(0.85, 1.25, 9)
    plot(vol, p(vol, T),'LineWidth', .7)
    hold on;
end

plot(vVec, pVec, 'k-', 'LineWidth', 1.3);
title("Van der Waals gas equation","Interpreter",'latex')
xlabel("$v$","Interpreter",'latex')
ylabel("$p$","Interpreter",'latex')
axis([0, 5, 0, 1.5])
hold off;
drawnow;

%%

%% Auxiliar codes

% % Code 20: Newton’s method for n-dimensional systems
% % Input: x0 - initial guess (column vector)
% %        tol - tolerance so that ||x_{k+1} - x_{k} || < tol
% %        itmax - max number of iterations
% %        fun - function’s name
% % Output:  XK - iterated
% %  resd: resulting residuals of iteration: ||F_k||
% %  it:   number of required iterations to satisfy tolerance
% function [XK,resd,it] = newtonn(x0,tol,itmax,fun)
%     xk = [x0];
%     resd = [norm(feval(fun,xk))];
%     XK = [x0]; 
%     it = 1;
% 
%     tolk = 1.0; 
%     n = length(x0);
% 
%     while it < itmax && tolk > tol
%         Fk = feval(fun, xk);
% 
%         DFk = jac(fun, xk); 
%         [P,L,U] = pplu(DFk);
% 
%         dxk = plusolve(L,U,P,-Fk);
% 
%         xk = xk + dxk; 
%         XK = [XK xk]; 
%         resd = [resd norm(Fk)];
%         tolk = norm(XK(:, end)-XK(:, end-1)); 
%         it = it + 1;
%     end
% end

% % Code 19: Computation of the Jacobian J
% % Input:   F(x) : R^m ---> R^n
% %          x : (m x 1)-vector ; F: (n x 1)-vector
% % Output: DF(x) (n x m) Jacobian matrix at x
% function DF = jac(F,x)
%     f1 = feval(F,x); 
%     n = length(f1); 
%     m = length(x);
% 
%     DF = zeros(n,m); 
%     H = sqrt(eps)*eye(m);
% 
%     for j = 1:m
%         f2 = feval(F,x+H(:,j)); 
%         DF(:,j) = (f2 - f1)/H(j,j);
%     end
% end

% % Code 13: PA = LU factorization (partial pivoting)
% % Input: A (non-singular square matrix)
% % Output: L (unit lower triangular matrix)
% %         U (upper triangular matrix)
% %         P (reordering vector)
% function [P, L, U] = pplu(A)
%     [m,n] = size(A);
% 
%     if m~=n
%            error('not square matrix'); 
%     end
% 
%     U = A;
%     L = eye(n);
% 
%     P = [1:n]';
% 
%     for k = 1:n-1
%         [~, imax] = max(abs(U(k:end,k)));
%         imax = imax+k-1;
%         i1 = [k, imax];
%         i2 = [imax, k];
% 
%         U(i1,:) = U(i2,:); % Column k will be column imax and column imax will be column k
%         P(k) = imax;
% 
%         L(i1,1:k-1) = L(i2, 1:k-1);
% 
%         for jj = [k+1:n]
%             L(jj, k) = U(jj, k)/U(k, k);
%             U(jj, k:n) = U(jj, k:n) - L(jj, k)*U(k,k:n);
%         end
%     end
% end

% % Code 14: PA = LU (Solver for Ax = b)
% % Input:    L (unit lower triangular matrix)
% %           U (upper triangular matrix)
% %           P (reordering vector)
% %           b (right-hand side)
% % Output:   solution x
% function x = plusolve(L, U, P, b)
%     n = length(b);
%     for k = 1:n-1
%         b([k P(k)]) = b([P(k) k]);
%     end
%     y = fs(L, b);
%     x = bs(U, y);
% end

% % Code 11: Forward Substitution for Lower Triangular Systems
% % Input:    L: Low Triangular non-singular square matrix
% %           b: column right-hand side
% % Output:   x: solution of Lx=b
% function x = fs(L, b)
%     x = 0*b;
%     n = length(b);
%     x(1) = b(1)/L(1,1);
% 
%     for ii = 2:n
%         x(ii) = (b(ii)-L(ii, 1:ii-1)*x(1:ii-1))/L(ii,ii);
%     end
% end

% % Code 12: Backward Substitution for Upper Triangular Systems
% % Input:    U: Upp. Triangular non-singular square matrix
% %           b: column right-hand side
% % Output:   x: solution of Ux=b
% function x = bs(U, b)
%     x = 0*b;
%     n = length(b);
%     x(n) = b(n)/U(n,n);
% 
%     for ii = n-1:-1:1
%         x(ii) = (b(ii)-U(ii, ii+1:n)*x(ii+1:n))/U(ii,ii);
%     end
% end