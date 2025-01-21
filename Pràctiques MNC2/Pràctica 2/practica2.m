%% Practical 2 - Cortés & García

clear all;
format long;

%%
% Define the function we are going to be working with. Also, define a mesh
% of points in which we will evaluate our interpolations.

mesh = linspace(-1,1,1000);
f = @(x)  exp(-20*(x+0.25).^2) + 0.25*sin(30*x).*( exp(-20*(x-0.25).^2));

%% (a) The Vandermonde matrix

%%
% As stated, we have $n+1$ nodes in which we will impose 
% $p(x_j)=f(x_j)=f_j$ for $j = \{0, 1, ..., n \}$. Also, we are looking 
% for a polynomial of the form $\Pi_m(x)=a_0 + a_1x + a_2x^2 + ... +
% a_mx^m$. Therefore, we'll have that for any $j$ between $0$ and $n$ (both
% included), we'll have an equation that looks like this:
%
% $$a_0 + a_1x_j + a_2x_j^2 + ... + a_mx_j^m = f_j$$
%
% Therefore, if we write all the equations in a matricial form where the
% vector $a=\{ a_i \}$ are the variables we want to solve the system for. 
% Let $F = \{ f_j \}$ be the vector of the function evaluated in all the 
% $n+1$ nodes and $V = \{ x_j^i \}$ the matrix formed by the powers (up 
% until $m$) of the $n+1$ nodes $x_j$.
% 
% We can see the coefficients that form $V$ look like the ones in the
% Vandermonde Matrix. The only difference is that in our system, the matrix
% will not be squared: it will have $n+1$ rows and $m \leq n$ columns.

%% (b) Equispaced nodes

cases = [14 7; 28 14; 28 20; 64 30];

for ii = 1:4
    n = cases(ii,1);
    m = cases(ii,2);
    
    % Create n equispaced nodes
    jj = 0:n;
    xj = -1+2*jj/n; 
    
    % Create the vandermonde matrix of the equispaced nodes
    V = fliplr(vander(xj));
    V = V(:, 1:m+1);
    
    % Solve the system to find the a_i coefficients of the polynomial interpolation 
    fx = f(xj)';
    UR = myqr(V);
    a = qrsolve(UR,fx);
    
    % Create the Vandermonde matrix for a mesh to evaluate the polynomial
    X = fliplr(vander(mesh));
    X = X(:, 1:m+1);
    
    % Plot f(x) vs the polynomical approximation    
    figure(1)
    subplot(2, 2, ii)
    fplot(f, [-1,1],'k','Linewidth',1)
    title('Least Squares Polynomial Fit', 'Interpreter', 'latex')
    subtitle('With equispaced nodes', 'Interpreter', 'latex')
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$f(x)$', 'Interpreter', 'latex')
    hold on;
    plot(mesh, X*a, 'Linewidth',0.8)    
    legend('$f(x)$',"n = " + num2str(n) + ", m =  " + num2str(m), 'Interpreter', 'latex')
    axis([-1 1 -0.2 1.2])
    hold off;
end 

%% 
% Now we can have a look at how the quality of the interpolation evolves
% when increasing the number of nodes used. We can clearly see that, in the
% central part of the interval, the interpolation gets closer to the actual
% function when $n$ is increased, which is good news. However, it can also
% be observed that when $n$ is takes a larger value, we eventually find the
% effect of Runge Phenomenon in the limits of the interval. This is not so
% good news, since, even when the approximation of the function may be good
% enough in the center, it will not be possible to say the same about the
% extremes.
%
% In MNC1 we already saw that a way of getting rid of this phenomena is
% changing our nodes and so, in the next section, Chebyshev nodes will be
% used.

%% (c) Chebyshev nodes

% We repeat the same procedure for Chevyshev nodes
for ii = 1:4
    n = cases(ii,1);
    m = cases(ii,2);
    
    % Create n Chebyshev nodes
    jj = 0:n;
    xj = cos(jj*pi/n); 
    
    % Create the vandermonde matrix of the Chebyshev nodes
    V = fliplr(vander(xj));
    V = V(:, 1:m+1);
    
    % Solve the system to find the a_i coefficients of the polynomial interpolation 
    fx = f(xj)';
    UR = myqr(V);
    a = qrsolve(UR,fx);
    
    % Create the Vandermonde matrix for a mesh to evaluate the polynomial
    X = fliplr(vander(mesh));
    X = X(:, 1:m+1);
    
    % Plot f(x) vs the polynomical approximation    
    figure(2)
    subplot(2, 2, ii)
    fplot(f, [-1,1],'k','Linewidth',1)
    title('Least Squares Polynomial Fit', 'Interpreter', 'latex')
    subtitle('With Chebyshev nodes', 'Interpreter', 'latex')
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$f(x)$', 'Interpreter', 'latex')
    hold on;
    plot(mesh, X*a, 'Linewidth',0.8)    
    legend('$f(x)$',"n = " + num2str(n) + ", m =  " + num2str(m), 'Interpreter', 'latex')
    hold off;    
end 

%% 
% In this case, we can also see that the interpolation gets better when $n$
% is increased, like when we used the equispaced nodes. Nevertheless, the
% effect of Runge Phenomena in the limits of the interval is not there this
% time. 
% 
% Depending on the number of nodes $n+1$ we are going to be using to 
% interpolate this function, we will probably prefer to choose one node 
% type or the other. For smaller $n$, when the Runge Phenomena doesn't have
% a great effect on the interpolation in the extremes, we will probably
% prefer the equispaced nodes, since they seem to give better results in
% this case. However, if we are going to be using a larger value for $n$,
% the results will be noticeably more accurate if Chebyshev nodes are used.

%% Auxiliar codes

%% 
% The auxiliar codes are commented because if not, Figure 2 did not appear
% when the script was published.

% % Code 15: QR - Householder factorization
% % Input:    A (n x m matrix with n >= m and rank(A) = m)
% % Output:   UR matrix containing (see text):
% %           1) Householder vectors u_j (low triang. part)
% %           2) R (upper triangular matrix)
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
% 
% % Code 16: QR - Solving Ax=b (requires bs.m)
% % Input:    1) UR matrix provided by Code 15, i.e. UR = myqr(A)
% %           2) b right-hand side
% % Output:   x solution minimizing ||Ax-b||
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
% 
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