%% Practical 1 - Cortés & García

format long; 
clear all;

%% Solving an electric circuit using LU-factorization
%
% Let us define the variables given in the statement of the problem.
r = 3;  % Resistance value
v = 1;  % Potential value

%%
% We can deduce that the matrix to solve the system is the given one by
% doing an analysis on the meshes and nodes of the circuit:
% 
% * Mesh 1: $-V + I_1*R + I_2*R + I_1*R=0 \rightarrow 2*I_1*R + I_2*R = V$
%
% * Mesh k (for $1<k<n$): $-I_{2k-2}*R + 2*I_{2k-1}*R + I_{2k}*R = 0$
%
% * Mesh n: $-I_{2n-2}*R + 3*I_{2n-1}*R = 0$
%
% * Upper-right node of the mesh k (for $1\leq k < n$): $I_{2k-1} =
% I_{2k}+I_{2k+1}$
% 
% This equations end up creating the matricial system given in the
% statement.

%% (a) Intensities

for n = 5:5:30
    A = zeros(2*n-1,2*n-1);
    A(1,1:2) = [2*r r];
    A(end, end-1:end) = [-r 3*r];

    % Build the matrix 
    for ii = 2:2*n-2
        if (mod(ii,2) == 0)
            A(ii, ii-1:ii+1) = [1 -1 -1];
        else
            A(ii, ii-1:ii+1) = [-r 2*r r];
        end
    end
    
    % Build the independent term b
    b = zeros(2*n-1,1);
    b(1) = v;
    
    % LU factorization with partial pivoting
    [P, L, U] = pplu(A);
    
    % Solve the system
    I = plusolve(L, U, P, b);
    
    % For n = 10 we show the results
    if(n == 10)
        disp(['The intensity I1 is ' num2str(I(1))]);
        disp(['The intensity I2 is ' num2str(I(2))]);
        disp(['The intensity I15 is ' num2str(I(15))]);
        
        figure(1) % Plot the intensity decay for n=10
        semilogy([1:2*n-1], I, "k-o")
        title('Intensity decay','interpreter', 'latex')
        xlabel('k','interpreter', 'latex')
        ylabel('$log(I_k)$','interpreter', 'latex')
    end
end

%%
% Figure 1 has been plotted with the intensities $I_k$ in a logarithmic
% scale. Besides, it can be seen that $log(I_k)$ seems to decay linearly 
% with $k$. This means that $log(I_k)=\alpha + \beta k$, which implies that
% $I_k = 10^{\alpha + \beta k} = ab^k$. So, it can be said that the 
% intensities have an exponential decay with the "number" of the mesh $k$.

%% (b) Equivalence resistance
R_eq = [];  % Vector storing equivalent resistances

for n = 2:1:16
    A = zeros(2*n-1,2*n-1);
    A(1,1:2) = [2*r r];
    A(end,end-1:end) = [-r 3*r];

    % Build the matrix 
    for ii = 2:2*n-2
        if (mod(ii,2) == 0)
            A(ii, ii-1:ii+1) = [1 -1 -1];
        else
            A(ii, ii-1:ii+1) = [-r 2*r r];
        end
    end
    
    % Build the independent term b
    b = zeros(2*n-1,1);
    b(1) = v;
    
    % LU factorization with partial pivoting
    [P, L, U] = pplu(A);
    
    % Solve the system
    I = plusolve(L, U, P, b);
    
    % Add new equivalent resistance to the vector
    R_eq = [R_eq v/I(1)];
   
end

%%
% If we have a look at the circuit we are dealing with, we can see that,
% for $n$ meshes, we can rewrite the equivalent resistance as:
%
% $R_e(n) = 2R + \frac{1}{\frac{1}{R} + \frac{1}{R_e(n-1)}}$
%
% Now, when $n\to \infty$, $R_e(n)=R_e(n-1)=R_e$. If we work out the 
% expression using this limit, we eventually get that $R_e^2-2RR_e-2R^2$, which
% has two possible solutions: $R_e=(1+\sqrt{3})R$ and $R_e=(1-\sqrt{3})R$.
% The second one gives a negative result, which has no phyisical meaning,
% and thus we will consider the first one the correct result, which is the
% value given in the statement.

figure(2)
plot(2:1:16, abs(R_eq - r*(1+sqrt(3))), "k-o")
title('Distance $|R_{eq}-R(1-\sqrt{3})|$ as a function of n','interpreter', 'latex')
xlabel('n','interpreter', 'latex')
ylabel('$|R_{eq}-R(1-\sqrt{3})|$','interpreter', 'latex')

figure(3)
semilogy(2:1:16, abs(R_eq - r*(1+sqrt(3))), "k-o")
title('Distance $|R_{eq}-R(1-\sqrt{3})|$ as a function of n','interpreter', 'latex')
xlabel('n','interpreter', 'latex')
ylabel('$log |R_{eq}-R(1-\sqrt{3})|$','interpreter', 'latex')

%%
% When we plot the absolut difference between the equivalent resistance 
% and $R(1+\sqrt{3})$ in Figure 2, we observe that the plot tends to 0, 
% confirming that the limit of the equivalent resistance tends to 
% $R(1+\sqrt{3})$.
%
% Besides, to get a better look at it, we have also plot it with the 
% y-axis in a logarithmic scale (Figure 3), we see that the last 
% iteration ($n=15$) does not appear in this case. We will then consider 
% the limit es reached with $n=15$, since we can see that MatLab will 
% not consider the logarithm of this distance as an actual number (as if 
% it were computing $log(0)$).
%
% Since the limit is considered to be reached with $n=15$, where MatLab 
% consideres that the distance to the mathematical limit equals 0, we can 
% conclude that the limit is surpassed in part (a), since we solved for 
% circuits with $n$ greater than 15.

%% Auxiliar codes

% Code 11: Forward Substitution for Lower Triangular Systems
% Input:    L: Low Triangular non-singular square matrix
%           b: column right-hand side
% Output:   x: solution of Lx=b
function x = fs(L, b)
    x = 0*b;
    n = length(b);
    x(1) = b(1)/L(1,1);

    for ii = 2:n
        x(ii) = (b(ii)-L(ii, 1:ii-1)*x(1:ii-1))/L(ii,ii);
    end
end

% Code 12: Backward Substitution for Upper Triangular Systems
% Input:    U: Upp. Triangular non-singular square matrix
%           b: column right-hand side
% Output:   x: solution of Ux=b
function x = bs(U, b)
    x = 0*b;
    n = length(b);
    x(n) = b(n)/U(n,n);

    for ii = n-1:-1:1
        x(ii) = (b(ii)-U(ii, ii+1:n)*x(ii+1:n))/U(ii,ii);
    end
end

% Code 13: PA = LU factorization (partial pivoting)
% Input: A (non-singular square matrix)
% Output: L (unit lower triangular matrix)
%         U (upper triangular matrix)
%         P (reordering vector)
function [P, L, U] = pplu(A)
    [m,n] = size(A);

    if m~=n
           error('not square matrix'); 
    end

    U = A;
    L = eye(n);

    P = [1:n]';

    for k = 1:n-1
        [~, imax] = max(abs(U(k:end,k)));
        imax = imax+k-1;
        i1 = [k, imax];
        i2 = [imax, k];
        
        U(i1,:) = U(i2,:); % Column k will be column imax and column imax will be column k
        P(k) = imax;

        L(i1,1:k-1) = L(i2, 1:k-1);

        for jj = [k+1:n]
            L(jj, k) = U(jj, k)/U(k, k);
            U(jj, k:n) = U(jj, k:n) - L(jj, k)*U(k,k:n);
        end
    end
end

% Code 14: PA = LU (Solver for Ax = b)
% Input:    L (unit lower triangular matrix)
%           U (upper triangular matrix)
%           P (reordering vector)
%           b (right-hand side)
% Output:   solution x
function x = plusolve(L, U, P, b)
    n = length(b);
    for k = 1:n-1
        b([k P(k)]) = b([P(k) k]);
    end
    y = fs(L, b);
    x = bs(U, y);
end