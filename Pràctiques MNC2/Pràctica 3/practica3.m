%% Practical 3 - Cortés & García

clear all;
format long;

%%
nn = 2:16;

%% (a) CGS

% Vector d'errors que guardara ||I_n - QtQ||
err_cgs = [];
for n = nn
    H = hilb(n);
    
    % Si n = 8, monotoritzem les normes dels vectors
    if n == 8
        [Q, R, S] = mycgs(H);
        
        fprintf("For n=%d \n", n)
        for kk = 1:n
            for jj = 1:kk-1
                fprintf("k = %d,    j=%d,   q_j*q_k = %d \n \n", kk, jj, S(jj, kk))
            end
        end

        figure(1)
        semilogy(1:n, diag(R), 'k', 'LineWidth', 1);
        title('$\tilde{q_k}$ norms', 'Interpreter', 'latex')
        xlabel('$\tilde{q_k}$','Interpreter', 'latex')
        ylabel('norm($\tilde{q_k}$)','Interpreter', 'latex')

    else
        [Q, ~, ~] = mycgs(H);
    end

    err_cgs = [err_cgs norm(eye(n)-Q'*Q)];
end

%%
% We can see that scalar products should be equally 0 (since the vectors
% are orthogonal, but they are not (they reach small values but not even
% close to 0). This will improve when we do the reorthonormalization. Also,
% we observe that the norms of the q_k vectors decrease exponentially as
% expected.

%% (b) GSR

err_gsr = [];
for n = nn
    H = hilb(n);
    
    if n == 8
        [Q, ~, S1, S2] = gsr(H);
        
        fprintf("For n=%d \n", n)
        for kk = 1:n
            for jj = 1:kk
                fprintf("k = %d,    j=%d \n", kk, jj)
                fprintf("After ortogonalization: q_j*q_k = %d \n", S1(jj, kk))
                fprintf("After reortogonalization: q_j*q_k = %d \n \n", S2(jj, kk))
            end
        end

    else
        [Q, ~, ~, ~] = gsr(H);
    end

    err_gsr = [err_gsr norm(eye(n)-Q'*Q)];
end

%%
% In the output of this section we can observe the values of the scalar
% products after the first orthogonalization and after the second
% reorthogonalization. We can see that first one reaches a small value but
% not zero. The second scalar product is nearly zero. This shows us the
% improvement when doing a second orthogonalization. 

%% (c) QR

err_qr = [];
for n = nn
    H = hilb(n);

    [Q, ~] = qr(H);

    err_qr = [err_qr norm(eye(n)-Q'*Q)];
end


figure(2)
semilogy(nn, err_cgs, 'LineWidth', 1)
hold on;
semilogy(nn, err_gsr, 'LineWidth', 1)
semilogy(nn, err_qr, 'LineWidth', 1)
hold off;
legend('CGS', 'GSR', 'QR',  'Interpreter', 'latex')
title('Error commited $\|I_n - Q^T*Q\|$', 'Interpreter', 'latex')
xlabel('n','Interpreter', 'latex')
ylabel('$\|I_n - Q^T*Q\|$','Interpreter', 'latex')

%%
% In this figure we represent the norm of $\|I_n - Q^T*Q\|$. This error
% increases as we increase the size of the matrix. 

%% Auxiliar codes

% Classical Gram--Schmidt (CGS) unstable algorithm
function [Q, R, S] = mycgs(A)
    [~, m] = size(A);
    Q = 0*A;
    R = zeros(m);
    S = R;

    R(1, 1) = norm(A(:, 1));
    Q(:, 1) = A(:, 1)/R(1, 1);

    for k = 2:m
        R(1:k-1, k) = Q(:, 1:k-1)'*A(:, k);
        
        Q(:, k) = A(:, k) - Q(:, 1:k-1)*R(1:k-1, k);
        
        S(1:k-1, k) = Q(:, 1:k-1)'*Q(:, k);
        
        R(k, k) = norm(Q(:, k));

        Q(:, k) = Q(:, k)/R(k, k);
    end
end

% Code 17: GSR (Gram--Schmidt-Reorthogonalized)
% Input:    1) A (n x m matrix -- n >= m)
% Output:   2) Q (n x m isometric matrix)
%              R (m x m upper triangular matrix)
function [Q, R, T1, T2] = gsr(A)
    [~, m] = size(A);
    Q = 0*A;
    R = zeros(m);
    T1 = R;
    T2 = R;

    R(1, 1) = norm(A(:, 1));
    Q(:, 1) = A(:, 1)/R(1, 1);

    for k = 2:m
        % Orthogonalization
        R(1:k-1, k) = Q(:, 1:k-1)'*A(:, k);
        Q(:, k) = A(:, k) - Q(:, 1:k-1)*R(1:k-1, k);
        S = Q(:, 1:k-1)'*Q(:, k);

        T1(1:k-1, k) = Q(:, 1:k-1)'*Q(:, k);
        
        % Reorthonormalization
        Q(:, k) = Q(:, k) - Q(:, 1:k-1)*S;
        R(k, k) = norm(Q(:, k));

        T2(1:k-1, k) = Q(:, 1:k-1)'*Q(:, k);

        if R(k, k) > 1e-18
            Q(:, k) = Q(:, k)/R(k, k);
        else
            ['Lin. Dep.']
        end
    end
end