% Code 17: GSR (Gram--Schmidt-Reorthogonalized)
% Input:    1) A (n x m matrix -- n >= m)
% Output:   2) Q (n x m isometric matrix)
%              R (m x m upper triangular matrix)
function [Q, R] = gsr(A)
    [~, m] = size(A);
    Q = 0*A;
    R = zeros(m);

    R(1, 1) = norm(A(:, 1));
    Q(:, 1) = A(:, 1)/R(1, 1);

    for k = 2:m
        % Orthogonalization
        R(1:k-1, k) = Q(:, 1:k-1)'*A(:, k);
        Q(:, k) = A(:, k) - Q(:, 1:k-1)*R(1:k-1, k);
        S = Q(:, 1:k-1)'*Q(:, k);
        
        % Reorthonormalization
        Q(:, k) = Q(:, k) - Q(:, 1:k-1)*S;
        R(k, k) = norm(Q(:, k));

        if R(k, k) > 1e-16
            Q(:, k) = Q(:, k)/R(k, k);
        else
            ['Lin. Dep.']
        end
    end
end