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