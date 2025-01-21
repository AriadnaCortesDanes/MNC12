% A = LU factorization (no pivoting)
function [L, U] = milu(A)
    [m,n]=size(A);

    if m~=n
           error('not square matrix'); 
    end

    U = A;
    L = eye(n);
    %['***************************   k = ' int2str(k) ' ***************]
    for ii = [k+1:n]
        L(ii, k) = U(ii, k)/U(k, k); % ii-th Multiplier m_ii,k

        U(ii, k:n) = U(ii, k:n) - L(ii, k)*U(k,k:n); % Update rows

    end
    % L(k+1:n, k)
    % U

