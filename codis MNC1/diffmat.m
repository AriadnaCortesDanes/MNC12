% Code 5A: Differentiation matrix
function D = diffmat(x)
    % Input: vector of nodes x= [x_0; x_1; . . .; x_n]
    % Output: differentiation matrix D
    N = length(x); lamb = zeros(N,1); D = zeros(N,N);
    for jj = 1:N
        lamb(jj) = 1/prod(x(jj) - x([1:jj-1 jj+1:N]));
    end
    for ii = 1:N
        for jj = 1:N
            if jj == ii
                D(ii,jj) = sum(1./(x(ii)-x([1:ii-1 ii+1:N])));
            else
                D(ii,jj) = lamb(jj)/((x(ii)-x(jj))*lamb(ii));
            end
        end
    end
end