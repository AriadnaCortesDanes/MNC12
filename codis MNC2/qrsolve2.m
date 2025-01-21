% Code 16: QR - Solving Ax=b (requires bs.m)
% Input:    1) Q, R matrix provided by Code 15, i.e. QR = mycgs(A)
%           2) b right-hand side
% Output:   x solution minimizing ||Ax-b||
function x = qrsolve2(Q, R, b)
    % QRx = b => Rb = Q'b since Q is orthogonal
    right_side = Q'*b;
    x = bs(R, right_side);
end