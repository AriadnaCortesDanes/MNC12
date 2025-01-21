% Code 5B: Chebyshev Differentiation matrix
% Input: n 
% Output: differentiation matrix D and Chebyshev nodes
function [D,x] = chebdiff(n)
x =cos([0:n]'*pi/n); d = [.5 ;ones(n-1,1);.5];
D = zeros(n+1,n+1); 
for ii = 0:n
 for jj = 0:n
  ir = ii + 1 ; jc = jj + 1;
  if ii == jj
   kk = [0:ii-1 ii+1:n]'; num = (-1).^kk.*d(kk+1) ;
   D(ir,jc) =((-1)^(ir)/d(ir))*sum(num./(x(ir)-x(kk+1)));
  else
   D(ir,jc) = d(jc)*(-1)^(ii+jj)/((x(ir)-x(jc))*d(ir));
  end
 end
end
    

   