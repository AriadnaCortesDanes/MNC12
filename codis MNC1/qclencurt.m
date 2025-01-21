% Input: a-b (low-up lim.) 
%        n (# nodes-1) 
%        fun (func. name) 
% Output: I_n(f)
function Icc = qclencurt(a, b, n, fun)
    l = [0:n]'; 
    k = [2:n]'; 
    x = cos(l*pi/n);
    
    w = cos(l*l'*pi/n) \ [2;0;(1+(-1).^k)./(1-k.^2)];
    
    z = a+.5*(b-a)*(x+1); 
    f = feval(fun,z);
    
    Icc =.5*(b-a)*w'*f;
end