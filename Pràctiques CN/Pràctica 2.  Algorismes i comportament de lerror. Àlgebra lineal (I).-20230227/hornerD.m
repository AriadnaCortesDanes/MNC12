function [ pa ] = hornerD(coeff,a,d)
%   HORNER avalua el polinomi de grau n
%   p(x) = coeff(1)x.^n + coeff(2)x.^(n-1) + ... + coeff(n-1)x + coeff(n)
%   en x=a pel mètode de horner fent ús de d dígits
% Falta implementar :) 
pa = coeff(1);
x=a;
for k=2:length(coeff)
    pa = pa.*x + coeff(k);
end
return
end