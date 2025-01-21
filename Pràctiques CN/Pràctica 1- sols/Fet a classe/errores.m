function [ea,er,dc,xs] = errores(x,y)
% exemple de funció
%   x valor exacte
%   y valor aproximat
ea = abs(x-y);
er = ea./abs(x);
erp = er*100;
dc = fix(-log10(2*ea));
xs = fix(-log10(2*er));
end