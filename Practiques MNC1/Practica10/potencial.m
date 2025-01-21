% Input: x,y,z (coordenadas) 
%        lamnda, a (parametros anilla)
% Output: V (potencial en el punto)
function V = potencial(x, y, z, lambda, a)
    % Definimos la función del diferencial de potencial
    dV = @(t) (lambda*a)./sqrt((x-a*cos(t)).^2 + (y-a*sin(t)).^2 + z.^2);
    
    V = -1*qclencurt(0, 2*pi, 20, dV);
end