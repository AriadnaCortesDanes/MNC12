% dibuixa el gradient d'una funcio de dues variables
% Input:
%      - la funcio f
%      - una malla per a x (vector fila)
%      - una malla per a y (vector columna)
% Output:
%      - els valors de la f en la malla 2d
function z = Gradient2d(f,x,y)
    z = x .* exp(-x.^2 - y.^2);
    [px,py] = gradient(z);
    contour(x,y,z)
    hold on
    quiver(x,y,px,py)
    hold off
end