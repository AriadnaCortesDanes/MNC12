%---------------------------------------
function [t, w ] = EulerV( f,a,b,h,alpha )
% Mètode d'Euler
% f una funcion @algo que retorna un vector columna [y1; y2; ...]
t=[a:h:b];
N=length(t);
w(1,:)=alpha;
    for i=1:N-1
        w(i+1,:)=w(i,:)+h*f(t(i),w(i,:))';
    end
end
%%%%