function  [w,t]  = EulerMod( f,a,b,h,alpha )
    % Mètode d'Euler Modificat o punto medio
    t=[a:h:b];
    N=length(t);
    w(1)=alpha;
    for i=1:N-1
        k1 = w(i)+0.5*h*f(t(i),w(i) );
        w(i+1)=w(i)+h*f(t(i)+ h/2,k1);
    end
end