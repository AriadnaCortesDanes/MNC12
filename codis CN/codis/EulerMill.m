function  [w,t]  = EulerMill( f,a,b,h,alpha )
% Mètode d'Euler Millorat o HEUN
    t=[a:h:b];
    N=length(t);
    w(1)=alpha;
    for i=1:N-1
        k1 = h*f(t(i),w(i) );
        k2 = h*f(t(i+1),w(i) +k1);
        w(i+1)=w(i)+0.5*(k1+k2);
    end
end