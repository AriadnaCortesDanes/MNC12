function [s,err] = myfunction(m)
% IN: m --> enter
% OUT: s --> sumatori 
%      err --> error
s=0;
for n = 0:m
    s = s + (-1/3)^n/(2*n+1);
end
s = s * sqrt(12);
err = abs(s-pi);
end