function [ M ] = axisrot( axis, t )
% [ M ] = axisrot( axis, t )
% AXISROT construct a 3x3 axis rotation matrix

fcos = cos(t);
fsin = sin(t);
fminuscos = 1 - fcos;
fx2 = axis(1)*axis(1);
fy2 = axis(2)*axis(2);
fz2 = axis(3)*axis(3);
fxym = axis(1)*axis(2)*fminuscos;
fxzm = axis(1)*axis(3)*fminuscos;
fyzm = axis(2)*axis(3)*fminuscos;
fxsin = axis(1)*fsin;
fysin = axis(2)*fsin;
fzsin = axis(3)*fsin;
M = [  fx2*fminuscos + fcos,  fxym-fzsin,  fxzm+fysin;
       fxym+fzsin,  fy2*fminuscos+fcos,    fyzm-fxsin;
       fxzm-fysin,  fyzm+fxsin,  fz2*fminuscos+fcos ];



%  matlab version...
%M = makehgtform('axisrotate',axis,t);
%M = M(1:3,1:3);

end
