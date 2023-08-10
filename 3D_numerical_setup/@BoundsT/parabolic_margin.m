

function [y,W]  = parabolic_margin(xa,xv,yv,x)

% Function to describe the linear margin of a terranes 
a = -yv./(xa-xv).^2;
y = a.*(x-xv).^2+yv; 
W = length_parabolic_margin(xa,xv,yv,x);
end

function [W] = length_parabolic_margin(xa,xv,yv)
a = -yv./(xa-xv).^2;
xb  = -xa; 
t  = (xa-xv);
t2 = (xb-xv); 
W1 = (2.*a.*t.*(4.*a.*a.*t.^2+1)+asinh(2.*a.*t))./(4.*a);
W2 = (2.*a.*t2.*(4.*a.*a.*t2.^2+1)+asinh(2.*a.*t2))./(4.*a);
W = W2-W1;
end