% ========================================================================%
% LaMEM grid generator                                                    %
%=========================================================================%
nel_total = 2048;
grid = [-1500.0,-400.0,400.0,1500.0];
nel_seg = [512,1024,512]; 
min_r   = (grid(3)-grid(2))./nel_seg(2);
mean_r   = (grid(2)-grid(1))./(nel_seg(1));
bias(1) = min_r./(2.*mean_r-min_r); 
max_r   = (2.*mean_r)./(1+bias(1));
bias(2) = 1.0; 
bias(3) = max_r./min_r; 
dx = zeros(nel_total,1); 
bias
for i=1:1:nel_total
if i<=nel_seg(1)
    if i ==1 
        dx(i) = max_r; 
    else
        dx(i)= max_r+((min_r-max_r)./(nel_seg(1))).*i;
    end
elseif i>nel_seg(1) && i <= nel_seg(2)+nel_seg(1)
        dx(i) = min_r; 
else
        dx(i)=min_r+((max_r-min_r)./(nel_seg(3))).*(i-sum(nel_seg(1:2)));
end
end
