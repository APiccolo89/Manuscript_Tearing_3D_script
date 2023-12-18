function [obj,Phase,Temp] = find_slab_mode_2(obj,A,Weak_Slab,Phase,Temp,Boundary,theta)
% linear part of the slab: should I compute?
sl    = 1    ; % Always YES, NO, only if we are considering weak zone, which has 90 deg theta
% Convert the angle in radians

x     = A.Xpart(:); % x part coordinate
z     = A.Zpart(:); % z part coordinate
y     = A.Ypart(:); % y part coordinate
data_boundary = obj.Boundary.(Boundary);
%Select the boundary = limits of the boundary
[loc1, ~] = find(tril(data_boundary{1} == data_boundary{1}', -1)); % find the main boundary.
if rem(loc1,2)==0
    ya = data_boundary{1}(1);
    yb = data_boundary{1}(3);
    xc = data_boundary{1}(2);
else
    ya = data_boundary{1}(2);
    yb = data_boundary{1}(4);
    xc = data_boundary{1}(3);

end
% Restrict our research to the point the t ARE WORTH TO LOOK. You do not
% find God in a Casino, but in church.
C      = [xc,0.0-obj.R(2)]; % center of the curvature
dl = obj.L0./obj.nseg;
% Spell out the surfaces
T = obj.Slab_surface.Top; % Top surface

M = obj.Slab_surface.MidS; % MidSurface

B = obj.Slab_surface.Bottom; % Bottom

WZ = obj.Slab_surface.WZ_surf;

% Select better the area of the chosen one
if ~strcmp(Weak_Slab,'Weak')
ind = x>=C(1) & x(:)<=max(T(:,1)) & z(:)>=min(B(:,2)) & z(:)<=1.0 & y>ya & y<yb;
else
ind = x>=C(1) & x(:)<=max(WZ(:,1)) & z(:)>=obj.Decoupling_depth(1) & y>ya & y<yb;
end


%Prepare the main vectors
d_slab = x.*nan;

l_slab = x.*nan;

continent = x.*nan;
%Prepare other variables

xs = x(ind==1);

zs = z(ind==1);

ds = d_slab(ind==1);

ls = l_slab(ind==1);

cs = l_slab(ind==1).*nan;
%free memory
x = [];
z = [];
y = [];

it = 1;

l = 0;


ind_decoupling_1 = find(M(:,2)>obj.Decoupling_depth,1);



CPU_AF = cputime;

while l<obj.L0
    % = Marcel suggested to use scatterinterpolant and inpolygon. Now the
    % routine is faster.
    %======================================================================

    ln = l+dl;

    itn = it+1;
    %[Spell out the coordinate of the small element]
    if ~strcmp(Weak_Slab,'Weak')
        Ap = [T(it,1),T(it,2)]; % coordinate point A

        Bp = [T(itn,1),T(itn,2)];% coordinate point B

        Cp = [B(it,1),B(it,2)]; %coordinate point C

        Dp = [B(itn,1),B(itn,2)]; % coordinate D;
    else
        Ap = [WZ(it,1),WZ(it,2)]; % coordinate point A

        Bp = [WZ(itn,1),WZ(itn,2)];% coordinate point B

        Cp = [T(it,1),T(it,2)]; %coordinate point C

        Dp = [T(itn,1),T(itn,2)]; % coordinate D;

    end


    % Create the polygon

    xv = [Ap(1),Bp(1),Dp(1),Cp(1)];

    zv = [Ap(2),Bp(2),Dp(2),Cp(2)];
    % Find the chosen one

    [ind_chosen,~] = inpolygon(xs,zs,xv,zv);

    %Set up the scatter interpolant classes
    F1  = scatteredInterpolant(xv',zv',[l l+dl l+dl l]');

    F2 = scatteredInterpolant(xv',zv',[0 0 -obj.D0 -obj.D0]');

    %Compute the length and d as l with linear interpolation
    ls(ind_chosen==1) = F1(double(xs(ind_chosen==1)),double(zs(ind_chosen==1)));

    ds(ind_chosen==1) = F2(double(xs(ind_chosen==1)),double(zs(ind_chosen==1)));

    if it == ind_decoupling_1
        l_dc = l;
    end

    it = itn;
    l = ln;
   
end
CPU_BF = cputime;

% Compute the poligon of the continental crust
disp(['Slab segments identification took ',num2str(CPU_BF-CPU_AF,2), 'seconds'])

if ~strcmp(Weak_Slab,'Weak')
if ~strcmp(obj.length_continent{2},'none')

    [L] = Compute_properties_along_function(ya,yb,obj.length_continent,A.Length_along(ind==1));
    Points = obj.Subducted_crust_L;

    Pl(1)     = Points(1);

    Pl(2)     = Points(3);

    Pl(3)     = 1.0;

    Pd(1)     = Points(2);

    Pd(2)     = Points(4);

    Pd(3)     = Points(6);

    [in,~] = inpolygon(ls./L(:),ds,Pl,Pd);

    cs(in==1)=1.0;
else
    Points = obj.Subducted_crust_L;

    Pl(1)     = Points(1);

    Pl(2)     = Points(3);

    Pl(3)     = Points(5);

    Pd(1)     = Points(2);

    Pd(2)     = Points(4);

    Pd(3)     = Points(6);

    [in,~] = inpolygon(ls,ds,Pl,Pd);

    cs(in==1)=1.0;

    
end
% Fill up the main portion 
d_slab(ind==1) = ds; 

l_slab(ind==1) = ls;

continent(ind==1) = cs; 

% Correct Phase and Temperature 

Temp(ind==1 & isnan(d_slab) & A.Zpart(:)<-obj.D0)= obj.Thermal_information.TP;

Phase(ind==1 & isnan(d_slab) & A.Zpart(:)<-obj.D0) = obj.Thermal_information.Ph_Ast;

% Find the first particles fullfilling the required dept

ind_decoupling = find(l_slab>=l_dc,1);

% Update 
obj.l_slab = reshape(l_slab,size(A.Xpart));

obj.d_slab = reshape(d_slab,size(A.Xpart));

obj.continent = reshape(continent,size(A.Xpart));

obj.Decoupling_depth(1) = obj.Decoupling_depth;

obj.Decoupling_depth(2) = ind_decoupling;


else
d_slab(ind==1) = ds; 
obj.d_slab = reshape(d_slab,size(A.Xpart));
end



bla = 0;

end































% 
% function [d] = find_distance_linear(P,p1,p2)
% 
% % Ok, wiki saves my day because I was too lazy:
% A  = abs((p2(1)-p1(1)).*(p1(2)-P(:,2))-(p1(1)-P(:,1)).*(p2(2)-p1(2)));
% B  = (p2(1)-p1(1)).^2+(p2(2)-p1(2)).^2;
% d  = A./sqrt(B);
% 
% end
% 
% function [l,ind_decoupling,Phase,Temp] = find_length_slab(obj,x,z,C,r,d,l,Phase,Temp,theta)
% slab = ~isnan(d)==1;
% % Compute the two point along the midplane of the slab
% % 1) After the curved part of the slab
% PA(1)=C(1)+r*sin(theta*pi/180);
% PA(2)=C(2)+r*cos(theta*pi/180);
% %Coordinate of the midsurface at the depth of slab  (bottom of the slab)
% PB(1)=PA(1)+obj.L0*cos(theta*pi/180);
% PB(2)=PA(2)-obj.L0*sin(theta*pi/180);
% % Compute the slope of the line:
% m=(PB(2)-PA(2))./(PB(1)-PA(1));
% %m = 1./m;
% % project all the x points into the line (from:
% % https://math.stackexchange.com/questions/62633/orthogonal-projection-of-a-point-onto-a-line)
% % Project all the points that belong to the slab into the mid surface
% xprojection = d.*0.0;
% zprojection = d.*0.0;
% b = PA(2)-PA(1).*m;
% xprojection(slab==1 & isnan(l)) = (x(slab==1 & isnan(l)) + m.*z(slab==1 & isnan(l))-m.*b)./(1+m.^2);
% zprojection(slab==1 &isnan(l)) = (m.*x(slab==1 & isnan(l)) + m^2.*z(slab==1 & isnan(l))+b)./(1+m.^2);
% % Compute the length of the slab along the midsurface
% l(slab==1 & isnan(l)) = r.*theta*pi/180+sqrt((xprojection(slab==1 & isnan(l))-PA(1)).^2+(zprojection(slab==1 & isnan(l))-PA(2)).^2);
% ind_decoupling = find(zprojection>=obj.Decoupling_depth,1);
% % Corret the damn Phase and Temperature
% Phase(z<PA(2)+m*(x-PA(1)) & isnan(d) & z<PA(2)& x>C(1) & ~isnan(Phase(:))) = obj.Thermal_information.Ph_Ast;
% Temp(z<PA(2)+m*(x-PA(1)) & isnan(d) & z<PA(2) & x>C(1) &~isnan(Phase(:))) = obj.Thermal_information.TP;
% end
% 
% function [z] = linear_equation_two_points(s,A,x)
% 
% z = s.*(x-A(1))+A(2);
% 
% end
% 
% function [length] = compute_real_length(s,M,x,z,l)
% 
% % project all the x points into the line (from:
% % https://math.stackexchange.com/questions/62633/orthogonal-projection-of-a-point-onto-a-line)
% % Project all the points that belong to the slab into the mid surface
% 
% prx = x.*0.0;
% 
% prz = z.*0.0;
% 
% b = M(2)-M(1).*s;
% 
% prx = (x + s.*z-s.*b)./(1+s.^2);
% 
% prz = (s.*x + s^2.*z+b)./(1+s.^2);
% % Compute distance from l along the midline
% 
% d = sqrt((prx-M(1)).^2+(prz-M(2)).^2);
% 
% length = l+d;
% 
% end