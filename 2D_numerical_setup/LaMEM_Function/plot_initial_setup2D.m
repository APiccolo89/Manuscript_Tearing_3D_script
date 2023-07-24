function plot_initial_setup2D(A,surf)
figure(1)
ratioX = max(A.Xpart(:,1,1))-min(A.Xpart(:,1,1));
ratioZ = max(A.Zpart(1,1,:))- min(A.Zpart(1,1,:));
r = [1,ratioZ./ratioX,1];
x = squeeze(A.Xpart(:,1,:));
z = squeeze(A.Zpart(:,1,:));
T = squeeze(A.Temp(:,1,:)); 
ph = squeeze(A.Phase(:,1,:)); 
T(ph==0)=nan;
ph(ph==0)=nan;
hold on 
plot(surf(1,:),surf(2,:),"Color",'r','LineStyle','--','LineWidth',1.2)
pcolor(x,z,T);
hold off
shading interp; 
pbaspect(r); 
xlabel('x, [km]', Interpreter='latex'); 
ylabel('z, [km]',Interpreter='latex'); 
title('Temperature Field [$^\circ$ C]',Interpreter='latex');
grid on; 
colormap(crameri('bilbao',13)); 
box on 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; 
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; 
c = colorbar;
caxis([0 1350])
c.Label.String = 'Temperature $[^\circ C]$';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter='latex';
print('Temperature','-dpng','-r0')


figure(2)
lv = 1:12;
hold on 
plot(surf(1,:),surf(2,:),"Color",'r','LineStyle','--','LineWidth',1.2)
contourf(x,z,ph,13); 

hold off
shading interp; 
pbaspect(r); 
xlabel('x, [km]', Interpreter='latex'); 
ylabel('z, [km]',Interpreter='latex'); 
title('Phase Field [$^\circ$ C]',Interpreter='latex');
grid on
caxis([1,13]); 
box on 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; 
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; 
xlim([-300 300]);
ylim([-100 10])
colormap(crameri('nuuk',13))
c = colorbar;%(['Upper Crust 1','Lower Crust 1', 'Lithospheric Mantle 1', ' Lithospheric Mantle 2', 'Upper Mantle','Oceanic Slab', 'Oceanic Crust', 'Weak Zone', 'Upper Crust 2', 'Lower Crust 2', 'Oceanic Sed.', 'Passive/Prism']);
c.Label.String = 'Phase ';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter='latex';
print('Phase','-dpng','-r0')

figure(3)
plot(surf(1,:),surf(2,:),"Color",'r','LineStyle','--','LineWidth',1.2)
xlim([-600 600])
bla = 0;
end