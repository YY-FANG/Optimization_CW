clear;

%% Rosenbrock Function
x = -2.5:.05:2.5;  y = -1:.05:6;  [X,Y] = meshgrid(x,y);
V = 100*(Y-X.^2).^2+(1-X).^2; 

%% Visualization 
figure(1)
surf(X,Y,V);  xlabel('x');  ylabel('y');  zlabel('v(x,y)'); 
title('Plot of Rosenbrock Function');

%% Level Sets
figure(2)
levels = [0.1,0.25,0.5,1,5,10,25,50];
contour(X,Y,V,levels,'LineWidth',1,'LineColor','#0072BD','ShowText','on');
xlabel('x');  ylabel('y');  title('Level Sets of Rosenbrock Function');