clear;

%% Rosenbrock Function
x = -2.5:.05:2.5;  y = -1:.05:6;  [X,Y] = meshgrid(x,y);
V = 100*(Y-X.^2).^2+(1-X).^2; 

%% BFGS Method with Armijo Line Search
x0 = [-0.75,1]';                                                           % initial point
max_k = 100;                                                               % maximum iteration number
sigma = 0.5;  gamma = 0.4;  a = 2;                                         % initial estimate of a=2 for line search
k = 0; H0 = eye(2);
xk = x0;  xk_1 = x0;
algorithm(:,1) = x0; 

jk = log((xk(1)-1)^2 + (xk(2)-1)^2);
Jk(:,1) = jk;

[g0,~] = fun_grad(x0);
[g0k_1,~] = fun_grad(x0);
g = g0;  gk_1 = g0k_1;
r0 = g0-g0k_1;  delta0 = x0-x0;
v0 = (r0'*H0*r0)^0.5*(delta0/(delta0'*r0)-(H0*r0)/(r0'*H0*r0)); vk = v0;   % Initialization 
                                                                  
while k < max_k

    if k == 0
        Hk = H0;
    else
        [g,~] = fun_grad(xk);
        [gk_1,~] = fun_grad(xk_1);
        rk = g-gk_1;  deltak = xk-xk_1;
        vk = (rk'*Hk*rk)^0.5*(deltak/(deltak'*rk)-Hk*rk/(rk'*Hk*rk));
        Hk = Hk+deltak*deltak'/(deltak'*rk)-Hk*(rk*rk')*Hk/(rk'*Hk*rk)+vk*vk';
    end
    
    dk = -Hk*g;
    
    alphak = armijo(xk, sigma, gamma, dk, a);                              % Armijo Line search
    xk_1 = xk;  xk = xk+alphak*dk;
    k = k+1;
    algorithm(:,k+1) = xk;
    
    jk = log((xk(1)-1)^2 + (xk(2)-1)^2);
    Jk(:,k+1) = jk;    

    if g == 0
        break;
    end

    
    k                                                                      % number of iteration
    xk                                                                     % iteration point                                                     % norm of the direction
    val = fun_obj(xk)                                                      % value of V at iteration point
end

%% Visualization
figure(1)
levels = [0.1,0.25,0.5,1,5,10,25,50];
contour(X,Y,V,levels,'LineWidth',1,'LineColor','#0072BD','ShowText','on');
hold on; plot(algorithm(1,:),algorithm(2,:),'r-.','LineWidth',1.5); 
xlabel('x');  ylabel('y');  title('BFGS Method with Armijo Line Search'); 
legend('Level sets','Behavior of BFGS method');

figure(2)
plot(Jk(1,:));
xlabel('Iteration number: k');  ylabel('Jk');  title('Speed of Convergence of BFGS Method with Armijo Line Search');