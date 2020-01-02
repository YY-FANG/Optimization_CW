clear;

%% Rosenbrock Function
x = -2.5:.05:2.5;  y = -1:.05:6;  [X,Y] = meshgrid(x,y);
V = 100*(Y-X.^2).^2+(1-X).^2; 

%% Gradient Method
x0 = [-0.75,1]';                                                           % initial point
max_k = 1115;                                                              % maximum iteration number $2339$
sigma = 0.5;  gamma = 0.4;  a = 2;                                         % initial estimate of a=2 for line search
k = 0;  epsilon = 1e-5;
xk = x0; algorithm(:,1) = x0; 

jk = log((xk(1)-1)^2 + (xk(2)-1)^2);
Jk(:,1) = jk;

while k < max_k
    [g,~] = fun_grad(xk);                                                  % compute the gradient
    dk = -g;                                                               % compute the search direction
    if norm(dk) < epsilon
        break; 
    end
    
    alphak = armijo(xk, sigma, gamma, dk, a);                              % Armijo Line search
    
    xk = xk+alphak*dk;
    k = k+1;
    algorithm(:,k+1) = xk;
    
    jk = log((xk(1)-1)^2 + (xk(2)-1)^2);
    Jk(:,k+1) = jk;    
    
    if g == 0
        break;
    end
    
    k                                                                      % number of iteration
    xk                                                                     % iteration point
    norm_dk = norm(dk)                                                     % norm of the direction
    val = fun_obj(xk)                                                      % value of V at iteration point
end

%% Visualization
figure(1)
levels = [0.1,0.25,0.5,1,5,10,25,50];
contour(X,Y,V,levels,'LineWidth',1,'LineColor','#0072BD','ShowText','on');
hold on; plot(algorithm(1,:),algorithm(2,:),'r-.','LineWidth',1.5); 
xlabel('x');  ylabel('y');  title('The Gradient Algorithm with Armijo Line Search'); 
legend('Level sets','Behavior of the gradient algorithm');

figure(2)
plot(Jk(1,:));
xlabel('Iteration number: k');  ylabel('Jk');  title('Speed of Convergence of The Gradient Method with Armijo Line Search');