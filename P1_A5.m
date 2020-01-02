clear;

%% Rosenbrock Function
x = -2.5:.05:2.5;  y = -1:.05:6;  [X,Y] = meshgrid(x,y);
V = 100*(Y-X.^2).^2+(1-X).^2; 

%% Polak-Ribiere Method with Armijo Line Search
x0 = [-0.75,1]';                                                           % initial point
max_k = 37;                                                                % maximum iteration number $89$
sigma = 0.5;  gamma = 0.35;  a = 1;                                        % initial estimate of a=2 for line search
k = 0; %epsilon = 1e-5;
xk = x0;  xk_1 = x0;
algorithm(:,1) = x0; 

jk = log((xk(1)-1)^2 + (xk(2)-1)^2);
Jk(:,1) = jk;

while k < max_k
    [g,~] = fun_grad(xk);                                                  % compute the gradient
    [gk_1,~] = fun_grad(xk_1);
    
%     if norm(dk) < epsilon
%         break; 
%     end    
    
    if k == 0
        dk = -g;
    else
        dk_1 = dk;
        dk = -g+((g'*(g-gk_1)*dk_1)/(norm(gk_1))^2);
    end
    
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
%     norm_dk = norm(dk)                                                     % norm of the direction    
    val = fun_obj(xk)                                                      % value of V at iteration point
end

%% Visualization
figure(1)
levels = [0.1,0.25,0.5,1,5,10,25,50];
contour(X,Y,V,levels,'LineWidth',1,'LineColor','#0072BD','ShowText','on');
hold on; plot(algorithm(1,:),algorithm(2,:),'r-.','LineWidth',1.5); 
xlabel('x');  ylabel('y');  title('Polak-Ribiere Method with Armijo Line Search'); 
legend('Level sets','Behavior of Polak-Ribiere method');

figure(2)
plot(Jk(1,:));
xlabel('Iteration number: k');  ylabel('Jk');  title('Speed of Convergence of Polak-Ribiere Method with Armijo Line Search');