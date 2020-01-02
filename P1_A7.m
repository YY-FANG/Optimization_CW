clear;

%% Rosenbrock Function
x = -2.5:.05:2.5;  y = -1:.05:6;  [X,Y] = meshgrid(x,y);
V = 100*(Y-X.^2).^2+(1-X).^2; 

%% Simplex Method
x0 = [-0.75,1]';                                                           % initial point
max_k = 1114;                                                              % maximum iteration number
k = 0;  e = 1e-25;  a = 1.9733; imp = 1;

f0 = fun_obj(x0);
stepsize = 0.5;
x1 = [x0(1)+stepsize;x0(2)];
x2 = [x0(1)+stepsize/2; x0(2)-sqrt(3)/2*stepsize];
points = [x0 x1 x2]; 

algorithm(:,1) = x0;  algorithm(:,2) = x1;  algorithm(:,3) = x2;

jk = log((points(1,1)-1)^2 + (points(2,1)-1)^2);
Jk(:,1) = jk;
                                                                  
while k < max_k
    average = (fun_obj(points(:,1))+fun_obj(points(:,2))+fun_obj(points(:,3)))/3;
    criterion = ((fun_obj(points(:,1))-average)^2+(fun_obj(points(:,2))-average)^2+(fun_obj(points(:,3))-average)^2)/3;
    
    F = [fun_obj(points(:,1));fun_obj(points(:,2));fun_obj(points(:,3))];
    [F_sort,indice] = sort(F);  F_max = F(indice(3));  xmax = indice(3);
    
    xc = (points(:,1)+points(:,2)+points(:,3))/3;    
    dk = xc-points(:,xmax); x_new = xc+a*dk;
    
    point = points;      
    
    if imp == 1
        gamma = 2; rho = 0.5;
        point(:,xmax) = x_new;
        F_im = [fun_obj(point(:,1));fun_obj(point(:,2));fun_obj(point(:,3))];
        [F_sort_im,indice_im] = sort(F_im);  
        F_max_im = F_im(indice_im(3));  xmax_im = indice_im(3);
        F_min_im = F_im(indice_im(1));  xmin_im = indice_im(1);
        if fun_obj(x_new) ~= F_max_im && fun_obj(x_new) ~= F_min_im
            points(:,xmax) = x_new;
        elseif fun_obj(x_new) == F_min_im
            xe = xc+gamma*(x_new-xc);
            if fun_obj(xe)<fun_obj(x_new)
                points(:,xmax) = xe;
            else
                points(:,xmax) = x_new;
            end
        elseif fun_obj(x_new) == F_max_im
            xn = x_new;
            while(1)
                xo = xc+rho*(xn-xc);
                if fun_obj(xo)<fun_obj(xn)
                    points(:,xmax) = xo;
                    break;
                else
                    xn = xo;
                end
            end
        end
    else
        points(:,xmax) = x_new;                                            % get rid of the maximum point 
    end
    
    k = k+1;
    algorithm(:,k+3) = points(:,xmax);
    
    jk = log((points(1,xmax)-1)^2 + (points(2,xmax)-1)^2);
    Jk(:,k+1) = jk;
    
    if criterion < e
        break;
    end

    k                                                                      % number of iteration
    xk = points(:,xmax)                                                    % iteration point                                                     % norm of the direction
    val = fun_obj(points(:,xmax))                                          % value of V at iteration point

end

%% Visualization
if imp == 1
    figure(1)
    levels = [0.1,0.25,0.5,1,5,10,25,50];
    contour(X,Y,V,levels,'LineWidth',1,'LineColor','#0072BD','ShowText','on');
    hold on; plot(algorithm(1,:),algorithm(2,:),'r-.','LineWidth',1.5); 
    xlabel('x');  ylabel('y');  title('Modified Simplex Method'); 
    legend('Level sets','Behavior of Modified Simplex method');

    figure(2)
    plot(Jk(1,:));
    xlabel('Iteration number: k');  ylabel('Jk');  title('Speed of Convergence of Modified Simplex Method');
else
    figure(1)
    levels = [0.1,0.25,0.5,1,5,10,25,50];
    contour(X,Y,V,levels,'LineWidth',1,'LineColor','#0072BD','ShowText','on');
    hold on; plot(algorithm(1,:),algorithm(2,:),'r-.','LineWidth',1.5); 
    xlabel('x');  ylabel('y');  title('Simplex Method'); 
    legend('Level sets','Behavior of Simplex method');

    figure(2)
    plot(Jk(1,:));
    xlabel('Iteration number: k');  ylabel('Jk');  title('Speed of Convergence of Simplex Method');  
end