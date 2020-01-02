function [g,g2] = fun_grad(x)                                                
% compute the gradient of function at the point x

syms X Y;

V = 100*(Y-X.^2).^2+(1-X).^2;  

grad = gradient(V);  grad2 = jacobian(grad,[X;Y]);

X = x(1);  Y = x(2);

g = eval(grad);  g2 = eval(grad2);

end

