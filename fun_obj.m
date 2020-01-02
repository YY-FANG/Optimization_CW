function V = fun_obj(x)                                                    
% compute the value of the function at point x

X = x(1);  Y = x(2);

V = 100*(Y-X.^2).^2+(1-X).^2; 

end