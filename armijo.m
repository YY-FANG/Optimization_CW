function alphak = armijo(xk, sigma, gamma, dk, a)                             
% Armijo line search method

assert(sigma >= 0 && sigma <= 1);  assert( gamma >= 0 && gamma <= 0.5 );

[g,~] = fun_grad(xk);                                                      % compute the gradient

j = 0; max_j = 10000;

while j <= max_j
    alpha = a*sigma^j;  x = xk+alpha*dk;
    if fun_obj(x)-fun_obj(xk)<=gamma*alpha*g'*dk
        alphak = alpha;  
        break;
    end
    j = j + 1;
end

end
