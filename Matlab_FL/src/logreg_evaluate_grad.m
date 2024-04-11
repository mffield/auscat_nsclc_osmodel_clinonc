function [gradient, hessian, L] = logreg_evaluate_grad(w,X,y,lambda)       

N = size(X,2);

gradient = [(ones(1,N) - 1./(1+exp(-w*X.*y))).*y*X'];% calculate gradient
A = (1./(1+exp(-w*X.*y))).*(1-1./(1+exp(-w*X.*y)));% calculate Hessian
D = (1./(1+exp(-w*X))).*(1-1./(1+exp(-w*X)));% calculate Hessian

% hessian = -(sum(((u*X).^2).*A));

hessian = -X*diag(A)*X' - lambda*eye();
% diag(inv(-hessian))
L  = -sum(log(1+exp(-w*X.*y(:)')));

end

