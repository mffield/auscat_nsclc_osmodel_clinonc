function [lik, sumError, prob] = test_logreg(w, X, y)

y=y(:)'; 

%calculate likelihood
temp2 = logsumexp((-w*X.*y)');
lik  = -sum(log(1+exp(-w*X.*y)));

y0 = y; y0(y0==-1)=0;
lik2 = sum(y0.*log(sigmoidfun(w*X))+(1-y0).*log(1-sigmoidfun(w*X)));

prob=1./(1+exp(-w*X));
y_=double(1./(1+exp(-w*X))<0.5); y_(y_==0)=-1;
sumError = sum((y(:)-y_(:))==0);
   
if any(isinf(lik))
    disp('inf')
    
end
end

