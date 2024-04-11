function [model, loglik, error, valloglik, valerror] = update_logreg_model(model)

Nboot = length(model.settings.randseed);
NModel = size(model.w,1)/Nboot;

for j=1:NModel
 
    gradient = zeros(size(model.client(1).info(j).gradient));
    H = zeros(size(model.client(1).info(j).hessian));
    loglik = zeros(size(model.client(1).info(j).loglik));
    valloglik = zeros(size(model.client(1).info(j).valloglik));
    error = zeros(size(model.client(1).info(j).loglik));
    valerror = zeros(size(model.client(1).info(j).valloglik));
    for i=1:length(model.client)
        gradient = gradient + model.client(i).info(j).gradient;
        H = H + model.client(i).info(j).hessian;
        loglik = loglik + model.client(i).info(j).loglik;
        valloglik = valloglik + model.client(i).info(j).valloglik;
    end
    
    if min(size(model.settings.varselect))==1
        varselect = logical(model.settings.varselect);
    else
        varselect = logical(model.settings.varselect(j,:));
    end
    
    for k=1:Nboot
        m_ind = ((j-1)*Nboot)+k; %current model index
        d = sum(varselect)+1;
        if d<size(model.w,2)
            model.w(m_ind,d+1:end)=0;model.g(m_ind,d+1:end)=0;model.u(m_ind,d+1:end)=0;
        end
        % get current parameters and previous gradient
        u=model.u(m_ind,[true true(1,sum(varselect))]);
        g_ = model.g(m_ind,[true true(1,sum(varselect))]);
        w_ = model.w(m_ind,[true true(1,sum(varselect))]);
        w = model.w(m_ind,[true true(1,sum(varselect))]);
        
        %combine likelihoods
        loglik(k) = loglik(k) - (model.settings.lambda/2)*(w*w');
        valloglik(k) = valloglik(k) - (model.settings.lambda/2)*(w*w');
        
        %     error(k) = sum(model.trainSumError(k,:));
        %     valerror(k) = sum(model.valSumError(k,:));
        if any(isnan(loglik(:)))
            disp('nan')
        end
        if ~model.convergence(m_ind)
            
            lambda_vec = [0; model.settings.lambda*ones(length(w)-1,1)]; %do not penalize bias term
            
            
            
            Hk = reshape(H(k,1:d^2),[d d]);
            g = gradient(k,[true true(1,sum(varselect))]) - lambda_vec(:)'.*w(:)';
            
            uTHu = (u*Hk*u') - model.settings.lambda*(u*u');
            
            %         g = sum(model.localGrad(:,k),length(size(model.localGrad))) - lambda_vec.*w;
            
            %         if model.w_fixed(1)
            %             g(1)=0; u(1)=0;
            %         elseif model.w_fixed(2)
            %             g(2)=0; u(2)=0;
            %         end
            %
            %         %combine Hessian terms
            %         uTHu = sum(model.localHess(k,:)) - model.lambda*(u'*u);
            
            %update w, u, g
            model.w(m_ind,1:d) = w_ - ((g*u')./uTHu)*u;
            model.u(m_ind,1:d) = g - u*model.kappa(k);
            model.kappa(m_ind) = g*(g-g_)'/(u*(g-g_)');
            model.g(m_ind,1:d) = g;
            
            if sum((model.w(m_ind,1:d)-w_).^2)<model.tol
                disp(['Model ' num2str(m_ind) ' converged within tolerance at iteration: ' num2str(model.itr)])
                model.convergence(m_ind)=true;
            end
        end
    end
    
end