function [UL, p] = KM_logrank(N1,D1,N2,D2, alpha, str1, str2)
UL=0;p=0;
if ~(N1(1)==0 || N2(1)==0)
    
    SC_km_D1 = D1;
    SC_km_D2 = D2;
    SC_km_D = SC_km_D1 + SC_km_D2;
    loc = SC_km_D>0; 
    SC_km_D = SC_km_D(loc);
    SC_km_D1 = SC_km_D1(loc);
    SC_km_D2 = SC_km_D2(loc);
    SC_km_N1 = N1; 
    SC_km_N2 = N2;
    SC_km_N = SC_km_N1 + SC_km_N2;
    SC_km_N = SC_km_N(loc);
    SC_km_N1 = SC_km_N1(loc);
    SC_km_N2 = SC_km_N2(loc);
    
    %Compute the difference between observed deaths for treatment 1 and
    %expected deaths in the hypothesis that the treatments are similar
    lr_stat = SC_km_D1 - (SC_km_N1).*(SC_km_D)./(SC_km_N);
    J = sum(lr_stat); UL=abs(J); %Log-rank statistic
    
    %Compute the contribute to the standard error
    std_err = prod([SC_km_N1 SC_km_N2 SC_km_D],2).*(SC_km_N-SC_km_D) ./ (SC_km_N.^2.*(SC_km_N-ones(size(SC_km_N,1),1)));
    std_err(isnan(std_err))=0;  %find if there is some NaN (i.e. 0/0)
    V=sum(std_err); SUL=sqrt(V); %Compute the total standard error
    K=J/V; HR=exp(K); HRci=[exp(K-1.956/SUL) exp(K+1.956/SUL)]; %NOTE: alpha should be considered here too!
    z=abs((UL-0.5)/SUL); %normalized UL with Yates' correction
    p=2*(1-0.5*erfc(-z/realsqrt(2))); %p-value

    if p < alpha
        disp(['For ' str1 ',' str2 ' group. Survival curves significantly different: p=' num2str(p)]);
    else
        disp(['For ' str1 ',' str2 ' group. Survival curves NOT significantly different: p=' num2str(p)]);
    end

else
    if N1(1)==0; disp([str1 ': truncated']); end
    if N2(1)==0; disp([str2 ': truncated']); end
end


end