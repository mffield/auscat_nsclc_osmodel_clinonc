function p = calculatePvalueProportion(n1, p1, n2, p2)

pmean=(n1*p1+n2*p2)/(n1+n2);
sd=sqrt(pmean*(1-pmean)*(1/n1+1/n2));
z=abs((p1-p2)/sd);
p = 2*(1-normcdf(z,0,1));

end