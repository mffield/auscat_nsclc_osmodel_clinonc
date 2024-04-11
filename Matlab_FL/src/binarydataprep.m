function [X, y] = binarydataprep(data,features,target, threshold)

x = data(:,features);
N = size(x,1);
X = [ones(N,1),table2array(x)]';

y = double(target>=threshold);
y(y==0)=-1;
y=y';

end