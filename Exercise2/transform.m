function [y_trans,s,R,t,sigma_FRA]=transform(x,y)

% transpose to vector, input expected as coloumns
x = x.';
y = y.';

mu_x=mean(x,2);
mu_y=mean(y,2);

x_=x-repmat(mu_x,1, length(y));
y_=y-repmat(mu_y,1, length(x));

% Singular value decomposition of H as stated on page 20 MIA course notes
H=x_*y_.';
[U,~,V]=svd(H);

xz=zeros(length(x),1);

% Using R definition
R=V*diag(det(V*U))*U.';

for i=1:size(x,2)
    xz(i)=x_(:,i).'*R*y_(:,i);
end

% s difinition
s=sum(xz)/(sum(sum(x_.^2,1),2));

% t definition
t=mu_y-s*R*mu_x;

y_trans=s*R*x+repmat(t,1, length(x));
sigma_FRA=mean(sum((y_trans-y).^2))/(2*length(x)); 
