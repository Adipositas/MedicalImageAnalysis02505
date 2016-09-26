function [mean, Co_Var, princ_ax11, princ_ax12, princ_ax21, princ_ax22] = MeanCovar(im)

Nx=size(im,2);
Ny=size(im,1);

x=repmat(1:Nx,1,Ny);
% Kronecker tensor product, taking all possible products between A and B
y=kron(1:Ny, ones(1,Nx));
vectorized=[x;y]';

% Mean of each axis
mean=sum(vectorized.*repmat(im(:),1, 2),1)/sum((im(:)));
vectorized_n=vectorized'-repmat(mean',1,size(vectorized,1));

% Co-variance
Co_Var=zeros(2,2);

for l=1:size(vectorized,1)
    Co_Var=Co_Var+(im(l)*vectorized_n(:,l)*vectorized_n(:,l)');
end

Co_Var=Co_Var/sum(im(:));

% Estimate principal axis
% Eigen vectors and values of co-var
[eig_vec,eig_val]=eig(Co_Var);
eig_val=diag(eig_val);
auxiliar=eig_vec;
eig_vec(1,:)=auxiliar(2,:);
eig_vec(2,:)=auxiliar(1,:);

princ_ax11=mean'-eig_vec(:,1)*sqrt(eig_val(1));
princ_ax12=mean'+eig_vec(:,1)*sqrt(eig_val(1));
princ_ax21=mean'-eig_vec(:,2)*sqrt(eig_val(2));
princ_ax22=mean'+eig_vec(:,2)*sqrt(eig_val(2));


