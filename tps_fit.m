function [alpha,beta,Kreg,K,P,SMat,y,Z,df,SMx] = tps_fit(pos,val,lambda)

c1=pos(:,1);
c2=pos(:,2);
n = length(pos);
K = zeros(n);
R= zeros(n);

% Calculating the radial basis function values.
for i=1:n
    for j=1:n
          U=sqrt((c1(i,1)-c1(j,1))^2+(c2(i,1)-c2(j,1))^2);
        if U==0
            K(i,j)=0;
        else 
            K(i,j)=U^2*log(U);
        end
    end
end

%Calculating Smoothing Matrix      
P=[ones(n,1),c1,c2]';
Kreg=K+lambda.*eye(n,n);
O=zeros(3);

SMat=[Kreg,P';P,O];
y=[val;0;0;0];
SMx=[K P']*inv(SMat);
df=trace(SMx(1:n,1:n));


%Extracting parameters
alpha=SMat\y;
beta=alpha(end-2:end);
alpha=alpha(1:end-3);


no_lines=256; %number of lines in grid, i.e. size of image 
x1 = linspace(1,256,no_lines);
x2 = linspace(1,256,no_lines);
Z=zeros(no_lines);
wU=zeros(n,1);

% Fitting the parameters to space
for i=1:no_lines
    for j=1:no_lines
        for k=1:n
            
              U=sqrt((x1(i)-c1(k,1))^2+(x2(j)-c2(k,1))^2);
            if U~=0
                wU(k)=alpha(k)*U^2*log(U);
            end
        end
        Z(i,j)=beta(1)+beta(2)*x1(i)+beta(3)*x2(j)+sum(wU);
    end

end



end          
            