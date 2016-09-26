
load('abdomen.mat');
load('inputpoints.mat');
pos_rev = pos(:,2:-1:1);
n = length(pos);

for i=1:n
    samples(i,1)=abdomen(floor(pos(i,2)),floor(pos(i,1)));
end

c1=pos_rev(:,1);
c2=pos_rev(:,2);

K = zeros(n);

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

lambda = 1390;
P=[ones(n,1),c1,c2]';
y=[samples;0;0;0];
L = [K + eye(n,n)*lambda, P'; P, zeros(3)];

SMx=[K P']*inv(L);
df=trace(SMx(1:n,1:n));
sol = L \ double(y);

beta=sol(end-2:end);
alpha=sol(1:end-3);

ress=256; %number of lines in grid, i.e. size of image 
x1 = linspace(1,256,ress);
x2 = linspace(1,256,ress);
ints=zeros(ress);
wU=zeros(n,1);

% Fitting the parameters to space
for i=1:ress
    for j=1:ress
        for k=1:n
            
              U=sqrt((x1(i)-c1(k,1))^2+(x2(j)-c2(k,1))^2);
            if U~=0
                wU(k)=alpha(k)*U^2*log(U);
            end
        end
        ints(i,j)=beta(1)+beta(2)*x1(i)+beta(3)*x2(j)+sum(wU);
    end

end
figure()
surf(double(roi).*ints)
title('Bias field');

sca1=max(max(double(roi.*abdomen)));
sca2=max(max(double(roi.*abdomen)./ints));
scale=sca1/sca2;

figure()
imshow((double(roi.*abdomen)./ints)*scale,[]) %Final corrected image
title('Final corrected image');
