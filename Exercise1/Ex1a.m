X = [1,1 ; 1,-1 ; -1,1 ; -1,-1];

y = [3 ; 1 ; 1; 2];
y_pad = [3; 1 ; 1; 2; 0 ; 0 ; 0];

P = [1, 1,1 ;1, 1, -1;1, -1, 1;1, -1, -1];
[m n] = size(P);
n = m;
K = [0, 4*log(2), 4*log(2), 8 * log(2*sqrt(2)); 
    4* log(2), 0, 8 *log(2*sqrt(2)), 4*log(2);
    4*log(2),8*log(2*sqrt(2)), 0, 4*log(2);
    8 * log(2*sqrt(2)), 4*log(2), 4*log(2), 0];

sam_lam =0:1000;
for n=1:size(sam_lam)
    lambda = sam_lam(n);
    L = [K + eye(m,n)*lambda, P; transpose(P), zeros(n-1,m-1)];

    sol = L \ y_pad;

    alpha = sol(1:m, :);
    beta = sol(m+1:m+3);

    grid_step = 0.1;
    x1 = -1.5:grid_step:1.5;
    x2 = -1.5:grid_step:1.5;
    [xx, yy] = meshgrid(x1, x2);
    Xgr =[reshape(xx,[],1) reshape(yy, [], 1)];

    Px = [ones(size(Xgr,1),1) Xgr(:,1) Xgr(:,2)].';

    Kx=zeros(size(Xgr(:,1),1),4);

    % Find distance from all points to reference points
    for k=1:size(Kx,1)
        c_a=(X(:,1) + 1i *X(:,2));
        c_b=(Xgr(k,1) + 1i *Xgr(k,2));

        b=repmat( c_b.', m, 1 );
        d = c_a - b;
        r=d.*conj( d );
        in_ = find(r >0);
        Kx(k,in_) = log(sqrt(r(in_))).*r(in_);
    end

    % Estimate grid values
    y_hat = [Kx Px.']*[alpha; beta];

    surf(xx, yy, reshape(y_hat, [size(xx,1),size(yy,1)]));
    str=sprintf('Thine Plate Spline \\lambda = %2.2f', lambda);
    title(str,'FontSize', 15);

    Slambda = [Kx Px.'].*inv([Kx + eye(size(Kx))*lambda, Px.'; Px, zeros(size(Kx,1)-1,m-1).']);
    freedom(n) = trace(Slambda);
end

plot(sam_lam,freedom,'LineWidth',2)
xlabel('\lambda','FontSize',12);
ylabel('df_l','FontSize',12);