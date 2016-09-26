%% Set-up
load('visiblehuman.mat');

subplot(2,2,1);
imshow(head_fresh);
title('CT scan');

subplot(2,2,2);
imshow(head_mri);
title('MR scan');

subplot(2,2,3);
imshow(head);
title('Color cryosection');

subplot(2,2,4);
imshow(head_frozen);
title('CT scan of frozen cadaver');
%% Task 2.1

%[CT_points] = cpselect(head, head_frozen);
% Save for easy restart
% save ref_points.mat ref_points
% save ref_points2.mat ref_points2
% save ref_points3.mat ref_points3
% save MR_points.mat MR_points;
% save CT_points.mat CT_points;

load('ref_points.mat');
load('ref_points2.mat');
load('ref_points3.mat');

load('MR_points.mat');
load('CT_points.mat');

mu_points = (ref_points + ref_points2) / 2;

FLE = (ref_points - mu_points).^2;
FLE_xy = sum(FLE) / length(FLE);

sigma_FLE =sum(sum((FLE)))/(length(FLE)*2);

figure();
hold on;
set(gca,'Ydir','reverse')
xlim([0 size(head,2)]);
ylim([0 size(head,1)]);
imagesc(head);

scatter(ref_points(:,1), ref_points(:,2),'filled', 'LineWidth', 2, 'markerfacecolor', 'red');
scatter(ref_points2(:,1), ref_points2(:,2),'filled' , 'LineWidth', 2, 'markerfacecolor', 'blue');

%% Task 2.3
% MR to CT
subplot(1,2,1);
imshow(head_frozen);
title('CT scan of frozen cadaver');

subplot(1,2,2);
imshow(head_mri);
title('MR scan of the cadaver');

% MR to CT
[y_trans_MR, s_MR, R_MR, t_MR, sigma_FRA_MR] = transform(MR_points, CT_points);

hold on;
scatter(y_trans_MR(1,:), y_trans_MR(2,:),'markerfacecolor','red');
scatter(CT_points(:,1), CT_points(:,2).','markerfacecolor','blue');
legend('show');
legend('Transformed MR','reference CT');

% CT to MR
[y_trans_CT, s_CT, R_CT, t_CT, sigma_FRA_CT] = transform(CT_points, MR_points);

hold on;
scatter(y_trans_CT(1,:), y_trans_CT(2,:),'markerfacecolor','blue');
scatter(MR_points(:,1), MR_points(:,2).','markerfacecolor','red');
legend('show');
legend('Transformed CT','reference MR');

%% Task 2.4
% Mr image transformd to CT

% Grid, CT size (we want to end with CT size)
x11 = linspace(1, size(head_frozen,2), size(head_frozen,2));
x22 = linspace(1, size(head_frozen,1), size(head_frozen,1));
[x1, x2] = meshgrid(x11, x22);

x1z = zeros(size(x1));
x2z = zeros(size(x2));
grid_points = [x1(:),x2(:)]';

% Transform with params
% Interpolate intensities, based on grid location
grid_transformed = s_CT*R_CT*grid_points+repmat(t_CT,1,size(grid_points,2));
int_ctspace = interp2(double(head_mri),grid_transformed(1,:),grid_transformed(2,:),'bilinear');

% Reshape to image size
x1_mr = reshape(grid_transformed(1,:), size(x1));
x2_mr = reshape(grid_transformed(2,:), size(x2));
int_ctspace = reshape(int_ctspace, size(head_frozen));

% Gather images
im_joint = zeros (size(head_frozen,1), size(head_frozen,2), 3);
im_joint(:,:,1) = head_frozen;
im_joint(:,:,2) = int_ctspace;

imshow(im_joint./255);
title('CT as red, MR as green');
