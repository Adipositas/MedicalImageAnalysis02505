clear all; close all;
load('visiblehuman.mat')

subplot(1,2,1);
imshow(head_frozen);
title('CT Scan');

subplot (1,2,2);
imshow(head_mri);
title('MR scan');

% Double cast, to ensure mean and Co-var calculations

head_frozen = double(head_frozen);
head_mri = double(head_mri);

[mean_mr, co_var_mr, princ_ax11mr, princ_ax12mr, princ_ax21mr, princ_ax22mr] = ...
    MeanCovar(head_mri);
[mean_ct, co_var_ct, princ_ax11ct, princ_ax12ct, princ_ax21ct, princ_ax22ct] = ...
    MeanCovar(head_frozen);

figure;
imagesc(head_mri);
colormap('Gray')
hold on;

plot([princ_ax11mr(1),princ_ax12mr(1)],[princ_ax12mr(2),princ_ax11mr(2)],'r','LineWidth',3);
plot([princ_ax21mr(1),princ_ax22mr(1)],[princ_ax22mr(2),princ_ax21mr(2)],'g','LineWidth',3);

title('MRI with principal axes')

figure;
imagesc(head_frozen);
colormap('Gray')
hold on;
plot([princ_ax11ct(1),princ_ax12ct(1)],[princ_ax12ct(2),princ_ax11ct(2)],'r','LineWidth',3);
plot([princ_ax21ct(1),princ_ax22ct(1)],[princ_ax22ct(2),princ_ax21ct(2)],'g','LineWidth',3);

title('CT with principal axes')

%% Transformation

Xtrans=[mean_ct(1),princ_ax11ct(1);mean_ct(2),princ_ax12ct(2)];
% Y2=[prinax_mri_11,prinax_mri_2]
ytran=[mean_mr(1), princ_ax11mr(1); mean_mr(2),princ_ax12mr(2)];

[t_trans,s,R,sp,sigma_FRA]=transform(Xtrans.',ytran.');

% Grid, CT size (we want to end with CT size)
x11 = linspace(1, size(head_frozen,2), size(head_frozen,2));
x22 = linspace(1, size(head_frozen,1), size(head_frozen,1));
[x1, x2] = meshgrid(x11, x22);

x1z = zeros(size(x1));
x2z = zeros(size(x2));
grid_points = [x1(:),x2(:)]';

% Transform with params
% Interpolate intensities, based on grid location
grid_transformed = s*R*grid_points+repmat(sp,1,size(grid_points,2));
% For some reason image is rotated 180
int_ctspace = imrotate(interp2(double(head_mri),grid_transformed(1,:),grid_transformed(2,:),'bilinear'),180);

% Reshape to image size
x1_mr = reshape(grid_transformed(1,:), size(x1));
x2_mr = reshape(grid_transformed(2,:), size(x2));
int_ctspace = reshape(int_ctspace, size(head_frozen));

% Gather images
im_joint = zeros (size(head_frozen,1), size(head_frozen,2), 3);
im_joint(:,:,1) = head_frozen;
im_joint(:,:,2) = int_ctspace;

image(im_joint./255);
title('CT as red, MR as green');

%% Optimizing position and rotation

im_jointCT=im_joint(:,:,2);
sp1=linspace(-10 ,10 ,100);
sp2=sp1;
for t1 =1:100
    for t2=1:100
        sp =[sp1(t1), sp2(t2)];
        head_shift = imtranslate(im_jointCT, sp , 'FillValues', 255);
        hist = histogram2(double(head_frozen(:))', double(head_shift(:)'), [0 256 64; 0 256 64]);
        sumh = sum(hist(:));
        hist = hist/sumh;
        % Using Philip M. Hanna, Wright State University Entropy.m
        mi(t1, t2) = Entropy(hist);
    end
end;
figure;
imagesc(mi); colorbar;
ax = gca;
ax.XAxis.TickLabelFormat = '%,.1f';
%xticks([-10:1:10]);
%% Find minimum
% 
[val, index]= min(mi(:));
[I_row, I_col] = ind2sub(size(mi),index);
sp_min =[-sp1(I_col), -sp2(I_row)];
%sp_min = [20,-30];

im_jointCT1=imtranslate(int_ctspace,sp_min,'FillValues',0) ;

im_joint(:,:,2) = im_jointCT1;

image(im_joint);
figure;
image(im_joint./255);