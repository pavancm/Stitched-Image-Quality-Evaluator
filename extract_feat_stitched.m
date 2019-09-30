function feat = extract_feat_stitched(stitched_im,patch_size)

%function to extract features of PIQE framework for stitched stitched_images
[s_x,s_y] = extract_patch(stitched_im,patch_size);

%Steerable Pyramid Decomposition
Nor = 6;    %number of orientations
scale = 2;  %number of scales
sig = 0.1;  %sigma for patch weighting non-linearity

patch_sel_s_h = [];patch_sel_s_v = [];

%feature extraction for stitched images
disp('Extracting features of Marginal Distribution for stitched image....');

parfor t = 1:length(s_x)
    block_im = rgb2gray(stitched_im(s_x(t,1):s_x(t,2),s_y(t,1):s_y(t,2),:));
    [pyr_s{t},pind_s{t}] = buildSFpyr(double(block_im),scale,Nor-1);
    
    %GSM based divisive normalization and fitting marginal distribution using GGD
    marginal_s{t} = marginal_est(pyr_s{t},pind_s{t},scale,Nor,1,1,3,3,50);
    
    %values for patch weighting
    %Horizontal neighbors
    
    du = zeros(13,13);  %13 quantization levels
    quant_im = floor(double(block_im)/21)+1;
    [glcm_constituent,SI] = graycomatrix(quant_im,'Offset',[0,1],'NumLevels',...
        max(quant_im(:)) - min(quant_im(:)) + 1,'GrayLimits',[min(quant_im(:)) max(quant_im(:))]);
    du(min(quant_im(:)):max(quant_im(:)),min(quant_im(:)):max(quant_im(:))) = glcm_constituent;
    stats = graycoprops(du,'Energy');       
    patch_sel_s_h(t) = 1 - stats.Energy;        %weight = 1 - energy
    patch_sel_s_h(t) = non_lin_fn(patch_sel_s_h(t),sig);  %non_linearity
    
    %Vertical Neighbors
    du = zeros(13,13);  %13 quantization levels
    [glcm_constituent,SI] = graycomatrix(quant_im,'Offset',[1,0],'NumLevels',...
        max(quant_im(:)) - min(quant_im(:)) + 1,'GrayLimits',[min(quant_im(:)) max(quant_im(:))]);
    du(min(quant_im(:)):max(quant_im(:)),min(quant_im(:)):max(quant_im(:))) = glcm_constituent;
    stats = graycoprops(du,'Energy');      
    patch_sel_s_v(t) = 1 - stats.Energy;        %weight = 1 - energy
    patch_sel_s_v(t) = non_lin_fn(patch_sel_s_v(t),sig);  %non_linearity
end

%patch weights
wx = patch_sel_s_h';wy = patch_sel_s_v';
wx = repmat(wx,1,12);wy = repmat(wy,1,12);

%steerable features weighted averaging
wt = wx';
g_shape = [];

parfor p = 1:length(marginal_s)
    if(~isempty(marginal_s{p}))
        g_shape(:,p) = marginal_s{p}(1:12);     %shape parameters of ggd
    else
        g_shape(:,p) = -1*ones(12,1);
    end
end

g_shape = g_shape.*wt;
ft_marginal_s = sum(g_shape,2)/sum(wt(1,:));

%patches with non-zero weights
ix = find(patch_sel_s_h>0);iy = find(patch_sel_s_v>0);

%Bivariste GMM fit
disp('Extracting features of Bivariate Distribution for stitched image....');
sc = scale - 1;
for orien = 1:Nor
    op = 2*orien - 1:2*orien;
    nband = (sc-1)*Nor+orien+1;
    parfor t = 1:min(length(ix),length(iy))
        aux_m = pyrBand(pyr_s{ix(t)}, pind_s{ix(t)}, nband);
        aux_h = find_neighbor(aux_m,1);  %Horizontal neighbors
        
        Xh = [aux_h{1}(:) - mean(aux_h{1}(:)), aux_h{2}(:) - mean(aux_h{2}(:))];
        joint_sx(t,op) = gmm_eig(Xh,4,2,35);
        
        aux_m = pyrBand(pyr_s{iy(t)}, pind_s{iy(t)}, nband);
        aux_v = find_neighbor(aux_m,3);  %Vertical neighbors
        
        Xv = [aux_v{1}(:) - mean(aux_v{1}(:)), aux_v{2}(:) - mean(aux_v{2}(:))];
        joint_sy(t,op) = gmm_eig(Xv,4,2,35);
    end
end

it = length(ix) >= length(iy);
switch it
    case 1
        for orien = 1:Nor
            op = 2*orien - 1:2*orien;
            nband = (sc-1)*Nor+orien+1;
            parfor t = length(iy)+1:length(ix)
                aux_m = pyrBand(pyr_s{ix(t)}, pind_s{ix(t)}, nband);
                aux_h = find_neighbor(aux_m,1);  %Horizontal neighbors
                
                Xh = [aux_h{1}(:) - mean(aux_h{1}(:)), aux_h{2}(:) - mean(aux_h{2}(:))];
                joint_sx(t,op) = gmm_eig(Xh,4,2,35);
            end
        end
    case 0
        for orien = 1:Nor
            op = 2*orien - 1:2*orien;
            nband = (sc-1)*Nor+orien+1;
            parfor t = length(ix)+1:length(iy)
                aux_m = pyrBand(pyr_s{iy(t)}, pind_s{iy(t)}, nband);
                aux_v = find_neighbor(aux_m,3);  %Vertical neighbors
                
                Xv = [aux_v{1}(:) - mean(aux_v{1}(:)), aux_v{2}(:) - mean(aux_v{2}(:))];
                joint_sy(t,op) = gmm_eig(Xv,4,2,35);
            end
        end
end

joint_s_avex = zeros(size(wx));joint_s_avey = zeros(size(wy));
joint_s_avex(ix,:) = joint_sx;joint_s_avey(ix,:) = joint_sy;

%square root of eigen values
joint_s_avex = sqrt(joint_s_avex.*(wx.^2));joint_s_avey = sqrt(joint_s_avey.*(wy.^2));
joint_featx = sum(joint_s_avex,1)'/sum(wx(:,1));joint_featy = sum(joint_s_avey,1)'/sum(wy(:,1));

feat = [ft_marginal_s;joint_featx;joint_featy];
