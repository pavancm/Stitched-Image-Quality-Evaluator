function feat = extract_feat_constituent(constituent_im,patch_size)

%function to extract features of PIQE framework for constituent images

for k = 1:length(constituent_im)
    [c_x{k},c_y{k}] = extract_patch(constituent_im{k},patch_size);
end

%Steerable Pyramid Decomposition
Nor = 6;    %number of orientations
scale = 2;  %number of scales
sig = 0.1;  %sigma for patch weighting non-linearity

patch_sel_c_h = [];patch_sel_c_v = [];

disp('Extracting features of Marginal Distribution for constituent images....');
for k = 1:length(constituent_im)
    sz = length(patch_sel_c_h);
    parfor t = 1:length(c_x{k})
        block_im = rgb2gray(constituent_im{k}(c_x{k}(t,1):c_x{k}(t,2),c_y{k}(t,1):c_y{k}(t,2),:));
        
        %Steerable pyramid decomposition
        [pyr_c{k,t},pind_c{k,t}] = buildSFpyr(double(block_im),scale,Nor-1);
        
        %GSM based divisive normalization and fitting marginal distribution using GGD
        marginal_c{sz+t} = marginal_est(pyr_c{k,t},pind_c{k,t},scale,Nor,1,1,3,3,50);
        
        %values for patch weighting
        %Horizontal neighbors
        
        du = zeros(13,13);  %13 quantization levels
        quant_im = floor(double(block_im)/21)+1;
        [glcm_constituent,SI] = graycomatrix(quant_im,'Offset',[0,1],'NumLevels',...
            max(quant_im(:)) - min(quant_im(:)) + 1,'GrayLimits',[min(quant_im(:)) max(quant_im(:))]);
        du(min(quant_im(:)):max(quant_im(:)),min(quant_im(:)):max(quant_im(:))) = glcm_constituent;
        stats = graycoprops(du,'Energy');
        patch_sel_c_h(sz+t) = 1 - stats.Energy;   %weight = 1 - energy
        patch_sel_c_h(sz+t) = non_lin_fn(patch_sel_c_h(sz+t),sig);  %non_linearity
        
        %Vertical Neighbors
        du = zeros(13,13);  %13 quantization levels
        [glcm_constituent,SI] = graycomatrix(quant_im,'Offset',[1,0],'NumLevels',...
            max(quant_im(:)) - min(quant_im(:)) + 1,'GrayLimits',[min(quant_im(:)) max(quant_im(:))]);
        du(min(quant_im(:)):max(quant_im(:)),min(quant_im(:)):max(quant_im(:))) = glcm_constituent;
        stats = graycoprops(du,'Energy');
        patch_sel_c_v(sz+t) = 1 - stats.Energy;   %weight = 1 - energy;
        patch_sel_c_v(sz+t) = non_lin_fn(patch_sel_c_v(sz+t),sig);  %non_linearity
    end
end

%concatenation
sz = length(c_x{1});
for k = 2:length(constituent_im)
    pyr_c(1,sz+1:sz+length(c_x{k})) = pyr_c(k,1:length(c_x{k}));
    pind_c(1,sz+1:sz+length(c_x{k})) = pind_c(k,1:length(c_x{k}));
    sz = size(pyr_c,2);
end
pyr_c(2:end,:) = [];pind_c(2:end,:) = [];

%patch weights
wx = patch_sel_c_h';wy = patch_sel_c_v';
wx = repmat(wx,1,12);wy = repmat(wy,1,12);

%steerable features weighted averaging
wt = wx';
g_shape = [];

parfor p = 1:length(marginal_c)
    if(~isempty(marginal_c{p}))
        g_shape(:,p) = marginal_c{p}(1:12);     %shape parameters of ggd
    else
        g_shape(:,p) = -1*ones(12,1);
    end
end

g_shape = g_shape.*wt;
ft_marginal_c = sum(g_shape,2)/sum(wt(1,:));

%patches with non-zero weights
ix = find(patch_sel_c_h>0);iy = find(patch_sel_c_v>0);

%Bivariste GMM fit
disp('Extracting features of Bivariate Distribution for constituent images....');
sc = scale-1;
for orien = 1:Nor
    op = 2*orien - 1:2*orien;
    nband = (sc-1)*Nor+orien+1;
    parfor t = 1:min(length(ix),length(iy))
        aux_m = pyrBand(pyr_c{ix(t)}, pind_c{ix(t)}, nband);
        aux_h = find_neighbor(aux_m,1);  %Horizontal neighbors
        
        Xh = [aux_h{1}(:) - mean(aux_h{1}(:)), aux_h{2}(:) - mean(aux_h{2}(:))];
        joint_cx(t,op) = gmm_eig(Xh,4,2,35);
        
        aux_m = pyrBand(pyr_c{iy(t)}, pind_c{iy(t)}, nband);
        aux_v = find_neighbor(aux_m,3);  %Vertical neighbors
        
        Xv = [aux_v{1}(:) - mean(aux_v{1}(:)), aux_v{2}(:) - mean(aux_v{2}(:))];
        joint_cy(t,op) = gmm_eig(Xv,4,2,35);
    end
end

it = length(ix) >= length(iy);
switch it
    case 1
        for orien = 1:Nor
            op = 2*orien - 1:2*orien;
            nband = (sc-1)*Nor+orien+1;
            parfor t = length(iy)+1:length(ix)
                aux_m = pyrBand(pyr_c{ix(t)}, pind_c{ix(t)}, nband);
                aux_h = find_neighbor(aux_m,1);  %Horizontal neighbors
                
                Xh = [aux_h{1}(:) - mean(aux_h{1}(:)), aux_h{2}(:) - mean(aux_h{2}(:))];
                joint_cx(t,op) = gmm_eig(Xh,4,2,35);
            end
        end
    case 0
        for orien = 1:Nor
            op = 2*orien - 1:2*orien;
            nband = (sc-1)*Nor+orien+1;
            parfor t = length(ix)+1:length(iy)
                aux_m = pyrBand(pyr_c{iy(t)}, pind_c{iy(t)}, nband);
                aux_v = find_neighbor(aux_m,3);  %Vertical neighbors
                
                Xv = [aux_v{1}(:) - mean(aux_v{1}(:)), aux_v{2}(:) - mean(aux_v{2}(:))];
                joint_cy(t,op) = gmm_eig(Xv,4,2,35);
            end
        end
end

joint_c_avex = zeros(size(wx));joint_c_avey = zeros(size(wy));
joint_c_avex(ix,:) = joint_cx;joint_c_avey(ix,:) = joint_cy;

%square root of eigen values
joint_c_avex = sqrt(joint_c_avex.*(wx.^2));joint_c_avey = sqrt(joint_c_avey.*(wy.^2));
joint_featx = sum(joint_c_avex,1)'/sum(wx(:,1));joint_featy = sum(joint_c_avey,1)'/sum(wy(:,1));
    
feat = [ft_marginal_c;joint_featx;joint_featy];