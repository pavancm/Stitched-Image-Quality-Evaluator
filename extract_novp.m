function im = extract_novp(im)

M     = 500;  % Number of hypotheses for RANSAC.
thr   = 0.1;  % RANSAC threshold.

%global functions
global fitfn resfn degenfn psize numpar
fitfn = 'homography_fit';
resfn = 'homography_res';
degenfn = 'homography_degen';
psize   = 4;
numpar  = 9;

[points, descriptors] = vl_sift(single(rgb2gray(im{1})),'PeakThresh', 0,'edgethresh',500);

for n=2:length(im)
    pointsPrev = points;
    descriptorsPrev = descriptors;
    
    [points, descriptors] = vl_sift(single(rgb2gray(im{n})),'PeakThresh', 0,'edgethresh',500);
    matches   = vl_ubcmatch(descriptorsPrev,descriptors);
    
    %data normalization
    data_orig = [pointsPrev(1:2,matches(1,:)) ; ones(1,size(matches,2)) ; points(1:2,matches(2,:)) ; ones(1,size(matches,2)) ];
    [dat_norm_img1,T1] = normalise2dpts(data_orig(1:3,:));
    [dat_norm_img2,T2] = normalise2dpts(data_orig(4:6,:));
    data_norm = [ dat_norm_img1 ; dat_norm_img2 ];
    
    %RANSAC
    [ ~,res,~,~ ] = multigsSampling(100,data_norm,M,10);
    con = sum(res<=thr);
    [ ~, maxinx ] = max(con);
    inliers = find(res(:,maxinx)<=thr);
    
    max_x = max(data_orig(4,inliers));
    im{n} = im{n}(:,round(max_x+1):end,:);
end