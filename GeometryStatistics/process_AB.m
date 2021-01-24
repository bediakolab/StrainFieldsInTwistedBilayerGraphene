% process_AB.m

% addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/Dataset26MaddieFit');
% load('DS26_Day2_TEAMI_objectdata.mat');

displacement_field = m4.DSC_fit_storage;
amplitude = displacement_field(:,:,1).^2 + displacement_field(:,:,2).^2;
[ RGB_color_stack, HSV_color_stack ] = getCustomDisplacementColor( displacement_field, [], [], 1, 1 );
S = HSV_color_stack(:,:,2);
V = HSV_color_stack(:,:,3);
mask1 = S > 0.3;
morphfilt = ~bwareaopen(~bwareaopen(bwmorph(bwmorph(bwmorph(bwmorph(mask1,'clean'),'spur',10),'clean'),'close'),100),100);
morphfilt = trimArray(morphfilt,0,mask1);
res = imgaussfilt(double(morphfilt),5);
mask1v2 = res > 0.8;
% mask1v2 = mask1v2 | edge(mask1v2);
mask2 = V < 0.7;
res2 = imgaussfilt(double(mask2),2);
mask2v2 = res2 > 0.6;
mask = mask1v2 | mask2v2;
figure;
imagesc(mask);
connstruct = bwconncomp(~mask);
nnz_storage = zeros(connstruct.NumObjects,1);
for i = 1:connstruct.NumObjects
    nnz_storage(i) = numel(connstruct.PixelIdxList{i});
end
% figure;
% histogram(nnz_storage);
conres = flipud(sortrows(horzcat(nnz_storage,(1:connstruct.NumObjects)')));
diffs = abs(diff(conres(:,1)));
avdiff = movmean(diffs,5);
[~,idx] = max(avdiff);
idx = idx + 10;
connected_regions = conres;
connected_regions(idx+1:end,:) = [];
connected_regions_map = connected_regions(:,2);

% pause(5);
% figure
combo = zeros(size(S));
ABcentroids = zeros(numel(connected_regions_map),2);
[y, x] = ndgrid(1:size(S, 1), 1:size(S, 2));
count = 1;
for j = connected_regions_map'
    thismat = zeros(size(S));
    thismat(connstruct.PixelIdxList{j}) = 1;
%     imagesc(emptymat);
%     pause(1)
    combo = thismat | combo;
    ABcentroids(count,:) = mean([x(logical(thismat)), y(logical(thismat))]);
    count = count + 1;
end

figure;
imagesc(combo);
hold on
scatter(ABcentroids(:,1),ABcentroids(:,2),'filled','r');

[overlay_figh,axh] = m4.makeCustomDisplacementColorPlot([],[],[],[],1);
hold on
scatter(ABcentroids(:,1),ABcentroids(:,2),200,'filled','k');

%% Process the AA regions
figure
imagesc(mask2v2);
set(gca,'yDir','normal');

fAA = imgaussfilt(double(mask2v2),10);
maxvals = imregionalmax(fAA);
lininds = find(maxvals);
[yv,xv] = ind2sub(size(fAA),lininds);
AA_centroids = horzcat(xv,yv);
hold on
scatter(AA_centroids(:,1),AA_centroids(:,2),100,'filled','r');
figure(overlay_figh);
scatter(AA_centroids(:,1),AA_centroids(:,2),50,'filled','g');

%% Getting saddle point connectivities

allowable_mask = bwmorph(bwneighborerode(morphfilt,5,20),'thin',5);
dists1 = bwdistgeodesic(allowable_mask,xv(5)',yv(5)');
dists1(isnan(dists1)) = inf;
dists2 = bwdistgeodesic(allowable_mask,xv(7)',yv(7)');
dists2(isnan(dists2)) = inf;
dists = dists1 + dists2;
dists = round(dists * 8) / 8;

paths = imregionalmin(dists);
line_path = bwmorph(paths,'thin',inf);
figure
imagesc(line_path);

%% New BW morphology method for processing the AA regions
% f1 = bwmorph(bwmorph(bwmorph(mask2,'spur',10),'clean'),'fill',3);
f1 = bwmorph(bwmorph(mask2,'clean'),'fill',3);
f2 = ~bwareaopen(~bwareaopen(f1,50),50,4);  % man
f2AA = imgaussfilt(double(f2),15);
maxvals = imregionalmax(f2AA);
lininds = find(maxvals);
[yv,xv] = ind2sub(size(f2AA),lininds);
AA_centroids_medmorph = horzcat(xv,yv);

f3 = bwmorph(f2,'shrink',inf);
lininds = find(f3);
[yv,xv] = ind2sub(size(fAA),lininds);
AA_centroids_allmorph = horzcat(xv,yv);

figure(overlay_figh);
scatter(AA_centroids_medmorph(:,1),AA_centroids_medmorph(:,2),40,'filled','y');
scatter(AA_centroids_allmorph(:,1),AA_centroids_allmorph(:,2),30,'filled','b');
legend('im','Org gauss filter method','medium morphology processing','all morphology processing')


%% Yet another method

amplitude_for_calc = amplitude;
amplitude_for_calc(~mask2) = 2;
mask2 = V < 0.7;
res3 = imgaussfilt(amplitude_for_calc,4);
minvals = imregionalmin(res3);
linindsn = find(minvals);
[yv,xv] = ind2sub(size(fAA),linindsn);
AA_centroids_useampinfo = horzcat(xv,yv);
figure(overlay_figh);
scatter(AA_centroids_useampinfo(:,1),AA_centroids_useampinfo(:,2),30,'filled','m');

%% Perform AA region connection via only saddle points of a particular class.
% New method 04/02/2020

SP_hsv = [0,1,1;
          0.33,1,1;
          0.66,1,1];
AB_hsv = [0,0,1;
          0.33,0,1;
          0.66,0,1];
% 
% [rgb2,hsv2] = getCustomDisplacementColor( displacement_field, SP_hsv, AB_hsv, 1, 1 );
% H = hsv2(:,:,1);  % contains the saddle point uniqueness info
% % mask1 here has thresholded the saturation to know that we are in an SP
% % now threshold the hue to figure out which SP

[~,AAfitparams] = m4.fitAAToCircle(0);
radius_multiplier = 3;
AA_lininds_mask = m4.getAACircleMask(radius_multiplier);
[fighn] = m4.makeCustomDisplacementColorPlot(SP_hsv,AB_hsv,0.3,9,1);
plot_flag = 1;
upsample_flag = 0;
SP_id = 1;
[ SP_line ] = registerUniqueSPLines( displacement_field, AA_lininds_mask, SP_id, plot_flag, fighn, upsample_flag );
SP_id = 2;
[ SP_line ] = registerUniqueSPLines( displacement_field, AA_lininds_mask, SP_id, plot_flag, fighn, upsample_flag );
SP_id = 3;
[ SP_line ] = registerUniqueSPLines( displacement_field, AA_lininds_mask, SP_id, plot_flag, fighn, upsample_flag );

% SP1mask = (H <= 1/6 | H > 5/6) & mask1; 
% SP2mask = (H <= 3/6 & H > 1/6) & mask1; 
% SP3mask = (H <= 5/6 & H > 3/6) & mask1; 
% SP1morphfilt = bwmorph(SP1mask,'clean');
% SP1morphfilt = ~bwareaopen(~SP1morphfilt,50);  % close in dark spots
% SP1morphfilt = bwareaopen(SP1morphfilt,20);  % clear out light spots
% SP1morphfilt(AA_lininds_mask) = SP1mask(AA_lininds_mask);  % Protect AA regions from the filter
% inclusive_thresh = 5;
% n = 1;
% 
% % Fatten within the AA mask
% niter = 10;
% for i = 1:niter
%     filteres = ~bwneighborerode( ~SP1morphfilt, 5, 1 );
%     SP1morphfilt(AA_lininds_mask) = filteres(AA_lininds_mask);
% end
% 
% 
% SP1morphfilt = ~bwneighborerode( ~SP1morphfilt, inclusive_thresh, n );
% SP1morphfilt = bwneighborerode( SP1morphfilt, inclusive_thresh, n );
% SP1morphfilt = ~bwneighborerode( ~SP1morphfilt, inclusive_thresh, n );
% SP1morphfilt = bwneighborerode( SP1morphfilt, inclusive_thresh, n );
% SP1morphfilt = ~bwneighborerode( ~SP1morphfilt, inclusive_thresh, n );
% SP1morphfilt = bwneighborerode( SP1morphfilt, inclusive_thresh, n );
% SP1morphfilt = ~bwneighborerode( ~SP1morphfilt, inclusive_thresh, n );
% SP1morphfilt = bwneighborerode( SP1morphfilt, inclusive_thresh, n );
% SP1morphfilt = ~bwneighborerode( ~SP1morphfilt, inclusive_thresh, n );
% SP1morphfilt = bwneighborerode( SP1morphfilt, inclusive_thresh, 7 );
% SP1morphfilt = bwmorph( SP1morphfilt, 'thin', inf );
% figure;imagesc(SP1morphfilt);
% [fighn] = m4.makeCustomDisplacementColorPlot(SP_hsv,AB_hsv,0.3,9,1);
% figure(fighn);
% hold on
% h = imagesc(cat(3,SP1morphfilt,SP1morphfilt,zeros(size(SP1morphfilt))));
% h.AlphaData = 0.5;

inclusive_thresh = 5;
n = 5;
% SP1morphfilt = bwneighborerode( SP1morphfilt, inclusive_thresh, n )
% SP1morphfilt = ~bwareaopen(~bwareaopen(bwmorph(bwmorph(bwmorph(bwmorph(SP1mask,'clean'),'spur',10),'clean'),'close'),100),100);
SP2morphfilt = ~bwareaopen(~bwareaopen(bwmorph(bwmorph(bwmorph(bwmorph(SP2mask,'clean'),'spur',10),'clean'),'close'),100),100);
SP3morphfilt = ~bwareaopen(~bwareaopen(bwmorph(bwmorph(bwmorph(bwmorph(SP3mask,'clean'),'spur',10),'clean'),'close'),100),100);

figure;imagesc(SP2morphfilt);
figure;imagesc(SP3morphfilt);



%% Connect all AA and store structural information
% Older, uses full SP network
if false
AA_centroids_to_use = AA_centroids; % For now, just go with the initial gaussian filter.
AA_centroids_sorted = sortrows(AA_centroids_to_use);
% Attempt 1 at SP registration.
% for each AA region

line_paths = cell(0,1);
nnz_linepaths = zeros(size(AA_centroids_sorted,1) - 1,1);
count = 1;
idx_to_use = 6;
AA_1_coords = AA_centroids_sorted(idx_to_use,:);
others = 1:size(AA_centroids_sorted,1);
others(others == idx_to_use) = [];
for i = others
    count = count + 1
    AA_2_coords = AA_centroids_sorted(i,:);
    [ line_path ] = registerSP( AA_1_coords, AA_2_coords, morphfilt );
    nnz_linepaths(count) = nnz(line_path);
    line_paths{count} = line_path;
end

aug = horzcat(nnz_linepaths,(1:numel(nnz_linepaths))');
augsort = sortrows(aug);
sort_map = augsort(:,2);
% TODO
% Need to come up with better ways of handling the image boundary, both
% from the standpoint of computing the mask (currently a mess) and also of
% specifying how many nearest neighbors the given AA region should have, to
% get connectivity.
if false
    for i = 1:numel(line_paths)
        [fighh,~] =m4.makeCustomDisplacementColorPlot([],[],[],[],1);
        figure(fighh);
        hold on
        spy(line_paths{sort_map(i)});
        pause(2);
    end
end
end



