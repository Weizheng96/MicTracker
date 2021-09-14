function [zMap, synId, fMap] = Synquant4Embryo_Paramater(imgIn, q, minSz, maxSz, flag_2d)
% cell detection using synQuant for 3D image
% INPUT:
% vid: is YXZ 3D image
% q: minIntensity
% minSz: the minimum size of valid region (default 100 for 3d, 20 for 2d)
% maxSz: the maximum size of valid region
% flag_2d: true means we process image one by one, default is false
% OUTPUT:
% zMap: the zscore map of detected regions
% synId: the id map of the detection regions
% fMap: foreground maps
% contact: ccwant@vt.edu 02/04/2020

if nargin < 5
    flag_2d = false;
end
if nargin < 3
    if flag_2d % TODO: 2d funciton has not been implemented
        minSz = 20;
        maxSz = 450;
    else
        minSz = 100;
        maxSz = 3000;
    end
end

if flag_2d % TODO: 2d funciton has not been implemented
    minfill = 0.2;
    maxWHRatio = 10;
else
    minfill = 0.0001;
    maxWHRatio = 100;
end

fMap =  imgIn > q.minIntensity;
maxVid = max(imgIn(:));
if maxVid>255
    imgIn = round(255 * imgIn / maxVid);
end
q = paraQ3D(1, 0.5); % number of channels (always 1); ratio used for noise estimation, small value==>small variance
% INPUTS of paraP3D
% 1. fdr threshold; 2. z-score threshold; 3. lowest intensity; 4. highest
% intensity; 5. minimum size; 6. largest size; 7. min fill; 8. max WHRatio
p = paraP3D(0, 0,0,255,minSz, maxSz, minfill, maxWHRatio);
vox_x = 2e-7; % useless indeed
% detection for one channel
[H1,W1,D1] = size(imgIn);
datx = zeros(D1,H1*W1,'uint8');
for ii=1:D1
    tmp = imgIn(:,:,ii)';
    datx(ii,:) = tmp(:);
end
det_res = ppsd3D(datx, W1, H1, vox_x, p,q);
zMap1 = det_res.ppsd_main.zMap;
synId1 = det_res.ppsd_main.kMap;

zMap = zeros(size(imgIn));
synId = zeros(size(imgIn));
for i=1:size(zMap,3)
    zMap(:,:,i) = zMap1(i,:,:);
    synId(:,:,i) = synId1(i,:,:);
end
synId(zMap<=1 | ~fMap) = 0;
% implay(mat2gray(zMap))
% temp=zMap(:,:,32);
% imagesc(temp)

% s = regionprops3(synId, {'VoxelIdxList'});
% cnt = 0;
% synId = zeros(size(imgIn));
% for i=1:numel(s.VoxelIdxList)
%     if length(s.VoxelIdxList{i}) >= minSz
%         cnt = cnt + 1;
%         synId(s.VoxelIdxList{i}) = cnt;
%     end
% end

s = regionprops(synId, {'PixelIdxList'});
cnt = 0;
synId = zeros(size(imgIn));
for i=1:numel(s)
    if length(s(i).PixelIdxList) >= minSz
        cnt = cnt + 1;
        synId(s(i).PixelIdxList) = cnt;
    end
end




% overlay_im = zeros([H1,W1,3,D1], 'uint8');
% for i=1:D1
%     overlay_im(:,:,1,i) = uint8((fMap(:,:,i)>0) * 100);
%     overlay_im(:,:,2,i) = uint8(imgIn(:,:,i));
% end
% svf = 'C:\Users\Congchao\Desktop\cell_detection_samples\crop_embryo_data_500x500x30x40\';
% tifwrite(overlay_im, fullfile(svf, 'overlay_synQuant'));