% function videoRegister3D(dataStack,TI,minSize,isParallel,ImName,slctYrg)

clear all; close all; clc
dataPath = '..\..\Data\Oct072019_upload\Sept062019 Timelapse\';
outPath = '..\Data\';
logPath = '..\logs\';

% isParallel = false;

fmFileNameList = dir([dataPath,'*_s4_t*.tif']);
ImName = fmFileNameList(1).name(1:end-7);

slctYrg = [351, 650];


%%%%%%   Input data
dataInfo = imfinfo([dataPath,fmFileNameList(1).name]);
nY = dataInfo(1).Height;
nX = dataInfo(1).Width;
nZ = length(dataInfo);
nT = length(fmFileNameList);
disp(['Size of raw data: ', num2str([nY,nX,nZ,nT])]);

nY = slctYrg(2) - slctYrg(1) + 1;
disp(['Size of selected field of view: ', num2str([nY,nX,nZ,nT])]);


dataStack = zeros(nY,nX,nZ,nT);
for t = 1:nT
    disp(['Reading... t', num2str(t)]);
    fmName = [ImName,'_t',num2str(t),'.tif'];
    for z = 1:nZ
        tempImg = imread([dataPath,fmName],'Index',z);
        tempImg = tempImg(351:650, :, :);
        tempImg = double(tempImg)/(2^dataInfo(1).BitDepth - 1);
        dataStack(:,:,z,t) = tempImg;
    end
end




%%  Normalize by histogram
binCntrs = (0.01:0.01:1)-0.005;
fmEcdfVec = zeros(nT,100);
fmEpdfVec = zeros(nT,100);
for t = 1:nT
    % t = 102;
    tempImg = dataStack(:,:,:,t);
    fmEpdf = hist(tempImg(:),binCntrs);
    fmEpdf = fmEpdf/sum(fmEpdf);
    fmEcdf = cumsum(fmEpdf);
    % figure; plot(fmEcdf);
    fmEcdfVec(t,:) = fmEcdf;
    fmEpdfVec(t,:) = fmEpdf;
end
% figure; plot(binCntrs, fmEcdfVec'); xlabel('Intensity after histogram normalization'); ylabel('ecdf');
% figure; plot(binCntrs, fmEpdfVec'); xlabel('Intensity after histogram normalization'); ylabel('epdf');
% xlim([0,0.06]);

avgEpdf = mean(fmEpdfVec);
% hold on; plot(binCntrs, avgEpdf','r-','LineWidth',3);

for t = 1:nT
    % t = 102;
    disp(t);
    tempImg = dataStack(:,:,:,t);
    % figure; imshow(tempImg(:,:,11));
    % figure; imshow(max(tempImg,[],3));
    tempImg2 = reshape(tempImg,[nY,nX*nZ]);
    % figure; imshow(tempImg2(:,10*nX+1:11*nX));
    tempImg2 = histeq(tempImg2, avgEpdf);
    tempImg = reshape(tempImg2,[nY,nX,nZ]);
    % figure; imshow(tempImg(:,:,11));
    dataStack(:,:,:,t) = tempImg;
end

minVal = min(dataStack(:));
maxVal = max(dataStack(:));
dataStack = (dataStack - minVal) / (maxVal - minVal);

% imgMaxProj = zeros(nY,nX,nT);
% for t = 1:nT
%     imgMaxProj(:,:,t) = max(dataStack(:,:,:,t),[],3);
% end
% % imgMaxProj = max(dataStack,[],3);
% % figure; imshow(imgMaxProj);
% implay(imgMaxProj);


% %%%%%  Y-Z stack of raw data
% YZimgStk = permute(imgStkOri,[1, 3, 2]);




%% Global motion estimation by MATLAB image registration toolbox
%%%% Try simplest model: average motion of the whole foreground
disp('Global image registration...');

% isParallel = false;

% tic
% [nY,nX,nZ,nT] = size(dataStack);

dxdydzGlbVec = zeros(nT,3);  % Displacement from frame 1
dxdydzGlbVec(1,:) = [0,0,0];

[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 300;

for t = 2:nT
    % t = 91;
    disp(['Global motion estimation... t: ', num2str(t)]);
    
    tic
    % imgFrom = dataStack(:,:,round(nFM/2));
    tRef = floor(t/10)*10;
    if (tRef==t)
        tRef = tRef - 10;
    end
    tRef = max(tRef,1);
    
    imgFrom = dataStack(:,:,:,tRef);
    imgTo = dataStack(:,:,:,t);
    % figure; imshow(max(imgFrom,[],3));
    % figure; imshow(max(permute(imgFrom,[1, 3, 2]),[],3));
    % figure; imshow(max(imgTo,[],3));
    % figure; imshow(max(permute(imgTo,[1, 3, 2]),[],3));

    % showImg = imfuse(max(imgFrom,[],3),max(imgTo,[],3),'falsecolor','ColorChannels',[2 1 2]);
    % showImg = imfuse(max(permute(imgFrom,[1, 3, 2]),[],3),max(permute(imgTo,[1, 3, 2]),[],3),'falsecolor','ColorChannels',[2 1 2]);
    % showImg = imfuse(max(permute(imgFrom,[3, 2, 1]),[],3),max(permute(imgTo,[3, 2, 1]),[],3),'falsecolor','ColorChannels',[2 1 2]);
    % figure; imshow(showImg);
    
    % [moving_reg, R_reg] = imregister(imgFrom,imgTo,'translation',optimizer,metric);
    tform = imregtform(imgFrom,imgTo,'translation',optimizer,metric);
    % moving_reg = imwarp(imgFrom,tform,'OutputView',imref3d(size(imgTo)));
    % showImg = imfuse(max(moving_reg,[],3),max(imgTo,[],3),'falsecolor','ColorChannels',[2 1 2]);
    % showImg = imfuse(max(permute(moving_reg,[1, 3, 2]),[],3),max(permute(imgTo,[1, 3, 2]),[],3),'falsecolor','ColorChannels',[2 1 2]);
    % showImg = imfuse(max(permute(moving_reg,[3, 2, 1]),[],3),max(permute(imgTo,[3, 2, 1]),[],3),'falsecolor','ColorChannels',[2 1 2]);
    % figure; imshow(showImg);
        
    dxdydzGlbVec(t,:) = tform.T(4,1:3) + dxdydzGlbVec(tRef,:);
    
    toc
    
end

% outFileName = [logPath,ImName,'_dxdydzGlbVec.mat'];
% save(outFileName, 'dxdydzGlbVec');

% if (isParallel)
%     delete(poolobj);
% end



%%% Align back to the "central" frame
maxAbsDis = ceil(max(abs(dxdydzGlbVec)));
dataStack_aligned = zeros(nY+maxAbsDis(2)*2, nX+maxAbsDis(1)*2, nZ+maxAbsDis(3)*2, nT);
paddedOriImg = zeros(nY+maxAbsDis(2)*2, nX+maxAbsDis(1)*2, nZ+maxAbsDis(3)*2);
itsctMask = true(size(paddedOriImg));
unionMask = false(size(paddedOriImg));

paddedOriMask = false(size(paddedOriImg));
lxLm = maxAbsDis(1)+1;
hxLm = nX+maxAbsDis(1);
lyLm = maxAbsDis(2)+1;
hyLm = nY+maxAbsDis(2);
lzLm = maxAbsDis(3)+1;
hzLm = nZ+maxAbsDis(3);

paddedOriMask(lyLm:hyLm, lxLm:hxLm, lzLm:hzLm) = true;

tic
for t = 1:nT
    disp(['Aligning... t: ', num2str(t)]);
    paddedOriImg(lyLm:hyLm, lxLm:hxLm, lzLm:hzLm) = dataStack(:,:,:,t);
        
    tform.T(4,1:3) = -dxdydzGlbVec(t,:);
    moving_reg = imwarp(paddedOriImg,tform,'OutputView',imref3d(size(paddedOriImg)));
    mask_aligned = imwarp(paddedOriMask,tform,'OutputView',imref3d(size(paddedOriImg)));
    % showImg = imfuse(max(moving_reg,[],3),max(paddedOriImg,[],3),'falsecolor','ColorChannels',[2 1 2]);
    % showImg = imfuse(max(permute(moving_reg,[1, 3, 2]),[],3),max(permute(imgTo,[1, 3, 2]),[],3),'falsecolor','ColorChannels',[2 1 2]);
    % showImg = imfuse(max(permute(moving_reg,[3, 2, 1]),[],3),max(permute(imgTo,[3, 2, 1]),[],3),'falsecolor','ColorChannels',[2 1 2]);
    % figure; imshow(showImg);
    % figure; imshow(sqrt(max(moving_reg,[],3)));
    
    dataStack_aligned(:, :, :, t) = moving_reg;
    
    unionMask = unionMask | mask_aligned;
    itsctMask = itsctMask & mask_aligned;
end
toc

lyLm = find(max(max(unionMask,[],3),[],2),1,'first');
hyLm = find(max(max(unionMask,[],3),[],2),1,'last');
lxLm = find(max(max(unionMask,[],3),[],1),1,'first');
hxLm = find(max(max(unionMask,[],3),[],1),1,'last');
lzLm = find(max(max(unionMask,[],1),[],2),1,'first');
hzLm = find(max(max(unionMask,[],1),[],2),1,'last');
disp([lyLm,hyLm,lxLm,hxLm,lzLm,hzLm]);

dataStack_aligned = dataStack_aligned(lyLm:hyLm, lxLm:hxLm, lzLm:hzLm, :);

[nYa,nXa,nZa,nTa] = size(dataStack_aligned);

% imgMaxProj = zeros(nYa,nXa,nTa);
% for t = 1:nTa
%     imgMaxProj(:,:,t) = max(dataStack_aligned(:,:,:,t),[],3);
% end
% % imgMaxProj = max(dataStack,[],3);
% % figure; imshow(imgMaxProj);
% implay(sqrt(imgMaxProj));
% 
% % outFileName = [logPath,'alignedAsUnion_',ImName,'_YXmaxProj.tif'];
% % for t = 1:nT
% %     imwrite(sqrt(imgMaxProj(:,:,t)), outFileName, 'WriteMode', 'append','Compression','none');
% % end
% 
% 
% imgMaxProj = zeros(nYa,nZa,nTa);
% for t = 1:nTa
%     imgMaxProj(:,:,t) = max(permute(dataStack_aligned(:,:,:,t),[1,3,2]),[],3);
% end
% % imgMaxProj = max(dataStack,[],3);
% % figure; imshow(sqrt(imgMaxProj(:,:,nTa)));
% implay(sqrt(imgMaxProj));
% 
% % outFileName = [logPath,'alignedAsUnion_',ImName,'_YZmaxProj.tif'];
% % for t = 1:nT
% %     imwrite(sqrt(imgMaxProj(:,:,t)), outFileName, 'WriteMode', 'append','Compression','none');
% % end





%%  Interpolate the z slices and save the processed data
% % % %  (3 times of original resolution, making z comparable to x and y)

% [nYa,nXa,nZa,nTa] = size(dataStack_aligned);

[yYgrid,xXgrid,zZgrid] = ndgrid((1:nYa)', (1:nXa)', (1:nZa)');
zZgrid = (zZgrid-1)*3 + 1;
[qYgrid,qXgrid,qZgrid] = ndgrid((1:nYa)', (1:nXa)', (1:max(zZgrid(:)))');

nZintplt = size(qZgrid,3);

for t = 1:nTa
    disp(['Interpolating z slices... t', num2str(t)]);
    
    imgStkOri = dataStack_aligned(:,:,:,t);
    
    FintrpGrid = griddedInterpolant(yYgrid, xXgrid, zZgrid, dataStack_aligned(:,:,:,t));
    imgStkIntrplt = FintrpGrid(qYgrid, qXgrid, qZgrid);
    % implay(sqrt(imgStkIntrplt));
    
    
    disp(['Saving pre-processed data... t', num2str(t)]);
    outFileName = [outPath,'prePrcssd_',ImName,'_t',num2str(t),'.tif'];
    for z = 1:nZintplt
%         imwrite(imgStkIntrplt(:,:,z), outFileName, 'WriteMode', 'append','Compression','none');
    end
end




%%  Save the pre-processed data as 3D frames
% 
% for t = 1:nTa
%     disp(['Saving pre-processed data... t', num2str(t)]);
%     outFileName = [outPath,'alignedAsUnion_',ImName,'_t',num2str(t),'.tif'];
%     for z = 1:nZa
% %         imwrite(dataStack_aligned(:,:,z,t), outFileName, 'WriteMode', 'append','Compression','none');
%     end
% end







%% %%%%%%%%%% Load the saved aligned data %%%%%%%%%
% dataInfo = imfinfo([logPath,'aligned_',ImName]);
% % dataInfo = imfinfo([logPath,'toyExample2.tif']);
% nY = dataInfo(1).Height;
% nX = dataInfo(1).Width;
% nZ = length(dataInfo);
% % nFM = 100;
% dataStack_aligned = zeros(nY,nX,nZ);
% for kFM = 1:nZ
%     tempImg = imread([logPath,'aligned_',ImName],'Index',kFM);
% %     tempImg = imread([logPath,'toyExample2.tif'],'Index',kFM);
%     tempImg = double(tempImg)/255;
%     dataStack_aligned(:,:,kFM) = tempImg;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


