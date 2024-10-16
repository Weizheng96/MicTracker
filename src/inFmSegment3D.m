function [RoiLst, pcXY3dMapFm] = inFmSegment3D(imgFm, PCsigma, minSz, maxSz)
% %  Intra-frame segmentation (3D) for a single frame
% %
% %  Input:
% %      imgFm      -   3d image, taking values in [0,1]
% %      PCsigma    -   sigma of principle curvature scores
% %      minSz, maxSz   -   The min/max size of expected detections
% %  Output:
% %      RoiLst     -   List of detections


[nnY,nnX,nnZ] = size(imgFm);


%% Add synQuant java path
% Pij = fullfile(synQuantPath, 'src_synquant\ij-1.52i.jar');
% javaaddpath(Pij);
% p1 = fullfile(synQuantPath,'src_synquant\commons-math3-3.6.1.jar');
% javaaddpath(p1);%
% p0 = fullfile(synQuantPath, 'src_synquant\SynQuantVid_v1.2.4.jar');
% javaaddpath(p0);%



%% Get foreground maps

%%%  Local maximum of edge contrast
img = max(imgFm,[],3);
contrastVec = zeros(254,1);
for thr = 1:254
    roiMap = img > thr/255;
    roiOuterEdgeMap = imdilate(roiMap,strel('disk',1)) & ~roiMap;
    % roiInnerEdgeMap = roiMap & ~imerode(roiMap,strel('disk',1));
    contrastVec(thr) = mean(img(roiMap)) - mean(img(roiOuterEdgeMap));
    % contrastVec(thr) = mean(img(roiInnerEdgeMap)) - mean(img(roiOuterEdgeMap));
end
% figure; plot(1:254, contrastVec);

%%% Smooth the contrast curve
windowSz = 3;
b = (1/windowSz)*ones(1,windowSz);
fltdContrastVec = filter(b,1,contrastVec);
% hold on; plot(1:254, fltdContrastVec);

%%% Find threshold: the first peak
slctThr = find((fltdContrastVec >= [0;fltdContrastVec(1:end-1)]) &...
    (fltdContrastVec >= [fltdContrastVec(2:end);0]), 1, 'first');
% myImOverlay_white(sqrt(img),img>slctThr/255,true,true,1);

FGmapStack = imgFm > slctThr/255;
FGmapStack = medfilt3(FGmapStack, [3,3,3]);
FGmapStack = imfill(FGmapStack,'holes');
seR = 3;
SE = false(seR*2+1, seR*2+1, seR*2+1);
SE(seR+1,seR+1,seR+1) = true;
temp = bwdist(SE);
SE = temp <= seR;
FGmapStack = imopen(FGmapStack,SE);
FGmapStack = imclose(FGmapStack,SE);
% implay(FGmapStack);
% myImageStackPrint(FGmapStack(:,:,1:3:end), 10, true);

%%%  Foreground map in Y-Z planes
% YZFGmapStack = permute(FGmapStack,[1, 3, 2]);
% % myImageStackPrint(YZFGmapStack(:,:,1:3:end), 10, true);





%% Get the principle curvature score maps: 3D volume
% tic
% PCsigma = 3;

% pcXY3dMapFm = principalCurvature3d(imgFm, PCsigma, FGmapStack);
% temp=imgFm*255;
% temS=temp(:,:,99);
pcXY3dMapFm = principalCurvature3d(imgFm, PCsigma);
pcXY3dMapFm(1:2,:,:)=0;
pcXY3dMapFm(end-1:end,:,:)=0;
pcXY3dMapFm(:,1:2,:)=0;
pcXY3dMapFm(:,end-1:end,:)=0;

% max(pcXY3dMapFm(:))*255
% temp=pcXY3dMapFm(:,:,99)*255;
% imshow(pcXY3dMapFm(:,:,21)<0);
% myImageStackPrint(1-sqrt(max(0,pcXY3dMapFm(:,:,1:3:nnZ)/max(pcXY3dMapFm(:)))), 10, true);

% tempMapStk3 = permute(pcXY3dMapFm,[1, 3, 2]);
% % implay(1-sqrt(max(0,tempMapStk3/max(tempMapStk3(:)))));

% toc

% [x,y,z]=ind2sub(size(pcXY3dMapFm),find(pcXY3dMapFm==max(pcXY3dMapFm(:))))


%% Get ROIs by thresholding PC score maps

%%% 3D version SynQuant
imgIn = 1 - sqrt(max(0, pcXY3dMapFm/max(pcXY3dMapFm(:))));
imgIn = round(imgIn*255);
% tem=max(pcXY3dMapFm/max(pcXY3dMapFm(:)),0);
% temp=tem(:,:,17);
q.minIntensity = 0;
% minSz = 10;
% maxSz = 10000;
% minSz = 5;
% maxSz = 500;
[~, RoiIdMap, ~] = Synquant4Embryo_Paramater(imgIn, q, minSz, maxSz);
% implay(RoiIdMap>0);
% implay(imgIn==255);
% implay(mat2gray(imgFm))
% addpath('/home/wei/Projects/Kenichi/ImageAnalysis/ImageAnalysis_ForKenichi/NewMethod/boyu');
% outIm = imdisplayWithROI3DMoreRandom(mat2gray(imgFm), RoiIdMap);
% imshow(outIm(:,:,:,21))
% implay(outIm)

RoiIdMap(~FGmapStack) = 0;
% RoiIdMap = bwlabeln(RoiIdMap > 0);


%%%%  Get the ROI list
% tempRoiLst = bwconncomp(RoiIdMap > 0);
% nROI = tempRoiLst.NumObjects;
% RoiLst = cell(nROI,1);
% for iROI = 1:nROI
%     RoiLst{iROI}.voxelIdxLst = tempRoiLst.PixelIdxList{iROI};
%     [yY, xX, zZ] = ind2sub([nnY,nnX,nnZ], RoiLst{iROI}.voxelIdxLst);
%     RoiLst{iROI}.coordinateRgs = [min(yY),max(yY), min(xX),max(xX), min(zZ),max(zZ)];
% end

RoiLst = bwconncomp(RoiIdMap > 0);




%%  Refine ROIs by PC gap
tic

nROI = RoiLst.NumObjects;
imgIn = 1 - sqrt(max(0, pcXY3dMapFm/max(pcXY3dMapFm(:)))); %gk - not needed
% implay(imgIn==1);
% implay(pcXY3dMapFm<=0);
% figure; imshow(RoiIdMap(:,:,58));

RoiLst.coordinateRgs = zeros(RoiLst.NumObjects,6);

for iROI = 1:nROI
    % iROI = 197;
    vxIdxVec = RoiLst.PixelIdxList{iROI};
    [yY, xX, zZ] = ind2sub([nnY,nnX,nnZ],vxIdxVec);
    RoiLst.coordinateRgs(iROI,:) = [min(yY),max(yY), min(xX),max(xX), min(zZ),max(zZ)];
end
    
roiBnMap = false(nnY,nnX,nnZ);
flag2RmVec = false(nROI,1);
for iROI = 1:nROI
    % iROI = 32;
%     disp(['Refine by PC gap... ROI ', num2str(iROI)]);
    
    vxIdxVec = RoiLst.PixelIdxList{iROI};
    yRg = RoiLst.coordinateRgs(iROI,1:2); %gk - remove this
    xRg = RoiLst.coordinateRgs(iROI,3:4);
    zRg = RoiLst.coordinateRgs(iROI,5:6);
    roiBnMap(:) = false;
    roiBnMap(vxIdxVec) = true;

    
    % %%%%  Visualization
    % myImageStackPrint(imgFm(:,:,1:3:end), 10, true);
    % myImageStackPrint(imgIn(:,:,1:3:end), 10, true);
    % myImageStackPrint(roiBnMap(:,:,1:3:end), 10, true);

    % roiBnYZPatch = permute(roiBnMap,[1, 3, 2]);
    % % implay(roiBnYZPatch);
    % % myImageStackPrint(roiBnYZPatch(:,:,1:3:end), 7, true);
    % pcMapYZPatch = permute(imgIn,[1, 3, 2]);
    % % myImageStackPrint(pcMapYZPatch(:,:,1:3:end), 7, true);

    
    
    %%%%  Refinement by princial curvature gap
    tempMap = imgIn;
    tempMap(~roiBnMap) = 0;
    % myImageStackPrint(tempMap, 7, true);
    
    lcMaxBnPatch = imregionalmax(tempMap, 26);
%     implay(lcMaxBnPatch);
%     implay(tempMap==1);
    
    % myImageStackPrint(lcMaxBnPatch, 7, true);
    
    %%% Remove tiny (volume<=10) noisy local maxima
    lcMaxLbPatch = bwlabeln(lcMaxBnPatch);
    % disp(max(lcMaxLbPatch(:)));
    [freqVec, valVec] = hist(lcMaxLbPatch(:),unique(lcMaxLbPatch(:)));
    freqVec(valVec==0) = [];
    valVec(valVec==0) = [];
    toRmLcMaxs = valVec(freqVec <= 10);
    for iLcM = 1:length(toRmLcMaxs)
        lcMaxBnPatch(lcMaxLbPatch == toRmLcMaxs(iLcM)) = false;
    end
    % disp(max(lcMaxLbPatch(:)));

    lcMaxLbPatch = bwlabeln(lcMaxBnPatch);
    
    
    % [showImg, resCmapLb] = myImOverlay3D_multi((lcMaxLbPatch>0)*1, [], lcMaxLbPatch, false, false, []);
    % % [showImg, resCmapLb] = myImOverlay3D_multi(imgIn, [], lcMaxLbPatch, false, true, resCmapLb);
    % % [showImg, resCmapLb] = myImOverlay3D_multi(imgFm, [], lcMaxLbPatch, false, true, resCmapLb);
    % myImageStackPrint(showImg, 7, true);
    
    
    %%%  Reshape the sub-ROIs by watershed
    scMap = 1 - imgIn;
    subRoiLbMap = splitROIbyWtrshd3D(roiBnMap, lcMaxLbPatch, scMap);
    

%     [showImg, resCmapLb] = myImOverlay3D_multi((subRoiLbMap(:,:,1:3:end)>0)*1, [], subRoiLbMap(:,:,1:3:end), false, false, []);
%     % [showImg, resCmapLb] = myImOverlay3D_multi(imgIn(:,:,1:3:end), [], subRoiLbMap(:,:,1:3:end), false, true, resCmapLb);
%     % [showImg, resCmapLb] = myImOverlay3D_multi(imgFm(:,:,1:3:end), [], subRoiLbMap(:,:,1:3:end), false, true, resCmapLb);
%     myImageStackPrint(showImg, 10, true);
%     showImgYZ = permute(showImg,[1, 4, 3, 2]);
%     myImageStackPrint(showImgYZ, 7, true);
    
    
    %%%  Modify the ROI map
    nSubRoi = max(subRoiLbMap(:));
    if nSubRoi > 1
        flag2RmVec(iROI) = true;
        subRoiLbMap(~roiBnMap) = 0;
        roiBnMap = subRoiLbMap > 0;
        tempRois = bwconncomp(roiBnMap);
        RoiLst.PixelIdxList(end+1:end+nSubRoi) = tempRois.PixelIdxList;
        tempVec = zeros(nSubRoi,6);
        for iSubRoi = 1:nSubRoi
            vxIdxVec = tempRois.PixelIdxList{iSubRoi};
            [yY, xX, zZ] = ind2sub([nnY,nnX,nnZ],vxIdxVec);
            tempVec(iSubRoi,:) = [min(yY),max(yY), min(xX),max(xX), min(zZ),max(zZ)];
        end
        RoiLst.coordinateRgs(end+1:end+nSubRoi, :) = tempVec;
    end
    
end
% nnz(flag2RmVec)
tempInd = find(flag2RmVec);
RoiLst.PixelIdxList(tempInd) = [];
RoiLst.coordinateRgs(tempInd,:) = [];
RoiLst.NumObjects = length(RoiLst.PixelIdxList);
 

disp(['By PC gap, ',num2str(nROI),' ROIs are devided into ', num2str(RoiLst.NumObjects)]);

toc

% im=zeros(size(imgIn));
% for i=1:RoiLst.NumObjects
%     im(RoiLst.PixelIdxList{i})=i;
% end
% 
% addpath('/home/wei/Projects/Kenichi/ImageAnalysis/ImageAnalysis_ForKenichi/NewMethod/boyu');
% outIm = imdisplayWithROI3DMoreRandom(mat2gray(imgFm), im);
% implay(outIm)






