
clear all; close all; clc
dataPath = '/home/gaurang/My Projects/micTracker/OPCtracking_software_3D/Data/';
logPath = '/home/gaurang/My Projects/micTracker/OPCtracking_software_3D/logs/';
outPath = '/home/gaurang/My Projects/micTracker/OPCtracking_software_3D/Results/';

% synQuantPath = '..\..\Seg3D_4Dissertation\src_code_cellSegment_cc\';
synQuantPath = './';

%%   Pre-processing
% % %  Run "preProcessing.m"



%%   Setting

%%%%    Image name: prefix of all files
fmFileNameList = dir([dataPath,'*_s4_t*.tif']);
ImName = fmFileNameList(1).name(11:end-7);


%%%%    Image size
dataInfo = imfinfo([dataPath,fmFileNameList(1).name]);
nYraw = dataInfo(1).Height;
nXraw = dataInfo(1).Width;
nZraw = length(dataInfo);
nTraw = length(fmFileNameList);
disp(['Size of raw data: ', num2str([nYraw,nXraw,nZraw,nTraw])]);


%%%%    Select field of view of interest
% yRg = [156,255];
yRg = [106,255];
xRg = [151,300];
% xRg = [696,795];
zRg = [1,nZraw];
tRg =  [61,80];

nY = yRg(2) - yRg(1) + 1;
nX = xRg(2) - xRg(1) + 1;
nZ = zRg(2) - zRg(1) + 1;
nT = tRg(2) - tRg(1) + 1;



%%%%    Parallel computing
% if isempty(gcp('nocreate'))
%     poolobj = parpool(15);
% end
% delete(poolobj);




%%   Segmentation I: per-frame

%%%%  Parameters
PCsigma = 3;
minSz = 50;
maxSz = 10000;


%%%%  Add synQuant java path
Pij = fullfile(synQuantPath, 'src_synquant/ij-1.52i.jar');
javaaddpath(Pij);
p1 = fullfile(synQuantPath,'src_synquant/commons-math3-3.6.1.jar');
javaaddpath(p1);%
p0 = fullfile(synQuantPath, 'src_synquant/SynQuantVid_v1.2.4.jar');
javaaddpath(p0);%


%%%%  Per-frame
cellLst = cell(nT,1);
imgStack = zeros(nY, nX, nZ, nT);
PCmapStack = zeros(nY, nX, nZ, nT);

imgFmRaw = zeros(nYraw,nXraw,nZraw);

for t = tRg(1):tRg(2)
    % t = 61;
    tInRg = t - tRg(1) + 1;
    
    disp(['Per-frame segmentation... t', num2str(t)]);
    fmName = ['prePrcssd_',ImName,'_t',num2str(t),'.tif'];
    for z = 1:nZ
        imgFmRaw(:,:,z) = imread([dataPath,fmName], 'Index', z);
    end
    imgFmRaw = imgFmRaw/(2^dataInfo(1).BitDepth - 1); %gk - refactor
    
    imgFm = imgFmRaw(yRg(1):yRg(2), xRg(1):xRg(2), zRg(1):zRg(2));
    imgStack(:,:,:,tInRg) = imgFm;
    % implay(sqrt(imgStack(:,:,:,tInRg)));
    % myImageStackPrint(sqrt(imgStack(:,:,1:3:end,tInRg)), 10, true);
    
    
    % [RoiLst, pcXY3dMapFm] = inFmSegment3D(imgStack(:,:,:,t), PCsigma, minSz, maxSz);
    [RoiLst, PCmapStack(:,:,:,tInRg)] = inFmSegment3D(imgFm, PCsigma, minSz, maxSz);
    
    % myImageStackPrint(1 - sqrt(max(0, pcXY3dMapFm(:,:,1:3:end)/max(pcXY3dMapFm(:)))), 10, true);
    
    
    % % %  Record the detections
    cellLst{tInRg}.NumObjects = RoiLst.NumObjects;
    cellLst{tInRg}.VoxelIdxList = RoiLst.PixelIdxList';
    
    % % %  Get average intensity and area of each cell
    cellLst{tInRg}.avgItstyVec = cellfun(@(x) mean(imgFm(x)), cellLst{tInRg}.VoxelIdxList);
    cellLst{tInRg}.areaVec = cellfun(@length, cellLst{tInRg}.VoxelIdxList);
    
    % % %  Get center position of each cell
    cellLst{tInRg}.ctrPt = zeros(cellLst{tInRg}.NumObjects, 3);
    cellLst{tInRg}.CoordinateRgs = zeros(cellLst{tInRg}.NumObjects, 6);
    for iCell = 1:cellLst{tInRg}.NumObjects
        [yY, xX, zZ] = ind2sub([nY,nX,nZ], cellLst{tInRg}.VoxelIdxList{iCell});
        cellLst{tInRg}.ctrPt(iCell,:) = mean([yY, xX, zZ],1);
        cellLst{tInRg}.CoordinateRgs(iCell,:) = [min(yY),max(yY), min(xX),max(xX), min(zZ),max(zZ)];
    end
end
% outputFileName = [logPath,ImName,'_inFmSegRes.mat'];
% imgRgs = [yRg; xRg; zRg; tRg];
% save(outputFileName, 'cellLst', 'PCmapStack', 'imgRgs');

invPCmapStack = 1 - sqrt(max(0, PCmapStack/max(PCmapStack(:))));
% myImageStackPrint(PCmapStack(:,:,1:3:end,tInRg), 10, true);
% myImageStackPrint(invPCmapStack(:,:,1:3:end,tInRg), 10, true);




%% Visualize Data 
%%%%%%%%%%%%%%%%%%%%%%%%%      Visualize Data      %%%%%%%%%%%%%%%%%%%%%%%%
% imgFm = imgStack(:,:,:,2);
% 
% outputFileName = [logPath,ImName(1:end-4),'_sqrtImg.tif'];
% for z = 1:nZ
%     imwrite(sqrt(imgFm(:,:,z)), outputFileName, 'WriteMode', 'append','Compression','none');
% end
% 
% %%%  2D max projection to Y-X plane
% imgMaxProj = max(imgFm,[],3);
% % figure; imshow(imgMaxProj);
% 
% %%%  Y-Z stack of raw data
% YZimgFm = permute(imgFm,[1, 3, 2]);
% % implay(sqrt(YZimgFm));
% % figure; imshow(YZimgFm(:,:,88));

%%%%%   Given detection(s)
nDtct = sum(cellfun(@(x) x.NumObjects, cellLst));
dtctOIs = 1:nDtct; %%% gk - refactor

%%%%%  Show the detections
traceVxLstLst = cell(length(dtctOIs), 1);
%
cellSubIdxsLUT = cell(nT,1);
for t = 1:nT
    cellSubIdxsLUT{t} = [t + zeros(cellLst{t}.NumObjects,1), (1:cellLst{t}.NumObjects)'];
end
cellSubIdxsLUT = cell2mat(cellSubIdxsLUT);
%
for ii = 1:length(dtctOIs) %%%gk refactor
    traceVxLstLst{ii} = cell(nT,1);
    iDtct = dtctOIs(ii);
    tt = cellSubIdxsLUT(iDtct,1);
    iSub = cellSubIdxsLUT(iDtct,2);
    traceVxLstLst{ii}{tt} = cellLst{tt}.VoxelIdxList{iSub};
end

imgSzNrg.size = [nY,nX,nZ,nT];
imgSzNrg.range = [1,nY; 1,nX; 1,nZ; 1,nT];

% [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, [], false, []);
% % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, (cellIdPatch>0)*1, false, []);
[printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, false, []);
% [printImgStack,printImgYZStack,~] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCpatchStack, true, cmapLb);
implay(printImgStack);
implay(printImgYZStack);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%  Parameter learning for key features
[cellLst,nullMtScPhat,altMtScPhat,apprPara, posParaStatic,effZinFmTwoEndIdx] = paraLearn(cellLst, imgStack);


    

%%  Segmentation II: Refinement by temporal consistency
cellLst_refined = segRefineByTempCons(cellLst, imgStack, invPCmapStack, nullMtScPhat, effZinFmTwoEndIdx);

% outputFileName = [logPath,ImName,'_tempConsRefinedSegRes.mat'];
% imgRgs = [yRg; xRg; zRg; tRg];
% save(outputFileName, 'cellLst_refined', 'imgRgs');







%%   Linking by tow-model MCF

cmpltOpcTrcLst = linkByIntegrtedMCF(cellLst_refined, imgStack, nullMtScPhat, altMtScPhat, apprPara, posParaStatic);



%% Visualization

%%%%%%%%%%%%%%%%%%%%       Visualization       %%%%%%%%%%%%%%%%%%%%

%%%%  Show motion-broken static traces with the significant detections
traceVxLstLst = cellfun(@(x) x.traceVxLst, cmpltOpcTrcLst, 'UniformOutput', false);

% [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, [], false, []);
% % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, (cellIdPatchStack*1, false, []);
% [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, false, []);
[printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, imgStack, false, []);
% [printImgStack,printImgYZStack,~] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, true, cmapLb);
implay(printImgStack);
implay(printImgYZStack);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




