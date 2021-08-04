function [cellLst,nullMtScPhat,altMtScPhat,apprPara,posParaStatic,effZinFmTwoEndIdx] = paraLearn(cellLst, imgStack)
% %  Learn the distribution parameters of features
% %  Input: 
% %     cellLst     -   List of detections
% %     imgStack    -   Raw image stack (3d+t)
% %  Ouput:
% %     dtctMtScVec     -   Vector of rapid-motion score of the detections
% %     nullMtScPhat    -   Parameters of rapid-motion score (null/non-rapid)
% %     altMtScPhat     -   Parameters of rapid-motion score (alternative/rapid)
% %     apprPara        -   Parameters of appearance features
% %     posParaStatic   -   Parameters of position features


disp('Learning the parameters...');

tic

[nnY, nnX, nnZ, nnT] = size(imgStack);


%%  Prepare: index LUTs

%%% The overall cell representer indices
cellIdxsLst = cellfun(@(x) (1:x.NumObjects)', cellLst, 'UniformOutput', false);
for t = 2:nnT
    cellIdxsLst{t}(:,1) = cellIdxsLst{t}(:,1) + cellIdxsLst{t-1}(end,1);
end



%%  Rapid motion: regional appearance change

%%%%%   Check the ineffective z-stacks
effZinFmTwoEndIdx = zeros(nnT, 2);
for tt = 1:nnT
    for zz = 1:3:nnZ
        imgPatchZ0 = imgStack(:,:,zz,tt);
        if var(imgPatchZ0(:)) ~= 0
            effZinFmTwoEndIdx(tt, 1) = zz;
            break;
        end
    end
    for zz = nnZ:-3:1
        imgPatchZ0 = imgStack(:,:,zz,tt);
        % inEffZlbs(zz) = var(imgPatchZ0(:)) == 0;
        if var(imgPatchZ0(:)) ~= 0
            effZinFmTwoEndIdx(tt, 2) = zz;
            break;
        end
    end
end



%%%%%   Get score for any detection
% dtctMtScVec = nan(nDtct, 2);
cellLst{1}.dtctMtScVec = nan(cellLst{1}.NumObjects, 2);

for tt = 1:(nnT-1)
    img1 = imgStack(:,:,:,tt);
    img2 = imgStack(:,:,:,tt+1);
    
    cellLst{tt+1}.dtctMtScVec = nan(cellLst{tt+1}.NumObjects, 2);
    
    effVxBnds1 = [(effZinFmTwoEndIdx(tt,1)-1)*nnY*nnX + 1, effZinFmTwoEndIdx(tt,2)*nnY*nnX];
    % disp(var(img1([1:(effVxBnds1(1)-1), (effVxBnds1(2)+1):end])));
    effVxBnds2 = [(effZinFmTwoEndIdx(tt+1,1)-1)*nnY*nnX + 1, effZinFmTwoEndIdx(tt+1,2)*nnY*nnX];
    
    for iSub = 1:cellLst{tt}.NumObjects
        vxVec = cellLst{tt}.VoxelIdxList{iSub};
        temp = vxVec >= effVxBnds2(1) & vxVec <= effVxBnds2(2);
        if nnz(temp) >= 0.5*length(vxVec)
            vxVec = vxVec(temp);
            % dtctMtScVec(cellIdxsLst{tt}(iSub), 2) = sqrt(mean((img1(vxVec)-img2(vxVec)).^2));
            cellLst{tt}.dtctMtScVec(iSub, 2) = sqrt(mean((img1(vxVec)-img2(vxVec)).^2));
        end
    end
    
    for iSub = 1:cellLst{tt+1}.NumObjects
        vxVec = cellLst{tt+1}.VoxelIdxList{iSub};
        temp = vxVec >= effVxBnds1(1) & vxVec <= effVxBnds1(2);
        if nnz(temp) >= 0.5*length(vxVec)
            vxVec = vxVec(temp);
            % dtctMtScVec(cellIdxsLst{tt+1}(iSub), 1) = sqrt(mean((img1(vxVec)-img2(vxVec)).^2));
            cellLst{tt+1}.dtctMtScVec(iSub, 1) = sqrt(mean((img1(vxVec)-img2(vxVec)).^2));
        end
    end
    
end

dtctMtScAllVec = cell2mat(cellfun(@(x) x.dtctMtScVec, cellLst, 'UniformOutput', false));

%%%%%   Fit null probability distribution model to the empirical distribution
data4Null = dtctMtScAllVec(~isnan(dtctMtScAllVec));
% [mu, sigma] = fitTruncGauss(data4Null);
% temp = randn(10000,1)*sigma + mu;
% mtScPhat = expfit(1-data4Null);
% temp = exprnd(mtScPhat,10000,1);
nullMtScPhat = gamfit(data4Null + 1e-10*rand(length(data4Null),1));
% temp = gamrnd(nullMtScPhat(1),nullMtScPhat(2),10000,1);

% figure; histogram(data4Null,100); xlim([0,1]);
% figure; histogram(temp,200); xlim([0,1]);



%%%%%   Fit alernative probability distribution model to the empirical distribution
data4alt = cell2mat(cellfun(@(x) x.avgItstyVec, cellLst, 'UniformOutput', false));
altMtScPhat = gamfit(data4alt + 1e-10*rand(length(data4alt),1));

% figure; histogram(data4Null,'BinWidth',0.01,'Normalization','probability'); xlim([0,1]);
% hold on; histogram(data4alt,'BinWidth',0.01,'Normalization','probability'); xlim([0,1]);
% figure; histogram(gamrnd(nullMtScPhat(1),nullMtScPhat(2),10000,1),'BinWidth',0.01,'Normalization','probability'); xlim([0,1]);
% hold on; histogram(gamrnd(altMtScPhat(1),altMtScPhat(2),10000,1),'BinWidth',0.01,'Normalization','probability'); xlim([0,1]);
% temp = 0:0.001:1;
% figure; plot(temp,gampdf(temp,nullMtScPhat(1),nullMtScPhat(2)),'LineWidth',2);
% hold on; plot(temp,gampdf(temp,altMtScPhat(1),altMtScPhat(2)),'LineWidth',2);








%%  Features for linking

%%%%%   Get the properties of spatially most overlapped detection pairs
clsstOvrlpRate = cell(nnT-1,1);
clsst2ndNbIncldRate = cell(nnT-1,1);
clsstIntDiff = cell(nnT-1,1);
clsstDsplcmt = cell(nnT-1,1);
clsstSizeChange = cell(nnT-1,1);
clsstSizeChangeRate = cell(nnT-1,1);

cellIdMap2 = zeros(nnY, nnX, nnZ);

for t = 1:nnT-1
    % t = 52;
    % disp(['iFM: ', num2str(t)]);
    % cellIdMap1 = cellIdMapStack(:,:,iFM);
    % cellIdMap2 = cellIdMapStack(:,:,t + 1);
    cellIdMap2(:) = 0;
    for iCnn = 1:cellLst{t+1}.NumObjects
        cellIdMap2(cellLst{t+1}.VoxelIdxList{iCnn}) = iCnn;
    end
    % figure; imshow(max(cellIdMap2,[],3));
    
    nCnt = cellLst{t}.NumObjects;
    ovrlpRate = zeros(nCnt,1);
    nbIncldRate = zeros(nCnt,1);
    intDiff = zeros(nCnt,1);
    dsplsmt = zeros(nCnt,3);
    sizeChange = zeros(nCnt,1);
    sizeChangeRate = zeros(nCnt,1);
    
    for iCnt0 = 1:nCnt
        % iCnt = 50;
        % figure; imshow(max(cellIdMap1 == iCnt,[],3));
        % figure; imshow(max(cellIdMap2 == 89,[],3));
        area1 = cellLst{t}.areaVec(iCnt0);
        temp = cellIdMap2(cellLst{t}.VoxelIdxList{iCnt0});
        unqVals = unique(temp);
        if length(unqVals)==1
            if unqVals == 0
                continue;
            else
                candIds = unqVals;
                ovlpAreas = area1;
            end
        else
            [ovlpAreas, candIds] = hist(temp, unique(temp));
        end
        ovlpAreas(candIds == 0) = [];
        candIds(candIds == 0) = [];
        
        candOvrlpRates = zeros(size(candIds));
        candIncldRates = zeros(size(candIds));
        for iCand = 1:length(candIds)
            area2 = cellLst{t+1}.areaVec(candIds(iCand));
            % candOvrlpRates(iCand) = ovlpAreas(iCand)/max(area1, area2);
            candOvrlpRates(iCand) = ovlpAreas(iCand)/(area1 + area2 - ovlpAreas(iCand));
            candIncldRates(iCand) = ovlpAreas(iCand)/area2;
        end
        [ovrlpRate(iCnt0), idx] = max(candOvrlpRates);
        candIncldRates(idx) = 0;
        nbIncldRate(iCnt0) = max(candIncldRates);
        clstNbId = candIds(idx);
        
        intDiff(iCnt0) = cellLst{t+1}.avgItstyVec(clstNbId) - cellLst{t}.avgItstyVec(iCnt0);        
        dsplsmt(iCnt0,:) = cellLst{t+1}.ctrPt(clstNbId,:) - cellLst{t}.ctrPt(iCnt0,:);        
        vlm1 = cellLst{t}.areaVec(iCnt0);
        vlm2 = cellLst{t+1}.areaVec(clstNbId);
        sizeChange(iCnt0) = vlm2 - vlm1;
        sizeChangeRate(iCnt0) = (vlm2 - vlm1) / ((vlm2 + vlm1)/2);        
        
    end
 
    clsstOvrlpRate{t} = ovrlpRate;
    clsst2ndNbIncldRate{t} = nbIncldRate;
    clsstIntDiff{t} = intDiff;
    clsstDsplcmt{t} = dsplsmt;
    clsstSizeChange{t} = sizeChange;
    clsstSizeChangeRate{t} = sizeChangeRate;
end
clsstOvrlpRate = cell2mat(clsstOvrlpRate);
clsst2ndNbIncldRate = cell2mat(clsst2ndNbIncldRate);
clsstIntDiff = cell2mat(clsstIntDiff);
clsstDsplcmt = cell2mat(clsstDsplcmt);
clsstSizeChange = cell2mat(clsstSizeChange);
clsstSizeChangeRate = cell2mat(clsstSizeChangeRate);

% % xlim([-20,20]); ylim([-20,20]); colormap jet;
% % figure; histogram(clsstIntDiff,100); xlim([-1,1]);
% % figure; histogram(clsstDsplcmt(:,3),100); xlim([-50,50]);
% % figure; histogram(sqrt(sum(clsstDsplcmt.^2,2)),100); xlim([0,50]);
% % figure; histogram(clsstOvrlpRate(clsstOvrlpRate>0),100);
% % figure; histogram(1-clsstOvrlpRate(clsstOvrlpRate>0),100);
% % figure; histogram(clsst2ndNbIncldRate(clsst2ndNbIncldRate>0),100);
% % figure; histogram(clsstSizeChange,400); xlim([-4000,4000]);
% % figure; histogram(clsstSizeChangeRate(:),500); xlim([-2,2]);




% % %  Estimate P(aj | ai, dt) and P(xj | xi, dt, non-OPC)
% [apprPara.intDiff_mu, apprPara.intDiff_sigma] = fitTruncGauss(clsstIntDiff);
[apprPara.intDiff_mu, apprPara.intDiff_sigma] = normfit(clsstIntDiff);
% figure; histogram(randn(10000,1)*apprPara.intDiff_sigma+apprPara.intDiff_mu, 100); xlim([-1,1]);
% [apprPara.sizeChange_mu, apprPara.sizeChange_sigma] = fitTruncGauss(clsstSizeChange);
[apprPara.sizeChange_mu, apprPara.sizeChange_sigma] = normfit(clsstSizeChange);
% figure; histogram(randn(10000,1)*apprPara.sizeChange_sigma+apprPara.sizeChange_mu, 100); xlim([-4000,4000]);
% [apprPara.sizeChangeRt_mu, apprPara.sizeChangeRt_sigma] = fitTruncGauss(clsstSizeChangeRate);
[apprPara.sizeChangeRt_mu, apprPara.sizeChangeRt_sigma] = normfit(clsstSizeChangeRate);
% figure; histogram(randn(10000,1)*apprPara.sizeChangeRt_sigma+apprPara.sizeChangeRt_mu, 100); xlim([-2,2]);


posParaStatic.DsplcmtMu = [0,0];
posParaStatic.DsplcmtSigma = zeros(3,3);
[posParaStatic.DsplcmtMu(1), posParaStatic.DsplcmtSigma(1,1)] = fitTruncGauss(clsstDsplcmt(:,1));
[posParaStatic.DsplcmtMu(2), posParaStatic.DsplcmtSigma(2,2)] = fitTruncGauss(clsstDsplcmt(:,2));
[posParaStatic.DsplcmtMu(3), posParaStatic.DsplcmtSigma(3,3)] = fitTruncGauss(clsstDsplcmt(:,3));
posParaStatic.DsplcmtSigma = posParaStatic.DsplcmtSigma.^2;
% temp = [randn(10000,1)*sqrt(posParaStatic.DsplcmtSigma(1,1))*3 + posParaStatic.DsplcmtMu(1),...
%     randn(10000,1)*sqrt(posParaStatic.DsplcmtSigma(2,2))*3 + posParaStatic.DsplcmtMu(2),...
%     randn(10000,1)*sqrt(posParaStatic.DsplcmtSigma(3,3))*3 + posParaStatic.DsplcmtMu(3)];
% figure; histogram(temp(:,1),100); xlim([-50,50]);
% figure; histogram(sqrt(sum(temp.^2,2)),100); xlim([0,50]);

clsstDist = sqrt(sum(clsstDsplcmt.^2,2));
posParaStatic.distPhat = gamfit(clsstDist(clsstDist>0));
% temp = gamrnd(posParaStatic.distPhat(1),posParaStatic.distPhat(2),10000,1);
% figure; histogram(temp(:,1),100); xlim([0,50]);



% [posParaStatic.ovrlpRtMu, posParaStatic.ovrlpRtSigma] = fitTruncGauss(clsstOvrlpRate(clsstOvrlpRate>0));
% % figure; histogram(randn(10000,1)*posParaStatic.ovrlpRtSigma + posParaStatic.ovrlpRtMu, 100); xlim([0,1]);
% posParaStatic.invOvrlpRtPhat = gamfit(1-clsstOvrlpRate(clsstOvrlpRate>0));
% % temp = gamrnd(posParaStatic.invOvrlpRtPhat(1),posParaStatic.invOvrlpRtPhat(2),10000,1);
posParaStatic.invOvrlpRtPhat = expfit(1-clsstOvrlpRate(clsstOvrlpRate>0));
% temp = exprnd(posParaStatic.invOvrlpRtPhat,10000,1);
% figure; histogram(temp,200); xlim([0,1]);

posParaStatic.nbIncldRtPhat = expfit(clsst2ndNbIncldRate(clsst2ndNbIncldRate>0));
% posParaStatic.nbIncldRtPhat = expfit(clsst2ndNbIncldRate);
% temp = exprnd(posParaStatic.nbIncldRtPhat,10000,1);
% figure; histogram(temp,200); xlim([0,1]);
% disp(exppdf(0, posParaStatic.nbIncldRtPhat));



% outputFileName = [outPath, videoName, '_learnedFeatures.mat'];
% save(outputFileName, 'apprPara', 'posParaStatic', 'mtScPhat', 'dtctMtScVec');

toc



