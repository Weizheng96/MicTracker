function [showImgSub1, showImgSub2] = visualizeFmSegDiff2D(iFM1,iFM2,resCellLst,flagShow,dataStack_aligned,cellIdMapStack)

% iFM1 = 1;
% iFM2 = 2;

img1 = dataStack_aligned(:,:,iFM1);
img2 = dataStack_aligned(:,:,iFM2);

[nRow,nCol] = size(img1);

cellIdMap1 = cellIdMapStack(:,:,iFM1);
cellIdMap2 = cellIdMapStack(:,:,iFM2);

nCnt1 = resCellLst{iFM1}.NumObjects;
nCnt2 = resCellLst{iFM2}.NumObjects;


%%%  Assign each detection in fm1 to the best correspondence in fm2
adjMatFm1toFm2 = false(nCnt1, nCnt2);
for iCnt1 = 1:nCnt1
    % figure; imshow(cellIdMap1 == iCnt1);
    % figure; imshow(cellIdMap2 == iCnt2);
    temp = cellIdMap2(resCellLst{iFM1}.PixelIdxList{iCnt1});
    unqVals = unique(temp);
    
    if length(unqVals) == 1 && unqVals == 0  % % If no overlap, consider distance
        distVec = sqrt(sum((ones(nCnt2,1)*resCellLst{iFM1}.ctrPt(iCnt1,:) - resCellLst{iFM2}.ctrPt).^2,2));
        [~,assgndIcnt2] = min(distVec);
        adjMatFm1toFm2(iCnt1,assgndIcnt2) = true;
        continue;
    end
    
    if length(unqVals) == 1  % % If only one overlap
        assgndIcnt2 = unqVals;
    else
        [ovlpAreas, candIds] = hist(temp, unique(temp));
        ovlpAreas(candIds == 0) = [];
        candIds(candIds == 0) = [];
        [~,tempI] = max(ovlpAreas);
        assgndIcnt2 = candIds(tempI);
    end
    adjMatFm1toFm2(iCnt1,assgndIcnt2) = true;    
end

%%%  Assign each detection in fm2 to the best correspondence in fm1
adjMatFm2fromFm1 = false(nCnt1, nCnt2);
for iCnt2 = 1:nCnt2
    % figure; imshow(cellIdMap1 == iCnt1);
    % figure; imshow(cellIdMap2 == iCnt2);
    temp = cellIdMap1(resCellLst{iFM2}.PixelIdxList{iCnt2});
    unqVals = unique(temp);
    
    if length(unqVals) == 1 && unqVals == 0  % % If no overlap, consider distance
        distVec = sqrt(sum((ones(nCnt1,1)*resCellLst{iFM2}.ctrPt(iCnt2,:) - resCellLst{iFM1}.ctrPt).^2,2));
        [~,assgndIcnt1] = min(distVec);
        adjMatFm2fromFm1(assgndIcnt1,iCnt2) = true;
        continue;
    end
    
    if length(unqVals) == 1  % % If only one overlap
        assgndIcnt1 = unqVals;
    else
        [ovlpAreas, candIds] = hist(temp, unique(temp));
        ovlpAreas(candIds == 0) = [];
        candIds(candIds == 0) = [];
        [~,tempI] = max(ovlpAreas);
        assgndIcnt1 = candIds(tempI);
    end
    adjMatFm2fromFm1(assgndIcnt1,iCnt2) = true;    
end


%%% Show any one-to-multiple cases
hugeDtctInFm1Lbs = sum(adjMatFm2fromFm1,2)>1;

%%% Show any multiple-to-one cases
hugeDtctInFm2Lbs = sum(adjMatFm1toFm2,1)>1;

%%% Check consistency of assignment in the two opposite directions
incnsstDtctInFm1Lbs = false(nCnt1,1);
for iCnt1 = 1:nCnt1
    if ~adjMatFm2fromFm1(iCnt1,adjMatFm1toFm2(iCnt1,:))
        incnsstDtctInFm1Lbs(iCnt1) = true;
    end
end
incnsstDtctInFm2Lbs = false(nCnt2,1);
for iCnt2 = 1:nCnt2
    if ~adjMatFm1toFm2(adjMatFm2fromFm1(:,iCnt2),iCnt2)
        incnsstDtctInFm2Lbs(iCnt2) = true;
    end
end


%%%  Show given detections (one-to-multiple)
showMap1 = zeros(nRow,nCol);
showMap2 = zeros(nRow,nCol);
for iCnt1 = 1:nCnt1
    if hugeDtctInFm1Lbs(iCnt1)
        showMap1(resCellLst{iFM1}.PixelIdxList{iCnt1}) = iCnt1;
        tempIds = find(adjMatFm2fromFm1(iCnt1,:));
        for iii = 1:length(tempIds)
            iCnt2 = tempIds(iii);
            showMap2(resCellLst{iFM2}.PixelIdxList{iCnt2}) = iCnt2;
        end
    end
end
lbmapIdxLst1 = {imdilate(cellIdMap1,strel('disk',1))-cellIdMap1 > 0, imdilate(showMap1,strel('disk',1))-showMap1 > 0};
% myImOverlay_multi_white(sqrt(img1), lbmapIdxLst1, false, true, 0);
lbmapIdxLst2 = {imdilate(cellIdMap2,strel('disk',1))-cellIdMap2 > 0, imdilate(showMap2,strel('disk',1))-showMap2 > 0};
% myImOverlay_multi_white(sqrt(img2), lbmapIdxLst2, false, true, 0);


%%%  Show given detections (multiple-to-one)
showMap1 = zeros(nRow,nCol);
showMap2 = zeros(nRow,nCol);
for iCnt2 = 1:nCnt2
    if hugeDtctInFm2Lbs(iCnt2)
        showMap2(resCellLst{iFM2}.PixelIdxList{iCnt2}) = iCnt2;
        tempIds = find(adjMatFm1toFm2(:,iCnt2));
        for iii = 1:length(tempIds)
            iCnt1 = tempIds(iii);
            showMap1(resCellLst{iFM1}.PixelIdxList{iCnt1}) = iCnt1;
        end
    end
end
lbmapIdxLst1{end+1} = imdilate(showMap1,strel('disk',1))-showMap1 > 0;
lbmapIdxLst2{end+1} = imdilate(showMap2,strel('disk',1))-showMap2 > 0;

showImgSub1 = myImOverlay_multi_white(sqrt(img1), lbmapIdxLst1, false, false, 0);
showImgSub2 = myImOverlay_multi_white(sqrt(img2), lbmapIdxLst2, false, false, 0);

if flagShow
    showImg = [showImgSub1; 1+zeros(3,nCol,3); showImgSub2];
    figure; imshow(showImg);
end




