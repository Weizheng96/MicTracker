
[nR,nC,nT] = size(dataStack_aligned);

TLlst = trackletLst;


%%% The tracklet time Look-Up-Table:[stFM, endFM, stCtrPt[r,c], endCtrPt[r,c], stAvgItsty, endAvgItsty]
TLLUT = cell2mat(cellfun(@(x) [x(1,1),x(end,1)], TLlst, 'UniformOutput', false));
temp = cell2mat(cellfun(@(x) resCellLst{x(1,1)}.ctrPt(x(1,2),:), TLlst, 'UniformOutput', false));
temp2 = cell2mat(cellfun(@(x) resCellLst{x(end,1)}.ctrPt(x(end,2),:), TLlst, 'UniformOutput', false));
TLLUT = [TLLUT, temp, temp2];
TLLUT = [TLLUT, cellfun(@(x) resCellLst{x(1,1)}.avgItstyVec(x(1,2)), TLlst), ...
    cellfun(@(x) resCellLst{x(end,1)}.avgItstyVec(x(end,2)), TLlst)];

%%% The starting/ending tracklet in each frame Look-Up-List
TLstEndLUL = cell(nT,1);
for tEnd = 1:nT
    TLstEndLUL{tEnd}.stTLs = find(TLLUT(:,1)==tEnd);
    TLstEndLUL{tEnd}.endTLs = find(TLLUT(:,2)==tEnd);
end



showImgStack = zeros(nRow*2+3,nCol,3,nFM-1);

for iFM1 = 1:nFM-1
    % iFM1 = 103;
    % iFM2 = 104;
    iFM2 = iFM1 + 1;
    
    % cellIdMap1 = cellIdMapStack(:,:,iFM1);
    % cellIdMap2 = cellIdMapStack(:,:,iFM2);
    
    TLOIlist1 = TLlst(TLstEndLUL{iFM1}.endTLs);
    TLOIlist2 = TLlst(TLstEndLUL{iFM2}.stTLs);
    
    showImg1 = zeros(nRow,nCol);
    for iTL = 1:length(TLOIlist1)
        iCnt = TLOIlist1{iTL}(end,2);
        showImg1(resCellLst{iFM1}.PixelIdxList{iCnt}) = iCnt;
    end
    showImg1 = imdilate(showImg1,strel('disk',1)) & ~showImg1;
    % figure; imshow(showImg1);
    % myImOverlay_white(sqrt(dataStack_aligned(:,:,iFM1)),showImg1>0,false,true,0);
    
    
    showImg2 = zeros(nRow,nCol);
    for iTL = 1:length(TLOIlist2)
        iCnt = TLOIlist2{iTL}(1,2);
        showImg2(resCellLst{iFM2}.PixelIdxList{iCnt}) = iCnt;
    end
    showImg2 = imdilate(showImg2,strel('disk',1)) & ~showImg2;
    % figure; imshow(showImg2);
    % myImOverlay_white(sqrt(dataStack_aligned(:,:,iFM2)),showImg2>0,false,true,0);
    
    
    showImgSub1 = myImOverlay_white(sqrt(dataStack_aligned(:,:,iFM1)), showImg1>0, false, false, 0);
    showImgSub2 = myImOverlay_white(sqrt(dataStack_aligned(:,:,iFM2)), showImg2>0, false, false, 0);
    showImg = [showImgSub1; 1+zeros(3,nCol,3); showImgSub2];
    % figure; imshow(showImg);
    
    showImgStack(:,:,:,iFM1) = showImg;
    
end
implay(showImgStack);

% outputFileName = [logPath,'trackletConstruction\',ImName(1:end-4),...
%     '_initialSeg_FG_FPrmd_tracklet_twoObjectiveSimSc_maxOvlp_TLEndsHeads.tif'];
% for kFM = 1:nFM-1
%     imwrite(showImgStack(:,:,:,kFM), outputFileName, 'WriteMode', 'append','Compression','none');
% end



