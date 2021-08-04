function showImgStack = showTraces_lines_withSeg(traceSet, dataStack_aligned, cellIdMapStack, isAll, tgtTraceIDs, isImgBkg, isWithSeg)

if ~isAll
    traceSet = traceSet(tgtTraceIDs);
end

nTrace = length(traceSet);
[nR, nC, nT] = size(dataStack_aligned);


%%%%   Map stack of segmentation results
if isWithSeg
    segLbMapStack = false(nR, nC, nT);
    for iFM = 1:nT
        segLbMapStack(:,:,iFM) = imdilate(cellIdMapStack(:,:,iFM),strel('disk',1)) - cellIdMapStack(:,:,iFM) > 0;
    end
    % figure; imshow(segLbMapStack(:,:,116));
end


% %%%%  Map of traces (centers)
% lbMapStack = zeros(nR, nC, nT);
% for iTrc = 1:nTrace
%     dtcts = traceSet{iTrc};
%     dtcts = round(dtcts);
%     dtctPosInds = sub2ind([nR,nC,nT], dtcts(:,1), dtcts(:,2), dtcts(:,3));
%     lbMapStack(dtctPosInds) = iTrc;
% end
% % implay(lbMapStack);


%%%%  Map of traces (lines)
lbMapStack = zeros(nR, nC, nT);
for iTrc = 1:nTrace
    % iTrc = 142;
    subVec = [];
    dtcts = traceSet{iTrc};
    dtcts = round(dtcts(:,1:3));
    lbMapStack(dtcts(1,1), dtcts(1,2), dtcts(1,3)) = iTrc;
    if size(dtcts,1) > 1
        for iiT = 2:size(dtcts,1)
            r1 = dtcts(iiT-1,1);
            c1 = dtcts(iiT-1,2);
            r2 = dtcts(iiT,  1);
            c2 = dtcts(iiT,  2);
            if r1==r2 && c1==c2
                rVec = r1;
                cVec = c1;
            elseif abs(r1-r2) >= abs(c1-c2)
                slope = (c2-c1)/(r2-r1);
                rVec = (min(r1,r2):max(r1,r2))';
                cVec = round(c1 + slope*(rVec - r1));
            else
                slope = (r2-r1)/(c2-c1);
                cVec = (min(c1,c2):max(c1,c2))';
                rVec = round(r1 + slope*(cVec - c1));
            end
            subVec = [subVec; [rVec, cVec]];
            stckIdxVec = sub2ind([nR,nC,nT], subVec(:,1), subVec(:,2), dtcts(iiT,3)+zeros(size(subVec,1),1));
            lbMapStack(stckIdxVec) = iTrc;
        end
    end
end
% implay(lbMapStack);


%%%%   Image Stack of Result Presentation
showImgStack = zeros(nR, nC, 3, nT);
cmap = hsv(nTrace);
[~,temp] = sort(rand(size(cmap,1),1));
cmap = cmap(temp,:);
cmapSeg = [1,1,1] * 0.8;

for iFM = 1:nT
    lbMap = lbMapStack(:,:,iFM);
    %     lbMap = imdilate(lbMap, [1,1,1; 1,1,1; 1,1,1]);
    if isImgBkg
        temp = sqrt(dataStack_aligned(:,:,iFM)) * 0.8;
    else
        temp = zeros(nR, nC) + 0.5;
    end
    lbs = lbMap(lbMap>0);
    lbMap = lbMap>0;
    for cl = 1:3
        if isWithSeg
            temp(segLbMapStack(:,:,iFM)) = cmapSeg(cl);
        end
        temp(lbMap) = cmap(lbs,cl);
        showImgStack(:,:,cl,iFM) = temp;
    end
end
% implay(showImgStack);
% figure; imshow(showImgStack(:,:,:,nFM));





