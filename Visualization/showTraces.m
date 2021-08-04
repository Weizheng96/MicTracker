function [showImgStack,lbMapStack] = showTraces(traceSet, dataStack_aligned, isAll, tgtTraceIDs)

if ~isAll
    traceSet = traceSet(tgtTraceIDs);
end

nTrace = length(traceSet);
[nR, nC, nT] = size(dataStack_aligned);

lbMapStack = zeros(nR, nC, nT);
for iTrc = 1:nTrace
    dtcts = traceSet{iTrc};
    dtcts = round(dtcts);
    dtctPosInds = sub2ind([nR,nC,nT], dtcts(:,1), dtcts(:,2), dtcts(:,3));
    lbMapStack(dtctPosInds) = iTrc;
end
% implay(lbMapStack);

showImgStack = zeros(nR, nC, 3, nT);
cmap = hsv(nTrace);
[~,temp] = sort(rand(size(cmap,1),1));
cmap = cmap(temp,:);
for iFM = 1:nT
    lbMap = lbMapStack(:,:,iFM);
    lbMap = imdilate(lbMap, [1,1,1; 1,1,1; 1,1,1]);
    temp = sqrt(dataStack_aligned(:,:,iFM));
    lbs = lbMap(lbMap>0);
    lbMap = lbMap>0;
    for cl = 1:3
        temp(lbMap) = cmap(lbs,cl);
        showImgStack(:,:,cl,iFM) = temp;
    end
end
% implay(showImgStack);






