function [lbMapStack,showImgStack] = visualizeTrace2D(tracePxLstLst, dataStack_aligned, cellIdMapStack, isContour, isImgBkg, isImgBkgSqrt, isWithSeg, isWithTrace)

%%%%   Inputs
% %   tracePxLstLst: a list of trace pixel lists
% %   (Each trace pixel list:  nT*1 cell array, cell i is the pixel list of the trace in frame i)


nTrace = length(tracePxLstLst);
[nR, nC, nT] = size(dataStack_aligned);


%%%%   Map stack of segmentation results
if isWithSeg
    segLbMapStack = false(nR, nC, nT);
    for iFM = 1:nT
        segLbMapStack(:,:,iFM) = imdilate(cellIdMapStack(:,:,iFM),strel('disk',1)) - cellIdMapStack(:,:,iFM) > 0;
    end
    % figure; imshow(segLbMapStack(:,:,116));
end


%%%  Map of regions in traces
lbMapStack = zeros(nR, nC, nT);
if isContour
    tempMap = false(nR, nC);
    for iTrc = 1:nTrace
        % iTrc = 142;
        disp(['Contour map: Trace No. ', num2str(iTrc)]);
        TOIs = find(~cellfun(@isempty, tracePxLstLst{iTrc}));
        for iiT = 1:length(TOIs)
            t = TOIs(iiT);
            tempMap(:) = false;
            tempMap(tracePxLstLst{iTrc}{t}) = true;
            tempMap = imdilate(tempMap, strel('disk',1)) & ~imerode(tempMap, strel('disk',1));
            % tempMap = tempMap & ~imerode(tempMap, strel('disk',1));
            [rR,cC] = find(tempMap);
            lbMapStack(sub2ind([nR,nC,nT], rR, cC, t*ones(length(rR),1))) = iTrc;
        end
    end
else
    for iTrc = 1:nTrace
        % iTrc = 142;
        disp(['ROI map: Trace No. ', num2str(iTrc)]);
        TOIs = find(~cellfun(@isempty, tracePxLstLst{iTrc}));
        for iiT = 1:length(TOIs)
            t = TOIs(iiT);
            [rR,cC] = ind2sub([nR,nC], tracePxLstLst{iTrc}{t});
            lbMapStack(sub2ind([nR,nC,nT], rR, cC, t*ones(size(rR)))) = iTrc;
        end
    end
end
% figure; imshow(lbMapStack(:,:,1));
% implay(lbMapStack);


%%%  Footprint map of detection center in traces
if isWithTrace
    for iTrc = 1:nTrace
        % iTrc = 142;
        disp(['Center footprint map: Trace No. ', num2str(iTrc)]);
        TOIs = find(~cellfun(@isempty, tracePxLstLst{iTrc}));
        
        dtcts = zeros(length(TOIs),3);
        for iiT = 1:length(TOIs)
            t = TOIs(iiT);
            [rVec,cVec] = ind2sub([nR,nC], tracePxLstLst{iTrc}{t});
            dtcts(iiT,:) = [round(mean([rVec,cVec],1)), t];
        end
        
        subVec = [];
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
end
% implay(lbMapStack);


%%%%   Image Stack of Result Presentation
showImgStack = zeros(nR, nC, 3, nT);
% cmap = jet(nTrace);
% [~,temp] = sort(rand(size(cmap,1),1));
% cmap = cmap(temp,:);
cmap = rand(nTrace,3)*0.85 + 0.15;

cmapSeg = [1,1,1] * 0.8;

for iFM = 1:nT
    lbMap = lbMapStack(:,:,iFM);
    %     lbMap = imdilate(lbMap, [1,1,1; 1,1,1; 1,1,1]);
    if isImgBkg
        % temp = sqrt(dataStack_aligned(:,:,iFM)) * 0.8;
        if isImgBkgSqrt
            temp = sqrt(dataStack_aligned(:,:,iFM));
        else
            temp = min(dataStack_aligned(:,:,iFM)*1.3, 1);
        end
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
% figure; imshow(showImgStack(:,:,:,24));





