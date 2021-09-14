function [numCellResVec, dtctIncldLst] = tempConsistSort(cellLst_cur, nullMtScPhat)



nnT = length(cellLst_cur);

%%  Prepare: single-frame segmentation map  &  index LUTs

nDtct = sum(cellfun(@(x) x.NumObjects, cellLst_cur));

%%% The overall cell representer indices
cellIdxsLst = cellfun(@(x) (1:x.NumObjects)', cellLst_cur, 'UniformOutput', false);
for t = 2:nnT
    cellIdxsLst{t}(:,1) = cellIdxsLst{t}(:,1) + cellIdxsLst{t-1}(end,1);
end

% cellSubIdxsLUT = cell(nnT,1);
% for t = 1:nnT
%     cellSubIdxsLUT{t} = [t + zeros(cellLst{t}.NumObjects,1), (1:cellLst{t}.NumObjects)'];
% end
% cellSubIdxsLUT = cell2mat(cellSubIdxsLUT);




%%  Infer the inclusion relationships among detections

dtctIncldLst = cell(nDtct,1);

for t0 = 1:nnT
    % t0 = 2;
    % disp(['Inferring the inclusion relationship among detections... t = ', num2str(t0)]);
        
    for iSubId0 = 1:cellLst_cur{t0}.NumObjects
        % iSubId0 = 10;
        iDtct0 = cellIdxsLst{t0}(iSubId0);
        dtct0.vxLst = cellLst_cur{t0}.VoxelIdxList{iSubId0};
        dtct0.mtSc = cellLst_cur{t0}.dtctMtScVec(iSubId0, :);
        dtctIncldLst{iDtct0}.lastFm = [];
        dtctIncldLst{iDtct0}.nextFm = [];
        
        %%%  Last frame
        if t0 > 1
            [childDtctSet, ~] = getChildSet(dtct0, cellLst_cur{t0-1}, 0, nullMtScPhat);
            dtctIncldLst{iDtct0}.lastFm = cellIdxsLst{t0-1}(childDtctSet);
        end
        %%%  Next frame
        if t0 < nnT
            [childDtctSet, ~] = getChildSet(dtct0, cellLst_cur{t0+1}, 1, nullMtScPhat);
            dtctIncldLst{iDtct0}.nextFm = cellIdxsLst{t0+1}(childDtctSet);
        end
    end
end




%%   Construct the "inclusion" network of detections
edgeVec = zeros(nDtct*10, 2);
countEdge = 0;
for iDtct0 = 1:nDtct
    % iDtct0 = 6000;
    if ~isempty(dtctIncldLst{iDtct0}.lastFm)
        temp = dtctIncldLst{iDtct0}.lastFm;
        edgeVec(countEdge+1 : countEdge+length(temp), :) = [iDtct0+zeros(length(temp),1), temp];
        countEdge = countEdge + length(temp);
    end
    if ~isempty(dtctIncldLst{iDtct0}.nextFm)
        temp = dtctIncldLst{iDtct0}.nextFm;
        edgeVec(countEdge+1 : countEdge+length(temp), :) = [iDtct0+zeros(length(temp),1), temp];
        countEdge = countEdge + length(temp);
    end
end
edgeVec(countEdge+1 : end, :) = [];


%%%  Directed graph of "inclusion" relationship among detections
GdtctIncld = digraph(edgeVec(:,1), edgeVec(:,2), [], nDtct);
% isdag(GdtctIncld)

% figure; plot(GdtctIncld, 'NodeLabel', (1:nDtct)',...
%     'layout', 'layered', 'AssignLayers','alap',...
%     'NodeColor', [0.7,0.7,0.7], 'MarkerSize', 10, 'EdgeColor', 'k', 'ArrowSize', 8);



%%  Objective function score

%%%  Method 2 of constructing Amat: with m2m relationships
% Amat = zeros(nDtct*3, nDtct);
% count = 0;
% for t = 1:nnT-1
%     tempRltdDtctIds = cell2mat(cellIdxsLst(t:t+1));
%     fmGraph = subgraph(GdtctIncld,tempRltdDtctIds);
%     % figure; plot(fmGraph, 'NodeLabel', tempRltdDtctIds,...
%     %     'layout', 'layered', 'AssignLayers','alap',...
%     %     'NodeColor', [0.7,0.7,0.7], 'MarkerSize', 10, 'EdgeColor', 'k', 'ArrowSize', 8);
%     
%     % % % Find the connected components just between two frames
%     cnnIdVec = conncomp(fmGraph, 'Type', 'weak')';
%     [binCnt, binID] = hist(cnnIdVec, unique(cnnIdVec));
%     binID(binCnt <= 1) = [];
%     
%     for iii = 1:length(binID)
%         iCnn = binID(iii);
%         candDtcts = tempRltdDtctIds(cnnIdVec == iCnn);
%         candDtctTs = cellSubIdxsLUT(candDtcts, 1);
%         Amat(count+1, candDtcts(candDtctTs==t)) = 1;
%         Amat(count+1, candDtcts(candDtctTs==t+1)) = -1;
%         count = count + 1;
%     end
% end
% Amat(count+1:end, :) = [];
% 
% 
% Nc = size(Amat,1);   % % Number of constraints




%%  Solve by MCC:  Construct the MCC flow network

nnNode = nDtct * 4 + 2;   % %  4 nodes for each detection: [pre, unit, extra, post]
source = nnNode - 1;
sink = nnNode;

infNum = nDtct * 100;

alphaFP = 3;
alphaFC = 1;
beta = 4;
% FPpnlty = 5;

strtIdx = ((1:nDtct)-1)'*4;

%%%%  Arcs: pre --> unit  (format: [<tail> <head> <capacity l.b.> <capacity u.b> <cost>])
arcSet1 = [strtIdx + 1, strtIdx + 2, zeros(nDtct,1), ones(nDtct,1), zeros(nDtct,1)];

%%%%  Arcs: unit --> post
arcSet2 = [strtIdx + 2, strtIdx + 4, zeros(nDtct,1), ones(nDtct,1), zeros(nDtct,1) - alphaFP];

%%%%  Arcs: pre --> extra
arcSet3 = [strtIdx + 1, strtIdx + 3, zeros(nDtct,1), zeros(nDtct,1) + infNum, zeros(nDtct,1)];

%%%%  Arcs: extra --> post
arcSet4 = [strtIdx + 3, strtIdx + 4, zeros(nDtct,1), zeros(nDtct,1) + infNum, zeros(nDtct,1) + alphaFC];

%%%%  Arcs: source --> pre (starting detections have special costs)
arcSet5 = [zeros(nDtct,1) + source, strtIdx + 1, zeros(nDtct,1), zeros(nDtct,1) + infNum,...
    zeros(nDtct,1) + beta];
% arcSet5(1:cellLst{1}.NumObjects, 5) = 0;
for iDtct = 1 : nDtct
    nbDtctIds = unique([predecessors(GdtctIncld, iDtct); successors(GdtctIncld, iDtct)]);
    if isempty(nbDtctIds) || (~isempty(nbDtctIds) && all(nbDtctIds > iDtct))
        arcSet5(iDtct, 5) = 0;
    end
end

%%%%  Arcs: post --> sink (ending detections have special costs)
arcSet6 = [strtIdx + 4, zeros(nDtct,1) + sink, zeros(nDtct,1), zeros(nDtct,1) + infNum,...
    zeros(nDtct,1) + beta];
% arcSet6((nDtct-cellLst{nnT}.NumObjects+1):nDtct, 5) = 0;
for iDtct = 1 : nDtct
    nbDtctIds = unique([predecessors(GdtctIncld, iDtct); successors(GdtctIncld, iDtct)]);
    if isempty(nbDtctIds) || (~isempty(nbDtctIds) && all(nbDtctIds < iDtct))
        arcSet6(iDtct, 5) = 0;
    end
end

%%%%  Arcs: transitions between frames
temp = GdtctIncld.Edges.EndNodes;
temp(temp(:,1)>temp(:,2), :) = [];
nTrans = size(temp,1);
arcSet7 = [strtIdx(temp) + ones(nTrans,1)*[4,1], zeros(nTrans,1), zeros(nTrans,1) + infNum,...
    zeros(nTrans,1)];

%%%%  Arc: sink --> source
% arcSet8 = [sink, source, 0, infNum, 0];

%%%%  The combined arc list
arcVec = [arcSet1; arcSet2; arcSet3; arcSet4; arcSet5; arcSet6; arcSet7];
% arcVec = [arcSet1; arcSet2; arcSet3; arcSet4; arcSet5; arcSet6; arcSet7; arcSet8];
nnArc = size(arcVec,1);


% Gmcf = digraph(arcVec(:,1), arcVec(:,2), arcVec(:,5), nnNode);
% isdag(GdtctIncld)
% figure; plot(Gmcf, 'NodeLabel', (1:nnNode)','layout', 'layered', 'AssignLayers','alap',...
%     'Sources', 1:dtctLst{1}.NumObjects, 'Sinks', ((nDtct-dtctLst{nnT}.NumObjects+1):nDtct)*4,...
%     'NodeColor', [0.7,0.7,0.7], 'MarkerSize', 10, 'EdgeColor', 'k', 'ArrowSize', 8);




%%  Solve by MCC:  by searching for best demand/supply

upbndSpplVal = nDtct;

tic

[minCC, optSpplVal, resFlowVec, spplValnCostVec] = myMCCbyCMF(arcVec, nnNode, nnArc, source, sink, upbndSpplVal);

toc

disp(['minCC:  ', num2str(minCC)]);
disp(['optSpplVal:  ', num2str(optSpplVal)]);
% disp('spplValnCostVec: ');
% disp(spplValnCostVec);



%%  Interprete the results into conclusions
dtctFcsdFlowMat = int64(zeros(nDtct,6));
for iArc = 1:size(resFlowVec,1)
    dtctIdEnds = int64(ceil(double(resFlowVec(iArc,1:2))/4));
    if dtctIdEnds(1) == dtctIdEnds(2)
        nodeSubEnds = resFlowVec(iArc,1:2) - (dtctIdEnds-1)*4;
        if(nodeSubEnds(1) == 1 && nodeSubEnds(2) == 2)
            dtctFcsdFlowMat(dtctIdEnds(1),2) = dtctFcsdFlowMat(dtctIdEnds(1),2) + resFlowVec(iArc,3);
        end
        if(nodeSubEnds(1) == 2 && nodeSubEnds(2) == 4)
            dtctFcsdFlowMat(dtctIdEnds(1),3) = dtctFcsdFlowMat(dtctIdEnds(1),3) + resFlowVec(iArc,3);
        end
        if(nodeSubEnds(1) == 1 && nodeSubEnds(2) == 3)
            dtctFcsdFlowMat(dtctIdEnds(1),4) = dtctFcsdFlowMat(dtctIdEnds(1),4) + resFlowVec(iArc,3);
        end
        if(nodeSubEnds(1) == 3 && nodeSubEnds(2) == 4)
            dtctFcsdFlowMat(dtctIdEnds(1),5) = dtctFcsdFlowMat(dtctIdEnds(1),5) + resFlowVec(iArc,3);
        end
    elseif resFlowVec(iArc,1) == source
        dtctFcsdFlowMat(dtctIdEnds(2),1) = dtctFcsdFlowMat(dtctIdEnds(2),1) + resFlowVec(iArc,3);
    elseif resFlowVec(iArc,2) == sink
        dtctFcsdFlowMat(dtctIdEnds(1),6) = dtctFcsdFlowMat(dtctIdEnds(1),6) + resFlowVec(iArc,3);
    end
end


numCellResVec = double(dtctFcsdFlowMat(:,2) + dtctFcsdFlowMat(:,4));


% disp(['Number of false positives: ',num2str(nnz(numCellResVec == 0)), ' out of ', num2str(nDtct)]);

% eqViolationTerm = sum(dtctFcsdFlowMat(cellLst{1}.NumObjects+1:end, 1)) + ...
%     sum(dtctFcsdFlowMat(1:(nDtct-cellLst{nnT}.NumObjects), 6));
% 
% oriSegErrorTerm =  sum(abs(dtctFcsdFlowMat(:, 2)+dtctFcsdFlowMat(:, 4)-1));



% %%%%  Show result in graph
% temp = GdtctIncld.Edges.EndNodes;
% temp(temp(:,1) > temp(:,2), :) = temp(temp(:,1) > temp(:,2), end:-1:1);
% temp = unique(temp, 'rows');
% 
% GdtctRes = digraph(temp(:,1), temp(:,2), [], nDtct);
% 
% tempNodeLbels = numCellResVec;
% % temp = num2str([(1:nDtct)', numCellResVec]);
% % tempNodeLbels = mat2cell(temp, ones(size(temp,1),1), size(temp,2));
% figure; plot(GdtctRes, 'NodeLabel', tempNodeLbels,'layout', 'layered', 'AssignLayers','alap',...
%     'Sources', 1:cellLst{1}.NumObjects, 'Sinks', (nDtct-cellLst{nnT}.NumObjects+1):nDtct,...
%     'NodeColor', [0.7,0.7,0.7], 'MarkerSize', 10, 'EdgeColor', 'k', 'ArrowSize', 8);





%%  Visualization

% dtctMtScAllVec = cell2mat(cellfun(@(x) x.dtctMtScVec, cellLst, 'UniformOutput', false));
% 
% %%%%%   Rapid motion related 
% dtctPvalAllVec = gamcdf(dtctMtScAllVec, nullMtScPhat(1), nullMtScPhat(2), 'upper');
% % dtctPvalVec = normcdf(dtctMtScVec, mu, sigma);
% 
% pThr = 0.005;
% rpdMtRltdDtcts = dtctPvalAllVec < pThr;
% % disp(nnz(rpdMtRltdDtcts));
% % disp(nnz(sum(rpdMtRltdDtcts,2)>0));
% % scThr = gaminv(1-pThr, nullMtScPhat(1), nullMtScPhat(2));
% % figure; histogram(dtctPvalVec(:), 'BinWidth',0.01,'Normalization','probability'); xlim([0,1]);
% 
% % dtctOIs = find(rpdMtRltdDtcts(:,1) | rpdMtRltdDtcts(:,2));
% dtctOIs = find(rpdMtRltdDtcts(:,2));
% % dtctOIs = find(sum(abs(rpdMtRltdDtctPostP1-rpdMtRltdDtctPostP1_old)>0.01, 2)>0);
% 
% 
% %%%%%   Isolated detections in the correspondence network
% dtctOIs = find(indegree(GdtctIncld) + outdegree(GdtctIncld) == 0);
%     
% 
% %%%%%   False positives'
% dtctOIs = find(numCellResVec == 0);
% dtctOIs = find(numCellResVec > 1);
% 
% 
% %%%%%   Given detection(s)
% dtctOIs = find(numCellResVec ~= 1 & indegree(GdtctIncld) + outdegree(GdtctIncld) == 0);
% 
% 
% %%%%%  Show the significant detections
% traceVxLstLst = cell(length(dtctOIs), 1);
% for ii = 1:length(dtctOIs)
%     traceVxLstLst{ii} = cell(nT,1);
%     iDtct = dtctOIs(ii);
%     tt = cellSubIdxsLUT(iDtct,1);
%     iSub = cellSubIdxsLUT(iDtct,2);
%     traceVxLstLst{ii}{tt} = cellLst{tt}.VoxelIdxList{iSub};
% end
% 
% [nnY, nnX, nnZ, nnT] = size(imgStack);
% 
% imgSzNrg.size = [nnY,nnX,nnZ,nnT];
% imgSzNrg.range = [1,nnY; 1,nnX; 1,nnZ; 1,nnT];
% 
% % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, [], false, []);
% % % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, (cellIdPatch>0)*1, false, []);
% [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, false, []);
% % [printImgStack,printImgYZStack,~] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCpatchStack, true, cmapLb);
% implay(printImgStack);
% implay(printImgYZStack);
% 
% 
% % tt = 75;
% % figure; imshow(cellIdPatchStack(:,:,(22-1)*3+1,tt));
% % iSub = 57;
% % iDtct = cellIdxsLst{tt}(iSub)
% % iTrc = find(cellfun(@(x) any(x.dtctSubs(:,1)==tt & x.dtctSubs(:,2)==iSub), traceLst))











