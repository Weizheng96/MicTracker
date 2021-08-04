function cmpltOpcTrcLst = linkByIntegrtedMCF(cellLst_cur, imgStack, nullMtScPhat, ...
    altMtScPhat, apprPara, posParaStatic)


[nnY, nnX, nnZ, nnT] = size(imgStack);


%%  Prepare: current segmentation map  &  index LUTs

% nDtct = sum(cellfun(@(x) x.NumObjects, cellLst_cur));

%%% The overall cell representer indices
cellIdxsLst = cellfun(@(x) (1:x.NumObjects)', cellLst_cur, 'UniformOutput', false);
for t = 2:nnT
    cellIdxsLst{t}(:,1) = cellIdxsLst{t}(:,1) + cellIdxsLst{t-1}(end,1);
end

cellSubIdxsLUT = cell(nnT,1);
for t = 1:nnT
    cellSubIdxsLUT{t} = [t + zeros(cellLst_cur{t}.NumObjects,1), (1:cellLst_cur{t}.NumObjects)'];
end
cellSubIdxsLUT = cell2mat(cellSubIdxsLUT);


%%  Parameter setting

%%%%  Appearance parameters
% apprPara.intDiff_sigma = apprPara.intDiff_sigma * 3;
% figure; histogram(randn(10000,1)*apprPara.intDiff_sigma*9+apprPara.intDiff_mu, 'BinWidth', 0.01, 'Normalization', 'pdf'); xlim([-1,1]);
% hold on; histogram(randn(10000,1)*apprPara.intDiff_sigma+apprPara.intDiff_mu, 'BinWidth', 0.01, 'Normalization', 'pdf'); xlim([-1,1]);


%%%%  Motion parameters: static
xStaticPara = posParaStatic;


%%%%  Position (motion) parameters: rapid moving
xRpMtPara = posParaStatic;
xRpMtPara.DsplcmtSigma = posParaStatic.DsplcmtSigma * 15 ^ 2;
temp = [randn(10000,1)*sqrt(xRpMtPara.DsplcmtSigma(1,1)) + xRpMtPara.DsplcmtMu(1),...
    randn(10000,1)*sqrt(xRpMtPara.DsplcmtSigma(2,2)) + xRpMtPara.DsplcmtMu(2),...
    randn(10000,1)*sqrt(xRpMtPara.DsplcmtSigma(3,3)) + xRpMtPara.DsplcmtMu(3)];
tempDist = sqrt(sum(temp.^2,2));
xRpMtPara.distPhat = gamfit(tempDist(tempDist>0));

% temp = gamrnd(xRpMtPara.distPhat(1),xRpMtPara.distPhat(2),10000,1);
% figure; histogram(temp(:,1),'BinWidth', 1, 'Normalization', 'pdf'); 
% temp = gamrnd(xStaticPara.distPhat(1),xStaticPara.distPhat(2),10000,1);
% hold on; histogram(temp(:,1),'BinWidth', 1, 'Normalization', 'pdf'); xlim([0,60]);


%%%%   Posterior probability of rapid moving
dtctMtScAllVec = cell2mat(cellfun(@(x) x.dtctMtScVec, cellLst_cur, 'UniformOutput', false));
% dtctMtScAllVec(isnan(dtctMtScAllVec)) = 0;
temp0 = gampdf(dtctMtScAllVec, nullMtScPhat(1), nullMtScPhat(2));
temp1 = gampdf(dtctMtScAllVec, altMtScPhat(1), altMtScPhat(2));
postPvec1 =  temp1 ./ (temp0 + temp1);
% postPvec1(isnan(postPvec1)) = 0;
% postPvec0 = 1 - postPvec1;
rpdMtRltdDtctPostP1 = postPvec1;
rpdMtRltdDtctPostP1(isnan(rpdMtRltdDtctPostP1)) = 0;
dtctPvalVec = gamcdf(dtctMtScAllVec, nullMtScPhat(1), nullMtScPhat(2), 'upper');
% figure; histogram(postPvec1(:),'BinWidth',0.01,'Normalization','probability'); xlim([0,1]);
% figure; histogram(postPvec1(~rpdMtRltdDtcts),'BinWidth',0.01,'Normalization','probability'); xlim([0,1]);
% hold on; histogram(postPvec1(rpdMtRltdDtcts),'BinWidth',0.01,'Normalization','probability'); xlim([0,1]);
% temp = 0:0.001:1;
% temp0 = gampdf(temp, nullMtScPhat(1), nullMtScPhat(2));
% temp1 = gampdf(temp, altMtScPhat(1), altMtScPhat(2));
% tempPostPvec1 =  temp1 ./ (temp0 + temp1);
% figure; plot(temp, tempPostPvec1,'LineWidth',2);





%%%%  Other auxiliary parameters
para.nnT = nnT;
para.nnY = nnY;
para.nnX = nnX;
para.nnZ = nnZ;

para.maxJump = 5;
para.FPdtctRate = 5.0000e-4;
% para.FPdtctRate = 0.0025;
% disp(-log(para.FPdtctRate/(1-para.FPdtctRate)));
% para.pLinkThr = 1.0000e-10;
% pCtrPos = mvnpdf(sqrt(diag(xRpMtPara.DsplcmtSigma)')*2, xRpMtPara.DsplcmtMu, xRpMtPara.DsplcmtSigma);
pCtrPos = gampdf(xRpMtPara.distPhat(1)*xRpMtPara.distPhat(2),xRpMtPara.distPhat(1),xRpMtPara.distPhat(2));
pIntst = normpdf(sqrt(apprPara.intDiff_sigma)*2, apprPara.intDiff_mu, apprPara.intDiff_sigma);
% pSize = normpdf(sqrt(aRpMtPara.sizeChange_sigma), aRpMtPara.sizeChange_mu, aRpMtPara.sizeChange_sigma);
pSize = normpdf(sqrt(apprPara.sizeChangeRt_sigma)*2 ,apprPara.sizeChangeRt_mu, apprPara.sizeChangeRt_sigma);
pOvlpRt = exppdf(1, xStaticPara.invOvrlpRtPhat);
pNbIncldRt = exppdf(1, xStaticPara.nbIncldRtPhat);
pT = exppdf(1+para.maxJump,(para.maxJump+1)/2);
para.pLinkThr = min(pIntst.*pSize.*pCtrPos.*pOvlpRt.*pNbIncldRt.*pT, para.FPdtctRate/(1-para.FPdtctRate));



%%  Integrated-motion-model MCF linking
%%%%  Construct the graph
[trackG, para] = buildTrackGraph_linkDetections_twoMtModel_3D(cellLst_cur,...
    rpdMtRltdDtctPostP1, apprPara, xStaticPara, xRpMtPara, para);


%%%%  Solve the graph
[~, trajectories] = sspTracker(trackG, para);


%%%%  Translate the result into traces
isIncldFPs = true;
traceDtctIdSet = traceTranslator_detection(trajectories, cellSubIdxsLUT, cellIdxsLst, isIncldFPs);


imgSzNrg.size = [nnY,nnX,nnZ,nnT];
imgSzNrg.range = [1,nnY; 1,nnX; 1,nnZ; 1,nnT];
traceLst = traceCreate_3D(traceDtctIdSet, cellLst_cur, imgSzNrg);

nTrace = length(traceLst);




% %%%%%%%%%%%%%%%%       Visualization of results      %%%%%%%%%%%%%%%%%%%
% %%%%  Show traceDtctIdSet (MCF result)
% traceVxLstLst = cell(length(traceDtctIdSet),1);
% trace0 = [];
% for iTrace = 1:length(traceDtctIdSet)
%     % iTrace = 203;
%     trace0.dtctSubs = traceDtctIdSet{iTrace};
%     trace0.traceVxLst = cell(nnT,1);
%     temp = mat2cell(trace0.dtctSubs, ones(size(trace0.dtctSubs,1),1), 2);
%     trace0.traceVxLst(trace0.dtctSubs(:,1)) = cellfun(@(x) cellLst_cur{x(1)}.VoxelIdxList{x(2)}, temp, 'UniformOutput', false);
%     traceVxLstLst{iTrace} = trace0.traceVxLst;
% end
% 
% 
% %%%%   Visualization: Show starting/ending detections in traces
% showTrLst = traceLst;
% traceVxLstLst = {cell(nnT,1); cell(nnT,1); cell(nnT,1)};
% for tt = 1:nnT
%     tempVec = cellfun(@(x) x.stT == tt && x.endT ~= tt, showTrLst);
%     traceVxLstLst{1}{tt} = cell2mat(cellfun(@(x) x.traceVxLst{tt}, showTrLst(tempVec), 'UniformOutput', false));
%     tempVec = cellfun(@(x) x.stT == tt && x.endT == tt, showTrLst);
%     traceVxLstLst{2}{tt} = cell2mat(cellfun(@(x) x.traceVxLst{tt}, showTrLst(tempVec), 'UniformOutput', false));
%     tempVec = cellfun(@(x) x.stT ~= tt && x.endT == tt, showTrLst);
%     traceVxLstLst{3}{tt} = cell2mat(cellfun(@(x) x.traceVxLst{tt}, showTrLst(tempVec), 'UniformOutput', false));
% end
% 
% 
% %%%%    Print
% % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, [], false, []);
% % % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, (cellIdPatchStack>0)*1, false, []);
% [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, false, []);
% % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, imgStack, false, []);
% % [printImgStack,printImgYZStack,~] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, true, cmapLb);
% implay(printImgStack);
% implay(printImgYZStack);
% % figure; imshow(printImgStack(:,:,:,25));
% 
% 
% % tt = 84;
% % figure; imshow(cellIdPatchStack(:,:,(35-1)*3+1,tt));
% % iSub = 2;
% % iDtct = cellIdxsLst{tt}(iSub)
% % iTrc = find(cellfun(@(x) any(x.dtctSubs(:,1)==tt & x.dtctSubs(:,2)==iSub), traceLst))
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%  Analyze the motion pattern of the composite traces, looking for the OPCs' basic traces
% opcWholeTrcLst
tic

trcMtnDistMat = nan(nTrace, nnT);
for iTrace = 1:nTrace
    trace0 = traceLst{iTrace};
    effctvTs = find(trace0.isDtctd);
    for jj = 1 : length(effctvTs)-1
        tt1 = effctvTs(jj);
        tt2 = effctvTs(jj + 1);
        trcMtnDistMat(iTrace, tt1:tt2-1) = ...
            sqrt(sum((trace0.cntrPtSub(tt2,:) - trace0.cntrPtSub(tt1,:)).^2)) / (tt2 - tt1);
    end
end
% figure; plot(trcMtnDistMat');

empNullVlctyVec = trcMtnDistMat(~isnan(trcMtnDistMat));
% figure; histogram(empNullVlctyVec, 100);



%%%%%   Fit probability distribution model to the empirical distribution
phat_vlcty = gamfit(empNullVlctyVec);
% temp = gamrnd(phat_vlcty(1),phat_vlcty(2),10000,1);
% figure; histogram(temp,100); xlim([0,60]);


%%%%%  Test: instantaneous speed
vlctyPvalMat = gamcdf(trcMtnDistMat, phat_vlcty(1), phat_vlcty(2), 'upper');

pThr = 0.001;
vlctySigLbMat = vlctyPvalMat < pThr;

% opcBasicTrcIds = find(sum(vlctySigLbMat, 2) > 0);
% opcBasicTrcLst = traceLst(opcBasicTrcIds);
% nOpcBasicTrc = length(opcBasicTrcLst);
% disp(nOpcBasicTrc);



%%%%%  Find OPCs as traces with appearance change
pThr = 0.001;

apprChgSigLbMat = false(nTrace, nnT);
for iTrace = 1:nTrace
    % iTrace = 389;
    trace0 = traceLst{iTrace};
    rltdDtctIdxs = zeros(size(trace0.dtctSubs,1),1);
    for iT = 1:size(trace0.dtctSubs,1)
        rltdDtctIdxs(iT) = cellIdxsLst{trace0.dtctSubs(iT,1)}(trace0.dtctSubs(iT,2));
    end
    rltdDtctPvalVec = dtctPvalVec(rltdDtctIdxs,:);
    apprChgSigLbMat(iTrace, trace0.dtctSubs(1:end-1,1)) = rltdDtctPvalVec(1:end-1,2) < pThr | rltdDtctPvalVec(2:end,1) < pThr;
end
% rpdMtRltdTrcLbs = sum(apprChgSigLbMat,2) > 0;

% opcBasicTrcIds = find(rpdMtRltdTrcLbs);
% opcBasicTrcIds = find(rpdMtRltdTrcLbs & sum(vlctySigLbMat, 2) > 0);
opcBasicTrcIds = find(sum(apprChgSigLbMat & vlctySigLbMat, 2) > 0);
opcBasicTrcLst = traceLst(opcBasicTrcIds);
nOpcBasicTrc = length(opcBasicTrcLst);
disp(nOpcBasicTrc);






% %%%%%%%%%%%%%%%%%%%%       Visualization       %%%%%%%%%%%%%%%%%%%%
% 
% %%%%  Show motion-broken static traces with the significant detections
% traceVxLstLst = cellfun(@(x) x.traceVxLst, opcBasicTrcLst, 'UniformOutput', false);
% 
% % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, [], false, []);
% % % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, (cellIdPatchStack>0)*1, false, []);
% [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, false, []);
% % [printImgStack,printImgYZStack,~] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, true, cmapLb);
% implay(printImgStack);
% implay(printImgYZStack);
% 
% 
% % tt = 53;
% % figure; imshow(cellIdPatchStack(:,:,(18-1)*3+1,tt));
% % iSub = 92;
% % iDtct = cellIdxsLst{tt}(iSub)
% % iTrc = find(cellfun(@(x) any(x.dtctSubs(:,1)==tt & x.dtctSubs(:,2)==iSub), traceLst))
% % iOpcBasicTrc = find(opcBasicTrcIds == iTrc)
% % find(apprChgSigLbMat(iTrc,:))
% % find(vlctySigLbMat(iTrc,:))
% % 
% % % speedVec = trcMtnDistMat(iTrc,:);
% % % disp(speedVec)
% % trace0 = traceLst{iTrc};
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%   For incomplete opc trace candidates, complete them or elongate them with related traces

%%%%%   Find all traces related with OPC basic traces
intpltdTrcLst = cell(nTrace, 1);
rltdTrcFoundLbs = false(nTrace, 1);
rltdTrcSetLst = cell(nTrace, 1);
% img1 = zeros(nnY,nnX,nnZ);
for iOpcBasicTrc = 1:nOpcBasicTrc
    % iOpcBasicTrc = 15;
    iTrc0 = opcBasicTrcIds(iOpcBasicTrc);
    if rltdTrcFoundLbs(iTrc0)
        continue;
    end
    candTrcIds = [iTrc0];
    grpIsDtctdVec = false(nnT,1);

    while ~isempty(candTrcIds)
        disp(['Current candTrcIds: ', num2str(candTrcIds')]);
        candTrcIds_new = [];
        for ii = 1:length(candTrcIds)
            iTrc = candTrcIds(ii);
            trace0 = traceLst{iTrc};
            if any(grpIsDtctdVec & trace0.isDtctd)
                continue;
            end
            grpIsDtctdVec = grpIsDtctdVec | trace0.isDtctd;
            
            %%%  Get the interpolated trace
            if ~isempty(intpltdTrcLst{iTrc})
                intpltdTrace0 = intpltdTrcLst{iTrc};
            else
                intpltdTrace0 = traceInterpolate_3D_trcVersion(trace0, imgSzNrg);
                disp(['Interpolating trace: iTrc = ', num2str(iTrc)]);
                intpltdTrcLst{iTrc} = intpltdTrace0;
            end
            
            %%%  Get the related traces
            tempVec = ~intpltdTrace0.isDtctd;
            tempVec([1:intpltdTrace0.stT, intpltdTrace0.endT+1:nnT]) = false;
            incmpltTs = find(tempVec);
            rltdTrcIdVec = nan(length(incmpltTs), 1);
            for iT = 1:length(incmpltTs)
                tt = incmpltTs(iT);
                vxLst = intpltdTrace0.traceVxLst{tt};
                candCounts = cellfun(@(x) length(intersect(vxLst, x)), cellLst_cur{tt}.VoxelIdxList);
                candIds = find(candCounts > 0);
                if isempty(candIds)
                    continue;
                elseif length(candIds)==1
                    iSub = candIds;
                else
                    counts = candCounts(candIds);
                    [~,jj] = max(counts);
                    iSub = candIds(jj);
                    vxLst2 = cellLst_cur{tt}.VoxelIdxList{iSub};
                    if length(intersect(vxLst,vxLst2))/length(union(vxLst,vxLst2)) < 0.5
                        continue;
                    end
                end
                jTrc = find(cellfun(@(x) any(x.dtctSubs(:,1)==tt & x.dtctSubs(:,2)==iSub), traceLst));
                rltdTrcIdVec(iT) = jTrc;
            end
            rltdTrcSetLst{iTrc}.incmpltTs = incmpltTs;
            rltdTrcSetLst{iTrc}.rltdTrcIdVec = rltdTrcIdVec;
            rltdTrcFoundLbs(iTrc) = true;
            candTrcIds_new = [candTrcIds_new; unique(rltdTrcIdVec(~isnan(rltdTrcIdVec)))];
        end
        candTrcIds = unique(candTrcIds_new);
        candTrcIds(rltdTrcFoundLbs(candTrcIds)) = [];
    end
end



%%%%%   Find sets of mutually related traces
trcOIids = find(rltdTrcFoundLbs);
trcAssctMatEdges = zeros(length(trcOIids)^2, 2);
countEdge = 0;
for ii = 1:length(trcOIids)
    % iTrc = 15;
    iTrc = trcOIids(ii);
    rltdTrcIdVec = rltdTrcSetLst{iTrc}.rltdTrcIdVec;
    candTrcIds = unique(rltdTrcIdVec(~isnan(rltdTrcIdVec)));
    trcAssctMatEdges(countEdge+1:countEdge+length(candTrcIds), :) = [iTrc+zeros(length(candTrcIds),1), candTrcIds];
    countEdge = countEdge + length(candTrcIds);
    trcAssctMatEdges(countEdge+1:countEdge+length(candTrcIds), :) = [candTrcIds, iTrc+zeros(length(candTrcIds),1)];
    countEdge = countEdge + length(candTrcIds);
end
trcAssctMatEdges(countEdge+1:end, :) = [];
trcAssctMatEdges = unique(trcAssctMatEdges, 'rows');
trcAssctMat = sparse(trcAssctMatEdges(:,1), trcAssctMatEdges(:,2), ones(size(trcAssctMatEdges,1),1), nTrace, nTrace);
GopcBasicRltdTrc = graph(trcAssctMat);
opcBasicTrcGrpIdVec = conncomp(GopcBasicRltdTrc);
opcBasicTrcGrpIdVec(~rltdTrcFoundLbs) = 0;
[grpFreq, grpIds] = hist(opcBasicTrcGrpIdVec, unique(opcBasicTrcGrpIdVec));
grpFreq(grpIds == 0) = [];
grpIds(grpIds == 0) = [];




% %%%%%%%%%%%%%%%%%%%%       Visualization       %%%%%%%%%%%%%%%%%%%%
% 
% %%%%  Show motion-broken static traces with the significant detections
% iGrp = 2;
% traceVxLstLst = cellfun(@(x) x.traceVxLst, traceLst(opcBasicTrcGrpIdVec==grpIds(iGrp)), 'UniformOutput', false);
% 
% % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, [], false, []);
% % % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, (cellIdPatchStack>0)*1, false, []);
% [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, false, []);
% % [printImgStack,printImgYZStack,~] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, true, cmapLb);
% implay(printImgStack);
% % implay(printImgYZStack);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%   Merge members in the same group
newOpcTrcLst = cell(length(grpIds),1);
newIntpltdOpcTrcLst = cell(length(grpIds),1);
for iiGrpOI = 1:length(grpIds)
    iGrp = grpIds(iiGrpOI);
    candTrcIds = find(opcBasicTrcGrpIdVec == iGrp);
    if length(candTrcIds) == 1
        newOpcTrcLst{iiGrpOI} = traceLst{candTrcIds};
        newIntpltdOpcTrcLst{iiGrpOI} = intpltdTrcLst{candTrcIds};
    else
        candTrcs = traceLst(candTrcIds);
        trace0 = candTrcs{1};
        trace0.stT = min(cellfun(@(x) x.stT, candTrcs));
        trace0.endT = max(cellfun(@(x) x.endT, candTrcs));
        trace0.dtctSubs = nan(1,2);
        for tt = trace0.stT:trace0.endT
            tempVec = cellfun(@(x) x.isDtctd(tt), candTrcs);
            trace0.isDtctd(tt) = any(tempVec);
            if ~trace0.isDtctd(tt)
                continue;
            end
            tempVec2 = cellfun(@(x) length(x.traceVxLst{tt}), candTrcs);
            tempVec2(~tempVec) = 0;
            [~,iC] = max(tempVec2);
            trace0.traceVxLst{tt} = candTrcs{iC}.traceVxLst{tt};
            trace0.posBoundsPerT(tt,:) = candTrcs{iC}.posBoundsPerT(tt,:);
            trace0.cntrPtSub(tt,:) = candTrcs{iC}.cntrPtSub(tt,:);
            trace0.dtctSubs(end+1,:) = candTrcs{iC}.dtctSubs(candTrcs{iC}.dtctSubs(:,1)==tt,:);
        end
        trace0.dtctSubs(1,:) = [];
        temp = [min(trace0.posBoundsPerT); max(trace0.posBoundsPerT)];
        trace0.posBounds = temp(logical([1,0,1,0,1,0; 0,1,0,1,0,1]))';
        
        newOpcTrcLst{iiGrpOI} = trace0;
        newIntpltdOpcTrcLst{iiGrpOI} = traceInterpolate_3D_trcVersion(trace0, imgSzNrg);
        disp(['Interpolating new trace for group: iiGrpOI = ', num2str(iiGrpOI)]);
    end    
end



% %%%%%%%%%%%%%%%%%%%%       Visualization       %%%%%%%%%%%%%%%%%%%%
% 
% %%%%  Show motion-broken static traces with the significant detections
% traceVxLstLst = cellfun(@(x) x.traceVxLst,newIntpltdOpcTrcLst, 'UniformOutput', false);
% 
% % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, [], false, []);
% % % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, (cellIdPatchStack>0)*1, false, []);
% [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, false, []);
% % [printImgStack,printImgYZStack,~] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, true, cmapLb);
% implay(printImgStack);
% implay(printImgYZStack);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%  Further exclude confounders using long-term (whole-trace) features

%%%%%%   Using largest displacement in any time window

%%%%  Get the null distribution
nullMaxDistVec = nan(nTrace, 1);
for iTrc = 1:nTrace
    trace0 = traceLst{iTrc};
    if trace0.endT == trace0.stT
        continue;
    end
    tempMaxDistVec = zeros(trace0.endT - trace0.stT, 1);
    for tDiff = 1:(trace0.endT - trace0.stT)
        tt2Vec = trace0.stT+tDiff : trace0.endT;
        tt1Vec = trace0.stT : trace0.endT-tDiff;
        tempDist = sqrt(sum((trace0.cntrPtSub(tt2Vec,:) - trace0.cntrPtSub(tt1Vec,:)).^2, 2));
        tempMaxDistVec(tDiff) = max(tempDist);
    end
    nullMaxDistVec(iTrc) = max(tempMaxDistVec);
end
% figure; histogram(nullMaxDistVec, 50);

maxDistPhat = gamfit(nullMaxDistVec(~isnan(nullMaxDistVec)));
% figure; histogram(gamrnd(maxDistPhat(1), maxDistPhat(2), 10000, 1), 100); xlim([0,50]);



%%%%  Test on the OPC trace candidates
nOPC = length(newIntpltdOpcTrcLst);
maxDistOpcVec = nan(nOPC, 1);
for iOpcTrc = 1:nOPC
    intpltdTrace0 = newIntpltdOpcTrcLst{iOpcTrc};
    if intpltdTrace0.endT == intpltdTrace0.stT
        continue;
    end
    tempMaxDistVec = zeros(intpltdTrace0.endT - intpltdTrace0.stT, 1);
    for tDiff = 1:(intpltdTrace0.endT - intpltdTrace0.stT)
        tt2Vec = intpltdTrace0.stT+tDiff : intpltdTrace0.endT;
        tt1Vec = intpltdTrace0.stT : intpltdTrace0.endT-tDiff;
        tempDist = sqrt(sum((intpltdTrace0.cntrPtSub(tt2Vec,:) - intpltdTrace0.cntrPtSub(tt1Vec,:)).^2, 2));
        tempMaxDistVec(tDiff) = max(tempDist);
    end
    maxDistOpcVec(iOpcTrc) = max(tempMaxDistVec);
end

maxDistPvalVec = gamcdf(maxDistOpcVec, maxDistPhat(1), maxDistPhat(2), 'upper');

pThr = 0.01;
maxDistSigLbVec = maxDistPvalVec < pThr;


cmpltOpcTrcLst = newIntpltdOpcTrcLst(maxDistSigLbVec);

nOPC = length(cmpltOpcTrcLst);
disp(['Final number of OPC traces: ', num2str(nOPC)]);


% %%%%%%%%%%%%%%%%%%%%       Visualization       %%%%%%%%%%%%%%%%%%%%
% 
% %%%%  Show motion-broken static traces with the significant detections
% traceVxLstLst = cellfun(@(x) x.traceVxLst, cmpltOpcTrcLst, 'UniformOutput', false);
% 
% % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, [], false, []);
% % % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, (cellIdPatchStack*1, false, []);
% % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, false, []);
% [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, imgStack, false, []);
% % [printImgStack,printImgYZStack,~] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, true, cmapLb);
% implay(printImgStack);
% implay(printImgYZStack);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 





