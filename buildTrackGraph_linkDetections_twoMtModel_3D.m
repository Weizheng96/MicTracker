function [trackG, para] = buildTrackGraph_linkDetections_twoMtModel_3D(cellLst,...
    rpdMtRltdDtctPostP1, aPara, xStaticPara, xRpMtPara, para)
% %  Input:
% %     cellLst  -   nFM*1 list of detection results
% %     rpdMtRltdDtctPostP1   -   detections' probability of "rapid motion related", nDtct*2 matrix
% %     aPara   - perameters w.r.t the appearance model
% %     xStaticPara   - perameters w.r.t the static-status motion model
% %     xRpMtPara   - perameters w.r.t the rapid-status motion model
% %     para	-   (structure: priors & hyper parameters for tracking)
% %  Output:
% %     trackG      -   digraph of tracking (min-cost flow model)

nnT = para.nnT;
% nnX = para.nnX;
% nnY = para.nnY;
% nnZ = para.nnZ;


%%%  The overall cell representer indices
cellIdxsLst = cellfun(@(x) (1:x.NumObjects)', cellLst, 'UniformOutput', false);
for t = 2:nnT
    cellIdxsLst{t}(:,1) = cellIdxsLst{t}(:,1) + cellIdxsLst{t-1}(end,1);
end
para.nCellRep = cellIdxsLst{nnT}(end);


%%%  Get the dummy scores
% pIntstRef = normpdf(aPara.intDiff_mu, aPara.intDiff_mu, aPara.intDiff_sigma);
% pSizeChangeRtRef = normpdf(aPara.sizeChangeRt_mu, aPara.sizeChangeRt_mu, aPara.sizeChangeRt_sigma);

% pDsplcmtRefRpMt = mvnpdf(xRpMtPara.DsplcmtMu, xRpMtPara.DsplcmtMu, xRpMtPara.DsplcmtSigma);
% pDsplcmtRefStatic = mvnpdf(xStaticPara.DsplcmtMu, xStaticPara.DsplcmtMu, xStaticPara.Sigma);
% pInvOvrlpRtRefRpMt = exppdf(xRpMtPara.invOvrlpRtPhat, xRpMtPara.invOvrlpRtPhat);
% pInvOvrlpRtRefStatic = exppdf(xStaticPara.invOvrlpRtPhat, xStaticPara.invOvrlpRtPhat);
pInvOvrlpRtRefStatic = 1;
% pNbIncldRtRefRpMt = exppdf(xRpMtPara.nbIncldRtPhat, xRpMtPara.nbIncldRtPhat);
% pNbIncldRtRefStatic = exppdf(xStaticPara.nbIncldRtPhat, xStaticPara.nbIncldRtPhat);
% pNbIncldRtRefStatic = exppdf(0, xStaticPara.nbIncldRtPhat);
pNbIncldRtRefStatic = 1;

pInvOvrlpRtStaticBase = 1 / expcdf(1, xStaticPara.invOvrlpRtPhat);
pNbIncldRtStaticBase = 1 / expcdf(1, xStaticPara.nbIncldRtPhat);



nNode = para.nCellRep * 2 + 2;  % % i-th cell representor has (2*i-1)-th & (2*i)-th nodes
source = nNode - 1;
sink = nNode;
para.source = source;
para.sink = sink;

arcTails = zeros(para.nCellRep*3, 1);
arcHeads = zeros(para.nCellRep*3, 1);
arcCosts = zeros(para.nCellRep*3, 1);


%%%  P(yi)
c_yi = log(para.FPdtctRate/(1-para.FPdtctRate));
arcTails(1 : para.nCellRep) = (1:2:nNode-2)';
arcHeads(1 : para.nCellRep) = (2:2:nNode-2)';
arcCosts(1 : para.nCellRep) = c_yi;


%%%  P_enter, P_exit
c_enter = - c_yi / 2 - 10^(-10);
c_exit = c_enter;

arcTails(para.nCellRep+1 : 2*para.nCellRep) = source;
arcHeads(para.nCellRep+1 : 2*para.nCellRep) = (1:2:nNode-2)';
arcCosts(para.nCellRep+1 : 2*para.nCellRep) = c_enter;

arcTails(2*para.nCellRep+1 : 3*para.nCellRep) = (2:2:nNode-2)';
arcHeads(2*para.nCellRep+1 : 3*para.nCellRep) = sink;
arcCosts(2*para.nCellRep+1 : 3*para.nCellRep) = c_exit;



%%%  P_xij
tempTails = zeros(nNode*20, 1);
tempHeads = zeros(nNode*20, 1);
tempCosts = zeros(nNode*20, 1);

countEdge = 0;

% disp(cellSubIdxsLUT(3517,:));
% disp(cellSubIdxsLUT(3568,:));
% disp(cellIdxsLst{t1}(49));
% disp(cellIdxsLst{t2}(56));

for t1 = 1:nnT-1
    % t1 = 52;
    disp(['t: ', num2str(t1)]);
    % cellIdMap1 = cellIdMapStack(:,:,t);
    
    nCnt = cellLst{t1}.NumObjects;
    
    for iCnt1 = 1:nCnt
        % iCnt1 = 2;
        pxVec1 = cellLst{t1}.VoxelIdxList{iCnt1};
        pos1 = cellLst{t1}.ctrPt(iCnt1,:);
        coord1 = cellLst{t1}.CoordinateRgs(iCnt1,:);
        intst1 = cellLst{t1}.avgItstyVec(iCnt1);
        vlm1 = cellLst{t1}.areaVec(iCnt1);
        % rpMtLb1 = rpdMtRltdDtcts(cellIdxsLst{t1}(iCnt1),2);
        rpMtP1 = rpdMtRltdDtctPostP1(cellIdxsLst{t1}(iCnt1),2);
        
        for kJump = 0:para.maxJump
            % kJump = 1;
            if (t1 + 1 + kJump) > nnT
                break;
            end
            t2 = t1 + 1 + kJump;
            
            % pCtrPos1 = mvnpdf(cellLst{t2}.ctrPt, pos1, xRpMtPara.DsplcmtSigma*(kJump+1));
            % pCtrPos0 = mvnpdf(cellLst{t2}.ctrPt, pos1, xStaticPara.DsplcmtSigma*(kJump+1));
            candDists = sqrt(sum((cellLst{t2}.ctrPt-ones(cellLst{t2}.NumObjects,1)*pos1).^2,2));
            pCtrPos1 = gampdf(candDists,xRpMtPara.distPhat(1),xRpMtPara.distPhat(2));
            pCtrPos0 = gampdf(candDists,xStaticPara.distPhat(1),xStaticPara.distPhat(2));
            rpMtEdgeP1 = rpMtP1 * rpdMtRltdDtctPostP1(cellIdxsLst{t2},1);
            pCtrPos = (1-rpMtEdgeP1).*pCtrPos0 + rpMtEdgeP1.*pCtrPos1;
            candIds = (1:cellLst{t2}.NumObjects)';
            candIds(pCtrPos < max(para.pLinkThr, para.FPdtctRate/(1-para.FPdtctRate))) = [];
            
            if isempty(candIds)
                continue;
            end
            
            
            %%%%%    Rapid motion related probabilities
            % candCtrPos = cellLst{t2}.ctrPt(candIds,:);
            % candDists = candDists(candIds);
            if kJump == 0
                pCtrPos1 = pCtrPos1(candIds);
                pOvlpRt1 = zeros(length(candIds),1) + pInvOvrlpRtRefStatic;
                pNbIncldRt1 = zeros(length(candIds),1) + pNbIncldRtRefStatic;
                
            else
                pCtrPos1 = zeros(length(candIds),1);
                pOvlpRt1 = zeros(length(candIds),1);
                pNbIncldRt1 = zeros(length(candIds),1);
            end
            
            
            
            
            %%%%%    Rapid motion UNRELATED (static) probabilities
            pCtrPos0 = pCtrPos0(candIds);
            
            %%%  Inherit limit: should be at least partially overlaping
            candCoord = cellLst{t2}.CoordinateRgs(candIds, :);
            cands2rm0 = candCoord(:,2) < coord1(1) | candCoord(:,1) > coord1(2) ...
                | candCoord(:,4) < coord1(3) | candCoord(:,3) > coord1(4) ...
                | candCoord(:,6) < coord1(5) | candCoord(:,5) > coord1(6);
            
            %%%  Key features: overlapRate
            candOvlpRts = zeros(length(candIds),1);
            candOvlps = zeros(length(candIds),1);
            pOvlpRt0 = zeros(length(candIds),1);
            for iiDtct2 = 1:length(candIds)
                % iiDtct2 = 19;
                if cands2rm0(iiDtct2)
                    continue;
                end
                dtctId2 = candIds(iiDtct2);
                pxVec2 = cellLst{t2}.VoxelIdxList{dtctId2};
                % candDtctSod(iiDtct2) = nnz(intersect(pxVec1,pxVec2))/ max(length(pxVec1),length(pxVec2));
                % candDtctSod(iiDtct2) = 1 - nnz(intersect(pxVec1,pxVec2))/ max(length(pxVec1),length(pxVec2));
                ovlp = nnz(intersect(pxVec1,pxVec2));
                candOvlpRts(iiDtct2) = ovlp / (length(pxVec1) + length(pxVec2) - ovlp);
                candOvlps(iiDtct2) = ovlp;
            end
            cands2rm0 = (candOvlpRts == 0);
            pOvlpRt0(~cands2rm0) = exppdf(1-candOvlpRts(~cands2rm0), xStaticPara.invOvrlpRtPhat)...
                * pInvOvrlpRtStaticBase;
            
            
            %%%  Key features: 2nd neighbor inclusion (both sides)
            candVlms = cellLst{t2}.areaVec(candIds);
            candIncldRts = zeros(length(candIds),1);
            pNbIncldRt0 = zeros(length(candIds),1);
            for iiDtct2 = 1:length(candIds)
                % iiDtctTemp = 2;
                if cands2rm0(iiDtct2)
                    continue;
                end
                iCnt2 = candIds(iiDtct2);
                pxVec2 = cellLst{t2}.VoxelIdxList{iCnt2};
                pos2 = cellLst{t2}.ctrPt(iCnt2,:);
                % intst2 = cellLst{t2}.avgItstyVec(iCnt2);
                coord2 = cellLst{t2}.CoordinateRgs(iCnt2,:);
                
                % % %  Check possible 2nd neighbor (positive time order)
                if length(candIds) == 1
                    candIncldRtPosDir = 0;
                else
                    candIncldRtPosDir = candOvlps./candVlms;
                    candIncldRtPosDir(iiDtct2) = 0;
                end
                
                % % %  Check possible 2nd neighbor (negative time order)
                candPctrPos1 = mvnpdf(cellLst{t1}.ctrPt, pos2, xStaticPara.DsplcmtSigma*(kJump+1));
                candIds1 = (1:cellLst{t1}.NumObjects)';
                candIds1(candPctrPos1 <= 0 | candIds1 == iCnt1) = [];
                candCoord1 = cellLst{t1}.CoordinateRgs(candIds1, :);
                cands2rm = candCoord1(:,2) < coord2(1) | candCoord1(:,1) > coord2(2) ...
                    | candCoord1(:,4) < coord2(3) | candCoord1(:,3) > coord2(4) ...
                    | candCoord1(:,6) < coord2(5) | candCoord1(:,5) > coord2(6);
                candIds1(cands2rm) = [];
                
                candIncldRtNegDir = zeros(length(candIds1),1);
                for iiDtct1 = 1:length(candIds1)
                    % iiDtct1 = 19;
                    candDtctId1 = candIds1(iiDtct1);
                    candPxVec1 = cellLst{t1}.VoxelIdxList{candDtctId1};
                    candOvlp = nnz(intersect(candPxVec1, pxVec2));
                    candIncldRtNegDir(iiDtct1) = candOvlp / length(candPxVec1);
                end
                candIds1(candIncldRtNegDir == 0) = [];
                candIncldRtNegDir(candIncldRtNegDir == 0) = [];
                
                candIncldRts(iiDtct2) = max([candIncldRtPosDir; candIncldRtNegDir]);
                % pNbIncldRt = exppdf(nbIncldRt, xNOPC.nbIncldRtPhat);
            end
            pNbIncldRt0(~cands2rm0) = exppdf(candIncldRts(~cands2rm0), xStaticPara.nbIncldRtPhat)...
                * pNbIncldRtStaticBase;
            
            
            
            %%%%%   Combined motion model
            rpMtEdgeP1 = rpMtP1 * rpdMtRltdDtctPostP1(cellIdxsLst{t2}(candIds),1);
            pSpace0 = pCtrPos0 .* pOvlpRt0 .* pNbIncldRt0;
            pSpace1 = pCtrPos1 .* pOvlpRt1 .* pNbIncldRt1;
            pSpace = (1-rpMtEdgeP1).*pSpace0 + rpMtEdgeP1.*pSpace1;

            
                        
            %%%  Appearance features: intensity, sizeChangeRate
            candDtctIntst = cellLst{t2}.avgItstyVec(candIds);
            pIntst = normpdf(candDtctIntst, intst1, aPara.intDiff_sigma);
            
            candDtctVlm = cellLst{t2}.areaVec(candIds);
            pSize = normpdf((candDtctVlm - vlm1)./((candDtctVlm + vlm1)/2) ,...
                aPara.sizeChangeRt_mu, aPara.sizeChangeRt_sigma);
            
            
            %%%  Time feature
            pT = exppdf(1+kJump,(para.maxJump+1)/2);
            
            %%%%%  Overall score
            % pLink = pIntst.*pSpace.*pT;
            pLink = pIntst.*pSize.*pSpace.*pT;
            % pLink = pIntst.^3.*pSize.*pSpace.*pT;
            % pLink = pIntst.*pSize.*pCtrPos.*pOvlpRt.*pNbIncldRt.*pT;
            % pLink = pIntst.*pSize.*pCtrPos.*pOvlpRt.*pNbIncldRt;
            % pLink = pIntst.*pCtrPos.*pOvlpRt.*pNbIncldRt.*pT;
            % boolActvLink = pLink > para.pLinkThr;
            boolActvLink = pLink >= max(para.pLinkThr, para.FPdtctRate/(1-para.FPdtctRate));
            
            
            
%             %%%%%%%  For debugging  %%%%%%%
%             % boolActvLink(candIds==92) = true;
%             disp([length(boolActvLink), nnz(boolActvLink)]);
%             showMat = [candIds(boolActvLink), pIntst(boolActvLink), pSize(boolActvLink),...
%                 pSpace0(boolActvLink),pSpace1(boolActvLink),pSpace(boolActvLink),...
%                 pT+zeros(nnz(boolActvLink),1), pLink(boolActvLink), -log(pLink(boolActvLink))];
%             % showIdx = (1:length(candIds))';
%             % showIdx = [find(candIds == 2), find(candIds == 61)];
%             % disp(showMat(showIdx,:));
%             disp('[candIds,   pIntst,   pSize,   pSpace0,   pSpace1,   pSpace,   pT,   pLink,   -log(pLink)]');
%             disp(showMat);
%             disp('[candIds,   rpMtP1,   rpMtP2,  rpMtEdgeP1,  pCtrPos0, pCtrPos1, pOvlpRt0, pOvlpRt1, pNbIncldRt0, pNbIncldRt1]')
%             disp([candIds(boolActvLink), rpMtP1+zeros(nnz(boolActvLink),1),...
%                 rpdMtRltdDtctPostP1(cellIdxsLst{t2}(candIds(boolActvLink)),1),...
%                 rpMtEdgeP1(boolActvLink), pCtrPos0(boolActvLink), pCtrPos1(boolActvLink),...
%                 pOvlpRt0(boolActvLink), pOvlpRt1(boolActvLink),...
%                 pNbIncldRt0(boolActvLink), pNbIncldRt1(boolActvLink)]);
%             % disp(candDists(candIds(boolActvLink)));
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            nNewEdge = nnz(boolActvLink);
            tempTails(countEdge+1: countEdge+nNewEdge) = cellIdxsLst{t1}(iCnt1) * 2;
            tempHeads(countEdge+1: countEdge+nNewEdge) = cellIdxsLst{t2}(candIds(boolActvLink))*2 - 1;
            tempCosts(countEdge+1: countEdge+nNewEdge) = -log(pLink(boolActvLink));
            % tempCosts(countEdge+1: countEdge+nNewEdge) = -log(candOvrlpRates * pT);
            countEdge = countEdge + nNewEdge;
        end
    end
end
% figure; histogram(tempCosts(1:countEdge),100);

arcTails = [arcTails; tempTails(1:countEdge)];
arcHeads = [arcHeads; tempHeads(1:countEdge)];
arcCosts = [arcCosts; tempCosts(1:countEdge)];

arcCosts = arcCosts + 1e-5*randn(length(arcCosts),1);

% clear tempTails; clear tempHeads; clear tempCosts;

trackG = digraph(arcTails,arcHeads,arcCosts);

disp('Graph building for linking is finished!');
disp(['Average links pointed out from one detection: ',num2str(countEdge / para.nCellRep)]);



