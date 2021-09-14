function [minCC, optSpplVal, resFlowVec, spplValnCostVec] = myMCCbyCMF(arcVec, nnNode, nnArc, source, sink, upbndSpplVal)

if ~exist('cs2_mex.mexw64','file')
    mex -setup;
    mex cs2_mex.c;
end


costMat = int64(zeros(nnNode));
costMat(sub2ind([nnNode,nnNode], arcVec(:,1), arcVec(:,2))) = arcVec(:,5);


% spplVal_l = 0;
% spplVal_u = nAllDtct;
spplValIdx_l = 1;
spplValIdx_u = 2;


spplValnCostVec = nan(2, upbndSpplVal);
spplValnCostVec(1,1:2) = [0, upbndSpplVal];
grdDecDirRight = nan(1, upbndSpplVal);

for iSppl = 1 : size(spplValnCostVec,2)
% for iSppl = 1 : 4
    % iSppl = 5;
    % spplVal = round(nAllDtct/nnT);
    
    spplVal_l = spplValnCostVec(1, spplValIdx_l);
    spplVal_u = spplValnCostVec(1, spplValIdx_u);

    if iSppl == 1
        spplVal = spplVal_l;
        costAllVal = 0;
        spplValnCostVec(:, iSppl) = [spplVal, costAllVal]';
        continue;
    elseif iSppl == 2
        spplVal = spplVal_u;
    elseif spplVal_l == spplVal_u - 1
        break;
    elseif spplValnCostVec(2, spplValIdx_l) == spplValnCostVec(2, spplValIdx_u)
        spplVal = spplVal_l + 1;
    else
        spplVal = round(0.5*(spplVal_l+spplVal_u));
    end
    spplVal = int64(spplVal);
    
    
    %%%%  The node supply list  (format: [<nodeId>, <supply>])
    % nodeVec = [];
    nodeVec = [source, spplVal; sink, -spplVal];
    
    
    %%%%  Solve the MCF given supply
    resFlowVec = myMCFportal(arcVec, nodeVec, nnNode, nnArc);
    
    
    %%%%  Get the objective function value (overall cost)
    resFlowMat = int64(zeros(nnNode));
    resFlowMat(sub2ind([nnNode,nnNode], resFlowVec(:,1), resFlowVec(:,2))) = resFlowVec(:,3);
    costAllVal = sum(sum(resFlowMat.*costMat));
    
    
    %%%%  Get the "gradient" direction
    nodeVec = [source, spplVal+1; sink, -(spplVal+1)];
    resFlowVec = myMCFportal(arcVec, nodeVec, nnNode, nnArc);
    resFlowMat = int64(zeros(nnNode));
    resFlowMat(sub2ind([nnNode,nnNode], resFlowVec(:,1), resFlowVec(:,2))) = resFlowVec(:,3);
    costAllVal_rnb = sum(sum(resFlowMat.*costMat));
    if costAllVal_rnb < costAllVal
        grdDecDirRight(iSppl) = true;
        spplValnCostVec(:, iSppl) = [spplVal+1, costAllVal_rnb]';
    else
        grdDecDirRight(iSppl) = false;
        spplValnCostVec(:, iSppl) = [spplVal, costAllVal]';
    end
    
    
    if iSppl <= 2
        continue;
    elseif costAllVal > max(spplValnCostVec(2, [spplValIdx_l,spplValIdx_u]))
        break;
    else
        if grdDecDirRight(iSppl)
            spplValIdx_l = iSppl;
        else
            spplValIdx_u = iSppl;
        end
    end
    
end

spplValnCostVec(:, isnan(spplValnCostVec(1,:))) = [];

[minCC, tempI] = min(spplValnCostVec(2,:));
optSpplVal = spplValnCostVec(1,tempI);


%%%%  Get the final results with optimal supply value
nodeVec = [source, optSpplVal; sink, -optSpplVal];
resFlowVec = myMCFportal(arcVec, nodeVec, nnNode, nnArc);





