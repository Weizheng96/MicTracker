function [childDtctSet, ovlpSc] = getChildSet(dtct0, tgtDtctLst, flagTempFw, nullMtScPhat)

% flagTempFw = 0;
% tgtDtctLst = cellLst{t0-1};
% flagTempFw = 1;
% tgtDtctLst = cellLst{t0+1};

childDtctSet = [];
ovlpSc = 0;


pMtScThr = 0.005;
if (flagTempFw == 1 && gamcdf(dtct0.mtSc(2), nullMtScPhat(1), nullMtScPhat(2), 'upper') < pMtScThr) ||...
        (flagTempFw == 0 && gamcdf(dtct0.mtSc(1), nullMtScPhat(1), nullMtScPhat(2), 'upper') < pMtScThr)
    return;
end


candDtctIds = find(cellfun(@(x) ~isempty(intersect(dtct0.vxLst, x)), tgtDtctLst.VoxelIdxList));
% candDtctIds = unique(tgtDtctIdMap(dtct0.vxLst));
candDtctIds(candDtctIds == 0) = [];

if isempty(candDtctIds)
    return;
end


%%%%   Find all o2m (/o2o) sets with dtct0 as the "one" side

vxVec0 = dtct0.vxLst;
vlm0 = length(vxVec0);


%%%   Check the best members ('multi' side)

candDtctVlm = tgtDtctLst.areaVec(candDtctIds);
candDtctOvlp = cellfun(@(x) nnz(intersect(vxVec0,x)), tgtDtctLst.VoxelIdxList(candDtctIds));
% candDtctUnion = candDtctVlm + vlm0 - candDtctOvlp;

%%%   Start from the full candidate neighbor set, reduce the set
rmCandLbVec = false(length(candDtctIds), 1);
ovlpSc = sum(candDtctOvlp) / (vlm0 + sum(candDtctVlm) - sum(candDtctOvlp));

for iNcand = 1:length(candDtctIds)
    tempOvlpScVec = zeros(length(candDtctIds),1);
    for iiCand = 1:length(candDtctIds)
        if rmCandLbVec(iiCand)
            continue;
        end
        tryRmCandLbVec = rmCandLbVec;
        tryRmCandLbVec(iiCand) = true;
        temp = sum(candDtctOvlp(~tryRmCandLbVec));
        tempOvlpScVec(iiCand) = temp / (vlm0 + sum(candDtctVlm(~tryRmCandLbVec)) - temp);
    end
    [tempOvlpSc,iBestCand] = max(tempOvlpScVec);
    if tempOvlpSc > ovlpSc
        rmCandLbVec(iBestCand) = true;
        ovlpSc = tempOvlpSc;
    else
        break;
    end
end
childDtctSet = candDtctIds(~rmCandLbVec);
% childDtctVlms = candDtctVlm(~rmCandLbVec);
% childDtctOvlps = candDtctOvlp(~rmCandLbVec);
% childDtctUnions = childDtctVlms + vlm1 - childDtctOvlps;

if ovlpSc < 0.2
    childDtctSet = [];
    % disp(ovlpSc);
end


if flagTempFw == 1
    pMtScVec = gamcdf(tgtDtctLst.dtctMtScVec(childDtctSet,1), nullMtScPhat(1), nullMtScPhat(2), 'upper');
else
    pMtScVec = gamcdf(tgtDtctLst.dtctMtScVec(childDtctSet,2), nullMtScPhat(1), nullMtScPhat(2), 'upper');
end
if any(pMtScVec < pMtScThr)
    childDtctSet = [];
end


% tgtVxs = cell2mat(tgtDtctLst.VoxelIdxList(childDtctSet));
% tgtApprVals = tgtApprMap(tgtVxs);
% pAppr = 2*cdf('norm', abs(median(tgtApprVals) - median(dtct0.apprValVec)), aPara.mu, aPara.sigma, 'upper');
% if pAppr < 0.05
%     childDtctSet = [];
% end

% figure; histogram(dtct0.apprValVec, 30); xlim([0,1]);
% hold on; histogram(tgtApprVals, 30);









