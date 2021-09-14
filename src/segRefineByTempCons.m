function cellLst_refined = segRefineByTempCons(cellLst, imgStack, invPCmapStack, nullMtScPhat, effZinFmTwoEndIdx)


nDtct = sum(cellfun(@(x) x.NumObjects, cellLst));


%%  Iterative refinement: outer loop

cellLst_cur = cellLst;
isChanged = true;
numCellResVec = zeros(nDtct,1) + nDtct;
iOutIter = 0;
while any(numCellResVec > 1) && isChanged
    iOutIter = iOutIter + 1;
    disp(['Outer loop #', num2str(iOutIter), ' ...']);
    %%%%   Infer the number of cells in each detection
    [numCellResVec, dtctIncldLst] = tempConsistSort(cellLst_cur, nullMtScPhat);
    
    %%%% Segmentation refinement operation according to the inference
    [cellLst_cur, isChanged] = segRefineOperations(cellLst_cur, numCellResVec, dtctIncldLst,...
        imgStack, invPCmapStack, effZinFmTwoEndIdx);
end

cellLst_refined = cellLst_cur;


%%  Visualize the results
% nDtct_res = sum(cellfun(@(x) x.NumObjects, cellLst_cur));
% cellIdxsLst_res = cellfun(@(x) (1:x.NumObjects)', cellLst_cur, 'UniformOutput', false);
% for t = 2:nnT
%     cellIdxsLst_res{t}(:,1) = cellIdxsLst_res{t}(:,1) + cellIdxsLst_res{t-1}(end,1);
% end
% cellSubIdxsLUT_res = cell(nnT,1);
% for t = 1:nnT
%     cellSubIdxsLUT_res{t} = [t + zeros(cellLst_cur{t}.NumObjects,1), (1:cellLst_cur{t}.NumObjects)'];
% end
% cellSubIdxsLUT_res = cell2mat(cellSubIdxsLUT_res);
% 
% 
% 
% %%%  All detections
% traceVxLstLst = cell(nDtct_res, 1);
% for iDtct = 1:nDtct_res
%     traceVxLstLst{iDtct} = cell(nnT,1);
%     tt = cellSubIdxsLUT_res(iDtct, 1);
%     iDtctSub = cellSubIdxsLUT_res(iDtct, 2);
%     traceVxLstLst{iDtct}{tt} = cellLst_cur{tt}.VoxelIdxList{iDtctSub};
% end
% 
% % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, [], false, []);
% % % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, (cellIdPatchStack>0)*1, false, []);
% [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, false, []);
% % [printImgStack,printImgYZStack,~] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, true, cmapLb);
% % implay(printImgStack);
% % implay(printImgYZStack);



