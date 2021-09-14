function traceLst = traceCreate_3D(traceDtctIdSet, cellLst, imgSzNrg)


nTrace = length(traceDtctIdSet);

nnY = imgSzNrg.size(1);
nnX = imgSzNrg.size(2);
nnZ = imgSzNrg.size(3);
nnT = imgSzNrg.size(4);

traceLst = cell(nTrace,1);

for iTrace = 1:nTrace
    % iTrace = 284;

    % t = 92;
    % cellIdPatch1 = cellIdPatchStack(:,:,:,t);
    % % myImageStackPrint(cellIdPatch1, 10, true);
    % iSubDtct1 = 287;
    % iTrace = find(cellfun(@(x) any(x(:,1)==t & x(:,2)==iSubDtct1),traceDtctIdSet));

    disp(['Interpolating iTrace: ', num2str(iTrace), ' in total ', num2str(nTrace)]);

    %%%%  Basic info of the selected trace
    trace0.dtctSubs = traceDtctIdSet{iTrace};
    trace0.stT = trace0.dtctSubs(1,1);
    trace0.endT = trace0.dtctSubs(end,1);
    trace0.isDtctd = false(nnT,1);
    trace0.isDtctd(trace0.dtctSubs(:,1)) = true;

    
    %%%%  Voxel lists of the selected trace (deteced part)
    trace0.traceVxLst = cell(nnT,1);
    temp = mat2cell(trace0.dtctSubs, ones(size(trace0.dtctSubs,1),1), 2);
    trace0.traceVxLst(trace0.dtctSubs(:,1)) = cellfun(@(x) cellLst{x(1)}.VoxelIdxList{x(2)}, temp, 'UniformOutput', false);

    
    traceLst{iTrace} = trace0;

end


%%%%%   Pre-calculate features of traces: spatial position
for iTrace = 1:nTrace
    % iTrace = 36;
    %     [rR,cC] = ind2sub([nRow,nCol], unique(cell2mat(intpltdTraceLst{iTrace}.traceVxLst)));
    %     intpltdTraceLst{iTrace}.posBounds = [min(rR), max(rR), min(cC), max(cC)];
    traceLst{iTrace}.posBoundsPerT = nan(nnT,6);
    traceLst{iTrace}.cntrPtSub = nan(nnT, 3);
    for t = traceLst{iTrace}.stT : traceLst{iTrace}.endT
        if traceLst{iTrace}.isDtctd(t)
            [yY,xX,zZ] = ind2sub([nnY,nnX,nnZ], traceLst{iTrace}.traceVxLst{t});
            traceLst{iTrace}.posBoundsPerT(t,:) = [min(yY), max(yY), min(xX), max(xX), min(zZ), max(zZ)];
            traceLst{iTrace}.cntrPtSub(t,:) = mean([yY,xX,zZ]);
        end
    end
    temp = [min(traceLst{iTrace}.posBoundsPerT); max(traceLst{iTrace}.posBoundsPerT)];
    traceLst{iTrace}.posBounds = temp(logical([1,0,1,0,1,0; 0,1,0,1,0,1]))';
end

