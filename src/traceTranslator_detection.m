function traceDtctIdSet = traceTranslator_detection(trajectories, cellSubIdxsLUT, cellIdxsLst, isIncldFPs)

traceDtctIdSet = cellfun(@(x) cellSubIdxsLUT(x,:), trajectories, 'UniformOutput', false);

if isIncldFPs
    cellIncldLbLst = cellfun(@(x) false(size(x)), cellIdxsLst, 'UniformOutput', false);
    for jTrc = 1:length(traceDtctIdSet)
        for jDtct = 1:size(traceDtctIdSet{jTrc},1)
            cellIncldLbLst{traceDtctIdSet{jTrc}(jDtct,1)}(traceDtctIdSet{jTrc}(jDtct,2)) = true;
        end
    end
    % sum(cellfun(@nnz, cellIncldLbLst))

    FPsToBeIncldd = cellfun(@(x) find(~x), cellIncldLbLst, 'UniformOutput', false);
    for t = 1:length(FPsToBeIncldd)
        FPsToBeIncldd{t} = [t + zeros(size(FPsToBeIncldd{t})), FPsToBeIncldd{t}];
    end
    FPsToBeIncldd = cell2mat(FPsToBeIncldd);
    FPsToBeIncldd = mat2cell(FPsToBeIncldd, ones(size(FPsToBeIncldd,1),1), 2);
    traceDtctIdSet = [traceDtctIdSet; FPsToBeIncldd];
end
