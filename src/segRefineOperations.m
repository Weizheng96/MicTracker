function [cellLst_res, isChanged] = segRefineOperations(cellLst_cur, ...
    numCellResVec, dtctIncldLst, imgStack, invPCmapStack, effZinFmTwoEndIdx)

% %   Segmentation refinement operation according to the inference


[nnY, nnX, nnZ, nnT] = size(imgStack);

imgSzNrg.size = [nnY,nnX,nnZ,nnT];
imgSzNrg.range = [1,nnY; 1,nnX; 1,nnZ; 1,nnT];



%%  Prepare: current segmentation map  &  index LUTs

nDtct = sum(cellfun(@(x) x.NumObjects, cellLst_cur));

%%% The overall cell representer indices
% cellIdxsLst = cellfun(@(x) (1:x.NumObjects)', cellLst_cur, 'UniformOutput', false);
% for t = 2:nnT
%     cellIdxsLst{t}(:,1) = cellIdxsLst{t}(:,1) + cellIdxsLst{t-1}(end,1);
% end

cellSubIdxsLUT = cell(nnT,1);
for t = 1:nnT
    cellSubIdxsLUT{t} = [t + zeros(cellLst_cur{t}.NumObjects,1), (1:cellLst_cur{t}.NumObjects)'];
end
cellSubIdxsLUT = cell2mat(cellSubIdxsLUT);




%%  Iterative operations

xRes = numCellResVec;


%%%%%    Refinement operation: splitting    %%%%%
dtcts2Split = find(xRes >= 2);
revIdLUT = nan(nDtct,1);
revIdLUT(dtcts2Split) = 1:length(dtcts2Split);

newCellVxLst = cell(length(dtcts2Split), 1);


%%%  Find out the detections of interest that initially have valid reference of splitting
vldRefDirFlags = false(length(dtcts2Split),2);
% potentVldRefDirFlags = false(length(dtcts2Split),2);
for ii = 1:length(dtcts2Split)
    % ii = 6;
    % ii = find(dtcts2Split == 198);
    iDtct = dtcts2Split(ii);
    t0 = cellSubIdxsLUT(iDtct, 1);
    numCellInDtct0 = xRes(iDtct);
    
    refCellNums = [sum(xRes(dtctIncldLst{iDtct}.lastFm)), sum(xRes(dtctIncldLst{iDtct}.nextFm))];
    
    potentVldRefDirFlag = [all(xRes(dtctIncldLst{iDtct}.lastFm) == 1), ...
        all(xRes(dtctIncldLst{iDtct}.nextFm) == 1)];
    tempVec = potentVldRefDirFlag & (refCellNums == numCellInDtct0);
    if any(tempVec)
        vldRefDirFlags(ii, :) = tempVec;
    else
        vldRefDirFlags(ii, :) = potentVldRefDirFlag & (refCellNums > 1);
    end
        
    %     % potentVldRefDirFlags(ii, :) = refCellNums == numCellInDtct0;
    %     potentVldRefDirFlags(ii, :) = refCellNums >= numCellInDtct0;
    %     if potentVldRefDirFlags(ii, 1) && all(xRes(dtctIncldLst{iDtct}.lastFm) == 1)
    %         vldRefDirFlags(ii, 1) = true;
    %     end
    %     if potentVldRefDirFlags(ii, 2) && all(xRes(dtctIncldLst{iDtct}.nextFm) == 1)
    %         vldRefDirFlags(ii, 2) = true;
    %     end
end

% %%%%%%%%%%%%%%   Debug   %%%%%%%%%%%%%%%
% disp(find(sum(vldRefDirFlags,2)>0 & ismember(dtcts2Split, dtctOI(tempFlagVec))));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic
toProcFlagVec = true(length(dtcts2Split),1);
tempMap = false(nnY,nnX,nnZ);
tempMap2 = zeros(nnY,nnX,nnZ);

oldNnz1 = -1;
newCellNumVec = cellfun(@length, newCellVxLst);
% toFrthrProcDtcts = find((newCellNumVec>0) & (newCellNumVec~=xRes(dtcts2Split)));
toFrthrProcDtcts = find(newCellNumVec~=xRes(dtcts2Split));

% %%%%%%%%%%%%%%%%%%%%%%%   Debug   %%%%%%%%%%%%%%%%%%%%%%%%
% count = 0;
% rcdVec = oldNnz1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while sum(abs(newCellNumVec-xRes(dtcts2Split))) ~= oldNnz1
    % % %  Outer loop: to ensure utilizing the best direction (s.t. not influenced by sequence of processing)
    oldNnz1 = sum(abs(newCellNumVec-xRes(dtcts2Split)));
    toProcFlagVec(toFrthrProcDtcts) = true;
    
%     %%%%%%%%%%%%%%%%%%%%%%%   Debug   %%%%%%%%%%%%%%%%%%%%%%%%
%     count = count + 1;
%     rcdVec(end+1) = oldNnz1;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oldNnz2 = -1;
    while any(toProcFlagVec) && nnz(toProcFlagVec)~=oldNnz2
        oldNnz2 = nnz(toProcFlagVec);
        
        for ii = 1:length(dtcts2Split)
            %     for ii = 1:80
            % ii = 82;
            % ii = find(dtcts2Split == 1500);
            if ~toProcFlagVec(ii) || ~any(vldRefDirFlags(ii,:))
                continue;
            end
            iDtct = dtcts2Split(ii);
            % disp(['Splitting... Detection #', num2str(iDtct)]);
            
            t0 = cellSubIdxsLUT(iDtct, 1);
            dtctSub0 = cellSubIdxsLUT(iDtct, 2);
            numCellInDtct0 = xRes(iDtct);
            
            scMap = 1 - invPCmapStack(:,:,:,t0);
            
            tempMap(:) = false;
            tempMap(cellLst_cur{t0}.VoxelIdxList{dtctSub0}) = true;
            roiBnMap = tempMap(imgSzNrg.range(1,1):imgSzNrg.range(1,2), ...
                imgSzNrg.range(2,1):imgSzNrg.range(2,2), imgSzNrg.range(3,1):imgSzNrg.range(3,2));
            
            %%%  Find reference
            % childrenVxLst = vldRefLst{ii};
            candRefDir = vldRefDirFlags(ii, :);
            candRefNum = [sum(xRes(dtctIncldLst{iDtct}.lastFm)), sum(xRes(dtctIncldLst{iDtct}.nextFm))];
            % candRefNum(candRefNum <= 1) = Inf;
            candRefNum(~candRefDir) = Inf;
            [~, refDir] = min(abs(xRes(iDtct) - candRefNum));
            if candRefNum(1) == candRefNum(2)
                revIds = revIdLUT(dtctIncldLst{iDtct}.lastFm);
                temp = newCellVxLst(revIds(~isnan(revIds)));
                temp2 = cellLst_cur{t0-1}.VoxelIdxList(cellSubIdxsLUT(dtctIncldLst{iDtct}.lastFm(isnan(revIds)), 2));
                childrenVxLst1 = [cat(1, temp{:}); temp2];
                vxVec1 = cellfun(@length, childrenVxLst1);
                
                revIds = revIdLUT(dtctIncldLst{iDtct}.nextFm);
                temp = newCellVxLst(revIds(~isnan(revIds)));
                temp2 = cellLst_cur{t0+1}.VoxelIdxList(cellSubIdxsLUT(dtctIncldLst{iDtct}.nextFm(isnan(revIds)), 2));
                childrenVxLst2 = [cat(1, temp{:}); temp2];
                vxVec2 = cellfun(@length, childrenVxLst2);
                
                if length(childrenVxLst1) == length(childrenVxLst2)
                    if min(vxVec1)/max(vxVec1) >= min(vxVec2)/max(vxVec2)
                        childrenVxLst = childrenVxLst1;
                    else
                        childrenVxLst = childrenVxLst2;
                    end
                else
                    if abs(xRes(iDtct) - length(childrenVxLst1)) <= abs(xRes(iDtct) - length(childrenVxLst2))
                        childrenVxLst = childrenVxLst1;
                    else
                        childrenVxLst = childrenVxLst2;
                    end
                end
            else
                if refDir == 1
                    revIds = revIdLUT(dtctIncldLst{iDtct}.lastFm);
                    temp = newCellVxLst(revIds(~isnan(revIds)));
                    temp2 = cellLst_cur{t0-1}.VoxelIdxList(cellSubIdxsLUT(dtctIncldLst{iDtct}.lastFm(isnan(revIds)), 2));
                    childrenVxLst = [cat(1, temp{:}); temp2];
                elseif refDir == 2
                    revIds = revIdLUT(dtctIncldLst{iDtct}.nextFm);
                    temp = newCellVxLst(revIds(~isnan(revIds)));
                    temp2 = cellLst_cur{t0+1}.VoxelIdxList(cellSubIdxsLUT(dtctIncldLst{iDtct}.nextFm(isnan(revIds)), 2));
                    childrenVxLst = [cat(1, temp{:}); temp2];
                end
            end
            if length(childrenVxLst) > xRes(iDtct)
                refDtctVxLst = cellfun(@(x) intersect(x, cellLst_cur{t0}.VoxelIdxList{dtctSub0}), childrenVxLst, 'UniformOutput', false);
                [~,bb] = sort(cellfun(@length, refDtctVxLst), 'descend');
                childrenVxLst = childrenVxLst(bb(1:xRes(iDtct)));
                % disp(['Not exact reference... Detection #', num2str(iDtct)]);
            end
            
            % refDtctVxLst = cellLst{cellSubIdxsLUT(children(1), 1)}.VoxelIdxList(cellSubIdxsLUT(children, 2));
            refVlmVec = cellfun(@length, childrenVxLst);
            refDtctVxLst = cellfun(@(x) intersect(x, cellLst_cur{t0}.VoxelIdxList{dtctSub0}), childrenVxLst, 'UniformOutput', false);
            if any(cellfun(@isempty, refDtctVxLst))
                % disp(['Not exact reference... Detection #', num2str(iDtct)]);
                refVlmVec = refVlmVec(~cellfun(@isempty, refDtctVxLst));
                refDtctVxLst = refDtctVxLst(~cellfun(@isempty, refDtctVxLst));
                if length(refDtctVxLst) <= 1
                    disp(['Error (Detection #',num2str(iDtct),'): lack of reference for splitting operation!']);
                    continue;
                end
            end
            seedVxLst = cell(length(refDtctVxLst), 1);
            for iChild = 1:length(refDtctVxLst)
                [yy,xx,zz] = ind2sub([nnY,nnX,nnZ], refDtctVxLst{iChild});
                seedVxLst{iChild} = sub2ind([nnY,nnX,nnZ], round(mean(yy)), round(mean(xx)), round(mean(zz)));
            end
            tempMap(:) = false;
            tempMap(cell2mat(seedVxLst)) = true;
            % tempMap(cell2mat(refDtctVxLst)) = true;
            seedLbMap = tempMap(imgSzNrg.range(1,1):imgSzNrg.range(1,2), ...
                imgSzNrg.range(2,1):imgSzNrg.range(2,2), imgSzNrg.range(3,1):imgSzNrg.range(3,2));
            
            subRoiLbMap = splitROIbyWtrshd3D(roiBnMap, seedLbMap, scMap);
            % myImageStackPrint(roiBnMap(:,:,1:3:nnZ), 10, true);
            % myImageStackPrint(seedLbMap(:,:,1:3:nnZ), 10, true);
            % myImageStackPrint(subRoiLbMap(:,:,1:3:nnZ), 10, true);
            
            
            %%%  Double check the result of splitting operation: round I
            redoFlag = false;
            if any(subRoiLbMap(seedLbMap) == 0)
                redoFlag = true;
                % disp(['Error (Detection #',num2str(iDtct),'): unreliable Watershed!']);
            else
                tempMap2(:) = false;
                tempMap2(imgSzNrg.range(1,1):imgSzNrg.range(1,2), imgSzNrg.range(2,1):imgSzNrg.range(2,2), ...
                    imgSzNrg.range(3,1):imgSzNrg.range(3,2)) = subRoiLbMap;
                subRoiVxLst = cell(length(refDtctVxLst), 1);
                for iChild = 1:length(refDtctVxLst)
                    % subRoiVxLst{iChild} = find(tempMap2 == iChild);
                    subRoiVxLst{iChild} = find(tempMap2 == tempMap2(seedVxLst{iChild}));
                end
                
                %%%  Double-check the result of Watershed
                subRoiVec = cellfun(@length, subRoiVxLst);
                tempVec = (subRoiVec/sum(subRoiVec)) ./ (refVlmVec/sum(refVlmVec));
                minDtctSize = 20;
                if any(subRoiVec < minDtctSize | subRoiVec > length(cellLst_cur{t0}.VoxelIdxList{dtctSub0})) ...
                        || min(tempVec) < 0.25
                    redoFlag = true;
                end
            end
            if redoFlag
                tempMap(:) = false;
                % tempMap(cell2mat(seedVxLst)) = true;
                tempMap(cell2mat(refDtctVxLst)) = true;
                seedLbMap = tempMap(imgSzNrg.range(1,1):imgSzNrg.range(1,2), ...
                    imgSzNrg.range(2,1):imgSzNrg.range(2,2), imgSzNrg.range(3,1):imgSzNrg.range(3,2));
                
                subRoiLbMap = splitROIbyWtrshd3D(roiBnMap, seedLbMap, scMap);
                % myImageStackPrint(subRoiLbMap(:,:,1:3:nnZ), 10, true);
                
                
                %%%  Double check the result of splitting operation: round II
                redoFlag2 = false;
                if any(subRoiLbMap(seedLbMap) == 0)
                    redoFlag2 = true;
                    % disp(['Error (Detection #',num2str(iDtct),'): unreliable Watershed!']);
                else
                    tempMap2(:) = false;
                    tempMap2(imgSzNrg.range(1,1):imgSzNrg.range(1,2), imgSzNrg.range(2,1):imgSzNrg.range(2,2), ...
                        imgSzNrg.range(3,1):imgSzNrg.range(3,2)) = subRoiLbMap;
                    subRoiVxLst = cell(length(refDtctVxLst), 1);
                    for iChild = 1:length(refDtctVxLst)
                        tempVec = find(tempMap2 == iChild);
                        [~,tempI] = max(cellfun(@(x) length(intersect(x,tempVec)), refDtctVxLst));
                        subRoiVxLst{tempI} = tempVec;
                        % subRoiVxLst{iChild} = find(tempMap2 == tempMap2(seedVxLst{iChild}));
                    end
                    
                    %%%  Double-check the result of Watershed
                    subRoiVec = cellfun(@length, subRoiVxLst);
                    tempVec = (subRoiVec/sum(subRoiVec)) ./ (refVlmVec/sum(refVlmVec));
                    minDtctSize = 20;
                    if any(subRoiVec < minDtctSize | subRoiVec > length(cellLst_cur{t0}.VoxelIdxList{dtctSub0})) ...
                            || min(tempVec) < 0.25
                        redoFlag2 = true;
                    end
                end
                if redoFlag2
                    % subRoiVxLst = refDtctVxLst;
                    subRoiVxLst = [];
                end
            end
            
            
            %%%  Record the result
            if isempty(subRoiVxLst)
                toProcFlagVec(ii) = false;
            else
                toProcFlagVec(ii) = false;
                newCellVxLst{ii} = subRoiVxLst;
                % newCellVxLst{t0}(end+1: end+length(subRoiVxLst)) = subRoiVxLst;
                
                % % %  Currently processed detection contributes to the "valid reference"
                nbDtcts = [dtctIncldLst{iDtct}.lastFm; dtctIncldLst{iDtct}.nextFm];
                nbDtcts = nbDtcts(~isnan(revIdLUT(nbDtcts)));
                for iNb = 1:length(nbDtcts)
                    nbDtct = nbDtcts(iNb);
                    nbLcIdx = revIdLUT(nbDtct);
                    if ~vldRefDirFlags(nbLcIdx, 1) && sum(xRes(dtctIncldLst{nbDtct}.lastFm)) > 1
                        revIds = revIdLUT(dtctIncldLst{nbDtct}.lastFm);
                        revIds = revIds(~isnan(revIds));
                        vldRefDirFlags(nbLcIdx, 1) = all(~cellfun(@isempty, newCellVxLst(revIds)));
                    end
                    if ~vldRefDirFlags(nbLcIdx, 2) && sum(xRes(dtctIncldLst{nbDtct}.nextFm)) > 1
                        revIds = revIdLUT(dtctIncldLst{nbDtct}.nextFm);
                        revIds = revIds(~isnan(revIds));
                        vldRefDirFlags(nbLcIdx, 2) = all(~cellfun(@isempty, newCellVxLst(revIds)));
                    end
                end
                %             %%%%%%%%%%%%%   Debug   %%%%%%%%%%%%%
                %             if any(vldRefDirFlags(108,:))
                %                 return;
                %             end
                %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end
    
    newCellNumVec = cellfun(@length, newCellVxLst);
    % toFrthrProcDtcts = find((newCellNumVec>0) & (newCellNumVec~=xRes(dtcts2Split)));
    toFrthrProcDtcts = find(newCellNumVec~=xRes(dtcts2Split));
end
toc





%%%%%    Refinement operation: merging (or neglecting) FPs    %%%%%
% tic
% 
% dtcts2Merge = find(xRes == 0);
% mergedFlagVec = false(length(dtcts2Merge),1);
% 
% for ii = 1:length(dtcts2Merge)
%     iDtct = dtcts2Merge(ii);
%     t0 = cellSubIdxsLUT(iDtct, 1);
%     iDtctSub0 = cellSubIdxsLUT(iDtct, 2);
%     % numCellInDtct0 = xRes(iDtct);
%     
%     nbDtctIds = [predecessors(GdtctIncld,iDtct); successors(GdtctIncld,iDtct)];
%     nbDtctIds(xRes(nbDtctIds) ~= 1) = [];
%     for jj = 1:length(nbDtctIds)
%         jDtct = nbDtctIds(jj);
%         if length(dtctIncldLst{jDtct}.lastFm) > 1 && ismember(iDtct, dtctIncldLst{jDtct}.lastFm)
%             toMergeDtcts = dtctIncldLst{jDtct}.lastFm;
%         elseif length(dtctIncldLst{jDtct}.nextFm) > 1 && ismember(iDtct, dtctIncldLst{jDtct}.nextFm)
%             toMergeDtcts = dtctIncldLst{jDtct}.nextFm;
%         else
%             continue;
%         end
%         
%         cellLst_res{t0}.VoxelIdxList{iDtctSub0} = ...
%             cell2mat(cellLst_res{t0}.VoxelIdxList(cellSubIdxsLUT(toMergeDtcts, 2)));
%         toMergeDtcts(toMergeDtcts == iDtct) = [];
%         for ll = 1:length(toMergeDtcts)
%             cellLst_res{t0}.VoxelIdxList{cellSubIdxsLUT(toMergeDtcts(ll), 2)} = [];
%         end
%         mergedFlagVec(ii) = true;
%         disp(['Merging... Detection #', num2str(iDtct), ' into ', num2str(toMergeDtcts)]);
%         break;
%     end
% end
% 
% toc



%%   Update the detection list with the new results
tic

cellLst_res = cellLst_cur;
if all(cellfun(@isempty, newCellVxLst))
    isChanged = false;
    return;
else
    isChanged = true;
    disp(['Modified ',num2str(nnz(~cellfun(@isempty, newCellVxLst))),' detections']);
end


%%%  Remove dtcts2Split from cellLst_res & Add newCellVxLst to cellLst_res
for ii = 1:length(dtcts2Split)
    if ~isempty(newCellVxLst{ii})
        iDtct = dtcts2Split(ii);
        t0 = cellSubIdxsLUT(iDtct, 1);
        dtctSub0 = cellSubIdxsLUT(iDtct, 2);
        cellLst_res{t0}.VoxelIdxList{dtctSub0} = [];
        
        cellLst_res{t0}.VoxelIdxList(end+1:end+length(newCellVxLst{ii})) = newCellVxLst{ii};
    end
end



%%%  Update all features in cellLst_res
for tt = 1:nnT
    cellLst_res{tt}.VoxelIdxList(cellfun(@isempty, cellLst_res{tt}.VoxelIdxList)) = [];
    cellLst_res{tt}.NumObjects = length(cellLst_res{tt}.VoxelIdxList);
    cellLst_res{tt}.CoordinateRgs = nan(cellLst_res{tt}.NumObjects, 6);
    cellLst_res{tt}.ctrPt = zeros(cellLst_res{tt}.NumObjects, 3);
    cellLst_res{tt}.dtctMtScVec = nan(cellLst_res{tt}.NumObjects, 2);
    for iDtctSub = 1:cellLst_res{tt}.NumObjects
        [yY, xX, zZ] = ind2sub([nnY,nnX,nnZ], cellLst_res{tt}.VoxelIdxList{iDtctSub});
        cellLst_res{tt}.CoordinateRgs(iDtctSub,:) = [min(yY),max(yY), min(xX),max(xX), min(zZ),max(zZ)];
        cellLst_res{tt}.ctrPt(iDtctSub,:) = mean([yY, xX, zZ], 1);
    end
    imgStkIntrplt = imgStack(:,:,:,tt);
    cellLst_res{tt}.avgItstyVec = cellfun(@(x) mean(imgStkIntrplt(x)), cellLst_res{tt}.VoxelIdxList);
    cellLst_res{tt}.areaVec = cellfun(@length, cellLst_res{tt}.VoxelIdxList);
end


for tt = 1:(nnT-1)
    img1 = imgStack(:,:,:,tt);
    img2 = imgStack(:,:,:,tt+1);
    
    cellLst_res{tt+1}.dtctMtScVec = nan(cellLst_res{tt+1}.NumObjects, 2);
    
    effVxBnds1 = [(effZinFmTwoEndIdx(tt,1)-1)*nnY*nnX + 1, effZinFmTwoEndIdx(tt,2)*nnY*nnX];
    % disp(var(img1([1:(effVxBnds1(1)-1), (effVxBnds1(2)+1):end])));
    effVxBnds2 = [(effZinFmTwoEndIdx(tt+1,1)-1)*nnY*nnX + 1, effZinFmTwoEndIdx(tt+1,2)*nnY*nnX];
    
    for iSub = 1:cellLst_res{tt}.NumObjects
        vxVec = cellLst_res{tt}.VoxelIdxList{iSub};
        temp = vxVec >= effVxBnds2(1) & vxVec <= effVxBnds2(2);
        if nnz(temp) >= 0.5*length(vxVec)
            vxVec = vxVec(temp);
            cellLst_res{tt}.dtctMtScVec(iSub, 2) = sqrt(mean((img1(vxVec)-img2(vxVec)).^2));
        end
    end
    
    for iSub = 1:cellLst_res{tt+1}.NumObjects
        vxVec = cellLst_res{tt+1}.VoxelIdxList{iSub};
        temp = vxVec >= effVxBnds1(1) & vxVec <= effVxBnds1(2);
        if nnz(temp) >= 0.5*length(vxVec)
            vxVec = vxVec(temp);
            cellLst_res{tt+1}.dtctMtScVec(iSub, 1) = sqrt(mean((img1(vxVec)-img2(vxVec)).^2));
        end
    end
end




toc



%%    Visualize the result
% nDtct_res = sum(cellfun(@(x) x.NumObjects, cellLst_res));
% cellIdxsLst_res = cellfun(@(x) (1:x.NumObjects)', cellLst_res, 'UniformOutput', false);
% for t = 2:nnT
%     cellIdxsLst_res{t}(:,1) = cellIdxsLst_res{t}(:,1) + cellIdxsLst_res{t-1}(end,1);
% end
% cellSubIdxsLUT_res = cell(nnT,1);
% for t = 1:nnT
%     cellSubIdxsLUT_res{t} = [t + zeros(cellLst_res{t}.NumObjects,1), (1:cellLst_res{t}.NumObjects)'];
% end
% cellSubIdxsLUT_res = cell2mat(cellSubIdxsLUT_res);
% 
% 
% %%%  All detections
% traceVxLstLst = cell(nDtct_res, 1);
% for iDtct = 1:nDtct_res
%     traceVxLstLst{iDtct} = cell(nnT,1);
%     tt = cellSubIdxsLUT_res(iDtct, 1);
%     iDtctSub = cellSubIdxsLUT_res(iDtct, 2);
%     traceVxLstLst{iDtct}{tt} = cellLst_res{tt}.VoxelIdxList{iDtctSub};
% end
% 
% %%%  New detections
% traceVxLstLst = cell(sum(cellfun(@length, newCellVxLst)), 1);
% tempCnt = 0;
% for ii = 1:length(dtcts2Split)
%     if ~isempty(newCellVxLst{ii})
%         iDtct = dtcts2Split(ii);
%         tt = cellSubIdxsLUT(iDtct, 1);
%         for jj = 1:length(newCellVxLst{ii})
%             tempCnt = tempCnt + 1;
%             traceVxLstLst{tempCnt} = cell(nnT,1);
%             traceVxLstLst{tempCnt}{tt} = newCellVxLst{ii}{jj};
%         end
%     end
% end
% 
% % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, [], false, []);
% % % [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, (cellIdPatchStack>0)*1, false, []);
% [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, false, []);
% % [printImgStack,printImgYZStack,~] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, invPCmapStack, true, cmapLb);
% % implay(printImgStack);
% % implay(printImgYZStack);
% 






