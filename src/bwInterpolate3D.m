function intrpltVxLstLst = bwInterpolate3D(imgSize, vxLst1, vxLst2, intrpltRatios)

% imgSize = [nY,nX,nZ];
% vxLst1 = trace0.traceVxLst{t1};
% vxLst2 = trace0.traceVxLst{t2};
% intrpltRatios = [0.25,0.5,0.75];


%%%%  Output initialization
intrpltVxLstLst = cell(length(intrpltRatios),1);


%%%%  Centering
vxSub1 = zeros(length(vxLst1), 3);
[vxSub1(:,1),vxSub1(:,2),vxSub1(:,3)] = ind2sub(imgSize, vxLst1);
cntr1 = mean(vxSub1);

vxSub2 = zeros(length(vxLst2), 3);
[vxSub2(:,1),vxSub2(:,2),vxSub2(:,3)] = ind2sub(imgSize, vxLst2);
cntr2 = mean(vxSub2);

midCntr = (cntr1 + cntr2) / 2;
vxSub1cntrd = round(vxSub1 + ones(length(vxLst1), 1) * (midCntr - cntr1));
vxSub2cntrd = round(vxSub2 + ones(length(vxLst2), 1) * (midCntr - cntr2));
temp = ones(length(vxLst1), 1) * imgSize;
temp2 = min(max(vxSub1cntrd(:), 1), temp(:));
vxSub1cntrd = reshape(temp2, length(vxLst1), 3);
temp = ones(length(vxLst2), 1) * imgSize;
temp2 = min(max(vxSub2cntrd(:), 1), temp(:));
vxSub2cntrd = reshape(temp2, length(vxLst2), 3);


intrpltCntrs = (1-intrpltRatios)'*cntr1 + intrpltRatios'*cntr2;


%%%%  Boundaries of the centered detections
SE = true(3,3,3);

ptchRgs = [max(min([vxSub1;vxSub2])' - 10, 1), min([max([vxSub1;vxSub2])' + 10, imgSize'],[],2)];

BW1 = false(imgSize);
BW1(sub2ind(imgSize,vxSub1cntrd(:,1),vxSub1cntrd(:,2),vxSub1cntrd(:,3))) = true;
% myImageStackPrint(BW1(ptchRgs(1,1):ptchRgs(1,2), ptchRgs(2,1):ptchRgs(2,2),...
%     ptchRgs(3,1):ptchRgs(3,2)), 10, true);

BW2 = false(imgSize);
BW2(sub2ind(imgSize,vxSub2cntrd(:,1),vxSub2cntrd(:,2),vxSub2cntrd(:,3))) = true;
% myImageStackPrint(BW2(ptchRgs(1,1):ptchRgs(1,2), ptchRgs(2,1):ptchRgs(2,2),...
%     ptchRgs(3,1):ptchRgs(3,2)), 10, true);


%%%%  Get boundaries
outerBndofUnionBW = BW1 | BW2;
outerBndofUnionBW = outerBndofUnionBW & ~imerode(outerBndofUnionBW,SE);
% myImageStackPrint(outerBndofUnionBW(ptchRgs(1,1):ptchRgs(1,2), ptchRgs(2,1):ptchRgs(2,2),...
%     ptchRgs(3,1):ptchRgs(3,2)), 10, true);

tempMap = false(imgSize+2);
tempMap(2:end-1, 2:end-1, 2:end-1) = BW1;
tempMap = tempMap & ~imerode(tempMap,SE);
BW1 = tempMap(2:end-1, 2:end-1, 2:end-1);

tempMap = false(imgSize+2);
tempMap(2:end-1, 2:end-1, 2:end-1) = BW2;
tempMap = tempMap & ~imerode(tempMap,SE);
BW2 = tempMap(2:end-1, 2:end-1, 2:end-1);

% BW1 = BW1 & ~imerode(BW1,SE);
% BW2 = BW2 & ~imerode(BW2,SE);

bdVxLst1cntrd = find(BW1);
[yy,xx,zz] = ind2sub(imgSize, bdVxLst1cntrd);
% bdVxSub1cntrd = [yy,xx,zz];
[azimuth1,elevation1,r1] = cart2sph(xx-midCntr(2),yy-midCntr(1),zz-midCntr(3));

bdVxLst2cntrd = find(BW2);
[yy,xx,zz] = ind2sub(imgSize, bdVxLst2cntrd);
% bdVxSub2cntrd = [yy,xx,zz];
[azimuth2,elevation2,r2] = cart2sph(xx-midCntr(2),yy-midCntr(1),zz-midCntr(3));

bdVxLstOutBndcntrd = find(outerBndofUnionBW);
[yy,xx,zz] = ind2sub(imgSize, bdVxLstOutBndcntrd);
% bdVxSubOutBndcntrd = [yy,xx,zz];
[azimuthOuter,elevationOuter,rOuter] = cart2sph(xx-midCntr(2),yy-midCntr(1),zz-midCntr(3));



%%%%  For each voxel on the outer boundary (of the union), find the interpolated points
azimuthMid = zeros(length(azimuthOuter),1);
elevationMid = zeros(length(azimuthOuter),1);
% rMid = zeros(length(azimuthOuter),1);
rMids = zeros(length(azimuthOuter),length(intrpltRatios));

for iVx = 1:length(azimuthOuter)
    % iVx = 228;
    [~, I1] = min((azimuth1 - azimuthOuter(iVx)).^2 + (elevation1 - elevationOuter(iVx)).^2);
    [~, I2] = min((azimuth2 - azimuthOuter(iVx)).^2 + (elevation2 - elevationOuter(iVx)).^2);
    %     disp([azimuth1(I1),elevation1(I1),r1(I1)]);
    %     disp([azimuth2(I2),elevation2(I2),r2(I2)]);
    %     disp([azimuthOuter(iVx),elevationOuter(iVx),rOuter(iVx)]);
    azimuthMid(iVx) = (azimuth1(I1) + azimuth2(I2)) / 2;
    elevationMid(iVx) = (elevation1(I1) + elevation2(I2)) / 2;
    % rMid(iVx) = (r1(I1) + r2(I2)) / 2;
    rMids(iVx,:) = r1(I1)*(1-intrpltRatios) + r2(I2)*intrpltRatios;
    %     disp([azimuthMid(iVx),elevationMid(iVx),rMid(iVx)]);
end
% [xx,yy,zz] = sph2cart(azimuthMid, elevationMid, rMid);
% bdVxSubMid = unique(round([yy,xx,zz] + ones(length(xx),1) * midCntr), 'rows');
% bdVxLstMid = sub2ind(imgSize, bdVxSubMid(:,1), bdVxSubMid(:,2), bdVxSubMid(:,3));


%%%%  Get the map for each (d-d1)/(d2-d1) ratio
for iRatio = 1:length(intrpltRatios)
    %%%  Fill all line segments between boundary points and the center
    azimuthFull = zeros(ceil(max(rMids(:, iRatio)))*length(azimuthMid), 1);
    elevationFull = zeros(length(azimuthFull),1);
    rFull = zeros(length(azimuthFull),1);
    count = 0;
    for iVx = 1:length(azimuthMid)
        idxs = (count+1) : (count+ceil(rMids(iVx, iRatio)));
        azimuthFull(idxs) = azimuthMid(iVx);
        elevationFull(idxs) = elevationMid(iVx);
        rFull(idxs(1:end-1)) = 1 : (ceil(rMids(iVx, iRatio))-1);
        rFull(idxs(end)) = rMids(iVx, iRatio);
        count = count + length(idxs);
    end
    azimuthFull((count+1):end) = [];
    elevationFull((count+1):end) = [];
    rFull((count+1):end) = [];
    [xx,yy,zz] = sph2cart(azimuthFull, elevationFull, rFull);
    vxSubMid = unique(round([yy,xx,zz] + ones(length(xx),1) * midCntr), 'rows');
    vxLstMid = sub2ind(imgSize, vxSubMid(:,1), vxSubMid(:,2), vxSubMid(:,3));
    
    
    %%%  Some modifications to the map
    BW = false(imgSize);
    BW(vxLstMid) = true;
    % seR = 1;
    % SE = false(seR*2+1, seR*2+1, seR*2+1);
    % SE(seR+1,seR+1,seR+1) = true;
    % temp = bwdist(SE);
    % SE = temp <= seR;
    BW = imclose(BW, true(5,5));
    % BW = imfill(BW,'holes');
    % BW = BW & ~imerode(BW,SE);
    % myImageStackPrint(BW(ptchRgs(1,1):ptchRgs(1,2), ptchRgs(2,1):ptchRgs(2,2),...
    %     ptchRgs(3,1):ptchRgs(3,2)), 10, true);
    
    %%%  Re-centered to interpolated centers
    vxLstMid = find(BW);
    [yy,xx,zz] = ind2sub(imgSize, vxLstMid);
    vxSubMid = round([yy,xx,zz] + ones(length(yy),1)*(intrpltCntrs(iRatio,:) - midCntr));
    vxLstMid = sub2ind(imgSize, vxSubMid(:,1), vxSubMid(:,2), vxSubMid(:,3));

    intrpltVxLstLst{iRatio} = vxLstMid;
end









