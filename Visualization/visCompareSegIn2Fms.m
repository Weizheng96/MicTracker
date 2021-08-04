function [showImg, showImgYZ] = visCompareSegIn2Fms(lbIdMap1, lbIdMap2)

if any(size(lbIdMap1) ~= size(lbIdMap2))
    disp('Visualization of between-frame segmentation comparison: Image sizes are incompatible.')
    return;
end

vals1 = unique(lbIdMap1(:));
vals1(vals1==0) = [];
for iVal = 1:length(vals1)
    lbIdMap1(lbIdMap1 == vals1(iVal)) = iVal;
end

vals2 = unique(lbIdMap2(:));
vals2(vals2==0) = [];
for iVal = 1:length(vals2)
    lbIdMap2(lbIdMap2 == vals2(iVal)) = iVal;
end

nRoi1 = max(lbIdMap1(:));
nRoi2 = max(lbIdMap2(:));

SE = reshape([0,0,0; 0,1,0; 0,0,0; 0,1,0; 1,1,1; 0,1,0; 0,0,0; 0,1,0; 0,0,0],3,3,3);
bndryMap1 = lbIdMap1 - imerode(lbIdMap1, SE);
bndryMap2 = lbIdMap2 - imerode(lbIdMap2, SE);
% myImageStackPrint(bndryMap2, 5, true);

showMatchMap = bndryMap1;
showMatchMap(bndryMap2>0) = bndryMap2(bndryMap2>0) + nRoi1;

temp1 = hsv(nRoi1*3+1);
temp2 = jet(nRoi2*3+1);
%     cmapLb = [temp1(1:nRoi1,:); temp2(nRoi2*2+1:end,:)];
cmapLb = [temp1(nRoi1+2 : nRoi1*2+1, :); temp2(nRoi2*2+2 : end, :)];

[showImg, ~] = myImOverlay3D_multi((showMatchMap>0)*1, [], showMatchMap, false, true, cmapLb);
% myImageStackPrint(showImg, 10, true);
showImgYZ = permute(showImg,[1, 4, 3, 2]);
% myImageStackPrint(showImgYZ, 10, true);
% % myImageStackPrint(roiIdPatch1, 5, true);
% % myImageStackPrint(roiIdPatch2, 5, true);

