function [showImg, cmapLb] = myImHeatMap(imShowIn)

% pcMapPatch = pcXY3dMapStk(yRg(1):yRg(2), xRg(1):xRg(2), zRg(1):zRg(2));
% pcMapPatch = 1 - sqrt((pcMapPatch-min(pcMapPatch(:)))/(max(pcMapPatch(:))-min(pcMapPatch(:))));
% imShowIn = pcMapPatch;

[nR, nC, nZ] = size(imShowIn);

% cmapLb = hsv(1000);
cmapLb = jet(1000);
% cmapLb(1,:) = [0,0,0];
lbmap = round(pcMapPatch*999) + 1;

showImg = zeros(nR, nC, nZ, 3);
tempMap = zeros(nR, nC, nZ);
for cl = 1:3
    tempMap(:) = cmapLb(lbmap(:),cl);
    showImg(:,:,:,cl) = tempMap;
end
showImg = permute(showImg,[1, 2, 4, 3]);

% myImageStackPrint(showImg, 7, true);





