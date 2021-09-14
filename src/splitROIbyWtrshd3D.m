function subRoiLbMap = splitROIbyWtrshd3D(roiBnMap, seedLbMap, scMap)

% roiBnMap = roiBnPatch;
% seedLbMap = curLbMap;
% scMap = 1- dist2BdMap/max(dist2BdMap(:));
% figure; imshow(scMap(:,:,10))

seedBnMap = seedLbMap > 0;

tempMap = scMap;
tempMap(~roiBnMap) = 1;
mdfdScMap = imimposemin(tempMap, seedBnMap);
% figure; imshow(mdfdScMap(:,:,10));

resLbMap = watershed(mdfdScMap);
resLbMap(~roiBnMap) = 0;
% figure; imshow(double(resLbMap(:,:,10))/double(max(resLbMap(:))));

subRoiLbMap = double(resLbMap);
% myImageStackPrint(subRoiLbMap, 7, true);
% implay(subRoiLbMap)