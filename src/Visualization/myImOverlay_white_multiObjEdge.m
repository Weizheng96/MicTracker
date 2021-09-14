function showImg = myImOverlay_white_multiObjEdge(img, lbmap, flagShow, markSize)

% cmap = jet(size(img,3));
cmapLb = [1, 0, 0];
cmapImg = [1, 1, 1];

SE = strel('diamond',markSize);

edgeMap = false(size(lbmap));
eleVec = unique(lbmap);
eleVec(eleVec==0) = [];
for ii = 1:length(eleVec)
    tempMap = lbmap==eleVec(ii);
    edgeMap = edgeMap | (imdilate(tempMap,SE) & ~tempMap);
end
% figure; imshow(edgeMap);

[nR, nC] = size(img);
showImg = zeros(nR, nC, 3);
for cl = 1:3
    temp2 = img * cmapImg(cl);
    temp2(edgeMap) = cmapLb(cl);
    showImg(:,:,cl) = temp2;
end

if flagShow
    figure; imshow(showImg);
end
