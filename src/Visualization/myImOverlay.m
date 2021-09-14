function showImg = myImOverlay(img, lbmap, flagEdgeTrans, flagShow, markSize)

% cmap = jet(size(img,3));
cmapLb = [1, 0, 0];
cmapImg = [0, 1, 0];

SE = strel('diamond',markSize);
if flagEdgeTrans
    lbmap = imdilate(lbmap,SE) & ~lbmap;
else
    lbmap = imdilate(lbmap,SE);
end

[nR, nC] = size(img);

showImg = zeros(nR, nC, 3);
for cl = 1:3
    temp2 = img * cmapImg(cl);
    temp2(lbmap) = cmapLb(cl);
    showImg(:,:,cl) = temp2;
end

if flagShow
    figure; imshow(showImg);
end
