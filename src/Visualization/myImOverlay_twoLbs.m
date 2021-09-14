function showImg = myImOverlay_twoLbs(img, lbmap1, lbmap2, flagEdgeTrans, flagShow, markSize)

% cmap = jet(size(img,3));
cmapImg = [1, 1, 1];
cmapLb1 = [0, 1, 1];
cmapLb2 = [1, 0, 0];


SE = strel('diamond',markSize);
if flagEdgeTrans
    lbmap1 = imdilate(lbmap1,SE) & ~lbmap1;
    lbmap2 = imdilate(lbmap2,SE) & ~lbmap2;
else
    lbmap1 = imdilate(lbmap1,SE);
    lbmap2 = imdilate(lbmap2,SE);
end

[nR, nC] = size(img);

showImg = zeros(nR, nC, 3);
for cl = 1:3
    temp2 = img * cmapImg(cl);
    temp2(lbmap1) = cmapLb1(cl);
    temp2(lbmap2) = cmapLb2(cl);
    showImg(:,:,cl) = temp2;
end

if flagShow
    figure; imshow(showImg);
end
