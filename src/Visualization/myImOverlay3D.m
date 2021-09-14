function showImg = myImOverlay3D(imgBg, lbmap, flagShow)

if max(imgBg(:)) > 1
    imgBg = imgBg / max(imgBg(:));
end

lbmap = lbmap > 0;

% cmap = jet(size(img,3));
cmapLb = [1, 0, 0];
cmapImg = [1, 1, 1];

[nR, nC, nZ] = size(imgBg);

showImg = zeros(nR, nC, nZ, 3);
for cl = 1:3
    temp2 = imgBg * cmapImg(cl);
    temp2(lbmap) = imgBg(lbmap) * cmapLb(cl);
    showImg(:,:,:,cl) = temp2;
end

showImg = permute(showImg,[1, 2, 4, 3]);

if flagShow
    implay(showImg);
end