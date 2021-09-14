function showImg = myImOverlay_multi(img, lbmapIdxLst, flagEdgeTrans, flagShow, markSize)

cmapLb = autumn(length(lbmapIdxLst));
[~,temp] = sort(rand(size(cmapLb,1),1));
cmapLb = cmapLb(temp,:);

% cmapLb = [1, 0, 0];
cmapImg = [0, 1, 0];

SE = strel('diamond',markSize);
lbmap = false(size(img));

[nR, nC] = size(img);

showImg = zeros(nR, nC, 3);
for cl = 1:3
    temp2 = img * cmapImg(cl);
    for iLb = 1:length(lbmapIdxLst)
        lbmap(:) = false;
        lbmap(lbmapIdxLst{iLb}) = true;
        if flagEdgeTrans
            lbmap = imdilate(lbmap,SE) & ~lbmap;
        else
            lbmap = imdilate(lbmap,SE);
        end
        temp2(lbmap) = cmapLb(iLb,cl);
    end
    showImg(:,:,cl) = temp2;
end

if flagShow
    figure; imshow(showImg);
end
