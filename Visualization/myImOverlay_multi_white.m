function showImg = myImOverlay_multi_white(img, lbmapIdxLst, flagEdgeTrans, ...
    flagShow, markSize, isGivenCmapLb, cmapLb)

if ~isGivenCmapLb
    % cmapLb = hsv(length(lbmapIdxLst));
    % cmapLb = jet(length(lbmapIdxLst));
    cmapLb = rand(length(lbmapIdxLst), 3);
    cmapLb = cmapLb./(max(cmapLb,[],2) * ones(1,3));
    cmapLb = cmapLb*0.8 + 0.2;
end

% [~,temp] = sort(rand(length(lbmapIdxLst),1));
% cmapLb = cmapLb(temp,:);

% cmapLb = [1, 0, 0];
cmapImg = [1, 1, 1];

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
