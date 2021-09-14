function [showImg, cmapLb] = myImOverlay3D_multi(imgBg, lbIdxLst, lbIdMap, flagShow, isGivenCmapLb, cmapLb)

if max(imgBg(:)) > 1
    imgBg = imgBg / max(imgBg(:));
end

if isempty(lbIdxLst)
    idVec = unique(lbIdMap(:));
    idVec(idVec == 0) = [];
    lbIdxLst = cell(length(idVec),1);
    for jId = 1:length(idVec)
        lbIdxLst{jId} = find(lbIdMap == idVec(jId));
    end
end

if ~isGivenCmapLb
    %     % cmapLb = hsv(length(lbIdxLst));
    %     % cmapLb = jet(length(lbIdxLst));
    %     cmapLb = rand(length(lbIdxLst), 3);
    %     cmapLb = cmapLb./(max(cmapLb,[],2) * ones(1,3));
    %     cmapLb = cmapLb*0.7 + 0.3;
    %     invalidClInds = sum(cmapLb>0.9,2)==3 | sum(cmapLb<0.1,2)==3;
    %     if any(invalidClInds)
    %         temp = rand(nnz(invalidClInds), 3);
    %         temp = temp./(max(temp,[],2) * ones(1,3));
    %         temp = temp*0.7 + 0.3;
    %         cmapLb(invalidClInds, :) = temp;
    %     end
    nCl = length(lbIdxLst);
    hueVec = (0.5:nCl-0.5)'/nCl + 0.1*randn(nCl,1);
    hueVec = max(0,min(1,hueVec));
    [~,temp] = sort(rand(nCl,1));
    hueVec = hueVec(temp);
    satVec = 0.25 + 0.5*rand(nCl,1);
    brValVec = 0.7 + 0.3*rand(nCl,1);
    cmapLb = hsv2rgb([hueVec, satVec, brValVec]);
end

cmapImg = [1, 1, 1];
% cmapImg = [1, 1, 1] * 0.8;

[nR, nC, nZ] = size(imgBg);

showImg = zeros(nR, nC, nZ, 3);
lbmap = false(size(imgBg));
for cl = 1:3
    temp2 = imgBg * cmapImg(cl);
    for jId = 1:length(lbIdxLst)
        lbmap(:) = false;
        lbmap(lbIdxLst{jId}) = true;
        % temp2(lbmap) = cmapLb(jId,cl);
        % temp2(lbmap) = imgBg(lbmap) * cmapLb(jId,cl);
        temp2(lbmap) = sqrt(imgBg(lbmap)) * cmapLb(jId,cl);
    end
    showImg(:,:,:,cl) = temp2;
end

showImg = permute(showImg,[1, 2, 4, 3]);

if flagShow
    implay(showImg);
end