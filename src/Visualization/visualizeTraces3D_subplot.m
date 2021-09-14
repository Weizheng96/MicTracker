function [printImgStack,printImgYZStack,cmapLb] = visualizeTraces3D_subplot(traceVxLstLst, imgSzNrg, bgImgStack, isCmapGiven, cmapLb)
% %  Visualize any given trace(s)
% %  Input:
% %     traceVxLstLst  --  a list of traces. Each element of traceVxLstLst is a 
% %                        nT-long list of voxel lists, representing a trace.
% %     imgSzNrg       --  imgSzNrg.size: 4*1 vector, size of 3D+T image stack
% %                        imgSzNrg.range: 4*2 array, selected range of visualization
% %     bgImgStack     --  background image for visualization overlay. If [], set the background as pure black

nTrc = length(traceVxLstLst);

nYptch = imgSzNrg.range(1,2) - imgSzNrg.range(1,1) + 1;
nXptch = imgSzNrg.range(2,2) - imgSzNrg.range(2,1) + 1;
nZptch = imgSzNrg.range(3,2) - imgSzNrg.range(3,1) + 1;
nTptch = imgSzNrg.range(4,2) - imgSzNrg.range(4,1) + 1;

nT = imgSzNrg.size(4);
nSubsInRow = 5;


%%%%  Color map of traces
if ~isCmapGiven
    nCl = length(traceVxLstLst);
    hueVec = (0.5:nCl-0.5)'/nCl + 0.1*randn(nCl,1);
    hueVec = max(0,min(1,hueVec));
    [~,temp] = sort(rand(nCl,1));
    hueVec = hueVec(temp);
    satVec = 0.25 + 0.5*rand(nCl,1);
    brValVec = 0.7 + 0.3*rand(nCl,1);
    cmapLb = hsv2rgb([hueVec, satVec, brValVec]);
end

%%%%  LUT of frame-trace
assMat = false(nTrc,nT);
for iTrc = 1:nTrc
    % iTrc = 42;
    assMat(iTrc, :) = ~cellfun(@isempty, traceVxLstLst{iTrc})';
end
% figure; imshow(assMat);


%%%%  Generate images for display
nullImgStack = zeros(nYptch, nXptch, nZptch, nTptch);
printImg = myImageStackPrint(nullImgStack(:,:,1:3:nZptch,1), nSubsInRow, false);
printYZImg = myImageStackPrint(permute(nullImgStack(:,1:3:nXptch,:,1),[1, 3, 2]), nSubsInRow, false);
% myImageStackPrint(bgImgStack(:,:,1:3:nZptch,1), 10, true);

printImgStack = zeros([size(printImg,1), size(printImg,2), 3, nTptch]);
printImgYZStack = zeros([size(printYZImg,1), size(printYZImg,2), 3, nTptch]);
% showImgStack = zeros(nYptch, nXptch, nZptch, 3, nTptch);
for t = imgSzNrg.range(4,1) : imgSzNrg.range(4,2)
    % t = 6;
    disp(['Printing patch... t == ', num2str(t)]);
    tPtch = t - imgSzNrg.range(4,1) + 1;
    slctCmapLb = cmapLb(assMat(:,t),:);
    lbIdxLst = cellfun(@(x) x{t}, traceVxLstLst(assMat(:,t)), 'UniformOutput', false);
    for iDtct = 1:length(lbIdxLst)
        [yY,xX,zZ] = ind2sub(imgSzNrg.size(1:3), lbIdxLst{iDtct});
        yY = yY - imgSzNrg.range(1,1) + 1;
        xX = xX - imgSzNrg.range(2,1) + 1;
        zZ = zZ - imgSzNrg.range(3,1) + 1;
        tempLb = yY<1 | yY>nYptch | xX<1 | xX>nXptch | zZ<1 | zZ>nZptch;
        lbIdxLst{iDtct} = sub2ind([nYptch, nXptch, nZptch], yY(~tempLb), xX(~tempLb), zZ(~tempLb));
    end
    if isempty(bgImgStack)
        imgBg = nullImgStack(:,:,:,tPtch);
        imgBg(cell2mat(lbIdxLst)) = 1;
    else
        imgBg = bgImgStack(:,:,:,tPtch);
    end
    [showImg, ~] = myImOverlay3D_multi(imgBg, lbIdxLst, [], false, true, slctCmapLb);
    % implay(showImg);
    % myImageStackPrint(showImg(:,:,:,1:3:nZptch), 10, true);

    % showImgStack(:,:,:,:,tPtch) = showImg;
    printImgStack(:,:,:,tPtch) = myImageStackPrint(showImg(:,:,:,1:3:nZptch), nSubsInRow, false);
    % figure; imshow(printImgStack(:,:,:,tPtch));
    
    printImgYZStack(:,:,:,tPtch) = myImageStackPrint(permute(showImg(:,1:3:nXptch,:,:),[1,4,3,2]), nSubsInRow, false);
    % figure; imshow(printImgYZStack(:,:,:,tPtch));
end



