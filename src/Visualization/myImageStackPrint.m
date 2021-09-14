function printImg = myImageStackPrint(showImgStackPatch, nSubsInRow, isShow)

if length(size(showImgStackPatch)) == 4
    [subR,subC,~,subT] = size(showImgStackPatch);
    imDim = 4;
elseif length(size(showImgStackPatch)) == 3
    [subR,subC,subT] = size(showImgStackPatch);
    imDim = 3;
else
    printImg = [];
    return;
end

blankWidth = 3;
% nSubsInRow = 8;

printImg = ones(subR*ceil(subT/nSubsInRow) + blankWidth*ceil(subT/nSubsInRow-1),...
    subC*nSubsInRow + blankWidth*(nSubsInRow-1), 3);

for t = 1:subT
    rSubIdx = ceil(t/nSubsInRow);
    cSubIdx = t - (rSubIdx-1)*nSubsInRow;
    blankShiftR = blankWidth*(rSubIdx-1);
    blankShiftC = blankWidth*(cSubIdx-1);

    if imDim == 4
        tempPatch = showImgStackPatch(:,:,:,t);
    else
        tempPatch = zeros(subR,subC,3);
        tempPatch(:,:,1) = showImgStackPatch(:,:,t);
        tempPatch(:,:,2) = showImgStackPatch(:,:,t);
        tempPatch(:,:,3) = showImgStackPatch(:,:,t);
    end
    printImg(subR*(rSubIdx-1)+blankShiftR+1 : subR*rSubIdx + blankShiftR,...
        subC*(cSubIdx-1)+blankShiftC+1:subC*cSubIdx+blankShiftC, :)...
        = tempPatch;
end

if isShow
    figure; imshow(printImg);
end

