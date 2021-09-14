ROIs=cellLst{1}.VoxelIdxList;
labelMap=zeros(size(imgFm));
for i=1:length(ROIs)
    labelMap(ROIs{i})=i;
end


orgIm3d=imgStack(:,:,:,1);
addpath("/home/wei/Projects/Kenichi/ImageAnalysis/ImageAnalysis_ForKenichi/method/boyu");
addpath("/home/wei/Projects/Kenichi/ImageAnalysis/ImageAnalysis_ForKenichi/NewMethod/cc_ImHandle");
outIm = imdisplayWithROI3DMoreRandom(orgIm3d, labelMap);

implay(outIm);
tifwrite(outIm,'1');