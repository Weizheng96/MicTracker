x1 = -5:5;
x2 = -5:5;
x3 = -5:5;
[X1,X2,X3] = meshgrid(x1,x2,x3);
X = [X1(:) X2(:) X3(:)];

outLabel=zeros(11,11,11);

y = mvnpdf(X,[-3 0 0],eye(3)*2);
y = reshape(y,length(x2),length(x1),length(x3));
outLabel=outLabel+y;

y = mvnpdf(X,[3 0 0],eye(3)*2);
y = reshape(y,length(x2),length(x1),length(x3));
outLabel=outLabel+y;

outLabel=outLabel/max(outLabel(:))*255;

implay(mat2gray(outLabel));

outLabel=uint8(outLabel);

ImName='1';

imwrite(outLabel(:,:,1),[ImName,'.tif']);
for i = 2:size(outLabel,3)
    imwrite(outLabel(:,:,i),[ImName,'.tif'],'WriteMode','append');
end

Im1 = imread("/home/wei/Downloads/sample_data_and_analysis_results/1.tif",1);
Im2 = imread("/home/wei/Projects/MicTrackerMatab/testData/1.tif",1);
Im3 = imread("/home/wei/Projects/MicTrackerMatab/testData2/1.tif",1);