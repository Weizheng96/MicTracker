for n=1:5

    Im1 = uint8(tifread("/home/wei/Projects/MicTrackerMatab/testData/"+n+".tif"));
    outLabel=Im1(151:250,201:400,40:139);

    ImName=num2str(n);

    imwrite(outLabel(:,:,1),[ImName,'.tif']);
    for i = 2:size(outLabel,3)
        imwrite(outLabel(:,:,i),[ImName,'.tif'],'WriteMode','append');
    end

end

