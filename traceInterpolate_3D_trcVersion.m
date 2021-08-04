function intpltdTrace = traceInterpolate_3D_trcVersion(trace0, imgSzNrg)


nnY = imgSzNrg.size(1);
nnX = imgSzNrg.size(2);
nnZ = imgSzNrg.size(3);
% nnT = imgSzNrg.size(4);



%%%%  Voxel lists of the selected trace (infer the undetected part)
%     BW = false(nnY,nnX,nnZ);
for iiT = 1:(size(trace0.dtctSubs,1)-1)
    % iiT = 3;
    t1 = trace0.dtctSubs(iiT,1);
    t2 = trace0.dtctSubs(iiT+1,1);
    % t1 = 3;
    % t2 = 5;
    
    if t2 == t1 + 1
        continue;
    end
    
    %%%  Infer the frames between t1 and t2
    tInfr = t1+1:t2-1;
    % tInfr = 4;
    % disp(['Inferring frame ', num2str(tInfr)]);
    intrpltRatios = (tInfr - t1) / (t2 - t1);
    intrpltVxLstLst = bwInterpolate3D([nnY,nnX,nnZ], trace0.traceVxLst{t1}, trace0.traceVxLst{t2}, intrpltRatios);
    trace0.traceVxLst(tInfr) = intrpltVxLstLst;
    
    
    %         %%%%%%    Visualize    %%%%%%
    %         % vxLstShow = trace0.traceVxLst{t2};
    %         vxLstShow = intrpltVxLstLst{2};
    %         BW = false([nnY,nnX,nnZ]);
    %         BW(vxLstShow) = true;
    %         vxSub1 = zeros(length(trace0.traceVxLst{t1}), 3);
    %         [vxSub1(:,1),vxSub1(:,2),vxSub1(:,3)] = ind2sub(imgSize, vxLst1);
    %         vxSub2 = zeros(length(trace0.traceVxLst{t2}), 3);
    %         [vxSub2(:,1),vxSub2(:,2),vxSub2(:,3)] = ind2sub(imgSize, vxLst2);
    %         ptchRgs = [max(min([vxSub1;vxSub2])' - 10, 1), min([max([vxSub1;vxSub2])' + 10, imgSize'],[],2)];
    %         myImageStackPrint(BW(ptchRgs(1,1):ptchRgs(1,2), ptchRgs(2,1):ptchRgs(2,2),...
    %             ptchRgs(3,1):ptchRgs(3,2)), 10, true);
    %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%   Pre-calculate features of traces: spatial position
    for iInfrT = 1:length(tInfr)
        tt = tInfr(iInfrT);
        [yY,xX,zZ] = ind2sub([nnY,nnX,nnZ], trace0.traceVxLst{tt});
        trace0.posBoundsPerT(tt,:) = [min(yY), max(yY), min(xX), max(xX), min(zZ), max(zZ)];
        trace0.cntrPtSub(tt,:) = [mean(yY), mean(xX), mean(zZ)];
    end
end

temp = [min(trace0.posBoundsPerT); max(trace0.posBoundsPerT)];
trace0.posBounds = temp(logical([1,0,1,0,1,0; 0,1,0,1,0,1]))';


intpltdTrace = trace0;



