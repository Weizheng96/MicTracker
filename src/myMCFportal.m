function resFlowVec = myMCFportal(arcVec, nodeVec, nnNode, nnArc)
% % Input:
% %     arcVec   -  nnArc*5 array, (format: [<tail> <head> <capacity l.b.> <capacity u.b> <cost>])
% %     nodeVec  -  nSourceSink*2 array, (format: [<nodeId>, <supply>])

%%%%  Write the graph into the log file
fileID = fopen('myCs2Imput.txt', 'w');
fprintf(fileID,'p min %ld %ld\n', nnNode, nnArc);
if ~isempty(nodeVec)
    tempMat = nodeVec;
    tempMat(:,1) = tempMat(:,1) - 1;  % % Different index for C
    fprintf(fileID,'n %ld %ld\n', tempMat');
end
if ~isempty(arcVec)
    tempMat = arcVec;
    tempMat(:,1:2) = tempMat(:,1:2) - 1;  % % Different index for C
    fprintf(fileID,'a %ld %ld %ld %ld %ld\n', tempMat');
end
fclose(fileID);



%%%%  Run the min cost flow/circulation solver
mex cs2_mex.c
cs2_mex();



%%%%  Read the results and interprete them into conclusions
fileID = fopen('myCs2Result.txt', 'r');
resFlowVec = fscanf(fileID,'%ld %ld %ld', [3, Inf]);
fclose(fileID);

resFlowVec = resFlowVec';
resFlowVec(:, 1:2) = resFlowVec(:, 1:2) + 1;   % % Back to MATLAB index





