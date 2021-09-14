function [res_G, trajectories] = sspTracker(trackG, para)
% successive shortest path for min-cost flow
vNum = para.nCellRep;
res_G = trackG;
pSet = cell(vNum,1);
numTracks = 0;
tic;
for i=1:vNum
    %     if i==940
    %         keyboard;
    %     end
    [P,d] = shortestpath(res_G, para.source, para.sink);
    if d>=-0
        break;
    end
    %     if sum(P==14082*2+1)>0
    %         %keyboard;
    %         pp = P(mod(P,2)==0);
    %         pp = pp/2;
    %     end
    fprintf('The %d paths found with cost: %.2f\n', i, d);
    % change weight as -1*weight
    %     for j=1:length(P)-1
    %         tmpIdx = find(res_G.Edges.EndNodes(:,1) == P(j) & res_G.Edges.EndNodes(:,2) == P(j+1)); % should only have one
    %         if length(tmpIdx) ~= 1
    %             error('There should be only one edge be changed.');
    %         end
    %         res_G.Edges.Weight(tmpIdx) = -1*res_G.Edges.Weight(tmpIdx);
    %     end
    
    % inverse edge
    edgeIdx = findedge(res_G, P(1:end-1), P(2:end));
    res_G = addedge(res_G, P(2:end), P(1:end-1), -1*res_G.Edges.Weight(edgeIdx));
    res_G = rmedge(res_G, P(1:end-1), P(2:end));

    numTracks = numTracks + 1;
    pSet{i} = P;
end
timeLapse = toc;

pSet(numTracks+1:end) = [];

fprintf('finish successive shortest path alg with %d paths, %.3fs!\n',numTracks, timeLapse);


% Recover the tracks
trajectories = cell(numTracks,1);    
stNodes = res_G.Edges.EndNodes(res_G.Edges.EndNodes(:,1)== para.sink,2);
% pathCost = zeros(numTracks,1);

for i = 1:numTracks
    disp(['Recovering tracks: ',num2str(i),'/',num2str(numTracks)]);
    head = stNodes(i);
%     pathCost(i) = pathCost(i)+res_G.Edges.Weight(findedge(res_G, para.t, head));
%     pathCost(i) = pathCost(i)+res_G.Edges.Weight(findedge(res_G, head, head-1));
    
    trajectories{i} = stNodes(i);
    lastHead = head;
    head = res_G.Edges.EndNodes(res_G.Edges.EndNodes(:,1) == head-1, 2);
    
    while head ~= para.source
        if length(head) ~= 1
            error('There should be only one edge.');
        end
%         pathCost(i) = pathCost(i)+res_G.Edges.Weight(findedge(res_G, lastHead-1, head));
%         pathCost(i) = pathCost(i)+res_G.Edges.Weight(findedge(res_G, head, head-1));
        trajectories{i} = [trajectories{i}, head];
        
        lastHead = head;
        head = res_G.Edges.EndNodes(res_G.Edges.EndNodes(:,1) == head-1, 2);
    end
%     pathCost(i) = pathCost(i)+res_G.Edges.Weight(findedge(res_G, lastHead-1, para.s));
    trajectories{i}  = trajectories{i}(end:-1:1)/2;
end

% pathCost = -pathCost;


% %%%%%%     For debugging     %%%%%%
% pathOI = [2403 923 924 945 946 2404];
% pathOI = [2404 826 825 806 805 2403];
% % pathOI = [2403 923 924 945 946 2404 826 825 806 805 2403];
% testCost = 0;
% for ii = 1:length(pathOI)-1
%     testCost = testCost + res_G.Edges.Weight(findedge(res_G, pathOI(ii), pathOI(ii+1)));
% end
% disp(testCost);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


