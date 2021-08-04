function eig_all = principalCurvature3d(imgIn, sigma, fmap)
% % Calculate 3D principle curvature on given foreground

%%%  get the hessian matrix
% [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz, sm_vid] = Hessian3D(imgIn,sigma);
[Dxx, Dyy, Dzz, Dxy, Dxz, Dyz, ~] = Hessian3D(imgIn,sigma);
if nargin < 3
    fmap = ones(size(imgIn)); % if we cal all the voxes
end

% if(sigma>1)
%     % Correct for scaling
%     c=(sigma^2);
%     Dxx = c*Dxx; Dxy = c*Dxy;
%     Dxz = c*Dxz; Dyy = c*Dyy;
%     Dyz = c*Dyz; Dzz = c*Dzz;
% end
% [h,w,z] = size(imgIn);

%%% test each connected component
eig_all = zeros(size(imgIn));

dir_y = zeros(size(imgIn));
dir_x = zeros(size(imgIn));
dir_z = zeros(size(imgIn));
% eig2 = zeros(size(synId));
% eig3 = zeros(size(synId));


s = regionprops(fmap, {'PixelIdxList'});
for i=1:numel(s)
    %disp(i);
    vox = s(i).PixelIdxList;
    xx = Dxx(vox); yy = Dyy(vox); zz = Dzz(vox);
    xy = Dxy(vox); xz = Dxz(vox); yz = Dyz(vox);
    
    C = zeros(numel(s(i).PixelIdxList),3);
    dir_xyz = zeros(numel(s(i).PixelIdxList),3);
    for j=1:numel(s(i).PixelIdxList)
    % parfor j=1:numel(s(i).PixelIdxList)
        MM = [xx(j), xy(j), xz(j);...
            xy(j), yy(j), yz(j);...
            xz(j), yz(j), zz(j)];
        [Evec,Eval] = eig(MM);
        dEval = diag(Eval);
        %[~,od] = sort(abs(dEval),'descend');
        %C(j,:) = dEval(od)';
        [c,od] = sort(dEval,'descend');
        C(j,:) = c';
        dir_xyz(j,:) = Evec(:, od(1))';
    end
    dir_x(vox) = dir_xyz(:,1);
    dir_y(vox) = dir_xyz(:,2);
    dir_z(vox) = dir_xyz(:,3);
    eig_all(vox) = C(:,1);
end

end