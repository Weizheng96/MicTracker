function [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz, F] = Hessian3D(Volume,Sigma)
%  This function Hessian3D filters the image with an Gaussian kernel
%  followed by calculation of 2nd order gradients, which aprroximates the
%  2nd order derivatives of the image.
% 
% [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(I,Sigma)
% 
% inputs,
%   I : The image volume, class preferable double or single
%   Sigma : The sigma of the gaussian kernel used. If sigma is zero
%           no gaussian filtering.
%
% outputs,
%   Dxx, Dyy, Dzz, Dxy, Dxz, Dyz: The 2nd derivatives
%
% Function is written by D.Kroon University of Twente (June 2009)
% defaults
if nargin < 2, Sigma = 1; end

if length(Sigma)==1
    if(Sigma>0)
        F=imgaussian(Volume,Sigma);
    else
        F=Volume;
    end
else
    if Sigma(end) == 0
        F = Volume * 0;
        for i=1:size(F,3)
            F(:,:,i) = imgaussian(Volume(:,:,i), Sigma(1));
        end
    elseif length(Sigma) == 2
        F = imgaussfilt3(Volume,[Sigma(1), Sigma(1), Sigma(2)]);
    elseif length(Sigma) == 3
        F = imgaussfilt3(Volume,[Sigma(1), Sigma(2), Sigma(3)]);
    else
        F=Volume;
    end
end
Dz=gradient3(F,'z');
% temp=Dz(20,20,100)*255
% temp=Dzz(10,10,99)*255
% (F(10,10,100)-F(10,10,99))*255
Dzz=(gradient3(Dz,'z'));
clear Dz;

Dy=gradient3(F,'y');
Dyy=(gradient3(Dy,'y'));
Dyz=(gradient3(Dy,'z'));
clear Dy;

Dx=gradient3(F,'x');
Dxx=(gradient3(Dx,'x'));
Dxy=(gradient3(Dx,'y'));
Dxz=(gradient3(Dx,'z'));
clear Dx;

function D = gradient3(F,option)
% This function does the same as the default matlab "gradient" function
% but with one direction at the time, less cpu and less memory usage.
%
% Example:
%
% Fx = gradient3(F,'x');
[k,l,m] = size(F);
D  = zeros(size(F),class(F)); 
switch lower(option)
case 'x'
    % Take forward differences on left and right edges
    D(1,:,:) = (F(2,:,:) - F(1,:,:));
    D(k,:,:) = (F(k,:,:) - F(k-1,:,:));
    % Take centered differences on interior points
    D(2:k-1,:,:) = (F(3:k,:,:)-F(1:k-2,:,:))/2;
case 'y'
    D(:,1,:) = (F(:,2,:) - F(:,1,:));
    D(:,l,:) = (F(:,l,:) - F(:,l-1,:));
    D(:,2:l-1,:) = (F(:,3:l,:)-F(:,1:l-2,:))/2;
case 'z'
    D(:,:,1) = (F(:,:,2) - F(:,:,1));
    D(:,:,m) = (F(:,:,m) - F(:,:,m-1));
    D(:,:,2:m-1) = (F(:,:,3:m)-F(:,:,1:m-2))/2;
otherwise
    disp('Unknown option')
end
        

