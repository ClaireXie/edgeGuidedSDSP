function rpatch = extractRotatedPatch(img, center, width, height, angle)
% To extract a rotated patch from image
%
% Input:
% img - image (gray scale or color)
% center - the pixel location of the center of rectangle 1x2 vector
% width - width of rectangle in pixels
% height - height of rectangle in pixels
% angle - rotation angle of rectangle in degrees
% output:
% rpatch = rotated patch
%
% Tian Cao, University of North Carolina at Chapel Hill (tiancao@cs.unc.edu)

theta = angle/180*pi;

% rotate coordinate
[x,y]   = meshgrid(1:size(img,1), 1:size(img,2));
rotatex = x(center(2)-floor(height/2):center(2)+floor(height/2), ...
    center(1)-floor(width/2):center(1)+floor(width/2));
rotatey = y(center(2)-floor(height/2):center(2)+floor(height/2), ...
    center(1)-floor(width/2):center(1)+floor(width/2));
coords   = [rotatex(:)'-center(1); rotatey(:)'-center(2)];
roatemat = [cos(theta) sin(theta);...
    -sin(theta) cos(theta)]; % generate rotation matrix
rotatedcoords = roatemat*coords;

if size(img, 3) == 1
    rpatch  = interp2(double(img), rotatedcoords(1,:)+center(2), ...
        rotatedcoords(2,:)+center(1), 'linear');
    rpatch  = reshape(rpatch, [height width]);
else
    rpatch(:,:,1)  = interp2(double(img(:,:,1)), rotatedcoords(1,:)+center(1), ...
        rotatedcoords(2,:)+center(2), 'linear');
    rpatch(:,:,2)  = interp2(double(img(:,:,2)), rotatedcoords(1,:)+center(1), ...
        rotatedcoords(2,:)+center(2), 'linear');
    rpatch(:,:,3)  = interp2(double(img(:,:,3)), rotatedcoords(1,:)+center(1), ...
        rotatedcoords(2,:)+center(2), 'linear');
    rpatch  = reshape(rpatch, [height width 3]);
end


