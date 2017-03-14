function [imresize3D] = my_resize(in, nx, ny, nz)
% This function will convert the input 3D inrix (im) into a resized 3D
% inrix witht he following dimensions:
% new dimensions: nx, ny, nz

[x y z]=ndgrid(linspace(1,size(in,1),nx),linspace(1,size(in,2),ny),linspace(1,size(in,3),nz));
imresize3D=interp3(in,y,x,z);

end