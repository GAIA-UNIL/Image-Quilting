function Surfacecut3D=DP_3D_surface_v2_function_Kashif_v2(e,dir,showslices)

% Function to find the 3D Surface mincut of lowest values through a 3D bloc

% e: matrix containing the bloc of continuous values where to find the surface cut
% dir: direction of the surface cut: 1=X, 2=Y, 3=Z
% showslices = 1 for viewing the slices in specified direction
%            = Other values for not viewing the slices in specified direction

% VARIABLES USED:
% e: error matrix
% e2: rotated error matrix
% ee: cumulative minimum error with nan crust
% E: cumulative minimum error without nan crust
% C3D: the surface cut of the 1D shortest path in the 3D block

[x y z]=size(e);

%% rotate block to account for direction of trace
if dir==3
    e2=zeros(y,z,x);
    for i=1:x
        e2(:,:,i)=squeeze(e(i,:,:));
    end
elseif dir==1
    e2=zeros(z,x,y);
    for i=1:y
%         e2(:,:,i)=squeeze(e(i,:,:));
        e2(:,:,i)=squeeze(e(:,i,:))';
%         e2(:,:,i)=squeeze(e(:,:,i));
    end
elseif dir==2
    e2=e;
end

% rotated dimensions
[x y z]=size(e2);

% %%
% figure(4);clf
% ViewGrid(e2)
% axis equal tight
% colormap gray
% view(-50,50) 

%% add crust of nan
nn=nan(x+2,y+2,z+2);
nn(2:end-1,2:end-1,2:end-1)=e2;
e2=nn;

%% compute cumulative min errors in 3D
ee=nan(size(e2));
ee(:,:,2)=e2(:,:,2);
for i=3:z+1
    for j=2:y+1
        for k=2:x+1
            patch2D=ee(k-1:k+1,j-1:j+1,i-1);
            ee(k,j,i)=e2(k,j,i)+min(patch2D(isfinite(patch2D)));
        end
    end
end
E=ee(2:end-1,2:end-1,2:end-1);

%% 2D DP for the last layer
C0 = mincut_corrected_kashif(E(:,:,end),2);

C3D=zeros(size(E));
C3D(:,:,end)=C0;

%% go back along Z checking only the potential locations that do not break the plane
%backtrace, which are around the divide of the previous Z slice
% to do this, we first do for each layer a mask which has higher values in all excluded locations (double the error value),
% then we do a 2D DP in this mask. Since the mincut cannot go through the
% high values, it will stick to the non-masked areas.

for k=z-1:-1:1
    
    mask=E(:,:,k)*2;
    for i=x:-1:1
        ind=find(C3D(i,:,k+1)==0);
        
        if ind==1
            mask(i,ind:ind+1)=E(i,ind:ind+1,k);
        elseif ind==y
            mask(i,ind-1:ind)=E(i,ind-1:ind,k);
        else
            mask(i,ind-1:1:ind+1)=E(i,ind-1:1:ind+1,k);
        end
    end
    C3D(:,:,k) = mincut_corrected_kashif(mask,2);
end

e2=e2(2:end-1,2:end-1,2:end-1);
Surfacecut3D=C3D;

%% back-rotate final trace
T2=zeros(size(e,1),size(e,2),size(e,3));
if dir==3
    for i=1:size(e,3)
        T2(:,:,i)=rot90(squeeze(Surfacecut3D(:,i,:)));
    end
elseif dir==1
    for i=1:size(e,3)
        T2(:,:,i)=squeeze(Surfacecut3D(i,:,:));
%         T2(:,:,i)=squeeze(Surfacecut3D(:,i,:));
%         T2(:,:,i)=squeeze(Surfacecut3D(:,:,i));
    end
elseif dir==2
    T2=Surfacecut3D;
end
Surfacecut3D=T2;

%% show slice per slice
if showslices==1
    for i=z:-1:1
        figure(1);clf;colormap gray
        
        subplot(2,1,1)
        imagesc(e2(:,:,i))
        colorbar
        
        subplot(2,1,2)
        imagesc(Surfacecut3D(:,:,i))
        colorbar
    end
end
