function [Trace3D,maxcost]=ConnecTrace_v2(e,dir,nbtraces)

% function to find the path of highest values in a 3D bloc

% e: matrix containing the bloc of continuous values where to find the trace
% dir: direction of the traces: 1=X, 2=Y, 3=Z
% nbtraces: number of traces to find. 1: trace of minimal resistance
%                                     >1: other traces in order of
%                                     increasing resistance.
%                                  max nb of traces: nb of pixels in a 2D plane

%VARIABLES USED:
% e: error matrix
% e2: rotated error matrix
% ee: cumulative minimum error with nans crust
% E: cumulative minimum error without nans crust
% Trace3D: the trace of the 1D shortest path in the 3D block
% maxcost: cost of first trace

%add infinitesimal noise to have only a single max each time
e=e+(rand(size(e))*(max(e(:)-min(e(:))))/100000);

%original dimensions
[x,y,z]=size(e);

%maxcost=0;

%% rotate block to account for direction of trace
if dir==1
    e2=zeros(y,z,x);
    for i=1:x     
        e2(:,:,i)=squeeze(e(i,:,:));
    end
elseif dir==2
    e2=zeros(x,z,y);
    for i=1:y
        e2(:,:,i)=squeeze(e(:,i,:));
    end
elseif dir==3
    e2=e;
end

%rotated dimensions
[x,y,z]=size(e2);

%% add crust of nan
nn=nan(x+2,y+2,z+2);
nn(2:end-1,2:end-1,2:end-1)=e2;
e2=nn;

% compute cumulative min errors in 3D
ee=nan(size(e2));
ee(:,:,2)=e2(:,:,2);
for i=3:z+1
    for j=2:y+1
        for k=2:x+1
            patch2D=ee(k-1:k+1,j-1:j+1,i-1); %patch 3x3 gives us all directions
            ee(k,j,i)=e2(k,j,i)+max(patch2D(isfinite(patch2D)));
        end
    end
end

%remove the crust
E=ee(2:end-1,2:end-1,2:end-1);

%% add crust of nan
nn=nan(x+2,y+2,z+2);
nn(2:end-1,2:end-1,2:end-1)=E;
E=nn;

Trace3D=zeros(size(E,1),size(E,2),size(E,3),nbtraces);
for t=1:nbtraces
    %find the end point of each trace
    m=E(:,:,z);
    s = sort(m(isfinite(m)),'descend');
    I=m==s(t);
    
    if t==1
        maxcost=s(t);
    end
    
    %backtrace
    Trace3D(:,:,end-1,t)=I;
    for k=z:-1:2
        ind=find(Trace3D(:,:,k+1,t)==1);
        [i,j] = ind2sub(size(I),ind);
        m=E(i-1:i+1,j-1:j+1,k);
        I=( E(:,:,k) == max(m(:)) );
        Trace3D(:,:,k,t)=I;
    end
end

%remove the crust
Trace3D=Trace3D(2:end-1,2:end-1,2:end-1,:);

%% back-rotate final trace

T2=zeros(size(e,1),size(e,2),size(e,3),nbtraces);
for t=1:nbtraces
    if dir==1
        for i=1:size(e,3)            
            T2(:,:,i,t)=fliplr(rot90(squeeze(Trace3D(:,i,:,t)),3));
        end
    elseif dir==2
        for i=1:size(e,3)
            T2(:,:,i,t)=squeeze(Trace3D(:,i,:,t));
        end
    elseif dir==3
        T2=Trace3D;
    end
end

Trace3D=T2;

