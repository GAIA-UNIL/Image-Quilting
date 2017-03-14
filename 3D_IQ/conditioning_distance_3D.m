function distances2 = conditioning_distance_3D(X,var_type,P,D,x,y,z,startI,endI,startJ,endJ,startK,endK,tilesize,distances1_scaled)

% Function to compute the distances for conditioning
% X = the input image
% P = Number of Conditioning points
% D = Conditioning values
% x = x coordinates of conditioning points
% y = y coordinates of conditioning points
% startI = starting position of each tile in i direction
% endI = end position of each tile in i direction
% startJ = starting position of each tile in j direction
% endJ = end position of each tile in j direction
% n = The number of tiles to be placed in the output image, in each dimension
% distances1_scaled = The scaled distance matrix computing the sum of squared distances from each tile to the overlap region

%tic
nb_data_in_patch = 0;
c=0;
Dc=[];
for k=1:P
    if (x(k) >= startI && x(k) <= endI && y(k) >= startJ && y(k) <= endJ && z(k) >= startK && z(k) <= endK)
        c=c+1;
        Dc(c,:) = D(k,:);
        nb_data_in_patch = nb_data_in_patch + 1;
    end;
end;
x0 = startI;
y0 = startJ;
z0 = startK;
if nb_data_in_patch>0
    xc=Dc(:,1);
    yc=Dc(:,2);
    zc=Dc(:,3);
    for kk = 1:nb_data_in_patch
        xc(kk) = xc(kk)-x0+1;
        yc(kk) = yc(kk)-y0+1;
        zc(kk) = zc(kk)-z0+1;
    end
    if var_type==1 % Represents catagorical variable
        for ii = 1:(size(X,1)-tilesize(1)+1)
            for jj = 1:(size(X,2)-tilesize(2)+1)
                for kk = 1:(size(X,3)-tilesize(3)+1)
                    ec=0;
                    for ll = 1:nb_data_in_patch
                        ec = ec + (X((xc(ll)+ii-1),(yc(ll)+jj-1),(zc(ll)+kk-1))~=Dc(ll,4));
                    end;
                    distances2(ii,jj,kk) = ec;
                end;
            end;
        end
    elseif (var_type==2) % Represents continuous variable
        for ii = 1:(size(X,1)-tilesize(1)+1)
            for jj = 1:(size(X,2)-tilesize(2)+1)
                for kk = 1:(size(X,3)-tilesize(3)+1)
                    ec=0;
                    for ll = 1:nb_data_in_patch
                        ec = ec + (X((xc(ll)+ii-1),(yc(ll)+jj-1),(zc(ll)+kk-1))-Dc(ll,4)).^2;
                    end;
                    distances2(ii,jj,kk) = ec;
                end;
            end;
        end
    end
    distances2 = distances2(1:size(distances1_scaled,1), 1:size(distances1_scaled,2) ,1:size(distances1_scaled,3));
    distances2 = sqrt(distances2);
    %distances2_scaled = (distances2)./(max(X(:))*nb_data_in_patch);  %%%HERE WHY DIVIDE BY max(X(:))???
    distances2 = distances2./nb_data_in_patch; 
else
    distances2 = zeros(size(distances1_scaled));
end
%toc