function [y_template distances2 s] = cond_dist_temp...
    (s,X,Y,Y_temp_store,pos,D,w,x,y,startI,endI,startJ,endJ,tilesize,distances1_scaled,overlap,nbreplicates)

% Function to compute the Template Splitting Conditioning distances

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
y_template=0;
nb_data_in_patch = 0;
c=0;
Dc=[];
for k=1:size(D,1)
    if (x(k) >= startI(end) && x(k) <= startI(end)+tilesize-1 &&...
            y(k) >= startJ(end) && y(k) <= startJ(end)+tilesize-1)
        c=c+1;
        Dc(c,:) = D(k,:);
        nb_data_in_patch = nb_data_in_patch + 1;
    end;
end;
distances2 = zeros(size(distances1_scaled));

if nb_data_in_patch == 1
    xc=Dc(:,1);
    yc=Dc(:,2);
    xc = xc-startI(end)+1;
    yc = yc-startJ(end)+1;
    for ii = 1:(size(X,1)-tilesize+1)
        for jj = 1:(size(X,2)-tilesize+1)
            ec = (X((xc+ii-1),(yc+jj-1))-Dc(1,3)).^2;
            distances2(ii,jj) = ec;
        end;
    end;
    
elseif nb_data_in_patch > 1
    xc=Dc(:,1);
    yc=Dc(:,2);
    for kk = 1:nb_data_in_patch
        xc(kk) = xc(kk)-startI(end)+1;
        yc(kk) = yc(kk)-startJ(end)+1;
    end
    for ii = 1:(size(X,1)-tilesize+1)
        for jj = 1:(size(X,2)-tilesize+1)
            ec=0;
            for kk = 1:nb_data_in_patch
                ec = ec + (X((xc(kk)+ii-1),(yc(kk)+jj-1))-Dc(kk,3)).^2;
            end;
            distances2(ii,jj) = ec;
        end;
    end;
    
    distances2 = distances2(1:size(distances1_scaled,1), 1:size(distances1_scaled,2));
    distances2 = sqrt(distances2);
    distances2 = distances2./nb_data_in_patch;
    
    cand = sort(distances2(:),'ascend');
    cand = cand(1:nbreplicates);
    
    if ( max(cand) > 0 )
        % if ( max(cand) > abs(0.05*(max(D(:))-min(D(:)))) )
        % if ( max(cand) > abs(0.1*(min(D(:,3)))) )
        % if ( max(cand) > 0.01 )
        Data{1,1} = Dc;
        [y_template s]= iq_multi_v7_temp(s,X,Y,Y_temp_store,pos,Data,w,ceil(tilesize/2+overlap/2),ceil(overlap/2),nbreplicates,1,startI,endI,startJ,endJ);
    end
else
    distances2 = zeros(size(distances1_scaled));
end
