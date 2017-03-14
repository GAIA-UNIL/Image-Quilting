function [y_template distances2_scaled s] = cond_dist_var_layer_temp...
    ( s,X,Y,Y_temp_store,pos,tilesize,D,w,w_v,startI,endI,startJ,endJ,distances1,overlap,nbreplicates)

% Function to compute the Template Splitting Conditioning distances considering all variables and TI layers

% X = the input image
% D = Conditioning values
% startI = starting position of each tile in i direction
% endI = end position of each tile in i direction
% startJ = starting position of each tile in j direction
% endJ = end position of each tile in j direction
% distances1 = The distance matrix computing the sum of squared distances from each tile to the overlap region

distances2 = zeros(size(X,1)-tilesize+1, size(X,2)-tilesize+1, size(X,3));
distances2_scaled=zeros(size(distances2,1), size(distances2,2), size(X,4));
for t=1:size(X,4) % No of TI layers
    for v=1:size(X,3)  % No of variables
%         Data=D{1,v};
        Data_check=sum(D{1,v},1);
        if (Data_check(1,3)~=0)
            x=D{1,v}(:,1); % x coordinates of conditioning points
            y=D{1,v}(:,2); % y coordinates of conditioning points
            [y_template dist_2 s] = cond_dist_temp(s,X(:,:,v,t),Y,Y_temp_store,pos,D{1,v},w,x,y,startI,endI,startJ,endJ,...
                tilesize,distances1,overlap,nbreplicates);
            distances2(:,:,v) = w_v(v)*dist_2(:,:,1);
        else
            y_template=0;
        end
    end
    distances2_scaled(:,:,t) = sum(distances2,3);
end
