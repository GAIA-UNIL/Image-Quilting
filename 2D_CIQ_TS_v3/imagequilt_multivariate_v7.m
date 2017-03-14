function Y = imagequilt_multivariate_v7(X, m, n, D, w, tilesize, overlap, nbreplicates, w_v,temp_split)
% Performs the Image Quilting algorithm on the input with Data Conditioning

% Inputs
% X:  The source image to be used in synthesis
% Tilesize: The dimensions of each square tile
% m:  The number of tiles to be placed in the output image, in X dimension
% n:  The number of tiles to be placed in the output image, in Y dimension
% D: Conditioning values
% w: Weightage of conditioning
% overlap: The amount of overlap to allow between pixels (def: 1/6 tilesize)
% nbreplicates: Number of replicas considered to find the best possible patch
% w_v: Weightage of different variables

%%
X = double(X);  % convert X to double precision X
% [P Q Nvar] = size (D);
Nvar = size(X,3);
Nlayer = size(X,4);

if( overlap >= tilesize )
    error('Overlap must be less than tilesize');
end;

destsize_x = m * tilesize - (m-1) * overlap;
destsize_y = n * tilesize - (n-1) * overlap;
Y = zeros(destsize_x, destsize_y, Nvar);   % Y = the size of output
pos = zeros(1,3);

nx=size(Y,1);

% % fig = figure('position',[100 100 850 600]);
% avifile='avifile.avi';
% 
% vidObj = VideoWriter(avifile);
% vidObj.FrameRate = 3;
% vidObj.Quality = 100;
% open(vidObj);
% video=0;

%% Determining starting and ending points of each tile in both i, j direction
for i=1:m,
    pos(1) = pos(1) + 1;
    pos(2) = 0;
    for j=1:n,
        pos
        Total = [m n]
        pos(2) = pos(2) + 1;
        startI(1) = (i-1)*tilesize - (i-1) * overlap + 1; % starting position of each tile in i direction
        startJ(1) = (j-1)*tilesize - (j-1) * overlap + 1; % end position of each tile in i direction
        endI(1) = startI + tilesize -1; % starting position of each tile in j direction
        endJ(1) = startJ + tilesize -1; % end position of each tile in j direction
        % Determining the distances from each tile to the overlap region
        distances1 = zeros( size(X,1)-tilesize+1, size(X,2)-tilesize+1, Nlayer); % initialize the distances1 matrix
        %         distances2 = zeros(size(X,1)-tilesize+1, size(X,2)-tilesize+1, Nvar); % initialize the distances2 matrix
        [y_template distances2_scaled] = conditioning_distance_var_layer_v2(X,Y,pos,tilesize,D,w,w_v,...
            startI,endI,startJ,endJ,distances1,overlap,nbreplicates,temp_split);
        
        if sum(y_template(:)~=0)~=0
            y_template = y_template(1:endI-startI+1, 1:endJ-startJ+1, 1:Nvar);
            Y(startI:endI, startJ:endJ, 1:Nvar) = y_template;
        else
            %%
            % Compute the sum of squared distances in j direction between X and Y(startI:endI, startJ:startJ+overlap-1, 1:3)
            % for each possible overlap of Y(startI:endI, startJ:startJ+overlap-1, 1:3) on X.
            if( j > 1 )
                z = zeros(size(X,1)-tilesize+1, size(X,2)-overlap+1, size(X,3));
                Z = zeros(size(X,1)-tilesize+1, size(X,2)-overlap+1, size(X,3), size(X,4));
                distances1 = ssd_v4( X, Y(startI:endI, startJ:startJ+(overlap-1), 1:Nvar), w_v, z, Z );
                distances1 = distances1(1:end, 1:end-tilesize+overlap, 1:Nlayer);   % considering only the overlapping region
            end
            
            %%
            % Compute the sum of squared distances in i direction between X and Y(startI:startI+overlap-1, startJ:endJ, 1:3)
            % for each possible overlap of Y(startI:startI+overlap-1, startJ:endJ, 1:3) on X.
            if( i > 1 )
                z = zeros(size(X,1)-overlap+1, size(X,2)-tilesize+1, size(X,3));
                Z = zeros(size(X,1)-overlap+1, size(X,2)-tilesize+1, size(X,3), size(X,4));
                Z = ssd_v4( X, Y(startI:startI+(overlap-1), startJ:endJ, 1:Nvar), w_v, z, Z );
                Z = Z(1:end-tilesize+overlap, 1:end, 1:Nlayer);   % considering only the overlapping region
                if( j > 1 )
                    distances1 = distances1 + Z;  % considering the overlapping regions in both i & j directions
                else
                    distances1 = Z; % Compute the distances of overlap for i > 1 & j = 1
                end;
            end;
            
            %%
            % If both i > 1 & j > 1, compute the distances of the overlap
            if( i > 1 && j > 1 )
                z = zeros(size(X,1)-overlap+1, size(X,2)-overlap+1, size(X,3));
                Z = zeros(size(X,1)-overlap+1, size(X,2)-overlap+1, size(X,3), size(X,4));
                Z = ssd_v4( X, Y(startI:startI+overlap-1, startJ:startJ+(overlap-1), 1:Nvar), w_v, z, Z );
                Z = Z(1:end-tilesize+overlap, 1:end-tilesize+overlap, 1:Nlayer);
                distances1 = distances1 - Z;
            end;

            %% Rescaling by the nb of pixels in overlap region
            distances1_scaled=distances1./(tilesize*overlap);
            
            %%
            if( i == 1 && j == 1 )
                sub(1) = round(rand*(size(X,1)-tilesize));
                sub(2) = round(rand*(size(X,2)-tilesize));
                %                 sub(3) = round(rand*(size(X,3)));
                sub(3) = 1;
                if (sub(1) <= 0)
                    sub(1) = 1;
                end
                if (sub(2) <= 0)
                    sub(2) = 1;
                end
                if (sub(3) <= 0)
                    sub(3) = 1;
                end
            else % For rest of the blocks
                % Compute the distances considering the weights
                distances = (1-w)*distances1_scaled + w*distances2_scaled;
                % sort distances
                [~,distances_index]=sort(distances(:),'ascend');
                % draw one in nbreplicates of the best distances
                idx = distances_index(ceil(rand(1)*nbreplicates));
                [sub(1), sub(2), sub(3)] = ind2sub(size(distances1), idx);
            end
            
            %%
            % Initialize the mask to all ones
            M = ones(tilesize, tilesize);
            
            % If we have a left overlap in j direction
            if( j > 1 )
                %Compute the SSD in the border region
                Err = ( X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+overlap-1, :, sub(3)) - Y(startI:endI, startJ:startJ+overlap-1, :) ).^2;
                
                % Calculate combined error for all the variables considering the weightage
                E = zeros(size(Err));
                for k=1:size(X,3), % No of variables
                    E(:,:,k) = w_v(k)*Err(:,:,k);
                end
                E = sum(E,3);
                
                %Compute the mincut array
                C = mincut_corrected_kashif(E.^2, 0);
                %Compute the mask and write to the destination
                M(1:end, 1:overlap) = double(C >= 0);
            end;
            
            %We have a top overlap
            if( i > 1 )
                %Compute the SSD in the border region
                Err = ( X(sub(1):sub(1)+overlap-1, sub(2):sub(2)+tilesize-1, :, sub(3)) - Y(startI:startI+overlap-1, startJ:endJ, :) ).^2;
                
                % Calculate combined error for all the variables considering the weightage
                E = zeros(size(Err));
                for k=1:size(X,3), % No of variables
                    E(:,:,k) = w_v(k)*Err(:,:,k); 
                end
                E = sum(E,3);
                
                %Compute the mincut array
                C = mincut_corrected_kashif(E.^2, 1);
                %Compute the mask and write to the destination
                M(1:overlap, 1:end) = M(1:overlap, 1:end) .* double(C >= 0);
            end;
            
            if( i == 1 && j == 1) % For the first block
                Y(startI:endI, startJ:endJ, 1:Nvar) = X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, 1:Nvar, sub(3));
            else
                %Write to the destination using the mask
                Y(startI:endI, startJ:endJ, :) = filtered_write(Y(startI:endI, startJ:endJ, :), ...
                    X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, :, sub(3)), M);
            end
        end
        
% %         if (pos(1)>=7 && pos(2)>=4)
%         figure(100);hold on
% %         subplot(1,2,1);hold on
%         imagesc(flipud(Y))
%         
%         scatter(D{1,1}(:,2) , nx-D{1,1}(:,1) , ones(size(D{1,1}(:,1)))*30 , D{1,1}(:,3));
%         plot(D{1,1}(:,2), nx-D{1,1}(:,1), 'r.');
%         axis equal tight
%         colormap gray
%         
%         if (pos(2)==5 && pos(1)==4)
%             % %         if (pos(2)<=4)
%             video=1;
%         end
% %         if (pos(2)==5 && pos(1)==2)
% %             % %         if (pos(2)<=4)
% %             video=1;
% %         end
%         if video==1
%             set(gca,'NextPlot','replacechildren');
%             currFrame = getframe;
%             writeVideo(vidObj,currFrame);
%         end
%         
%         'stop';
%         %         end
        
    end
end

% close(vidObj);

function A = filtered_write(A, B, M)
Nvar=size(A,3);
for i = 1:Nvar
    A(:, :, i) = A(:,:,i) .* (M == 0) + B(:,:,i) .* (M == 1);
end;
