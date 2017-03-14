function [Y_temp s]= iq_multi_v7_temp(s,X,Y,Y_temp_store,pos,D,w,tilesize,overlap,nbreplicates,w_v,startI,endI,startJ,endJ)
% Performs the Image Quilting algorithm with Tempalte Splitting Conditioning

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
s = s + 1;
X = double(X);  % convert X to double precision X
% [P Q Nvar] = size (D);
Nvar = size(X,3);
Nlayer = size(X,4);

if( overlap >= tilesize )
    error('Overlap must be less than tilesize');
end;

destsize_x = 2 * tilesize - overlap;
destsize_y = 2 * tilesize - overlap;

Y_temp = zeros(destsize_x, destsize_y, Nvar);
% nx=size(Y_temp,1);
pos(s) = 0;
y_template = 0;

% figure(101);clf;
% avifile_1='avifile_1.avi';
% 
% vidObj_1 = VideoWriter(avifile_1);
% vidObj_1.FrameRate = 2;
% vidObj_1.Quality = 100;
% open(vidObj_1);

%% Determining starting and ending points of each tile in both i, j direction
for i=1:2,
    for j=1:2,
        pos(s) = pos(s) + 1;
        startI(s-1) = startI(s-2) + (i-1)*tilesize - (i-1) * overlap; % starting position of each tile in i direction
        startJ(s-1) = startJ(s-2) + (j-1)*tilesize - (j-1) * overlap; % end position of each tile in i direction
        endI(s-1) = startI(s-1) + tilesize -1; % starting position of each tile in j direction
        endJ(s-1) = startJ(s-1) + tilesize -1; % end position of each tile in j direction
        % Determining the distances from each tile to the overlap region
        distances1 = zeros( size(X,1)-tilesize+1, size(X,2)-tilesize+1, Nlayer); % initialize the distances1 matrix
        if (tilesize ~= 2 )
            [y_template distances2_scaled s] = cond_dist_var_layer_temp...
                (s,X,Y,Y_temp_store,pos,tilesize,D,w,w_v,startI,endI,startJ,endJ,distances1,overlap,nbreplicates);
        end
        
        % Adjustments with Y for last row/column filling
        gg = endI(s-1) - size (Y,1);
        if gg > 0
            Y(size(Y,1)+1:endI(s-1), :) = 0;
        end
        
        g = endJ(s-1) - size (Y,2);
        if g > 0
            Y(:, size(Y,2)+1:endJ(s-1)) = 0;
        end
        
        %% For 4th template
        if (pos(s) == 4)
            % Left overlap with Y_template
            z = zeros(size(X,1)-tilesize+1, size(X,2)-overlap+1, size(X,3));
            Z = zeros(size(X,1)-tilesize+1, size(X,2)-overlap+1, size(X,3), size(X,4));
            distances1 = ssd_v4( X, Y_temp(tilesize-overlap+1:2*tilesize-overlap, startJ(s-1)-startJ(s-2)+1:startJ(s-1)-startJ(s-2)+overlap, 1:Nvar), w_v, z, Z );
            distances1 = distances1(1:end, 1:end-tilesize+overlap, 1:Nlayer);   % considering only the overlapping region
            
            % Top overlap with Y_template
            z = zeros(size(X,1)-overlap+1, size(X,2)-tilesize+1, size(X,3));
            Z = zeros(size(X,1)-overlap+1, size(X,2)-tilesize+1, size(X,3), size(X,4));
            Z = ssd_v4( X, Y_temp(startI(s-1)-startI(s-2)+1:startI(s-1)-startI(s-2)+overlap, end-tilesize+1:end, 1:Nvar), w_v, z, Z );
            Z = Z(1:end-tilesize+overlap, 1:end, 1:Nlayer);   % considering only the overlapping region
            distances1 = distances1 + Z;  % considering the overlapping regions in both i & j directions
            
            % Common overlap deduction with Y_template
            z = zeros(size(X,1)-overlap+1, size(X,2)-overlap+1, size(X,3));
            Z = zeros(size(X,1)-overlap+1, size(X,2)-overlap+1, size(X,3), size(X,4));
            Z = ssd_v4( X, Y_temp(startI(s-1)-startI(s-2)+1:startI(s-1)-startI(s-2)+overlap, startJ(s-1)-startJ(s-2)+1:startJ(s-1)-startJ(s-2)+overlap, 1:Nvar), w_v, z, Z );
            Z = Z(1:end-tilesize+overlap, 1:end-tilesize+overlap, 1:Nlayer);
            distances1 = distances1 - Z;
        end
        
        %% For 3rd template
        if (pos(s) == 3)
            if (pos(s-1) > 1 && s > 3)
                % Left overlap with previous Y_temp
                Y_temp_store{1,s-1} = Y_temp_store{1,s-1}(1:endI(s-2)-startI(s-2)+1, 1:endI(s-2)-startI(s-2)+1, 1:Nvar);
                Y11 = Y_temp_store{1,s-1}(end-tilesize+1:end, end-2*overlap+1:end, 1:Nvar);
                z = zeros(size(X,1)-tilesize+1, size(X,2)-2*overlap+1, size(X,3));
                Z = zeros(size(X,1)-tilesize+1, size(X,2)-2*overlap+1, size(X,3), size(X,4));
                distances11 = ssd_v4( X, Y11(1:tilesize, 1:2*overlap, 1:Nvar), w_v, z, Z );
                distances11 = distances11(1:end, 1:end-tilesize+2*overlap, 1:Nlayer);   % considering only the overlapping region
            else
                distances11 = distances1;
            end
            
            % Top overlap with Y_template
            z = zeros(size(X,1)-overlap+1, size(X,2)-tilesize+1, size(X,3));
            Z = zeros(size(X,1)-overlap+1, size(X,2)-tilesize+1, size(X,3), size(X,4));
            distances1 = ssd_v4( X, Y_temp(startI(s-1)-startI(s-2)+1:startI(s-1)-startI(s-2)+overlap, 1:tilesize, 1:Nvar), w_v, z, Z );
            distances1 = distances1(1:end-tilesize+overlap, 1:end, 1:Nlayer);   % considering only the overlapping region
            distances1 = distances11 + distances1;
            
            if (pos(2) > 1 && s == 3) || (pos(2) > 1 && pos(s-1) == 1 && s > 3)
                % Left overlap with Y
%                 if endI(s-1) > size (Y,1)
%                     Y(size(Y,1):endI(s-1), :) = 0;
%                 end
                Y1 = Y(startI(s-1):endI(s-1), startJ(s-1):endJ(s-1), Nvar);
                z = zeros(size(X,1)-tilesize+1, size(X,2)-2*overlap+1, size(X,3));
                Z = zeros(size(X,1)-tilesize+1, size(X,2)-2*overlap+1, size(X,3), size(X,4));
                Z = ssd_v4( X, Y1(1:tilesize, 1:2*overlap, 1:Nvar), w_v, z, Z );
                Z = Z(1:end, 1:end-tilesize+2*overlap, 1:Nlayer);   % considering only the overlapping region
                distances1 = distances1 + Z;  % considering the overlapping regions in both i & j directions
%                 if size (Y,1) ~= size (Y,2)
%                     Y(size(Y,1)+1:endI(s-1), :) = [];
%                 end
            end
            
            if (pos(1) > 1)
                % Top overlap with Y
                Y1 = Y(startI(s-1):endI(s-1), startJ(s-1):endJ(s-1), Nvar);
                z = zeros(size(X,1)-2*overlap+1, size(X,2)-tilesize+1, size(X,3));
                Z = zeros(size(X,1)-2*overlap+1, size(X,2)-tilesize+1, size(X,3), size(X,4));
                Z = ssd_v4( X, Y1(1:2*overlap, 1:tilesize, 1:Nvar), w_v, z, Z );
                Z = Z(1:end-tilesize+2*overlap, 1:end, 1:Nlayer);   % considering only the overlapping region
%                 if (pos(s-1) == 1)
%                     distances1 = Z;
%                 else % (pos(2) > 1)
                    distances1 = distances1 + Z;
%                 end
            end
        end
        
        %% For 2nd template
        if (pos(s) == 2)
            % Left overlap with Y_template
            z = zeros(size(X,1)-tilesize+1, size(X,2)-overlap+1, size(X,3));
            Z = zeros(size(X,1)-tilesize+1, size(X,2)-overlap+1, size(X,3), size(X,4));
            distances1 = ssd_v4( X, Y_temp(1:tilesize, startJ(s-1)-startJ(s-2)+1:startJ(s-1)-startJ(s-2)+overlap, 1:Nvar), w_v, z, Z );
            distances1 = distances1(1:end, 1:end-tilesize+overlap, 1:Nlayer);   % considering only the overlapping region
            
            if (pos(1) > 1)
                % Top overlap with Y
%                 g = endJ(s-1) - size (Y,2);
%                 if g > 1
%                     Y(:, size(Y,2)+1:endJ(s-1)) = 0;
%                 end
                Y1 = Y(startI(s-1):endI(s-1), startJ(s-1):endJ(s-1), Nvar);
                z = zeros(size(X,1)-2*overlap+1, size(X,2)-tilesize+1, size(X,3));
                Z = zeros(size(X,1)-2*overlap+1, size(X,2)-tilesize+1, size(X,3), size(X,4));
                Z = ssd_v4( X, Y1(1:2*overlap, 1:tilesize, 1:Nvar), w_v, z, Z );
                Z = Z(1:end-tilesize+2*overlap, 1:end, 1:Nlayer);   % considering only the overlapping region
                distances1 = distances1 + Z;  % considering the overlapping regions in both i & j directions
%                 if size (Y,1) ~= size (Y,2)
%                     Y(:, end-g+1:end) = [];
%                 end
            end
        end
        
        %% For 1st template
        if (pos(s) == 1)
%             gg = endI(s-1) - size (Y,1);
%             if gg > 1
%                 Y(size(Y,1)+1:endI(s-1), :) = 0;
%             end
            if (pos(s-1) > 1 && s > 3)
                % Left overlap with previous Y_temp
                Y_temp_store{1,s-1} = Y_temp_store{1,s-1}(1:endI(s-2)-startI(s-2)+1, 1:endI(s-2)-startI(s-2)+1, 1:Nvar);
                Y11 = Y_temp_store{1,s-1}(1:tilesize, end-2*overlap+1:end, 1:Nvar);
                z = zeros(size(X,1)-tilesize+1, size(X,2)-2*overlap+1, size(X,3));
                Z = zeros(size(X,1)-tilesize+1, size(X,2)-2*overlap+1, size(X,3), size(X,4));
                distances11 = ssd_v4( X, Y11(1:tilesize, 1:2*overlap, 1:Nvar), w_v, z, Z );
                distances11 = distances11(1:end, 1:end-tilesize+2*overlap, 1:Nlayer);   % considering only the overlapping region
            else
                distances11 = distances1;
            end
            
            if (pos(2) > 1 && s == 3) || (pos(2) > 1 && pos(s-1) == 1 && s > 3) || ( pos(2) > 1 && pos(s-1) == 1 && s > 3)
                % Left overlap with Y
                Y1 = Y(startI(s-1):endI(s-1), startJ(s-1):endJ(s-1), Nvar);
                z = zeros(size(X,1)-tilesize+1, size(X,2)-2*overlap+1, size(X,3));
                Z = zeros(size(X,1)-tilesize+1, size(X,2)-2*overlap+1, size(X,3), size(X,4));
                distances1 = ssd_v4( X, Y1(1:tilesize, 1:2*overlap, 1:Nvar), w_v, z, Z );
                distances1 = distances1(1:end, 1:end-tilesize+2*overlap, 1:Nlayer);   % considering only the overlapping region
                distances1 = distances11 + distances1;
            else
                distances1 = distances11;
            end
            if (pos(1) > 1)
                % Top overlap with Y
                Y1 = Y(startI(s-1):endI(s-1), startJ(s-1):endJ(s-1), Nvar);
                z = zeros(size(X,1)-2*overlap+1, size(X,2)-tilesize+1, size(X,3));
                Z = zeros(size(X,1)-2*overlap+1, size(X,2)-tilesize+1, size(X,3), size(X,4));
                Z = ssd_v4( X, Y1(1:2*overlap, 1:tilesize, 1:Nvar), w_v, z, Z );
                Z = Z(1:end-tilesize+2*overlap, 1:end, 1:Nlayer);   % considering only the overlapping region
%                 if (pos(s-1) == 1)
%                     distances1 = Z;
%                 else % (pos(2) > 1)
                    distances1 = distances1 + Z;
%                 end
            end
            %             if size (Y,1) > size (Y,2)
            %                 Y(end-gg+1:end, :) = [];
            %             end
        end
        
        % Re-adjustments with Y for last row/column filling
        if size (Y,1) < size (Y,2)
            Y(:, end-g+1:end) = [];
            % Y(size(Y,1)+1:endI(s-1), :) = [];
        end
        
        if size (Y,1) > size (Y,2)
            Y(end-gg+1:end, :) = [];
        end
        
        %%
%         %         if (pos(1)>=7 && pos(2)>=4)
%                 figure(4);clf;hold on
%                 subplot(1,2,1);
%                 imagesc(flipud(distances1))
%                 axis equal tight
%                 subplot(1,2,2);
%                 imagesc(flipud(distances2_scaled));pos
%                 axis equal tight
%                 %                 subplot(2,2,3);
%                 %                 imagesc(flipud(Y1))
%                 %                 axis equal tight
%         %         subplot(2,2,4);
%         %         imagesc(distances2_scaled);pos
%         %         axis equal tight
%                 'stop';
%         %         end
        
        %%
        if y_template~=0
            if (pos(s) == 1)
                Y_temp(1:tilesize, 1:tilesize, 1:Nvar) = y_template(1:tilesize, 1:tilesize, 1:Nvar);
            elseif (pos(s) == 2)
                Y_temp(1:tilesize, tilesize-overlap+1:2*tilesize-overlap, 1:Nvar) = y_template(1:tilesize, 1:tilesize, 1:Nvar);
            elseif (pos(s) == 3)
                Y_temp(tilesize-overlap+1:2*tilesize-overlap, 1:tilesize, 1:Nvar) = y_template(1:tilesize, 1:tilesize, 1:Nvar);
            else
                Y_temp(tilesize-overlap+1:2*tilesize-overlap, tilesize-overlap+1:2*tilesize-overlap, 1:Nvar) = y_template(1:tilesize, 1:tilesize, 1:Nvar);
            end
        else
            
            % Rescaling by the nb of pixels in overlap region
            distances1_scaled=distances1./(tilesize*overlap); %%%HERE
            
            %%
            if (tilesize == 2 )
                distances2_scaled = zeros(size(distances1_scaled)); % initialize the distances2 matrix
            end
            
            % Compute the distances considering the weights
            distances = (1-w)*distances1_scaled + w*distances2_scaled;
            % sort distances
            [~,distances_index]=sort(distances(:),'ascend');
            % draw one in nbreplicates of the best distances
            idx = distances_index(ceil(rand(1)*nbreplicates));
            [sub(1), sub(2), sub(3)] = ind2sub(size(distances), idx);
            
            %%
            % Initialize the mask to all ones
            M = ones(tilesize, tilesize);
            
            %% For 4th template
            if (pos(s) == 4)
                % Left overlap with Y_template
                % Compute the SSD in the border region
                Err = ( X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+overlap-1, :, sub(3))...
                    - Y_temp(tilesize-overlap+1:2*tilesize-overlap, startJ(s-1)-startJ(s-2)+1:startJ(s-1)-startJ(s-2)+overlap, :) ).^2;
                % Calculate combined error for all the variables considering the weightage
                E = zeros(size(Err));
                for k=1:size(X,3), % No of variables
                    E(:,:,k) = w_v(k)*Err(:,:,k);
                end
                E = sum(E,3);
                
                %Compute the mincut array
                C = mincut_corrected_kashif(E, 0);
                %Compute the mask and write to the destination
                M(1:end, 1:overlap) = double(C >= 0);
                %                 M(1:end, 1:2*overlap) = double(C >= 0);
                
                % Top overlap with Y_template
                Err = ( X(sub(1):sub(1)+overlap-1, sub(2):sub(2)+tilesize-1, :, sub(3))...
                    - Y_temp(startI(s-1)-startI(s-2)+1:startI(s-1)-startI(s-2)+overlap, tilesize-overlap+1:2*tilesize-overlap, :) ).^2;
                
                % Calculate combined error for all the variables considering the weightage
                E = zeros(size(Err));
                for k=1:size(X,3), % No of variables
                    E(:,:,k) = w_v(k)*Err(:,:,k);
                end
                E = sum(E,3);
                
                %Compute the mincut array
                C = mincut_corrected_kashif(E, 1);
                %Compute the mask and write to the destination
                M(1:overlap, 1:end) = M(1:overlap, 1:end) .* double(C >= 0);
                
                %%
                Y_temp(end-tilesize+1:end, end-tilesize+1:end, :) = filtered_write(Y_temp(end-tilesize+1:end, end-tilesize+1:end, :), ...
                    X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, :, sub(3)), M);
            end
            
            %% For 3rd template
            if (pos(s) == 3)
                % Top overlap with Y_template
                Err = ( X(sub(1):sub(1)+overlap-1, sub(2):sub(2)+tilesize-1, :, sub(3))...
                    - Y_temp(startI(s-1)-startI(s-2)+1:startI(s-1)-startI(s-2)+overlap, 1:tilesize, :) ).^2;
                % Calculate combined error for all the variables considering the weightage
                E = zeros(size(Err));
                for k=1:size(X,3), % No of variables
                    E(:,:,k) = w_v(k)*Err(:,:,k);
                end
                E = sum(E,3);
                
                %Compute the mincut array
                C = mincut_corrected_kashif(E, 1);
                %Compute the mask and write to the destination
                M(1:overlap, 1:end) = double(C >= 0);
                
                if (pos(s-1) > 1 && s > 3)
                    % Left overlap with previous Y_temp
                    % Compute the SSD in the border region
                    Err = ( X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+2*overlap-1, :, sub(3))...
                        - Y11(1:tilesize, 1:2*overlap, :) ).^2;
                    % Calculate combined error for all the variables considering the weightage
                    E = zeros(size(Err));
                    for k=1:size(X,3), % No of variables
                        E(:,:,k) = w_v(k)*Err(:,:,k);
                    end
                    E = sum(E,3);
                    
                    %Compute the mincut array
                    C = mincut_corrected_kashif(E, 0);
                    %Compute the mask and write to the destination
                    M(1:end, 1:2*overlap) = M(1:end, 1:2*overlap) .* double(C >= 0);
                end
                
                if (pos(2) > 1 && s == 3)
                    % Left overlap with Y
                    % Compute the SSD in the border region
                    Err = ( X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+2*overlap-1, :, sub(3))...
                        - Y1(end-tilesize+1:end, 1:2*overlap, :) ).^2;
                    % Calculate combined error for all the variables considering the weightage
                    E = zeros(size(Err));
                    for k=1:size(X,3), % No of variables
                        E(:,:,k) = w_v(k)*Err(:,:,k);
                    end
                    E = sum(E,3);
                    
                    %Compute the mincut array
                    C = mincut_corrected_kashif(E, 0);
                    %Compute the mask and write to the destination
                    M(1:end, 1:2*overlap) = M(1:end, 1:2*overlap) .* double(C >= 0);
                end
                
                %%
                if ( pos(1) > 1 )
                    % Top overlap with Y_template and Left overlap with Y
                    if s == 3  || (pos(s-1) == 1 && s > 3)
                        Y_crop = Y1(end-tilesize+1:end, 1:tilesize, :);
                    else
                        Y_crop = zeros(tilesize,tilesize,Nvar);
                        Y_crop(1:tilesize, 1:2*overlap, :) = Y11;
                    end
                    Y_crop(1:overlap, 1:tilesize, :) = Y_temp(startI(s-1)-startI(s-2)+1:startI(s-1)-startI(s-2)+overlap, 1:tilesize, :);
                    
                    Y_temp(tilesize-overlap+1:2*tilesize-overlap, 1:tilesize, :) = filtered_write(Y_crop(1:tilesize, 1:tilesize, :), ...
                        X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, :, sub(3)), M);
                else % ( pos(2) == 1 )
                    % Top overlap with Y_template
                    Y_temp(tilesize-overlap+1:2*tilesize-overlap, 1:tilesize, :) = filtered_write(Y_temp(tilesize-overlap+1:2*tilesize-overlap, 1:tilesize, :), ...
                        X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, :, sub(3)), M);
                end
            end
            
            %% For 2nd template
            if (pos(s) == 2)
                % Left overlap with Y_template
                Err = ( X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+overlap-1, :, sub(3))...
                    - Y_temp(1:tilesize, startJ(s-1)-startJ(s-2)+1:startJ(s-1)-startJ(s-2)+overlap, :) ).^2;
                
                % Calculate combined error for all the variables considering the weightage
                E = zeros(size(Err));
                for k=1:size(X,3), % No of variables
                    E(:,:,k) = w_v(k)*Err(:,:,k);
                end
                E = sum(E,3);
                
                %Compute the mincut array
                C = mincut_corrected_kashif(E, 1);
                %Compute the mask and write to the destination
                M(1:end, 1:overlap) = double(C >= 0);
                
                if (pos(1) > 1)
                    % Top overlap with Y
                    % Compute the SSD in the border region
                    Err = ( X(sub(1):sub(1)+2*overlap-1, sub(2):sub(2)+tilesize-1, :, sub(3))...
                        - Y1(1:2*overlap, end-tilesize+1:end, :) ).^2;
                    
                    % Calculate combined error for all the variables considering the weightage
                    E = zeros(size(Err));
                    for k=1:size(X,3), % No of variables
                        E(:,:,k) = w_v(k)*Err(:,:,k);
                    end
                    E = sum(E,3);
                    
                    %Compute the mincut array
                    C = mincut_corrected_kashif(E, 0);
                    %Compute the mask and write to the destination
                    M(1:2*overlap, 1:end) = M(1:2*overlap, 1:end) .* double(C >= 0);
                end
                
                %%
                if ( pos(1) > 1 )
                    % Left overlap with Y_template and Top overlap with Y
                    Y_crop = Y1(1:tilesize, end-tilesize+1:end, :);
                    Y_crop(1:tilesize, 1:overlap, :) = Y_temp(1:tilesize, startJ(s-1)-startJ(s-2)+1:startJ(s-1)-startJ(s-2)+overlap, :);
                    Y_temp(1:tilesize, tilesize-overlap+1:2*tilesize-overlap, :) = filtered_write(Y_crop(1:tilesize, 1:tilesize, :), ...
                        X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, :, sub(3)), M);
                else % ( pos(1) == 1 )
                    % Left overlap with Y_template
                    Y_temp(1:tilesize, tilesize-overlap+1:2*tilesize-overlap, :) = filtered_write(Y_temp(1:tilesize, tilesize-overlap+1:2*tilesize-overlap, :), ...
                        X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, :, sub(3)), M);
                end
            end
            
            %% For 1st template
            if (pos(s) == 1)
                %                 if (pos(2) > 1)  % (pos(1) == 1 && pos(2) > 1)
                if (pos(s-1) > 1 && s > 3)
                    % Left overlap with previous Y_temp
                    % Compute the SSD in the border region
                    Err = ( X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+2*overlap-1, :, sub(3))...
                        - Y11(1:tilesize, 1:2*overlap, :) ).^2;
                    % Calculate combined error for all the variables considering the weightage
                    E = zeros(size(Err));
                    for k=1:size(X,3), % No of variables
                        E(:,:,k) = w_v(k)*Err(:,:,k);
                    end
                    E = sum(E,3);
                    
                    %Compute the mincut array
                    C = mincut_corrected_kashif(E, 0);
                    %Compute the mask and write to the destination
                    %                     M(1:end, 1:overlap) = double(C >= 0);
                    M(1:end, 1:2*overlap) = double(C >= 0);
                end
                
                if (pos(2) > 1 && s == 3)  % (pos(1) == 1 && pos(2) > 1)
                    %                 if (pos(s-1) > 1)
                    % Left overlap with Y
                    % Compute the SSD in the border region
                    Err = ( X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+2*overlap-1, :, sub(3))...
                        - Y1(1:tilesize, 1:2*overlap, :) ).^2;
                    % Calculate combined error for all the variables considering the weightage
                    E = zeros(size(Err));
                    for k=1:size(X,3), % No of variables
                        E(:,:,k) = w_v(k)*Err(:,:,k);
                    end
                    E = sum(E,3);
                    
                    %Compute the mincut array
                    C = mincut_corrected_kashif(E, 0);
                    %Compute the mask and write to the destination
                    M(1:end, 1:2*overlap) = double(C >= 0);
                end
                
                if (pos(1) > 1) % (pos(2) == 1 && pos(1) > 1)
                    %                 if (pos(s-2) > 1)
                    % Top overlap with Y
                    % Compute the SSD in the border region
                    Err = ( X(sub(1):sub(1)+2*overlap-1, sub(2):sub(2)+tilesize-1, :, sub(3))...
                        - Y1(1:2*overlap, 1:tilesize, :) ).^2;
                    
                    % Calculate combined error for all the variables considering the weightage
                    E = zeros(size(Err));
                    for k=1:size(X,3), % No of variables
                        E(:,:,k) = w_v(k)*Err(:,:,k);
                    end
                    E = sum(E,3);
                    
                    %Compute the mincut array
                    C = mincut_corrected_kashif(E, 0);
                    %Compute the mask and write to the destination
                    if (pos(2) == 1)
                        M(1:2*overlap, 1:end) = double(C >= 0);
                    else % pos(2) > 1
                        M(1:2*overlap, 1:end) = M(1:2*overlap, 1:end) .* double(C >= 0);
                    end
                end
                
                %%
                if ( pos(2) == 1 && pos(1) == 1 )
                    % No overlap for the first block
                    Y_temp(1:tilesize, 1:tilesize, 1:Nvar) = ...
                        X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, 1:Nvar, sub(3));
                elseif ( pos(2) > 1 && pos(1) == 1 && s == 3) || ( pos(2) > 1 && pos(s-1) == 1 && s > 3)
                    % Left overlap with Y
                    Y_temp(1:tilesize, 1:tilesize, :) = filtered_write(Y1(1:tilesize, 1:tilesize,...
                        :), X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, :, sub(3)), M);
%                 elseif ( pos(2) > 1 && pos(1) == 1 && s > 3)
%                     % Left overlap with Y
%                     Y_crop = zeros(tilesize,tilesize,Nvar);
%                     Y_crop(1:tilesize, 1:2*overlap, :) = Y11;
%                     Y_temp(1:tilesize, 1:tilesize, :) = filtered_write(Y_crop(1:tilesize, 1:tilesize,...
%                         :), X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, :, sub(3)), M);
                elseif ( pos(1) > 1 && pos(2) == 1 )
                    % Top overlap with Y
                    Y_temp(1:tilesize, 1:tilesize, :) = filtered_write(Y1(1:tilesize, 1:tilesize,...
                        :), X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, :, sub(3)), M);
                else % ( pos(1) > 1 && pos(2) > 1 )
                    if (s == 3)  || (pos(s-1) == 1 && s > 3)
                        Y_crop = Y1(1:tilesize, 1:tilesize, :);
                    else
                        Y_crop = zeros(tilesize,tilesize,Nvar);
                        Y_crop(1:tilesize, 1:2*overlap, :) = Y11;
                    end
                    Y_temp(1:tilesize, 1:tilesize, :) = filtered_write(Y_crop(1:tilesize, 1:tilesize, :), ...
                        X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, :, sub(3)), M);
                end
            end
        end
        Y_temp_store{1,s} = Y_temp;
        
%         %         if (pos(1)>=7 && pos(2)>=4)
%         figure(101);hold on
%         colormap gray
% %         subplot(1,2,2)
%         imagesc(flipud(Y_temp)); %pos
%         nx=size(Y_temp,1);
%         scatter(D{1,1}(:,2)-startJ(1)+1 , nx-(D{1,1}(:,1)-startI(1)) , ones(size(D{1,1}(:,1)))*20 , D{1,1}(:,3));
%         plot(D{1,1}(:,2)-startJ(1)+1 , nx-(D{1,1}(:,1)-startI(1)) ,'r.');
% %         scatter3(D{1,1}(:,2)-startJ(1)+1,D{1,1}(:,1)-startI(1)+1,ones(size(D{1,1}(:,1)))+200,100,D{1,1}(:,3),'filled','MarkerEdgeColor','r')
% %         axis([0 32 0 32])
%         axis equal tight
%         'stop';
%         
%         set(gca,'NextPlot','replacechildren');
%         currFrame_1 = getframe;
%         writeVideo(vidObj_1,currFrame_1);
        %         end
        
    end
end
s = s - 1;
% close(vidObj_1);

if (tilesize == 2 )
    for ll=1:size(D{1,1},1)
        Y_temp(D{1,1}(ll,1)-startI(end-1)+1,D{1,1}(ll,2)-startJ(end-1)+1)=D{1,1}(ll,3);
    end
end

function A = filtered_write(A, B, M)
Nvar=size(A,3);
for i = 1:Nvar
    A(:, :, i) = A(:,:,i) .* (M == 0) + B(:,:,i) .* (M == 1);
end;
