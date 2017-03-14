function Y = imagequilt_conditional_3D_v13(X, var_type, m, n, o, D, w, tilesize, overlap, nbreplicates, w_v, temp_split ,do_cut)
% Performs the Efros/Freeman Unconditional Image quilting algorithm on the input TI

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
% Nvar = size(X,3);
% Nlayer = size(X,4);
x_cord=D(:,1); % x coordinates of conditioning points
y_cord=D(:,2); % y coordinates of conditioning points
z_cord=D(:,3); % y coordinates of conditioning points
P = size (D,1); % Number of Conditioning points

if( overlap(1) >= tilesize(1) || overlap(2) >= tilesize(2) || overlap(3) >= tilesize(3))
    error('Overlap must be less than tilesize');
end;

destsize_x = m * tilesize(1) - (m-1) * overlap(1);
destsize_y = n * tilesize(2) - (n-1) * overlap(2);
destsize_z = o * tilesize(3) - (o-1) * overlap(3);
Y = zeros(destsize_x, destsize_y, destsize_z);   % Y = the size of output
% Y = zeros(destsize_x, destsize_y, destsize_z, Nvar);

%% Determining starting and ending points of each tile in both i, j direction
for i=1:m,
    for j=1:n,
        for k=1:o,
            pos = [i,j,k]
            Total = [m,n,o]
            if (pos(1) == 1 && pos(2) >= 3 && pos(3) >= 9)
                'stop';
            end
            startI = (i-1)*tilesize(1) - (i-1) * overlap(1) + 1; % starting position of each tile in x direction
            startJ = (j-1)*tilesize(2) - (j-1) * overlap(2) + 1; % end position of each tile in y direction
            startK = (k-1)*tilesize(3) - (k-1) * overlap(3) + 1; % end position of each tile in z direction
            endI = startI + tilesize(1) -1; % starting position of each tile in x direction
            endJ = startJ + tilesize(2) -1; % end position of each tile in y direction
            endK = startK + tilesize(3) -1; % end position of each tile in z direction
            
            % Determining the distances from each tile to the overlap region
            distances1 = zeros( size(X,1)-tilesize(1)+1, size(X,2)-tilesize(2)+1, size(X,3)-tilesize(3)+1); % initialize the distances1 matrix
            distances2 = zeros( size(X,1)-tilesize(1)+1, size(X,2)-tilesize(2)+1, size(X,3)-tilesize(3)+1); % initialize the distances1 matrix
            
            %%
            % Compute the sum of squared distances in z direction between X and Y(startI:endI, startJ:startJ+overlap-1, 1:3)
            % for each possible overlap of Y(startI:endI, startJ:startJ+overlap-1, 1:3) on X.
            if( k > 1 )
                z = zeros(size(X,1)-tilesize(1), size(X,2)-tilesize(2), size(X,3)-overlap(3), size(X,4));
                %                 Z = zeros(size(X,1)-tilesize+1, size(X,2)-overlap+1, size(X,3), size(X,4));
                distances1 = ssd_v4_3D( X, Y(startI:endI, startJ:endJ, startK:startK+(overlap(3)-1)));
                distances1 = distances1(ceil((size(X,1)-size(z,1))/2):end-floor((size(X,1)-size(z,1))/2), ceil((size(X,2)-size(z,2))/2):end-floor((size(X,2)-size(z,2))/2)...
                    , ceil((size(X,3)-size(z,3))/2):end-floor((size(X,3)-size(z,3))/2));
                distances1 = distances1(1:end, 1:end, 1:end-tilesize(3)+overlap(3));   % considering only the overlapping region
                distances2 = conditioning_distance_3D(X,var_type,P,D,x_cord,y_cord,z_cord,startI,endI,startJ,endJ,startK,endK,tilesize,distances1);
            end
            
            %%
            % Compute the sum of squared distances in j direction between X and Y(startI:endI, startJ:startJ+overlap-1, 1:3)
            % for each possible overlap of Y(startI:endI, startJ:startJ+overlap-1, 1:3) on X.
            if( j > 1 )
                z = zeros(size(X,1)-tilesize(1), size(X,2)-overlap(2), size(X,3)-tilesize(3), size(X,4));
                %                 Z = zeros(size(X,1)-tilesize+1, size(X,2)-overlap+1, size(X,3), size(X,4));
                Z = ssd_v4_3D( X, Y(startI:endI, startJ:startJ+(overlap(2)-1), startK:endK));
                %                 Z = Z(1:end, 1:end-tilesize+overlap, 1:end);   % considering only the overlapping region
                Z = Z(ceil((size(X,1)-size(z,1))/2):end-floor((size(X,1)-size(z,1))/2), ceil((size(X,2)-size(z,2))/2):end-floor((size(X,2)-size(z,2))/2)...
                    , ceil((size(X,3)-size(z,3))/2):end-floor((size(X,3)-size(z,3))/2));
                Z = Z(1:end, 1:end-tilesize(2)+overlap(2), 1:end);
                %                if( k > 1 )
                distances1 = distances1 + Z;  % considering the overlapping regions in both i & j directions
                distances2 = conditioning_distance_3D(X,var_type,P,D,x_cord,y_cord,z_cord,startI,endI,startJ,endJ,startK,endK,tilesize,distances1);
                %             else
                %                    distances1 = Z; % Compute the distances of overlap for i > 1 & j = 1
                %                end
            end
            
            %%
            % Compute the sum of squared distances in i direction between X and Y(startI:startI+overlap-1, startJ:endJ, 1:3)
            % for each possible overlap of Y(startI:startI+overlap-1, startJ:endJ, 1:3) on X.
            if( i > 1 )
                z = zeros(size(X,1)-overlap(1), size(X,2)-tilesize(2), size(X,3)-tilesize(3), size(X,4));
                %                 Z = zeros(size(X,1)-overlap+1, size(X,2)-tilesize+1, size(X,3), size(X,4));
                Z = ssd_v4_3D( X, Y(startI:startI+(overlap(1)-1), startJ:endJ, startK:endK));
                %                 Z = Z(1:end-tilesize+overlap, 1:end, 1:end);   % considering only the overlapping region
                Z = Z(ceil((size(X,1)-size(z,1))/2):end-floor((size(X,1)-size(z,1))/2), ceil((size(X,2)-size(z,2))/2):end-floor((size(X,2)-size(z,2))/2)...
                    , ceil((size(X,3)-size(z,3))/2):end-floor((size(X,3)-size(z,3))/2));
                Z = Z(1:end-tilesize(1)+overlap(1), 1:end, 1:end);
 %               if( j > 1 )
                    distances1 = distances1 + Z;  % considering the overlapping regions in both i & j directions
                    distances2 = conditioning_distance_3D(X,var_type,P,D,x_cord,y_cord,z_cord,startI,endI,startJ,endJ,startK,endK,tilesize,distances1);
 %               else
 %                   distances1 = Z; % Compute the distances of overlap for i > 1 & j = 1
 %               end
            end
            
             %%
            % If both i > 1 & j > 1, compute the distances of the overlap
            if( i > 1 && j > 1 )
                z = zeros(size(X,1)-overlap(1), size(X,2)-overlap(2), size(X,3)-tilesize(3), size(X,4));
                %                 Z = zeros(size(X,1)-overlap+1, size(X,2)-overlap+1, size(X,3), size(X,4));
                Z = ssd_v4_3D( X, Y(startI:startI+(overlap(1)-1), startJ:startJ+(overlap(2)-1), startK:endK));
                %                 Z = Z(1:end-tilesize+overlap, 1:end-tilesize+overlap, 1:end);
                Z = Z(ceil((size(X,1)-size(z,1))/2):end-floor((size(X,1)-size(z,1))/2), ceil((size(X,2)-size(z,2))/2):end-floor((size(X,2)-size(z,2))/2)...
                    , ceil((size(X,3)-size(z,3))/2):end-floor((size(X,3)-size(z,3))/2));
                Z = Z(1:end-tilesize(1)+overlap(1), 1:end-tilesize(2)+overlap(2), 1:end);
                distances1 = distances1 - Z;
                distances2 = conditioning_distance_3D(X,var_type,P,D,x_cord,y_cord,z_cord,startI,endI,startJ,endJ,startK,endK,tilesize,distances1);
            end
            
            % If both i > 1 & j > 1, compute the distances of the overlap
            if( j > 1 && k > 1 )
                z = zeros(size(X,1)-tilesize(1), size(X,2)-overlap(2), size(X,3)-overlap(3), size(X,4));
                %                 Z = zeros(size(X,1)-overlap+1, size(X,2)-overlap+1, size(X,3), size(X,4));
                Z = ssd_v4_3D( X, Y(startI:endI, startJ:startJ+(overlap(2)-1), startK:startK+(overlap(3)-1)));
                %                 Z = Z(1:end, 1:end-tilesize+overlap, 1:end-tilesize+overlap);
                Z = Z(ceil((size(X,1)-size(z,1))/2):end-floor((size(X,1)-size(z,1))/2), ceil((size(X,2)-size(z,2))/2):end-floor((size(X,2)-size(z,2))/2)...
                    , ceil((size(X,3)-size(z,3))/2):end-floor((size(X,3)-size(z,3))/2));
                Z = Z(1:end, 1:end-tilesize(2)+overlap(2), 1:end-tilesize(3)+overlap(3));
                distances1 = distances1 - Z;
                distances2 = conditioning_distance_3D(X,var_type,P,D,x_cord,y_cord,z_cord,startI,endI,startJ,endJ,startK,endK,tilesize,distances1);
            end
            
            % If both i > 1 & j > 1, compute the distances of the overlap
            if( k > 1 && i > 1 )
                z = zeros(size(X,1)-overlap(1), size(X,2)-tilesize(2), size(X,3)-overlap(3), size(X,4));
                %                 Z = zeros(size(X,1)-overlap+1, size(X,2)-overlap+1, size(X,3), size(X,4));
                Z = ssd_v4_3D( X, Y(startI:startI+(overlap(1)-1), startJ:endJ, startK:startK+(overlap(3)-1)));
                %                 Z = Z(1:end-tilesize+overlap, 1:end, 1:end-tilesize+overlap);
                Z = Z(ceil((size(X,1)-size(z,1))/2):end-floor((size(X,1)-size(z,1))/2), ceil((size(X,2)-size(z,2))/2):end-floor((size(X,2)-size(z,2))/2)...
                    , ceil((size(X,3)-size(z,3))/2):end-floor((size(X,3)-size(z,3))/2));
                Z = Z(1:end-tilesize(1)+overlap(1), 1:end, 1:end-tilesize(3)+overlap(3));
                distances1 = distances1 - Z;
                distances2 = conditioning_distance_3D(X,var_type,P,D,x_cord,y_cord,z_cord,startI,endI,startJ,endJ,startK,endK,tilesize,distances1);
            end
            
            % If both i > 1 & j > 1, compute the distances of the overlap
            if( i > 1 && j > 1 && k > 1 )
                z = zeros(size(X,1)-overlap(1), size(X,2)-overlap(2), size(X,3)-overlap(3), size(X,4));
                %                 Z = zeros(size(X,1)-overlap+1, size(X,2)-overlap+1, size(X,3), size(X,4));
                Z = ssd_v4_3D( X, Y(startI:startI+(overlap(1)-1), startJ:startJ+(overlap(2)-1), startK:startK+(overlap(3)-1)));
                %                 Z = Z(1:end-tilesize+overlap, 1:end-tilesize+overlap, 1:end-tilesize+overlap);
                Z = Z(ceil((size(X,1)-size(z,1))/2):end-floor((size(X,1)-size(z,1))/2), ceil((size(X,2)-size(z,2))/2):end-floor((size(X,2)-size(z,2))/2)...
                    , ceil((size(X,3)-size(z,3))/2):end-floor((size(X,3)-size(z,3))/2));
                Z = Z(1:end-tilesize(1)+overlap(1), 1:end-tilesize(2)+overlap(2), 1:end-tilesize(3)+overlap(3));
                distances1 = distances1 + Z;
                distances2 = conditioning_distance_3D(X,var_type,P,D,x_cord,y_cord,z_cord,startI,endI,startJ,endJ,startK,endK,tilesize,distances1);
            end
            
            %%
            if( i == 1 && j == 1 && k == 1 && sum(distances2(:)) == 0)
                sub(1) = ceil(rand*(size(X,1)-tilesize(1)));
                sub(2) = ceil(rand*(size(X,2)-tilesize(2)));
                sub(3) = ceil(rand*(size(X,3)-tilesize(3)));
                %                 sub(3) = round(rand*size(X,4));
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
                distances = (1-w)*distances1 + w*distances2;
                [~,distances_index]=sort(distances(:),'ascend');
                % draw one in nbreplicates of the best distances
                idx = distances_index(ceil(rand(1)*nbreplicates));
                [sub(1), sub(2), sub(3)] = ind2sub(size(distances), idx);
            end
            
            
            %%
            if do_cut==1
                % Initialize the mask to all ones
                M = ones(tilesize(1), tilesize(2), tilesize(3));
                
%                 allC=zeros(tilesize, tilesize, tilesize);
                % If we have a end overlap in Z direction
                if( k > 1 )
                    %Compute the SSD in the border region
                    Err = ( X(sub(1):sub(1)+tilesize(1)-1, sub(2):sub(2)+tilesize(2)-1, sub(3):sub(3)+overlap(3)-1, :) - Y(startI:endI, startJ:endJ, startK:startK+overlap(3)-1, :) ).^2;
%                     Err = ( X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+overlap-1, sub(3):sub(3)+tilesize-1, :) - Y(startI:endI, startJ:startJ+overlap-1, startK:endK, :) ).^2;
%                     Err = ( X(sub(1):sub(1)+overlap-1, sub(2):sub(2)+tilesize-1, sub(3):sub(3)+tilesize-1, :) - Y(startI:startI+overlap-1, startJ:endJ, startK:endK, :) ).^2;
                  
                    %                 % Calculate combined error for all the variables considering the weightage
                    %                 E = zeros(size(Err));
                    %                 for k=1:size(X,4), % No of variables
                    %                     E(:,:,:,k) = w_v(k)*Err(:,:,:,k);
                    %                 end
                    %                 E = sum(E,3);
                    
                    %Compute the mincut array
                    %                 C = mincut_corrected_kashif(Err.^2, 0);
                    C = DP_3D_surface_v2_function_Kashif_v2(Err,3,0);
%                     allC(C(1:end, 1:end, 1:overlap)==0)=1;
                    %Compute the mask and write to the destination
                    M(1:end, 1:end, 1:overlap(3)) = double(C >= 0);
%                     M(1:end, 1:overlap, 1:end) = double(C >= 0);
%                     M(1:overlap, 1:end, 1:end) = double(C >= 0);                    
                    
%                     figure(31);clf;
%                     ViewGrid(C);
%                     view(-50,-50)
%                     figure(32);clf;
%                     ViewGrid(M);
%                     view(-50,-50)
                end
                
                % If we have a left overlap in Y direction
                if( j > 1 )
                    %Compute the SSD in the border region
                    Err = ( X(sub(1):sub(1)+tilesize(1)-1, sub(2):sub(2)+overlap(2)-1, sub(3):sub(3)+tilesize(3)-1, :) - Y(startI:endI, startJ:startJ+overlap(2)-1, startK:endK, :) ).^2;
%                     Err = ( X(sub(1):sub(1)+overlap-1, sub(2):sub(2)+tilesize-1, sub(3):sub(3)+tilesize-1, :) - Y(startI:startI+overlap-1, startJ:endJ, startK:endK, :) ).^2;
                    
                    %                 % Calculate combined error for all the variables considering the weightage
                    %                 E = zeros(size(Err));
                    %                 for k=1:size(X,3), % No of variables
                    %                     E(:,:,k) = w_v(k)*Err(:,:,k); %    ???????HERE is this okay???????
                    %                 end
                    %                 E = sum(E,3); %    ???????HERE is this okay???????
                    
                    %Compute the mincut array
                    %                 C = mincut_corrected_kashif(Err.^2, 0);
                    C = DP_3D_surface_v2_function_Kashif_v2(Err,2,0);
%                     allC(C(1:end, 1:overlap, 1:end)==0)=1;
                    %Compute the mask and write to the destination
%                    if ( k == 1)
                        %                     M(1:overlap, 1:end, 1:end) = M(1:overlap, 1:end, 1:end) .* double(C >= 0);
%                        M(1:end, 1:overlap, 1:end) = double(C >= 0);
                    %elseif ( k > 1)
                        M(1:end, 1:overlap(2), 1:end) = M(1:end, 1:overlap(2), 1:end) .* double(C >= 0);
                    %end
                    
%                     figure(31);clf;
%                     ViewGrid(C);
%                     view(-50,-50)
%                     figure(32);clf;
%                     ViewGrid(M);
%                     view(-50,-50)
                end
                
                % We have a top overlap
                if( i > 1 )
                    %Compute the SSD in the border region
%                     Err = ( X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+overlap-1, sub(3):sub(3)+tilesize-1, :) - Y(startI:endI, startJ:startJ+overlap-1, startK:endK, :) ).^2;
                  Err = ( X(sub(1):sub(1)+overlap(1)-1, sub(2):sub(2)+tilesize(2)-1, sub(3):sub(3)+tilesize(3)-1, :) - Y(startI:startI+overlap(1)-1, startJ:endJ, startK:endK, :) ).^2;
                    
                    %                 % Calculate combined error for all the variables considering the weightage
                    %                 E = zeros(size(Err));
                    %                 for k=1:size(X,3), % No of variables
                    %                     E(:,:,k) = w_v(k)*Err(:,:,k); 
                    %                 end
                    %                 E = sum(E,3); 
                    
                    %Compute the mincut array
                    C = DP_3D_surface_v2_function_Kashif_v2(Err,1,0);
%                     if k > 1
%                         C = DP_3D_surface_v2_function_Kashif(Err,2,0);
%                     elseif j >= 1
%                         C = DP_3D_surface_v2_function_Kashif(Err,3,0);
%                     end
                    %                     allC(C(1:overlap, 1:end, 1:end)==0)=1;
                    %Compute the mask and write to the destination
 %                   if ( k == 1 && j == 1)
                        %                         M(1:end, 1:overlap, 1:end) = M(1:end, 1:overlap, 1:end) .* double(C >= 0);
 %                       M(1:overlap, 1:end, 1:end) = double(C >= 0);
 %                   elseif ( k > 1 || j > 1) %                  ??????? HERE is this okay ???????
                        M(1:overlap(1), 1:end, 1:end) = M(1:overlap(1), 1:end, 1:end) .* double(C >= 0);
 %                   end
                    
%                     figure(31);clf;
%                     ViewGrid(C);
%                     view(-50,-50)
%                     figure(32);clf;
%                     ViewGrid(M);
%                     view(-50,-50)
                end
            end
            
            if( i == 1 && j == 1 && k == 1) % For the first block
                Y(startI:endI, startJ:endJ, startK:endK) = X(sub(1):sub(1)+tilesize(1)-1, sub(2):sub(2)+tilesize(2)-1, sub(3):sub(3)+tilesize(3)-1);
%                 allC=M;
            end
            
            %% Write to the destination using the mask
            %             if( k > 1 )
            %                 for l=startK:1:endK
            %                     Y(startI:endI, startJ:endJ, l) = X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, sub(3)+l-startK);
            %                 end
            %             end
            %
            %             if( j > 1 )
            %                 for l=startJ:1:endJ
            %                     Y(startI:endI, l, startK:endK) = X(sub(1):sub(1)+tilesize-1, sub(2)+l-startJ, sub(3):sub(3)+tilesize-1);
            %                 end
            %             end
            %
            %             if( i > 1 )
            %                 for l=startI:1:endI
            %                     Y(l, startJ:endJ, startK:endK) = X(sub(1)+l-startI, sub(2):sub(2)+tilesize-1, sub(3):sub(3)+tilesize-1);
            %                 end
            %             end
            
            
%             Cutsmap(startI:endI, startJ:endJ, startK:1:endK)=allC;
            
            %% Write to the destination using the mask
            if do_cut==1
                if( k > 1 )
                    for l=startK:1:endK
                        Y(startI:endI, startJ:endJ, l) = filtered_write(Y(startI:endI, startJ:endJ, l), ...
                            X(sub(1):sub(1)+tilesize(1)-1, sub(2):sub(2)+tilesize(2)-1, sub(3)+l-startK), M(:,:,l-startK+1));
                    end
                end
                
                if( j > 1 )
                    for l=startJ:1:endJ
                        Y(startI:endI, l, startK:endK) = filtered_write(Y(startI:endI, l, startK:endK), ...
                            X(sub(1):sub(1)+tilesize(1)-1, sub(2)+l-startJ, sub(3):sub(3)+tilesize(3)-1), M(:,l-startJ+1,:));
                    end
                end
                
                if( i > 1 )
                    for l=startI:1:endI
                        Y(l, startJ:endJ, startK:endK) = filtered_write(Y(l, startJ:endJ, startK:endK), ...
                            X(sub(1)+l-startI, sub(2):sub(2)+tilesize(2)-1, sub(3):sub(3)+tilesize(3)-1), M(l-startI+1,:,:));
                    end
                end
                
            else
                Y(startI:endI, startJ:endJ, startK:endK) = X(sub(1):sub(1)+tilesize(1)-1, sub(2):sub(2)+tilesize(2)-1, sub(3):sub(3)+tilesize(3)-1);
            end
            
%                         %%
%                         figure(100);clf
%                         ViewGrid(distances1)
%                         colorbar
%                         figure(101);clf;hold on
%                         subplot(2,2,1)
%                         ViewGrid(Y(startI:endI, startJ:endJ, startK:endK))
%                         view(-50,-50)
%                         subplot(2,2,2)
%                         ViewGrid(X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, sub(3):sub(3)+tilesize-1))
%                         view(-50,-50)
%                         %             subplot(2,2,3)
%                         %             ViewGrid(M)
%                         subplot(2,2,4)
%                         ViewGrid(Y)
%                         view(-50,-50)
        end
    end
end

function A = filtered_write(A, B, M)
% Nvar=size(A,4);
% for i = 1:Nvar
A(:,:,:) = A(:,:,:) .* (M == 0) + B(:,:,:) .* (M == 1);
% end;
