clear; home
clc;

% This code can be used for both unconditional and conditional image 
% quilting simulations with Template Spliting (TS)

% Authors: Kashif Mahmud and Gregoire Mariethoz, 2013

% List of Inputs:
% 1. X1:  The source image to be used in synthesis (TI)
% 2. p and q: The required size of output realizations in X and Y dimensions
% 3. tilesize: The dimensions of each square tile
% 4. overlap: The amount of overlap to allow between pixels (def: 1/6 tilesize)
% 4. nbreplicates: Number of replicas considered to find the best possible patch
% 5. c: No of conditioning points (randomly taken from TI) or D: Conditioning hard data
% 6. w: Weightage of conditioning
% 7. w_v: VARIABLES' WEIGHTS

%% 1st INPUT: Read the input training image (TI)
X1 = LoadGrid('ti_lena_crop.SGEMS');

%%%%%%%%%%% BELOW ADJUST THE PARAMETERS %%%%%%%%%%%
tilesize = 30; % Optimal size for the patch size is between 1/7 and 1/8 of the TI
% (To get better realizations tilesize should be minimum)
tilesize_vary = 0; % Randomize patch size (max +-10% of patch size)
overlap = 7; % Optimal size for the overlap region is between 1/3 and 1/4 of the tilesize
nbreplicates = 3; % Standard Number of replicas = 5 to 10 to avoid "Verbatim Copy"
cond = 1; % Define whether Conditional or Unconditional simulation (0 or 1)
c = 50; % No of conditioning points (Don’t need to be defined if you load your Data file in the code)
w = 0.99; % Conditioning weight (Range [0 1])
w_v = 1; % Each variable's weight (For Multivariable IQ analysis)
% This is univariate IQ version
nbrealz = 1; % Number of realizations to be produced
temp_split = 1; % temp_split = 1 Perform template splitting; temp_split = 0 Do not perform template splitting
do_cut = 1; % do_cut = 1 Perform optimum cutting; do_cut = 0 Do not perform optimum cutting

show_TI = 1; % Showing the TI
show_realizations = 1; % Showing all the simulated realizations
show_etype_stdev_maps = 1; % Showing the E-type and Standard deviation maps for all variables
show_histogram = 1; % Showing the histograms
show_variogram = 1; % Showing the variograms
show_data = 1; % Showing the hard data
show_connectivity = 0; % Showing the connectivity in both X and Y directions for multi-facies TI & Realizations


%% Define the size of final realizations
% q = 256; % Size of the final realizations in X dimension
% p = 256; % Size of the final realizations in Y dimension

% Otherwise realization size is similar to TI:
[p q Nvar] = size (X1);

%% Calculate m, n from tilesize and overlap
% m: The number of tiles to be placed in the output image, in X dimension
% n: The number of tiles to be placed in the output image, in Y dimension
m = (p-overlap)/(tilesize-overlap);
m = ceil(m);
n = (q-overlap)/(tilesize-overlap);
n = ceil(n);

%% Consider the full TI
X = X1;

% % OR consider the TI having similar size of realizations
% X = X1(1:p,1:q);

%% Unif transform of TI %%%%%%%%%%%%%%
% UTbins=100;
% UTtail=0;
% [X,z,F]=UnifSTransform(X,UTbins,UTtail);
% 
% X1=(round(X1*10))/10;
% % X=(round(X*10))/10;
% % X=(round(X*100));
%% Check whether tile size is smaller or larger than the size of TI
if( tilesize > p || tilesize > q)
    error('Size of TI must be bigger than the Tile size');
end

%% Produce data set from TI with randomized locations or read given data set
if cond > 0
    %     fprintf('Produce data set from a unconditional realization with randomized locations:\n');
    %     X_data = imagequilt_Unconditional_v12(X, m, n, tilesize, overlap, 10, w_v,do_cut);
    %     X_data = X_data(1:p,1:q);
    %     Data = Produce_data(c, X_data);
    
    fprintf('Produce data set from TI with randomized locations:\n');
    Data = Produce_data(c, X); % Produce data set from TI with randomized locations
    D{1,1} = Data;
    
    %% Uniform transform data
    %%%%%%%% NOT NECESSARY SINCE THE DATA COME FROM A TRANSFORMED TI ALREADY,
    %%%%%%%% HOWEVER SHOULD BE DONE WITH REAL DATA
    %D{1,1}(:,3)=UnifSTransform_data(D{1,1}(:,3),UTtail,F,z)
    
end

%% Run the IQ function to get realizations

Y=zeros(p,q,nbrealz); % preallocate - 3D array to store all realizations in the 3rd dimension

if cond > 0
    fprintf('Required time for Conditioning Image Quilting:\n');
else
    fprintf('Required time for Unconditioning Image Quilting:\n');
end

I = 0;
% loop all realizations
for i=1:nbrealz
    fprintf(['Output realization #',num2str(i),': ']);
    tic
    if tilesize_vary > 0
        I = I + 1;
        if I <= tilesize_vary
            tilesize_new = tilesize + I;
        else
            I = I - tilesize_vary;
            tilesize_new = tilesize;
        end
    else
        tilesize_new = tilesize;
    end
    
    if cond > 0
        % generate realization at the size needed for the tiles arrangement, possibly too large
        Y1 = imagequilt_multivariate_v7(X, m, n, D, w, tilesize_new, overlap, nbreplicates, 1, temp_split);
    else
        Y1 = imagequilt_Unconditional_v12(X, m, n, tilesize_new, overlap, nbreplicates, w_v,do_cut);
    end
    Y(:,:,i) = Y1(1:p,1:q); % crop the realization to its true dimensions and store in 3D array Y
    toc
end

%% back-transforming TI, simuls and conditioning data  %%%%%%%
% X=UnifSTransform_inv(z,F,X);
% Y=UnifSTransform_inv(z,F,Y);
% D{1,1}(:,3)=UnifSTransform_inv(z,F,D{1,1}(:,3));
% 
% X=(round(X*10))/10;
% Y=(round(Y*10))/10;
% D{1,1}(:,3)=(round(D{1,1}(:,3)*10))/10;

%% Find out the simulated satisfied conditioning points in all realizations
if cond > 0
    x=D{1,1}(:,1); % x coordinates of conditioning points
    y=D{1,1}(:,2); % y coordinates of conditioning points
    Satisfied_conditioning_points = zeros(c,1,nbrealz);
    for i=1:nbrealz
        for k = 1:c
            if x(k) <= size (Y(:,:,i),1) && y(k) <= size (Y(:,:,i),2)
                Satisfied_conditioning_points(k,1,i) = Y(x(k),y(k),i);
            else
                fprintf('Given data at location (%d, %d) is outside of final realization.\n', x(k), y(k));
            end
        end
    end
end

%% Showing the input training image
originx=0;
originy=0;
stepx=1;
stepy=1;
if show_TI == 1
    nx_ti=size(X,1);
    ny_ti=size(X,2);
    [xx_ti,yy_ti]=meshgrid(originy:stepy:(ny_ti-1)*stepy , originx:stepx:(nx_ti-1)*stepx);
    
    figure(1);clf;colormap default;hold on
    surf(xx_ti,yy_ti,flipud(double(X)))
    view(0,90)
    shading flat
    axis equal tight
    colormap default
end

%% Showing the output realizations
originx=0;
originy=0;
stepx=1;
stepy=1;
nx=size(Y(:,:,i),1);
ny=size(Y(:,:,i),2);
[xx,yy]=meshgrid(originy:stepy:(ny-1)*stepy , originx:stepx:(nx-1)*stepx);
if show_realizations == 1
    if nbrealz > 5
        for i=1:5
            figure(1+i);clf;colormap gray;hold on
            surf(xx,yy,flipud(double(Y(:,:,i))))
            view(0,90)
            shading flat
            axis equal tight
            colormap default
        end
    else
        for i=1:nbrealz
            figure(1+i);clf;colormap gray;hold on
            surf(xx,yy,flipud(double(Y(:,:,i))))
            view(0,90)
            shading flat
            axis equal tight
            colormap default
            colorbar
        end
    end
end

%% Plotting the Satisfied conditioning points vs Given field data
if cond > 0
    figure(2+nbrealz);clf;hold on
    meanerror=zeros(1,nbrealz);
    stderror=zeros(1,nbrealz);
    for i=1:nbrealz
        plot(D{1,1}(:,3),Satisfied_conditioning_points(:,:,i),'.')
        meanerror(i)=mean(Satisfied_conditioning_points(:,:,i) - D{1,1}(:,3));
        stderror(i)=std(Satisfied_conditioning_points(:,:,i) - D{1,1}(:,3));
    end
    plot([0,max(D{1,1}(:,3))],[0,max(D{1,1}(:,3))],'r')
    disp(['Mean conditioning error over all realizations: ',num2str(mean(meanerror))])
    disp(['Mean std error over all realizations: ',num2str(mean(stderror))])
    axis equal square
    axis tight
    xlabel('Conditioning Data (Grayscale)')
    ylabel('Simulated Data (Grayscale)')
end

%% Showing the E-type and standard error maps
if show_etype_stdev_maps == 1
    meanrealz=mean(Y,3);
    stdrealz=std(Y,[],3);
    figure(3+nbrealz);clf;hold on
    surf(xx,yy,flipud(double(meanrealz)));
    if cond > 0
        scatter3(D{1,1}(:,2),nx-D{1,1}(:,1),ones(size(D{1,1}(:,1)))+200,20,D{1,1}(:,3),'o','MarkerEdgeColor','r')
    end
    view(0,90)
    shading flat
    axis equal tight
    colormap default
    colorbar
    title('E-type map');
    
    figure(4+nbrealz);clf;hold on
    surf(xx,yy,flipud(double(stdrealz)));
    if cond > 0
        plot3(D{1,1}(:,2),nx-D{1,1}(:,1),ones(size(D{1,1}(:,1)))+200,'or');
    end
    view(0,90)
    shading flat
    axis equal tight
    colormap default
    colorbar
    title('Standard daviation map');
end

%% Showing the Variogram of TI and output realizations
if show_variogram==1
    % Histogram transformation to match the TI
    %         for i=1:nbrealz
    %             Y1=Y(:,:,i);
    %             Y(:,:,i)=( ( (Y1 - mean(Y1(:))) / std(Y1(:)) ) * std(X(:))) + mean(X(:));
    %         end
    
    % Variogram of TI and output realizations
    tic
    figure(5+nbrealz);clf;hold on
    [xxti,yyti]=meshgrid(1:size(X,2),1:size(X,1));
    pts_ti=Grid2Pts(X,xxti,yyti,zeros(size(xxti)));
    sampling_ti = round(size(pts_ti,1)/1000);
    [meanhti,gammahti]=Variogram(pts_ti,20,60,sampling_ti,0);
    pts = zeros(p*q,1+c,nbrealz);
    for i=1:nbrealz
        pts=Grid2Pts(Y(:,:,i),xx,yy,zeros(size(xx)));
        sampling = round(size(pts,1)/1000);
        [meanh,gammah]=Variogram(pts,20,60,sampling,0);
        plot(meanh,gammah,'k--')
    end
    plot(meanhti,gammahti,'r','LineWidth',2)
    
    axis equal square
    axis tight
    xlabel('Distance')
    ylabel('Variance')
    fprintf('To draw the Variogram: ');
    toc
end

%% Showing the Histogram of TI and output realizations
if show_histogram==1
    figure(6+nbrealz);clf;hold on
    [Nti,Xti]=hist(X);
    Nti=Nti/sum(Nti);
    
    for i=1:nbrealz
        Y1=Y(:,:,i);
        [N1,X1]=hist(Y1(:));
        N1=N1/sum(N1);
        [X1,N1] = smoothLine(X1,N1,10);
        plot(X1,N1,'k--');
    end
    
    [Xti,Nti] = smoothLine(Xti,Nti,10);
    plot(Xti,Nti,'r','LineWidth',2);
    
    xlabel('log K')
    ylabel('Probability')
    axis equal square
    axis tight
end

%% Plotting the locations of field Data
if cond > 0
    if show_data==1
        figure(7+nbrealz);clf;hold on
        scatter(D{1,1}(:,2),nx-D{1,1}(:,1),30,D{1,1}(:,3),'o','filled','MarkerEdgeColor','r')
        axis equal tight
        colormap default
        colorbar
    end
end

%% Showing the connectivity in both X and Y directions
if show_connectivity==1
    figure(8+nbrealz);clf;hold on
    Ctr = ConnectFct(X,1,0);
    for i=1:nbrealz
        C1re = ConnectFct(Y(:,:,i),1,0);
        plot(C1re.x(:,1),'r--');
    end
    plot(Ctr.x(:,1),'k','LineWidth',2);
    axis equal square
    axis tight
    
    figure(9+nbrealz);clf;hold on
    for i=1:nbrealz
        C1re = ConnectFct(Y(:,:,i),1,0);
        plot(C1re.x(:,2),'g--');
    end
    plot(Ctr.x(:,2),'k','LineWidth',2);
    axis equal square
    axis tight
    
    
    figure(10+nbrealz);clf;hold on
    for i=1:nbrealz
        C1re = ConnectFct(Y(:,:,i),1,0);
        plot(C1re.y(:,1),'r--');
    end
    plot(Ctr.y(:,1),'k','LineWidth',2);
    axis equal square
    axis tight
    
    figure(11+nbrealz);clf;hold on
    for i=1:nbrealz
        C1re = ConnectFct(Y(:,:,i),1,0);
        plot(C1re.y(:,2),'g--');
    end
    plot(Ctr.y(:,2),'k','LineWidth',2);
    axis equal square
    axis tight
end