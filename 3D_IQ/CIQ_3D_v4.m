clear; home

% RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));
% profile on

% List of Inputs:
% 1. X1:  The source image to be used in synthesis (TI)
% 2. p and q: The required size of output realizations in X and Y dimensions
% 3. tilesize: The dimensions of each square tile
% 4. overlap: The amount of overlap to allow between pixels (def: 1/6 tilesize)
% 4. nbreplicates: Number of replicas considered to find the best possible patch
% 5. c: No of conditioning points (randomly taken from TI) or D: Conditioning hard data
% 6. w: Weightage of conditioning

%% Read the input training image (TI)
% X = LoadGrid('TI1_WCA.SGEMS');
% X = X(1:50, 1:40, 1:30);

% X = LoadGrid('TI3d-c1.SGEMS');

X = LoadGrid('paola_bloc.SGEMS');
X = X(51:150, 51:130, 1:45);

% X = LoadGrid('ti_ifp_stationary_part_6f_101x101x536.SGEMS');
% X = MGSimulFFT(200,200,60,5,5,2,100,50,30);
% X1 = load('3D_cont');
% X = X1.imin;
% X = MGSimulFFT(100,100,30,5,5,2,50,30,30);

%% Showing the input training image
X=(round(X*100))/100;
figure(1);clf;
ViewGrid(X)
colormap default
axis equal tight
title('TI');
view(30,30)

%% Define the size of final realizations
% p = 250; % Size of the final realizations in I,X dimension
% q = 80; % Size of the final realizations in J,Y dimension
% r = 30; % Size of the final realizations in K,Z dimension

% SG = LoadGrid('SG.SGEMS');
% [p q r ~] = size (SG);

% Otherwise realization size is similar to TI:
[p q r] = size (X);

tilesize = [30 20 20]; % Optimal size for the patch size is between 1/7 and 1/8 of the TI

% Tilesize can be assigned as size of TI/8 if not given as input
% tilesize(1) = round(p/7);
% tilesize(2) = round(q/7);
% tilesize(3) = round(r/7);

overlap = [7 5 5]; %Optimal size for the overlap region is between 1/3 and 1/4 of the patch size

% Overlap can be assigned as tilesize/4 in case of not given as input
% overlap(1) = round(tilesize(1)/4);
% overlap(2) = round(tilesize(2)/4);
% overlap(3) = round(tilesize(3)/4);

for i=1:3
    if overlap(i) < 2;
        overlap(i) = 2; % Assign minimum value of overlap as 2
    end
end

nbreplicates = 3; % Standard Number of replicas = 5 to avoid "Verbatim Copy"
% % Number of replicas can be assigned as 5 in case of not given as input
do_cut = 1; % do_cut = 1 Perform optimum cutting; do_cut = 0 Do not perform optimum cutting
c = 0; % No of conditioning points (Generate from a simulation with randomized location)
% If one has to read field data set then put any value of c greater than 0
w = 0; % Weightage of conditioning

show_TI = 1; % Showing the TI
show_realizations = 1; % Showing all the simulated realizations
show_histogram = 0; % Showing the histograms
show_variogram = 0; % Showing the variograms
show_data = 1; % Showing the hard data
show_connectivity = 0; % Showing the connectivity in both X and Y directions
w_v = 1; % Each variable's weight (For Multivariable IQ analysis)
% This is an univariate IQ version, Multivariate version will follow soon
temp_split = 0; % Template splitting scheme will follow soon
nbrealz = 1; % Number of realizations to be produced
var_type = 2; % catagorical = 1 or continuous = 2 

%% Calculate m, n from tilesize and overlap
% m: The number of tiles to be placed in the output image, in X dimension
% n: The number of tiles to be placed in the output image, in Y dimension
% o: The number of tiles to be placed in the output image, in Z dimension
m = (p-overlap(1))/(tilesize(1)-overlap(1));
m = ceil(m);
n = (q-overlap(2))/(tilesize(2)-overlap(2));
n = ceil(n);
o = (r-overlap(3))/(tilesize(3)-overlap(3));
o = ceil(o);

%% Consider the full TI
% X = X1;

% % OR consider the TI having similar size of realizations
% X = X1(1:p,1:q);

%% Check whether tile size is smaller or larger than the size of TI
if( tilesize(1) > p || tilesize(2) > q || tilesize(3) > r)
    error('Size of TI must be bigger than the Tile size');
end

%% Produce data set from TI with randomized locations or read given data set
if c > 0
%         X_data = imagequilt_Unconditional_3D_v13(X, m, n, o, tilesize, overlap, 3, do_cut);
%         X_data = X_data(1:p,1:q,1:r);
    %     WriteGrid(X_data,'IQ_3D_realzation_WCA_10.SGEMS', 'SGEMS');
%         X_data = LoadGrid('IQ_3D_realzation_WCA_10.SGEMS');
        X_data = LoadGrid('IQ_3D_realzation_Paola_50.SGEMS');
        D = SampleGrid(X_data,c);
    %     D = SampleGrid(X,c);
    %     D(:,4)=(round(D(:,4)*100))/100;
    
%     Data = LoadPts('Well_data.SGEMS');
%     D = Data;
%     j = 1;
%     for i = 1:size(Data,1)
%         if Data(i,3) < 35 || Data(i,3) > (35 + r)
%             D(j,:) = [];
%         else
%             j = j + 1;
%         end
%     end

    % D(:,3) = D(:,3) + 1;
    c = size(D,1);
end

%% Plotting the locations of field Data
if c > 0
    if show_data==1
        figure(7+nbrealz);clf;colormap gray;hold on
%         scatter(D{1,1}(:,1),D{1,1}(:,2),30,D{1,1}(:,3),'o','filled')
%         nx=size(Y,1);
        scatter3(D(:,2),p-D(:,1),D(:,3),30,D(:,4),'o','filled','MarkerEdgeColor','r')
%         scatter(D{1,1}(:,2),nx-D{1,1}(:,1),ones(size(D{1,1}(:,1)))+200,40,D{1,1}(:,3),'filled','MarkerEdgeColor','r')
%         view(0,90)
%         shading flat
        axis equal tight
        colormap default
        colorbar 
        view(30,30)
    end
end

%% Run the IQ function to get realizations

Y=zeros(p,q,r,nbrealz); % preallocate - 3D array to store all realizations in the 3rd dimension

if c > 0
    fprintf('Required time for Conditioning Image Quilting:\n');
else
    fprintf('Required time for Unconditioning Image Quilting:\n');
end

% loop all realizations
for i=1:nbrealz
    fprintf(['Output realization #',num2str(i),': ']);
    tic
    if c > 0
        % generate realization at the size needed for the tiles arrangement, possibly too large
        Y1 = imagequilt_conditional_3D_v13(X, var_type, m, n, o, D, w, tilesize, overlap, nbreplicates, w_v, temp_split ,do_cut);
    else
        Y1 = imagequilt_Unconditional_3D_v13(X, m, n, o, tilesize, overlap, nbreplicates,do_cut);
    end
    Y(:,:,:,i) = Y1(1:p,end-q+1:end,1:r);
    toc
end

% WriteGrid(X,'TI.SGEMS', 'SGEMS');
% WriteGrid(Y(:,:,:,1),'Realzation_1.SGEMS', 'SGEMS');
% WriteGrid(Y(:,:,:,2),'Realzation_2.SGEMS', 'SGEMS');
% WriteGrid(Y(:,:,:,3),'Realzation_3.SGEMS', 'SGEMS');

% %% MDS plot with ANODI
% DistMtrx = calculate3DModelVar_MPH(Y,X,10);
% 
% show3DModelVar2(DisMtrx,X,Y);

%% Find out the simulated satisfied conditioning points in all realizations
if c > 0
    x=D(:,1); % x coordinates of conditioning points
    y=D(:,2); % y coordinates of conditioning points
    z=D(:,3); % y coordinates of conditioning points
    Satisfied_conditioning_points = zeros(c,1,nbrealz);
    for i=1:nbrealz
        for k = 1:c
            if x(k) <= size (Y(:,:,:,i),1) && y(k) <= size (Y(:,:,:,i),2) && z(k) <= size (Y(:,:,:,i),3)
                Satisfied_conditioning_points(k,1,i) = Y(x(k),y(k),z(k),i);
            else
                fprintf('Given data at location (%d, %d) is outside of final realization.\n', x(k), y(k), z(k));
            end
        end
    end
end

%% Showing the input training image
cmin = min(X(:));
cmax = max(X(:));
if show_TI == 1
    figure(1);clf;
    ViewGrid(X)
    colormap default
    axis equal tight
    colorbar
    caxis([cmin, cmax])
    view(30,30)
    title('TI');
end

%% Showing the output realizations
if show_realizations == 1
    for i=1:nbrealz
        figure(1+i);clf;colormap gray;hold on
        ViewGrid(Y(:,:,:,i))
        colormap default
        colorbar
        caxis([cmin, cmax])
        axis equal tight
        view(30,30)
        title('Realization');
    end
end

%% Plotting the Satisfied conditioning points vs Given field data
if c > 0
    figure(2+nbrealz);clf;hold on
    meanerror=zeros(1,nbrealz);
    stderror=zeros(1,nbrealz);
    for i=1:nbrealz
        plot(D(:,4),Satisfied_conditioning_points(:,:,i),'o')
        meanerror(i)=mean(Satisfied_conditioning_points(:,:,i) - D(:,4));
        stderror(i)=std(Satisfied_conditioning_points(:,:,i) - D(:,4));
    end
    plot([min(D(:,4)),max(D(:,4))],[min(D(:,4)),max(D(:,4))],'r')
    disp(['Mean conditioning error over all realizations: ',num2str(mean(meanerror))])
    disp(['Mean std error over all realizations: ',num2str(mean(stderror))])
    axis equal square
    axis tight
    xlabel('Conditioning Data (Grayscale)')
    ylabel('Simulated Data (Grayscale)')
end

% %% Showing the E-type and standard error maps   
% meanrealz=mean(Y,3);
% stdrealz=std(Y,[],3);
% figure(3+nbrealz);clf;hold on
% %     subplot(2,1,1);hold on
% surf(xx,yy,flipud(double(meanrealz)));
% if c > 0
% %     scatter(D{1,1}(:,2),nx-D{1,1}(:,1),3,D{1,1}(:,3),'o','filled')
%     scatter3(D{1,1}(:,2),nx-D{1,1}(:,1),ones(size(D{1,1}(:,1)))+200,60,D{1,1}(:,3),'MarkerEdgeColor','r')
% %     plot3(D(:,2),nx-D(:,1),ones(size(D(:,1)))+200,'or');
% end
% % view(0,90)
% % shading flat
% axis equal tight
% colormap default
% colorbar
% title('E-type map');
% 
% %     subplot(2,1,2);hold on
% figure(4+nbrealz);clf;hold on
% surf(xx,yy,flipud(double(stdrealz)));
% if c > 0
% %     scatter(D{1,1}(:,2),nx-D{1,1}(:,1),3,D{1,1}(:,3),'MarkerEdgeColor','r')
%     scatter3(D{1,1}(:,2),nx-D{1,1}(:,1),ones(size(D{1,1}(:,1)))+200,60,D{1,1}(:,3),'MarkerEdgeColor','r')
% %     plot3(D(:,2),nx-D(:,1),ones(size(D(:,1)))+200,'or');
% end
% view(0,90)
% shading flat
% axis equal tight
% colormap default
% colorbar
% title('Standard daviation map');


% %% Showing the Variogram of TI and output realizations
% if show_variogram==1
% % Variogram of TI and output realizations
%     tic
%     figure(5+nbrealz);clf;hold on
%     [xxti,yyti]=meshgrid(1:size(X,2),1:size(X,1));
%     pts_ti=Grid2Pts(X,xxti,yyti,zeros(size(xxti)));
%     sampling_ti = round(size(pts_ti,1)/1000);
%     [meanhti,gammahti]=Variogram(pts_ti,20,60,sampling_ti,0);
%     pts = zeros(p*q,1+c,nbrealz);
%     for i=1:nbrealz
%         pts=Grid2Pts(Y(:,:,i),xx,yy,zeros(size(xx)));
%         sampling = round(size(pts,1)/1000);
%         [meanh,gammah]=Variogram(pts,20,60,sampling,0);
%         plot(meanh,gammah,'k--')
%     end
%     plot(meanhti,gammahti,'r','LineWidth',2)
%     
%     axis equal square
%     axis tight
% %     ylim([0 1200])
%     xlabel('Distance')
%     ylabel('Variance')
%     fprintf('To draw the Variogram: ');
%     toc
% end

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
    [xxti,yyti,zzti]=meshgrid(1:size(X,3),1:size(X,2),1:size(X,1));
    pts_ti=Grid2Pts(X,xxti,yyti,zzti);
    sampling_ti = round(size(pts_ti,1)/10000);
    [meanhti,gammahti]=Variogram(pts_ti,20,50,sampling_ti,0);
    pts = zeros(p*q,1+c,nbrealz);
    for i=1:nbrealz
        [xx,yy,zz]=meshgrid(1:size(Y(:,:,:,i),3),1:size(Y(:,:,:,i),2),1:size(Y(:,:,:,i),1));
        pts=Grid2Pts(Y(:,:,:,i),xx,yy,zz);
        sampling = round(size(pts,1)/10000);
        [meanh,gammah]=Variogram(pts,20,50,sampling,0);
        plot(meanh,gammah,'k--')
    end
    plot(meanhti,gammahti,'r','LineWidth',2)
    
    axis equal square
    axis tight
    ylim([0 2])
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

% %% Plotting the locations of field Data
% if c > 0
%     if show_data==1
%         figure(7+nbrealz);clf;colormap gray;hold on
% %         scatter(D{1,1}(:,1),D{1,1}(:,2),30,D{1,1}(:,3),'o','filled')
%         nx=size(Y,1);
%         scatter3(D(:,2),nx-D(:,1),D(:,3),30,D(:,4),'o','filled','MarkerEdgeColor','r')
% %         scatter(D{1,1}(:,2),nx-D{1,1}(:,1),ones(size(D{1,1}(:,1)))+200,40,D{1,1}(:,3),'filled','MarkerEdgeColor','r')
% %         view(0,90)
% %         shading flat
%         axis equal tight
%         colormap default
%         colorbar 
%         view(30,30)
%     end
% end

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

% profile off

% WriteGrid(X,'TI.SGEMS', 'SGEMS');
% WriteGrid(Y(:,:,:,1),'Realzation_1.SGEMS', 'SGEMS');
% WriteGrid(Y(:,:,:,2),'Realzation_2.SGEMS', 'SGEMS');
% WriteGrid(Y(:,:,:,3),'Realzation_3.SGEMS', 'SGEMS');
% WriteGrid(Y(:,:,:,1),'Realzation_4.SGEMS', 'SGEMS');
% WriteGrid(Y(:,:,:,1),'Realzation_5.SGEMS', 'SGEMS');
% WritePts(D,'Data_3D.SGEMS', 'SGEMS');
