%function Z = ssd(X, Y)
%Computes the sum of squared distances between X and Y for each possible
% overlap of Y on X.  Y is thus smaller than X
%
%Inputs:
%   X - larger image
%   Y - smaller image
%
%Outputs:
%   Each pixel of Z contains the ssd for Y overlaid on X at that pixel

function Z = ssd_v4_3D(X, Y)


K = ones(size(Y,1), size(Y,2), size(Y,3));

% for t=1:size(X,4), % No of TI layers 
% for k=1:size(K,4), % No of variables
% A = X(:,:,:);
% B = Y(:,:,:);

% a2 = imfilter(K, A.^2, 'symmetric', 'full');
% ab = imfilter(B, A, 'symmetric', 'full').*2;
% a2 = imfilter(X.^2, K);
a2 = filter3_v1(K, X.^2);
% a22 = convn(K, X.^2, 'full');
% a22 = a22(ceil((size(a22,1)-size(X,1)+1)/2):end-floor((size(a22,1)-size(X,1)+1)/2), ceil((size(a22,2)-size(X,2)+1)/2):end-floor((size(a22,2)-size(X,2)+1)/2)...
%                     , ceil((size(a22,3)-size(X,3)+1)/2):end-floor((size(a22,3)-size(X,3)+1)/2));
                
% figure(1);clf
% subplot(1,2,1)
% ViewGrid(a2)
% colorbar
% subplot(1,2,2)
% ViewGrid(a22)
% colorbar

% a22 = a22(1:size(X,1), 1:size(X,2), 1:size(X,3));

                
% figure(2);clf
% subplot(1,2,1)
% ViewGrid(a2)
% colorbar
% subplot(1,2,2)
% ViewGrid(a22)
% colorbar

b2 = sum(sum(sum(Y.^2)));
% ab = imfilter(X, Y).*2;
ab = filter3_v1(Y, X).*2;
% ab2 = convn(Y, X, 'full').*2;

% figure(3);clf
% subplot(1,2,1)
% ViewGrid(ab)
% colorbar
% subplot(1,2,2)
% ViewGrid(ab2)
% colorbar

% ab2 = ab2(1:size(X,1), 1:size(X,2), 1:size(X,3));
% ab2 = ab2(ceil((size(ab2,1)-size(X,1)+1)/2):end-floor((size(ab2,1)-size(X,1)+1)/2), ceil((size(ab2,2)-size(X,2)+1)/2):end-floor((size(ab2,2)-size(X,2)+1)/2)...
%                     , ceil((size(ab2,3)-size(X,3)+1)/2):end-floor((size(ab2,3)-size(X,3)+1)/2));
% figure(4);clf
% subplot(1,2,1)
% ViewGrid(ab)
% colorbar
% subplot(1,2,2)
% ViewGrid(ab2)
% colorbar

% a2 = a2(1:size(z,1),1:size(z,2),1:size(z,3));
% ab = ab(1:size(z,1),1:size(z,2),1:size(z,3));

% for i=1:size(K,3)
%     %             z(:,:,i,k) = ((a2(:,:,i) - ab(:,:,i)) + b2(:,:,i));
%     z(:,:,i) = ((a2(:,:,i) - ab(:,:,i)) + b2(:,:,i));
% end
%         z(:,:,k) = w_v(k)*((a2 - ab) + b2);
% end;
%     Z(:,:,t) = sum(z,3); 
% end
% z(:,:,:) = (a2(:,:,:) - ab(:,:,:)) + b2(1,1,:);

z = ((a2 - ab) + b2);
% z2 = ((a22 - ab2) + b2);

% z = z(1:size(z,1),1:size(z,2),1:size(z,3));

Z=max(z,0);
Z=sqrt(Z);

% Z2=max(z2,0);
% Z2=sqrt(Z2);

% figure(5);clf
% subplot(1,2,1)
% ViewGrid(Z)
% colorbar
% subplot(1,2,2)
% ViewGrid(Z2)
% colorbar

% figure(101);clf;hold on
% subplot(2,2,1)
% ViewGrid(A)
% colorbar
% subplot(2,2,2)
% ViewGrid(B)
% colorbar
% subplot(2,2,4)
% ViewGrid(Z);
% colorbar
% 'stop';
