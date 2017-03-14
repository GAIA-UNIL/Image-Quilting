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

function Z = ssd_v4(X, Y, w_v, z, Z)

K = ones(size(Y,1), size(Y,2));

for t=1:size(X,4), % No of TI layers    ???????HERE is this okay???????
    for k=1:size(X,3), % No of variables
        A = X(:,:,k,t);
        B = Y(:,:,k);
        
        a2 = filter2(K, A.^2, 'valid');
        b2 = sum(sum(B.^2));
        ab = filter2(B, A, 'valid').*2;
        z(:,:,k) = w_v(k)*((a2 - ab) + b2); 
    end;
    Z(:,:,t) = sum(z,3); %    ???????HERE is this okay???????
end

Z=max(Z,0); %    ???????HERE is this okay???????
Z=sqrt(Z);

