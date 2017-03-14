function y = filter3_v1(b,x,shape)
%FILTER2 Three-dimensional digital filter.
%   Y = FILTER3(B,X) filters the data in X with the 3-D FIR
%   filter in the matrix B.  The result, Y, is computed
%   using 3-D correlation and is the same size as X.
%
%   Y = FILTER3(B,X,'shape') returns Y computed via 3-D
%   correlation with size specified by 'shape':
%     'same'  - (default) returns the central part of the
%               correlation that is the same size as X.
%     'valid' - returns only those parts of the correlation
%               that are computed without the zero-padded
%               edges, size(Y) < size(X).
%     'full'  - returns the full 3-D correlation,
%               size(Y) > size(X).
%
%   FILTER3 uses CONVN to do most of the work.  3-D correlation
%   is related to 3-D convolution by a 180 degree rotation of the
%   filter matrix.
%
%   Class support for inputs B,X:
%      float: double, single
%
%   See also FILTER, CONV2.

%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 5.13.4.5 $  $Date: 2010/09/02 13:35:21 $

error(nargchk(2,3,nargin,'struct'));
if nargin<4, shape = 'same'; end

if (~isa(b,'float')),  b = double(b); end % True if object is a given class
if (~isa(x,'float')),  x = double(x); end

code = [shape,' ']; code = code(1);
if isempty(find(code=='svf', 1)) % Check the size specification value 'shape': same/valid/full
    error(message('MATLAB:filter2:InvalidParam'));
end


% [mx,nx,ox] = size(x);
stencil = zeros(size(b,1),size(b,2),size(b,3));
size_b = [size(b,1) size(b,2) size(b,3)];

%% Rotate b by 2*90 degree
% if (min(size_b) == size(b,3))
%     for i=1:size(b,3)
%         stencil(:,:,i) = rot90(b(:,:,i),2); % Rotate 2*90 degree: Why??????
%     end
% end
if (min(size_b) == size(b,3)) % If we have K directional overlap
   stencil = rot90_3D(b,3,2); % Rotate 2*90 degree: Why??????
end
% if (min(size_b) == size(b,2))
%     for i=1:size(b,2)
%         stencil(:,i,:) = rot90(b(:,i,:),2); % Rotate 2*90 degree: Why??????
%     end
% end
if (min(size_b) == size(b,2)) % If we have J directional overlap
   stencil = rot90_3D(b,2,2); % Rotate 2*90 degree: Why??????
end
% if (min(size_b) == size(b,1))
%     for i=1:size(b,1)
%         stencil(i,:,:) = rot90(b(i,:,:),2); % Rotate 2*90 degree: Why??????
%     end
% end
if (min(size_b) == size(b,1)) % If we have I directional overlap
   stencil = rot90_3D(b,1,2); % Rotate 2*90 degree: Why??????
end
% [ms,ns,os] = size(stencil);
y = convn(x,stencil,shape);

% %%
% % 1-D stencil?
% % if (ms == 1)
% %   y = convn(1,stencil,x,shape);     % Not sure??????
% % elseif (ns == 1)
% %   y = convn(stencil,1,x,shape);     % Not sure??????
% % elseif (os == 1)
% %   y = convn(stencil,1,x,shape);     % Not sure??????
% % else
% if (ms*ns*os > mx*nx*ox)
%     % The filter is bigger than the input.  This is a nontypical
%     % case, and it may be counterproductive to check the
%     % separability of the stencil.
%     if (min(size_b) == size(b,3)) % If we have K directional overlap
%         for i=1:size(b,3)
%             y(:,:,i) = conv2(x(:,:,i),stencil(:,:,i),shape);     % Not sure??????
%         end
%     end
%     if (min(size_b) == size(b,2)) % If we have J directional overlap
%         for i=1:size(b,2)
%             y(:,i,:) = conv2(x(:,i,:),stencil(:,i,:),shape);     % Not sure??????
%         end
%     end
%     if (min(size_b) == size(b,1)) % If we have I directional overlap
%         for i=1:size(b,1)
%             y(i,:,:) = conv2(x(i,:,:),stencil(i,:,:),shape);     % Not sure??????
%         end
%     end
% else
%     separable = false;
%     if all(isfinite(stencil(:)))
%         % Check rank (separability) of stencil
%         u = zeros(size(b,1),size(b,2),size(b,3));
%         S = zeros(size(b,1),size(b,2),size(b,3));
%         v = zeros(size(b,1),size(b,2),size(b,3));
%         if (min(size_b) == size(b,3)) % If we have K directional overlap
%             for i=1:size(b,3)
%                 [u(:,:,i),S(:,:,i),v(:,:,i)] = svd(stencil(:,:,i)); 
%                 % Singular value decomposition: svd(X) produces a diagonal matrix S, of the same
%                 %   dimension as X and with nonnegative diagonal elements in decreasing order,
%                 %   and unitary matrices U and V so that X = U*S*V'
%             end
%         end
%         %         if (min(size_b) == size(b,3))
%         %             R=get_SVD1(b);
%         %         end
%         
%         
%         if (min(size_b) == size(b,2)) % If we have J directional overlap
%             for i=1:size(b,2)
%                 [u(:,i,:),S(:,i,:),v(:,i,:)] = svd(stencil(:,i,:));
%             end
%         end
%         if (min(size_b) == size(b,1)) % If we have I directional overlap
%             for i=1:size(b,1)
%                 [u(i,:,:),S(i,:,:),v(i,:,:)] = svd(stencil(i,:,:));
%             end
%         end
%         
%         s = zeros(size(b,1),1,size(b,3));
%         if (min(size_b) == size(b,3)) % If we have K directional overlap
%             for i=1:size(b,3)
%                 s(:,:,i) = diag(S(:,:,i));
%             end
%         end
%         if (min(size_b) == size(b,2)) % If we have J directional overlap
%             for i=1:size(b,2)
%                 s(:,i,:) = diag(S(:,i,:));
%             end
%         end
%         if (min(size_b) == size(b,1)) % If we have I directional overlap
%             for i=1:size(b,1)
%                 s(i,:,:) = diag(S(i,:,:));
%             end
%         end
%         
%         tol = length(stencil) * eps(max(s(:)));
%         rank = sum(s > tol);
%         separable = (rank ==1);
%     end
%     %% If we have K directional overlap
%     if (min(size_b) == size(b,3)) 
%         if separable
%         % Separable stencil
%             y = zeros(size(x,1), size(x,2), size(b,3));
%             for i=1:size(b,3)
%                 hcol = u(:,1,i) * sqrt(s(1));
%                 hrow = conj(v(:,1,i)) * sqrt(s(1));
%                 if (all(all(all((round(stencil) == stencil)))) && all(all(all((round(x) == x)))))
%                     % Output should be integer
% %                     y = round(convn(hcol, hrow, x, shape));
%                     y(:,:,i) = round(conv2(hcol, hrow, x(:,:,i), shape));
%                 else
% %                     y = conv2(hcol, hrow, x, shape);
%                     y(:,:,i) = conv2(hcol, hrow, x(:,:,i), shape);     % Not sure??????
%                 end
%             end
%         else
% %             for i=1:size(b,3)
%                 % Nonseparable stencil
%                 y = convn(x,stencil,shape);
% %                 y(:,:,i) = convn(x(:,:,i),stencil(:,:,i),shape);     % Not sure??????
% %             end
%         end
%     end
%     %% If we have J directional overlap
%     if (min(size_b) == size(b,2))
%         if separable
%         % Separable stencil
%         y = zeros(size(x,1), size(b,2), size(x,3));
%             for i=1:size(b,2)
%                 hcol = u(:,i,1) * sqrt(s(1));
%                 hrow = conj(v(:,i,1)) * sqrt(s(1));
%                 if (all(all(all((round(stencil) == stencil)))) && all(all(all((round(x) == x)))))
%                     % Output should be integer
%                     y(:,i,:) = round(conv2(hcol, hrow, x(:,i,:), shape));
%                 else
%                     y(:,i,:) = conv2(hcol, hrow, x(:,i,:), shape);     % Not sure??????
%                 end
%             end
%         else
% %             for i=1:size(b,2)
%                 % Nonseparable stencil
%                 y = convn(x,stencil,shape);
% %                 y(:,i,:) = conv2(x(:,i,:),stencil(:,i,:),shape);     % Not sure??????
% %             end
%         end
%     end
%     %% If we have I directional overlap
%     if (min(size_b) == size(b,1))
%         if separable
%         % Separable stencil
%         y = zeros(size(b,1), size(x,2), size(x,3));
%             for i=1:size(b,1)
%                 hcol = u(:,1,i) * sqrt(s(1));
%                 hrow = conj(v(:,1,i)) * sqrt(s(1));
%                 if (all(all(all((round(stencil) == stencil)))) && all(all(all((round(x) == x)))))
%                     % Output should be integer
%                     y(:,:,i) = round(conv2(hcol, hrow, x(:,:,i), shape));
%                 else
%                     y(:,:,i) = conv2(hcol, hrow, x(:,:,i), shape);     % Not sure??????
%                 end
%             end
%         else
% %             for i=1:size(b,3)
%                 % Nonseparable stencil
%                 y = convn(x,stencil,shape);
% %                 y(:,:,i) = conv2(x(:,:,i),stencil(:,:,i),shape);     % Not sure??????
% %             end
%         end
%     end
%     
%     
% end
% end
