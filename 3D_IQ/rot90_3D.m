function M = rot90_3D(M, Dim, NumRot)
% Extends the rot90 command to 3-Dimensions.
% The input matrix M is assumed to be a 3D matrix. This matrix is
% rotated by 90degrees (or 180 or 270) about any of the 3 axes.
% Dim is the dimension to rotate normal to.
% NumRot is the number of 90-degree rotations.
% M must be specified as input matrix.
% Dim should be an integer in the range 1..3 for 3 dimensions.
% if Dim is not specified, then Dim=1. 
% If NumRot not specified then NumRot=1.

aSize = size(M);

switch Dim
  case 1
      switch NumRot
        case 1
            X = permute(M, [1, 3, 2]);
            M = X(:, aSize(3):-1:1, :);
        case 2
            M = M(:, aSize(2):-1:1, aSize(3):-1:1);
        case 3
            M = permute(M(:, aSize(2):-1:1, :), [1, 3, 2]);
      end
      
  case 2
      switch NumRot
        case 1
            X = permute(M, [3, 2, 1]);
            M = X(aSize(3):-1:1, :, :);
        case 2
            M = M(aSize(1):-1:1, :, aSize(3):-1:1);
        case 3
            M = permute(M(aSize(1):-1:1, :, :), [3, 2, 1]);
      end
      
  case 3
      switch NumRot
        case 1
            X = permute(M, [2, 1, 3]);
            M = X(aSize(2):-1:1, :, :);
        case 2
            M = M(aSize(1):-1:1, aSize(2):-1:1, :);
        case 3
            M = permute(M(aSize(1):-1:1, :, :), [2, 1, 3]);
      end
      
  otherwise
      error('Dim must be 1, 2 or 3');
end

return;