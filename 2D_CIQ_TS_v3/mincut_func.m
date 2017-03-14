% -----------------------------------------------------------------------------
% This function can be used for finding the minimum boundary

% Author: Pejman Tahmasebi
% Reference: Tahmasebi, P., Hezarkhani, A., Sahimi, M., 2012. 
% Multiple-point geostatistical modeling based on the cross-correlation
% functions, Journal of Computational Geosciences, 
% doi:10.1007/s10596-012 9287-1

%-----------------------------------------------------------------------------------

function [ M ] = mincut_func( Target, imout, T, OL, i, j )

%% Input Parameters
% - Target: one side of pattern (selected pattern)
% - imout: another side of pattern (previous pattern) 
% - T: Template size
% - OL: Overlap size
% - i: i-th row of output 
% - j: j-th column of output

%% Output Parameters
% - M: minimum error boundary

%----------------------------------------------------------------------------------

      M = ones(T, T);
      
      if (j>1)
      E = ( Target(1:1+T-1, 1:1+OL-1) - imout(1:1+T-1, 1:1+OL-1) ).^2;
      C = mincut(E, 0);
      M(1:end, 1:OL) = double(C >= 0);
      end

      if (i>1)
      E = ( Target(1:1+OL-1, 1:1+T-1) - imout(1:1+OL-1, 1:1+T-1) ).^2;          
      C = mincut(E, 1);
      M(1:OL, 1:end) = M(1:OL, 1:end) .* double(C >= 0);
      end;

end

