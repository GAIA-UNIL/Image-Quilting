function [ M ] = mincut_func( Target, imout, T, OL, i, j )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

      M = ones(T, T);
      
      if (j>1)
      % Left OL
      %Compute the error in the border region
      E = ( Target(1:1+T-1, 1:1+OL-1) - imout(1:1+T-1, 1:1+OL-1) ).^2;
                 
      %Compute the mincut array
      C = mincut(E, 0);
                 
      %Compute the mask and write to the destination
      M(1:end, 1:OL) = double(C >= 0);
      end

      if (i>1)
      %We have a bottom OL
      %Compute the SSD in the border region
      E = ( Target(1:1+OL-1, 1:1+T-1) - imout(1:1+OL-1, 1:1+T-1) ).^2;
                 
      %Compute the mincut array
      C = mincut(E, 1);
                 
      %Compute the mask and write to the destination
      M(1:OL, 1:end) = M(1:OL, 1:end) .* double(C >= 0);
      end;

end

