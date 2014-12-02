function [stack, stackinfo] = stackRead(stackpath)
%
% [stack, stackinfo] = stackRead(stackpath)
%
% Wrapper function for Fran�ois Nedelec's tiffread.m
% This works on STK files as well as multipage TIFFs.
% Outputs:
%
%   stack     : 3D data array
%   stackinfo : Any other information included in original file
%
% Francois Aguet, 01/2010
%
% Copyright (C) 2012 LCCB 
%
% This file is part of QFSM.
% 
% QFSM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% QFSM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with QFSM.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

if (nargin == 0 || isempty(stackpath))
    stackinfo = tiffread();
else
    stackinfo = tiffread(stackpath);
end

stack = cat(3,stackinfo(:).data);
stackinfo = rmfield(stackinfo, 'data');