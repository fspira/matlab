function statFunc = regiondata(varargin)
% FUNC = regiondata(LABELED_IMAGE, OPTIONS)
%   Create a function that returns data a la REGIONPROPS for a given index
%   in LABELED_IMAGE.  This saves a great deal on memory and fragmentation,
%   as REGIONPROP stat structs are not created for every index and are
%   instead viewed one at a time.  OPTIONS are REGIONPROPS like options.
%
%   PolygonCorners is also added, which are the corners of the polygon that
%   make up a given image (usually used in convex hull calculation).
%
%   To get all the indices, use unique(lbl(:)), and take out the zero.
%        ie
%      ids = unique(labeledImage);
%      ids = ids(ids ~= 0);
%
%   See regionprops for more info.
%

%   ** BASED ON WORK BELONGING TO MATHWORKS **
%   ** THIS IS THEIR COPYRIGHT NOTICE FROM REGIONPROPS **
%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision.4.2.3 $  $Date: 2006/05/24 03:32:56 $

  officialStats = {'Area'
                   'Centroid'
                   'BoundingBox'
                   'SubarrayIdx'
                   'MajorAxisLength'
                   'MinorAxisLength'
                   'Eccentricity'
                   'Orientation'
                   'ConvexHull'
                   'ConvexImage'
                   'ConvexArea'
                   'Image'
                   'FilledImage'
                   'FilledArea'
                   'EulerNumber'
                   'Extrema'
                   'EquivDiameter'
                   'Solidity'
                   'Extent'
                   'PixelIdxList'
                   'PixelList'
                   'Perimeter'
                   'PolygonCorners'
                   };

%  tempStats = {'PerimeterCornerPixelList'};
%  allStats = [officialStats; tempStats];
  allStats = officialStats;

  [L, requestedStats] = ParseInputs(officialStats, varargin{:});

  if ndims(L) > 2
    % Remove stats that aren't supported for N-D input and issue
    % warning messages as appropriate.
    requestedStats = PreprocessRequestedStats(requestedStats);
  end

  if isempty(requestedStats)
    eid = sprintf('Images:%s:noPropertiesWereSelected',mfilename);
    msg = 'No input properties';
    error(eid,'%s',msg);
  end

  if (isempty(L))
    numObjs = 0;
  else
    numObjs = round(double(max(L(:))));
  end

  % Initialize the stats structure array.
  numStats = length(allStats);
  empties = cell(numStats, 1);
  zz = cell(numStats, 1);
  for z = 1:numStats
    zz{z} = 0;
  end

  function stats = ComputeStats(idx)
    if idx > numObjs
      error('RegionData:Index', 'Index %d is out of range.  %d total objects.', idx, numObjs);
    end

    computedStats = cell2struct(zz, allStats, 1);
    stats = cell2struct(empties, allStats, 1);

    computedStats.PixelIdxList = 1;
    stats.PixelIdxList = find(L == idx);

    stats.objectIdx = idx;

    for k = 1:length(requestedStats)
      switch requestedStats{k}
        case {'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Eccentricity'}
          [stats, computedStats] = ComputeEllipseParams(L, stats, computedStats);
        otherwise
          f = str2func(['Compute', requestedStats{k}]);
          [stats, computedStats] = f(L, stats, computedStats);
      end
    end
  end

  statFunc = @ComputeStats;
end

function [stats, computedStats] = ComputePixelIdxList(L, stats, computedStats)
  if ~computedStats.PixelIdxList
    % This should just plain never happen - but just in case we mess up
    % with what's above us, let's leave it as a safety.
    computedStats.PixelIdxList = 1;
    stats.PixelIdxList = find(L == stats.objectIdx);
  end
end

function [stats, computedStats] = ComputeArea(L, stats, computedStats)
%   The area is defined to be the number of pixels belonging to
%   the region.
  if ~computedStats.Area
    computedStats.Area = 1;
    stats.Area = numel(stats.PixelIdxList);
  end
end

function [stats, computedStats] = ComputeEquivDiameter(L, stats, computedStats)
%   Computes the diameter of the circle that has the same area as
%   the region.
%   Ref: Russ, The Image Processing Handbook, 2nd ed, 1994, page
%   511.

  if ~computedStats.EquivDiameter
    computedStats.EquivDiameter = 1;
    
    if ndims(L) > 2
      NoNDSupport('EquivDiameter');
      return
    end

    [stats, computedStats] = ComputeArea(L, stats, computedStats);

    factor = 2/sqrt(pi);
    stats.EquivDiameter = factor * sqrt(stats.Area);
  end
end

function [stats, computedStats] = ComputeFilledImage(L, stats, computedStats)
%   Uses imfill to fill holes in the region.

  if ~computedStats.FilledImage
    computedStats.FilledImage = 1;
    
    [stats, computedStats] = ComputeImage(L, stats, computedStats);
    
    conn = conndef(ndims(L),'minimal');
    stats.FilledImage = imfill(stats.Image, conn, 'holes');
  end
end

function [stats, computedStats] = ComputeConvexArea(L, stats, computedStats)
%   Computes the number of "on" pixels in ConvexImage.

  if ~computedStats.ConvexArea
    computedStats.ConvexArea = 1;
    
    if ndims(L) > 2
      NoNDSupport('ConvexArea');
      return
    end
    
    [stats, computedStats] = ComputeConvexImage(L, stats, computedStats);
    stats.ConvexArea = sum(stats.ConvexImage(:));
  end
end

function [stats, computedStats] = ComputeFilledArea(L, stats, computedStats)
%   Computes the number of "on" pixels in FilledImage.

  if ~computedStats.FilledArea
    computedStats.FilledArea = 1;
    
    [stats, computedStats] = ComputeFilledImage(L,stats,computedStats);

    stats.FilledArea = sum(stats.FilledImage(:));
  end
end


function [stats, computedStats] = ComputeConvexImage(L, stats, computedStats)
%   Uses ROIPOLY to fill in the convex hull.

  if ~computedStats.ConvexImage
    computedStats.ConvexImage = 1;
    
    if ndims(L) > 2
      NoNDSupport('ConvexImage');
      return
    end
    [stats, computedStats] = ComputeConvexHull(L, stats, computedStats);
    [stats, computedStats] = ComputeBoundingBox(L, stats, computedStats);
    
    M = stats.BoundingBox(4);
    N = stats.BoundingBox(3);
    hull = stats.ConvexHull;
    if (isempty(hull))
      stats.ConvexImage = false(M,N);
    else
      firstRow = stats.BoundingBox(2) + 0.5;
      firstCol = stats.BoundingBox(1) + 0.5;
      r = hull(:,2) - firstRow + 1;
      c = hull(:,1) - firstCol + 1;
      stats.ConvexImage = roipoly(M, N, c, r);
    end
  end
end

function [stats, computedStats] = ComputeCentroid(L, stats, computedStats)
%   [mean(r) mean(c)]

  if ~computedStats.Centroid
    computedStats.Centroid = 1;
    
    [stats, computedStats] = ComputePixelList(L, stats, computedStats);

    % Save the warning state and disable warnings to prevent divide-by-zero
    % warnings.
    warning off MATLAB:divideByZero;
    stats.Centroid = mean(stats.PixelList, 1);
    
    % Restore the warning state.
    warning on MATLAB:divideByZero;
  end
end

function [stats, computedStats] = ComputeEulerNumber(L, stats, computedStats)
%   Calls BWEULER on 'Image' using 8-connectivity

  if ~computedStats.EulerNumber
    computedStats.EulerNumber = 1;
    
    if ndims(L) > 2
      NoNDSupport('EulerNumber');
      return
    end
    
    [stats, computedStats] = ComputeImage(L, stats, computedStats);
    
    stats.EulerNumber = bweuler(stats.Image,8);
  end
end

function [stats, computedStats] = ComputeExtrema(L, stats, computedStats)
%   A 8-by-2 array; each row contains the x and y spatial
%   coordinates for these extrema:  leftmost-top, rightmost-top,
%   topmost-right, bottommost-right, rightmost-bottom, leftmost-bottom,
%   bottommost-left, topmost-left. 
%   reference: Haralick and Shapiro, Computer and Robot Vision
%   vol I, Addison-Wesley 1992, pp. 62-64.

  if ~computedStats.Extrema
    computedStats.Extrema = 1;
    
    if ndims(L) > 2
      NoNDSupport('Extrema');
      return
    end
    
    [stats, computedStats] = ComputePixelList(L, stats, computedStats);
    
    pixelList = stats.PixelList;
    if (isempty(pixelList))
      stats.Extrema = zeros(8,2) + 0.5;
    else
      r = pixelList(:,2);
      c = pixelList(:,1);

      minR = min(r);
      maxR = max(r);
      minC = min(c);
      maxC = max(c);

      minRSet = r==minR;
      maxRSet = r==maxR;
      minCSet = c==minC;
      maxCSet = c==maxC;

      % Points 1 and 2 are on the top row.
      r1 = minR;
      r2 = minR;
      % Find the minimum and maximum column coordinates for
      % top-row pixels.
      tmp = c(minRSet);
      c1 = min(tmp);
      c2 = max(tmp);

      % Points 3 and 4 are on the right column.
      % Find the minimum and maximum row coordinates for
      % right-column pixels.
      tmp = r(maxCSet);
      r3 = min(tmp);
      r4 = max(tmp);
      c3 = maxC;
      c4 = maxC;

      % Points 5 and 6 are on the bottom row.
      r5 = maxR;
      r6 = maxR;
      % Find the minimum and maximum column coordinates for
      % bottom-row pixels.
      tmp = c(maxRSet);
      c5 = max(tmp);
      c6 = min(tmp);

      % Points 7 and 8 are on the left column.
      % Find the minimum and maximum row coordinates for
      % left-column pixels.
      tmp = r(minCSet);
      r7 = max(tmp);
      r8 = min(tmp);
      c7 = minC;
      c8 = minC;

      stats.Extrema = [c1-0.5 r1-0.5
                       c2+0.5 r2-0.5
                       c3+0.5 r3-0.5
                       c4+0.5 r4+0.5
                       c5+0.5 r5+0.5
                       c6-0.5 r6+0.5
                       c7-0.5 r7+0.5
                       c8-0.5 r8-0.5];
    end
  end
end
  
function [stats, computedStats] = ComputeBoundingBox(L, stats, computedStats)
%   [minC minR width height]; minC and minR end in .5.

  if ~computedStats.BoundingBox
    computedStats.BoundingBox = 1;
    
    [stats, computedStats] = ComputePixelList(L, stats, computedStats);
      
    num_dims = ndims(L);
    
    list = stats.PixelList;
    if (isempty(list))
      stats.BoundingBox = [0.5*ones(1,num_dims) zeros(1,num_dims)];
    else
      min_corner = min(list,[],1) - 0.5;
      max_corner = max(list,[],1) + 0.5;
      stats.BoundingBox = [min_corner (max_corner - min_corner)];
    end
  end
end

function [stats, computedStats] = ComputeSubarrayIdx(L, stats, computedStats)
%   Find a cell-array containing indices so that L(idx{:}) extracts the
%   elements of L inside the bounding box.

  if ~computedStats.SubarrayIdx
    computedStats.SubarrayIdx = 1;
    
    [stats, computedStats] = ComputeBoundingBox(L, stats, computedStats);
    num_dims = ndims(L);
    idx = cell(1,num_dims);

    boundingBox = stats.BoundingBox;
    left = boundingBox(1:(end/2));
    right = boundingBox((1+end/2):end);
    left = left(1,[2 1 3:end]);
    right = right(1,[2 1 3:end]);
    for p = 1:num_dims
      first = left(p) + 0.5;
      last = first + right(p) - 1;
      idx{p} = first:last;
    end
    stats.SubarrayIdx = idx;
  end
end


function [stats, computedStats] = ComputeEllipseParams(L, stats, ...
                                                    computedStats)  
%   Find the ellipse that has the same normalized second central moments as the
%   region.  Compute the axes lengths, orientation, and eccentricity of the
%   ellipse.  Ref: Haralick and Shapiro, Computer and Robot Vision vol I,
%   Addison-Wesley 1992, Appendix A.


  if ~(computedStats.MajorAxisLength && computedStats.MinorAxisLength && ...
       computedStats.Orientation && computedStats.Eccentricity)
    computedStats.MajorAxisLength = 1;
    computedStats.MinorAxisLength = 1;
    computedStats.Eccentricity = 1;
    computedStats.Orientation = 1;
    
    if ndims(L) > 2
      NoNDSupport({'MajorAxisLength', 'MinorAxisLength', ...
                   'Eccentricity', 'Orientation'});
      return
    end
    
    [stats, computedStats] = ComputePixelList(L, stats, computedStats);
    [stats, computedStats] = ComputeCentroid(L, stats, computedStats);

    % Disable divide-by-zero warning
    warning off MATLAB:divideByZero;
    
    list = stats.PixelList;
    if (isempty(list))
      stats.MajorAxisLength = 0;
      stats.MinorAxisLength = 0;
      stats.Eccentricity = 0;
      stats.Orientation = 0;

    else
      % Assign X and Y variables so that we're measuring orientation
      % counterclockwise from the horizontal axis.

      xbar = stats.Centroid(1);
      ybar = stats.Centroid(2);

      x = list(:,1) - xbar;
      y = -(list(:,2) - ybar); % This is negative for the 
                               % orientation calculation (measured in the
                               % counter-clockwise direction).

      N = length(x);

      % Calculate normalized second central moments for the region. 1/12 is
      % the normalized second central moment of a pixel with unit length.
      uxx = sum(x.^2)/N + 1/12; 
      uyy = sum(y.^2)/N + 1/12;
      uxy = sum(x.*y)/N;

      % Calculate major axis length, minor axis length, and eccentricity.
      common = sqrt((uxx - uyy)^2 + 4*uxy^2);
      stats.MajorAxisLength = 2*sqrt(2)*sqrt(uxx + uyy + common);
      stats.MinorAxisLength = 2*sqrt(2)*sqrt(uxx + uyy - common);
      stats.Eccentricity = 2*sqrt((stats.MajorAxisLength/2)^2 - ...
                                     (stats.MinorAxisLength/2)^2) / ...
          stats.MajorAxisLength;

      % Calculate orientation.
      if (uyy > uxx)
        num = uyy - uxx + sqrt((uyy - uxx)^2 + 4*uxy^2);
        den = 2*uxy;
      else
        num = 2*uxy;
        den = uxx - uyy + sqrt((uxx - uyy)^2 + 4*uxy^2);
      end
      if (num == 0) && (den == 0)
        stats.Orientation = 0;
      else
        stats.Orientation = (180/pi) * atan(num/den);
      end
    end
  end

  % Restore warning state.
  warning on MATLAB:divideByZero;
end
  
%%%
%%% ComputeSolidity
%%%
function [stats, computedStats] = ComputeSolidity(L, stats, computedStats)
%   Area / ConvexArea

  if ~computedStats.Solidity
    computedStats.Solidity = 1;
    
    if ndims(L) > 2
      NoNDSupport('Solidity');
      return
    end
    
    [stats, computedStats] = ComputeArea(L, stats, computedStats);
    [stats, computedStats] = ComputeConvexArea(L, stats, computedStats);
    
    if (stats.ConvexArea == 0)
      stats.Solidity = NaN;
    else
      stats.Solidity = stats.Area / stats.ConvexArea;
    end
  end
end

%%%
%%% ComputeExtent
%%%
function [stats, computedStats] = ComputeExtent(L, stats, computedStats)
%   Area / (BoundingBox(3) * BoundingBox(4))

  if ~computedStats.Extent
    computedStats.Extent = 1;
    
    if ndims(L) > 2
      NoNDSupport('Extent');
      return
    end
    
    [stats, computedStats] = ComputeArea(L, stats, computedStats);
    [stats, computedStats] = ComputeBoundingBox(L, stats, computedStats);
    
    if (stats.Area == 0)
      stats.Extent = NaN;
    else
      stats.Extent = stats.Area / prod(stats.BoundingBox(3:4));
    end
  end
end

%%%
%%% ComputeImage
%%%
function [stats, computedStats] = ComputeImage(L, stats, computedStats)
%   Binary image containing "on" pixels corresponding to pixels
%   belonging to the region.  The size of the image corresponds
%   to the size of the bounding box for each region.

  if ~computedStats.Image
    computedStats.Image = 1;

    [stats, computedStats] = ComputeSubarrayIdx(L, stats, computedStats);

    subarray = L(stats.SubarrayIdx{:});
    if ~isempty(subarray)
      stats.Image = (subarray == stats.objectIdx);
    else
      stats.Image = logical(subarray);
    end
  end
end


%%%
%%% ComputePixelList
%%%
function [stats, computedStats] = ComputePixelList(L, stats, computedStats)
%   A P-by-2 matrix, where P is the number of pixels belonging to
%   the region.  Each row contains the row and column
%   coordinates of a pixel.

  if ~computedStats.PixelList
    computedStats.PixelList = 1;
    
    % Convert the linear indices to subscripts and store
    % the results in the pixel list.  Reverse the order of the first
    % two subscripts to form x-y order.
    In = cell(1,ndims(L));

    if ~isempty(stats.PixelIdxList)
      [In{:}] = ind2sub(size(L), stats.PixelIdxList);
      stats.PixelList = [In{:}];
      stats.PixelList = stats.PixelList(:,[2 1 3:end]);
    else
      stats.PixelList = zeros(0,ndims(L));
    end
  end
end

%%%
%%% ComputePolygonCorners
%%%
function [stats, computedStats] = ComputePolygonCorners(L, ...
                                  stats, computedStats)
  %   Find the pixels on the perimeter of the region; make a list
  %   of the coordinates of their corners; sort and remove
  %   duplicates.

  if ~computedStats.PolygonCorners
    computedStats.PolygonCorners = 1;
    
    if ndims(L) > 2
      NoNDSupport('PolygonCorners');
      return
    end
    
    [stats, computedStats] = ComputeImage(L, stats, computedStats);
    [stats, computedStats] = ComputeBoundingBox(L, stats, computedStats);

    perimImage = bwmorph(stats.Image, 'perim8');
    firstRow = stats.BoundingBox(2) + 0.5;
    firstCol = stats.BoundingBox(1) + 0.5;
    [r,c] = find(perimImage);
    % Force rectangular empties.
    r = r(:) + firstRow - 1;
    c = c(:) + firstCol - 1;
    rr = [r-.5 ; r    ; r+.5 ; r   ];
    cc = [c    ; c+.5 ; c    ; c-.5];
    stats.PolygonCorners = [cc rr];
  end
end

%%%
%%% ComputeConvexHull
%%%
function [stats, computedStats] = ComputeConvexHull(L, stats, computedStats)
%   A P-by-2 array representing the convex hull of the region.
%   The first column contains row coordinates; the second column
%   contains column coordinates.  The resulting polygon goes
%   through pixel corners, not pixel centers.

  if ~computedStats.ConvexHull
    computedStats.ConvexHull = 1;
    
    if ndims(L) > 2
      NoNDSupport('ConvexHull');
      return
    end
    
    [stats, computedStats] = ComputePolygonCorners(L, stats, ...
                                                   computedStats);
    [stats, computedStats] = ComputeBoundingBox(L, stats, computedStats);

    list = stats.PolygonCorners;
    if (isempty(list))
      stats.ConvexHull = zeros(0,2);
    else
      rr = list(:,2);
      cc = list(:,1);
      hullIdx = convhull(rr, cc);
      stats.ConvexHull = list(hullIdx,:);
    end
  end
end

%%%
%%% ComputePerimeter
%%%
function [stats, computedStats] = ComputePerimeter(L, stats, computedStats)
  if ~computedStats.Perimeter
    computedStats.Perimeter = 1;
    
    if ndims(L) > 2
      NoNDSupport('ComputePerimeter');
      return
    end
    
    B = regionboundariesmex(double(L),8);

    boundary = B{stats.objectIdx};
    delta = diff(boundary).^2;
    stats.Perimeter = sum(sqrt(sum(delta, 2)));
  end
end

%%%
%%% ParseInputs
%%%
function [L,reqStats] = ParseInputs(officialStats, varargin)

  L = [];
  reqStats = [];

  if (length(varargin) < 1)
    eid = sprintf('Images:%s:tooFewInputs',mfilename);
    msg = 'Too few input arguments.';
    error(eid,'%s',msg);
  end

  L = varargin{1};

  if islogical(L)
    eid = 'Images:regionprops:binaryInput';
    msg1 = 'Use bwlabel(BW) or double(BW) convert binary image to ';
    msg2 = 'a label matrix before calling regionprops.';
    msg = sprintf('%s\n%s',msg1,msg2);
    error(eid, '%s', msg);
  end

  iptcheckinput(L, {'numeric'}, {'real', 'integer', 'nonnegative'}, ...
                mfilename, 'L', 1);

  list = varargin(2:end);
  if (~isempty(list) && ~iscell(list{1}) && strcmpi(list{1}, 'all'))
    reqStats = officialStats;
    reqStatsIdx = 1:length(officialStats);
    
  elseif (isempty(list) || (~iscell(list{1}) && strcmpi(list{1},'basic')))
    % Default list
    reqStats = {'Area'
                'Centroid'
                'BoundingBox'};
  else
    
    if (iscell(list{1}))
      list = list{1};
    end
    list = list(:);

    officialStatsL = lower(officialStats);
    
    reqStatsIdx = [];
    eid = sprintf('Images:%s:invalidMeasurement',mfilename);
    for k = 1:length(list)
      if (~ischar(list{k}))
        msg = sprintf('This measurement is not a string: "%d".', list{k});
        error(eid,'%s',msg);
      end
      
      idx = strmatch(lower(list{k}), officialStatsL);
      if (isempty(idx))
        msg = sprintf('Unknown measurement: "%s".', list{k});
        error(eid,'%s',msg);
        
      elseif (length(idx) > 1)
        msg = sprintf('Ambiguous measurement: "%s".', list{k});
        error(eid,'%s',msg);
        
      else
        reqStatsIdx = [reqStatsIdx; idx];
      end
    end
    
    reqStats = officialStats(reqStatsIdx);
  end
end

%%%
%%% NoNDSupport
%%%
function NoNDSupport(str)
%   Issue a warning message about lack of N-D support for a given
%   measurement or measurements.
  
  wid = sprintf('Images:%s:measurementNotForN-D',mfilename);

  if iscell(str)
    warn_str = sprintf('%s: %s ', ...
                       'These measurements are not supported if ndims(L) > 2.', ...
                       sprintf('%s ', str{:}));
  else
    warn_str = sprintf('%s: %s', ...
                       'This measurement is not supported if ndims(L) > 2.', ...
                       str);
  end

  warning(wid,'%s',warn_str);
end

%%%
%%% PreprocessRequestedStats
%%%
function requestedStats = PreprocessRequestedStats(requestedStats)
%   Remove any requested stats that are not supported for N-D input
%   and issue an appropriate warning.

  no_nd_measurements = {'MajorAxisLength'
                      'MinorAxisLength'
                      'Eccentricity'
                      'Orientation'
                      'ConvexHull'
                      'ConvexImage'
                      'ConvexArea'
                      'ConvexCentroid'
                      'ConvexAspect'
                      'EulerNumber'
                      'Extrema'
                      'EquivDiameter'
                      'Solidity'
                      'Extent'
                      'Perimeter'};

  bad_stats = find(ismember(requestedStats, no_nd_measurements));
  if ~isempty(bad_stats)
    NoNDSupport(requestedStats(bad_stats));
  end

  requestedStats(bad_stats) = [];
end