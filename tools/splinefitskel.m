function [skelfit, skeldist, rowI, colI, rowseed, colseed, skelI]=splinefitskel(skel,P,K,varargin)
% SPLINEFITSKEL Compute a spline approximation to a 2D skeleton.
%   [skelfit skeldist rowI colI]=SPLINEFITSKEL(skel,P,K,<WeightMat>,<refI>,<PixSize>,<DebugPlot>)
%   fits a 2-dimensional spline of single variable to describe a spatial
%   curve. The variable is the contour length of the curve measured with
%   respect to a reference point and the 2-dimensions of the spline are the
%   row indices and column indices (estimated to a fraction of the pixel).
%
%   INPUTS:
%   skel: Binary image representing the skeleton to be fit.
%   P: Pieces (of equal length) in the spline fit.
%   K: Order of B-spline used for approximation.
%
%   OPTIONAL INPUTS:
%   WeightMat: Matrix of the size same as the skel that describes relative
%   weights to be assigned to the points on skel for fitting the curve. If
%   all points on skel carry same weight, use double(skel) as WeightMat.
%   The points off-the-skeleton in WeightMat do not matter.
%   REFPOINT: 2-element vector [R C] with row index (R) and column index
%   (C) of the reference point from which the distance on the skeleton is
%   measured.
%   PixSize: Do the fit assuming certain pixel size. Useful for fitting
%   data in the units in which they are acquired.
%   DebugPlot: Number of figure on which the debug info (overlay of the
%   weight matrix/skeleton and the fit are plotted).
%
%   OUTPUTS:
%   skelfit: Structure that describes the fit.
%   skeldist: Equispaced vector that represents the contour distant along
%   the skeleton with respect to the reference point.
%   rowI:   Row-indices (can be fraction of pixel) of the smooth curve at
%   distances containted in skeldist.
%   colI:   Column-indices (can be fraction of pixel) of the smooth curve
%   at distances contained in skeldist.
% 
%   See included script that analyzes experimental microscopy data
%   for usage example.
%   
%   VERSION HISTORY:
%   Ver 1: April 20, 2012 Shalin Mehta (Marine Biological Laboratory). 
%   shalin.mehta@gmail.com
%       First implementation, sanity-check of optional arguments postponed
%       to next version.

%%% Set-default values for optional arguments, no tests done at this
%%% point.

args.weightMat=double(skel);
args.refpoint=NaN;
args.PixSize=1;
args.DistFromRef=true;
args.DebugDisp=0;
args.dosplinefit=true;

args=parsepropval(args,varargin{:});

weightMat=args.weightMat;
refpoint=args.refpoint;
PixSize=args.PixSize;
DistFromRef=args.DistFromRef;
DebugDisp=args.DebugDisp;


%%% Find the endpoint of the skeleton closest to the reference.
%--------------
EPskel=bwmorph(skel,'endpoints',4);

[EProw, EPcol]=find(EPskel);
% if(length(EProw)~=2)
%     error('The skeleton does not describe a single spatial curve.');
% end
if(isnan(refpoint) || any(refpoint == 0))% If reference point is not valid, then assume that the left most point is the reference.
   refpoint(1)=EProw(1);
   refpoint(2)=EPcol(1);
end

EPdist=sqrt( (EProw-refpoint(1)).^2 + (EPcol-refpoint(2)).^2 );
[refdist, minidx]=min(EPdist);
rowseed=EProw(minidx); colseed=EPcol(minidx);

%%% Obtain pixels sorted according to the distance from the reference
%%% point.
%--------------

% Map of countour distance with respect to the seed point: geodesic
% distance transform constrains the path that can be used to compute the
% distance to the skeleton. The distances on pixels that are not part of
% the skeleton are set to NaN.

skeldistmap=bwdistgeodesic(skel,colseed,rowseed,'quasi-euclidean');

if(DistFromRef)
    skeldistmap=skeldistmap+refdist; % Measure distance from the reference point.
end

% Sort the distances and find corresponding indices of the skeletal pixels. 
[skeldist, skelI]=sort(skeldistmap(:));
%Find end of skeleton as the first sorted pixel at which distance is set to NaN.
EndofSkel=find(isnan(skeldist),1,'first'); 
skeldist=skeldist(1:EndofSkel-1);
skelI=skelI(1:EndofSkel-1);
[rowI, colI]=ind2sub(size(skel),skelI);

 
% Convert all distances in the units of the pixel (e.g., micron) starting
% with zero.
skeldist=skeldist*PixSize;
rowI=rowI*PixSize; 
colI=colI*PixSize;


if(args.dosplinefit)

    %%% Think of row indices and column indices of the skeleton as a vector
    %%% function of the skeleton's contour length and find a vector-spline of
    %%% specified K that approximates it the best.
    %--------------

    %Read the weights from the weight map.
    weights=weightMat(skelI); 
    % Obtain the knot sequence for least-squares approximation, the sequence
    % will be refined after the first approximation.
    uniform=linspace(min(skeldist),max(skeldist),2+P); 

    %NOTE: If multiplicity of internal knots is increased, the algorithm allows
    %for higher derivatives to be non-zero around knots, thereby allowing
    %discontinuities at the knots. This is the reason, why the function drops
    %to zero outside of defined support.
    knots=augknt(uniform,K);

    % Do the spline approximation.
    skelfit=spap2(knots,K,skeldist,[rowI colI]',weights);

    %% newknt useful for distributing errors evenly across the segments.
    % But not yet implemented for vector-valued functions. 
    % skelfit=spap2(newknt(skelfit),K,skeldist,[rIpix cIpix]',weights);

    % Obtain new row, and column positions that sample the
    % contour of the spline fitted to the skeleton. These are output.
    skeldist=linspace(min(skeldist),max(skeldist),length(skeldist));
    rowcolI=fnval(skelfit,skeldist);
    rowI=rowcolI(1,:)'; colI=rowcolI(2,:)';
else
    skelfit=NaN;
end

%%% Check if the spline 'fits snuggly' to higher weight
if(DebugDisp)
    figure(DebugDisp);
    xaxis=PixSize*(1:size(skel,2)-1);
    yaxis=PixSize*(1:size(skel,1)-1);
    imagesc(xaxis,yaxis,weightMat,[0 150]); axis equal; axis xy; colormap gray;
    hold on;
    plot(colI,rowI,'x-r','LineWidth',1);
    plot(colI,rowI,'x-b','LineWidth',1);
    hold off;
end

end