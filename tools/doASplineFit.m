function [splineFitOutline, tmp]= doASplineFit(segmentedOutline,tmp)

%%%%% if whole stack shall be fitted us a loop
lauf = 1;
tmp(:,:) = 0;
%img = S1;
%tmp = S1;

runningSplineA_x =  segmentedOutline(:,1);
runningSplineA_y =  segmentedOutline(:,2);

splineEnd = length(runningSplineA_x);
 
 knots(:,1:splineEnd) = [runningSplineA_x runningSplineA_y]';
%x = knots(1, :);
%y = knots(2, :);
%areaOfPolygon = polyarea(x,y);
finerSpacing = 1:0.01:splineEnd;
originalSpacing = 1:length(knots);

splineXY = spline(originalSpacing, knots, finerSpacing);

A = round(splineXY);
splineXY = A';



%imshow(tmp)
%hold on
tmpSpline = tmp;

X1 = splineXY(:,1);
Y1 = splineXY(:,2);



for subrun = 1:length(X1)
    tmpSpline(Y1(subrun),X1(subrun)) = 125;

end
imshow(tmpSpline,[])

tmp= tmpSpline;


%%%%%%%%% This section re-initializes the index of the BW image. This is
%%%%%%%%% important because the spline function uses a left to right
%%%%%%%%% indexing

tmp = bwmorph(tmp,'bridge',8);

B = bwboundaries(tmp,8);



 BLength = length(B);


[m n p] = size(B);

for subrun = 1:m
    
    BLength(subrun) = length(B{subrun});
    
end

B_Select = max(BLength);
B_Select = find(B_Select == BLength);
B = B{B_Select};
SplineXStore_1{lauf} = B(:,1);
SplineYStore_1{lauf} =  B(:,2);

splineFitOutline = [SplineXStore_1{lauf},SplineYStore_1{lauf}];

