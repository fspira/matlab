%%%%%%%% Function to detect curvature of the membrane

%[] = doDetermineCurvature()
[mm nn pp] = size(cNew1Ref);

for lauf = 1:nn
    close all
    
runningSplineA_x = cNew1Ref{lauf};
runningSplineA_y = rNew1Ref{lauf};

runningSplineB_x = cNew2Ref{lauf};
runningSplineB_y = rNew2Ref{lauf};

%%%% identifies the central region of the spline
center =   redMax1(lauf)
tmp1 = greenStack(:,:,1);
tmp1(:,:) = 0;

tmp1(rNewMid1(lauf), cNewMid1(lauf)) =255;

tmp = greenStack(:,:,1);
tmp(:,:) = 0;
imshow(tmp1)
hold on
line(runningSplineA_x,runningSplineA_y)

plot(cNewMid1(lauf),rNewMid1(lauf),'x')
plot(cNewMid2(lauf),rNewMid2(lauf),'x')
pause(0.1)



%%%% Calculate intersection points between ingressing furrow
%%%% First exrapolate line between marked center points of the ingressing
%%%% furrow
imshow(tmp)
xp = [cNewMid1(lauf) cNewMid2(lauf)];
yp = [rNewMid1(lauf) rNewMid2(lauf)];


    a = (yp(2)-yp(1)) / (xp(2)-xp(1));
    b = yp(1)-a*xp(1);

  %%%% Determine size of the active window
  xlims = xlim(gca);
    ylims = ylim(gca);

%%%% Calculate line
y = xlims*a+b;

hold on
plot(xp(1),yp(1),'x')
plot(xp(2),yp(2),'x')
line( xlims, y );

for subrun = 1:512
    
    y1(subrun) = subrun*a+b;

end

y1shape = y1(round(xlims(1))):round(xlims(2))
    
%%%% store values of the line connecting the center points of the
%%%% ingressing furrow
%%%% These lines eliminate negative values
%yCat = [y1' [round(xlims(1)):round(xlims(2))]']

yCat = [y1' [1:512]']

yCatTmp = yCat(:,1:2) < xlims(2)
yCat = yCat .* yCatTmp;
yCat = [round(yCat(:,1)),round(yCat(:,2))];

yCatNorm = (yCat > 0);
yCat = yCatNorm .* yCat;
catZero = find(yCat(:,1) >0);

yCat = yCat(catZero,:);
%catNew = yCat(catZero)

%flag =  isempty(catZero)
%if flag == 0 

 %   yCat= yCat(1:catZero(1)-1,1:2);
 %   else
%end

%%%% To obtain values for intersecting points geometrically draw the shape
%%%% and the line connecting two points assigning values of 125 to each
%%%% shape. After addition of both shapes pixel values having a value
%%%% greate than 125 are intersecting pixels.
%%%% Calculate flank 1
tmp(:,:) = 0;
tmp1 = tmp;
tmp1(:,:) = 0;
knots = [];

%%%%% This section fits a spline to the individual data points
 splineEnd = length(runningSplineA_x);
 
 knots(:,1:splineEnd) = [runningSplineA_x runningSplineA_y]';
%x = knots(1, :);
%y = knots(2, :);
%areaOfPolygon = polyarea(x,y);
finerSpacing = 1:0.01:splineEnd;
originalSpacing = 1:splineEnd;

splineXY = spline(originalSpacing, knots, finerSpacing);

A = round(splineXY);
splineXY = A';

%%%% The following tow lines will reduce the number of knots within the
%%%% spline. 

%%[u,id1,id2] = unique(A(:,1:2),'rows');

%%splineXY = u';


%plot(knots(1, :), knots(2, :), 'ro',...
%splineXY(1, :), splineXY(2, :), 'b+-', 'LineWidth', 2, 'MarkerSize',16);

imshow(tmp)
hold on
tmpSpline = tmp;

X1 = splineXY(:,1);
Y1 = splineXY(:,2);



for subrun = 1:length(X1)
    tmpSpline(Y1(subrun),X1(subrun)) = 125;

end
imshow(tmpSpline)

tmp= tmpSpline;

%%%%%%%%% This section re-initializes the index of the BW image. This is
%%%%%%%%% important because the spline function uses a left to right
%%%%%%%%% indexing


B = bwboundaries(tmp,8);

length(B);
clear BLength

[m n p] = size(B)

for subrun = 1:m
    
    BLength(subrun) = length(B{subrun});
    
end

B_Select = max(BLength);
B_Select = find(B_Select == BLength);
B = B{B_Select};
SplineXStore_1{lauf} = B(:,1);
SplineYStore_1{lauf} =  B(:,2);

tmp(:,:) = 0;

for subrun  = 1:length(B(:,1))
    tmp(B(subrun,1),B(subrun,2)) = 125;
end

clear B

for subrun = 1:length(yCat)

    tmp1((round(yCat(subrun,1))),round(yCat(subrun,2))) = 125;
    
end



tmp1 = bwmorph(tmp1,'bridge',8);
tmp1 = im2bw(tmp1);
tmp1 = tmp1.*125;

tmp2 = tmp + tmp1;
imshow(tmp2)

[rIntersect1,cIntersect1] = find(tmp2 > 130);

flag = isempty(rIntersect1)

if flag == 1
    
    se1 = fspecial('gaussian',2,1);
    tmp1 = imfilter(tmp1,se1);
     tmp1 = im2bw(tmp1);
    tmp1 = tmp1.*125;
    tmp2 = tmp + tmp1;
    
    
[rIntersect1,cIntersect1] = find(tmp2 > 130);

else
end



%%%% use central value of line
[mi ni pi] = size(rIntersect1);
if mi >1
    rIntersect1 = rIntersect1(round(mi/2));
    cIntersect1 = cIntersect1(round(mi/2));
    else
end
close all
%%%% calculate flank 2 %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp(:,:) = 0;


%%%%% This section fits a spline to the individual data points
 splineEnd = length(runningSplineB_x);
 
 knots(:,1:splineEnd) = [runningSplineB_x runningSplineB_y]';
%x = knots(1, :);
%y = knots(2, :);
%areaOfPolygon = polyarea(x,y);
finerSpacing = 1:0.01:splineEnd;
originalSpacing = 1:splineEnd;

splineXY = spline(originalSpacing, knots, finerSpacing);
%plot(knots(1, :), knots(2, :), 'ro',...
%splineXY(1, :), splineXY(2, :), 'b+-', 'LineWidth', 2, 'MarkerSize',16);


A = round(splineXY)';
splineXY= A;
%%%% The following tow lines will reduce the number of knots within the
%%%% spline. 
%[u,id1,id2] = unique(A(:,1:2),'rows');

%splineXY = u';


imshow(tmp)
hold on
tmpSpline = tmp;


X1 = splineXY(:,1);
Y1 = splineXY(:,2);



for subrun = 1:length(X1)
    tmpSpline(Y1(subrun),X1(subrun)) = 125;

end
imshow(tmpSpline)

tmp= tmpSpline;

%%%%%%%%% This section re-initializes the index of the BW image. This is
%%%%%%%%% important because the spline function uses a left to right
%%%%%%%%% indexing

B = bwboundaries(tmp,8);


clear BLength

[m n p] = size(B)
for subrun = 1:m
    
    BLength(subrun) = length(B{subrun});
    
end

B_Select = max(BLength);
B_Select = find(B_Select == BLength);
B = B{B_Select};

SplineXStore_2{lauf} =  B(:,1);
SplineYStore_2{lauf} =  B(:,2);

tmp(:,:) = 0;

for subrun  = 1:length(B(:,1))
    tmp(B(subrun,1),B(subrun,2)) = 125;
end

clear B
%se = strel('ball',3,3,8); % radius of 4
%BW1 = imclose(tmp1,se);



imshow(tmp,[])



tmp2 = tmp + tmp1;
imshow(tmp2)

[rIntersect2,cIntersect2] = find(tmp2 > 135);

flag = isempty(rIntersect2)

if flag == 1
    
    se1 = fspecial('gaussian',2,1);
    tmp1 = imfilter(tmp1,se1);
    tmp1 = im2bw(tmp1);
    tmp1 = tmp1.*125;
    tmp2 = tmp + tmp1;
    
[rIntersect2,cIntersect2] = find(tmp2 > 135);

else
end

[mi ni pi] = size(rIntersect2);
if mi >1
    rIntersect2 = rIntersect2(round(mi/2));
    cIntersect2 = cIntersect2(round(mi/2));
    else
end
close all

cInterStore1(lauf) = cIntersect1(1);
rInterStore1(lauf) = rIntersect1(1);

cInterStore2(lauf) = cIntersect2(1);
rInterStore2(lauf) = rIntersect2(1);

close all

end

for lauf = 1:nn
tmp1(:,:)=0;
imshow(tmp1)
hold on
plot(cNew1Ref{lauf},rNew1Ref{lauf})
plot(cNew2Ref{lauf},rNew2Ref{lauf})
plot(cInterStore1(lauf),rInterStore1(lauf),'xr')
plot(cInterStore2(lauf),rInterStore2(lauf),'xr')
plot(cNewMid1(lauf),rNewMid1(lauf),'xy')
plot(cNewMid2(lauf),rNewMid2(lauf),'xy')
pause(0.1)
lauf
end

for lauf = 1:nn
%%%% set sliding window for ellipse fitting
slidingWindow = 25;



%%%% crop the spline to a size as given by the sliding window

InterSplineA_x = SplineXStore_1{lauf};
InterSplineA_y = SplineYStore_1{lauf};



%%%%%%%% Find the center pixel


centerA = find(cInterStore1(lauf) == InterSplineA_x & rInterStore1(lauf) == InterSplineA_y)


centerA = find(cInterStore1(lauf) == test(:,2) & rInterStore1(lauf) == test(:,1))


splineA_x = InterSplineA_x(centerA-(slidingWindow)/2:centerA+(slidingWindow)/2)
splineA_y = InterSplineA_y(centerA-(slidingWindow)/2:centerA+(slidingWindow)/2)

tmp = greenStack(:,:,1);
tmp(:,:) = 0;


%tmp(rNew1Ref{lauf}, cNew1Ref{lauf})=255;
%spline(rNew1Ref{lauf}', cNew1Ref{lauf}');

%plot(splineA_x,splineA_y,'-')

for sublauf = 1:length(InterSplineA_y)

    tmp(test(sublauf,1),test(sublauf,2)) = 255;
    imshow(tmp)
    pause(0.1)
end




%line(splineA_x, splineA_y)

imBinary = im2bw(tmp);
imshow(imBinary)
cc = regionprops(imBinary,'all');
  close all
%B = bwboundaries(imBinary,8); 
%  boundary = B{1};


  
k = 1;
phi = linspace(0,2*pi,50);
cosphi = cos(phi);
sinphi = sin(phi);
    
    xbar = cc(k).Centroid(1);
    ybar = cc(k).Centroid(2);

    a = cc(k).MajorAxisLength/2;
    b = cc(k).MinorAxisLength/2;

    theta = pi*(cc(k).Orientation/180);
    
    R = [ cos(theta)   sin(theta)
         -sin(theta)   cos(theta)];

    xy = [a*cosphi; b*sinphi];
    xy = R*xy;

    x = xy(1,:) + xbar;
    y = xy(2,:) + ybar;

    imshow(imBinary,[])
    hold on
    plot(x,y,'r','LineWidth',2);
  
    x_EllipseStore{lauf} = x;
    y_EllipseStore{lauf} = y;
    
    splineA_xStore{lauf} = splineA_x;
    splineA_yStore{lauf} = splineA_y;
    
   % imshow(img,[])
   % hold on
   % plot(x_EllipseStore{lauf},y_EllipseStore{lauf},'r','LineWidth',2);
    
    Curvature(lauf) = cc(k).MajorAxisLength / cc(k).MinorAxisLength
close all
end
    

cc.Centroid(1,2)-(tan((cc.Orientation+90)/360*2*pi)*-1)




L2 = greenStack(:,:,1);
L2(:,:) = 0;
for lauf = 1:25
    
    
  imshow(imBinary,[])
    hold on
 line(splineA_xStore{lauf}, splineA_yStore{lauf})
  % plot(flank1Store{:lauf},flank1Store{lauf}, 'b', 'LineWidth', 2)
   
  plot(x_EllipseStore{lauf},y_EllipseStore{lauf},'r','LineWidth',2);
  
    plot( cNew1Ref{lauf}, rNew1Ref{lauf}, 'g', 'LineWidth', 2)
    pause(0.1)
     %title(['MMLCII-B Furrow Region:' tifFilename],'FontSize', 20);
   % mov(lauf) = getframe(gcf);
     hold off
end


close(gcf)

%# save as AVI file, and open it using system video player
movie2avi(mov, 'Ellipse.avi', 'compression','None', 'fps',10);
 
cNewMid1(lauf), rNewMid1(lauf)

spline(cNew1Ref{lauf}', rNew1Ref{lauf}')
