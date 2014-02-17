function [yCat] = doInterpolateLineThroughPoints(xp, yp, tmp)

imshow(tmp)
%xp = [cNewMid1(lauf) cNewMid2(lauf)];
%yp = [rNewMid1(lauf) rNewMid2(lauf)];


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

%%%%% if number of points is too sparse, interpolate points
if length(yCat) < 200
    
    
     splineEnd = length(yCat);
 
    knots(:,1:splineEnd) = [yCat(:,1) yCat(:,2)]';
    finerSpacing = 1:0.01:splineEnd;
    originalSpacing = 1:splineEnd;

    splineXY = spline(originalSpacing, knots, finerSpacing);

    A = round(splineXY);
    splineXY = A';

    tmpSpline = tmp;

    X1 = splineXY(:,1);
    Y1 = splineXY(:,2);
    


    for subrun = 1:length(X1)
        tmpSpline(X1(subrun),Y1(subrun)) = 125;

    end
    tmp1 = tmpSpline;
    imshow(tmpSpline,[]);
    clear yCat
    yCat = [X1,Y1];
    
end