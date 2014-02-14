function tmp1 = doDrawLongAxis(ccc,BW2,flipMarker)

%setLength = 73; %%% Arbitrary value

    %plot(ccc.Extrema(:,1),ccc.Extrema(:,2),'x');
    %plot(ccc.Centroid(1),ccc.Centroid(:,2),'xr');

    %%%%%% plot endpoints of ellipse long axis. split chromatin at the long
    %%%%%% axis

    
   %imshow(BW2)
   %hold on
   % xp(1) = ccc(1).Centroid(1,1)+setLength;
   % yp(1) = ccc(1).Centroid(1,2)-(tan((ccc(1).Orientation)/380*2*pi)*+1);
    
    
   % dx = chromatin2(1) - chromatin1(1);
   % dy = chromatin2(2) - chromatin1(2);
    
   
   % plot(x,y,'db');
    
    
    %      xp(2) = ccc(1).Centroid(1,1)-setLength;
    %      yp(2) = ccc(1).Centroid(1,2)+(tan((ccc(1).Orientation)/360*2*pi)*+1);
    %  plot(xp(2),yp(2),'dr');
    
      
   % flipMarker determines whether ellipse was flipped. If so long and
   % short axis are flipped and intersection points have to be flipped as
   % well.
      
      
  % if flipMarker == 0
        [xp(1) yp(1)] = doSelectEllipsePoint(ccc, 0,1)
   
        [xp(2) yp(2)] = doSelectEllipsePoint(ccc, 180,1)
  % else
  %     [xp(1) yp(1)] = doSelectEllipsePoint(ccc, -90,1)
   
  %     [xp(2) yp(2)] = doSelectEllipsePoint(ccc, 90,1)
       
  % end
       
   imshow(BW2)
   hold on
   
    plot(xp(2),yp(2),'ro')
    plot(xp(1),yp(1),'bo')
   
   
    a = (yp(2)-yp(1)) / (xp(2)-xp(1));
    b = yp(1)-a*xp(1);
    imshow(BW2);
  %%%% Determine size of the active window
    xlims = xlim(gca);
    ylims = ylim(gca);

%%%% Calculate long axis through the entire active window
y = xlims*a+b;
line( xlims, y );

[mmm nnn ppp] = size(BW2);

%%%% Calculate individual points of the long axis
    for subrun = 1:mmm
        
    
        y1(subrun) = subrun*a+b;

    end

    shapeStart = ceil(xlims(1));
    shapeStop = floor(xlims(2));
    
y1shape = y1(shapeStart:shapeStop)
    

%%%% format the entries, round entries and eliminate negative values

yCat = [y1' [1:mmm]'];

yCatTmp = yCat(:,1:2) < xlims(2);
yCat = yCat .* yCatTmp;
yCat = [round(yCat(:,1)),round(yCat(:,2))];

yCatNorm = (yCat > 0);
yCat = yCatNorm .* yCat;
catZero = find(yCat(:,1) >0);

yCat = yCat(catZero,:);
tmp1 = BW2;

%%%%%% plot long axis into the images


if length(yCat) <= 2
    
    
     close all
     imshow(BW2)
     hold on
     plot(xp(1),yp(1),'xr')
     plot( xp(2),yp(2),'xr')
     hLine = imline(gca,[xp(1),xp(2)], [yp(1),yp(2)]); 
     binaryImage1 = hLine.createMask();
     
     tmp1 = uint8(BW2) + uint8( binaryImage1);
     
else
    
    for subrun = 1:length(yCat)

      tmp1((round(yCat(subrun,1))),round(yCat(subrun,2))) = 125;

    end


    
     
   
    
end


%%%%%% it may happen that the line is not a line, but unconnected dots. In
%%%%%% this case a spline is interpolated to connect the pixels

if length(yCat) < 200
    
    
     splineEnd = length(yCat);
 
    knots(:,1:splineEnd) = [yCat(:,1) yCat(:,2)]';
    finerSpacing = 1:0.01:splineEnd;
    originalSpacing = 1:splineEnd;

    splineXY = spline(originalSpacing, knots, finerSpacing);

    A = round(splineXY);
    splineXY = A';

    tmpSpline = BW2;

    X1 = splineXY(:,1);
    Y1 = splineXY(:,2);
    


    for subrun = 1:length(X1)
        tmpSpline(X1(subrun),Y1(subrun)) = 125;

    end
    tmp1 = tmpSpline;
    imshow(tmpSpline);
end