function tmp1 = doConnectTwoPoints(xp,yp,BW2) 

     close all
     imshow(BW2)
     hold on
     plot(xp(1),yp(1),'xr')
     plot( xp(2),yp(2),'xr')
     hLine = imline(gca,[xp(1),xp(2)], [yp(1),yp(2)]); 
     binaryImage1 = hLine.createMask();
     
     tmp1 = uint8(BW2) + uint8( binaryImage1);
     close all