function [ ratioInterp_A ] = interpolatePixelPerpendicularToLine(inputImage, Fun_coords, lineLength, angleFunc)
    

        x(2) =  Fun_coords(1,1) + lineLength * cosd(angleFunc(1)-90);
        y(2) =  Fun_coords(2,1) + lineLength * sind(angleFunc(1)-90);


           lineLength = -3;

        x(3) =  Fun_coords(1,1) + lineLength * cosd(angleFunc(1)-90);
        y(3) =  Fun_coords(2,1) + lineLength * sind(angleFunc(1)-90);


        
        %%%%
        %%%% Interpolated line length
        %%%%
       
        InterpolatedPixel_I = sqrt((Fun_coords(1,1) - x(2))^2+ (Fun_coords(2,1) - y(2))^2);
        InterpolatedPixel_II = sqrt((Fun_coords(1,1) - x(3))^2+ (Fun_coords(2,1) - y(3))^2);
        
        
        %%%%
        %%%% interpolates the intensities along the line
        %%%%
      Vq_I = interp2(inputImage(:,:,1),linspace(Fun_coords(2,1),y(3),InterpolatedPixel_I),linspace(Fun_coords(1,1),x(3),InterpolatedPixel_I));
      Vq_II = interp2(inputImage(:,:,1),linspace(Fun_coords(2,1),y(2),InterpolatedPixel_II),linspace(Fun_coords(1,1),x(2),InterpolatedPixel_I));
       
      
      ratioInterp_A(1) = (mean(Vq_I) + mean(Vq_II))/2;
       
        

end

