function [ ratioInterp_Store ] = interpolatePixelPerpendicularToLine(inputImage, Fun_coords, lineLength, angleFunc)
    
for subrun = 1:length(Fun_coords)
    

        lineLengthP = lineLength;
        lineLengthN = 0 -  lineLength;


        x(2) =  Fun_coords(1,subrun) + lineLengthP * cosd(angleFunc(subrun)-90);
        y(2) =  Fun_coords(2,subrun) + lineLengthP * sind(angleFunc(subrun)-90);



        x(3) =  Fun_coords(1,subrun) + lineLengthN * cosd(angleFunc(subrun)-90);
        y(3) =  Fun_coords(2,subrun) + lineLengthN * sind(angleFunc(subrun)-90);


        
        %%%%
        %%%% Interpolated line length
        %%%%
       
        InterpolatedPixel_I = sqrt((Fun_coords(1,subrun) - x(2))^2+ (Fun_coords(2,subrun) - y(2))^2);
        InterpolatedPixel_II = sqrt((Fun_coords(1,subrun) - x(3))^2+ (Fun_coords(2,subrun) - y(3))^2);
        
        
        %%%%
        %%%% interpolates the intensities along the line
        %%%%
      Vq_I = interp2(inputImage(:,:,1),linspace(Fun_coords(1,subrun),x(3),InterpolatedPixel_I),linspace(Fun_coords(2,subrun),y(3),InterpolatedPixel_I));
      Vq_II = interp2(inputImage(:,:,1),linspace(Fun_coords(1,subrun),x(2),InterpolatedPixel_II),linspace(Fun_coords(2,subrun),y(2),InterpolatedPixel_II));
       
      
      ratioInterp_A(1) = (mean(Vq_I) + mean(Vq_II))/2;
      
      ratioInterp_Store(subrun) = ratioInterp_A;
     

% 
% imshow(inputImage,[])
% hold on
% plot(linspace(Fun_coords(1,subrun),x(3),InterpolatedPixel_I),linspace(Fun_coords(2,subrun),y(3),InterpolatedPixel_I),'r')
%  plot(   linspace(Fun_coords(1,subrun),x(2),InterpolatedPixel_II),linspace(Fun_coords(2,subrun),y(2),InterpolatedPixel_II),'g')
%       
      
      

subrun
       
end    

end

