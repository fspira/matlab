%%% aniso struct importer


function [chromatinDistance_center,chromatinDistance_pole,aniso_vertical_center, aniso_vertical_pole,...
   ratio_vertical_center,ratio_vertical_pole, predictedTime_vertical,fileName_vertical, aniso_horizontal_center,...
    aniso_horizontal_pole,ratio_horizontal_center,ratio_horizontal_pole,...
     predictedTime_horizontal, fileName_horizontal]  = doImmportAnisData(s)
 


load('ratioAnisoParameters.mat')
f = fieldnames(s);


increment = 1;
  counterHorizontal = 1;
  counterVertical = 1;

        for lauf = 1:length(f)/2
            
            if s.(f{increment}).orientation.orientation > 45 && s.(f{increment}).orientation.orientation > 135

                    chromatinDistance_center(counterVertical)= s(1).(f{increment}).chromatin_distance.chromatin_distance
                    chromatinDistance_pole(counterVertical)= s(1).(f{increment+1}).chromatin_distance.chromatin_distance

                    aniso_vertical_center(counterVertical) = s(1).(f{increment}).averageRatioD_Vertical.averageRatioD_Vertical  
     
                    aniso_vertical_pole(counterVertical) = s(1).(f{increment+1}).averageRatioD_Vertical.averageRatioD_Vertical  
                  

                    ratio_vertical_center(counterVertical) = s(1).(f{increment}).averageRatio_Vertical.averageRatio_Vertical 
               

                    ratio_vertical_pole(counterVertical) = s(1).(f{increment+1}).averageRatio_Vertical.averageRatio_Vertical
                

                    predictedTime_vertical(counterVertical) = s.(f{increment+1}).predicted_time.predicted_time
                    fileName_vertical{counterVertical} = f{increment};

                    increment = increment+2
                    counterVertical = counterVertical+1
                    
            else
                
                  aniso_horizontal_center(counterHorizontal) = s(1).(f{increment}).averageRatioD_Horizontal.averageRatioD_Horizontal
                  aniso_horizontal_pole(counterHorizontal) = s(1).(f{increment+1}).averageRatioD_Horizontal.averageRatioD_Horizontal
                  ratio_horizontal_center(counterHorizontal) = s(1).(f{increment}).averageRatio_Horizontal.averageRatio_Horizontal
                  ratio_horizontal_pole(counterHorizontal) = s(1).(f{increment+1}).averageRatio_Horizontal.averageRatio_Horizontal
                  
                    predictedTime_horizontal(counterHorizontal) = s.(f{increment+1}).predicted_time.predicted_time

                   fileName_horizontal{counterHorizontal} = f{increment};
                    increment = increment+2  
                    counterHorizontal = counterHorizontal+1
         
                    
            end
        end

        test = 1
          
            