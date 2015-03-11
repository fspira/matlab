%%%%%% Plot angle distribution from orientaionJ

%%%% Convert angle from degree to rad
%%%% 


orientationData(:,1) = angleInRadians
rose(orientationData(:,1),orientationData(:,2))


tiffr