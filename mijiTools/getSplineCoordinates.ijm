run("Fit Spline", "straighten"); 
  getSelectionCoordinates(x, y); 
  x2=0; y2=0; distance=0; 
  for (i=0; i<x.length; i++) { 
      if (i>0) { 
         dx = x[i] - x[i-1]; 
         dy = y[i] - y[i-1]; 
         distance = sqrt(dx*dx+dy*dy); 
      } 
      print(i, x[i], y[i], distance); 
  } 