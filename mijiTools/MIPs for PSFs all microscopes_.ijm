S=fromCharCode(92);

//_______________________Format date______________________________
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
Dayadjust = "";
if (dayOfMonth<10) {Dayadjust = "0";}
month=month+1;
Monthadjust = "";
if (month<10) {Monthadjust = "0";}
today = ""+year+"-"+Monthadjust+month+"-"+Dayadjust+dayOfMonth;

//_______________________Test location of user FMI vs. not FMI___________

pathTest=S+S+"imagestore"+S+"FAIM"+S+"Maintenance_All microscopes"+S+"TestFMI_DoNotMove.txt";
Loc = "";
if (File.exists(pathTest) == 1) Loc="FMI";

//_______________________Retrieve last saved Info________________________

	microscope="microscope";
	MA=100;
	NA=1.4;
	xyVoxel=70;
	zVoxel=200;
	date=today;
	ImageName = getInfo("image.filename");

Directory = getDirectory("imagej");
path = Directory+"LogFileImageJ.txt";
if (File.exists(path)==1) {
    Info = File.openAsString(path);
    if (Info!="") {
        firstK=substring(Info, 0, 1);
        if (firstK!="~") {
            index1=indexOf(Info, ";")+1;
            index2=indexOf(Info, ";", index1)+1;
            index3=indexOf(Info, ";", index2)+1;
            index4=indexOf(Info, ";", index3)+1;
            index5=indexOf(Info, ";", index4)+1;

            microscope=substring(Info, 0, index1-1);
            MA=substring(Info, index1, index2-1);
            NA=substring(Info, index2, index3-1);
            xyVoxel=substring(Info, index3, index4-1);
            zVoxel=substring(Info, index4, index5-1);
            date=substring(Info, index5, index5+10);
        }
    }
}

// _______________________Get info on the setup_______________________
Unit = newArray("uM", "nm");
getVoxelSize(xyVoxel_Info, height, zVoxel_Info, unit_Info);

//if (xyVoxel!=0 && xyVoxel_Info!=1) 
xyVoxel=xyVoxel_Info;
//if (zVoxel!=0 && zVoxel_Info!=1) 
zVoxel=zVoxel_Info;

Dialog.create("Image information");
  Dialog.addString("Microscope Name:", microscope);
  Dialog.addNumber("Magnification:", MA);
  Dialog.addString("NA:", NA);
  Dialog.addNumber("Pixel size :", xyVoxel);
  Dialog.addChoice("Unit", Unit, Unit[0]);
  Dialog.addNumber("Distance between stacks :", zVoxel)
  Dialog.addChoice("Unit", Unit, Unit[0]);
  Dialog.addString("Date:", date);
Dialog.show();
  microscope = Dialog.getString();
  MA = Dialog.getNumber();
  NA = Dialog.getString();
  xyVoxel = Dialog.getNumber();
  UnitStackXY = Dialog.getChoice();
  zVoxel=Dialog.getNumber();
  UnitStackZ = Dialog.getChoice();
  date = Dialog.getString();

if (UnitStackXY==Unit[0]) xyVoxel=xyVoxel*1000;
if (UnitStackZ==Unit[0]) zVoxel=zVoxel*1000;


//________________________Save Info on setup____________________________

Info = microscope+";"+MA+";"+NA+";"+xyVoxel+";"+zVoxel+";"+date;
File.saveString(Info, path);

//________________________Change image properties________________________

ImageName = getInfo("image.filename");
if (Loc == "FMI") ImageName = "";
MainName=ImageName+"_"+date+"_"+microscope+"_"+MA+"x_"+NA;
setVoxelSize(xyVoxel, xyVoxel, zVoxel, "nm");
rename("Stack");

// _______________Select the slice with highest intensity and crop_________________

selectWindow("Stack");
run("Duplicate...", "title=Stack duplicate range=1-100");
run("32-bit");
run("Z Project...", "start=1 stop=100 projection=[Average Intensity]");
Max=0;
width=getWidth();
height=getHeight();
for (i=1; i<width+1; i++) {
	for (j=1; j<height+1; j++) {
		PixelValue=getPixel(i,j);
		if (PixelValue>Max) {
			Max=PixelValue;
			x2=i; y2=j;
		}}}
close(); close();
Maxstack=0; Minstack=65500; OptSlice=0; Max=0;
for (i=1; i<nSlices+1; i++) {
	setSlice(i);
	getStatistics(area, mean, min, max, std, histogram);
	PixelValue=getPixel(x2,y2);
	if (PixelValue>Max) {
		Max=PixelValue;
		Maxstack=max;
		Minstack=min;
		OptSlice=i;
		}}

selectWindow("Stack");
ROIsize=round(15000/xyVoxel);
halfROIsize = round(ROIsize/2);
makeRectangle(x2-halfROIsize, y2-halfROIsize, ROIsize, ROIsize);
run("Crop");

	ROIsizeBG = round(ROIsize/10);
	makeRectangle(ROIsizeBG, ROIsizeBG, ROIsizeBG, ROIsizeBG);
	getStatistics(area, mean, min, max, std, histogram);
	run("Select None");
	run("Subtract...", "stack value="+mean);

run("Duplicate...", "title=Stack duplicate range=1-100");
run("32-bit");
run("Z Project...", "start=1 stop=100 projection=[Average Intensity]");

Max=0;
width=getWidth();
height=getHeight();

for (i=1; i<width+1; i++) {
	for (j=1; j<height+1; j++) {
		PixelValue=getPixel(i,j);
		if (PixelValue>Max) {
			Max=PixelValue;
			x2=i; y2=j;
		}}}
close(); close();

Maxstack=0; Minstack=65500; OptSlice=0; Max=0;
for (i=1; i<nSlices+1; i++) {
	setSlice(i);
	getStatistics(area, mean, min, max, std, histogram);
	PixelValue=getPixel(x2,y2);
	if (PixelValue>Max) {
		Max=PixelValue;
		Maxstack=max;
		Minstack=min;
		OptSlice=i;
		}}

//__________Redimension stack and set OptSlice to plane 50_________

while (OptSlice+50>nSlices) {
	setSlice(nSlices);
	run("Add Slice");
	}
while (OptSlice+50<nSlices) {
	setSlice(nSlices);
	run("Delete Slice");
	}
while (nSlices>100) {
	setSlice(1);
	run("Delete Slice");
	}
while (nSlices<100) {
	setSlice(1);
	run("Add Slice");
	}

OptSlice=50;

// _______________________Projections_______________________

selectWindow("Stack");
setSlice(OptSlice);
makeLine(0, y2, ROIsize, y2);
run("Reslice [/]...", "input="+zVoxel+" output="+zVoxel+" slice=1");
rename ("xProj");
	ROIsizeBG = round(ROIsize/10);
	makeRectangle(ROIsizeBG, ROIsizeBG, ROIsizeBG, ROIsizeBG);
	getStatistics(area, mean, min, max, std, histogram);
	run("Select None");
	run("Subtract...", "value="+min);
H=getHeight();

selectWindow("Stack");
setSlice(OptSlice);
makeLine(x2, 0, x2, ROIsize);
run("Reslice [/]...", "input="+zVoxel+" output="+zVoxel+" slice=1 rotate");
rename ("yProj");
	ROIsizeBG = round(ROIsize/10);
	makeRectangle(ROIsizeBG, ROIsizeBG, ROIsizeBG, ROIsizeBG);
	getStatistics(area, mean, min, max, std, histogram);
	run("Select None");
	run("Subtract...", "value="+min);

selectWindow("Stack");
run("Z Project...", "start=1 stop=100"+" projection=[Max Intensity]");
	ROIsizeBG = round(ROIsize/10);
	makeRectangle(ROIsizeBG, ROIsizeBG, ROIsizeBG, ROIsizeBG);
	getStatistics(area, mean, min, max, std, histogram);
	run("Select None");
	run("Subtract...", "value="+min);
rename("Project");

ProjectWidth=ROIsize+H;
run("Canvas Size...", "width="+ProjectWidth+" height="+ProjectWidth+" position=Top-Left zero");
selectWindow("yProj");
run("Canvas Size...", "width="+ProjectWidth+" height="+ProjectWidth+" position=Top-Right zero");
selectWindow("xProj");
run("Canvas Size...", "width="+ProjectWidth+" height="+ProjectWidth+" position=Bottom-Left zero");
run("Image Calculator...", "image1=Project operation=Add image2=xProj");
selectWindow("Project");
run("Image Calculator...", "image1=Project operation=Add image2=yProj");
run("Size...", "width=550 height=550 constrain interpolation");

selectWindow("xProj");
close();
selectWindow("yProj");
close();

selectWindow("Project");
run("32-bit");
run("Square Root");
getStatistics(area, mean, min, max, std, histogram);
setMinAndMax(mean, max);
run("Invert");
run("LUTforPSFs");
run("8-bit");
run("RGB Color");


// _______________________FWHM axial_______________________ ;

selectWindow("Stack");

zProfileX = newArray(nSlices); zProfileY = newArray(nSlices);
for (i=0; i<nSlices; i++) {
	setSlice(i+1);
	zProfileX[i] = i;
	zProfileY[i] = getPixel(x2, y2);
	}

Fit.doFit("Gaussian", zProfileX, zProfileY); a=Fit.p(0); b=Fit.p(1); c=Fit.p(2); d=Fit.p(3);

// _______Plot

Amplitude = 40;
XplotLatReal = newArray(Amplitude); YplotLatReal = newArray(Amplitude);
MaxGraph = 0;
for (i=0; i<Amplitude; i++) {
	XplotLatReal[i] = (i-Amplitude/2)*zVoxel;
	YplotLatReal[i] = zProfileY[OptSlice-Amplitude/2+i];
	if (YplotLatReal[i]>=MaxGraph) {MaxGraph = YplotLatReal[i];}
	}
XplotLatFit = newArray(Amplitude*4);
YplotLatFit = newArray(Amplitude*4);
Ymin=66000;
Ymax=0;
for (i=0; i<Amplitude*4; i++) {
	XplotLatFit[i] = (i/4-Amplitude/2)*zVoxel;
	X = OptSlice-Amplitude/2+i/4;
	YplotLatFit[i] =  a + (b-a)*exp(-(X-c)*(X-c)/(2*d*d));
	if (YplotLatFit[i]>=MaxGraph) {MaxGraph = YplotLatFit[i];}
	if (Ymin>YplotLatFit[i]) {Ymin=YplotLatFit[i];}
	if (Ymax<YplotLatFit[i]) {Ymax=YplotLatFit[i];}
	}

HM=(Ymax-Ymin)/2; k=-2*d*d*log((HM-a)/(b-a)); FWHMa=2*zVoxel*sqrt(k);

Plot.create("FWHM axial", "Z", "Intensity", XplotLatFit, YplotLatFit);
Plot.setLimits(-4000, 4000, 0, MaxGraph*1.1);
Plot.add("circles", XplotLatReal, YplotLatReal);
Text="FWHM axial ="+d2s(FWHMa,0)+"nm";
Plot.addText(Text, 0, 0);
Plot.show();
run("Canvas Size...", "width=528 height=510 position=Bottom-Center zero");


// _______________________FWHM lateral_______________________ ;

selectWindow("Stack");
setSlice(OptSlice);

x = newArray(-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8);
y = newArray(17);
for (i=0; i<17; i++) {
	y[i] = getPixel(x2-8+i,y2);
	}

Fit.doFit("Gaussian", x, y); a=Fit.p(0); b=Fit.p(1); c=Fit.p(2); d=Fit.p(3);

// _______Plot

XplotLatReal = newArray(17); YplotLatReal = newArray(17);
Ymin=66000;
Ymax=0;
MaxGraph = 0;
for (i=0; i<17; i++) {
	XplotLatReal[i] = (i-8)*xyVoxel;
	YplotLatReal[i] = y[i];
	if (y[i]>=MaxGraph) {MaxGraph = y[i];}
	}
XplotLatFit = newArray(65); YplotLatFit = newArray(65);
for (i=0; i<65; i++) {
	XplotLatFit[i] = (i/4-8)*xyVoxel;
	X = i/4-8;
	YplotLatFit[i] =  a + (b-a)*exp(-(X-c)*(X-c)/(2*d*d));
	if (YplotLatFit[i]>=MaxGraph) {MaxGraph = YplotLatFit[i];}
	if (Ymin>YplotLatFit[i]) {Ymin=YplotLatFit[i];}
	if (Ymax<YplotLatFit[i]) {Ymax=YplotLatFit[i];}
	}

HM=(Ymax-Ymin)/2; k=-2*d*d*log((HM-a)/(b-a)); FWHMl=2*xyVoxel*sqrt(k);

Plot.create("FWHM lateral", "X", "Intensity", XplotLatFit, YplotLatFit);
Plot.setLimits(-8*xyVoxel, 8*xyVoxel, 0, MaxGraph*1.1);
Plot.add("circles", XplotLatReal, YplotLatReal);
//Text="FWHM lateral ="+d2s(FWHMl,0)+"nm";
Plot.addText("FWHM lateral ="+d2s(FWHMl,0)+"nm", 0, 0);
Plot.show();
run("Canvas Size...", "width=528 height=510 position=Top-Center zero");
run("Image Calculator...", "image1=[FWHM lateral] operation=Add image2=[FWHM axial]");
run("Canvas Size...", "width=550 height=550 position=Center");
rename("FWHM");
run("RGB Color");
selectWindow("FWHM axial");
close();
selectWindow("Stack");
close();

MainName=MainName+" FWHMa="+d2s(FWHMa,0)+"nm FWHMl="+d2s(FWHMl,0)+"nm";
run("Images to Stack");
rename(MainName);

setFont("Arial", 14, " antialiased");
run("Colors...", "foreground=red background=white selection=yellow");
setSlice(1);
drawString(date, 250, 260);
String1="FWHM lateral = "+d2s(FWHMl,0)+"nm"; String2="FWHM axial = "+d2s(FWHMa,0)+"nm";
drawString(String1, 250, 280); drawString(String2, 250, 300);
run("Colors...", "foreground=blue background=white selection=yellow");
setFont("Arial", 14, " antialiased");
drawString("Theoretical values", 250, 320);
if (NA=="1.4") {drawString("NA1.4 Oil FWHMl 188nm - FWHMa 696nm", 250, 340);}
if (NA=="1.32") {drawString("NA1.32 Oil FWHMl 199nm - FWHMa 781nm", 250, 340);}
if (NA=="1.3") {
	drawString("NA1.3 Oil FWHMl 202nm - FWHMa 803nm", 250, 340);
	drawString("NA1.3 Glyc FWHMl 202nm - FWHMa 771nm", 250, 360);
	}
if (NA=="0.95") {drawString("NA0.95 Air FWHMl 276nm - FWHMa 996nm", 250, 340);}
if (NA=="0.8") {
	drawString("NA0.8 Glyc FWHMl 328nm - FWHMa 2038nm", 250, 340);
	drawString("NA0.8 Air FWHMl 383nm - FWHMa 1403nm", 250, 360);
	}
if (NA=="0.75") {drawString("NA0.75 Air FWHMl 350nm - FWHMa 1598nm", 250, 340);}
if (NA=="1.45") {drawString("NA1.45 Oil FWHMl 181nm - FWHMa 649nm", 250, 340);}
if (NA!="0.75" && NA!="0.8" && NA!="1.3" && NA!="1.4" && NA!="1.45" && NA!="0.95" && NA!="1.32") {
	drawString("Values not determined yet", 250, 340);
	}

//__________Save FWHM values_____________________
InfoFWHM=ImageName+" "+date+" "+microscope+" "+MA+"x "+NA+" "+d2s(FWHMa,0)+" "+d2s(FWHMl,0);

if (Loc == "FMI") {
	path2 = S+S+"imagestore"+S+"FAIM"+S+"Maintenance_All microscopes"+S+"FWHMValues.txt";
} else {
	path = getDirectory("Choose a Directory");
	path2=path+"FWHMValues.txt";
}
if (File.exists(path2) == 1) {
	File.append(InfoFWHM, path2);
} else {
	headers = "ImageName date microscope MA NA FWHMa FWHMl";
	File.append(headers, path2);
	File.append(InfoFWHM, path2);
}

//________Save MIP __________________________

if (Loc == "FMI") {
path=S+S+"imagestore"+S+"FAIM"+S+"Maintenance_All microscopes"+S+MainName;
PathTest=path+".tif";
i=2;
while (File.exists(PathTest)==1) {
    path=path+"_2";
    PathTest=path+".tif";
    }

	saveAs("tiff", path);
} else {
path=path+MainName;
PathTest=path+".tif";
while (File.exists(PathTest)==1) {
    path=path+"_2";
    PathTest=path+".tif";
    }
	saveAs("tiff", path);
}








