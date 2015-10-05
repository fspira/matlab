import os
from os.path import join, basename, splitext, isfile, dirname
from glob import glob
import re
import datetime

from ij import IJ

from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from loci.formats import ImageReader
from loci.formats import MetadataTools
from loci.formats import ImageWriter
from loci.formats.tiff import TiffParser
from loci.formats.tiff import TiffSaver
from loci.common import RandomAccessInputStream
from loci.common import RandomAccessOutputStream
from ome.xml.model.enums import DimensionOrder
from ome.xml.model.enums import PixelType
from ome.xml.model.primitives import PositiveInteger


restr1 = ".*_DE_2_W(?P<well>\d+)_P(?P<site>\d+)_T(?P<iteration>\d+)"
restr2 = ".*_%s_W(?P<well2>\d+)_P(?P<site2>\d+)_T(?P<time>\d+).lsm"
filepattern = "%(expid)s_%(treatment)s_Mark%(well)s_Tl%(timelapse)s_Movie%(movie)03d.tif"


def setupDialog(n=1):
    """Opens a generic dialog, n is the number of experimental
    conditions."""

    gd = GenericDialog("Cat & Rename")
    gd.addStringField("Search for", "TR1_2")
    gd.addMessage("")
    gd.addStringField("Experiment ID", "Unnamed")
    gd.addStringField("Timelapse", "0")
    gd.addMessage("Experimental Conditions")

    for i in xrange(1, n+1, 1):
        gd.addStringField("W%04d" %i, "None", 20)
    gd.showDialog()

    pattern = gd.getNextString()

    description = {
        "expid": gd.getNextString(),
    	"timelapse": gd.getNextString(),
	}

    conditions = {}
    for i in xrange(1, n+1, 1):
        conditions['%04d' %i] = gd.getNextString()

    if gd.wasCanceled():
        return (None, None, None)
    else:
        return description, conditions, pattern


def getTimePoint(reader, omeMeta):
    """ Extract timeStamp from file """
    time = datetime.datetime.strptime(str(omeMeta.getImageAcquisitionDate(0)), "%Y-%m-%dT%H:%M:%S")
    time = [reader.getSeriesMetadataValue("TimeStamp0"),
            reader.getSeriesMetadataValue("TimeStamp #1"),
            omeMeta.getPlaneDeltaT(0,0)]
    time = [x for x in time if x is not None]
    return time[0]


def setUpXml(ome, image, files):
    """setup Xml standard for concatenated file"""
    ome.setImageID("Image:0", 0)
    ome.setPixelsID("Pixels:0", 0)
    ome.setPixelsDimensionOrder(DimensionOrder.XYCZT,0)

    if image.getBitDepth() == 8:
        pixels = PixelType.UINT8
    if image.getBitDepth() == 12:
        pixels = PixelType.UINT12
    if image.getBitDepth() == 16:
        pixels = PixelType.UINT16

    ome.setPixelsType(pixels, 0)
    ome.setPixelsSizeX(PositiveInteger(image.getWidth()), 0)
    ome.setPixelsSizeY(PositiveInteger(image.getHeight()), 0)
    ome.setPixelsSizeZ(PositiveInteger(image.getNSlices()), 0)
    ome.setPixelsSizeT(PositiveInteger(len(files)), 0)
    ome.setPixelsSizeC(PositiveInteger(image.getNChannels()), 0)
    return ome

def concatenateImagePlus(files, outfile):
    """Concatenate images contained in files and save in outfile"""

    options = ImporterOptions()
    options.setId(files[0])
    options.setVirtual(1)
    options.setOpenAllSeries(1)
    options.setQuiet(1)
    images = BF.openImagePlus(options)
    imageG = images[0]
    nrPositions = len(images)
    options.setOpenAllSeries(0)

    nslices = imageG.getNSlices()
    nframes = len(files)
    nchannels = imageG.getNChannels()
    luts = imageG.getLuts()

    for i in range(0, nrPositions):
        concatImgPlus = IJ.createHyperStack(
            "ConcatFile", imageG.getWidth(), imageG.getHeight(),
            imageG.getNChannels(), imageG.getNSlices(),
            len(files), imageG.getBitDepth())


        concatStack = ImageStack(imageG.getWidth(), imageG.getHeight())
        IJ.showStatus("Concatenating files")
        for file_ in files:
            try:
                print '...', basename(file_)
                options.setSeriesOn(i, 1)
                options.setId(file_)
                image = BF.openImagePlus(options)[0]
                imageStack = image.getImageStack()

                sliceNr = imageStack.getSize()
                for j in range(1, sliceNr+1):
                    concatStack.addSlice(imageStack.getProcessor(j))
                image.close()
                options.setSeriesOn(i, 0)
            except Exception, e:
                IJ.log("ERROR")
                IJ.log(file_ + str(e))
                raise
            IJ.showProgress(files.index(file_), len(files))

        concatImgPlus.setStack(concatStack, nchannels, nslices, nframes)
        concatImgPlus.setCalibration(image.getCalibration())
        concatImgPlus.setOpenAsHyperStack(True)
        concatImgPlus.setLuts(luts)
        concatImgPlus.close()
        IJ.saveAs(concatImgPlus, "Tiff",  outfile)

    return outfile


def processMovie(root, files, outfile):
    """Concatenate images and write ome.tiff file.
    If image contains already multiple time points just copy the image"""

    files.sort()

    options = ImporterOptions()
    options.setId(files[0])
    options.setVirtual(1)

    image = BF.openImagePlus(options)
    image = image[0]
    if image.getNFrames() > 1:
        msg = ("%s Contains multiple time points. Can only concatenate"
               " single time points!" %files[0])
        raise RuntimeError(msg)
        image.close()

    reader = ImageReader()
    reader.setMetadataStore(MetadataTools.createOMEXMLMetadata())
    reader.setId(files[0])
    timeInfo = []
    omeOut = reader.getMetadataStore()
    omeOut = setUpXml(omeOut, image, files)
    reader.close()
    image.close()
    itime = 0

    for fileName in files:
        omeMeta = MetadataTools.createOMEXMLMetadata()
        reader.setMetadataStore(omeMeta)
        reader.setId(fileName)
        timeInfo.append(getTimePoint(reader, omeMeta))

        nrImages = reader.getImageCount()
        for i in range(0, reader.getImageCount()):
            try:
                dT = round(timeInfo[files.index(fileName)]-timeInfo[0],2)
            except:
                dT = (timeInfo[files.index(fileName)]-timeInfo[0]).seconds
            omeOut.setPlaneDeltaT(dT, 0, i + itime*nrImages)
            omeOut.setPlanePositionX(omeOut.getPlanePositionX(0,i), 0, i + itime*nrImages)
            omeOut.setPlanePositionY(omeOut.getPlanePositionY(0,i), 0, i + itime*nrImages)
            omeOut.setPlanePositionZ(omeOut.getPlanePositionZ(0,i), 0, i + itime*nrImages)
            omeOut.setPlaneTheC(omeOut.getPlaneTheC(0,i), 0, i + itime*nrImages)
            omeOut.setPlaneTheT(omeOut.getPlaneTheT(0,i), 0, i + itime*nrImages)
            omeOut.setPlaneTheZ(omeOut.getPlaneTheZ(0,i), 0, i + itime*nrImages)
        itime = itime + 1
        reader.close()

        IJ.showProgress(files.index(fileName), len(files))

    try:
        omeOut.setPixelsTimeIncrement(float(dT/(len(files)-1)), 0)
    except:
        omeOut.setPixelsTimeIncrement(0, 0)

    if len(files) <= 1:
        raise RuntimeError('Found only one file. Nothing to concatenate')

    outfile = concatenateImagePlus(files, outfile)
    filein = RandomAccessInputStream(outfile)
    fileout = RandomAccessOutputStream(outfile)
    saver = TiffSaver(fileout, outfile)
    saver.overwriteComment(filein, omeOut.dumpXML())
    fileout.close()
    filein.close()


def findAndCat(source_dir, target_dir):
    movies = dict()

    for root, dirs, files in os.walk(indir, topdown = True):
        if root == indir:

            n = len([d for d in dirs if re.match(".*_W(?P<mrk>\d+)", d)])
            description, conditions, pattern = setupDialog(n)
            if description is None:
                raise SystemExit()

            regex2 = re.compile(restr2 %pattern)
            regex1 = re.compile(restr1)

        for i, file_ in enumerate(files):
            ifile = join(root, file_)
            ifile2 = ifile.split(os.sep)
            wellpath = os.sep.join(ifile2[:-2])
            filepath = os.sep.join(ifile2[-2:])

            # grid regex
            match1 = regex1.match(wellpath)
            match2 = regex2.match(filepath)
            if match1 is not None and match2 is not None:
                wdict = match1.groupdict()
                fdict = match2.groupdict()

                if int(fdict['time']) == 1:
                    wdict.update({'movie': len(movies),
                                  'treatment': conditions[wdict["well"]]})
                    wdict.update(description)

                    key = dirname(ifile)
                    if movies.has_key(key):
                        raise RuntimeError('Key not unique')
                    movies[key] = wdict

                    files_ = glob(join(root, "*.lsm"))
                    outfile = join(target_dir, filepattern %wdict)
                    print "New file:" , basename(outfile)
                    processMovie(root, files_, outfile)

    return len(movies)

if __name__ == "__main__":
    indir = IJ.getDirectory("Choose the input directory")
    if indir is None:
        raise SystemExit
    outdir = IJ.getDirectory("Choose the output directory")
    if outdir is None:
        raise SystemExit

    # indir = "C:\Users\hoefler\Desktop\mergescript\data"
    # outdir = "C:\Users\hoefler\Desktop\mergescript\output"

    nmovies = findAndCat(indir, outdir)
    print "Concatenated %d movies" %nmovies
