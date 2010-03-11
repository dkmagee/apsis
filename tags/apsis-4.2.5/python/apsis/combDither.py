#!/usr/bin/env python

# $Id: combDither.py,v 1.47 2006/08/29 06:50:33 jpb Exp $
# ---------------------------------------------------------------------

__version__      = '$Revision: 1.47 $ '[11:-3]
__version_date__ = '$Date: 2006/08/29 06:50:33 $ '[7:-3]
__author__       = 'J Blakeslee <jpb@pha.jhu.edu>'

import os,string,glob,math
import numpy
import pyfits,matutil
import xydrizzle
import xmlUtil,fUtil,pyblot,augmask
import astrometer
from   pUtil import ptime
from   msg   import pMessage
from   sys   import version
pyversion = version
from   pyraf import iraf

# get the version of drizzle by pretending to run it
tmplist = iraf.drizzle("None","None",Stdout=1)
drversion = "unknown"
for line in tmplist:
    wordlist = line.split()
    for __ii in range(len(wordlist)):
        word = wordlist[__ii]
        if word.lower() == 'version':
            drversion = wordlist[__ii+1]
            break
        del word
    del wordlist
del line,tmplist

try:
    pydriz_version = xydrizzle.__version__
except:
    pydriz_version = xydrizzle.version


# Good bits marked by a flag value in the DQ array are
#
#    0  good pixel
#    2  data replaced by fill value
#   64  pre-existing hot pixel
#  128  bias level pixel
#  256  saturation (full well or a-to-d)
# 2048  a-to-d saturation*
# 4096  reserved (possibly for Multidrizzle CR rejected pixels)
# 8192  cosmic ray rejected during cr-split image combination

# Recent changes to DQ bit flags values has indicated that a DQ pixel value of
# 128 should not now be considered "good" as that value, a bias level pixel, now
# indicates bad columns as well.
# 128 bias level pixel (bad columns); was in BPIXTAB, will be in superbias DQ

#_goodBits_ = 2+64+128+256+2048+4096+8192

_goodBits_ = 2+64+256+2048+4096+8192   # good-pixel summed value

#_goodBits_ = 8578                     # CALACS good-pixel summed value (see run_drizzle docstring)

_rnVarSuppFac_ = 1.38                  # boost N*rn^2 by this factor for RMS images
_exVarPerSec_  = 0.02222               # add this much extra var per sec for RMS ims

class drizzleImage:
    """Some notes on the combDither module and its implementation of xydrizzle.

    combDither currently contains one class, drizzleImage, which defines
    methods for drizzling dithered image associations to a final images,
    including cosmic ray rejection, etc.  The original intention was to
    create other classes with methods for combining the images using other
    routines (e.g., SWARP, or something based on python numerical arrays).

    The basic procedure for processing a dataset through drizzle:

    drob = drizzleImage(obs,al.MatchDict) - create instance of this class
    drob.run_all()          -   run the drizzle/blot/drizzle sequence
    drob.makeFlagImage()    -   create Sextractor flag image for each asn
    drob.makeRmsImage()     -   create Sextractor rms image for each asn
    drob.fixAstrometry(obs) -   fix the astrometry of the images if possible.

    [The fixAstrometry method must receive a DataSet object (obs) and will generally
    be called by the apsis pipeline.  This behaviour is not internally consistent
    but the method is a bit of hack to removed the need for the pipeline to
    deal with a new object (an instance of astrometer.gscMatchup) and all the
    attendant messaging for that.  N.B. This method should not be called internally.]

    Alternatively, the user can replace run_all() with run_drizzle() if no
    cosmic ray rejection is desired.  This module no longer makes its own
    association tables, but just uses the original ones that are given as
    inputs to the pipeline.

    The run_all method does three things:

    1. run_drizzle() drizzles images onto separate output images
    2. runs pyblot (median stacks, blots, derivs, cr_rej) which
        returns the list of cr masks for 2nd drizzling
    3. run_drizzle() to drizzles images all together using CR masks
        from pyblot

    run_drizzle() is the heart of the drizzleImage class.  It runs pydrizzle
    on all of the association tables in the asnTableList and takes a number
    of optional parameters (outshape=None, deltmp=1, units='counts',
    bits=_goodBits_, separDriz=0, crmasks=None, exptimePyDrKludge=0).
    Here's what the params are for

    outshape:  x,y sizes of output drizzled image (def: determined 'otf')
    deltmp:    delete temporary (masks, coeffs) files? (def: 1 [yes])
    units:     units of output pixel values (def: counts)
    bits:      sum of all possible good-pixel values in dq array
          ('None' means don't use masks).
          The default value bits=8450 is the sum of the pix values
          that we accept as being good, bits = 0+2+256+8192 = 8450,
          which includes 'good', 'replaced', 'saturated',and
          'CR-repaired', but NOTE THAT THESE VALUES ARE PARTICULAR 
          TO CALACS! [ISR ACS-99-03]

    The remainder of the params relate to pyblot setup and results:
    crmasks:            List of (mask,nzap) tuples from pyblot
    separDriz:          Drizzle onto separate output images?

    run_drizzle loops over all associations and runs pydrizzle on them.  
    It measures the shift 'zeropoint' based on the mean of all the measured
    shifts and then sets the relative shifts accordingly, ensuring that all
    associations will have the same absolute shift zeropoint, as well as the
    same image size.  This *seems* to be ok, although the output WCS may then
    disagree with pydrizzle's idea of the WCS of the output product, but drizzle
    itself constructs the WCS image in the output product from the input and
    whatever image headers it gets.

    For each association, run_drizzle creates the pydrizzle object, then
    goes through the parlist (list of dictionaries holding input drizzle
    parameters).  It modifies these parlists in various way, e.g., it makes
    the output images different if 'separDriz' has been set.  It also sets
    the x,y shifts and rotations to whatever is in the MatchDict (thus,
    making delta shifts in the asn tables no longer necessary).  The list of
    modified parlists is saved as a class attribute.  It then runs drizzle.

    After the 2nd drizzle (separdriz = 0), the code checks if:
    (1) the number of science extensions in the multi-extension FITS (mef)
    file is greater than 1 and units='counts' (the default); if so, it
    divides the count levels by the number of extensions.
    (2) if the pydrizzle version is less than 1.4; if so and Nsci>1, then it
    divides the output exposure times by the number of science extensions.

    If it's the final drizzle, then it also goes through all the input images
    and sums the alignSky values in the headers and writes the sum to the
    science image header.  That's it.  The makeFlagImage and makeRmsImage
    methods are pretty straightforward.

    The drizzle kernel, for the *final drizzle only* can now be set when
    run_drizzle or run_all is called. The drizzle-default linear ('square')
    kernel is the fastest, and tends to smooth over bad pixels better,
    but the psf is not as tight and the noise correlation is much worse.

    The damped sinc interpolant ('lanczos3') is sometimes called the
    optimal kernel, though beware that you may get holes (negative
    pixels) in the output.  These holes can happen anywhere, though more
    commonly near stars, etc.  Because of the masking of bad pixels and
    columns (based on the DQ array), occasionally a pixel in the output
    image will have very little data from the input image, and by a
    quirk of the interpolation kernel with its negative sidelobes, will
    end up with a negative value.  Richard Hook said to check the
    drizzle weight image to make sure that those negative pixels are
    given very little weight; if they have low weight, just don't use
    them, but if their weight is high, then there might be a problem.
    I checked several of these pixels, and found they were given low
    weight, about 1/10 of the median, so it seems ok.  Of course, if
    there are multiple dithered images, these negative pixels get filled
    in quickly by other pixels with much higher weight.  For the linear
    kernel you don't get this effect because there are no negative
    sidelobes.

    One caveat: run_drizzle overwrites shifts, rotations, and output product
    names in the pydrizzle parlist rather than setting them through methods.
    This was historically necessary for the shifts (and *still* the most
    straightforward way), and it's still necessary for drizzling onto
    separate outputs for blotting.  So, I wouldn't recommend changing it,
    but we have to be careful in case later pydrizzles change the parlist
    elements or their meanings.

    """

    def __init__(self, obs, MatchDict, alignSky='ALIGNSKY',hdrGain=0, suppInstVar=1, crlower=None,
                 smartstack=0, notrim=0, padfac=None, outsize=None, outshift=None,origscale=None,noContext=None,
                 maskFile=None,dfilts=(None,None)):
        
        self.modName     = string.split(string.split(str(self))[0],'.')[0][1:]
        self.root        = obs.newobspath     # root path of the observation dir
        self.obsName     = obs.newobs
        self.obsAsnDict  = obs.asnDict
        self.obsAlign    = obs.newalign
        self.obsPars     = obs.newpar
        self.obsFits     = obs.newfits        # path to fits files
        self.obsFitsList = obs.fitslist
        self.messagedir  = obs.messagedir     # where the module message will go.
        self.logfile     = obs.logfile
        self.allShifts   = obs.allShifts      # list of tuples; 2nd element is dictionary of shifts
        self.alignSky    = alignSky           # header keyword for subtracted Sky's
        self.detector    = obs.detector       # detector matters for drizzle shifts
        self.refcdmat    = obs.refcdmat       # dictionary with final cd matrix
        self.refotherkeys = obs.refotherkeys  # dictionary with other keys relating to the wcs solution
        self.edgeMaskName= obs.edgeMaskName = 'Edgemask.fits'
        self.reflogfile_add = obs.reflogfile_add
        self.errorList   = []
        self.inputList   = []
        self.outputList  = {}
        self.hdrGain     = hdrGain
        self.crlower     = crlower
        self.smartstack  = smartstack
        self.notrim      = notrim
        self.padfac      = padfac
        self.outsize     = outsize
        self.outshift    = outshift
        self.maskFile    = maskFile
        self.origscale   = origscale
        self.noContext   = noContext
        self.dfilts      = dfilts
        self.suppInstVar = suppInstVar
        self.MatchDict   = MatchDict          # MatchDict, from the alignImage object is
                                              # a dictionary of dict's for each image
                                              
        # lists of the various image types . . .
        self.sciImageList     = obs.sciImageList     = []   
        self.contextImageList = obs.contextImageList = []
        self.weightImageList  = obs.weightImageList  = []
        self.flagImageList    = obs.flagImageList    = []
        self.rmsImageList     = obs.rmsImageList     = []
        self.removeList = []
        self.logfile.write('Instantiating drizzleImage object for observation: '+self.obsName)


        if not self.origscale:
            if obs.detector =='WFC':
                self.origscale = 0.05
            elif obs.detector == 'HRC':
                self.origscale = 0.025
            elif obs.detector == 'SBC':
                self.origscale = 0.025
            else:
                errtxt = "Detector type "+ obs.detector + " origscale unknown.\n"
                self.errorList.append((self.modName,errtxt))
                self.logfile.write(errtxt)
                raise Exception,errtxt
                
        # set up the inputList for mkMsg() method
        for i in self.obsFitsList:
            self.inputList.append(i)
            self.NumSci = 1
            _firstN = 0
            for im in self.MatchDict.keys():
                if not _firstN:
                    _firstN = self.MatchDict[im]['NumSci']
                else:
                    if self.MatchDict[im]['NumSci'] != _firstN:
                        self.NumSci = 0
                        errtxt = 'WARNING MatchDict purports that images to be combined'+\
                                 ' have different numbers of sci extensions. \n'+\
                                 ' CANNOT fix drizzle scale factor problem in this case!!!!'
                        self.logfile.write(errtxt)
                        self.errorList.append((self.modName,errtxt))
                        break
            if self.NumSci: self.NumSci = _firstN
            self.logfile.write('setting iraf min_lenuserarea=640000')
            iraf.set(min_lenuserarea=640000)


    def run_all(self, clean_up=0, userKeepBits=None, units='counts',
                asecpix=None, pixfrac=None, kernel=None,
                noRej=None):
        
        """ Run_all()
        This is it, the whole enchilada: drizzle, blot, drizzle again.
        """
        if noRej:
            # then just drizzle once
            self.logfile.write('Skipping cosmic ray rejection in combDither!')
            if not kernel:
                ## kernel = 'square' # default to lanczos3, Feb 25 02, anderson
                kernel = "lanczos3"
            
            if userKeepBits:
                print "Using user supplied good bits value of",userKeepBits
                self.run_drizzle(deltmp=clean_up, bits=userKeepBits, separDriz=0,
                                 crmasks=None, units=units, asecpix=asecpix, pixfrac=pixfrac,
                                 kernel=kernel, delwght=0)
            else:
                print "Using default good bits value of",_goodBits_
                self.run_drizzle(deltmp=clean_up, bits=_goodBits_, separDriz=0,
                                 crmasks=None, units=units, asecpix=asecpix, pixfrac=pixfrac,
                                 kernel=kernel, delwght=0)
            
            self.logfile.write('Drizzling complete.')
            return

        # otherwise, do the full cycle        
        self.logfile.write('Starting initial drizzle process...')
        
        if userKeepBits:
            print "using user supplied good bits value of",userKeepBits
            self.run_drizzle(deltmp=0, bits=userKeepBits, separDriz=1, crmasks=None,
                             units='counts', delwght=clean_up)
        else:
            print "using default good bits value of",_goodBits_
            self.run_drizzle(deltmp=0, bits=_goodBits_, separDriz=1, crmasks=None,
                             units='counts', delwght=clean_up)
        
        print "Creating pyblot object."
        self.logfile.write('getting blotter object.')        
        if self.crlower:
            self.logfile.write('crlower flag set; will use tighter driz_cr SNR thresholds.')
        PyBlOb = pyblot.blotter(self.shortparlists,self.parlists,self.obsFits,self.logfile, skyKey=self.alignSky,hdrGain=self.hdrGain,\
                                crlower=self.crlower, clean_up=clean_up, imNsci=self.NumSci)

        self.logfile.write('Blotting...')
        try:
            crmasks = PyBlOb.runblots()
            # elements of crmasks dict are now (maskname,nzero) tuples!
        except Exception,err:
            errtxt = "ERROR: runblots method failed.\n\t"+str(err)
            self.errorList.append((self.modName,errtxt))
            self.logfile.write(errtxt)
            raise Exception,err
        self.logfile.write('blotting complete.')

        # get the crmasks listed in the outputList
        for key in crmasks.keys():
            self.outputList[crmasks[key][0]] = [key]   # can't get at any predecessor info here.

        if not self.notrim:
            self.logfile.write('Assembling Edge mask...')
            auger = augmask.augmask(self.augimdict,self.edgeMaskName,logfile=self.logfile,clean_up=clean_up)
            auger.make()
            # test to see if edge mask has been made.
            if os.path.isfile(os.path.join(self.obsFits,self.edgeMaskName)):
                self.outputList[self.edgeMaskName] = []
        else:
            self.logfile.write('Skipping Edge mask construction.')

        self.logfile.write('Starting second drizzle process...')
        if asecpix:
            self.logfile.write('  using requested output scale: '+str(asecpix)+' asec/pix.')
        if pixfrac:
            self.logfile.write('  using requested drizzle pixfrac: '+str(pixfrac))
        if not kernel:
            ## kernel = 'square' # default to lanczos3, Feb 25 02, anderson
            kernel = "lanczos3"

        if userKeepBits:
            self.logfile.write('  using bits parameter = '+str(userKeepBits))
            self.run_drizzle(deltmp=clean_up, bits=userKeepBits,separDriz=0,crmasks=crmasks,units=units,\
                         asecpix=asecpix, pixfrac=pixfrac, kernel=kernel, delwght=0)
        else:
            self.logfile.write('  using bits parameter = '+str(_goodBits_))
            self.run_drizzle(deltmp=clean_up, bits=_goodBits_,separDriz=0,crmasks=crmasks,units=units,\
                         asecpix=asecpix, pixfrac=pixfrac, kernel=kernel, delwght=0)

        if clean_up:
            self.logfile.write('Removing 1st pass output files.')
            self._clean_driz()
        else:
            self.logfile.write('Keeping 1st pass drizzle output files.')

        self.logfile.write('Full drizzle cycle complete.')
        del PyBlOb
        return

    def _clean_driz(self):
        curdir = os.getcwd()
        os.chdir(self.obsFits)
        for file in self.removeList:
            try:
                os.remove(file)
                self.logfile.write("removed "+file)
            except Exception,err:
                self.logfile.write(str(err))
                print err
        os.chdir(curdir)
        return

    def _get_refAsn(self):
        "choose the reference association for the output image size"
        refAsn = None
        maxshift = -9999.
        
        for tab in self.obsAsnDict.keys():
            xshifts = []
            yshifts = []
            for im in self.obsAsnDict[tab]:
                xshifts.append(self.MatchDict[im]['xpix_shift'])
                yshifts.append(self.MatchDict[im]['ypix_shift'])
            xshifts.sort()
            yshifts.sort()
            _maxi = max((xshifts[-1] - xshifts[0]),(yshifts[-1] - yshifts[0]))
            if _maxi  > maxshift:
                maxshift = _maxi
                refAsn = tab
            del _maxi,xshifts,yshifts

        del tab
        if not refAsn:
            self.logfile.write('Error: ref Asn could not be determined.')
            raise RuntimeError,'Bug in _get_refAsn method.'
        self.logfile.write(refAsn+(' chosen as reference asn table.  maxshift: %.4f'%maxshift))
        print refAsn+' chosen as reference asn table.  maxshift: %.4f'%maxshift
        return refAsn

    def run_drizzle(self,outshape=None, deltmp=1, units='counts',
                    bits=_goodBits_, separDriz=0, crmasks=None,
                    asecpix=None, pixfrac=None, kernel='square', delwght=0):
        
        """ run_drizzle():
        Run drizzle on all asn tables in the obsAsnDict.keys()
           Input parameters (all optional) include:
             outshape: x,y sizes of output drizzled image (def: determined 'otf')
             deltmp:   delete temporary (masks, coeffs) files (def: 1 [yes])
             units:    units of output pixel values (def: counts)
             bits:     sum of possible good-pixel values in dq array (None => don't use masks).
                       The default bits=8578 is the sum of the pix values that we accept
                       as being good, bits = 0+2+256+8192+128 = 8578, which includes
                       'good', 'replaced', 'saturated', 'CR-repaired' and 'bias-level pixel' but
                       NOTE THAT THESE VALUES ARE PARTICULAR TO CALACS! [ISR ACS-99-03]
             The remainder of the params relate to pyblot setup and results:
             crmasks:      List of masks from pyblot
             separDriz:    Drizzle onto separate output images
             asecpix:      Scale in arsec/pix of final drizzled images.
             pixfrac:      Drizzle 'pixfrac' or 'dropsize' parameter.
             kernel:       interpolation kernel used by drizzle.
        """
        
        curdir = os.getcwd()
        self.logfile.write('Moving to directory '+self.obsFits+' for drizzling.')
        os.chdir(self.obsFits)       # move into the FITS dir.
        self.logfile.write('Using '+kernel+' drizzle kernel.')
        # will pick the output shape on-the-fly later if none specified
        self.outshape = outshape
        del outshape
        if bits==None:
            self.logfile.write('will drizzle without masks')
            if crmasks:
                warntxt='Warning: sent crmask dictionary but specified bits=None'
                self.errorList.append((self.modName,warntxt))
                self.logfile.write(warntxt)
        else:
            self.logfile.write('Will drizzle with masks; summed goodpix vals = '+str(bits))
            print '  Will drizzle with masks; summed goodpix vals =',bits
        # use drizzle pixfrac of 1 if not specified
        if not pixfrac:
            self.pixfrac = 1
        else:
            self.pixfrac = pixfrac
            self.logfile.write('setting pixfrac = '+str(pixfrac))
        self.asecpix = asecpix
        self.kernel  = kernel
        del pixfrac,asecpix,kernel

        # This bunch of thrashing is done to prevent multiple calls of this
        # method appending redundant or incorrect file names to the image
        # lists.  It is done this way so that the obs object passed to the
        # constructor also gets updated.  Just the setting lists to zero,
        # i.e.  self.sciImageList = [], and then appending items does not
        # result in the obs.sciImageList being updated as well.  
        while self.sciImageList:
            del self.sciImageList[0]
        while self.contextImageList:
            del self.contextImageList[0]
        while self.weightImageList:
            del self.weightImageList[0]

        # important to nuke any lingering inmasks!
        pyblot.junkem("*inmask*")
        self.parlists = []
        self.shortparlists = []
        self.wcslist  = []
        iraf.flpr()
        iraf.flpr()

        ### all this aug stuff is for making the mask for trimming the edges
        if separDriz:
            self.augimdict = {}
            for tab in self.obsAsnDict.keys():
                self.augimdict[tab] = {}
                self.augimdict[tab]['maskname'] = tab.split('.')[0][:-3]+'augmask.fits'
                self.augimdict[tab]['ctxlist'] = []
        ####################################################################

        for tab in self.obsAsnDict.keys():
            if self.outshape == None:
                # choose the reference Asn!
                self.refAsn = self._get_refAsn()
                REFtab = self.refAsn
                #
                if self.asecpix:
                    _skyf = xydrizzle.SkyField()
                    _skyf.set(psize=self.asecpix) 
                    PyDrOb = xydrizzle.PyDrizzle(REFtab, bits=None, field=_skyf, pixfrac=self.pixfrac, kernel=self.kernel,
                                                 filter1=self.dfilts[0], filter2=self.dfilts[1])
                    del _skyf
                else:
                    PyDrOb = xydrizzle.PyDrizzle(REFtab, bits=None, pixfrac=self.pixfrac, kernel=self.kernel,
                                                 filter1=self.dfilts[0], filter2=self.dfilts[1])

                # OK, we've just created a PyDrizzle object to serve
                # as a prototype for later ones, but not be used itself...
                # pad the axis sizes a little 
                _maxAng = 0
                for imkey in self.MatchDict.keys():
                    delta45ang = 45 - abs(45 - abs(divmod(self.MatchDict[imkey]['angle'],90)[1]))
                    if delta45ang > _maxAng:
                        _maxAng = delta45ang
                if(_maxAng > 5):
                    _padx,_pady = self._rotateRect(_maxAng,PyDrOb.parlist[0]['outnx'],PyDrOb.parlist[0]['outny'])
                    if self.padfac:
                        _padx *= self.padfac
                        _pady *= self.padfac
                    self.logfile.write("padding x,y by factors "+str(_padx)+" "+str(_pady)+" including rotation.")
                else:
                    if self.padfac:
                        _padx = _pady = self.padfac
                    else:
                        _padx = _pady = 1.025

                # if self.outsize:
                #     _outNx,_outNy = self.outsize
                # else:
                #     _outNx = int(_padx * PyDrOb.parlist[0]['outnx'])
                #     _outNy = int(_pady * PyDrOb.parlist[0]['outny'])
                # self.logfile.write("Output image size will be %d x %d"%(_outNx,_outNy))
                # self.outshape  = (_outNx,_outNy)
                # del _maxAng,delta45ang,_padx,_pady
                #######################################
                self.protoWCS = PyDrOb.observation.product.geometry.wcs
                
                # check to see if scale of protoWCS is same as specified asecpix
                if self.asecpix:
                    if abs(self.protoWCS.pscale - self.asecpix) > 1e-9:
                        errtxt = 'Error: self.protoWCS.pscale: '+str(self.protoWCS.pscale)+\
                                 '  !=  self.asecpix: '+str(self.asecpix)
                        self.logfile.write(errtxt)
                        raise Exception,errtxt
                else:
                    self.asecpix = self.protoWCS.pscale
                
                if self.outsize:
                    _outNx = int(self.outsize[0] / self.asecpix + 0.5 + 2*separDriz)
                    _outNy = int(self.outsize[1] / self.asecpix + 0.5 + 2*separDriz)
                else:
                    _outNx = int(_padx * PyDrOb.parlist[0]['outnx'])
                    _outNy = int(_pady * PyDrOb.parlist[0]['outny'])
                self.logfile.write("Output image size will be %d x %d"%(_outNx,_outNy))
                self.outshape  = (_outNx,_outNy)
                del _maxAng,delta45ang,_padx,_pady
                # self.ra_coord  = PyDrOb.observation.product.geometry.wcs.crval1
                # self.dec_coord = PyDrOb.observation.product.geometry.wcs.crval2
                # self.crpix1  = PyDrOb.observation.product.geometry.wcs.crpix1
                # self.crpix2  = PyDrOb.observation.product.geometry.wcs.crpix2
                ## Set what we can using the SkyField object.
                _skyf = xydrizzle.SkyField()
                _skyf.set(shape=self.outshape, psize=self.asecpix, \
                          ra=self.protoWCS.crval1, dec=self.protoWCS.crval2) 

                # now, remake the PyDrOb object, just to be sure of the right
                # zeropoint shifts for an output image of this size (this SkyField)
                del PyDrOb, self.protoWCS
                PyDrOb = xydrizzle.PyDrizzle(REFtab, bits=None, field=_skyf, pixfrac=self.pixfrac, kernel=self.kernel,
                                             filter1=self.dfilts[0], filter2=self.dfilts[1])
                self.protoWCS = PyDrOb.observation.product.geometry.wcs
                # paranoid check
                if abs(self.protoWCS.pscale - self.asecpix) > 1e-9:
                    errtxt = 'Error: self.protoWCS.pscale: '+str(self.protoWCS.pscale)+\
                             '  !=  self.asecpix: '+str(self.asecpix)
                    self.logfile.write(errtxt)
                    raise Exception,errtxt
            
                # calculate zeropoints of the x,y shifts, using *all* images in MatchDict
                _dX_EMP = []
                _dY_EMP = []
                ###for imdict in PyDrOb.parlist:
                for imkey in self.MatchDict.keys():
                    x_MsdShift = self.MatchDict[imkey]['xarc_shift'] / self.asecpix
                    y_MsdShift = self.MatchDict[imkey]['yarc_shift'] / self.asecpix
                    _dX_EMP.append(x_MsdShift)
                    _dY_EMP.append(y_MsdShift)
                del x_MsdShift, y_MsdShift
                self.logfile.write('dx_EMP '+ str(_dX_EMP))
                self.logfile.write('dy_EMP '+ str(_dY_EMP))

                ## this centers the refimage in the output image:
                # self.xZPT = 0
                # self.yZPT = 0
                if self.outshift:
                    self.xZPT = 1.0*int(0.5 + (self.outshift[0] / self.asecpix))
                    self.yZPT = 1.0*int(0.5 + (self.outshift[1] / self.asecpix))
                    self.logfile.write('using specified arcsec shift zeropoints '+\
                                       str(self.outshift[0])+' '+str(self.outshift[1]))
                else:
                    # this centers the avg image position:
                    self.xZPT = -1.0*int((max(_dX_EMP) + min(_dX_EMP))/2.0)   # preserves fractional
                    self.yZPT = -1.0*int((max(_dY_EMP) + min(_dY_EMP))/2.0)   # part of image shifts
                self.logfile.write('selected pixshift zeropoints '+str(self.xZPT)+' '+str(self.yZPT))
                del PyDrOb, _dX_EMP, _dY_EMP, REFtab

            _skyf.exptime = None
            PyDrOb = xydrizzle.PyDrizzle(tab, field=_skyf, clean=deltmp, bits=bits, units=units,\
                                         pixfrac=self.pixfrac, kernel=self.kernel,
                                         filter1=self.dfilts[0], filter2=self.dfilts[1])

            # verify all subsequent images get same pscale
            if abs(PyDrOb.observation.product.geometry.wcs.pscale - self.asecpix) > 1e-9:
                errtxt='Error: '+tab+' PyDrOb wcs.pscale: '+str(PyDrOb.observation.product.geometry.wcs.pscale)\
                        + '  !=  self.asecpix: '+str(self.asecpix)
                self.logfile.write(errtxt)
                raise Exception,errtxt

            # Now mask out specific regions if specified from the command line.
            # Usually, this is to explicitly remove figure 8's or other such
            # artifacts from the F850LP exposures...  

            if self.maskFile:
                maskreg_cmd = 'maskreg %s %s' % (self.maskFile,self.obsFits)
                self.logfile.write(maskreg_cmd)
                os.system(maskreg_cmd)

            # OK, PyDrOb should now have the desired shape and CRVAL's, and I'd have
            # thought we have to make sure the Drizzle Product has a consistent CRPIX
            # for this modified shifts, but drizzle seems to take care of it.
            # Now have to update the individual shifts accordingly...
            # Also, we need to change output product names if separDriz=1
            # and include crmasks if they've been sent.
            if separDriz and (self.MatchDict[self.MatchDict.keys()[0]]['NumSci'] == 2):
                if self.smartstack:
                    output = open('medianimages_input','w')
                    S = '%i\n' % (len(PyDrOb.parlist))
                    output.write(S)        
                    for imdict in PyDrOb.parlist:
                        imkey = string.split(imdict['data'],'[')[0]
                        x_MsdShift = self.MatchDict[imkey]['xarc_shift'] / self.asecpix
                        y_MsdShift = self.MatchDict[imkey]['yarc_shift'] / self.asecpix
                        xsh = (x_MsdShift + self.xZPT + _outNx/2)*self.asecpix
                        ysh = (y_MsdShift + self.yZPT + _outNy/2)*self.asecpix
                        S = '%s %g %g %g\n' % (imdict['data'],xsh,ysh,-(self.MatchDict[imkey]['angle']))
                        output.write(S)
                    output.close()
                    medianfilter_cmd = "whmedian medianimages_input medianimages %g %g 10" % (int(float(_outNx)*self.asecpix+0.9),int(float(_outNy)*self.asecpix+0.9))
                    self.logfile.write(medianfilter_cmd)
                    os.system(medianfilter_cmd)
                else:
                    miout = open('medianimages', 'w')
                    for imdict in PyDrOb.parlist:
                        miout.write('%s\n' % imdict['data'])
                    miout.close()
                        

            for imdict in PyDrOb.parlist:
                imkey = string.split(imdict['data'],'[')[0]
                x_MsdShift = self.MatchDict[imkey]['xarc_shift'] / self.asecpix
                y_MsdShift = self.MatchDict[imkey]['yarc_shift'] / self.asecpix

                imdict['xsh'] = x_MsdShift + self.xZPT
                imdict['ysh'] = y_MsdShift + self.yZPT
                ### changing to negative sign for drizzle convention (12/Apr/2002)
                imdict['rot'] = -(self.MatchDict[imkey]['angle'])

                tmpstring = string.split(string.split(imdict['data'],'[')[1],']')[0]
                Xname = (string.split(tmpstring,',')[0]).upper()
                Xver  = int(string.split(tmpstring,',')[1])
                basefits   = string.split(imkey,'.')[0]+"_"+Xname+"_"+str(Xver)+'.fits'
                if separDriz:
                    imdict['outdata']    = '_dz_'+basefits
                    imdict['outcontext'] = '_cx_'+basefits
                    imdict['outweight']  = '_wt_'+basefits
                    self.augimdict[tab]['ctxlist'].append(imdict['outcontext'])

                    # have to make sure exptime is right for a single image
                    _oldtexp = str(imdict['texptime'])
                    imdict['texptime']   = fUtil.getKeyVal(imdict['data'].split('[')[0],'EXPTIME')
                    self.logfile.write('changing '+imdict['outdata']+' parlist texptime from '+\
                                       _oldtexp+' to '+str(imdict['texptime']))
                    del _oldtexp
                    self.removeList.append(imdict['outdata'])
                    self.removeList.append(imdict['outcontext'])
                    # 'outweight' now done separately...
                    # self.removeList.append(imdict['outweight']) 
                    # be sure to remove the files if they're already there!
                    if os.path.isfile(imdict['outdata']):
                        os.remove(imdict['outdata'])
                    if os.path.isfile(imdict['outcontext']):
                        os.remove(imdict['outcontext'])
                    if os.path.isfile(imdict['outweight']):
                        os.remove(imdict['outweight'])
                if crmasks:
                    if bits:
                        # then multiply the pixel in_masks by crmasks
                        iraf.unlearn(iraf.imcalc)
                        iraf.imcalc.pixtype = 'short'
                        imcalclist = imdict['in_mask']+','+crmasks[basefits][0]
                        print imdict['in_mask'],' = ',imdict['in_mask']+' * '+crmasks[basefits][0]
                        iraf.imcalc(imcalclist,"_tmpfile.fits","im1 * im2")

                        orig_mask = 'Orig_'+imdict['in_mask']
                        os.rename(imdict['in_mask'], orig_mask)
                        os.rename("_tmpfile.fits", imdict['in_mask'])
                        self.removeList.append(orig_mask)

                        nCRzap = crmasks[basefits][1]
                        self.logfile.write(str(nCRzap)+" cosmic ray pixels masked by pyblot in "+basefits)
                        del imcalclist, orig_mask
                        iraf.flpr('imcalc')
                        iraf.flpr()
                    else:
                        # then we just copy crmask names to imdict
                        if imdict['in_mask']:  # this error should never happen
                            self.errorList.append((self.modName,"ERROR: bits=None, but in_mask exists!"))
                            raise Exception,"ERROR: bits=None, but in_mask exists!"
                        imdict['in_mask'] = crmasks[basefits][0]
                        self.logfile.write('Setting '+imdict['data']+' pyrizzle in_mask to '+imdict['in_mask'])
                        print 'Setting',imdict['data'],'pyrizzle in_mask to',imdict['in_mask']
                        nCRzap = crmasks[basefits][1]
                        self.logfile.write(str(nCRzap)+" cosmic ray pixels masked by pyblot in "+basefits)
                del imkey,tmpstring,Xname,Xver,basefits

            parlist = PyDrOb.parlist[:]
            self.parlists.append(parlist)
            del imdict
            if not separDriz:
                if self.noContext:
                    for junkdict in PyDrOb.parlist:
                        junkdict['outcontext'] = ''
                    del junkdict

                # then the first dict in list tells the whole data output story
                imdict = PyDrOb.parlist[0]
                self.logfile.write("All set for "+imdict['outdata']+' texptime: '+str(imdict['texptime']))

                # again, delete any files w/the output names so drizzle isn't confused
                if os.path.isfile(imdict['outdata']):
                    os.remove(imdict['outdata'])
                if not self.noContext:
                    if os.path.isfile(imdict['outcontext']):
                        os.remove(imdict['outcontext'])
                if os.path.isfile(imdict['outweight']):
                    os.remove(imdict['outweight'])
            if separDriz and (self.MatchDict[self.MatchDict.keys()[0]]['NumSci'] == 2):
                input = open('medianimages')
                L = input.readlines()
                i = 0
                for l in L:
                      L[i] = l[:-1]
                      i = i + 1
                print L
                input.close()
                RemoveList = []
                for imdict in PyDrOb.parlist:
                    print imdict['data']
                    if not (imdict['data'] in L):
                        RemoveList.append(imdict)
                    for imdict in RemoveList:
                        PyDrOb.parlist.remove(imdict)
                        print "Removing... %s" % (imdict['data'])

            if separDriz:
                self.shortparlists.append(PyDrOb.parlist)

            print "Beginning pydrizzle run method . . ."
            self.logfile.write('output pscale:  '+str(round(1e8*self.asecpix)/1e8))
            PyDrOb.run(save=1, build=0)   # never 'build'

            # delete 1st pass wght files, unless told otherwise
            if delwght and separDriz:
                for __tmpdict in PyDrOb.parlist:
                    if os.path.isfile(__tmpdict['outweight']):
                        os.remove(__tmpdict['outweight'])
                del __tmpdict

            self.wcslist.append(PyDrOb.observation.product.geometry.wcs)
            # drizzling all done.
            # Note: the WCS's of the the drizzle products change when the parlist shifts
            # change.  In particular, the CRPIX's of the images produced by drizzle 
            # with each other, and not with the WCS's that get appended to wcslist!
            # This is because drizzle sets the output image WCS, not xydrizzle.

            # Do some final fixes to drizzle products
            # if self.NumSci > 1 and float(xydrizzle.__version__[:3]) < 1.4:
            # break these conditions down....

            if not separDriz:
                # DANGER! this is a kludge to get the right count levels & ET's in
                # the simple drizzle products
                if self.NumSci > 1:
                    if units == 'counts':
                        # then drizzle doesn't get count levels right; divide by NumSci
                        _wtxt = "WARNING: dividing drizzled count levels by NumSci = "+str(self.NumSci)
                        self.logfile.write(_wtxt)
                        self.errorList.append((self.modName,_wtxt))
                        # divide images by NumSci
                        iraf.unlearn(iraf.imcalc)
                        operation = "im1 / "+str(float(self.NumSci))
                        iraf.imcalc(imdict['outdata'], imdict['outdata'], operation)
                        # only the 'outdata' should be divided by NumSci -- 
                        # not 'outweight' (exposure time image) or 'outcontext'
                        del operation,_wtxt
                        iraf.flpr('imcalc')
                        iraf.flpr()

                    if float(pydriz_version[:3]) < 1.4:
                        # then drizzle gets all exposure times wrong; divide by NumSci
                        _wtxt = "WARNING: dividing drizzled exptimes by NumSci = "+str(self.NumSci)
                        self.logfile.write(_wtxt)
                        self.errorList.append((self.modName,_wtxt))
                        expTimeWrong = fUtil.getKeyVal(imdict['outdata'],'EXPTIME')
                        expTimeVal = [('EXPTIME',expTimeWrong/float(self.NumSci))]
                        fUtil.fixHeader(os.path.join(self.obsFits,imdict['outdata']),expTimeVal)
                        fUtil.fixHeader(os.path.join(self.obsFits,imdict['outweight']),expTimeVal)
                        if not self.noContext:
                            fUtil.fixHeader(os.path.join(self.obsFits,imdict['outcontext']),expTimeVal)
                        del expTimeWrong, expTimeVal, _wtxt

                    else:
                        # then drizzle gets ctx and wgt exposure times wrong...
                        _wtxt = "WARNING: setting context and weight exptimes to same as sci image."
                        self.logfile.write(_wtxt)
                        self.errorList.append((self.modName,_wtxt))
                        expTimeVal = [('EXPTIME',fUtil.getKeyVal(imdict['outdata'],'EXPTIME'))]
                        # not sci image! ctx and wgt only!
                        fUtil.fixHeader(os.path.join(self.obsFits,imdict['outweight']),expTimeVal)
                        if not self.noContext:
                            fUtil.fixHeader(os.path.join(self.obsFits,imdict['outcontext']),expTimeVal)
                        del expTimeVal, _wtxt

                # Sum subtracted sky's and write result to drizzled image header
                # if the there is a single driz im, and self.alignSky non-Null
                if self.alignSky:
                    totNcomb,totSkySub = self._sumSubSkyNcombine(PyDrOb.parlist,doSky=1)
                    fUtil.fixHeader(os.path.join(os.getcwd(),PyDrOb.parlist[0]['outdata']),\
                                    [(self.alignSky,totSkySub)])
                    self.logfile.write(self.alignSky+" = "+str(totSkySub)+" set in "+\
                                       PyDrOb.parlist[0]['outdata']+" header.")
                else:
                    totNcomb,totSkySub = self._sumSubSkyNcombine(PyDrOb.parlist,doSky=0)
                # either way, we now have totNcomb
                fUtil.fixHeader(os.path.join(os.getcwd(),PyDrOb.parlist[0]['outdata']),\
                                [("NCOMBINE",totNcomb)])
                self.logfile.write("NCOMBINE = "+str(totNcomb)+" set in "+\
                                   PyDrOb.parlist[0]['outdata']+" header.")
                del totNcomb,totSkySub

            del PyDrOb

        # OK, loop over associations done!  All drizzling finished!
        # list the output for the mkMsg method
        if not separDriz:
            for list in self.parlists:
                pred_images = []
                cur_im  = ''
                next_im = ''
                sci_im = list[0]['outdata']
                if not self.noContext:
                    cxt_im = list[0]['outcontext']
                    self.contextImageList.append(cxt_im)
                wgt_im = list[0]['outweight']
                self.sciImageList.append(sci_im)
                self.weightImageList.append(wgt_im)

            for dict in list:
                next_im = dict['data'][:-7]
                if next_im != cur_im:
                    cur_im = next_im
                    pred_images.append(cur_im)
            self.outputList[sci_im] = pred_images
            if not self.noContext:
                self.outputList[cxt_im] = pred_images
            self.outputList[wgt_im] = pred_images

            # assemble all the WCS information into a list of tuples
            self.logfile.write("Making consistent WCS Info List . . .")
            wcsInfoList = self._makeWCS(template=self.sciImageList[0])

            # Now we want to fix up the headers of the final drizzle products.
            # This is done by the fUtil ufunc below.
            mod_string = self.modName + ", v"+ __version__
            date = ptime()  # ptime from pUtil
            for im in self.sciImageList:
                # with all the problems with pyfits/numpy
                # we will only try to fix the header of the file
                # warn and pass if exception
                try:    fUtil.fixDrzHeader2(im,mod_string,date,addCards=wcsInfoList)
                except Exception, err: 
                    self.logfile.write("ERROR: Call to fixDrzHeader failed on science image, "+im+\
                                       "\n\t\t DANGER!  WCS probably WRONG!  DO NOT ARCHIVE!  DANGER!\n")
                    self.logfile.write(" . . . is this a new version of pyfits . . . ?")
                    self.errorList.append((self.modName,"Call to fixDrzHeader failed on science image, "+im))
                    self.errorList.append((self.modName,"WCS of "+im+" likely WRONG!  DO NOT ARCHIVE!"))
                    self.errorList.append((self.modName,str(err)))
                    continue

            for im in self.contextImageList:
                try: 
                    fUtil.fixDrzHeader2(im,mod_string,date,file_type="CTX",addCards=wcsInfoList)
                except Exception, err: 
                    self.logfile.write("ERROR: Call to fixDrzHeader failed on context image, "+im)
                    self.logfile.write("Probabaly another flakey version of pyfits....")
                    self.errorList.append((self.modName,"Call to fixDrzHeader failed on context image, "+im))
                    self.errorList.append((self.modName,str(err)))
                    continue

            for im in self.weightImageList:
                try:
                    fUtil.fixDrzHeader2(im,mod_string,date,file_type="WGT",addCards=wcsInfoList)
                except Exception, err: 
                    self.logfile.write("ERROR: Call to fixDrzHeader failed on weight image, "+im)
                    self.logfile.write("Probabaly another flakey version of pyfits....")
                    self.errorList.append((self.modName,"Call to fixDrzHeader failed on weight image, "+im))
                    self.errorList.append((self.modName,str(err)))
                    continue

            # get the medriz_?_.fits images into the outputList. not archiveable.
            medriz_list = glob.glob("medriz_*.fits")
            if medriz_list:
                for file in medriz_list:
                    self.outputList[file] = []   # can't determine predecessors at this point but who cares?

        os.chdir(curdir)
        self.logfile.write('returning to directory '+curdir)
        return

    def makeFlagImage(self):
        """turns the weight images produced by drizzle into flag images that
            SExtractor will use.
        """
        if not self.weightImageList:
            errtxt="No Weight Images present."
            self.errorList.append((self.modName,errtxt))
            raise Exception, errtxt
            # reset flag image list
        while self.flagImageList:
            del self.flagImageList[0]

        curdir = os.getcwd()
        os.chdir(self.obsFits)
        for im in self.weightImageList:
            try:
                wgtfits = pyfits.open(im)
            except Exception,err:
                self.errorList.append((self.modName,str(err)))
                raise Exception,err

            if len(wgtfits) > 1:
                self.errorList.append((self.modName,"image file is not simple fits."+im))
                raise Exception,"image file is not simple fits."+im

            # build flag image name
            flgfile = im.split("_drz")[0]+'_FLAG.fits'
            self.flagImageList.append(flgfile)
            self.outputList[flgfile] = [im]
            
            # create and initialize the new pyfits object
            flgfits = pyfits.HDUList()
            flgfits.append(pyfits.PrimaryHDU())
            try:
                del flgfits[0].header.ascard["EXTEND"]
            except KeyError:
                pass
            flgfits[0].header = wgtfits[0].header
            flgfits[0].data = numpy.logical_not(wgtfits[0].data).astype(numpy.int16)
            wgtfits.close()
            flgfits[0].header.update('BITPIX',16)            
            flgfits[0].header.update('FILENAME',flgfile)
            flgfits[0].header.update('FILETYPE','FLG')

            # close (write out) the flag image
            flgfits.writeto(flgfile)
            
            self.logfile.write('Made flag image '+flgfile)
            del wgtfits, flgfits

        os.chdir(curdir)
        return

    def makeRmsImage(self):
        """turns the drizzled weight images into RMS images that
            SExtractor will use.
        """
        self.logfile.write("starting makeRmsImage . . .")
        if not self.weightImageList:
            errtxt="No Weight Images present."
            self.errorList.append((self.modName,errtxt))
            raise Exception, errtxt

        # reset rms image list
        while self.rmsImageList:
            del self.rmsImageList[0]

        curdir = os.getcwd()
        os.chdir(self.obsFits)
        for im_wgt in self.weightImageList:
            im_sci = im_wgt[:-12]+'.fits'
            if im_sci not in self.sciImageList:
                errtxt = 'makeRmsImage: '+im_sci+' not in sciImageList[]!'
                self.errorList.append((self.modName,errtxt))
                self.logfile.write(errtxt)
            try:
                wgtfits = pyfits.open(im_wgt)
                scifits = pyfits.open(im_sci)
            except:
                self.errorList.append((self.modName,"Cannot make a FITS object out of file "+im_wgt))
                raise Exception,"Cannot make a FITS object out of file "+im_wgt
            if len(wgtfits) > 1 or len(scifits) > 1:
                self.errorList.append((self.modName,"image file is not simple fits."))
                raise Exception,"image file is not simple fits."

            # build rms image name and open as a new file.
            rmsfile = im_wgt.split("_drz")[0]+'_RMS.fits'
            self.rmsImageList.append(rmsfile)
            self.outputList[rmsfile] = [im_wgt]
            
            # make new fits obj and copy WGT/SCI hdr/data to RMS image initially
            rmsfitsobj = pyfits.HDUList()
            rmsfitsobj.append(pyfits.PrimaryHDU())
            try:
                del rmsfitsobj[0].header.ascard["EXTEND"]
            except KeyError:
                pass
            rmsfitsobj[0].header = wgtfits[0].header
            rmsfitsobj[0].data   = scifits[0].data
            if os.path.isfile(rmsfile):
                os.remove(rmsfile)
            rmsfitsobj.writeto(rmsfile)
            del rmsfitsobj

            # reopen the rms image for editing.
            rmsfits = pyfits.open(rmsfile,'update')

            # ratio of default to specified output scales
            area_ratio = (self.asecpix / self.origscale)**2
            if abs(1-area_ratio) < 1e-4: area_ratio = 1
            self.logfile.write('Using area_ratio = %.6f in makeRmsImage' %(area_ratio))

            skyval  = scifits[0].header.get('ALIGNSKY')  # this rescaled below
            exptime = scifits[0].header.get('EXPTIME')
            Ncombed = scifits[0].header.get('NCOMBINE')
            if Ncombed == None:
                errtxt='Error: NCOMBINE not in '+im_sci+' header. _sumSubSkyNcombine() not run?'
                self.logfile.write(errtxt)
                self.errorList.append((self.modName,errtxt))
                raise Exception,errtxt
        
            gain,rn = pyblot._gain_rn(scifits, self.logfile, ext=0)
            # if not told to use header gain, then use 1.0 (data in electrons)
            if not self.hdrGain:
                gain = 1.0
            self.logfile.write(im_sci+":  gain,rn = "+str(gain)+","+str(rn)+\
                               "  NCOMBINE = "+str(Ncombed)+"  EXPTIME = "+str(exptime))
            if not exptime:
                raise Exception,"No EXPTIME in "+im_sci
            if (skyval == None or skyval < 1):
                warntxt = 'WARNING: found ALIGNSKY of '+str(skyval)+' in '+im_sci+\
                          ' : RMS image may be in WRONG!'
                self.logfile.write(warntxt)
                self.errorList.append((self.modName,warntxt))
                del warntxt
                if skyval == None: skyval=0

            skyval *= area_ratio
            
            # 1. construct variance from sky, sci[], wght[]
            # 2. clip zeros/negatives and infinities to values that will work w/sqrt(),
            #     and so that sqrt(val) = 2e30 > 1e30, => zero weight in SExtractor;
            #     have to work in Float64 to avoid Inf's.
            # 3. take sqrt() and cast as float32
            # 4. tidy header
            # 5. write it out
        
            readVariance = Ncombed*(rn/gain)*(rn/gain)
            self.logfile.write("total read variance = "+str(readVariance)+" for "+im_sci)
            if self.suppInstVar:
                # supplement factor for reference bais subtraction, etc
                # extra var per sec for dark subtraction, repaired cosmic rays, etc.
                totInstVar =  (_rnVarSuppFac_ * readVariance) + (_exVarPerSec_ * exptime)
                self.logfile.write("adjusted instrumental variance = "+str(totInstVar)+" for "+im_sci)
            else:
                totInstVar = readVariance

            totInstVar *= area_ratio
        
            # maybe doing arithmetic in two steps will help conserve memory...
            newDat  = ((skyval + scifits[0].data.astype(numpy.float64))/gain + totInstVar) * (exptime * area_ratio)
            # newDat[] is now variance *in counts* times maximum expTime; divide by expTime map...
            newDat /= wgtfits[0].data
            scifits.close()
            wgtfits.close()
            del scifits, wgtfits, im_wgt, im_sci, readVariance, totInstVar, area_ratio
        
            ## now fix up problem values...
            newDat = numpy.where(numpy.logical_or(numpy.greater_equal(newDat,1e60),\
                                                      numpy.less_equal(newDat,0.)),4e60,newDat)
            rmsfits[0].data = numpy.sqrt(newDat).astype(numpy.float32)

            # a few token updates to the header, then write it out
            rmsfits[0].header.update('FILENAME',rmsfile)
            rmsfits[0].header.update('FILETYPE','RMS')
            rmsfits.close()
            self.logfile.write('Made rms image '+rmsfile)
            del newDat, rmsfile, rmsfits

        os.chdir(curdir)
        return
    
    def fixAstrometry(self,obs):
        """Method to call the astrometer module's stuff and do the astrometric corrections.
        The astrometer constructor needs to receive the obs (i.e. DataSet) object.
        """

        print "Now correcting astrometric zeropoint..."
        astrom=astrometer.gscMatchup(obs)
        
        try:
            rval = astrom.findAstromCorrs()
        except astrometer.WebQueryError,err:
            warntxt = "Caught a WebQueryError. Astrometric matchup not successful."
            print warntxt
            self.logfile.write(warntxt)
            self.logfile.write(str(err))
            self.errorList.append((self.modName,warntxt))
            self.errorList.append((self.modName,str(err)))
            raise astrometer.WebQueryError,err
        
        if not rval:
            print "Astrometric matchup successful."
            self.logfile.write("Astrometric matchup successful.")
            self.logfile.write("Applying corrections.")
            astrom.applyCorrs()
        return
    

    def writeXml(self):
        """mark up products as xml."""

        curdir = os.getcwd()
        os.chdir(self.obsFits)

        if self.sciImageList:
            for im in self.sciImageList:
                file = xmlUtil.markupImage(im,dataset=self.obsName)
                if file not in self.outputList.keys():
                   self.outputList[file] = [im]
        if self.contextImageList:
            for im in self.contextImageList:
                file = xmlUtil.markupImage(im,dataset=self.obsName)
                if file not in self.outputList.keys():
                    self.outputList[file] = [im]
        if self.weightImageList:
            for im in self.weightImageList:
                file = xmlUtil.markupImage(im,dataset=self.obsName)
                if file not in self.outputList.keys():
                    self.outputList[file] = [im]
        if self.flagImageList:
            for im in self.flagImageList:
                file = xmlUtil.markupImage(im,dataset=self.obsName)
                if file not in self.outputList.keys():
                    self.outputList[file] = [im]
        if self.rmsImageList:
            for im in self.rmsImageList:
                file = xmlUtil.markupImage(im,dataset=self.obsName)
                if file not in self.outputList.keys():
                    self.outputList[file] = [im]
        os.chdir(curdir)
        return


    def mkMsg(self):
        """create and write module level message for this class.
        Most of this is just compiling the info. meta in a dictionary
        of lists where each list is a list of tuples describing the
        tag lines for the particular section of the message.  This tuple 
        format conforms to that used by the xmlMessage class which is
        modeled on basic python argument passing, i.e. (key,*value,**attr).
            outscale,pixfrac now converted to string - GRM following Magee
        """
        self.meta = {}
        self.meta['module']= []
        self.meta['meta']  = []
        self.meta['input'] = []
        self.meta['output']= []
        self.meta['errorlist'] = []

        self.meta['module'].append(('module','name='+self.modName,'version='+__version__,'dataset='+self.obsName))
        self.meta['module'].append(('root',self.root))
        self.meta['meta'].append(('meta',))
        self.meta['meta'].append(('configuration',))
        self.meta['meta'].append(('parameter','name=outscale',str(self.asecpix)))
        self.meta['meta'].append(('parameter','name=pixfrac',str(self.pixfrac)))
        self.meta['meta'].append(('parameter','name=kernel',self.kernel))
        self.meta['meta'].append(('depend',))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','python'))
        self.meta['meta'].append(('version',pyversion.split()[0]))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','xydrizzle'))
        self.meta['meta'].append(('version',pydriz_version))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','drizzle'))
        self.meta['meta'].append(('version',drversion))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','pyfits'))
        self.meta['meta'].append(('version',pyfits.__version__.split()[0]))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','numpy'))
        self.meta['meta'].append(('version',numpy.__version__))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','pyblot'))
        self.meta['meta'].append(('version',pyblot.__version__))        
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','astrometer'))
        self.meta['meta'].append(('version',astrometer.__version__))
                                 
        if self.errorList:
            self.meta['errorlist'].append(('errorlist',))
            for pkg,err in self.errorList:
                self.meta['errorlist'].append(('erroritem',err,'frompkg='+pkg))

        # input section
        self.meta['input'].append(('input',))
        for f in self.inputList:
            if string.find(f,"_asn") == -1:
                self.meta['input'].append(('file','type=image/x-fits'))
                self.meta['input'].append(('name',os.path.join("Images",f)))
            else:
                self.meta['input'].append(('file','type=image/x-fits'))
                self.meta['input'].append(('name',os.path.join("Images",f)))

        # output section
        if self.outputList:
            self.meta['output'].append(('output',))
        for f in self.outputList.keys():
            if string.find(f,".xml") == -1:
                self.meta['output'].append(('file','type=image/x-fits'))
                self.meta['output'].append(('name',os.path.join("Images",f)))
                for pred in  self.outputList[f]:
                    self.meta['output'].append(('predecessor',os.path.join("Images",pred)))
            else:
                self.meta['output'].append(('file','type=text/xml'))
                self.meta['output'].append(('name',os.path.join("Images",f)))
                for pred in  self.outputList[f]:
                    self.meta['output'].append(('predecessor',os.path.join("Images",pred)))
        

        # pass this dictionary to the class pMessage...
        msgFile = os.path.join(self.messagedir,self.modName+"_module.xml")
        mmsg = pMessage(self.meta)
        mmsg.writeMsg(msgFile)
        return


#------------------------------------ private stuff -------------------------------------#


    def _makeWCS(self, template=None):
        """ Make a list of tuples containing all the WCS and closely
        related info to be written to each output image header.
        """
        if template:  tempfits = pyfits.open(template)
        else:         tempfits = pyfits.open(self.sciImageList[0])
        temphdr = tempfits[0].header

        dec0 = temphdr['CRVAL2']
        wcstuples = [('CRPIX1', temphdr['CRPIX1']),
                     ('CRPIX2', temphdr['CRPIX2']),
                     ('CRVAL1', temphdr['CRVAL1']),
                     ('CRVAL2', temphdr['CRVAL2'])]
        tempfits.close()
        del temphdr,tempfits

        for L in self.reflogfile_add:
            self.logfile.write(L)
            
        if self.refcdmat['EXTREF']:
        	cdmat = self.refcdmat

        else:
        	PA_final = self.refotherkeys['PA_FINAL']
        	cdmat = matutil.makeCDdict(PA_final, self.asecpix)

        wcstuples.append(('CD1_1', round(1e12*cdmat['CD1_1'])/1.0e12))
        wcstuples.append(('CD1_2', round(1e12*cdmat['CD1_2'])/1.0e12))
        wcstuples.append(('CD2_1', round(1e12*cdmat['CD2_1'])/1.0e12))
        wcstuples.append(('CD2_2', round(1e12*cdmat['CD2_2'])/1.0e12))          

        for key in self.refotherkeys.keys():
            wcstuples.append((key, self.refotherkeys[key]))
    
        return wcstuples


    def _sumSubSkyNcombine(self,parlist, doSky):
        """ Sum the sky values in the input images and put the
        total in the output drizzled image header.
        """
        imSkyDict = {}
        outData = parlist[0]['outdata']        
        totNcomb = 0

        for parDict in parlist:
            if parDict['outdata'] != outData:
                raise Exception,"Images in parlist must have same outdata for _sumSubSky()."
            mefIm  = string.split(parDict['data'],'[')[0]
            NameVerString = string.split(string.split(parDict['data'],'[')[1],']')[0]
            Xname = string.split(NameVerString,',')[0]
            Xver  = int(string.split(NameVerString,',')[1])

            # get the sky and ncombine values
            ff  = pyfits.open(mefIm)
            if doSky:
                newSky = ff[Xname,Xver].header.get(self.alignSky) 
            icomb  = ff[Xname,Xver].header.get("NCOMBINE")
            if icomb == None or icomb < 1:
                self.logfile.write("WARNING: did not find NCOMBINE in "+mefIm+" ("+Xname+","+str(Xver)+").  Assuming N=1.")
                totNcomb += 1
            else:
                totNcomb += icomb
            ff.close()
            del ff,icomb

            if not doSky:
                continue
            
            if newSky == None:
                newSky = 0.0
                errtxt = "Warning: "+self.alignSky+" param not found in "+\
                         mefIm+" "+NameVerString
                self.errorList.append((self.modName,errtxt))
                self.logfile.write(errtxt)

            # append this extension sky to a list for this image
            try:
                imSkyDict[mefIm].append(newSky)
            except:
                imSkyDict[mefIm] = []
                imSkyDict[mefIm].append(newSky)

            del mefIm,NameVerString,Xname,Xver,newSky

        # end of list over input images (parDict's).
        # just normalize and return Ncombine if not doing sky
        totNcomb /= self.NumSci
        if not doSky:
            return (totNcomb,None)

        # ok, now we can just loop over mefIms, average the ext Sky's
        # for each mefIm and sum these averages for a total sky

        skyTot = 0.0
        for mefIm in imSkyDict.keys():
            imSky = 0.0
            skyList = imSkyDict[mefIm]
            for extSky in skyList:
                imSky += extSky
            imSky = imSky / len(skyList)
            skyTot = skyTot + imSky

        return (totNcomb,skyTot)


    def _rotateRect(self,_maxAng,nx,ny):
        """ rotate a rectangle and return ratio of new x,y sizes to old"""
        theta = _maxAng * math.pi/180.
        a     = ny/2.0
        b     = nx/2.0
        c     = math.sqrt(a*a + b*b)
        alpha = math.atan(b/a)

        halfwidth  = c*math.sin(alpha + theta)
        halfheight = c*math.cos(alpha - theta)

        xfac = halfwidth/b
        yfac = halfheight/a

        return (xfac,yfac)

    def printParlist(self, parlist):
        """ takes a PyDrizzle parlist and prints out the contents"""
        for i in range(len(parlist)):
            print 'Member ',i,':'
            for key in parlist[i].keys():
                print ('  %-10s :   ' % key), parlist[i][key]
            print '\n'
        return

    def help(self):
        print self.__doc__
        print self.run_drizzle.__doc__
        print self.run_all.__doc__
        return
    
