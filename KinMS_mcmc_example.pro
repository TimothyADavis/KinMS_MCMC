function mkkinms_model,param,_extra=_extra,filename=filename,ret_inclouds=ret_inclouds,regen=regen,ret_gassigma=ret_gassigma
    COMMON KinMS_COMMON, fdata,obspars
  ;;;;;;;;;;;;;;;
  ;;; Make your model here. Takes the inputs from param (and the
  ;;; common block) and
  ;;; returns the model cube.
  ;;;;;;;;;;;;;;;
    !EXCEPT = 0
    x=obspars.vrad                                                    ;;; radius vector for sbprofile
    fx=gaussian(x,[1,param[7],param[8]])                                    ;;; Ring morphology
    x1=obspars.vrad                                                   ;;; radius vector for velocity  
    vel=interpol([0,1,1,1],[0.01,1,10,200],x1)*param[6]               ;;; impose a flat velocity profile with a fitted normalisation
    KinMS,obspars.xsize,obspars.ysize,obspars.vsize,obspars.dx,obspars.dy,obspars.dv,obspars.beamsize,param[2],velrad=x1,velprof=vel,nsamps=obspars.nsamps,cubeout=fsim,posang=param[1],intflux=param[0],gassigma=1.,_extra=_extra,ra=obspars.ra,dec=obspars.dec,phasecen=[param[3],param[4]],voffset=param[5],filename=filename,inclouds=inclouds,sbrad=x,sbprof=fx,/fixseed ;;; run the model.
    return,fsim
 end
pro KinMS_mcmc_example
  COMMON KinMS_COMMON, fdata,obspars
  ;;;;;;;;;;;;;;;
  ;;; This is a test procedure for KinMS_mcmc, setting up everything
  ;;; for fitting the molecular ring in NGC4324, based on CARMA observations
  ;;; of Alatalo et al., 2013
  ;;;;;;;;;;;;;;;  
  

  ;;; Here we define the radius vector.
  vrad=(findgen(64))
  ;;;

  ;;; Here we set the default model parameters (and ranges). Set the
  ;;; minimum and maximum to the same as the guess to NOT fit this parameter
  intflux= 30        ;; Best fitting total flux
  minintflux= 1.             ;; lower bound total flux
  maxintflux= 50.            ;; upper bound total flux
  posang= 250.              ;; Best fit posang.
  minposang= 200.           ;; Min posang.
  maxposang= 265.           ;; Max posang.
  inc= 75.                   ;; degrees
  mininc=50.                ;; Min inc
  maxinc=89.                ;; Max inc
  centx=0.0                 ;; Best fit x-pos for kinematic centre
  mincentx=-5.0             ;; min cent x
  maxcentx=5.0              ;; max cent x
  centy=0.0               ;; Best fit y-pos for kinematic centre
  mincenty=-5.0             ;; min cent y
  maxcenty=5.0              ;; max cent y
  voffset= 0.0               ;; Best fit velocity centroid
  minvoffset=-20.0          ;; min velocity centroid
  maxvoffset=+20.0          ;; max velocity centroid
  vflat =  100.41546             ;; vflat
  min_vflat=10              ;; Lower range vflat
  max_vflat=300             ;; Upper range vflat
  r_ring=10.0               ;; ring radius
  minr_ring=10.              ;; Lower ring radius
  maxr_ring=25.              ;; Upper ring radius
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

rrad=3.
minrrad=0.
maxrrad=6.

  
  ;;;; Data cube operations ;;;;;
  fdata=readfits("NGC4324.fits",hdr,/silent) ; read in datacube to fit-  2^n sizes are MUCH quicker.
  fdata=fdata[50-32:50+31,50-32:50+31,*]
  RMS=stdev(fdata[*,*,0:1]) 
  sfdata=size(fdata)                   ; calculate the size of the resultant datacube
  fdata[where(finite(fdata) eq 0)]=0.0 ;; set NANs to zero
  ;;;;
  
  ;;;; setup the observation parameters, read from datacube header
  beam=[sxpar(hdr,'BMAJ')*3600.,sxpar(hdr,'BMIN')*3600.,sxpar(hdr,'BPA')] ; setup beam
  dv=sxpar(hdr,'CDELT3')/1e3
  vtot=dv*sfdata[3]
  dx=abs(sxpar(hdr,'CDELT1')*3600.)
  dy=abs(sxpar(hdr,'CDELT2')*3600.)
  xtot=dx*sfdata[1]
  ytot=dy*sfdata[2]
  ;;;;

  

  ;; Names of the fitted parameters
  names=["intflux","posang","inc",'centx','centy','voffset',"vflat","r-ring","rrad"]
  ;;; best guesses for the parameters
  guesses=[intflux,posang,inc,centx,centy,voffset,vflat,r_ring,rrad]
  ;;; minimum values for the parameters. Set this and maxpar to the same as the guess to not fit this parameter
  minpar=[minintflux,minposang,mininc,mincentx,mincenty,minvoffset,min_vflat,minr_ring,minrrad]
  ;;; maximum values for the parameters. Set this and minpar to the same as the guess to not fit this parameter
  maxpar=[maxintflux,maxposang,maxinc,maxcentx,maxcenty,maxvoffset,max_vflat,maxr_ring,maxrrad]
  ;;; How precisely do you need to know a give parameter? Increasing
  ;;; these numbers lets the chain coverge more easily. If set to
  ;;; zero, then 1% of the current best fit value is used.
  precision=[0.3,1.0,1.0,dx/3.,dy/3.,dv/3.,dv/3.,dx/3.,dx/3.]


  ;;;; Now set up the structure for the fitting
  param = { name :       names, $
            value :      guesses, $
            max :        maxpar, $ 
            min :        minpar, $
            knob :       replicate(1.0,n_elements(names)), $
            precision :       precision, $
            changeable : (minpar ne maxpar) }
  
;;;; and the structure to contain the observational parameters
  obspars={xsize:xtot,ysize:ytot,vsize:vtot,dx:dx,dy:dy,dv:dv,beamsize:beam,nsamps:1e4,ra:0.d,dec:0.d,vrad:vrad,rms:RMS,intscat:1.d}
;;;;

;;;; generate first guess model to display ;;;;
  !p.multi=[0,1,1]
  fsim= mkkinms_model(param.value)
  lookatplot_mkplot,fdata,fsim,param.value,obspars,xtitle='Position (")',ytitle=textoidl('Velocity (km s^{-1})'),charsize=1.5
  print,"Initial guess model and data (purposefully bad)"
  print,"[press any key to continue]."
  meep=""
  read,meep
;;;;

;;;; Call the MCMC routines - FDATA AND OBSPARS ARE PASSED IN COMMON BLOCK
  best=KinMS_mcmc(param,finaloutput="NGC4324_bestfit",iters=30000,outputll=outputll,silent=silent)
  print,"Reduced Chi-sqr:",(abs(outputll[0]*2.)+n_elements(fdata))/float(n_elements(fdata)-total(param.changeable))
  meep=""
  read,meep
;;;;

;;;; generate best model to display ;;;;
  fsim= mkkinms_model(best)
  print,"Best fit model and data"
  lookatplot_mkplot,fdata,fsim,param.value,obspars,xtitle='Position (")',ytitle=textoidl('Velocity (km s^{-1})'),charsize=1.5
  meep=""
  read,meep
  mk_chanmap,fdata,fsim,hdr,rms,[5,4],rmsfac=3,chans2do=[1,20]
  print,"Best fit channel maps"
;;;;
end

