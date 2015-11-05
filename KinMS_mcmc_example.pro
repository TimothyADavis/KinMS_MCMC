function mkkinms_model,param,_extra=_extra,filename=filename,ret_inclouds=ret_inclouds,regen=regen,ret_gassigma=ret_gassigma
    COMMON KinMS_COMMON, fdata,obspars
  ;;;;;;;;;;;;;;;
  ;;; Make your model here. Takes the inputs from obspars, param, and
  ;;; returns the model cube.
  ;;;;;;;;;;;;;;;
    !EXCEPT = 0
  svrad=n_elements(obspars.vrad)                 ;;; size of the velocity vector
  x=obspars.vrad                                 ;;; radius vector for sbprofile
  fx=gaussian(x,[1,param[(7+(svrad)+1)],3])      ;;; Ring morphology
  x1=obspars.vrad                                ;;; radius vector for velocity  
  vel=param[7:(7+(svrad)-1)]                     ;;; get the velocity profile from the inputs
  vel*=param[(7+(svrad))]                        ;;; multiply this by the fitted normalisation
  KinMS,obspars.xsize,obspars.ysize,obspars.vsize,obspars.dx,obspars.dy,obspars.dv,obspars.beamsize,param[2],velrad=x1,velprof=vel,nsamps=obspars.nsamps,cubeout=fsim,posang=param[1],intflux=param[0],gassigma=param[3],_extra=_extra,ra=obspars.ra,dec=obspars.dec,phasecen=[param[4],param[5]],voffset=param[6],filename=filename,inclouds=inclouds,sbrad=x,sbprof=fx,/fixseed ;;; run the model.
  return,fsim
end
FUNCTION mkkinms_likelihood,param,model,_EXTRA=ex
  COMMON KinMS_COMMON, fdata,obspars
  ;;;;;;;;;;;;;;;
  ;;; This function calls the model, then evaulates the log likelihood.
  ;;;;;;;;;;;;;;;
   !EXCEPT = 0
  modout = call_FUNCTION(model,param,_EXTRA=ex)
   w=where(finite(fdata) and finite(modout),nvalid)
  IF nvalid GT 1 THEN BEGIN      
     chisq = total(((fdata[w]-modout[w])^2)/((obspars.rms)^2))
     loglike=0.5*chisq
  ENDIF ELSE BEGIN
     loglike = 1e31
  ENDELSE 
  IF loglike EQ 1e31 THEN stop
  return,loglike
END
pro KinMS_mcmc_example
  COMMON KinMS_COMMON, fdata,obspars
  ;;;;;;;;;;;;;;;
  ;;; This is a test procedure for KinMS_mcmc, setting up everything
  ;;; for fitting the molecular ring in NGC4324, based on CARMA observations
  ;;; of Alatalo et al., 2012
  ;;;;;;;;;;;;;;;  
  

  ;;; Here we define the velocity curve, but normalised to one. We then
  ;;; rescale later with a single paramter for vmax. One could also
  ;;; fit the whole curve using the formalism below- but this is
  ;;; enough here.
  vrad=(findgen(64))
  vcirc=interpol([0,1,1,1],[0.01,1,10,200],vrad)
  ;;;

  ;;; Here we set the default model parameters (and ranges). Set the
  ;;; minimum and maximum to the same as the guess to NOT fit this parameter
  intflux= 27.2        ;; Best fitting total flux
  minintflux= 27.2     ;; lower bound total flux
  maxintflux= 27.2     ;; upper bound total flux
  gassigma=8.0         ;; Best fitting gas sigma
  mingassigma=8.0      ;; Min sigma
  maxgassigma=8.0     ;; Max sigma
  posang= 235.0        ;; Best fit posang.
  minposang= 235.      ;; Min posang.
  maxposang= 235.      ;; Max posang.
  inc=65.              ;; degrees
  mininc=50.           ;; Min inc
  maxinc=89.           ;; Max inc
  voffset=0.0          ;; Best fit velocity centroid
  minvoffset=0.0       ;; min velocity centroid
  maxvoffset=0.0       ;; max velocity centroid
  centx=0.0            ;; Best fit x-pos for kinematic centre
  mincentx=0.0         ;; min cent x
  maxcentx=0.0         ;; max cent x
  centy=0.0            ;; Best fit y-pos for kinematic centre
  mincenty=0.0         ;; min cent y
  maxcenty=0.0         ;; max cent y
  vflat =  162.0       ;; vflat
  min_vflat=10         ;; Lower range vflat
  max_vflat=300        ;; Upper range vflat
  r_ring=20.0          ;; ring radius
  minr_ring=2.         ;; Lower ring radius
  maxr_ring=40.        ;; Upper ring radius
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  
  ;;;; Data cube operations ;;;;;
  fdata=readfits("NGC4324.fits",hdr,/silent) ; read in datacube to fit-  2^n sizes are MUCH quicker.
  RMS=robust_sigma(fdata[*,*,0:2])
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



  names=["intflux","posang","inc","gassigma",'centx','centy','voffset',replicate('Vel',n_elements(vrad)),"vflat","r_ring"]
  ;;; best guesses for the parameters
  guesses=[intflux,posang,inc,gassigma,centx,centy,voffset,vcirc,vflat,r_ring]
  ;;;
  ;;; minimum values for the parameters. Set this and maxpar to the same as the guess to not fit this parameter
  minpar=[minintflux,minposang,mininc,mingassigma,mincentx,mincenty,minvoffset,vcirc,min_vflat,minr_ring]
  ;;;
  ;;; maximum values for the parameters. Set this and minpar to the same as the guess to not fit this parameter
  maxpar=[maxintflux,maxposang,maxinc,maxgassigma,maxcentx,maxcenty,maxvoffset,vcirc,max_vflat,maxr_ring]
  ;;;


  ;;;; Now set up the structure for the fitting
  param = { name :       names, $
            value :      guesses, $
            max :        maxpar, $
            min :        minpar, $
            knob :       replicate(1.0,n_elements(names)), $
           knobchange :       replicate(1.0,n_elements(names)), $
           changeable : (minpar ne maxpar) }
  ;;;;
  ;;;; and the structure to contain the observational parameters
  obspars={xsize:xtot,ysize:ytot,vsize:vtot,dx:dx,dy:dy,dv:dv,beamsize:beam,nsamps:5e5,ra:0.d,dec:0.d,vrad:vrad,rms:RMS,intscat:1.d}
  ;;;;

   
;;;; generate best guess model and display it ;;;;
  !p.multi=[0,1,1]

  t=systime(/sec)
  for i=0,19 do fsim= mkkinms_model(param.value)
  print,"model generation time (new):",(systime(/sec)-t)/20.
  
  lookatplot_mkplot,fdata,fsim,param.value,obspars,xtitle='Position (")',ytitle=textoidl('Velocity (km s^{-1})'),charsize=1.5 ;,/eps
  print,"Initial guess model and data"
  print,"[press any key to continue]."
  meep=""
  read,meep
;;;;

;;;; Call the MCMC routines - FDATA AND OBSPARS ARE PASSED IN COMMON BLOCK
  
  best=KinMS_mcmc(param,finaloutput="NGC4324_bestfit",iters=100)
;;;;
end

