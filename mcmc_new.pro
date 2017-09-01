;##############################################################################
;
; Copyright (C) 2015, Timothy A. Davis
; E-mail: DavisT -at- cardiff.ac.uk
;
; Updated versions of the software are available from github:
; github.com/TimothyADavis/KinMS_MCMC
; 
; If you have found this software useful for your research,
; I would appreciate an acknowledgment to the use of the
; "KINematic Molecular Simulation (KinMS; Davis et al., 2013)
; MCMC routines of Davis (in prep).".
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;##############################################################################
;+
; NAME:
;    MCMC_new
;
; PURPOSE:
;    Construct MCMC chains for parameter fitting, using an adaptive
;    step size MCMC algorithm with gibbs sampling. Tested for use with
;    the KinMS routines.
;
; CALLING SEQUENCE:

; mcmc_new,param,iters,numatonce,burn,nchains,[model=,likelihood=,plotchains=,silent=,covar=,log=,outputll=,outputvalue=]
;
; REQUIRED INPUTS:
;    PARAM       IDL structure containing the fitting details.
;                The structre requires the following fields:
;                param = { name       : ["Param1",Param2"], $  
;                        value      : [1.0D,2.0D], $   ;; Starting guesses
;                        max        : [10.0D,10.0D], $ ;; Maximum value
;                        min        : [0.0D,0.0D], $   ;; Mimimum value
;                        precision  : [0.1D,0.1D],$    ;; Convergance threshold
;                        changeable : [1,1] }          ;; Fit this?
;    ITERS      Number of iterations for the final results chain. Also
;               the number of iterations allowed to find covergance
;    NUMATONCE  How often to check for convergence.
;    BURN       Length of the burn in phase
;    NCHAINS    How many independant chains to try in order to find
;               global minima.
;   
; OPTIONAL INPUTS:
;    MODEL        Name of the function that generates the model you wish
;                 to fit. Must accept an array of guesses in the form of
;                 param.value. DEFAULT: "model"
;    LIKELIHOOD   Name of the function that generates the likelihood
;                 value for your model. Must accept param.value and 
;                 the name of the model function. Must return the
;                 log likelihood * (-1). Can also accept other
;                 variables, including COVAR and LOG.
;    PLOTCHAINS   If set, plot the chains as they converge. Turn off for speed.
;    SILENT       Supress information messages.
;    COVAR        If set, calculate Chi_2 using the inverse of the covariance
;                 matrix, which must be included in the
;                 "KinMS_COMMON_covar" common block and be called
;                 invmult. If not set caculate the Chi_2 assuming all variables are
;                 independant.
;    LOG          If set, the LIKELIHOOD function returns log
;                 likihoods. If not set, assumes linear likelihoods passed.
;
;  OPTIONAL OUTPUTS:
;    OUTPUTLL     Returns the negative of the likihood of each point
;                 in the final chain, after the burn in.
;    OUTPUTVALUE  Returns the value of each step in the final chain
;                 after the burn in. Used to get the parameter PDFs.
;   
; SIDE EFFECTS:
;    None.
;
; RESTRICTIONS:
;    Requires IDL 5.0 or higher (square bracket array
;    syntax). Requires the latest IDL astrolib.
;    
;
; USAGE EXAMPLE:
;    A simple usage example is given in the KinMS_MCMC_test.pro file
;    
;
; MODIFICATION HISTORY:
; DavisT@cardiff.ac.uk
; http://www.astro.cardiff.ac.uk/pub/Tim.Davis/
; 30/11/2015. v0.5 - Cleaned up code and documentation ready for
;                    internal use.  
;
;##############################################################################


function mcmc_new_stepper,values,change,knob,min,max,seed
  ;; Makes a step in the requested dimensions, within the specified range
  nchange=n_elements(change)
  newvals= double(values)
  acceptable=0
  cnt=0
  while acceptable eq 0 do begin   
     newvals[change]= values[change]+(randomn(seed,nchange,/double)*knob[change])
     if total(newvals[change] ge min[change]) eq nchange and total(newvals[change] le max[change]) eq nchange then acceptable++
     cnt++
     if cnt gt 1000 then stop,"Cannot find a good step. Check input parameters"
  endwhile
  return,double(newvals)
end

function mcmc_new_take_step,values,bestll,change,knob,min,max,seed,newll,accept,likelihood,model,holdknob=holdknob,tried=tried,log=log,covar=covar
  ;; Take a step in the MCMC chain, evalute the LL and decide if the
  ;; step is accepted or not

  ;;;;;; adaptive step parameter. Make the difference in these
  ;;;;;; parameters larger or smaller to change chain mixing and
  ;;;;;; convergence rate
  dec=0.95
  inc=1.05
  ;;;;;;;
  if n_elements(log) eq 0 then log=1
  accept=1
  newval=mcmc_new_stepper(values,change,knob,min,max,seed)
  newll = call_FUNCTION(likelihood,newval,model,covar=covar,_extra=ex)
 ; print,newll,bestll
  if keyword_set(log) then begin
    ;; Fit in log likihood. (chi2_new - chi2_best)*0.5
     ratio=  newll-bestll
     randacc=alog(randomu(seed,/double))
  endif else begin
     ;; Fit in likelihood. exp(chi2_new*0.5)/exp(chi2_best*0.5)
     ratio= newll/bestll
     randacc=randomu(seed,/double)
  endelse
  tried=newval
 ; print,newll,bestll
  if randacc gt ratio then begin
     ;; dont accept if the likelihood ratio is less than the random
     ;; acceptance threshold. Always accepts if the guess is better
     accept=0
     newval = values
     newll = bestll
     if not keyword_set(holdknob) then knob[change]*=dec 
  endif else if not keyword_set(holdknob) then knob[change]*=inc
  return,newval
end
  

FUNCTION mcmc_new_likelihood,param,model,covar=covar,_EXTRA=ex
  COMMON KinMS_COMMON, fdata,obspars
  COMMON KinMS_COMMON_covar,invmult
  ;;;;;;;;;;;;;;;
  ;;; This function calls the model, then evaluates the likelihood.
  ;;;;;;;;;;;;;;;
  !EXCEPT = 0
  modout = call_FUNCTION(model,param,_EXTRA=ex,ret_vel=ret_vel)
     w=where(finite(fdata) and finite(modout),nvalid)
     IF nvalid ge 1 THEN BEGIN
        if keyword_set(covar) then begin
        ;;; Use inverse covariance matrix in common block to calculate chi2 
           s=size(modout)
           chiconv=0.0
           delta=reform(((fdata[w]-modout[w])))           
           for j=0,s[3]-1 do begin
              chiconv+=MATRIX_MULTIPLY(MATRIX_MULTIPLY(delta[(s[1]*s[2]*(j-1)):(s[1]*s[2]*j)-1],invmult, /ATRANSPOSE),delta[(s[1]*s[2]*(j-1)):(s[1]*s[2]*j)-1])
           endfor
        endif else chiconv = (total((((fdata[w]-modout[w])*(fdata[w]-modout[w]))/((obspars.rms)^2)),/double)) ;; compute chi2 directly assuming no covariance.
        like= -0.5*(chiconv-n_elements(w))
     ENDIF ELSE BEGIN
        like = -1e31
     ENDELSE
  return,(like)
END





function mcmc_new_chain,param,likelihood,model,iter,seed,outputll=outputll,knob=knob,plotchains=plotchains,burn=burn,holdknob=holdknob,acceptrate=acceptrate,covar=covar,log=log
  COMMON KinMS_COMMON, fdata,obspars

  ;; Run an MCMC chain with gibbs sampling and adaptive step refinement.
  
  max = double(param.max)
  min = double(param.min)
  change = where(param.changeable EQ 1,nchange)
  seed=systime(/sec)
  gibbs_change=floor(randomu(seed,iter)*n_elements(change))
  values=double(param.value)
  bestll=call_FUNCTION(likelihood,values,model,_EXTRA=ex,log=log,covar=covar)
  outputvals=dblarr(n_elements(values),iter)
  outputll=dblarr(iter)
  accepted=dblarr(iter)
  outputvals[*,0]=values
  outputll[0]=bestll



;;;; do plotting if required
  if keyword_set(plotchains) then begin
     windows=n_elements(change)
     size=ceil(sqrt(float(windows)))
     if size eq windows then size2=1 else size2=size
     psave=replicate({p:!p,x:!x,y:!y},windows)
     !p.multi=[0,size2,size]
     for i=0,windows-1 do begin
        if max[change[i]]-min[change[i]] gt 1000 then ylog=1 else ylog=0
        cgplot,[0],[param.value[change[i]]],yrange=[min[change[i]],max[change[i]]],xrange=[-1,iter+1],psym=1,/xstyle,/ystyle,ylog=ylog,xtitle=param.name[change[i]],charsize=1.5
        psave[i]={p:!p,x:!x,y:!y}
     endfor
  endif
;;;;;
  n_accept=0L
  for i=1,iter-1 do begin
     outputvals[*,i]=mcmc_new_take_step(outputvals[*,i-1],outputll[i-1],change[gibbs_change[i]],knob,min,max,seed,newll,accept,likelihood,model,holdknob=holdknob,tried=tried,covar=covar,log=log)
     outputll[i]=newll
     accepted[i]=accept
     if keyword_set(plotchains) then begin
        for k=0,n_elements(change)-1 do begin
           wind2look=gibbs_change[i]
           !p=psave[wind2look].p
           !x=psave[wind2look].x
           !y=psave[wind2look].y
           cgplot,[i],tried[change[gibbs_change[i]]],psym=1+accept,/overplot
        endfor
     endif
     n_accept+=accept
     if keyword_set(holdknob) then print, $
      strjoin(strarr(22)+string(byte(8)),''), $
      'KinMS_MCMC: ',100L*i/iter,' percent', $
      format= '($,A,A12,I2.2,A)'
  endfor
if keyword_set(holdknob) then print, strjoin(strarr(22)+string(byte(8)),'')+'KinMS_MCMC: done        ';;;

acceptrate=float(n_accept)/iter
if n_elements(burn) then begin
   outputvals= outputvals[*,burn:-1]
   outputll= outputll[burn:-1]
   accepted=accepted[burn:-1]
endif
outputvals=outputvals[*,where(accepted)]
outputll=outputll[where(accepted)]
return,outputvals
end

pro mcmc_new,param,iters,numatonce,burn,nchains,plotchains=plotchains,silent=silent,outputll=outputll,outputvalue=outputvalue,model=model,likelihood=likelihood,covar=covar,log=log

;;; An adaptive step MCMC utalizing gibbs sampling
  
  COMMON KinMS_COMMON, fdata,obspars

  if not keyword_set(model) then model="model"
  if not keyword_set(likelihood) then likelihood="mcmc_new_likelihood"
  targetrate= 0.25 ;; need the acceptance rate to be higher than this or dont accept
  
  t=systime(/sec)

  verybestparam=param
  verybestknob=1.0
  verybestll=-1e31

  for chainno=0,nchains-1 do begin
      if not keyword_set(silent) then print,"Doing chain",chainno
     count=0
     converged=0
     oldmean=param.value*0.0
     newmean=0.0
     outputvals=param.value
     knob=fltarr(n_elements(outputvals))+0.5*(param.max-param.min)
     while (count lt iters/numatonce) and (converged eq 0) do begin
        param.value=outputvals[*,-1]
        outputvals= mcmc_new_chain(param,likelihood,model,numatonce,seed,outputll=outputll,knob=knob,plotchains=plotchains,acceptrate=acceptrate,covar=covar,log=log)
        if acceptrate*numatonce gt 1 then newmean=mean(outputvals,dim=2) else newmean=oldmean
        if count eq 0 then converged=0 else begin
           test=abs(newmean-oldmean) lt param.precision
           if total(test) eq n_elements(test) and (acceptrate ge targetrate) then begin
              converged=1      
            if not keyword_set(silent) then print,"Chain converged: LL:",max(outputll)," - Accept rate:",acceptrate
           endif else begin
              if not keyword_set(silent) then begin
                 print,"     Chain has not converged - Accept rate:",acceptrate
                 if total(test) eq n_elements(test) then print,"Target rate not reached" else Print,"     Still varying:",param.name[where(test eq 0)]
              endif
           endelse
        endelse
        oldmean=newmean
        count++
     endwhile
      if not keyword_set(silent) then if converged eq 0 then print,"     Chain did not converge: saving last values"

     if (max(outputll) gt verybestll) or chainno eq 0 then begin
        if not keyword_set(silent) then print,"Best chain so far!"
        verybestparam=param
        verybestknob=knob
        verybestll=max(outputll)
     endif
     
     
  endfor
   
     
  if not keyword_set(silent) then print,"Starting final chain"
  param=verybestparam
  knob=verybestknob
  print,"verybestparam",verybestparam.value
  outputvalue= mcmc_new_chain(param,likelihood,model,iters,seed,outputll=outputll,knob=knob,plot=0,/holdknob,burn=burn,covar=covar,log=log)
  print,"verybestparam final",verybestparam.value
  sindx= reverse(sort(outputll))
  outputvalue= outputvalue[*,sindx]
  outputll= outputll[sindx]

  print,"time:", systime(/sec)-t
end


