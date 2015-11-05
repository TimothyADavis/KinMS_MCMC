
function KinMS_mcmc,param,_extra=_extra,finaloutput=finaloutput,iters=iters
  COMMON KinMS_COMMON, fdata,obspars
  !EXCEPT = 0
  ;;;;;;;;;;;;;;;
  ;;; This functions calls the MCMC and sets all required parameters
  ;;;;;;;;;;;;;;;  
  
  ;;;;;;;; changeable parameters - rest should be OK as they are
 if n_elements(iters) eq 0 then iters=3000l  ;;; make this bigger to do more samples (get better chi-sqr contours). Smaller for testing.
  plotit=1      ;;; Display plots. Good for testing, somewhat slower for real runs. 
  ;;;;;;;;
  print,"Iterations for final run: ",iters
  if total(param.value gt param.max)+total(param.value lt param.min) ne 0 then stop,"Initial Guesses out of range"
  if plotit then begin
     set_plot, 'x'
     device,decomposed=0
  endif
  ;print,iters
  ;;;;;Call MCMC
  
  mcmc,param,outputvalue=outputvalue,model='mkkinms_model',likelihood='mkkinms_likelihood',outputll=outputll, numatonce = 50l*total(param.changeable),random_start = 1l,numchains = 1l,numspinups = 300l,ranacc=-10,iter=iters,plot=plotit
 meep=""
 read,meep
 fsim=mkkinms_model(outputvalue[*,0],filename=finaloutput)
 finalparams=param
 finalparams.value=outputvalue[*,0]
 save,param,outputvalue,obspars,filename=finaloutput+"_parsave.sav"
 mkchiplots,finaloutput+"_parsave.sav"
 return,outputvalue[*,0]
END
