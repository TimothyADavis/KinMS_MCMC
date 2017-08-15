
function KinMS_mcmc,param,_extra=_extra,iters=iters,outputll=outputll,finaloutput=finaloutput,silent=silent
  COMMON KinMS_COMMON, fdata,obspars
  !EXCEPT = 0
  ;;;;;;;;;;;;;;;
  ;;; This functions calls the MCMC and sets all required parameters
  ;;;;;;;;;;;;;;;  
  
  ;;;;;;;; changeable parameters - rest should be OK as they are
  if n_elements(iters) eq 0 then iters=3000l  ;;; make this bigger to do more samples (get better chi-sqr contours). Smaller for testing.
  plotit=1                                    ;;; Display plots. Good for testing, somewhat slower for real runs. 
  ;;;;;;;;
  print,"Iterations for final run: ",iters

  ;;;;;;;; Return an error if the initial guesses are outside of the allowed range
  testmax=param.value gt param.max
  testmin=param.value lt param.min
  if total(testmax)+total(testmin) ne 0 then begin
     print,"Initial Guesses out of range: ",param.name[where(testmax eq 1 or testmin eq 1)]
     stop,"Halting"    
  endif
  ;;;;;;;;;

  
  if plotit then begin
     set_plot, 'x'
     device,decomposed=0
  endif
  
  ;;;;;Call MCMC

  mcmc_new,param,iters,25l*total(param.changeable),300,1l,model="mkkinms_model",likelihood="mcmc_new_likelihood",plotchains=plotit,silent=silent,outputll=outputll,outputvalue=outputvalue

  save,param,outputvalue,obspars,outputll,filename=finaloutput+"_parsave.sav"
  ;;;;;;
  meep=""
  read,meep
  fsim=mkkinms_model(outputvalue[*,0],filename=finaloutput)
  mkchiplots,finaloutput+"_parsave.sav"
  return,outputvalue[*,0]
END
