function mcmc_new_stepper,values,change,knob,min,max,seed
  nchange=n_elements(change)
   newvals= values
   acceptable=0
   while acceptable eq 0 do begin   
       newvals[change]= values[change]+(randomn(seed,nchange,/double)*knob[change])
       if total(newvals[change] ge min[change]) eq nchange and total(newvals[change] le max[change]) eq nchange then acceptable++
    endwhile
   return,double(newvals)
end

function mcmc_new_take_step,values,bestll,change,knob,min,max,seed,newll,accept,likelihood,model,holdknob=holdknob,tried=tried
  dec=0.99
  inc=1.01
  log=1
  accept=1
  newval=mcmc_new_stepper(values,change,knob,min,max,seed)
  newll = (1.0)*call_FUNCTION(likelihood,newval,model,_extra=ex)
  if keyword_set(log) then begin
     ratio=  newll-bestll
     randacc=alog(randomu(seed2,/double))
  endif else begin
     ratio= newll/bestll
     randacc=randomu(seed2,/double)
  endelse
;  stop
  tried=newval
  if (newll le bestll) AND (randacc gt ratio) then begin
     ;; dont accept
     accept=0
     newval = values
     newll = bestll
     if not keyword_set(holdknob) then knob[change]*=dec 
  endif else if not keyword_set(holdknob) then knob[change]*=inc          ;< 1.
  ;stop
  return,newval
end
  


FUNCTION mcmc_new_likelihood,param,model,_EXTRA=ex
  COMMON KinMS_COMMON, fdata,obspars
  ;;;;;;;;;;;;;;;
  ;;; This function calls the model, then evaulates the likelihood.
  ;;;;;;;;;;;;;;;
  !EXCEPT = 0
  modout = call_FUNCTION(model,param,_EXTRA=ex)
  w=where(finite(fdata) and finite(modout),nvalid)
  IF nvalid ge 1 THEN BEGIN
     chisq = total((((fdata[w]-modout[w])*(fdata[w]-modout[w]))/((obspars.rms)^2)),/double)
     like= -0.5*(chisq-n_elements(fdata[w]))
  ENDIF ELSE BEGIN
     like = 1e31
  ENDELSE 
  IF like EQ 1e31 THEN stop
  return,(like)
END


function mcmc_new_chain,param,likelihood,model,iter,outputll=outputll,knob=knob,plotchains=plotchains,burn=burn,holdknob=holdknob,acceptrate=acceptrate
  COMMON KinMS_COMMON, fdata,obspars
  max = double(param.max)
  min = double(param.min)
  change = where(param.changeable EQ 1,nchange)
  seed=systime(/sec)
  values=double(param.value)
  bestll=call_FUNCTION(likelihood,values,model,_EXTRA=ex)
  outputvals=fltarr(nchange,iter)
  outputll=fltarr(iter)
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
        cgplot,[0],[param.value[change[i]]],yrange=[min[change[i]],max[change[i]]],xrange=[-1,iter+1],psym=1,/xstyle,/ystyle,ylog=ylog,xtitle=param.name[change[i]],charsize=1.5;,/xlog
        psave[i]={p:!p,x:!x,y:!y}
     endfor
  endif
;;;;;
  n_accept=0L
  for i=1,iter-1 do begin
     outputvals[*,i]=mcmc_new_take_step(outputvals[*,i-1],outputll[i-1],change,knob,min,max,seed,newll,accept,likelihood,model,holdknob=holdknob,tried=tried)

                             
     outputll[i]=newll
     ;; print,knob
     if keyword_set(plotchains) then begin
        for k=0,n_elements(change)-1 do begin
           wind2look=change[k]
           !p=psave[wind2look].p
           !x=psave[wind2look].x
           !y=psave[wind2look].y
           cgplot,[i],tried[change[k]],psym=1+accept,/overplot
         ;  wait,0.1
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
  endif
  return,outputvals
end

pro mcmc_new,param,iters,numatonce,burn,nchains,plotchains=plotchains,silent=silent,outputll=outputll,outputvalue=outputvalue,model=model,likelihood=likelihood
  COMMON KinMS_COMMON, fdata,obspars

  if not keyword_set(model) then model="model"
  if not keyword_set(likelihood) then likelihood="mcmc_new_likelihood"
  
  
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
        outputvals= mcmc_new_chain(param,likelihood,model,numatonce,outputll=outputll,knob=knob,plotchains=plotchains,acceptrate=acceptrate)
        newmean=mean(outputvals,dim=2)   
        if count eq 0 then converged=0 else begin
           test=abs(newmean-oldmean) lt param.precision
           if total(test) eq n_elements(test) then begin
              converged=1      
            if not keyword_set(silent) then print,"Chain converged: LL:",max(outputll)," - Accept rate:",acceptrate
           endif else begin
              if not keyword_set(silent) then begin
                 print,"     Chain has not converged - Accept rate:",acceptrate
                 Print,"     Still varying:",param.name[where(test eq 0)]
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
  outputvalue= mcmc_new_chain(param,likelihood,model,iters,outputll=outputll,knob=knob,plot=0,/holdknob,burn=burn)
      
  sindx= reverse(sort(outputll))
  outputvalue= outputvalue[*,sindx]
  outputll= outputll[sindx]

  print,"time:", systime(/sec)-t
end


