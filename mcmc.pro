PRO mcmc,param,$
         numatonce=numatonce,random_start=random_start,numchains=numchains,numspinups=numspinups,$
         iter=iter,model=model,likelihood=likelihood,$
         outputll=outputll,outputvalue=outputvalue,outputy=outputy,_EXTRA = ex,fast=fast,medium=medium,quiet=quiet,superfast=superfast,plotchains=plotchains
COMMON KinMS_COMMON, fdata,obspars

  
;MCMC parameter estimator
;based on sipnet
;Note: this can work for multiple data types (likelihood function
;would do the math)
;SEE BELOW FOR AN EXAMPLE FOR HOW TO USE IN IDL

;---REQUIRED INPUTS---
;x is the input data for the model (no requirements on format except
;for what likelihood and model expect)
;y is the output data to compare (same as above)

;param is a structure  with the following properties
;param.name is an arry of name
;param.value is parameter value (initial guess) array
;param.max max value array
;param.min min value array
;param.knob is knob array (not used in this version, can be all zero)
;param.changeable is whether it should be fixed (0) or estimated (1)

;---OPTIONAL KEYWORDS---

;model is the name of the model function (string) - default is "mcmc_testmodel"
;likelihood is the name of the likelihood function - it should call
;   model and compute the likelihood - default is "mcmc_likelihood"
;_EXTRA can be used to pass extra keywords to likelihood or model

;numatonce - how often to check for chain convergence (default is 10000)
;random_start - start at init guess if 0, else randomly within prior (default is 0)
;numchains - number of chains (default is 10)
;numspinups - how much burn in on final iteration (default is 125000)
;iter - max number of iterations both for chains and final (default is 375000)
;/fast, /medium, and /superfast are keywords with different default settings for the above settings - all are faster than default (see below)
;/quiet turns off all printing of messages

;---OUTPUTS---
;outputll is the likelihood of output values, best value is first
;outputvalues are accepted parameter values (same format as param.value)
;outputy is the model run with the best outputvalues parameter set (same format as y)

  A_STAR = 0.1
  THRESH = 0.02 
  DEC = 0.95
  INC = 1/DEC
  add_fraction = 1.0


;set defaults
  IF n_elements(numatonce) EQ 0 THEN numatonce = 10000l
  IF n_elements(random_start) EQ 0 THEN random_start = 0l
  IF n_elements(numchains) EQ 0 THEN numchains = 10l
  IF n_elements(numspinups) EQ 0 THEN numspinups = 125000l
  IF n_elements(iter) EQ 0 THEN iter = 375000l
  IF n_elements(model) EQ 0 THEN model = 'mcmc_testmodel'
  IF n_elements(likelihood) EQ 0 THEN likelihood = 'mcmc_likelihood'

  max = double(param.max)
  min = double(param.min)
  range = max-min
  change = where(param.changeable EQ 1,nchange)
   
  verybestll = double(-1e31)

;build some chains
  FOR c = 0,numchains-1 DO BEGIN

    IF ~keyword_set(quiet) THEN print,'Starting Chain ',c+1

    ;reset values
    value = double(param.value)
    IF (random_start EQ 1) AND (c GT 0) THEN value[change] = min[change] + (range[change] * randomu(systime(/sec),nchange))
    knob = double(param.knob)
    converged = 0
    steps = 0l
    seed = systime(/sec)
    oldvalue = value
    bestvalue = value
    ll = (-1.0) * call_FUNCTION(likelihood,value,model,_extra=ex)
    bestll = ll
    ll_old = bestll

;go through the chain until convergence

    means=fltarr(ceil(iter/numatonce) > 10,n_elements(value))
    stdevs=means
    whilecounter=0

    WHILE (converged EQ 0) && (steps LT iter) DO BEGIN 
      ichgs = change[long(randomu(seed,numatonce,/double)*nchange)]
      tune = randomu(seed,numatonce,/double)-0.5
      ran_accept = alog(randomu(seed,numatonce,/double)) ;; uniform between 0 and 1, but ln'ed
      accepted_vals=fltarr(numatonce,n_elements(value))

      if keyword_set(plotchains) then begin
         windows=n_elements(change)
         size=ceil(sqrt(float(windows)))
         if size eq windows then size2=1 else size2=size
         psave=replicate({p:!p,x:!x,y:!y},windows)
         !p.multi=[0,size2,size]
         for i=0,windows-1 do begin
            if max[change[i]]-min[change[i]] gt 1000 then ylog=1 else ylog=0
            cgplot,[0],[value[change[i]]],yrange=[min[change[i]],max[change[i]]],xrange=[-1,numatonce+1],psym=1,/xstyle,/ystyle,ylog=ylog,xtitle=param.name[change[i]]
            psave[i]={p:!p,x:!x,y:!y}
         endfor
      endif

      yes = 0l
      
      FOR k = 0l,numatonce-1l DO BEGIN
;randomly pick a parameter to change
        accept = 1
        ichg = ichgs[k]
        oldval = value[ichg] 
        newval = (knob[ichg] * range[ichg] * tune[k])+oldval
        IF (newval GT max[ichg]) OR (newval LT min[ichg]) THEN accept = 0
        
;run the model and calculate the likelihood
        IF accept EQ 1 THEN BEGIN
          value[ichg] = newval
          
          if keyword_set(plotchains) then begin
             wind2look=where(ichg eq change)
             !p=psave[wind2look].p
             !x=psave[wind2look].x
             !y=psave[wind2look].y
             cgplot,[k],[newval],psym=1,/overplot
          endif

          ll = (-1.0) * call_FUNCTION(likelihood,value,model,_extra=ex)
          IF (ll LE ll_old) && (ran_accept[k] ge ll-ll_old) THEN accept = 0
          IF (accept EQ 1) && (ll GT bestll) THEN BEGIN 
            bestll = ll
            bestvalue = value
          ENDIF
        ENDIF

       
;keep track of accepted parameter sets, tune knob
;       print,"ll ratio",ll-ll_old
        IF accept EQ 1 THEN BEGIN 
          ll_old = ll
          yes++
          knob[ichg]=(knob[ichg]*(1+((INC-1)*param.knobchange[ichg])))<(1.0)   
          if keyword_set(plotchains) then begin
             wind2look=where(ichg eq change)
                !p=psave[wind2look].p
                !x=psave[wind2look].x
                !y=psave[wind2look].y
                cgplot,[k],[newval],psym=2,/overplot;,/NOERASE
          endif
          accepted_vals[k,ichg]=value[ichg]
         ; stop
        ENDIF ELSE BEGIN
          value[ichg] = oldval
          knob[ichg] = (knob[ichg]*(1-((1-DEC)*param.knobchange[ichg])))
       ENDELSE
 ENDFOR

;check for convergence of this chain
      steps+=numatonce
     
      size=ceil(sqrt(float(n_elements(change))))
         !p.multi=[0,size,size]
         
     if yes gt 1 then begin
      for iternum=0,n_elements(change)-1 do begin
         w=where(accepted_vals[*,change[iternum]] ne 0.0)
         moms=moment(accepted_vals[w,change[iternum]],/nan)
         means[whilecounter,change[iternum]]=moms[0]
         stdevs[whilecounter,change[iternum]]=moms[1]
         binsize=((max(accepted_vals[w,change[iternum]])-min(accepted_vals[w,change[iternum]]))/float(n_elements(w)))*5.
           if binsize ne 0.0 then if keyword_set(plotchains) then plothist,accepted_vals[w,change[iternum]],bin=binsize else print,"converged to precision"
        endfor
      !p.multi=[0,1,1]
   endif
     IF ~keyword_set(quiet) THEN print,'  iteration ',steps,' accept ',(float(yes)/numatonce*100.),'% llmax ',bestll
     
     if whilecounter eq 0 then begin
        yes=0l
     endif else begin
                                ; stop
        
        
        errorvalue=abs(2*sqrt((stdevs[whilecounter,change]^2)+(stdevs[whilecounter-1,change]^2)))
        errorvalue = abs(((errorvalue/means[whilecounter,change]) > 0.01)*means[whilecounter,change])
        
        IF min(abs(means[whilecounter,change]-means[whilecounter-1,change]) le errorvalue) AND yes gt 3 then begin
           print,"Chain converged" 
        endif else begin
           print,"Chain has not converged" ; [enter to continue]"
         endelse
         IF min(abs(means[whilecounter,change]-means[whilecounter-1,change]) le errorvalue) AND yes gt 3 then begin
            ;; let the chain converge only if ALL the parameters that
            ;; are free agree within their own errors (or 1%, if larger), and at least 3
            ;; points are accepted.
            converged = 1
            IF ~keyword_set(quiet) THEN print,'  Chain ',c+1,' Converged LL: ',ll
            IF ~keyword_set(quiet) THEN print,'  Values: ',value
            IF bestll GE verybestll THEN BEGIN
               IF ~keyword_set(quiet) THEN print,'    And it is the best chain so far!'
               verybestll = bestll
               verybestvalue = bestvalue
               endvalue = value
               endll = ll
               endknob = knob
            ENDIF
         ENDIF ELSE BEGIN
            IF ~keyword_set(quiet) THEN print,'  Chain ',c+1,' not yet converged'
            yes = 0l
         ENDELSE 
      endelse
      whilecounter++

   ENDWHILE 
    
    IF converged EQ 0 THEN IF ~keyword_set(quiet) THEN print,'  Chain ',c+1,' did not converge'
    
  ENDFOR

;start at end of best chain
time=systime(/sec)

  
  IF n_elements(endvalue) EQ 0 THEN BEGIN
    IF ~keyword_set(quiet) THEN print,'No chains converged, try changing mcmc iterations'
    IF ~keyword_set(quiet) THEN print,'Starting from best value'
    endvalue = bestvalue
    endll = bestll
    endknob = knob
;    stop
  ENDIF 

  value = endvalue
  ll_old = endll
  knob = endknob
  seed = systime(/sec)
  bestvalue = value
  bestll = endll
  ichgs = change[long(randomu(seed,iter,/double)*nchange)]
  tune = randomu(seed,iter,/double)-0.5
  ran_accept = alog(randomu(seed,iter)) ;; uniform between 0 and 1, but ln'ed
  yes = 0l
  yes2 = yes

  outputll = fltarr(iter)
  outputvalue = fltarr(n_elements(value),iter)  

  FOR k = 0l,iter-1l DO BEGIN 
    IF k MOD numatonce EQ 0 THEN IF ~keyword_set(quiet) THEN print,'Final iteration ',k,' accepted ',yes2,' saved ',yes
;randomly pick a parameter to change
    accept = 1
    ichg = ichgs[k]
    oldval = value[ichg]
    newval = (knob[ichg] * range[ichg] * tune[k])+oldval
    IF (newval GT max[ichg]) OR (newval LT min[ichg]) THEN accept = 0

;run the model and calculate the likelihood
    IF accept EQ 1 THEN BEGIN
      value[ichg] = newval
      ll = (-1.0) * call_FUNCTION(likelihood,value,model,_extra=ex)

      IF (ll LE ll_old) && (ran_accept[k] GE ll-ll_old) THEN accept = 0
      IF (accept EQ 1) && (ll GT bestll) THEN BEGIN 
        bestll = ll
        bestvalue = value
      ENDIF
    ENDIF

;output values
    IF accept EQ 1 THEN BEGIN
      yes2++
      ll_old = ll
;if past numspinups, then start saving vals
      IF k GE numspinups THEN BEGIN 
        outputll[yes] = ll
        outputvalue[*,yes] = value
        yes++
      ENDIF 
    ENDIF ELSE BEGIN
      value[ichg] = oldval
    ENDELSE 

  ENDFOR

;create the history of accepted values
  IF yes EQ 0 THEN BEGIN 
    outputvalue = bestvalue
    outputll = bestll
  ENDIF ELSE BEGIN 
    outputvalue = [[bestvalue],[outputvalue[*,0:(yes-1)]]]
    outputll = [bestll,outputll[0:(yes-1)]]
    srt = reverse(sort(outputll))
    outputvalue = outputvalue[*,srt]
    outputll = outputll[srt]
  ENDELSE 

;output values if outputy is there
  IF arg_present(outputy) THEN BEGIN
    dummy = call_FUNCTION(likelihood,outputvalue[*,0],model,_extra=ex,modout=outputy)
  ENDIF 

  IF ~keyword_set(quiet) THEN print,'MCMC complete '
  IF ~keyword_set(quiet) THEN print,'Best LL: ',outputll[0]
  IF ~keyword_set(quiet) THEN print,'Values: ',outputvalue[*,0]
  print,"time taken",systime(/sec)-time
END
