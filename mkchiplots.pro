

function like_est,bhmasses,likelihoods,thresh=thresh
  if not keyword_set(thresh) then thresh=0.68
  lev= findchicont(likelihoods,thresh)
  w=where(likelihoods ge lev[0])
  error=[min(bhmasses[w]),max(bhmasses[w])]
  return,error
end

function mkhisto,plotval,binvar=binvar,_extra=_extra
  if not keyword_set(binvar) then autobin=1 else autobin=0
  if keyword_set(binvar) then binsize=((max(plotval)-min(plotval))/binvar[0]) else binsize=0
  mom=moment(plotval)
  ctload,0
  plothist,plotval,xhist,yhist,bin=binsize,axiscolor=0,color=0,background=255,ytickname=replicate(' ',10),xmargin=[4,4],xticks=4,peak=1,yrange=[0,1.1],/ystyle,/xstyle,autobin=autobin,charsize=2,_extra=_extra,yticks=4
  sig1_errs=like_est(xhist,yhist,thresh=0.66)
  sig3_errs=like_est(xhist,yhist,thresh=0.99)
  bestfit=median(plotval)
  orderofmag=float(floor(alog10(abs(bestfit))))
  if orderofmag gt 1.0 then begin
     ending=") x10^"+strcompress(floor(orderofmag),/rem)
     print,"Best fit (1sig) = (",bestfit/10^orderofmag,"-",(bestfit-sig1_errs[0])/10^orderofmag,"+",(sig1_errs[1]-bestfit)/10^orderofmag,ending,format='(A20,F7.2,A1,F4.2,A1,F4.2,A7)'
     print,"Best fit (3sig) = (",bestfit/10^orderofmag,"-",(bestfit-sig3_errs[0])/10^orderofmag,"+",(sig3_errs[1]-bestfit)/10^orderofmag,ending,format='(A20,F7.2,A1,F4.2,A1,F4.2,A7)'
  endif else begin
     ending=")       "
     orderofmag=0.0
     print,"Best fit (1sig) = (",bestfit/10^orderofmag,"-",(bestfit-sig1_errs[0])/10^orderofmag,"+",(sig1_errs[1]-bestfit)/10^orderofmag,ending,format='(A20,F7.2,A1,F4.2,A1,F4.2,A7)'
     print,"Best fit (3sig) = (",bestfit/10^orderofmag,"-",(bestfit-sig3_errs[0])/10^orderofmag,"+",(sig3_errs[1]-bestfit)/10^orderofmag,ending,format='(A20,F7.2,A1,F4.2,A1,F4.2,A7)'
  endelse
  
  w=where(xhist ge sig1_errs[0] and xhist le sig1_errs[1])
  
  
  xfill=[xhist[w[0]]-binsize/2.d]
  yfill=[0d]
  for j=0,n_elements(w)-1 do begin
     xfill=[xfill,xhist[w[j]]-binsize/2.d,xhist[w[j]]+binsize/2.d]
     yfill=[yfill,yhist[w[j]],yhist[w[j]]]
  endfor
  xfill=[xfill,xhist[w[-1]]+binsize/2.d]
  yfill=[yfill,0d]
  ctload,0
  for j=0,n_elements(w)-1 do polyfill,xfill,yfill,color=200
  plot,xhist,yhist,psym=10,color=0,background=255,ytickname=replicate(' ',10),xmargin=[4,4],xticks=4,yrange=[0,1.1],/ystyle,/xstyle,charsize=2,_extra=_extra,/noerase,yticks=4
  ctload,39
end

function findchicont,hist,fraction_wanted 
  order = reverse(sort(hist))
  sum = total(hist)
  fract = dblarr(n_elements(hist))
  for i=0, n_elements(hist) - 1 do begin
     fract[i] = total(hist[order[0:i]])/sum
  endfor
  loc = where(abs(fract - fraction_wanted) eq min(abs(fract - fraction_wanted)))
  level = hist[order[loc]]
  return,level
end

pro mkchiplots,savfile,eps=eps,ranges=ranges,names=names,_extra=_extra,tickinterval=tickinterval,binvar=binvar,tickformat=tickformat
  if not keyword_set(binvar) then binvar=[15.,15]
  ctload,39
  
  restore,savfile
  params=param
  changes=where(params.changeable eq 1)
  
  schange=sort(params.name[changes])
  if n_elements(names) ne 0 then params.name[changes]=names
  changes=changes[schange]
  
  if keyword_set(eps) then begin
     set_plot, 'ps'
     file=savfile
     file+="_conts.eps"
     device, filename=file, /color, /encap, xsize=12, ysize=8 ,/inches
  endif                       
  !p.multi=[0,n_elements(changes),n_elements(changes)]
  posstep=0.85/(n_elements(changes))
  for i=0,n_elements(changes)-1 do begin
     plotval1=reform(outputvalue[changes[i],*])
     binsize1=((max(plotval1)-min(plotval1))/binvar[1])
     first=0
     for k=n_elements(changes)-1,0,-1 do begin
        if k lt i then plot,[0],[0],xstyle=4,ystyle=4 else begin
           position=[0.1+(i)*posstep,(0.90-(k)*posstep),0.1+(i+1)*posstep,0.90-(k-1)*posstep]
           plotval2=reform(outputvalue[changes[k],*])
           binsize2=((max(plotval2)-min(plotval2))/binvar[1])
           
           if n_elements(ranges) gt 0 then begin
              yrange=ranges[*,changes[k]]
              xrange=ranges[*,changes[i]]
           endif else begin
                   xrange=[min(plotval1)-binsize1*5.,max(plotval1)+binsize1*5.]
                   yrange=[min(plotval2)-binsize2*5.,max(plotval2)+binsize2*5.]
                endelse
                
                if n_elements(tickinterval) gt 0 then begin
                   ytickinterval=tickinterval[changes[k]]
                   xtickinterval=tickinterval[changes[i]]
                endif else begin
                   ytickinterval=0
                   xtickinterval=0
                endelse
                

                
                if k eq i then begin
                   if first eq 0 then begin
                      xtitle=textoidl(params.name[changes[i]])
                      undefine,xtickname
                   endif else begin
                      undefine,xtitle
                      
                   endelse
                   xtickname=replicate(' ',10)
                   jam=mkhisto(plotval1,position=position,binvar=binvar[0],xrange=xrange,xtickinterval=xtickinterval,xtickname=xtickname)

                   if first eq 0 then begin
                      numticks=4.
                      ypos = Replicate(!Y.Window[0] - 0.01, numticks+1)
                      xpos = position[0] + (position[2] - position[0]) * Findgen(numticks + 1) / numticks
                      labels = xrange[0] + (xrange[1] - xrange[0]) * Findgen(numticks + 1) / numticks
                      fudge=CONVERT_COORD([0,!D.Y_CH_SIZE],/device,/to_norm)/2.
                      fudge=reform(fudge[0])
                      FOR j=1, numticks-1 DO XYOutS, xpos[j], ypos[j]-fudge, strcompress(string(labels[j],format='(F10.4)'),/rem), Alignment=0.0, Orientation=-90, /Normal,color=0
                      XYOutS, position[0] + (position[2] - position[0])/2., ypos[0]-0.1, xtitle, Alignment=0.5, Orientation=0, /Normal,color=0
                   endif

                   
                endif else begin
                

          
             hist=hist_2d(plotval1,plotval2,bin1=binsize1,bin2=binsize2,min1=min(plotval1),max1=max(plotval1),min2=min(plotval2),max2=max(plotval2))
             s=size(hist)
             ctload,0,/rev

             if first eq 0 then begin
                xtitle=textoidl(params.name[changes[i]])
                undefine,xtickname
             endif else begin
                undefine,xtitle
                
             endelse
             xtickname=replicate(' ',10)
             
             if i eq 0 then begin
                ytitle=textoidl(params.name[changes[k]])
                undefine,ytickname
             endif else begin
                undefine,ytitle
                undefine,ytickformat 
             endelse
               ytickname=replicate(' ',10)
               
               
               lev2do=[0,([findchicont(hist,0.95),findchicont(hist,0.68)]>0.01)]


               
               cgcontour,hist,(findgen(s[1])*binsize1)+min(plotval1),(findgen(s[2])*binsize2)+min(plotval2),/fill,levels=lev2do,color=255,background=0,charsize=2,/closed,xmargin=xmargin,ymargin=ymargin,xtickname=xtickname,ytickname=ytickname,position=position,xrange=xrange,xtickinterval=xtickinterval,yrange=yrange,ytickinterval=ytickinterval,/xstyle,/ystyle,_extra=_extra,ytickformat=ytickformat,xticks=4,yticks=4 
               fudge=CONVERT_COORD([!D.Y_CH_SIZE,!D.Y_CH_SIZE],/device,/to_norm)/2.
               
               
               if first eq 0 then begin
                  numticks=4.
                  ypos = Replicate(!Y.Window[0] - 0.01, numticks+1)
                  xpos = position[0] + (position[2] - position[0]) * Findgen(numticks + 1) / numticks
                  labels = xrange[0] + (xrange[1] - xrange[0]) * Findgen(numticks + 1) / numticks
                  if n_elements(tickformat) gt 0 then xtickformat=tickformat[changes[i]]
                  FOR j=1, numticks-1 DO XYOutS, xpos[j]-fudge[0], ypos[j], strcompress(string(labels[j],format=xtickformat),/rem), Alignment=0.0, Orientation=-90, /Normal,color=255
                  XYOutS, position[0] + (position[2] - position[0])/2., ypos[0]-0.1, xtitle, Alignment=0.5, Orientation=0, /Normal,color=255
               endif
               if i eq 0 then begin
                  numticks=4.
                  xpos = Replicate(!x.Window[0] - 0.005, numticks+1)
                  ypos = position[1] + (position[3] - position[1]) * Findgen(numticks + 1) / numticks
                  labels = yrange[0] + (yrange[1] - yrange[0]) * Findgen(numticks + 1) / numticks
                  if n_elements(tickformat) gt 0 then ytickformatz=tickformat[changes[k]]
                  FOR j=1, numticks-1 DO XYOutS, xpos[j], ypos[j]-fudge[0], strcompress(string(labels[j],format=ytickformatz),/rem), Alignment=1.0, Orientation=0, /Normal,color=255
                  XYOutS, xpos[0]-0.065, position[1] + (position[3] - position[1])/2., ytitle, Alignment=0.5, Orientation=90, /Normal,color=255
               endif
               first++
            endelse
             endelse
     endfor
  endfor
  if keyword_set(eps) then begin
     device, /close
     spawn,"epstopdf "+file
     spawn,"epstool --bbox --copy "+file+" la"+file+" >/dev/null"
     spawn,'epstopdf la'+file
     spawn,'mv '+STR_REPLACE('la'+file,"eps","pdf")+" "+STR_REPLACE(file,"eps","pdf")
     spawn,'mv la'+file+" "+file
     set_plot, 'x'
  endif
  
  !p.multi=[0,1,1]
end
