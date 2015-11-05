
pro mk_chanmap,f,fsim,hdrs,rms,nplots,rmsfac=rmsfac,_extra=_extra,eps=eps,file=file,chans2do=chans2do,vsys=vsys
 if not keyword_set(eps) then begin
   device, decomposed=0
   !p.thick=1
   !x.thick=1
   !y.thick=1
   !p.charthick=1
   ;charthick=1
endif else begin
   !p.thick=5
   !x.thick=5
   !y.thick=5
   ;charthick=5
   !p.charthick=5
   set_plot, 'ps'
  if not keyword_set(file) then file="mk_chanmap"
     file+=".eps"
   device, filename=file, /color, /encap, xsize=16, ysize=10 ,/inches
endelse
ctload,3,/rev
 
  cube=f
  
   if not keyword_set(rmsfac) then rmsfac=3.
  ;cube[where(cube le rms)]=0.0
 

  s=size(cube)
  if not keyword_set(chans2do) then chans2do=[0,s[3]-1]
  
  x1=(findgen(s[1])-s[1]/2.)*abs(sxpar(hdrs,'cdelt1')*3600.)
  v1=((findgen(s[3])-sxpar(hdrs,'crpix3'))*(sxpar(hdrs,'cdelt3')/1e3)+sxpar(hdrs,'crval3')/1e3);-vsys+7.5
  
  plots2do=chans2do[1]-chans2do[0]+1
  gofromx=range(0.1,1,nplots[0]+1)
  gofromy=range(0.1,1,nplots[1]+1)

  xpos=gofromx[0:nplots[0]-1]
  ypos=reverse(gofromy[0:nplots[1]-1])
  
  dxp=gofromx[1]-gofromx[0]
  dyp=gofromy[1]-gofromy[0]
  
!p.multi=[0,nplots[0],nplots[1]]
  for i=chans2do[0],chans2do[1] do begin
     j=i-chans2do[0]
     mom0=f[*,*,i]
     mom0sim=fsim[*,*,i]
     ctload,0,/rev;,clip=[200,255]
     levels=range(rms,rms*rmsfac,5)
     xtickname=replicate(' ',10)
     ytickname=replicate(' ',10)
     
     if j mod nplots[0] eq 0 then ytickname="" ;else 
      if floor(j/nplots[0]) eq nplots[1]-1 then xtickname="" ;else 

      contour,mom0,x1,x1,levels=levels,background=0,color=255,/fill,_extra=_extra,charsize=2.5,/xstyle,/ystyle,noclip=0,xtickname=xtickname,ytickname=ytickname,xmargin=[0,0],ymargin=[0,0],position=[xpos[j mod nplots[0]],ypos[floor(j/nplots[0])],xpos[j mod nplots[0]]+dxp,ypos[floor(j/nplots[0])]+dyp],/nodata

      ctload,0,/rev,clip=[200,255]
      
      contour,mom0,x1,x1,levels=levels,background=0,color=255,/fill,/overplot;,noclip=0,CLIP=[!x.crange[0]+0.1,!y.crange[0]+0.1,!x.crange[0]-0.1,!y.crange[0]-0.1]
      ctload,0,/rev
      contour,mom0,x1,x1,levels=levels,background=0,color=255,/fill,_extra=_extra,charsize=2.5,/xstyle,/ystyle,noclip=0,xtickname=xtickname,ytickname=ytickname,xmargin=[0,0],ymargin=[0,0],position=[xpos[j mod nplots[0]],ypos[floor(j/nplots[0])],xpos[j mod nplots[0]]+dxp,ypos[floor(j/nplots[0])]+dyp],/nodata,/noerase
      
     clrlevels=range(10,250,chans2do[1]-chans2do[0]+1,/open)
     sauron_colormap

     if rms*rmsfac lt max(mom0,/nan) then begin
        levels=range(rms*rmsfac,max(mom0,/nan),3)

        ;clrlevels2=range(clrlevels[j],clrlevels[j]+20,3)
        contour,mom0,x1,x1,levels=levels[0],/fill,/overplot,_extra=_extra,color=clrlevels[j]
       ; contour,mom0,x1,x1,levels=levels[1],/fill,/overplot,_extra=_extra,color=clrlevels2[1]
       ; contour,mom0,x1,x1,levels=levels[2],/fill,/overplot,_extra=_extra,color=clrlevels2[2]
       ; stop
        ctload,0
        contour,mom0,x1,x1,levels=levels,color=100,/overplot,_extra=_extra
     endif

    ; levels=range(rms,max(mom0,/nan),5)
     ctload,39
     if keyword_set(eps) then contour,mom0sim,x1,x1,levels=levels,color=0,/over,thick=8. else contour,mom0sim,x1,x1,levels=levels,color=0,/over,thick=2.
     
   ;  xyouts,-12,8,string(v1[i],format='(I5)'),color=0,charsize=1.5
     endfor
  xyouts,0.5,0.05,'RA offset (")',color=0,charsize=1.5,/norm
 xyouts,0.05,0.5,'DEC offset (")',color=0,orient=90,charsize=1.5,/norm

  

 ; ellipse, 0.5*sxpar(hdrs,'bmaj')*3600., 0.5*sxpar(hdrs,'bmin')*3600., sxpar(hdrs,'bpa'), 0,360, !x.crange[0]+(sxpar(hdrs,'bmaj')*3600.*1.3), !y.crange[0]+(sxpar(hdrs,'bmaj')*3600.*1.3)
  
  !p.multi=[0,1,1]
   if keyword_set(eps) then begin
   device, /close
   spawn,"epstool --bbox --copy "+file+" la"+file
   spawn,'mv la'+file+" "+file
   spawn,"epstopdf "+file
   set_plot, 'x'
endif
 ;  stop
end
