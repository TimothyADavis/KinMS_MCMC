pro lookatplot_mkplot,fdata,fsim,params,obspars,_extra=_extra,eps=eps
  
   r=params
   s=size(fsim)
   fdataclip=fdata
   fdataclip[where(fdata lt 3*obspars.rms)]=0.

   
   ;;; Create rough moment zero and one maps and a major axis PVD from
   ;;; the simulation and data
   mom0sim=total(fsim,3)
   mom0data=total(fdataclip,3)
   v1=(findgen(s[3])-(s[3]/2.))*obspars.dv
   x1=(findgen(s[1])-(s[1]/2.))*obspars.cellsize
   mom1sim=mom0sim*0.0
   mom1data=mom0data*0.0
   for i=0,s[1]-1 do begin
      for j=0,s[2]-1 do begin
         mom1sim[i,j]=total(v1*fsim[i,j,*])/total(fsim[i,j,*])
         mom1data[i,j]=total(v1*fdataclip[i,j,*])/total(fdataclip[i,j,*])
      endfor
   endfor
   mom1sim[where(mom0sim lt 0.1*max(mom0sim))] =-10000
   mom1data[where(mom0data lt 0.1*max(mom0data))] =-10000
   pvdcubedata=fdataclip*0.0
   pvdcubesim=fsim*0.0
   for i=0,s[3]-1 do pvdcubedata[*,*,i]=rot(fdataclip[*,*,i],r[1]-90,/interp)
   for i=0,s[3]-1 do pvdcubesim[*,*,i]=rot(fsim[*,*,i],r[1]-90,/interp)
   pvddata=reform(total(pvdcubedata[*,(s[2]/2.)-2.:(s[2]/2.)+2,*],2))
   pvdsim=reform(total(pvdcubesim[*,(s[2]/2.)-2.:(s[2]/2.)+2,*],2))
;;;
   
   ;;; Create the plot to your x-device
   device,decomposed=0,RETAIN=2
   loadct,39,/silent
   !p.multi=[0,2,2]
   contour,mom0data,x1,x1,levels=max(mom0data)*((findgen(9)+1)/10.),xrange=[-25,25],yrange=[-25,25],xtitle="RA (arcsec)",ytitle="DEC (arcsec)",/xstyle,/ystyle,background=255,color=0
   contour,mom0sim,x1,x1,/overplot,color=250,nlevels=10,levels=max(mom0data,/nan)*((findgen(9)+1)/10.)
   al_legend,["Moment Zero"],/top,/right,box=0,textcolor=0
   al_legend,["Data","Model"],/bottom,/left,box=0,textcolor=0,color=[250,0],linestyle=[0,0] 
   contour,pvddata,x1,v1,levels=max(pvddata)*((findgen(9)+1)/10.),xrange=[-30,30],yrange=[-200,200],/xstyle,/ystyle,xtitle="Position (arcsec)",ytitle="Velocity (km/s)",background=255,color=0
   contour,pvdsim,x1,v1,levels=max(pvddata)*((findgen(9)+1)/10.),/overplot,color=250
   al_legend,["PVD"],/top,/left,box=0,textcolor=0
   contour,mom1data,x1,x1,/fill,levels=v1,xrange=[-25,25],yrange=[-25,25],/cell,color=0,xtitle="RA (arcsec)",ytitle="DEC (arcsec)",xmargin=[15,0],/xstyle,/ystyle
   al_legend,["Moment One (data)"],/top,/right,box=0,textcolor=0
   contour,mom1sim,x1,x1,/fill,levels=v1,xrange=[-25,25],yrange=[-25,25],/cell,color=0,xtitle="RA (arcsec)",ytickname=replicate(' ',10),xmargin=[0,15],/xstyle,/ystyle
   al_legend,["Moment One (model)"],/top,/right,box=0,textcolor=0
   !p.multi=[0,1,1]
;;;
end
