;---------------------------------------------;
; DEPENDENT TASKS                             ;
;---------------------------------------------;

function prank,array,p
  nn = n_elements(array)
  ip =  long(float(p)*nn/100.)
  ii = sort(array)
  return,array(ii(ip))
end

pro vline, val,_extra=extra, min=min, max=max
  if !y.type eq 1 then yrange = 10^!y.crange else yrange = !y.crange
  nv = n_elements(val)
  for i = 0, nv-1 do oplot,fltarr(2)+val[i],yrange,_extra=extra
end

;---------------------------------------------;
; MAIN ROUTINE                                ;
;---------------------------------------------;

pro mospec1d,file,small=small,large=large,stddev=stddev,bootflag=bootflag,part=part,nomodel=nomodel,cpoly=cpoly,fire=fire,twocomp=twocomp,ascii=ascii

if n_params() lt 1 then begin
   print,'Syntax - MOSPEC1D, file, [/small, /large, /formal, /bootflag, /part, /nomodel, cpoly=]'
   print,'         Default is to use the full emission line catalog, cpoly=1, and the standard'
   print,'            error on the mean as the error spectrum for stacks'
   return
endif

resolve_all,/quiet,/continue_on_error,skip_routines='rsex'

if keyword_set(small) then small=1 else small=0
if keyword_set(large) then large=1 else large=0
if keyword_set(stddev) then formal=0 else formal=1
if keyword_set(bootflag) then bootflag=1 else bootflag=0
if keyword_set(part) then part=1 else part=0
if keyword_set(nomodel) then nomodel=1 else nomodel=0
if keyword_set(fire) then fire=1 else fire=0
if n_elements(cpoly) eq 0 then cpoly=1
if keyword_set(twocomp) then twocomp=1 else twocomp=0
if keyword_set(ascii) then ascii=1 else ascii=0

if strpos(file,'lris') ge 0 then lrisflag = 1 else lrisflag = 0

calibdir = getenv('MOSPEC_CALIB')
if strmid(calibdir,strlen(calibdir)-1) ne '/' then calibdir += '/'
if file_search('lines.dat') ne '' then readcol,'lines.dat',lines,airwav,vacwav,format='A,D,D',/silent 
if file_search('lines.dat') eq '' then begin
   if not part then readcol,calibdir+'rest_optical_emlines_vac.dat',lines,airwav,vacwav,format='A,D,D',/silent
   if part then readcol,calibdir+'rest_optical_lines_vac_part.dat',lines,airwav,vacwav,format='A,D,D',/silent
endif
linewav = vacwav

print,'Reading in '+file
if ~ascii then begin
   im1 = readfits(file,h1,exten=0,/silent)
   if fire then im1 /= 10.
   dim = size(im1,/dimen)
;;if dim[1] lt 2 or dim[1] gt 3 then print,'1D spectrum not in correct format'
   if file_test((strsplit(file,'.',/extract))[0]+'.stddev.fits') then stackflag = 1 else stackflag = 0
   if stackflag then begin
      stddev = readfits((strsplit(file,'.',/extract))[0]+'.stddev.fits',h2,exten=0,/silent)
      nspec = readfits((strsplit(file,'.',/extract))[0]+'.nspec.fits',h2,exten=0,/silent) & nspec = nspec[*,0]
      frej = readfits((strsplit(file,'.',/extract))[0]+'.frej.fits',h2,exten=0,/silent) & frej = frej[*,0]   
   endif
   
   object = strtrim(sxpar(h1,'OBJECT',/silent))
   field = (strsplit(object,'-',/extract))[0]
   unit = sxpar(h1,'BUNIT')
   filter = strtrim(sxpar(h1,'FILTER'))

   lambda = dblarr(dim[0])
   refpix = sxpar(h1,'CRPIX1')
   lam0 = sxpar(h1,'CRVAL1')
   delta = sxpar(h1,'CD1_1')
   if not keyword_set(delta) and fire then delta = sxpar(h1,'CDELT1')
   for n=long(0),long(n_elements(lambda))-1 do lambda[n]=lam0+(n-(refpix-1))*delta
   if fire then lambda = 10^lambda
;; lambda /= (1+sxpar(h1,'REDSHIFT'))
;; im1 *= (1+sxpar(h1,'REDSHIFT'))
   spec = im1[*,0]
   if n_elements(dim) eq 2 then errspec = im1[*,1] else errspec = spec*0
   if fire then begin
      errspec = readfits((strsplit(file,'_',/extract))[0]+'_E.fits',tmp,exten=0,/silent)
      ;; errspec = 1./sqrt(errspec)
      errspec /= 10.
      ;; badindex = where(~finite(spec) or ~finite(errspec))
      ;; spec[badindex] = 0
      ;; errspec[badindex] = 0
   endif
   if stackflag and formal then errspec = im1[*,1]
   if stackflag and not formal then errspec = stddev[*,0]
endif else begin
   stackflag = 0
   field = ''
   object = ''
   readcol,file,lambda,spec,errspec,format='D,D,D'
   nspec = replicate(1,n_elements(lambda))
endelse

specmask = 'specmask.dat'
if file_test(calibdir+specmask) then begin
   readcol,calibdir+specmask,masklam1,masklam2,format='F,F',/silent
   tempindex=fltarr(n_elements(lambda))
   for j=0, n_elements(masklam1)-1 do begin
      index=where(lambda ge masklam1[j] and lambda le masklam2[j])
      tempindex[index]=1.0
   endfor
   tempindex[where(lambda lt 1220 or lambda gt 1700)] = 1.0
   maskindex = where(tempindex eq 1.)
endif
if not file_test(calibdir+specmask) then begin
   tempindex=fltarr(n_elements(lambda))
   tempindex[where(lambda lt 1220 or lambda gt 1700)] = 1.0
   maskindex = where(tempindex eq 1.)
endif
maskspec = spec
maskspec[maskindex] = 0.0/0
goodindex = where(finite(maskspec))

lineflag=fltarr(n_elements(lambda))
for n=0,n_elements(linewav)-1 do begin
   index=where(abs(lambda-linewav[n]) lt 10.0)
   lineflag[index]=1.0
endfor

modelflag = 0
if not stackflag then begin
   seddir = getenv("MOSPEC_SED")
   sedfile = seddir+'bc03calz/'+strlowcase(field)+'/bestfit/bestfit.'+object+'.csf_agegt50.dat'
   if file_test(sedfile) and not nomodel then begin
      modelflag = 1
      print,'Reading in bestfit.'+object+'.csf_agegt50.dat'
      readcol,sedfile,modellambda,modelspec,format='D,D',/silent
      modelspec /= (1.0e-17)
      readcol,seddir+'inphot/'+strlowcase(field)+'_inphot.dat',sed_object,sed_redshift,format='A,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,D',/silent
      index = (where(sed_object eq object))[0]
      if index eq -1 then begin
         print,'No redshift found in SED inphot file'
         modelflag = 0
      endif else begin
         modellambda /= (1.+sed_redshift[index])
         modelspec_spl = spl_init(modellambda,modelspec,/double)
         modelspec_int = spl_interp(modellambda,modelspec,modelspec_spl,lambda,/double)
      endelse
   endif
endif

if lrisflag then begin
   models = file_search(calibdir+'bpass.v2p1.spectra/bpass.v2p1.300bin.*.t8.0.fnu.fits')
   chisq = fltarr(n_elements(models))
   best_ebmv = fltarr(n_elements(models))
   best_scale = fltarr(n_elements(models))
   for m = 0,n_elements(models)-1 do begin
      modelspec = readfits(models[m],mhead,exten=0,/silent)
      modellambda = dblarr(n_elements(modelspec))
      modelrefpix = sxpar(mhead,'CRPIX1')
      modellam0 = sxpar(mhead,'CRVAL1')
      modeldelta = sxpar(mhead,'CD1_1')
      for n=long(0),long(n_elements(modellambda))-1 do modellambda[n]=modellam0+(n-(modelrefpix-1))*modeldelta
      ;; modelspec_spl = spl_init(modellambda,modelspec,/double)
      ;; modelspec_int = spl_interp(modellambda,modelspec,modelspec_spl,lambda,/double)

      ebmv = findgen(51)/100+0.1
      scale = findgen(21)/100+0.9
      chi2 = dblarr(n_elements(ebmv),n_elements(scale))
      for e=0,n_elements(ebmv)-1 do begin
         for s=0,n_elements(scale)-1 do begin
            redmodelspec = modelspec
            for n=0,n_elements(redmodelspec)-1 do redmodelspec[n] *= 10^(-kmosdef(modellambda[n])*ebmv[e]/2.5)
            redmodel_spl = spl_init(modellambda,redmodelspec,/double)
            redmodel_int = spl_interp(modellambda,redmodelspec,redmodel_spl,lambda,/double)
            matchscale = median(spec[goodindex],/even)/median(redmodel_int[goodindex],/even)
            redmodel_int *= matchscale*scale[s]
            errspec = intarr(n_elements(errspec))+1
            chi2[e,s] = total((redmodel_int[goodindex]-spec[goodindex])^2./errspec[goodindex])
            ;; plot,lambda,spec,psym=10
            ;; oplot,lambda,redmodel_int,color=fsc_color('red')
            ;; wait,0.1
         endfor
      endfor
      array_indices = array_indices(chi2,where(chi2 eq min(chi2)))
      print,models[m],' ',min(chi2),ebmv[array_indices[0]] ;,minmax(ebmv[where(chi2 lt 20)])
      best_ebmv[m] = ebmv[array_indices[0]]
      best_scale[m] = scale[array_indices[1]]
      chisq[m] = min(chi2)
   endfor

   m = (where(chisq eq min(chisq)))[0]
   
   modelspec = readfits(models[m],mhead,exten=0,/silent)
   modellambda = dblarr(n_elements(modelspec))
   modelrefpix = sxpar(mhead,'CRPIX1')
   modellam0 = sxpar(mhead,'CRVAL1')
   modeldelta = sxpar(mhead,'CD1_1')
   for n=long(0),long(n_elements(modellambda))-1 do modellambda[n]=modellam0+(n-(modelrefpix-1))*modeldelta
   
   redmodelspec = modelspec
   for n=0,n_elements(redmodelspec)-1 do redmodelspec[n] *= 10^(-kmosdef(modellambda[n])*best_ebmv[m]/2.5)
   redmodel_spl = spl_init(modellambda,redmodelspec,/double)
   redmodel_int = spl_interp(modellambda,redmodelspec,redmodel_spl,lambda,/double)
   matchscale = median(spec[goodindex],/even)/median(redmodel_int[goodindex],/even)
   redmodel_int *= matchscale*best_scale[m]
   modelflag = 1
   
endif
   
xrange = [min(lambda),max(lambda)]
yrange = [0,1.1*max(spec[where(finite(spec))])] & yrange[0] = -1*yrange[1]/3.
yrange = prank(spec,[1,99])

dimensions = get_screen_size(resolution=resolution)
zoom = (dimensions[0]-226.)/float(n_elements(lambda))
if small then zoom *= 1/2. else if not large then zoom *= 2/3.
if large and not small then zoom = zoom
window,4,xsize=150+float(n_elements(lambda))*zoom,ysize=500,title='mospec1D',retain=2
device,retain=2
!p.charsize = 2.0
!p.thick = 1
!x.thick = 2
!y.thick = 2
!p.color=fsc_color('white')
!p.font = -1
!x.margin = [10,8]

if stackflag then begin
   xtitle='!6Rest wavelength (!6!sA!r!u!9 %!6!n)'  
   plot,lambda,nspec,xrange=xrange,yrange=[0,1.1*max(nspec)],xstyle=1,ystyle=1,psym=10,color=fsc_color('white'),yticks=1,ytickname=replicate(' ',2),yminor=1
   oplot,lambda,nspec,psym=10,color=fsc_color('turquoise')
   axis,yaxis=1,yminor=1,ytitle='Number of objects in stack',yrange=[0,1.1*max(nspec)],/ys
   plot,[min(xrange),max(xrange)],[0,0],xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,xtitle=xtitle,title=file,psym=10,color=fsc_color('white'),yticks=1,ytickname=replicate(' ',2),/noerase
   axis,yaxis=0,ytitle=textoidl("F_{\lambda} (10^{-17} erg/s/cm^2/!6!sA!r!u!9 %!6!n)"),ystyle=1
endif else begin
   xtitle='!6Observed wavelength (!6!sA!r!u!9 %!6!n)'  
   plot,[min(xrange),max(xrange)],[0,0],xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,xtitle=xtitle,ytitle=textoidl("F_{\lambda} (10^{-17} erg/s/cm^2/!6!sA!r!u!9 %!6!n)"),title=file,psym=10,color=fsc_color('white')
endelse
oplot,lambda,spec,psym=10,color=fsc_color('white')
oplot,lambda,errspec,psym=10,color=fsc_color('red')

lineid = 0
sval = 1.0
fitcheck = 0
mask = 0
redshift = 0.0
if stackflag then begin
   fixredshift = 0.0
   zflag = 1
endif else zflag = 0
vflag = 0

key = ''
while key ne 'q' do begin
      
   key = get_kbrd()

   ;;;;;
   if key eq ':' then begin
      print,':',format='(A,$)'
      key = get_kbrd()
      if key eq 'z' then begin
         print,key,format='(A,$)'
         zflag = 1
         read,fixredshift,prompt=''
         key = ''
      endif
      if key eq 'v' then begin
         print,key,format='(A,$)'
         vflag = 1
         read,fixvelocity,prompt=''
         key = ''
      endif
      if key eq 'c' then begin
         print,key,format='(A,$)'
         zflag = 0
         vflag = 0
         key = get_kbrd()
         print,''
         print,'Redshift and line width will no longer be fixed during fitting'
         key = ''
      endif
   endif
      
   ;;;;; 
   if key eq ' ' then begin
      cursor,x,y,/data,/nowait
      mindiff = min(abs(lambda-x),index)
      print,x,spec[index]
   endif
     
   ;;;;;  
   if key eq 'a' then begin
      cursor,x1,y1,/data,/nowait
      xyouts,0.15,0.82,"Press 'a' again",/normal,color=fsc_color('white')
      key = get_kbrd()
      cursor,x2,y2,/data,/nowait
      xrange = minmax([x1,x2])
      yrange = minmax([y1,y2])
      if x1 eq x2 and y1 eq y2 then begin
         xrange = minmax(lambda)
         yrange = [0,1.1*max(spec[where(finite(spec))])] & yrange[0] = -1*yrange[1]/3.
      endif
      key = 'r'
   endif

   ;;;;;
   if key eq 'b' then begin
      if not bootflag then begin
         bootflag = 1 
         print,'Error on line-fitting parameters will be estimated using bootstrap resampling'
      endif else begin
         bootflag = 0
         print,'Errors on line-fitting parameters will be calculated formally'
      endelse
   endif
     
   ;;;;; 
   if key eq 'c' then begin
      spec = im1[*,0]
      if formal or not stackflag then errspec = im1[*,1] else errspec = stddev[*,0]
      lineid = 0
      sval = 1.0
      fitcheck = 0
      zflag = 0
      vflag = 0
      xrange = [min(lambda),max(lambda)]
      yrange = [-0.35*max(spec[where(finite(spec))]),1.1*max(spec[where(finite(spec))])]
      key = 'r'
   endif

   ;;;;;
   if key eq 'f' then begin
      oldredshift = redshift
      cursor,fx,fy,/data,/nowait
      mindiff = min(abs(lambda-fix(fx)),startindex)
      startline = ''
      read,startline,prompt='Rest wavelength of line: '

      if stackflag then begin

         fitrange = ''
         read,fitrange,prompt='Range of rest wavelength to fit: '
         fitrange = fix(strsplit(fitrange,' ',/extract))
         
         if not nomodel then begin
            modelflag = 1
            ;; ebmv = sxpar(h1,'EBMV')
            ;;    modelspec = readfits(calibdir+'bpass_m004_t8.fnu.fits',modhead,exten=0,/silent)
            ;;    modelspec = readfits(calibdir+'t8_csf.fits',modhead,exten=0,/silent)
            ;;    modellambda = dblarr((size(modelspec))[1])
            ;;    m_refpix = sxpar(modhead,'CRPIX1')
            ;;    m_lam0 = sxpar(modhead,'CRVAL1')
            ;;    m_delta = sxpar(modhead,'CD1_1')
            ;;    for n=long(0),long(n_elements(modellambda))-1 do modellambda[n]=m_lam0+(n-(m_refpix-1))*m_delta
            readcol,calibdir+'tau0.0.10150',modellambda,modelspec,format='D,D',/silent
            airtovac,modellambda
            ;; modelspec *= (2.99792458e18)/(modellambda*modellambda)
            modelspec = modelspec[where(modellambda ge 2000 and modellambda le 8000)]
            modellambda = modellambda[where(modellambda ge 2000 and modellambda le 8000)]
            
            ebmv = findgen(301)/100
            chi2 = dblarr(n_elements(ebmv))
            for e=0,n_elements(ebmv)-1 do begin
               redmodelspec = modelspec
               for n=0,n_elements(redmodelspec)-1 do redmodelspec[n] = modelspec[n]*10^(-k_lambda(modellambda[n],/calzetti)*ebmv[e]/2.5)
               redmodel_spl = spl_init(modellambda,redmodelspec,/double)
               redmodel_int = spl_interp(modellambda,redmodelspec,redmodel_spl,lambda,/double)
               scale = median(spec[where(~lineflag)],/even)/median(redmodel_int[where(~lineflag)],/even)
               redmodel_int *= scale
               chi2[e] = total((redmodel_int[where(~lineflag and nspec gt ceil(0.5*max(nspec)))]-spec[where(~lineflag)])^2./errspec[where(~lineflag and nspec gt ceil(0.5*max(nspec)))])
            endfor
            print,ebmv[where(chi2 eq min(chi2))]
            
            redmodelspec = modelspec
            for n=0,n_elements(redmodelspec)-1 do redmodelspec[n] = modelspec[n]*10^(-k_lambda(modellambda[n],/calzetti)*ebmv[where(chi2 eq min(chi2))]/2.5)
            redmodel_spl = spl_init(modellambda,redmodelspec,/double)
            redmodel_int = spl_interp(modellambda,redmodelspec,redmodel_spl,lambda,/double)
            scale = median(spec[where(~lineflag)],/even)/median(redmodel_int[where(~lineflag)],/even)
            redmodel_int *= scale
            redmodelspec *= scale
            modelspec = redmodelspec
            modelspec_int = redmodel_int
         endif
         
         if modelflag then zeroindex = where(abs(spec) gt 0 and errspec gt 0 and finite(errspec) and finite(spec) and lambda ge fitrange[0] and lambda le fitrange[1] and lambda ge min(modellambda) and lambda le max(modellambda))
         if not modelflag then zeroindex = where(abs(spec) gt 0 and errspec gt 0 and finite(errspec) and finite(spec) and lambda ge fitrange[0] and lambda le fitrange[1])
      endif else begin
         zeroindex = where(abs(spec) gt 0 and errspec gt 0 and finite(errspec) and finite(spec))
         if fire then begin
            fitrange = ''
            read,fitrange,prompt='Range of rest wavelength to fit: '
            fitrange = fix(strsplit(fitrange,' ',/extract))
            zeroindex = where(abs(spec) gt 0 and errspec gt 0 and finite(errspec) and finite(spec) and lambda/(lambda[startindex]/startline) ge fitrange[0] and lambda/(lambda[startindex]/startline) le fitrange[1])
         endif
      endelse

      if not modelflag or not stackflag then begin
         line_fit_mc,lambda[zeroindex],spec[zeroindex],errspec[zeroindex],lambda[startindex],spec[startindex],startline,result=result,perror=perror,covar=covar,redxi2=redxi2,zflag=zflag,vflag=vflag,fixredshift=fixredshift,fixvelocity=fixvelocity,bootflag=bootflag,bootfile=file,stackflag=stackflag,full=full,cpoly=cpoly,twocomp=twocomp
         cfit = poly(lambda,result[0:cpoly])
         redshift = result[cpoly+1]
      endif
      
      if modelflag and stackflag then begin
         ;; line_fit_mc,lambda[zeroindex],spec[zeroindex],errspec[zeroindex],lambda[startindex],spec[startindex],startline,zflag=zflag,vflag=vflag,fixredshift=fixredshift,fixvelocity=fixvelocity,contmed=med_cont,cpoly=cpoly,stackflag=stackflag,/silent,full=full
         ;; modelspec_int *= med_cont/median(redmodelspec[where(modellambda ge min(lambda) and modellambda le max(lambda))])
         ;; modelspec *= med_cont/median(redmodelspec[where(modellambda ge min(lambda) and modellambda le max(lambda))])
         ;; if cpoly gt 0 then cpolyorig = cpoly
         cpoly = -1
         cfit = modelspec_int
         
         line_fit_mc,lambda[zeroindex],spec[zeroindex],errspec[zeroindex],lambda[startindex],spec[startindex],startline,result=result,perror=perror,covar=covar,redxi2=redxi2,zflag=zflag,vflag=vflag,fixredshift=fixredshift,fixvelocity=fixvelocity,contmed=med_cont,sedcont=cfit[zeroindex],cpoly=cpoly,stackflag=stackflag,bootfile=file,bootflag=bootflag,full=full,twocomp=twocomp
         redshift = result[cpoly+1]
      endif
      yfit = cfit
      if not twocomp then for n=0,(n_elements(result)-(cpoly+3))/2.-1 do yfit += result[2*n+cpoly+3]*exp(-(lambda-(1+result[cpoly+1])*result[2*n+cpoly+4])^2/(2*(result[cpoly+2]*(1+result[cpoly+1])*result[2*n+cpoly+4])^2))
      if twocomp then begin
         yfit1 = 0.0
         for n=0,(n_elements(result)-(cpoly+5))/4.-1 do yfit1 += result[2*n+cpoly+3]*exp(-(lambda-(1+result[cpoly+1])*result[2*n+cpoly+4])^2/(2*(result[cpoly+2]*(1+result[cpoly+1])*result[2*n+cpoly+4])^2))
         tmpn = n
         yfit2 = 0.0
         for n=tmpn,(n_elements(result)-(cpoly+5))/2.-1 do yfit2 += result[2*n+cpoly+3]*exp(-(lambda-(1+result[-1])*result[2*n+cpoly+4])^2/(2*(result[-2]*(1+result[-1])*result[2*n+cpoly+4])^2))
         yfit = cfit+yfit1+yfit2
      endif
      if keyword_set(cpolyorig) then cpoly = cpolyorig

      lambda *= (oldredshift+1.0)/(redshift+1.0)
      xrange *= (oldredshift+1.0)/(redshift+1.0)
      fitcheck = 1
      lineid = 1
      key = 'r'
   endif

   ;;;;;
   if key eq 'g' then begin
      if not keyword_set(startline) then print,"Please compute fit to spectrum using 'f' first"
      openw,1,'mospec1d.log',/append
      printf,1,systime()
      printf,1,file
      if not modelflag or not stackflag then line_fit_mc,lambda[zeroindex],spec[zeroindex],errspec[zeroindex],lambda[startindex],spec[startindex],startline,result=result,perror=perror,covar=covar,redxi2=redxi2,zflag=zflag,vflag=vflag,fixredshift=fixredshift,fixvelocity=fixvelocity,bootflag=bootflag,bootfile=file,stackflag=stackflag,full=full,cpoly=cpoly,twocomp=twocomp,/print
      if modelflag and stackflag then line_fit_mc,lambda[zeroindex],spec[zeroindex],errspec[zeroindex],lambda[startindex],spec[startindex],startline,result=result,perror=perror,covar=covar,redxi2=redxi2,zflag=zflag,vflag=vflag,fixredshift=fixredshift,fixvelocity=fixvelocity,contmed=med_cont,sedcont=cfit[zeroindex],cpoly=-1,stackflag=stackflag,bootfile=file,bootflag=bootflag,full=full,twocomp=twocomp,/print
      printf,1,''
      close,1
      print,'Results of line fitting printed to mospec1d.log'
   endif

   ;;;;; 
   if key eq 'h' then begin
      if not fitcheck then fitcheck = 1 else fitcheck = 0
      key = 'r'
   endif
   
   ;;;;;
   if key eq 'i' then begin
      if not lineid then lineid = 1 else lineid = 0
      key = 'r'        
   endif

   ;;;;;
   if key eq 'p' then begin
      tmpxrange = ''
      read,tmpxrange,prompt='Wavelength range to plot (Angstrom): '
      if min(tmpxrange) lt min(lambda) or max(tmpxrange) gt max(lambda) then read,tmpxrange,prompt='Invalid wavelength range. New range: '
      tmpyrange = 1.0
      read,tmpyrange,prompt='Maxmimum flux to plot (1e-17 erg/s/cm^2/Ang): '
      label = ''
      read,label,prompt='Plot label: '
      tmpyrange = [-1*tmpyrange[0]/3.,1.1*tmpyrange[0]]
      tmpxrange = fix(strsplit(tmpxrange[sort(tmpxrange)],' ',/extract))
      ;; mindiff = min(abs(lambda-tmpxrange[0]),tmpindex1)
      ;; mindiff = min(abs(lambda-tmpxrange[1]),tmpindex2) 
      ;; tmpxrange = [lambda[tmpindex1],lambda[tmpindex2]]

      set_plot,'PS'
      !p.charsize = 1
      !p.font = 1
      !p.charthick = 4
      !p.thick = 2
      !x.thick = 4
      !y.thick = 4
      xsize = 8.8;6.5*2.54/3
      ;; ysize = 9*2.54/3
      plot_params,xsize,ysize,position
      filename = (strsplit(file,'.',/extract))[0]+'.eps'
      device,filename=filename,xsize=xsize,ysize=ysize,/color,/encapsulated,font_size=12,/tt_font,decomposed=0
      angstrom = '!3' + STRING(197B) + '!X'
      if redshift gt 0 then xtitle='Observed wavelength ('+angstrom+')'
      if redshift eq 0 then xtitle='Rest wavelength ('+angstrom+')'

      if tmpxrange[1]-tmpxrange[0] lt 200 then xtickinterval=10 else xtickinterval=100
      plot,lambda,spec,xrange=tmpxrange,yrange=[-1*tmpyrange[1]/20.,tmpyrange[1]],xstyle=1,ystyle=1,psym=10,color=fsc_color('black'),xtitle=xtitle,ytitle=textoidl("F_{\lambda} (10^{-17} erg/s/cm^2/"+angstrom+")"),xminor=4,yminor=4,/normal,xmargin=[7,2],ymargin=[4,1],position=position,xtickinterval=xtickinterval
      oplot,!X.CRANGE,[0,0],thick=1
      if not mask then oplot,lambda,errspec,psym=10,color=fsc_color('red'),thick=1
      if modelflag then oplot,lambda,redmodel_int,psym=10,color=fsc_color('green'),thick=2
      
      xyouts,0.91*(!X.CRANGE[1]-!X.CRANGE[0])+!X.CRANGE[0],0.8*(!Y.CRANGE[1]-!Y.CRANGE[0])+!Y.CRANGE[0],label,/data,color=fsc_color('black'),align=1      
      ;xyouts,0.11*(!X.CRANGE[1]-!X.CRANGE[0])+!X.CRANGE[0],0.8*(!Y.CRANGE[1]-!Y.CRANGE[0])+!Y.CRANGE[0],label,/data,color=fsc_color('black')
      
      if fitcheck eq 1 then begin
         oplot,lambda,yfit,color=fsc_color('sky blue'),thick=4
         if keyword_set(yfit1) then oplot,lambda,yfit1,color=fsc_color('magenta'),linestyle=2
         if keyword_set(yfit2) then oplot,lambda,yfit2,color=fsc_color('orange'),linestyle=2
      endif
      if modelflag and fitcheck then oplot,lambda,cfit,psym=10,color=fsc_color('lime green'),thick=4
      
      if lineid then begin
         for n=0,n_elements(lines)-1 do begin
            if ((linewav[n] gt !X.CRANGE[0]) and (linewav[n] lt !X.CRANGE[1])) then begin
               vline,linewav[n],linestyle=1,color=fsc_color('dark green')
               printline = lines[n]
               if printline eq 'H-alpha' then printline = textoidl('H\alpha')
               if printline eq 'H-beta' then printline = textoidl('H\beta')
               if printline eq 'H-eps' then printline = textoidl('H\varepsilon')
               if printline eq '[OII]' and lines[n+1] eq '[OII]' then printline = ''
               if printline eq '[NeIII]' and lines[n+1] eq 'H-eps' then begin
                  xyouts,linewav[n]-10,0.85*!Y.CRANGE[1],printline,/data,color=fsc_color('dark green'),orient=270,charsize=0.9
               endif else begin
                  xyouts,linewav[n]+5,0.85*!Y.CRANGE[1],printline,/data,color=fsc_color('dark green'),orient=270,charsize=0.9
               endelse
            endif
         endfor
      endif
      device,/close
      print,'Postscript figure written to '+(strsplit(file,'.',/extract))[0]+'.eps'
      spawn,'open '+(strsplit(file,'.',/extract))[0]+'.eps'
      set_plot,'X'
      !p.charsize = 2
      !p.thick = 1
      !x.thick = 2
      !y.thick = 2
      !p.color=fsc_color('white')
      !p.font = -1
      !p.charthick = 1
      key = ''
   endif

   ;;;;; 
   if key eq 's' then begin
      xyouts,0.15,0.82,'Choose a smoothing value between 1 and 9',/normal,color=fsc_color('white')
      xyouts,0.15,0.77,'A value of 1 returns the original spectrum',/normal,color=fsc_color('white')
      sval = fix(1)
      sval = get_kbrd()
      spec = gauss_smooth(im1[*,0],width=fix(sval),/nan)
      key = 'r'
   endif

   ;;;;; 
   if key eq 'u' then begin
      if not mask then mask = 1 else mask = 0
      key = 'r'
   endif
   
   ;;;;; 
   if key eq 'z' then begin
      oldredshift = redshift
      read,redshift,prompt='Redshift: '
      lambda *= (oldredshift+1.0)/(redshift+1.0)
      xrange *= (oldredshift+1.0)/(redshift+1.0)
      lineid = 1
      key = 'r'
   endif

   ;;;;; 
   if key eq ',' then begin
      yrange *= 2.
      yrange[0] = -1*yrange[1]/3.
      key = 'r'
   endif
   
     ;;;;; 
   if key eq '.' then begin
      yrange /= 2.
      yrange[0] = -yrange[1]/3.
      key = 'r'
   endif
     
   ;;;;; 
   if key eq 'r' then begin
      erase      
      plot,[min(xrange),max(xrange)],[0,0],xrange=xrange,yrange=yrange,xstyle=1,xtitle=xtitle,ytitle=textoidl("F_{\lambda} (10^{-17} erg/s/cm^2/!6!sA!r!u!9 %!6!n)"),title=file,psym=10,color=fsc_color('white')
      oplot,lambda,spec,psym=10,color=fsc_color('white')
      if not mask then oplot,lambda,errspec,psym=10,color=fsc_color('red')

      !p.thick = 2
      if fitcheck eq 1 then begin
         oplot,lambda,yfit,color=fsc_color('orange'),thick=2
         if keyword_set(yfit1) then oplot,lambda,yfit1,color=fsc_color('magenta'),linestyle=2
         if keyword_set(yfit2) then oplot,lambda,yfit2,color=fsc_color('cyan'),linestyle=2
         oplot,lambda,cfit,color=fsc_color('black'),thick=2
         oplot,lambda,cfit,color=fsc_color('yellow'),thick=2
         if stackflag then vline,fitrange,color=fsc_color('blue')
      endif
      !p.thick = 1
      
      if lineid then begin
         for n=0,n_elements(lines)-1 do begin
            vline,linewav[n],linestyle=4,color=fsc_color('green')
            if ((linewav[n] gt !X.CRANGE[0]) and (linewav[n] lt !X.CRANGE[1])) then $
               xyouts,linewav[n]+1,yrange[0],lines[n],/data,color=fsc_color('green'),orient=270,align=1
         endfor
      endif

      if modelflag then oplot,lambda,redmodel_int,psym=10,color=fsc_color('green'),thick=2
      
      !p.color = fsc_color('white')
   endif
   
   if key eq 'x' then stop
   
   ;;;;; 
   if key eq '?' then begin
      print,"spacebar - print current location of cursor to command line"
      print,"a - zoom; press once in lower left corner of desired region and once again in the upper right corner; press twice in the same spot to unzoom"
      print,"b - turn bootstrapping error estimates on and off"
      print,"c - clear all changes to spectrum and replot"
      print,"d - toggle whether or not the fit is displayed"
      print,"f - calculate simultaneous fit to the continuum and strong emission lines"
      print,"g - write line fluxes and fit parameters to mospec.log"
      print,"i - identify emission lines using the default list or lines.dat"
      print,"p - prints figure of current object to file"
      print,"q - quit"
      print,"r - redraw the spectrum with current zoom, smoothing, and line list"
      print,"s - smooth spectrum (must choose a value from 1-9)"
      print,"x - pauses the IDL routine; to restart program, type '.c'"
      print,"? - help"
   endif
   
endwhile

if key eq 'q' then return

end
