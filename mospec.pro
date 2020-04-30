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

pro mospec,file,all=all,small=small,large=large,nocal=nocal,nomodel=nomodel,cpoly=cpoly,skip=skip,twocomp=twocomp

if n_params() lt 1 and not keyword_set(all) then begin
   print,'Syntax - MOSPEC, file, [/all, /small, /large, /nocal, /nomodel, cpoly=]'
   return
endif

resolve_all,/quiet,/continue_on_error,skip='rsex'

if keyword_set(all) then all=1 else all=0
if keyword_set(small) then small=1 else small=0
if keyword_set(large) then large=1 else large=0
if keyword_set(nocal) then nocal=1 else nocal=0
if keyword_set(nomodel) then nomodel=1 else nomodel=0
if n_elements(cpoly) eq 0 then cpoly = 1
if keyword_set(skip) then skip=1 else skip=0
if keyword_set(twocomp) then twocomp=1 else twocomp=0

common moscount,countnum,mosfilter
if not keyword_set(mosfilter) then mosfilter = ''

calibdir = getenv('MOSPEC_CALIB')
if strmid(calibdir,n_elements(calibdir)-1) ne '/' then calibdir += '/'
if file_search('lines.dat') ne '' then readcol,'lines.dat',lines,airwav,vacwav,format='A,D',/silent else readcol,calibdir+'rest_optical_lines_vac_part.dat',lines,airwav,vacwav,format='A,D,D',/silent
linewav=vacwav
hlines = linewav[where(strpos(lines,'H') eq 0 and strpos(lines,'e') ne 1)]

if all then begin
   filter = ''
   read,filter,prompt='Filter: '
   if mosfilter ne filter then countnum = 0
   mosfilter = filter
   files = file_search('*_'+strupcase(filter)+'_wtavg.fits')
   if files[0] eq '' then files = file_search('*_'+strupcase(filter)+'_*flx.fits',/fold_case) else files = [files,file_search('*_'+strupcase(filter)+'_*flx.fits',/fold_case)]
   if files[0] eq '' then files = file_search('*_'+strupcase(filter)+'_*eps.fits',/fold_case) else files = [files,file_search('*_'+strupcase(filter)+'_*eps.fits',/fold_case)]
endif else begin
   countnum = 0
   if strmid(file,3,/reverse) eq 'fits' then files=file else files = [file_search('*-'+file+'_*wtavg.fits'),file_search('*'+file+'_*flx.fits'),file_search('*'+file+'*_eps.fits')]
   if strlowcase(file) eq 'wtavg' then files = [file_search('*-*_*wtavg.fits')]
endelse 
files = files[where(files ne '')]
if files[0] eq '' then message,'Cannot find any matching files',/reset

for i=countnum,n_elements(files)-1 do begin
   countnum = i

   file = files[i]
   im1 = readfits(file,h1,exten=0,/silent)
   object = strtrim(sxpar(h1,'OBJECT'))
   filter = strtrim(sxpar(h1,'FILTER'))
   serendip = ''
   if strpos(file,'wtavg') eq -1 then outfile1d = strmid(file,0,strpos(file,'_',/reverse_search))+serendip+'_1d.fits' $
   else outfile1d = object+serendip+'_'+filter+'_1d.fits'
   if file_test('../onedspec/'+outfile1d) and skip then continue

   if all then print,''
   counter = '--- '+strtrim(string(countnum+1),1)+'/'+strtrim(string(n_elements(files)),1)+' ---'
   if all then print,counter
   print,'Reading in '+file
   dim = size(im1,/dimen)
   if file_search(strmid(file,0,strpos(file,'fits'))+'sig.fits') ne '' then begin
      im2 = readfits(strmid(file,0,strpos(file,'fits'))+'sig.fits',h2,exten=0,/silent) 
      print,'Reading in '+strmid(file,0,strpos(file,'fits'))+'sig.fits'
   endif
   if file_search(strmid(file,0,strpos(file,'eps'))+'ivar.fits') ne '' then begin
      im2 = readfits(strmid(file,0,strpos(file,'eps'))+'ivar.fits',h2,exten=0,/silent)
      print,'Reading in '+strmid(file,0,strpos(file,'eps'))+'ivar.fits'
      im2 = 1.0/im2
      im2 = sqrt(im2)
   endif
   if file_search(strmid(file,0,strpos(file,'eps'))+'sig.fits') ne '' then begin
      im2 = readfits(strmid(file,0,strpos(file,'eps'))+'sig.fits',h2,exten=0,/silent) 
      print,'Reading in '+strmid(file,0,strpos(file,'eps'))+'sig.fits'
   endif
   if not keyword_set(im2) then message,'Cannot find valid error or inverse variance spectrum'
   index = where(~finite(im2) and abs(im1) gt 0)
   if index[0] ne -1 then print,'WARNING: non-finite values found in 2D error spectrum or inverse variance map';,' ',n_elements(index)
   im2[index] = 0.0
   index = where(~finite(im1))
   im1[index] = 0.0
   im2[index] = 0.0

   scale = im1[where(finite(im1))]
   scale = prank(scale,[3,97])
   loadct,0
   ;; tvlct,r,g,b,/get
   ;; tvlct,reverse(r),reverse(g),reverse(b)

   field = (strsplit(object,'-',/extract))[0]
   unit = strtrim(sxpar(h1,'BUNIT'))
   header_ap1 = float(sxpar(h1,'MOSAP1'))
   header_ap2 = float(sxpar(h1,'MOSAP2'))
   if strpos(file,'wtavg') ge 0 then logmask = 'wtavg' else logmask = sxpar(h1,'MASKNAME')+'_'+strupcase(filter) ;strmid(file,0,strpos(file,(strsplit(object,'-',/extract))[1])-1)
   if nocal or (unit eq 'FLUX' or unit eq 'ergs/cm2/s/A') then cal = replicate(1.0,dim[0])
   if not nocal and unit eq 'electron/second' or unit eq 'ELECTRONS/S' then begin
      if filter eq 'K' then readcol,calibdir+'Kcal.comb.final.dat',callam,calflam,format='D,D',/silent
      if filter eq 'H' then readcol,calibdir+'Hcal_nov13.dat',callam,calflam,format='D,D',/silent
      if filter eq 'J' then readcol,calibdir+'Jcal_nov2013.dat',callam,calflam,format='D,D',/silent
      if filter eq 'Y' then readcol,calibdir+'Yspec.cal.drp.psp.dat',callam,calflam,format='D,D',/silent
      calfnu = callam^2*calflam/(3e18)
      cal = calflam/(1.0e-17)
   endif
   if not keyword_set(cal) then message,"Data units ('"+unit+"') not recognized"

   if header_ap1 gt 0 and header_ap2 gt 0 then line = [header_ap1,header_ap2]
   if header_ap1 eq 0 or header_ap2 eq 0 or header_ap1 eq header_ap2 then line = [fix(abs(sxpar(h1,'CRVAL2')))-3,fix(abs(sxpar(h1,'CRVAL2')))+3]
   boxcar = 1
   gaussian = 0
   empirical = 0
   test1 = total(im1[*,line[0]:line[1]],2)*cal
   test2 = im2*im2
   test2 = sqrt(total(test2[*,line[0]:line[1]],2))*cal

   obslambda = dblarr(dim[0])
   refpix = sxpar(h1,'CRPIX1')
   lam0 = sxpar(h1,'CRVAL1')
   delta = sxpar(h1,'CD1_1')
   for n=long(0),long(n_elements(obslambda))-1 do obslambda[n]=lam0+(n-(refpix-1))*delta

   if file_test(calibdir+'oh_lines_'+filter+'.txt') then begin
      readcol,calibdir+'oh_lines_'+filter+'.txt',skylam,f1,f2,f3,f4,f5, format='F,F,F,A,F,F,X',/silent
      tempindex=fltarr(n_elements(obslambda))
      temp2index=fltarr(n_elements(obslambda))
      for j=0, n_elements(skylam)-1 do begin
         index=where(abs(obslambda-skylam[j]) lt 6.0)
         tempindex[index]=1.0
         index=where(abs(obslambda-skylam[j]) lt 10.0)
         temp2index[index]=1.0
      endfor
      skyindex=where(tempindex eq 1.0)
      smoothskyindex=where(temp2index eq 1.0)
   endif
   if not file_test(calibdir+'oh_lines_'+filter+'.txt') then begin
      skyindex = []
      smoothskyindex = []
   endif
   
   objspec = test1
   objspec[skyindex] = 0.0/0
   
   lambda = obslambda
   spec = test1
   errspec = test2
   xrange = [min(lambda),max(lambda)]
   yrange = [0,1.1*max(test1[where(finite(test1))])] & yrange[0] = -1*yrange[1]/3.
   rowrange = [0,dim[1]]

   modelflag = 0
   seddir = getenv("MOSPEC_SED")
   sedfile = seddir+'bc03calz/'+strlowcase(field)+'/bestfit/bestfit.'+object+'.csf_agegt50.dat'
   if file_test(sedfile) and not nomodel then begin
      modelflag = 1
      print,'Reading in bestfit.'+object+'.csf_agegt50.dat'
      readcol,sedfile,modellambda,modelspec,format='D,D',/silent
      airtovac,modellambda
      modelspec /= (1.0e-17)
      photfile = seddir+'inphot/'+strlowcase(field)+'_inphot.dat'
      openr,1,photfile
      sed_object = strarr(file_lines(photfile))
      sed_redshift = dblarr(file_lines(photfile))
      for l=0,file_lines(photfile)-1 do begin
         tmps = ''
         readf,1,tmps
         tmps = strsplit(tmps,' ',/extract)
         sed_object[l] = tmps[0]
         sed_redshift[l] = tmps[-1]
      endfor
      close,1
      index = (where(sed_object eq object))[0]
      if index eq -1 then begin
         print,'No redshift found in SED inphot file'
         modelflag = 0
      endif else modelz = sed_redshift[index]
   endif
   kbss = read_kbss()
   kbss = kbss[where(kbss.field+'-'+kbss.obj eq object,foundobj)]
   if modelflag and foundobj and kbss.lmstar gt 0 then begin
      age_string = string(kbss.age/1000)
      age_string = strtrim(strmid(age_string,0,strpos(age_string,'00')),1)
      sedfile_highres = file_search(seddir+'highres_ssp/tau0/tau*.'+age_string+'*')
      if n_elements(sedfile_highres) gt 1 then begin
         age_highres = strmid(sedfile_highres,strpos(sedfile_highres[0],'.',strpos(sedfile_highres[0],'tau'))+1)
         mindiff = min(abs(double(age_highres)-kbss.age/1000),match)
         sedfile_highres = sedfile_highres[match]
      endif
      print,'Loading matching high-resolution model'
      readcol,sedfile_highres,highreslambda,highresspec,format='D,D',/silent
      highresspec = highresspec[where(highreslambda ge 3000 and highreslambda le 10000)]
      highreslambda = highreslambda[where(highreslambda ge 3000 and highreslambda le 10000)]
      airtovac,highreslambda
      for n=0,n_elements(highreslambda)-1 do highresspec[n] = highresspec[n]*10^(-k_lambda(highreslambda[n],/calzetti,/silent)*kbss.ebmv_sed/2.5)
      highreslambda *= (1.+modelz)
      med_highres = median(highresspec[where(highreslambda ge 5000*(1.+modelz) and highreslambda le 5700*(1.+modelz))],/even)
      med_lowres = median(modelspec[where(modellambda ge 5000*(1.+modelz) and modellambda le 5700*(1.+modelz))],/even)
      highresspec *= med_lowres/med_highres
      modellambda = highreslambda
      modelspec = highresspec
   endif
   
   if small then zoom = 1/2. else zoom = 2/3.
   if large and not small then zoom = 1.
   window,4,xsize=150+dim[0]*zoom,ysize=500+dim[1]*zoom,title='mospec',retain=2
   device,retain=2
   !p.charsize = 2.0
   !p.thick = 1
   !x.thick = 2
   !y.thick = 2
   !p.color=cgcolor('white')
   !p.font = -1
   twodspec = im1
  
   mask = 0
   lineid = 0
   redshift = 0.0
   sval = 1.0
   fitcheck = 0
   zflag = 0
   vflag = 0

   xtitle='!6Observed wavelength (!6!sA!r!u!9 %!6!n)'  
   plot,[0],[0],xrange=xrange,yrange=rowrange,position=[115,65,115+dim[0]*zoom,65+dim[1]*zoom],/device,xstyle=1,xticks=1,xtickname=replicate(' ',2),ystyle=1,yticks=1,ytickname=replicate(' ',2),color=cgcolor('white')
   tvimage,bytscl(twodspec,scale[0],scale[1]),/overplot
   plot,[0],[0],xrange=xrange,yrange=rowrange,position=[115,65,115+dim[0]*zoom,65+dim[1]*zoom],xtickformat='(I)',/device,xstyle=1,ystyle=1,yticks=1,ytickname=replicate(' ',2),xtitle=xtitle,color=cgcolor('white'),/noerase
   oplot,[0,xrange[1]],[line[0],line[0]],color=cgcolor('green')
   oplot,[0,xrange[1]],[line[1],line[1]],color=cgcolor('green')
   plot,[min(xrange),max(xrange)],[0,0],xrange=xrange,yrange=yrange,xstyle=1,title=object+serendip+' - '+file,ytitle=textoidl("F_{\lambda} (10^{-17} erg/s/cm^2/!6!sA!r!u!9 %!6!n)"),xticks=1,xtickname=replicate(' ',2),position=[115,65+dim[1]*zoom,115+dim[0]*zoom,450+dim[1]*zoom],/device,/noerase,psym=10,color=cgcolor('white')
   oplot,lambda,spec,psym=10,color=cgcolor('white')
   oplot,lambda,objspec,psym=10,color=cgcolor('white')
   oplot,lambda,errspec,psym=10,color=cgcolor('red')
   if file_test('../onedspec/'+outfile1d) then xyouts,115+dim[0]*zoom,460+dim[1]*zoom,'1D ',/align,/device,color=cgcolor('orange')
   
   if modelflag then begin
      print,'Best guess for redshift as used in SED-fitting: ',string(modelz,format='(F5.3)')
      print,'WARNING: do NOT fit using model continuum if you choose a different redshift!'
   endif
   key = ''
   while key ne 'q' and key ne 'n' and key ne 'b' do begin
      
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
        if key eq 's' then begin
           print,key,format='(A,$)'
           read,serendip,prompt=''
           serendip = strtrim(serendip,1)
           if serendip ne '' then object = field+'-'+serendip else object = sxpar(h1,'OBJECT')
           print,'Now measuring spectrum for object '+object
           if strpos(file,'wtavg') eq -1 then begin
              if serendip ne '' then outfile1d = logmask+'_'+serendip+'_1d.fits' else outfile1d = strmid(file,0,strpos(file,'_',/reverse_search))+'_1d.fits'
           endif else outfile1d = object+'_'+filter+'_1d.fits'
           if serendip ne '' then begin
              modelflag = 0
              print,'Disabling SED model continuum'
           endif
           if serendip eq '' and not nomodel then begin
              modelflag = 1
              print,'Re-enabling SED model continuum'
           endif
           key = 'c'
        endif
     endif

     ;;;;; 
     if key eq ' ' then begin
        cursor,x,y,/data,/nowait
        mindiff = min(abs(lambda-x),index)
        if y gt yrange[0] then print,x,spec[index]
        key = ''
        if y le yrange[0] then begin
           tvimage,bytscl(twodspec,scale[0],scale[1]),/overplot
           plot,[0],[0],xrange=xrange,yrange=rowrange,position=[115,65,115+dim[0]*zoom,65+dim[1]*zoom],xtickformat='(I)',/device,xstyle=1,ystyle=1,yticks=1,ytickname=replicate(' ',2),xtitle=xtitle,color=cgcolor('white'),/noerase
           cursor,x,y,/data,/nowait
           print,x,y
           key = 'r'
        endif
     endif
     
     ;;;;; 
     if key eq 'a' then begin
        cursor,x1,y1,/data,/nowait
        xyouts,0.11,0.85,"Press 'a' again",/normal,color=cgcolor('white')
        key = get_kbrd()
        cursor,x2,y2,/data,/nowait
        xrange = minmax([x1,x2])
        yrange = minmax([y1,y2])
        mindiff = min(abs(lambda-xrange[0]),index1)
        mindiff = min(abs(lambda-xrange[1]),index2)
        index = mean(line)
        rowrange = [max([fix(index-((index2-index1)*(dim[1])/dim[0])/2.),0]),min([fix(index+((index2-index1)*(dim[1])/dim[0])/2.),dim[1]-1])]
        twodspec = congrid(im1[index1:index2,rowrange[0]:rowrange[1]],dim[0],dim[1])

        if abs(x1-x2) le 1.0 and abs(y1-y2) le 0.01 then begin
           xrange = [min(obslambda),max(obslambda)]/(redshift+1.0)
           yrange = [0,1.1*max(test1[where(finite(test1))])] & yrange[0] = -1*yrange[1]/3.
           rowrange = [0,dim[1]]
           twodspec = im1[*,0:dim[1]-1]
        endif

        scale = twodspec[where(finite(twodspec))]
        scale = prank(scale,[3,97])

        key = 'r'
     endif
     
     ;;;;; 
     if key eq 'c' then begin
        twodspec = im1
        line = [fix(abs(sxpar(h1,'CRVAL2')))-3,fix(abs(sxpar(h1,'CRVAL2')))+3]
        test1 = total(im1[*,line[0]:line[1]],2)*cal
        test2 = im2*im2
        test2 = sqrt(total(test2[*,line[0]:line[1]],2))*cal
        objspec = test1
        objspec[skyindex] = 0.0/0
        lambda = obslambda
        spec = test1
        errspec = test2
        xrange = minmax(lambda)
        yrange = [0,1.1*max(test1[where(finite(test1))])] & yrange[0] = -1*yrange[1]/3.
        rowrange = [0,dim[1]]
        mask = 0
        lineid = 0
        redshift = 0.0
        sval = 1.0
        fitcheck = 0
        zflag = 0
        vflag = 0
        key = 'r'
     endif

     ;;;;; 
     if key eq 'd' then begin
        tmpimage = im1
        row = indgen(dim[1])
        flux = fltarr(n_elements(row))
        for n=0,n_elements(tmpimage[0,*])-1 do begin
           range = 50
           resistant_mean,tmpimage[fix((50-range/2)/100.*dim[0]):fix((50+range/2)/100.*dim[0]),n],3.,fluxavg,/double
           flux[n] = fluxavg
        endfor
        index = where(finite(flux))
        row = row[index]
        flux = flux[index]
        
        guess = row[where(flux eq max(flux[where(row ge line[0] and row le line[1])]))]
        expr = 'P[1]+P[2]*exp(-(X-P[3])^2/(2*P[0]^2))+P[4]*exp(-(X-P[5])^2/(2*P[8]^2))+P[6]*exp(-(X-P[7])^2/(2*P[9]^2))'
        start = [1.5,0.D,max(flux),guess,min(flux),guess-17,min(flux),guess+17,1.5,1.5]
        result = mpfitexpr(expr,row[where(row ge line[0]-1.5*17 and row le line[1]+1.5*17)],flux[where(row ge line[0]-1.5*17 and row le line[1]+1.5*17)],flux*0+1.,start,yfit=yfit,/quiet,status=status,pi=pi,errmsg=errmsg)
        if status le 0 then message,errmsg
        seeing = 2*sqrt(2*alog(2))*result[0]*.18
        maskname = sxpar(h1,'MASK1')
        if string(maskname) eq string(long(0)) then maskname = sxpar(h1,'MASKNAME')

        if filter eq 'Y' then slitcor_fitlers = ['Y']
        if filter eq 'J' then slitcor_filters = ['J','J2','J3']
        if filter eq 'H' then slitcor_filters = ['H','H1','H2']
        if filter eq 'K' then slitcor_filters = ['Ks']

        for n=0,n_elements(slitcor_filters)-1 do begin
           readcol,calibdir+'mosfire_'+slitcor_filters[n]+'.txt',filtlam,filtthru,format='D,D',/silent
           filtthru = filtthru[sort(filtlam)]
           filtlam = filtlam[sort(filtlam)]*10000
           filtthru_spl = spl_init(filtlam,filtthru,/double)
           index = where(obslambda ge min(filtlam) and obslambda le max(filtlam) and abs(spec) gt 0)
           filtthru_interp = spl_interp(filtlam,filtthru,filtthru_spl,obslambda[index],/double)
           filt_wgts = filtthru_interp/total(filtthru_interp)

           starlambda = obslambda[index]
           starerr_fnu = (errspec[index]*1.0e-17)*lambda[index]*lambda[index]/(2.99792458e18)
           starspec_fnu = (spec[index]*1.0e-17)*lambda[index]*lambda[index]/(2.99792458e18)

           y = (spec*1.0e-17)*lambda*lambda/(2.99792458e18)/(1.e-26)
           yerr = (errspec*1.0e-17)*lambda*lambda/(2.99792458e18)/(1.e-26)
           plot,obslambda,y,psym=10,color=cgcolor('white'),xrange=!X.CRANGE,/xs,yrange=[-1./3,1.7]*max(y),xtitle=xtitle,ytitle='mJy'
           oplot,!X.CRANGE,[0,0]
           oplot,obslambda,yerr,psym=10,color=cgcolor('red')
           oplot,filtlam,filtthru/max(filtthru)*0.9*max(y),color=cgcolor('green')
           xyouts,0.18,0.85,"Seeing estimated to be "+strtrim(string(seeing,format='(F0.2)'),1)+"'' on "+maskname,/normal,color=cgcolor('white')

           starspec_fnu_wgt = starspec_fnu*filt_wgts
           resistant_mean,starspec_fnu,3.,specmean,/double,goodvec=sigclip
           star_fnu = total(starspec_fnu_wgt[sigclip])
;;            star_fnu_var = total(filt_wgts[sigclip]*(starspec_fnu[sigclip]-star_fnu)^2.)/(1.-total(filt_wgts[sigclip]^2.)) ;; unbiased estimate of sample variance (http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance)
           star_fnu_var = total(filt_wgts[sigclip]^2.*starerr_fnu[sigclip]^2.) ;; variance of the weighted mean (http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Statistical_properties)
           mag_ab = -2.5*alog10(star_fnu)-48.6
           mag_ab_err = 2.5*(sqrt(star_fnu_var)/star_fnu)/alog(10.) ;; remove sqrt(n_elements(sigclip) when using variance of the mean
           xyouts,0.18,0.80,textoidl("<f_{\nu}> = "+strtrim(string(star_fnu,format='(e11.3)'),1)+" erg/s/cm/Hz     m_{AB} = "+strtrim(string(mag_ab,format='(F0.3)'),1)),/normal,color=cgcolor('white')
           oplot,!X.CRANGE,[star_fnu,star_fnu]/(1.e-26)
           oplot,!X.CRANGE,replicate(median(starspec_fnu),2)/(1.e-26),color=cgcolor('magenta')

           aperture = fix(line[1]-line[0]+1)
           slitcor_star = 1./(erf(0.247487/(result[0]*0.18))*erf(aperture/(2.*sqrt(2.)*result[0])))
           mag_cor = -2.5*alog10(star_fnu*slitcor_star)-48.6
           xyouts,0.11,0.75,slitcor_filters[n],/normal,color=cgcolor('green'),charsize=6
           xyouts,0.18,0.75,textoidl("Geometric slit correction: "+strtrim(string(slitcor_star,format='(F0.3)'),1)+"    m_{cor} = "+strtrim(string(mag_cor,format='(F0.3)'),1)),/normal,color=cgcolor('white')
           xyouts,0.18,0.70,"Press any key to continue",/normal,color=cgcolor('white')

           logfile = '../slitstars.'+strupcase(field)
           openw,1,logfile,/append
           printf,1,logmask,slitcor_filters[n],object,seeing,mag_ab,mag_cor,slitcor_star,aperture,format='(%"%s %s %s %f %f %f %f %d")'
           close,1
           key = get_kbrd()

           ;; xyouts,0.09,0.75,slitcor_filters[n],/normal,color=cgcolor('black'),charsize=6
           ;; xyouts,0.18,0.80,textoidl("<f_{\nu}> = "+strtrim(string(star_fnu,format='(e11.3)'),1)+" erg/s/cm/Hz     m_{AB} = "+strtrim(string(mag_ab,format='(F0.3)'),1)),/normal,color=cgcolor('black')
           ;; xyouts,0.18,0.75,textoidl("Geometric slit correction: "+strtrim(string(slitcor_star,format='(F0.3)'),1)+"    m_{cor} = "+strtrim(string(mag_cor,format='(F0.3)'),1)),/normal,color=cgcolor('black')
           ;; xyouts,0.18,0.70,"Press any key to continue",/normal,color=cgcolor('black')

        endfor
        print,'Seeing and slit star data written to slitstars.'+strupcase(field)
        key = 'r'        
     endif

     ;;;;;
     if key eq 'e' then begin
        write1d = 0
        if file_test('../onedspec/'+outfile1d) then begin
           printcheck = ''
           read,printcheck,prompt="There is already a 1d extracted spectrum for this object. Continue [y/n]? "
           if strlowcase(printcheck) eq 'y' then write1d = 1 else print,'No 1d spectrum written to file or results printed to mospec.'+logmask
        endif else write1d = 1
        
        if write1d then begin
           if not file_test('../onedspec') then spawn,'mkdir ../onedspec'
           h1d = h1
           sxaddpar,h1d,'MOSAP1',line[0]
           sxaddpar,h1d,'MOSAP2',line[1]
           if boxcar then sxaddpar,h1,'MOSPROF','boxcar'
           if gaussian then sxaddpar,h1,'MOSPROF','gaussian'
           if empirical then sxaddpar,h1,'MOSPROF','empirical'        
           if keyword_set(slitcor_2d) then sxaddpar,h1d,'SCOR_2D',slitcor_2d,' 2D slitcor'
           if keyword_set(slitcor_sed) then sxaddpar,h1d,'SCOR_SED',slitcor_sed,' SED slitcor'
           if keyword_set(twodsigma_arcsec) then sxaddpar,h1d,'LINESIZE',twodsigma_arcsec,' sigma in arcsec'
           if keyword_set(redshift) then sxaddpar,h1d,'REDSHIFT',redshift           
           sxaddpar,h1d,'WAT0_001',' system=equispec'
           sxaddpar,h1d,'WAT1_001',' wtype=linear label=Wavelength units=angstroms'
           sxaddpar,h1d,'WAT2_001',' wtype=linear'
           sxaddpar,h1d,'BUNIT','ergs/cm2/s/A'
           sxaddpar,h1d,'DISPAXIS','1'
           sxaddpar,h1d,'CTYPE1','LINEAR'
           sxaddpar,h1d,'CTYPE2','LINEAR',after='CTYPE1'
           sxaddpar,h1d,'CD2_2','0.0'
           sxaddpar,h1d,'DC-FLAG','0'
           cdelt = sxpar(h1,'CD1_1')
           sxaddpar,h1d,'CDELT1',cdelt
           sxaddpar,h1d,'APNUM1','1 0 object spectrum'
           sxaddpar,h1d,'APNUM2','2 1 error spectrum'
           sxaddpar,h1d,'APNUM3','3 2 abs. dev. of spectrum from fit (empty)'
           totspec = [[spec],[errspec],[replicate(0.0,n_elements(spec))]]
           writefits,'../onedspec/'+outfile1d,totspec,h1d
           print,'1D spectrum written as '+outfile1d
        endif
        key = ''
     endif
        
     ;;;;;
     if key eq 'f' then begin
        oldredshift = redshift
        cursor,fx,fy,/data,/nowait
        mindiff = min(abs(lambda-fix(fx)),startindex)
        startline = ''
        read,startline,prompt='Rest wavelength of line: '
        zeroindex = where(abs(spec) gt 0 and errspec gt 0 and finite(errspec) and finite(spec))

        bcor = 0
        if not modelflag then begin
           line_fit_mc,obslambda[zeroindex],spec[zeroindex],errspec[zeroindex],obslambda[startindex],spec[startindex],startline,result=result,perror=perror,covar=covar,redxi2=redxi2,zflag=zflag,vflag=vflag,fixredshift=fixredshift,fixvelocity=fixvelocity,cpoly=cpoly
           redshift = result[cpoly+1]
           cfit = poly(obslambda,result[0:cpoly])
        endif

        slitcor_sed = 0.0
        if modelflag then begin
           if cpoly gt 0 then cpolyorig = cpoly
           line_fit_mc,obslambda[zeroindex],spec[zeroindex],errspec[zeroindex],obslambda[startindex],spec[startindex],startline,result=result,perror=perror,covar=covar,zflag=zflag,vflag=vflag,fixredshift=fixredshift,fixvelocity=fixvelocity,contmed=contmed,/silent
           redshift = result[cpoly+1]
           zshift = (1.+redshift)/(1.+modelz)

           modellambda *= zshift
           modelspec_spl = spl_init(modellambda,modelspec,/double)
           modelspec_interp = spl_interp(modellambda,modelspec,modelspec_spl,obslambda,/double)
           modelspec_nobalmer = modelspec_interp
           nearbalmer = intarr(n_elements(obslambda))
           for n=0,n_elements(hlines)-1 do begin
              if hlines[n] lt 4000 then index = where(abs(obslambda/(1.+redshift)-hlines[n]) lt 25.0)
              if hlines[n] ge 4000 then index = where(abs(obslambda/(1.+redshift)-hlines[n]) lt 50.0)
              if index[0] ne -1 then begin
                 nearbalmer[index] = 1.0
                 if index[0] eq 0 or obslambda(index[-1]) eq max(obslambda) then begin
                    slope = 0
                    if index[0] eq 0 then modelspec_nobalmer[index] = slope*(obslambda[index]-obslambda[index[-1]+1])+modelspec_interp[index[-1]+1] else modelspec_nobalmer[index] = slope*(obslambda[index]-obslambda[index[0]-1])+modelspec_interp[index[0]-1]
                    
                 endif else begin
                    slope = (modelspec_interp[index[-1]+1]-modelspec_interp[index[0]-1])/(obslambda[index[-1]+1]-obslambda[index[0]-1])
                    modelspec_nobalmer[index] = slope*(obslambda[index]-obslambda[index[-1]+1])+modelspec_interp[index[-1]+1]
                 endelse
              endif
           endfor
           ;; index = where(~nearbalmer)
           ;; plot,modellambda,modelspec/slitcor_sed,psym=10,xrange=xrange,/xs
           ;; oplot,obslambda,modelspec_interp/slitcor_sed,psym=10,color=cgcolor('cyan')
           ;; oplot,obslambda,modelspec_nobalmer/slitcor_sed,psym=10,color=cgcolor('magenta')
           ;; stop

           modelmed = median(modelspec_nobalmer[zeroindex],/even)
           slitcor_sed = modelmed/contmed
           if slitcor_sed gt 0 then begin
              cpoly = -1
              cfit = modelspec_nobalmer/slitcor_sed          
              line_fit_mc,obslambda[zeroindex],spec[zeroindex],errspec[zeroindex],obslambda[startindex],spec[startindex],startline,result=result,perror=perror,covar=covar,redxi2=redxi2,zflag=zflag,vflag=vflag,fixredshift=fixredshift,fixvelocity=fixvelocity,sedcont=cfit[zeroindex],cpoly=cpoly,/silent
              ha_flag = 0 & hb_flag = 0
              mindiff = min(abs(result-6564),ha_index) & if(abs(result[ha_index]-6564) le 1.5) then ha_flag = 1
              mindiff = min(abs(result-4863),hb_index) & if(abs(result[hb_index]-4863) le 1.5) then hb_flag = 1
              ha_noabs = result[ha_index-1]*result[cpoly+2]*result[ha_index]*(1+result[cpoly+1])*sqrt(2*!DPI)*(1.0e-17)
              hb_noabs = result[hb_index-1]*result[cpoly+2]*result[hb_index]*(1+result[cpoly+1])*sqrt(2*!DPI)*(1.0e-17)

              cfit = modelspec_interp/slitcor_sed          
              line_fit_mc,obslambda[zeroindex],spec[zeroindex],errspec[zeroindex],obslambda[startindex],spec[startindex],startline,result=result,perror=perror,covar=covar,redxi2=redxi2,zflag=zflag,vflag=vflag,fixredshift=fixredshift,fixvelocity=fixvelocity,sedcont=cfit[zeroindex],cpoly=cpoly
              redshift = result[cpoly+1]
              ha_sed = result[ha_index-1]*result[cpoly+2]*result[ha_index]*(1+result[cpoly+1])*sqrt(2*!DPI)*(1.0e-17)
              hb_sed = result[hb_index-1]*result[cpoly+2]*result[hb_index]*(1+result[cpoly+1])*sqrt(2*!DPI)*(1.0e-17)
              if ha_flag then bcor = ha_sed/ha_noabs
              if hb_flag then bcor = hb_sed/hb_noabs

              modellambda /= (1.+redshift)
              modelspec_cor = modelspec/slitcor_sed
           endif else if modelflag and slitcor_sed le 0 then begin
              modelflag = 0
              cpoly = cpolyorig
              line_fit_mc,obslambda[zeroindex],spec[zeroindex],errspec[zeroindex],obslambda[startindex],spec[startindex],startline,result=result,perror=perror,covar=covar,redxi2=redxi2,zflag=zflag,vflag=vflag,fixredshift=fixredshift,fixvelocity=fixvelocity,cpoly=cpoly
              redshift = result[cpoly+1]
              cfit = poly(obslambda,result[0:cpoly])
           endif
        endif
        yfit = cfit
        for n=0,(n_elements(result)-(cpoly+3))/2.-1 do yfit += result[2*n+cpoly+3]*exp(-(obslambda-(1+result[cpoly+1])*result[2*n+cpoly+4])^2/(2*(result[cpoly+2]*(1+result[cpoly+1])*result[2*n+cpoly+4])^2))
        if keyword_set(cpolyorig) then cpoly = cpolyorig

        slitcor_2d = 0.0
        twodsigma_arcsec = 0.0
        twodfit = 0
        if not zflag and not vflag then begin
           mindiff = min(abs(obslambda-(startline*(1+redshift))),index)
           tmpimage = im1[index-25:index+25,*]
           tmperr = im2[index-25:index+25,*]
           tmpdim = size(tmpimage,/dim)
           twodindex = where(~finite(tmpimage) or ~finite(tmperr))
           tmpimage[twodindex] = 0.0
           tmperr[twodindex] = 0.0
           estimates = [0,max(tmpimage),5,2.5,25,mean(line),0]
           zfit = mpfit2dpeak(tmpimage,a,estimates=estimates,error=tmperr,perror=a_snr,/tilt)
           a_snr = a/a_snr
           if a_snr[2] ge 10. and a_snr[3] ge 10. then begin
              twodfit = 1
              tilt = abs(1./cos(a[6]))
              twodsigma_arcsec = 0.18*a[3]*tilt
              aperture = fix(line[1]-line[0]+1)
              slitcor_2d = 1./(erf(0.247487/twodsigma_arcsec)*erf(aperture*0.18/(2.*sqrt(2.)*twodsigma_arcsec)))
              window,2,xsize=3*tmpdim[0],ysize=3*tmpdim[1],title='2D fit'
              plot,[0],[0],/noerase,position=[0,0,3*tmpdim],/device,xrange=[0,tmpdim[0]],yrange=[0,tmpdim[1]],/xs,/ys
              tvimage,bytscl(tmpimage,scale[0],scale[1])
              plotsym,0,1,/fill
              oplot,[a[4]],[a[5]],psym=8,color=cgcolor('red')
              if a[1] ge 0 then contour,zfit,/overplot,c_thick=2,c_colors=cgcolor(['blue','green']),levels=a[1]*[0.1,0.5]
              sxaddpar,h1,'LINEXPOS',(index-25)+a[4]
              sxaddpar,h1,'LINEYPOS',a[5]
              if not keyword_set(serendip) or serendip eq '' then writefits,file,im1,h1
              wset,4
           endif
        endif
        if not twodfit then print,'Tilt: 0.00 degrees CCW   SED Slitcor: '+strtrim(string(slitcor_sed,format='(F0.2)'),1)+'   2D Slitcor: '+strtrim(string(slitcor_2d,format='(F0.2)'),1)+'   Balmer cor: '+strtrim(string(bcor,format='(F0.2)'),1) $
        else if twodfit then print,'Tilt: '+strtrim(string(a[6]*180/!DPI,format='(F0.2)'),1)+' degrees CCW   SED Slitcor: '+strtrim(string(slitcor_sed,format='(F0.2)'),1)+'   2D Slitcor: '+strtrim(string(slitcor_2d,format='(F0.2)'),1)+'   Balmer cor: '+strtrim(string(bcor,format='(F0.2)'),1)

        dev_spec = spec-yfit
        abs_dev = abs(dev_spec)
        mad_val = median(abs_dev[where(abs(spec) gt 0)],/even)
        mask_scale = abs(errspec/spec)
        mask_scale[where(mask_scale ge 1.0)] = 1.0
        mask_index = where(abs_dev*mask_scale lt 2.*mad_val or (obslambda ge 22500 and abs(spec) gt 0))
        mask_spec = spec & mask_spec[mask_index] = 0.0/0.0

        lambda = obslambda/(redshift+1.0)
        xrange *= (oldredshift+1.0)/(redshift+1.0)
        fitcheck = 1
        lineid = 1
        key = 'r'
     endif

     ;;;;;
     if key eq 'g' then begin
        if not keyword_set(startline) then print,"Please compute fit to spectrum using 'f' first"
        if keyword_set(startline) then begin
           write1d = 0
           if file_test('../onedspec/'+outfile1d) then begin
              printcheck = ''
              read,printcheck,prompt="There is already a 1d extracted spectrum for this object. Continue [y/n]? "
              if strlowcase(printcheck) eq 'y' then write1d = 1 else print,'No 1d spectrum written to file or results printed to mospec.'+logmask
           endif else write1d = 1
           
           if write1d then begin
              sxaddpar,h1,'MOSAP1',line[0]
              sxaddpar,h1,'MOSAP2',line[1]
              if boxcar then sxaddpar,h1,'MOSPROF','boxcar'
              if gaussian then sxaddpar,h1,'MOSPROF','gaussian'
              if empirical then sxaddpar,h1,'MOSPROF','empirical'
              sxaddpar,h1,'REDSHIFT',redshift           
              
              if filter eq 'Y' then slitcor_fitlers = ['Y']
              if filter eq 'J' then slitcor_filters = ['J','J2','J3']
              if filter eq 'H' then slitcor_filters = ['H1','H2']
              if filter eq 'K' then slitcor_filters = ['Ks']
              
              for n=0,n_elements(slitcor_filters)-1 do begin
                 readcol,calibdir+'mosfire_'+slitcor_filters[n]+'.txt',filtlam,filtthru,format='D,D',/silent
                 filtthru = filtthru[sort(filtlam)]
                 filtlam = filtlam[sort(filtlam)]*10000
                 filtthru_spl = spl_init(filtlam,filtthru,/double)
                 index = where(obslambda ge min(filtlam) and obslambda le max(filtlam) and abs(spec) gt 0 and ~finite(mask_spec))
                 filtthru_interp = spl_interp(filtlam,filtthru,filtthru_spl,obslambda[index],/double)
                 filt_wgts = filtthru_interp/total(filtthru_interp)
                 
                 galerr_fnu = (errspec[index]*1.0e-17)*lambda[index]*lambda[index]/(2.99792458e18)
                 galspec_fnu = (spec[index]*1.0e-17)*lambda[index]*lambda[index]/(2.99792458e18)
                 gal_fnu = total(galspec_fnu*filt_wgts)
                 gal_fnu_var = total(filt_wgts^2.*galerr_fnu^2.)
                 mag_ab = -2.5*alog10(gal_fnu)-48.6
                 mag_ab_err = 2.5*(sqrt(gal_fnu_var)/gal_fnu)/alog(10.)
                 sxaddpar,h1,strupcase(' MAG_'+slitcor_filters[n]),mag_ab
                 sxaddpar,h1,strupcase(' EMAG_'+slitcor_filters[n]),mag_ab_err
              endfor
              
              if not file_test('../onedspec') then spawn,'mkdir ../onedspec'
              h1d = h1
              sxaddpar,h1d,'OBJECT',object
              sxaddpar,h1d,'SCOR_2D',slitcor_2d,' 2D slitcor'
              sxaddpar,h1d,'SCOR_SED',slitcor_sed,' SED slitcor'
              sxaddpar,h1d,'LINESIZE',twodsigma_arcsec
              sxaddpar,h1d,'WAT0_001',' system=equispec'
              sxaddpar,h1d,'WAT1_001',' wtype=linear label=Wavelength units=angstroms'
              sxaddpar,h1d,'WAT2_001',' wtype=linear'
              sxaddpar,h1d,'BUNIT','ergs/cm2/s/A'
              sxaddpar,h1d,'DISPAXIS','1'
              sxaddpar,h1d,'CTYPE1','LINEAR'
              sxaddpar,h1d,'CTYPE2','LINEAR',after='CTYPE1'
              sxaddpar,h1d,'CD2_2','0.0'
              sxaddpar,h1d,'DC-FLAG','0'
              cdelt = sxpar(h1,'CD1_1')
              sxaddpar,h1d,'CDELT1',cdelt
              sxaddpar,h1d,'APNUM1','1 0 object spectrum'
              sxaddpar,h1d,'APNUM2','2 1 error spectrum'
              sxaddpar,h1d,'APNUM3','3 2 deviation of spectrum from fit'
              totspec = [[spec],[errspec],[dev_spec]]
              
              close,1
              if not file_test('../logs') then spawn,'mkdir ../logs'
              logfile = '../logs/mospec.'+logmask
              openw,1,logfile,/append
              printf,1,systime()
              printf,1,outfile1d,' ',fix(line)
;;               if not modelflag or (modelflag and slitcor_sed le 0) then line_fit_mc,obslambda[zeroindex],spec[zeroindex],errspec[zeroindex],obslambda[startindex],spec[startindex],startline,result=result,perror=perror,covar=covar,redxi2=redxi2,/print
              if modelflag and slitcor_sed gt 0 then line_fit_mc,obslambda[zeroindex],spec[zeroindex],errspec[zeroindex],obslambda[startindex],spec[startindex],startline,result=result,perror=perror,covar=covar,redxi2=redxi2,sedcont=cfit[zeroindex],cpoly=-1,/print else line_fit_mc,obslambda[zeroindex],spec[zeroindex],errspec[zeroindex],obslambda[startindex],spec[startindex],startline,result=result,perror=perror,covar=covar,redxi2=redxi2,cpoly=cpoly,/print
              if not twodfit then printf,1,'Tilt: 0.00 degrees CCW   SED Slitcor: '+strtrim(string(slitcor_sed,format='(F0.2)'),1)+'   2D Slitcor: '+strtrim(string(slitcor_2d,format='(F0.2)'),1)+'   Balmer cor: '+strtrim(string(bcor,format='(F0.2)'),1) $
              else if twodfit then printf,1,'Tilt: '+strtrim(string(a[6]*180/!DPI,format='(F0.2)'),1)+' degrees CCW   SED Slitcor: '+strtrim(string(slitcor_sed,format='(F0.2)'),1)+'   2D Slitcor: '+strtrim(string(slitcor_2d,format='(F0.2)'),1)+'   Balmer cor: '+strtrim(string(bcor,format='(F0.2)'),1)

              printf,1,''
              close,1
              print,'Results of line fitting printed to mospec.'+logmask
              
              writefits,'../onedspec/'+outfile1d,totspec,h1d
              if not keyword_set(serendip) or serendip eq '' then writefits,file,im1,h1
              print,'1d spectrum written as '+outfile1d           
              
           endif             
        endif
        key = ''
     endif

     ;;;;; 
     if key eq 'h' then begin
        if not fitcheck then fitcheck = 1 else fitcheck = 0
        key = 'r'
     endif

     ;;;;;
     if key eq 'i' then begin
        oldlineid = lineid
        if not lineid and redshift gt 0 then lineid = 1 else lineid = 0
        if redshift eq 0 then print,'Please provide a redshift before identifying lines'
        if lineid ne oldlineid then key = 'r'
     endif

     ;;;;;
     if key eq 'j' then begin
        cursor,jx,jy,/data,/nowait
        mindiff = min(abs(lambda-fix(jx)),index)
        a = [spec[index],obslambda[index],0]
        gauss = gaussfit(obslambda[index-25:index+25],spec[index-25:index+25],a,nterms=3,measure_errors=errspec[index-25:index+25])
        oplot,obslambda[index-25:index+25],gauss,color=cgcolor('sky blue')
        if not keyword_set(startline) then startline = ''
        read,startline,prompt='Rest wavelength of line: '
        mindiff = min(abs(linewav-startline),index) 
        wave = linewav[index]
        oldredshift = redshift
        redshift = a[1]/wave-1.0
        lambda = obslambda/(redshift+1.0)
        xrange = xrange*(oldredshift+1.0)/(redshift+1.0)
        lineid = 1
        key = 'r'
     endif
     
     ;;;;; 
     if key eq 'k' then begin
        if sxpar(h1,'LINEXPOS') gt 0 and sxpar(h1,'LINEYPOS') gt 0 then begin
           sxdelpar,h1,'LINEXPOS'
           sxdelpar,h1,'LINEYPOS'
           writefits,file,im1,h1
           print,'Emission line position removed from spectrum header'
        endif else print,'No emission line position recorded'
        key = ''
     endif

     ;;;;; 
     if key eq 'l' then begin
        line = ''
        read,line,prompt='Line range: '
        line = fix(strsplit(line,' ',count=lcount,/extract))
        while max(line) gt dim[1] do read,line,prompt='Out of bounds. New line range: '
        if lcount eq 1 then read,line,prompt='Please specify a beginning and ending row. New line range: '
        if lcount eq 2 then begin
           test1 = total(im1[*,line[0]:line[1]],2)*cal
           test2 = im2*im2
           test2 = sqrt(total(test2[*,line[0]:line[1]],2))*cal
        endif else message,'Invalid line range'
        test1 = gauss_smooth(test1,width=fix(sval),/nan)
        yrange = [0,1.1*max(test1[where(finite(test1))])] & yrange[0] = -1*yrange[1]/3.
        spec = test1
        objspec = test1
        errspec = test2
        if sval gt 1 then objspec[smoothskyindex] = 0.0/0 else objspec[skyindex] = 0.0/0
        if not keyword_set(serendip) or serendip eq '' then begin
           sxaddpar,h1,'MOSAP1',line[0]
           sxaddpar,h1,'MOSAP2',line[1]
           writefits,file,im1,h1
           print,'New aperture written to header'           
        endif
        key = 'r'
     endif

     ;;;;;
     if key eq 'm' then begin
        nmask = sxpar(h1,'NCOMBINE')
        if nmask gt 1 then begin
           for n=0,nmask-1 do begin
              if n eq 0 then masks = strtrim(sxpar(h1,'MASK'+strtrim(string(n+1),1)))
              if n gt 0 then masks += ', '+strtrim(sxpar(h1,'MASK'+strtrim(string(n+1),1)))
           endfor
           print,'2D stack made from: ',masks
           print,'Total exposure time (in seconds): ',sxpar(h1,'MOSITIME')
        endif else begin
           masks = sxpar(h1,'MASKNAME')
           print,'Shifted spectrum from: ',masks
           print,'Total exposure time (in seconds): ',sxpar(h1,'MOSITIME')
        endelse
     endif

     ;;;;;
     if key eq 'o' then begin
        cursor,ox,oy,/data,/nowait
        mindiff = min(abs(lambda-(ox-10)),lindex)
        mindiff = min(abs(lambda-(ox+10)),uindex)
        tmpimage = im1[lindex:uindex,*]
        row = indgen(dim[1])
        flux = fltarr(n_elements(row))        
        for n=0,n_elements(tmpimage[0,*])-1 do begin
           resistant_mean,tmpimage[*,n],3.,fluxavg,/double
           flux[n] = fluxavg
        endfor

        index = where(finite(flux))
        row = row[index]
        flux = flux[index]
        expr = 'P[1]+P[2]*exp(-(X-P[3])^2/(2*P[0]^2))+P[4]*exp(-(X-P[5])^2/(2*P[0]^2))+P[6]*exp(-(X-P[7])^2/(2*P[0]^2))'
        guess = mean(line);row[where(flux eq max(flux))]
        start = [5,0.D,max(flux),guess,min(flux),guess-17,min(flux),guess+17]
        result = mpfitexpr(expr,row,flux,flux*0+1.,start,yfit=yfit,/quiet,status=status,errmsg=errmsg)
        if status le 0 then message,errmsg
        profile1 = result[1]+result[2]*exp(-(row-result[3])^2/(2*result[0]^2)) & profile1[where(yfit le result[1])] = 0.0
        profile2 = flux & profile2[where(yfit le result[1])] = 0.0 & profile2[where(profile2 lt 0)] = 0.0
        plot,row,flux,title=object+'-'+file,xtitle='Row',ytitle='Average flux per pixel',yticks=3,ytickname=replicate(' ',4),psym=10,yrange=[-0.3*max(flux),1.5*max([profile1,profile2])],/ys
        oplot,!X.CRANGE,[0,0],linestyle=0
        oplot,row[where(profile1 gt 0)],profile1[where(profile1 gt 0)],thick=2,psym=10,color=cgcolor('green')
        oplot,row[where(profile2 gt 0)],profile2[where(profile2 gt 0)],thick=2,psym=10,color=cgcolor('magenta')

        boxcar = 0
        gaussian = 0
        empirical = 1
        xyouts,0.11,0.85,"Magenta = empirical profile, green = best-fit Gaussian profile",/normal,color=cgcolor('white')

        escape = 0
        while ~escape do begin
           if gaussian then begin
              xyouts,0.11,0.80,"Accept empirical optimal extraction weighting function? [y/n]",/normal,color=cgcolor('black')
              xyouts,0.11,0.75,"Press 'o' again to choose Gaussian profile",/normal,color=cgcolor('black')
              xyouts,0.11,0.80,"Accept Gaussian optimal extraction weighting function? [y/n]",/normal,color=cgcolor('white')
              xyouts,0.11,0.75,"Press 'o' again to choose empirical profile",/normal,color=cgcolor('white')
           endif else begin
              xyouts,0.11,0.80,"Accept Gaussian optimal extraction weighting function? [y/n]",/normal,color=cgcolor('black')
              xyouts,0.11,0.75,"Press 'o' again to choose empirical profile",/normal,color=cgcolor('black')
              xyouts,0.11,0.80,"Accept empirical optimal extraction weighting function? [y/n]",/normal,color=cgcolor('white')
              xyouts,0.11,0.75,"Press 'o' again to choose Gaussian profile",/normal,color=cgcolor('white')
           endelse
           
           key = get_kbrd()
           
           if key eq 'o' then begin
              gaussian = 1-gaussian
              empirical = 1-empirical
           endif
            
           if key eq 'y' then begin
              if gaussian then profile = replicate(1,dim[0]) # profile1/total(profile1)
              if empirical then profile = replicate(1,dim[0]) # profile2/total(profile2)
              tmpline = minmax(row[where(profile[0,*] gt 0)])
              
              ivar = 1./(im2*im2)
              ivar[where(~finite(ivar))] = 0.0
              test1 = total(im1*profile*ivar,2)/total(profile*profile*ivar,2)*cal
              test2 = sqrt(total(profile,2)/total(profile*profile*ivar,2))*cal
              
              spec = test1
              objspec = test1
              errspec = test2
              if sval gt 1 then objspec[smoothskyindex] = 0.0/0 else objspec[skyindex] = 0.0/0
              if gaussian then print,'1D spectrum determined using optimal extraction with a Gaussian profile'
              if empirical then print,'1D spectrum determined using optimal extraction with an empirical profile'
              line = tmpline
              if not keyword_set(serendip) or serendip eq '' then begin
                 sxaddpar,h1,'MOSAP1',line[0]
                 sxaddpar,h1,'MOSAP2',line[1]
                 writefits,file,im1,h1
                 print,'New aperture written to header'           
              endif
              escape = 1
           endif
           if key eq 'n' then begin
              print,'No new 1D spectrum extracted'
              boxcar = 1
              gaussian = 0
              empirical = 0
              escape = 1
           endif
        endwhile
     
        key = 'r'
     endif

     ;;;;;
     if key eq 'p' then begin
        tmpxrange = ''
        read,tmpxrange,prompt='Wavelength range to plot (Angstrom): '
        if min(tmpxrange) lt min(lambda) or max(tmpxrange) gt max(lambda) then read,tmpxrange,prompt='Invalid wavelength range. New range: '
        tmpyrange = 1.0
        read,tmpyrange,prompt='Maxmimum flux to plot (1e-17 erg/s/cm^2/Ang): '
        tmpyrange = [-1*tmpyrange[0]/3.,1.1*tmpyrange[0]]
        tmpxrange = fix(strsplit(tmpxrange[sort(tmpxrange)],' ',/extract))
        mindiff = min(abs(lambda-tmpxrange[0]),tmpindex1)
        mindiff = min(abs(lambda-tmpxrange[1]),tmpindex2) 
        tmpxrange = [lambda[tmpindex1],lambda[tmpindex2]]

        xsize = 8.8
        xmargin = 1.5
        ymargin = 1.1
        wall = 0.3
        xplotsize = xsize-(xmargin+wall)
        yplotsize = xplotsize/10.
        ysize = ymargin+yplotsize+(4*yplotsize)+wall

        ypixels = yplotsize*(tmpindex2-tmpindex1)/xplotsize
        if ypixels lt dim[1] and fix(mean(line)+ypixels/2.) lt dim[1] and fix(mean(line)-ypixels/2.) ge 0 then tmptwodspec = im1[tmpindex1:tmpindex2,fix(mean(line)-ypixels/2.):fix(mean(line)+ypixels/2.)]
        if ypixels ge dim[1] or fix(mean(line)+ypixels/2.) gt dim[1] or fix(mean(line)-ypixels/2.) lt 0 then tmptwodspec = im1[tmpindex1:tmpindex2,*]
        
        set_plot,'PS'
        !p.charsize = 1
        !p.font = 1
        !p.charthick = 4
        !p.thick = 4
        !x.thick = 4
        !y.thick = 4
        device,filename=object+'_'+filter+'.eps',xsize=xsize,ysize=ysize,/color,/encapsulated,font_size=12,/tt_font,decomposed=0
        angstrom = '!3' + STRING(197B) + '!X'
        if redshift eq 0 then xtitle='Observed wavelength ('+angstrom+')'
        if redshift gt 0 then xtitle='Rest wavelength ('+angstrom+')'

        plot,[0],[0],xrange=tmpxrange,yrange=rowrange,position=[xmargin/xsize,ymargin/ysize,(xsize-wall)/xsize,(ymargin+yplotsize)/ysize],/normal,xstyle=1,xticks=1,xtickname=replicate(' ',2),ystyle=1,yticks=1,ytickname=replicate(' ',2),color=cgcolor('black')
        tvimage,bytscl(tmptwodspec,scale[0],scale[1]),/overplot
        plot,[0],[0],xrange=tmpxrange,yrange=rowrange,position=[xmargin/xsize,ymargin/ysize,(xsize-wall)/xsize,(ymargin+yplotsize)/ysize],xtickformat='(I)',/normal,xstyle=1,ystyle=1,yticks=1,ytickname=replicate(' ',2),xtitle=xtitle,color=cgcolor('black'),/noerase
        plot,[min(tmpxrange),max(tmpxrange)],[0,0],xrange=tmpxrange,yrange=tmpyrange,xstyle=1,ystyle=1,xticks=1,xtickname=replicate(' ',2),position=[xmargin/xsize,(ymargin+yplotsize)/ysize,(xsize-wall)/xsize,(ysize-wall)/ysize],/normal,/noerase,psym=10,color=cgcolor('black'),ytitle=textoidl("F_{\lambda} (10^{-17} erg/s/cm^2/"+angstrom+")"),xminor=4,yminor=4
        oplot,lambda,spec,psym=10,color=cgcolor('dark grey'),thick=1
        oplot,lambda,objspec,psym=10,color=cgcolor('black'),thick=3
        if not mask then oplot,lambda,errspec+!Y.CRANGE[0],psym=10,color=cgcolor('red'),thick=1
        
        ;xyouts,0.09*(!X.CRANGE[1]-!X.CRANGE[0])+!X.CRANGE[0],0.8*(!Y.CRANGE[1]-!Y.CRANGE[0])+!Y.CRANGE[0],filter+'-band',/data,color=cgcolor('black')
        if filter eq 'K' then xyouts,0.91*(!X.CRANGE[1]-!X.CRANGE[0])+!X.CRANGE[0],0.8*(!Y.CRANGE[1]-!Y.CRANGE[0])+!Y.CRANGE[0],object,/data,color=cgcolor('black'),align=1
        if filter eq 'K' and redshift gt 0 then xyouts,0.91*(!X.CRANGE[1]-!X.CRANGE[0])+!X.CRANGE[0],0.65*(!Y.CRANGE[1]-!Y.CRANGE[0])+!Y.CRANGE[0],'z = '+strtrim(string(redshift,format='(F0.4)'),1),/data,color=cgcolor('black'),align=1
        
        if fitcheck eq 1 then begin
           ;; oplot,lambda,yfit_poly,color=cgcolor('magenta')
           oplot,lambda,yfit,color=cgcolor('steel blue'),thick=2
        endif
        
        if redshift gt 0 and lineid then begin
           for n=0,n_elements(lines)-1 do begin
              if ((linewav[n] gt !X.CRANGE[0]) and (linewav[n] lt !X.CRANGE[1])) then begin
                 vline,linewav[n],linestyle=1,color=cgcolor('dark green')
                 printline = lines[n]
                 if printline eq 'H-alpha' then printline = textoidl('H\alpha')
                 if printline eq 'H-beta' then printline = textoidl('H\beta')
                 xyouts,linewav[n]+1,0.85*!Y.CRANGE[0],printline,/data,color=cgcolor('dark green'),orient=270,align=1,charsize=0.75
              endif
           endfor
        endif
        device,/close
        print,'Postscript figure written to '+object+'_'+filter+'.eps'
        spawn,'open '+object+'_'+filter+'.eps'
        set_plot,'X'
        !p.charsize = 2
        !p.thick = 1
        !x.thick = 2
        !y.thick = 2
        !p.color=cgcolor('white')
        !p.font = -1
        !p.charthick = 1
        key = ''
     endif
     
     ;;;;; 
     if key eq 's' then begin
        xyouts,0.11,0.85,'Choose a smoothing value between 1 and 9',/normal,color=cgcolor('white')
        xyouts,0.11,0.80,'A value of 1 returns the original spectrum',/normal,color=cgcolor('white')
        sval = fix(1)
        sval = get_kbrd()
        test1 = total(im1[*,line[0]:line[1]],2)*cal
        test2 = im2*im2
        test2 = sqrt(total(test2[*,line[0]:line[1]],2))*cal
        test1 = gauss_smooth(test1,width=fix(sval),/nan)
        test2 = gauss_smooth(test2,wdith=fix(sval),/nan)
        spec = test1
        objspec = test1
        if sval gt 1 then objspec[smoothskyindex] = 0.0/0 else objspec[skyindex] = 0.0/0
        key = 'r'
     endif
     
     ;;;;; 
     if key eq 't' then begin
        cursor,tx,ty,/data,/nowait
        mindiff = min(abs(lambda-(tx-10)),lindex)
        mindiff = min(abs(lambda-(tx+10)),uindex)
        tmpimage = im1[lindex:uindex,*]
        row = indgen(dim[1])
        flux = fltarr(n_elements(row))        
        for n=0,n_elements(tmpimage[0,*])-1 do begin
           resistant_mean,tmpimage[*,n],3.,fluxavg,/double
           flux[n] = fluxavg
        endfor
        index = where(finite(flux))
        plot,row[index],flux[index],title=object+'-'+file,xtitle='Row',ytitle='Average flux per pixel',yticks=3,ytickname=replicate(' ',4),psym=10,yrange=[-0.3*max(flux),1.1*max(flux)]
        oplot,!X.CRANGE,[0,0],linestyle=0

        print,abs(sxpar(h1,'CRVAL2')),format='(%"Expected location of object: %i")'
        print,row[where(flux eq max(flux))],format='(%"Best guess for object location: %i")'
        escape = 0
        while not escape do begin
           vline,line,linestyle=2
           xyouts,0.11,0.85,"Press spacebar twice to mark limits of extraction aperture",/normal,color=cgcolor('white')
           xyouts,0.11,0.80,"Press 'q' to return without re-extracting 1D spectrum",/normal,color=cgcolor('white')
           key = get_kbrd()
           if key ne 'q' then begin
              if key eq ' ' then cursor,tx1,ty1,/data,/nowait
              vline,tx1,linestyle=1,color=cgcolor('white')
              key = get_kbrd()
              if key eq ' ' then cursor,tx2,ty2,/data,/nowait
              vline,tx2,linestyle=1,color=cgcolor('white')
              xyouts,0.11,0.85,"Press spacebar twice to mark limits of extraction aperture",/normal,color=cgcolor('black')
              xyouts,0.11,0.80,"Press 'q' to return without selecting aperture",/normal,color=cgcolor('black')
              xyouts,0.11,0.85,"Accept choice? [y/n]",/normal,color=cgcolor('white')
              xyouts,0.11,0.80,"Press 'q' to return without selecting aperture",/normal,color=cgcolor('white')
              key = get_kbrd()
           endif else escape = 1
           if key ne 'y' and key ne 'n' then escape = 1
           if key eq 'n' then begin
              vline,tx1,linestyle=1,color=cgcolor('black')
              vline,tx2,linestyle=1,color=cgcolor('black')
              xyouts,0.11,0.85,"Accept choice? [y/n]",/normal,color=cgcolor('black')
           endif
           if key eq 'y' then begin
              line = [round(tx1),round(tx2)] & line = minmax(line)
              print,line[0],line[1],format='(%"Selected line region: %i %i")'
              test1 = total(im1[*,line[0]:line[1]],2)*cal
              test2 = im2*im2
              test2 = sqrt(total(test2[*,line[0]:line[1]],2))*cal
              test1 = gauss_smooth(test1,width=fix(sval),/nan)
              yrange = [0,1.1*max(test1[where(finite(test1))])] & yrange[0] = -1*yrange[1]/3.
              spec = test1
              objspec = test1
              errspec = test2
              if sval gt 1 then objspec[smoothskyindex] = 0.0/0 else objspec[skyindex] = 0.0/0
              sxaddpar,h1,'MOSAP1',line[0]
              sxaddpar,h1,'MOSAP2',line[1]
              if not keyword_set(serendip) or serendip eq '' then begin
                 writefits,file,im1,h1
                 print,'New aperture written to header'
              endif
              escape = 1
           endif
        endwhile
        key = 'r'
     endif
     
     ;;;;; 
     if key eq 'u' then begin
        if not mask then mask = 1 else mask = 0
        key = 'r'
     endif
     
     ;;;;; 
     if key eq 'v' then begin
        if not keyword_set(serendip) or serendip eq '' then begin
           plot,[0],[0],xrange=xrange,yrange=rowrange,position=[115,65,115+dim[0]*zoom,65+dim[1]*zoom],xtickformat='(I)',/device,xstyle=1,ystyle=1,yticks=1,ytickname=replicate(' ',2),color=cgcolor('white'),/noerase
           cursor,x,y,/data,/nowait
           x = (x*(1.+redshift)-lam0)/delta+(refpix-1)
           print,'Write following coordinates to header [y/n]? ',x,y
           key = get_kbrd()
           if key eq 'y' then begin
              sxaddpar,h1,'LINEXPOS',x
              sxaddpar,h1,'LINEYPOS',y
              writefits,file,im1,h1
              print,'Selected position written to header'
           endif else if key eq 'n' then print,'No position written to header'
        endif else print,'Cannot record emission line position for serendipitous sources'
        key = ''
     endif

     ;;;;;
     if key eq 'w' then begin
        changed = 0
        if not keyword_set(storez) then storez = 0.0
        if redshift eq 0.0 and storez eq 0.0 then print,'Please identify a redshift first'
        if redshift gt 0.0 and not changed then begin
           storez = redshift
           redshift = 0.0
           lambda = obslambda
           xrange = xrange*(storez+1.0)
           changed = 1
        endif
        if storez gt 0.0 and not changed then begin
           redshift = storez
           lambda = obslambda/(redshift+1.0)
           xrange = xrange/(redshift+1.0)
           changed = 1
        endif
        key = 'r'
     endif

     ;;;;; 
     if key eq 'z' then begin
        oldredshift = redshift
        read,redshift,prompt='Redshift: '
        lambda = obslambda/(redshift+1.0)
        xrange = xrange*(oldredshift+1.0)/(redshift+1.0)
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
        yrange[0] = -1*yrange[1]/3.
        key = 'r'
     endif

     ;;;;;
     if key eq '/' then begin
        if not keyword_set(yfit) then print,"Please fit spectrum with 'f' first" else begin
           plotsym,0,1.0,/fill        
           oplot,lambda,mask_spec,psym=8,color=cgcolor('red')
           xyouts,0.11,0.85,"Red points will be masked in stack_1d",/normal,color=cgcolor('white')
           xyouts,0.11,0.80,"Press any key to continue",/normal,color=cgcolor('white')
           key = get_kbrd()
           key = 'r'
        endelse
     endif
     
     ;;;;; 
     if key eq 'r' then begin
        erase
        if redshift eq 0 then xtitle='Observed wavelength (!6!sA!r!u!9 %!6!n)'
        if redshift gt 0 then xtitle='Rest wavelength (!6!sA!r!u!9 %!6!n)'

        plot,[0],[0],xrange=xrange,yrange=rowrange,position=[115,65,115+dim[0]*zoom,65+dim[1]*zoom],/device,xstyle=1,xticks=1,xtickname=replicate(' ',2),ystyle=1,yticks=1,ytickname=replicate(' ',2),color=cgcolor('white')
        tvimage,bytscl(twodspec,scale[0],scale[1]),/overplot
        plot,[0],[0],xrange=xrange,yrange=rowrange,position=[115,65,115+dim[0]*zoom,65+dim[1]*zoom],xtickformat='(I)',/device,xstyle=1,ystyle=1,yticks=1,ytickname=replicate(' ',2),xtitle=xtitle,color=cgcolor('white'),/noerase
        oplot,[0,xrange[1]],[line[0],line[0]],color=cgcolor('green')
        oplot,[0,xrange[1]],[line[1],line[1]],color=cgcolor('green')
        if lineid and redshift gt 0 then for n=0,n_elements(lines)-1 do vline,linewav[n],linestyle=4,color=cgcolor('green')
        plot,[min(xrange),max(xrange)],[0,0],xrange=xrange,yrange=yrange,xstyle=1,title=object+' - '+file,ytitle=textoidl("F_{\lambda} (10^{-17} erg/s/cm^2/!6!sA!r!u!9 %!6!n)"),xticks=1,xtickname=replicate(' ',2),position=[115,65+dim[1]*zoom,115+dim[0]*zoom,450+dim[1]*zoom],/device,/noerase,psym=10,color=cgcolor('white')
        if mask then oplot,lambda,spec,psym=10,color=cgcolor('dark grey') else oplot,lambda,spec,psym=10,color=cgcolor('white')
        oplot,lambda,objspec,psym=10,color=cgcolor('white')
        if not mask then oplot,lambda,errspec,psym=10,color=cgcolor('red')
        if file_test('../onedspec/'+outfile1d) then xyouts,115+dim[0]*zoom,460+dim[1]*zoom,'1D ',/align,/device,color=cgcolor('orange')
        
        if redshift gt 0 then xyouts,0.95,0.85,'z = '+strtrim(string(redshift,format='(F0.4)'),1),/normal,alignment=1.0,color=cgcolor('white')

        if lineid and redshift gt 0 then begin
           for n=0,n_elements(lines)-1 do begin
              vline,linewav[n],linestyle=4,color=cgcolor('green')
              if (linewav[n] gt !X.CRANGE[0] and linewav[n] lt !X.CRANGE[1]) then $
                 xyouts,linewav[n]+1,yrange[0],lines[n],/data,color=cgcolor('green'),orient=270,align=1
           endfor
        endif

        !p.thick = 2
        if fitcheck eq 1 then begin
           oplot,lambda,yfit,color=cgcolor('yellow')
;;            oplot,lambda,yfit_poly,color=cgcolor('orange')
;;            if min(lambda) lt 4363 then begin
;;               te_level = (0.1)/(result[cpoly+2]*sqrt(2*!DPI)) ;in units of 1e-17
;;               oplot,lambda,cfit+te_level,color=cgcolor('white')
;;            endif
        endif
        !p.thick = 1
        
        if modelflag and fitcheck then begin
           oplot,modellambda,modelspec,psym=10,color=cgcolor('cyan')
           if keyword_set(modelspec_cor) then oplot,modellambda,modelspec_cor,psym=10,color=cgcolor('magenta')
        endif

        !p.color = cgcolor('white')
        key = ''
     endif
     
     if key eq 'x' then stop
          
     ;;;;; 
     if key eq '?' then begin
        print,""
        print,"--- AVAILABLE COMMANDS ---"
        print,"spacebar - print current location of cursor to command line"
        print,"a - zoom: press once in lower left corner of desired region and once again in the upper right corner; press twice in the same spot to unzoom"
        print,"b - move to previous object (only works when /all or multiple spectra have been called)"
        print,"c - clear all changes to spectrum and redraw with original aperture guess"
        print,"d - determine seeing on mask"
        print,"e - extract and save current 1d spectrum without absolute deviation array"
        print,"f - calculate simultaneous fit to the continuum and strong emission lines"
        print,"g - write line fluxes and fit parameters to mospec.log"
        ;print,"h - "
        print,"i - identify emission lines using the default list or lines.dat (requires a redshift)"
        print,"j - mark emission line and provide its rest wavelength to return redshift"
        print,"k - remove emission line location from spectrum header"
        print,"l - choose a range of lines to define the aperture (e.g., '60 70')"
        print,"m - print the MASKS header keyword, if available"
        print,"n - move to next object (only works when /all or multiple spectra have been called)"
        print,"o - optimally extract spectrum using spectrum near cursor to define the weighting profile"
        print,"p - prints figure of current object to file"
        print,"q - quit"
        print,"r - redraw the spectrum with current zoom, smoothing, masking, and line list"
        print,"s - smooth spectrum (must choose a value from 1-9)"
        print,"t - view the trace selection window"
        print,"u - mask or unmask noisy regions"
        print,"v - mark cursor position as location of emission line center"
        print,"w - switch between rest and observed wavelength"
        print,"x - pauses the IDL routine; to restart program, type '.c'"
        ;print,"y - "
        print,"z - apply a redshift solution to the spectrum"
        print,", - zoom out by a factor 2 in flux"
        print,". - zoom in by a factor 2 in flux"
        print,"/ - show points that would be masked in stack_1d"
        print,"? - help"
        print,""
     endif
        
  endwhile
   if key eq 'b' then i -= 2
   if i lt -1 then i = -1
   if key eq 'q' then return
endfor
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; CHANGE LOG ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; 01/25/2013
; + changed 't' extraction to be 20 pixels along wavelength direction
; + changed redshift error to be reported in km/s
; + removed dependence on gaussian() in line_fit.pro
; + changed 'e' to extract 1d spectrum using both even and odd numbers
;   of pixels in the aperture 
; 
; 01/31/2013
; + startup 1d spectrum now has an extraction aperture of 7 pixels, 
;   centered on the best guess location 
; + aperture selection prints best guess for location and selected
;   line region
; + continuum fit excludes first and last 100 Angstroms of spectrum
; + flux and line width errors include term from covariance matrix
; + can run /all in a directory with both single and combined 2d
;   spectra and measure both at the same time
;
; 02/05/2013
; + writing fit values to log with 'g' extracts a 1d spectrum
; + entries in mospec.log contain filename of associated 1d spectrum
; + can fix redshift and line width using :z and :v (clear using :c)
; + rest-frame equivalent widths get printed to log with fluxes
; + pressing 'c' key now clears fit and returns current object to 
;   default display
; + line_fit checks SNR of data versus SNR from fit and displays
;   warning if the two are very discrepant
; + can begin working on a serendipitous object by typing ":s <>", 
;   where <> is some letter/string; will return 2d spectrum to 
;   default view
; + command-line counter for #/total when using the /all option
; + fixed 'b' option so it will not try to go past the first object
;
; 02/13/2013
; + fixed filename search string to accept products of the DRP that
;   haven't been combined or renamed yet
; + fixed redshift errors in km/s (were too high by (1+z)^2)
; + line flux errors are always positive
;
; 04/24/2013
; + changed line width reporting to include instrumental profile
;
; 05/27/2013
; + the /nocal option no longer looks for a cal vector and will run
;   without having the necessary *.psp.dat files
; + you can provide a fixed line width, and it will add back in the
;   instrumental profile before it does the fit
; + printing in line_fit ('g') now skips doing the fit again and just
;   reports the parameters from the previous fit ('f')
;
; 06/11/2013
; + when you fit the spectrum, it overplots the fit as well as the 
;   line identifications
; + when you define the aperture (with either 'l' or 't'), mospec 
;   saves it to the header of the 2D spectrum (in ROW and APERTURE)
; + extracting a spectrum (with either 'e' or 'g') now writes the
;   redshift to the header keyword REDSHIFT
;
; 07/17/2013
; + condition for 'sigflag' in line_fit.pro changed to reflect
;   what has been used in Chuck's measurements
;
; 07/19/2013
; + line width correction for instrumental resolution has been 
;   been turned off (commented out in line_fit.pro, so you can
;   add it back in, if you want)
; 
; 09/24/2013
; + the counter in mospec is now saved as a common block variable,
;   so you can resume with the object you were working on when
;   mospec crashes (the counter is reset if you .reset_session,
;   however)
;
; 10/02/2013
; + added new option: 'm' prints out the MASKS header keyword
;
; 10/24/2013
; + line_fit.pro now constrains the ratio of the nebular [OIII] lines
;   to be 2.98 (note that this makes the error on the weaker line
;   incorrect, since mpfit doesn't calculate errors the same way for
;   tied parameters
;
; 11/05/2013
; + changed the ratio of [OIII] 5007/4959 to 3.01, to be consistent
;   with the output from cloudy models on ramekin
;
; 12/09/2013
; + added line at zero flux for trace-finding interface
; + problem with rounding during trace-finding solved
; + direct-integration fluxes are now always printed
;
; 01/22/2014
; + switched over to using line_fit_mc.pro exclusively
;
; 02/24/2014
; + continuum is fit and then fit again simultaneously with the 
;   lines
; + the best-fit SED model spectrum for each galaxy can now be used
;   as the continuum in the fit
;
; 07/25/2014
; + ratio of nebular [OIII] lines and nebular [NII] lines is fixed
;   to 3.0 and errors for both lines should reflect the fact that 
;   they are tied variables
; + slit star diagnostic ('d') now uses only the middle 50% of the
;   spectrum to calculate the synthetic magnitude
; 
; 07/28/2014
; + error in the calculation of slitloss corrections stemming from 
;   a missing factor of two in the argument for the term accounting
;   for a finite extraction aperture has been fixed; results in a
;   ~4-5% 
;
; 08/11/2014
; + slit star magnitudes are now calculated by converting F_lambda to
;   F_nu, rejecting points >3-sigma away from the median, and taking
;   the mean, weighted by the filter transmission curve
;
; 10/20/2014 
; + update to mospec and line_fit_mc complete, see details below:
;; 1) The 2D fitting portion of the code only tries to display and write information to the header if the profile is well-defined (>10 sigma in the measurement of the X and Y Gaussian sigmas). It is also zoomed in greatly, so you can see what's happening with the 2D fit.
;; 2) Masking will only occur for regions of the spectrum that deviate significantly from the fit AND are at wavelengths less than 22,500 AA. This will prevent large portions of the thermal end of K from being excluded unnecessarily. Use the '/' key to show the points that would be flagged. This has been changed in stack_1d as well.
;; 3) The continuum polynomial fit is now not weighted by the error array, which makes it less likely to be absurdly steep.
;; 4) The optimal extraction option uses a collapsed portion of the spectrum instead of a 2D fit to a line, which both makes it more robust and means it can be used for continuum objects (like the stars!). Quick checks (on continuum bright sources) show that it is doing what it should: providing essentially the same spectrum with lower errors.
;; 5) There was a mistake in the formula for calculating the line width error, which is now fixed; it actually reduces the reported error.
;; 6) For objects with a model continuum, mospec now calculates the ratio of H-alpha or H-beta fluxes as measured with and without the SED continuum, in order to quantify the contribution from Balmer absorption features (this isn't currently being printed to log just yet...)
;; 7) Most references to file names and locations have been changed to be less rigid, and mospec should now work in pretty much any location where there are *wtavg.fits, *flx.fits, and *eps.fits files, so the onus of making sure that only 2D spectra from a given field will have results written to file is on us, rather than the code (this has never been a problem anyway, and I like the greater flexibility). The slitstars.log and mospec.log files are written one directory up. 
;; 8) A synthetic magnitude is recorded in the header of spectra when 1d spectrum is extracted with 'g'
;
; 04/14/2015
; + the stretch on the 2D spectrum has been changed from 1/99% to 3/97%. In addition, it will rescale when you zoom in, which helps when there is a bright star on the slit with the faint emission line object
;
; 05/12/2015
; + updated fitting routine to have a fixed velocity (in km/s) instead
; of a fixed width (in Ang); this required a change to line_fit_mc.pro
; and mospec.pro (and mospec1d.pro)
; + smoothing is now done using a Gaussian smoothing algorithm instead
; of a boxcar
;
; 07/15/2015
; + updated [OII] and [SII] measurements to include the most recent
; atomic data (references in line_fit_mc.pro)
; + errors adjusted on [SII] and [OII] ratios pegged at the limit to
; include equivalent fractional errors on the ratio as the measured
; line (previously errors on pegged ratios were 0)
;
; 10/23/2015
; + mospec opens *eps.fits files again, after not being able to
; appropriately parse filenames and headers
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

