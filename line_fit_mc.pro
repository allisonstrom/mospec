function boot_spec,infile,range,seed

  readcol,infile,file,format='A',/silent
  
  randomindex = fix(randomu(seed,n_elements(file))*n_elements(file))
  
  openw,101,'tmpspec.dat'
  for n=0,n_elements(file)-1 do printf,101,file[randomindex[n]]
  close,101
  stack_1d,'tmpspec.dat',/mean,/sigclip,/mask,/silent
  spawn,'/bin/rm tmpspec.dat'
  
  file = 'tmpspec.spec.fits'
  im1 = readfits(file,h1,exten=0,/silent)
  dim = size(im1,/dimen)
  stddev = readfits((strsplit(file,'.',/extract))[0]+'.stddev.fits',h2,exten=0,/silent)
  nspec = readfits((strsplit(file,'.',/extract))[0]+'.nspec.fits',h2,exten=0,/silent) & nspec = nspec[*,0]
  spawn,'/bin/rm tmpspec*fits'
  
  lambda = dblarr(dim[0])
  refpix = sxpar(h1,'CRPIX1')
  lam0 = sxpar(h1,'CRVAL1')
  delta = sxpar(h1,'CD1_1')
  for n=long(0),long(n_elements(lambda))-1 do lambda[n]=lam0+(n-(refpix-1))*delta
  spec = im1[*,0]
  errspec = stddev[*,0]
  
  zeroindex = where(abs(spec) gt 0 and errspec gt 0 and finite(errspec) and nspec gt 0.5*max(nspec) and lambda ge min(range) and lambda le max(range))
  
  return,[[lambda[zeroindex]],[spec[zeroindex]],[errspec[zeroindex]]]
  
end


;;;;;;;;;;;;;;;;;;;;;;

pro line_fit_mc,lambda,flux,err,startwav,startpeak,startline,full=full,result=result,perror=perror,covar=covar,redshift=redshift,redxi2=redxi2,print=print,zflag=zflag,vflag=vflag,fixredshift=fixredshift,fixvelocity=fixvelocity,bootflag=bootflag,bootfile=bootfile,stackflag=stackflag,sedcont=modelspec2,silent=silent,contmed=contmed,cpoly=cpoly,refit=refit,twocomp=twocomp,empirical=empirical

scale = 1.0

if keyword_set(full) then full=1 else full=0
if keyword_set(print) then print=1 else print=0
if keyword_set(zflag) then zflag=1 else zflag=0
if keyword_set(vflag) then vflag=1 else vflag=0
if keyword_set(bootflag) then bootflag=1 else bootflag=0
if keyword_set(silent) then silent=1 else silent=0
if keyword_set(stackflag) then stackflag=1 else stackflag=0
if n_elements(cpoly) eq 0 then cpoly = 1
if keyword_set(refit) then refit=1 else refit=0
if keyword_set(twocomp) then twocomp=1 else twocomp=0
if keyword_set(empirical) then empirical=1 else empirical=0

;; if filter eq 'Y' then inst_sig = 0.0 ;;;FILL IN CORRECT VALUE
;; if filter eq 'J' then inst_sig = 37.7
;; if filter eq 'H' then inst_sig = 35.2
;; if filter eq 'K' then inst_sig = 34.5
inst_sig = 35.0

if zflag then redshift = fixredshift
if vflag then linewidth = fixvelocity

;load local line list or general vacuum line list
calibdir = getenv('MOSPEC_CALIB')
if strmid(calibdir,n_elements(calibdir)-1) ne '/' then calibdir += '/'
if file_search('lines.dat') ne '' then readcol,'lines.dat',lines,airwav,vacwav,format='A,D',/silent 
if file_search('lines.dat') eq '' then begin
   if full then readcol,calibdir+'rest_optical_emlines_vac.dat',lines,airwav,vacwav,format='A,D,D',/silent
   if not full then readcol,calibdir+'rest_optical_lines_vac_part.dat',lines,airwav,vacwav,format='A,D,D',/silent
endif
linewav=vacwav
lines = lines[sort(linewav)]
linewav = linewav[sort(linewav)]
mindiff = min(abs(linewav-startline),index)
startline = linewav[index]

modelflag = 0
if keyword_set(modelspec2) then modelflag = 1

if print then begin
   wave = []
   line = []
   for n=0,(n_elements(result)-(cpoly+3))/2.-1 do begin
      wave = [wave,result[2*n+4+cpoly]]
      mindiff = min(abs(linewav-wave[n]),index)
      line = [line,lines[index]]
   endfor
   redshift = result[cpoly+1]
endif

if not print then begin

   if bootflag then begin
      bootfile = (strsplit(bootfile,'.',/extract))[0]+'.dat'
      nruns=1000
      seeds = fix(randomu(9876543,nruns)*nruns)
      mcchain = []
      openw,26,(strsplit(bootfile,'.',/extract))[0]+'.chain'
   endif else nruns=1
   origlambda = lambda
   origflux = flux
   origerr = err

   for run=0,nruns-1 do begin
      if run gt 0 then print,run
      if run lt nruns-1 then begin
         newspec = boot_spec(bootfile,minmax(origlambda),seeds[run])
         lambda=newspec[*,0]
         flux=newspec[*,1]
         err=newspec[*,2]
      endif
      if run eq nruns-1 then begin
         lambda=origlambda
         flux=origflux
         err=origerr
      endif

;robust fit to single line
      if not refit then begin
         mindiff = min(abs(lambda-(startwav-100)),lindex)
         mindiff = min(abs(lambda-(startwav+100)),uindex)
         lambda2 = lambda[lindex:uindex]
         flux2 = flux[lindex:uindex]
         err2 = err[lindex:uindex]
         expr = 'P[0]+P[1]*exp(-(X-P[2])^2/(2*(P[3]*P[2])^2))'
         start = [0.D,startpeak/scale,startwav,80/3.e5]
         pi = replicate({fixed:0,limited:[0,0],limits:[0.D,0.D],tied:''},4)
         if zflag then start[2] = (1+redshift)*startline
         if vflag then start[3] = linewidth/(2.99792458e5)
         result = mpfitexpr(expr,lambda2,flux2/scale,err2/scale,start,yfit=yfit2,/quiet,status=status,errmsg=errmsg,parinfo=pi,perror=perror2)
         linesize = abs(5.0*result[2]*result[3])
         pixels = where(abs(lambda2-result[2]) lt linesize)
         profile = flux2[pixels]-result[0] & profile /= total(profile)
         if status le 0 then message,errmsg
         if not vflag and (result[3] lt inst_sig/3.e5 or result[3] gt 1) then result[3] = start[3]
         yfit2 *= scale
         ;constant,max,mean,sigma
      endif else begin
         redshift = fixredshift
         if median(flux) lt 0.01 then peakflux = 0.01 else peakflux = median(flux)
         if not vflag then result = [0.D,peakflux,0,80/3.e5]
         if vflag then result = [0.D,peakflux,0,linewidth/(2.99792458e5)]
      endelse
         
;make guess at redshift
      if not zflag and not refit then redshift = result[2]/startline-1

;locate regions near common strong lines
      tempindex = intarr(n_elements(lambda))
      for  n=0,n_elements(lines)-1 do begin
         index = where(abs(lambda-(1+redshift)*linewav[n]) lt 5*result[3])
         tempindex[index] = 1
      endfor
      cindex = where(tempindex eq 0 and lambda gt min(lambda)+5 and lambda lt max(lambda)-5)
      contmed = median(flux[cindex])
      lineindex = where(tempindex eq 1)
      
;fit continuum
      if not modelflag then begin
         clambda = lambda[cindex]
         cflux = flux[cindex]
         cerr = err[cindex]*0+1.
         for n=cpoly,cpoly do begin
            expr = 'POLY(X,P[0:'+strtrim(string(n),1)+'])'
            start = replicate(0.D,n+1)
            continuum = mpfitexpr(expr,clambda,cflux/scale,cerr/scale,start,yfit=cfit,/quiet,dof=dof,bestnorm=xi,perror=cperror,covar=covar,status=status,errmsg=errmsg)
            if status le 0 then message,errmsg
            cfit *= scale
         end
      endif

;simultaneous fit to multiple lines
      line = lines[where(linewav gt min(origlambda/(redshift+1.))+10 and linewav lt max(origlambda/(redshift+1.))-5)]
      wave = linewav[where(linewav gt min(origlambda/(redshift+1.))+10 and linewav lt max(origlambda/(redshift+1.))-5)]
      if twocomp then begin
         line = [line,line]
         wave = [wave,wave]
      endif

      if not modelflag then begin
         expr = 'POLY(X,P[0:'+strtrim(string(cpoly),1)+'])+'
         pi = replicate({fixed:0,limited:[0,0],limits:[0.D,0.D],tied:''},2*n_elements(line)+3+cpoly)
         for n=0,cpoly do pi[n].fixed = 1
         start = [continuum[0:cpoly],redshift,result[3]]
      endif else begin
         expr = ''
         pi = replicate({fixed:0,limited:[0,0],limits:[0.D,0.D],tied:''},2*n_elements(line)+2)
         start = [redshift,result[3]]
         flux -= modelspec2
      endelse
      if zflag then begin
         start[cpoly+1] = redshift
         if redshift gt 0 then begin
            ;; pi[cpoly+1].limited = [1,1]
            ;; pi[cpoly+1].limits = [start[cpoly+1]-(1+start[cpoly+1])*2e-4,start[cpoly+1]+(1+start[cpoly+1])*2e-4]
            pi[cpoly+1].fixed = 1
         endif else pi[cpoly+1].fixed = 1
      endif
      if vflag then begin
         start[cpoly+2] = linewidth/(2.99792458e5)
         pi[cpoly+2].limited = [1,1]
         pi[cpoly+2].limits = [0.9,1.1]*start[cpoly+2]
      endif else begin
         pi[cpoly+2].limited = [1,0]
         pi[cpoly+2].limits = [inst_sig/(2.99792458e5),0]
      endelse         

      ; fixing the ratio of the nebular [OIII] lines
      line_index = where(line eq '[OIII]')
      o3_index = []
      if n_elements(line_index) ge 2 then begin
         for n=0,n_elements(line_index)-1 do begin
            if abs(wave[line_index[n]]-4959) le 2 then o3_index = [o3_index,line_index[n]]
            if abs(wave[line_index[n]]-5007) le 2 then o3_index = [o3_index,line_index[n]]
         endfor
      o3_index = 2*o3_index+3+cpoly
      endif 
      if not twocomp and n_elements(o3_index) eq 2 then for n=0,n_elements(o3_index)-1 do pi[o3_index[0]].tied = 'P['+strtrim(string(o3_index[1]),1)+']/3.0*P['+strtrim(string(o3_index[1]+1),1)+']/P['+strtrim(string(o3_index[0]+1),1)+']'
      if twocomp and keyword_set(o3_index) and (n_elements(o3_index) mod 2) eq 0 then begin
         total_o3_index = o3_index
         o3_index = total_o3_index[0:1]
         for n=0,n_elements(o3_index)-1 do pi[o3_index[0]].tied = 'P['+strtrim(string(o3_index[1]),1)+']/3.0*P['+strtrim(string(o3_index[1]+1),1)+']/P['+strtrim(string(o3_index[0]+1),1)+']'
         o3_index = total_o3_index[-2:-1]
         for n=0,n_elements(o3_index)-1 do pi[o3_index[0]].tied = 'P['+strtrim(string(o3_index[1]),1)+']/3.0*P['+strtrim(string(o3_index[1]+1),1)+']/P['+strtrim(string(o3_index[0]+1),1)+']'
      endif
      
      ; fixing the ratio of the [NII] lines
      line_index = where(line eq '[NII]')
      n2_index = []
      if n_elements(line_index) ge 2 then begin
         for n=0,n_elements(line_index)-1 do begin
            if abs(wave[line_index[n]]-6548) le 2 then n2_index = [n2_index,line_index[n]]
            if abs(wave[line_index[n]]-6584) le 2 then n2_index = [n2_index,line_index[n]]
         endfor
      n2_index = 2*n2_index+3+cpoly
      endif 
      if not twocomp and n_elements(n2_index) eq 2 then for n=0,n_elements(n2_index)-1 do pi[n2_index[0]].tied = 'P['+strtrim(string(n2_index[1]),1)+']/3.0*P['+strtrim(string(n2_index[1]+1),1)+']/P['+strtrim(string(n2_index[0]+1),1)+']'
      if twocomp and keyword_set(n2_index) and (n_elements(n2_index) mod 2) eq 0 then begin
         total_n2_index = n2_index
         n2_index = total_n2_index[0:1]
         for n=0,n_elements(n2_index)-1 do pi[n2_index[0]].tied = 'P['+strtrim(string(n2_index[1]),1)+']/3.0*P['+strtrim(string(n2_index[1]+1),1)+']/P['+strtrim(string(n2_index[0]+1),1)+']'
         n2_index = total_n2_index[-2:-1]
         for n=0,n_elements(n2_index)-1 do pi[n2_index[0]].tied = 'P['+strtrim(string(n2_index[1]),1)+']/3.0*P['+strtrim(string(n2_index[1]+1),1)+']/P['+strtrim(string(n2_index[0]+1),1)+']'
      endif
      
      ; identifying [OII] lines, restricting ratio to physical range
      line_index = where(line eq '[OII]')
      o2_index = []
      o2_flag = 0
      if n_elements(line_index) ge 2 then begin
         for n=0,n_elements(line_index)-1 do begin
            if abs(wave[line_index[n]]-3726) le 1.5 then o2_index = [o2_index,line_index[n]]
            if abs(wave[line_index[n]]-3729) le 1.5 then o2_index = [o2_index,line_index[n]]
         endfor
      endif
      if not twocomp and n_elements(o2_index) eq 2 then begin
         o2_index = 2*o2_index+3+cpoly
         pi[o2_index[1]].limited = [1,1]
         ;; pi[o2_index[1]].limits = [0.35,1.5]
         pi[o2_index[1]].limits = [6./4*(3.382e-5+7.416e-6)/(2.209e-5+1.414e-4),0.834/0.554]*wave[(o2_index[0]-(cpoly+3))/2]/wave[(o2_index[1]-(cpoly+3))/2]
         pi[o2_index[1]].limits = [0.8*pi[o2_index[1]].limits[0],1.2*pi[o2_index[1]].limits[1]]
         ;; Froese Fischer & Tachiev (2004), Kisielius et al. (2009) @ 10kK
      endif
      if twocomp and keyword_set(o2_index) and (n_elements(o2_index) mod 2) eq 0 then begin
         total_o2_index = 2*o2_index+3+cpoly
         o2_index = total_o2_index[0:1]
         pi[o2_index[1]].limited = [1,1]
         pi[o2_index[1]].limits = [6./4*(3.382e-5+7.416e-6)/(2.209e-5+1.414e-4),0.834/0.554]*wave[(o2_index[0]-(cpoly+3))/2]/wave[(o2_index[1]-(cpoly+3))/2]
         pi[o2_index[1]].limits = [0.8*pi[o2_index[1]].limits[0],1.2*pi[o2_index[1]].limits[1]]
         o2_index = total_o2_index[-2:-1]
         pi[o2_index[1]].limited = [1,1]
         pi[o2_index[1]].limits = [6./4*(3.382e-5+7.416e-6)/(2.209e-5+1.414e-4),0.834/0.554]*wave[(o2_index[0]-(cpoly+3))/2]/wave[(o2_index[1]-(cpoly+3))/2]
         pi[o2_index[1]].limits = [0.8*pi[o2_index[1]].limits[0],1.2*pi[o2_index[1]].limits[1]]
      endif
      if (n_elements(o2_index) mod 2) gt 0 or not keyword_set(o2_index) then o2_index = [-1,-1]
      
      ; identifying [SII] lines, restricting ratio to physical range
      line_index = where(line eq '[SII]')
      s2_index = []
      s2_flag = 0
      if n_elements(line_index) ge 2 then begin
         for n=0,n_elements(line_index)-1 do begin
            if abs(wave[line_index[n]]-6716) le 3 then s2_index = [s2_index,line_index[n]]
            if abs(wave[line_index[n]]-6731) le 3 then s2_index = [s2_index,line_index[n]]
         endfor
      endif
      if n_elements(s2_index) eq 2 then begin
         s2_index = 2*s2_index+3+cpoly
         pi[s2_index[0]].limited = [1,1]
         ;; pi[s2_index[0]].limits = [0.43,1.5]
         pi[s2_index[0]].limits = [6./4*(2.66e-4/8.95e-4),3.83/2.56]*wave[(s2_index[1]-(cpoly+3))/2]/wave[(s2_index[0]-(cpoly+3))/2]
         pi[s2_index[0]].limits = [0.8*pi[s2_index[0]].limits[0],1.2*pi[s2_index[0]].limits[1]]
         ;; Mendoza & Bautista (2014), Tayal & Zatsarinny (2010) @ 10kK
      endif else s2_index = [-1,-1]

      
      if not twocomp then for n=0,n_elements(line)-1 do begin
         if (2*n+3+cpoly eq s2_index[0] and n_elements(s2_index) eq 2) then begin
            expr = expr+'P['+strtrim(string(2*(n+1)+3+cpoly),1)+']*P['+strtrim(string(2*n+3+cpoly),1)+']*exp(-(X-(1+P['+strtrim(string(1+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2/(2*(P['+strtrim(string(2+cpoly),1)+']*(1+P['+strtrim(string(1+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2))'
            start = [start,1.0,wave[n]]
         endif else if (2*n+3+cpoly eq o2_index[1] and n_elements(o2_index) eq 2) then begin
            expr = expr+'P['+strtrim(string(2*(n-1)+3+cpoly),1)+']*P['+strtrim(string(2*n+3+cpoly),1)+']*exp(-(X-(1+P['+strtrim(string(1+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2/(2*(P['+strtrim(string(2+cpoly),1)+']*(1+P['+strtrim(string(1+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2))'
            start = [start,1.0,wave[n]]
         endif else begin
            expr = expr+'P['+strtrim(string(2*n+3+cpoly),1)+']*exp(-(X-(1+P['+strtrim(string(1+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2/(2*(P['+strtrim(string(2+cpoly),1)+']*(1+P['+strtrim(string(1+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2))'
            start = [start,result[1],wave[n]]
         endelse
         if n lt n_elements(line)-1 then expr += '+'
         pi[2*n+4+cpoly].fixed = 1
;;          pi[2*n+3+cpoly].limited = [1,0]
;;          pi[2*n+3+cpoly].limits = [0.0,0.0]
      endfor
      
      if twocomp then for n=0,n_elements(line)-1 do begin
         pi[2*n+3+cpoly].limited = [1,0]
         pi[2*n+3+cpoly].limits = [0.0,0.0]
         if n lt n_elements(line)/2 then begin
            if keyword_set(total_o2_index) then o2_index = total_o2_index[0:1]
            if (2*n+3+cpoly eq s2_index[0] and n_elements(s2_index) eq 2) then begin
               expr = expr+'P['+strtrim(string(2*(n+1)+3+cpoly),1)+']*P['+strtrim(string(2*n+3+cpoly),1)+']*exp(-(X-(1+P['+strtrim(string(1+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2/(2*(P['+strtrim(string(2+cpoly),1)+']*(1+P['+strtrim(string(1+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2))'
               start = [start,1.0,wave[n]]
            endif else if (2*n+3+cpoly eq o2_index[1] and n_elements(o2_index) eq 2) then begin
               expr = expr+'P['+strtrim(string(2*(n-1)+3+cpoly),1)+']*P['+strtrim(string(2*n+3+cpoly),1)+']*exp(-(X-(1+P['+strtrim(string(1+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2/(2*(P['+strtrim(string(2+cpoly),1)+']*(1+P['+strtrim(string(1+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2))'
               start = [start,1.0,wave[n]]
            endif else begin
               expr = expr+'P['+strtrim(string(2*n+3+cpoly),1)+']*exp(-(X-(1+P['+strtrim(string(1+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2/(2*(P['+strtrim(string(2+cpoly),1)+']*(1+P['+strtrim(string(1+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2))'
               start = [start,result[1],wave[n]]
            endelse
         endif else begin
            if keyword_set(total_o2_index) then o2_index = total_o2_index[-2:-1]
            if (2*n+3+cpoly eq s2_index[0] and n_elements(s2_index) eq 2) then begin
               expr = expr+'P['+strtrim(string(2*(n+1)+3+cpoly),1)+']*P['+strtrim(string(2*n+3+cpoly),1)+']*exp(-(X-(1+P['+strtrim(string(2*n_elements(line)+4+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2/(2*(P['+strtrim(string(2*n_elements(line)-1+4+cpoly),1)+']*(1+P['+strtrim(string(2*n_elements(line)+4+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2))'
               start = [start,1.0,wave[n]]
            endif else if (2*n+3+cpoly eq o2_index[1] and n_elements(o2_index) eq 2) then begin
               expr = expr+'P['+strtrim(string(2*(n-1)+3+cpoly),1)+']*P['+strtrim(string(2*n+3+cpoly),1)+']*exp(-(X-(1+P['+strtrim(string(2*n_elements(line)+4+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2/(2*(P['+strtrim(string(2*n_elements(line)-1+4+cpoly),1)+']*(1+P['+strtrim(string(2*n_elements(line)+4+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2))'
               start = [start,1.0,wave[n]]
            endif else begin
               expr = expr+'P['+strtrim(string(2*n+3+cpoly),1)+']*exp(-(X-(1+P['+strtrim(string(2*n_elements(line)+4+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2/(2*(P['+strtrim(string(2*n_elements(line)-1+4+cpoly),1)+']*(1+P['+strtrim(string(2*n_elements(line)+4+cpoly),1)+'])*P['+strtrim(string(2*n+4+cpoly),1)+'])^2))'
               start = [start,result[1],wave[n]]
            endelse
         endelse
         if n lt n_elements(line)-1 then expr += '+'
         pi[2*n+4+cpoly].fixed = 1
      endfor
      if twocomp then begin
         if keyword_set(total_o2_index) then o2_index = total_o2_index[0:1]
         start = [start,5*start[2+cpoly],0.0]
         pi = [pi,pi[2+cpoly],pi[1+cpoly]]
         pi[-1].fixed = 0
      endif
      
      lambda3 = lambda;[lineindex]
      flux3 = flux;[lineindex]
      err3 = err  ;[lineindex]
      result = mpfitexpr(expr,lambda3,flux3/scale,err3/scale,start,parinfo=pi,yfit=yfit,bestnorm=xi,dof=dof,status=status,perror=perror,covar=covar,niter=niter,/quiet,errmsg=errmsg,nfree=nfree)
      if status le 0 then message,errmsg
      if where(~finite(covar)) eq -1 then begin
         if not modelflag then perror[0:cpoly] = cperror[0:cpoly]
         if n_elements(o3_index) eq 2 then begin
            perror[o3_index[0]] = perror[o3_index[1]]/3.0*result[o3_index[1]+1]/result[o3_index[0]+1]
            covar[o3_index[0],*] = covar[o3_index[1],*]/3.0*result[o3_index[1]+1]/result[o3_index[0]+1]
            covar[*,o3_index[0]] = covar[*,o3_index[1]]/3.0*result[o3_index[1]+1]/result[o3_index[0]+1]
         endif
         if n_elements(n2_index) eq 2 then begin
            perror[n2_index[0]] = perror[n2_index[1]]/3.0*result[n2_index[1]+1]/result[n2_index[0]+1]
            covar[n2_index[0],*] = covar[n2_index[1],*]/3.0*result[n2_index[1]+1]/result[n2_index[0]+1]
            covar[*,n2_index[0]] = covar[*,n2_index[1]]/3.0*result[n2_index[1]+1]/result[n2_index[0]+1]
         endif
         if n_elements(s2_index) eq 2 and pi[s2_index[0]].limits[0] gt 0 then begin
            if perror[s2_index[0]] eq 0 then s2_flag = 1
            if s2_flag and not silent then print,'WARNING: [SII] ratio pegged at '+strtrim(string(result[s2_index[0]],format='(F0.2)'),1)
            if s2_flag then perror[s2_index[0]] = perror[s2_index[1]]/result[s2_index[1]]*result[s2_index[0]]
            perror[s2_index[0]] = result[s2_index[0]]*result[s2_index[1]]*sqrt((perror[s2_index[0]]/result[s2_index[0]])^2.+(perror[s2_index[1]]/result[s2_index[1]])^2.+2*covar[s2_index[0],s2_index[1]]/(result[s2_index[0]]*result[s2_index[1]]))
            new_covar = result[s2_index[0]]*covar[s2_index[1],cpoly+2]+result[s2_index[1]]*covar[cpoly+2,s2_index[0]]
            covar[s2_index[0],cpoly+2] = new_covar
            covar[cpoly+2,s2_index[0]] = new_covar
            new_covar = result[s2_index[0]]*covar[s2_index[1],cpoly+1]+result[s2_index[1]]*covar[cpoly+1,s2_index[0]]
            covar[s2_index[0],cpoly+1] = new_covar
            covar[cpoly+1,s2_index[0]] = new_covar
            result[s2_index[0]] = result[s2_index[0]]*result[s2_index[1]]
            if s2_flag then perror[s2_index[1]] = perror[s2_index[0]]/result[s2_index[0]]*result[s2_index[1]]
         endif
         if n_elements(o2_index) eq 2 and pi[o2_index[1]].limits[0] gt 0 then begin
            if perror[o2_index[1]] eq 0 then o2_flag = 1
            if o2_flag and not silent then print,'WARNING: [OII] ratio pegged at '+strtrim(string(result[o2_index[1]],format='(F0.2)'),1)
            if o2_flag then perror[o2_index[1]] = perror[o2_index[0]]/result[o2_index[0]]*result[o2_index[1]]
            perror[o2_index[1]] = result[o2_index[0]]*result[o2_index[1]]*sqrt((perror[o2_index[0]]/result[o2_index[0]])^2.+(perror[o2_index[1]]/result[o2_index[1]])^2.+2*covar[o2_index[0],o2_index[1]]/(result[o2_index[0]]*result[o2_index[1]]))
            new_covar = result[o2_index[1]]*covar[o2_index[0],cpoly+2]+result[o2_index[0]]*covar[cpoly+2,o2_index[1]]
            covar[o2_index[1],cpoly+2] = new_covar
            covar[cpoly+2,o2_index[1]] = new_covar
            new_covar = result[o2_index[1]]*covar[o2_index[0],cpoly+1]+result[o2_index[0]]*covar[cpoly+1,o2_index[1]]
            covar[o2_index[1],cpoly+1] = new_covar
            covar[cpoly+1,o2_index[1]] = new_covar            
            result[o2_index[1]] = result[o2_index[0]]*result[o2_index[1]]
            if keyword_set(total_o2_index) then result[total_o2_index[-1]] = result[total_o2_index[-2]]*result[total_o2_index[-1]]
            if o2_flag then perror[o2_index[0]] = perror[o2_index[1]]/result[o2_index[1]]*result[o2_index[0]]
         endif
         redshift = result[cpoly+1]
         xfit = lambda3
         yfit *= scale
         if not modelflag then begin
            result[0:cpoly] *= scale
            perror[0:cpoly] *= scale
         endif
         for n=0,n_elements(line)-1 do begin 
            result[2*n+cpoly+3] *= scale
            perror[2*n+cpoly+3] *= scale
            covar[*,2*n+cpoly+3] *= scale
            covar[2*n+cpoly+3,*] *= scale
         endfor
         redxi2 = xi/dof
         
         if bootflag then begin
            linewidth = result[cpoly+2]*(2.99792458e5)/((1+result[cpoly+1])*result[cpoly+4])
            tmpmcchain = [linewidth]
            for l=0,(n_elements(result)-(cpoly+3))/2.-1 do tmpmcchain = [tmpmcchain,result[2*l+cpoly+3]*result[cpoly+2]*sqrt(2*!DPI)*(1.0e-17)]
            printf,26,tmpmcchain,format='(%"%e %e %e %e %e %e")'
            mcchain = [[mcchain],[tmpmcchain]]
         endif
      endif else run -= 1
   endfor
   close,26
endif

if bootflag then begin
   for n=0,(n_elements(result)-(cpoly+3))/2. do begin
      index = where(mcchain[n,*] le 0, subzero)
      f = 1-float(subzero)/n_elements(mcchain[n,*])
      interval = conf_int(mcchain[n,*],/minrange)
      print,mean(interval),interval[1]-mean(interval),2*mean(interval)/(interval[1]-interval[0]),f,sqrt(2.)*inverf(f)
   endfor
endif

if modelflag then begin
   continuum = modelspec2
   if not print then flux += modelspec2
endif else continuum = poly(lambda,result[0:cpoly])

linewidth = result[cpoly+2]*(2.99792458e5)
lineerror = perror[cpoly+2]*(2.99792458e5)
redshifterror = perror[cpoly+1]*(2.99792458e5);/(1+result[cpoly+1])
format1 = '(%"Redshift: %0.4f  z_error: %-0.1f km/s  Line width: %-0.2f +/- %-0.2f km/s  xi2: %-0.2f")'
format2 = '(%"%-10s %8.3f  %9.2e +/- %-9.2e  %8.3f")'
if zflag then format1 = '(%"Redshift*: %6.4f  z_error: %-0.1f km/s  Line width: %-0.2f +/- %-0.2f km/s  xi2: %-0.2f")'
if vflag then format1 = '(%"Redshift: %6.4f  z_error: %-0.1f km/s  Line width*: %-0.2f +/- %-0.2f km/s  xi2: %-0.2f")'
if zflag and vflag then format1 = '(%"Redshift*: %6.4f  z_error: %-0.1f km/s  Line width*: %-0.2f +/- %-0.2f km/s  xi2: %-0.2f")'

;; fluxes = []
;; for n=0,(n_elements(result)-(cpoly+3))/2.-1 do fluxes = [fluxes,result[2*n+cpoly+3]]
;; strongindex = (where(fluxes eq max(fluxes)))[0]
;; if stackflag then begin
;;    strongline = line[strongindex]
;;    strongwave = wave[strongindex]
;;    index = where(lambda ge strongwave-10 and lambda le strongwave+10)
;;    x = lambda[index]
;;    y = ((flux-continuum)/err)[index]
;;    devindex = where(y le 1.0 and x le strongwave)
;;    linerange = max(x[devindex])
;;    devindex = where(y le 1.0 and x ge strongwave)
;;    linerange = [linerange,min(x[devindex])]
;; ;;    sortindex = sort(x[devindex]-strongwave)
;; ;;    linerange = [x[devindex[sortindex[0]]]]
;; ;;    sortindex = sort(strongwave-x[devindex])
;; ;;    linerange = [linerange,x[devindex[sortindex[-1]]]]
;;    linesize = (linerange[1]-linerange[0])/2.
;; endif else linesize = abs(3.0*result[cpoly+2]*startline*(1+redshift))
;; if linesize lt 3*(lambda[1]-lambda[0]) then linesize = 3*(lambda[1]-lambda[0])

sigprint = 0
if print then printf,1,redshift,redshifterror,linewidth,lineerror,redxi2,format=format1
if not print and not silent then print,'LINE       LAMBDA     FLUX                         SNR'
;; for n=0,(n_elements(result)-(cpoly+5))/2.-1 do begin
for n=0,n_elements(line)-1 do begin
   sigflag = 0

   lineflux = result[2*n+cpoly+3]*result[cpoly+2]*result[2*n+cpoly+4]*(1+result[cpoly+1])*sqrt(2*!DPI)
   fluxerror = abs(lineflux)*sqrt((perror[2*n+cpoly+3]/result[2*n+cpoly+3])^2.+(perror[cpoly+2]/result[cpoly+2])^2.+(perror[cpoly+1]/(1+result[cpoly+1]))^2.+2*covar[2*n+cpoly+3,cpoly+2]/(result[2*n+cpoly+3]*result[cpoly+2])+2*covar[2*n+cpoly+3,cpoly+1]/(result[2*n+cpoly+3]*(1+result[cpoly+1]))+2*covar(cpoly+2,cpoly+1)/(result[cpoly+2]*(1+result[cpoly+1])))

   linesize = abs(3.0*result[cpoly+2]*result[2*n+cpoly+4]*(1+redshift))
   pixels = where(abs(lambda-(1+result[cpoly+1])*result[2*n+cpoly+4]) lt linesize)
   lineflux2 = int_tabulated(lambda[pixels],flux[pixels]-continuum[pixels])
;;    funct = result[2*n+cpoly+3]*exp(-(lambda-(1+result[cpoly+1])*result[2*n+cpoly+4])^2/(2*result[cpoly+2]^2))
;;    wgts = funct/lineflux
   dlambda = (max(lambda[pixels])-min(lambda[pixels]))/n_elements(pixels)
   fluxerror2 = sqrt(total(err[pixels]*err[pixels]))*dlambda
;;    fluxerror2 = sqrt(total(err*err*wgts))*n_elements(pixels)/(max(lambda[pixels])-min(lambda[pixels]))
;;    oplot,lambda[pixels],flux[pixels],color=fsc_color('orange'),psym=10
;;    stop
   
   if abs(lineflux2/fluxerror2) gt abs(1.25*lineflux/fluxerror) and (lineflux2/fluxerror2) gt 5 then begin
      sigflag = 1
      sigprint = 1
   endif
   if print then begin
      printf,1,line[n],wave[n],lineflux,fluxerror,abs(lineflux/fluxerror),format=format2
      printf,1,line[n],wave[n],lineflux2,fluxerror2,abs(lineflux2/fluxerror2),format=format2
   endif
   if not print and not silent then begin
      print,line[n],wave[n],lineflux,fluxerror,abs(lineflux/fluxerror),format=format2
      if sigprint then print,'',wave[n],lineflux2,fluxerror2,abs(lineflux2/fluxerror2),format=format2
   endif
endfor
if not print and not silent then print,redshift,redshifterror,linewidth,lineerror,redxi2,format=format1
if sigprint and not print and not silent then print,'SOME FEATURES MAY BE MORE SIGNIFICANT THAN LINE FITTING RESULTS INDICATE'

;; if not print and not stackflag and (modelflag or not modelflag) then begin
;;    set_plot,'X'
;;    window,1,xsize=1000,ysize=950
;;    multiplot,[2,2],xgap=0,/square
;;    scale = 5

;;    mindiff = min(abs(result-startline),line_index) & line_index -= 1

;;    index0 = line_index
;;    index1 = cpoly+1
;;    tmpcovar = [[covar[index0,index0],covar[index0,index1]],[covar[index1,index0],covar[index1,index1]]]
;;    range0 = [result[index0]-scale*sqrt(covar[index0,index0]),result[index0]+scale*sqrt(covar[index0,index0])]
;;    range1 = [result[index1]-scale*sqrt(covar[index1,index1]),result[index1]+scale*sqrt(covar[index1,index1])]
;;    parm0 = dindgen(101)/100.*(range0[1]-range0[0])+range0[0]
;;    parm1 = dindgen(101)/100.*(range1[1]-range1[0])+range1[0]
;;    cost = dblarr([n_elements(parm0),n_elements(parm1)])
;;    for i=0,n_elements(parm0)-1 do begin
;;       for j=0,n_elements(parm1)-1 do begin
;;          a = [[result[index0]-parm0[i]],[result[index1]-parm1[j]]]
;;          cost[i,j] = (transpose(a)##invert(tmpcovar))##a
;;       endfor
;;    endfor
;;    levels = [chisqr_cvf(0.32,nfree),chisqr_cvf(0.05,nfree)]
;;    contour,cost,parm0,parm1,c_labels=replicate(1,3),levels=levels,c_annotation=textoidl(['68%','95%']),/xstyle,/ystyle,xrange=range0,yrange=range1,ytitle=textoidl('z'),charsize=1.5,xminor=1,yminor=1,c_charsize=1.5
;;    oplot,[result[index0]],[result[index1]],psym=1
;;    multiplot
;;    multiplot
   
;;    index0 = line_index
;;    index1 = cpoly+2
;;    tmpcovar = [[covar[index0,index0],covar[index0,index1]],[covar[index1,index0],covar[index1,index1]]]
;;    range0 = [result[index0]-scale*sqrt(covar[index0,index0]),result[index0]+scale*sqrt(covar[index0,index0])]
;;    range1 = [result[index1]-scale*sqrt(covar[index1,index1]),result[index1]+scale*sqrt(covar[index1,index1])]
;;    parm0 = dindgen(101)/100.*(range0[1]-range0[0])+range0[0]
;;    parm1 = dindgen(101)/100.*(range1[1]-range1[0])+range1[0]
;;    cost = dblarr([n_elements(parm0),n_elements(parm1)])
;;    for i=0,n_elements(parm0)-1 do begin
;;       for j=0,n_elements(parm1)-1 do begin
;;          a = [[result[index0]-parm0[i]],[result[index1]-parm1[j]]]
;;          cost[i,j] = (transpose(a)##invert(tmpcovar))##a
;;       endfor
;;    endfor
;;    levels = [chisqr_cvf(0.32,nfree),chisqr_cvf(0.05,nfree)]
;;    contour,cost,parm0,parm1,c_labels=replicate(1,3),levels=levels,c_annotation=textoidl(['68%','95%']),/xstyle,/ystyle,xrange=range0,yrange=range1,xtitle=textoidl('n_{startline}'),ytitle=textoidl('v/c'),charsize=1.5,xminor=1,yminor=1,c_charsize=1.5
;;    oplot,[result[index0]],[result[index1]],psym=1
;;    multiplot

;;    index0 = cpoly+1
;;    index1 = cpoly+2
;;    tmpcovar = [[covar[index0,index0],covar[index0,index1]],[covar[index1,index0],covar[index1,index1]]]
;;    range0 = [result[index0]-scale*sqrt(covar[index0,index0]),result[index0]+scale*sqrt(covar[index0,index0])]
;;    range1 = [result[index1]-scale*sqrt(covar[index1,index1]),result[index1]+scale*sqrt(covar[index1,index1])]
;;    parm0 = dindgen(101)/100.*(range0[1]-range0[0])+range0[0]
;;    parm1 = dindgen(101)/100.*(range1[1]-range1[0])+range1[0]
;;    cost = dblarr([n_elements(parm0),n_elements(parm1)])
;;    for i=0,n_elements(parm0)-1 do begin
;;       for j=0,n_elements(parm1)-1 do begin
;;          a = [[result[index0]-parm0[i]],[result[index1]-parm1[j]]]
;;          cost[i,j] = (transpose(a)##invert(tmpcovar))##a
;;       endfor
;;    endfor
;;    levels = [chisqr_cvf(0.32,nfree),chisqr_cvf(0.05,nfree)]
;;    contour,cost,parm0,parm1,c_labels=replicate(1,3),levels=levels,c_annotation=textoidl(['68%','95%']),/xstyle,/ystyle,xrange=range0,yrange=range1,xtitle=textoidl('z'),charsize=1.5,xminor=1,yminor=1,c_charsize=1.5
;;    oplot,[result[index0]],[result[index1]],psym=1
;;    multiplot,/reset

;;    print,'Press any key to close window and continue'
;;    key = ''
;;    key = get_kbrd()
;;    wdelete,1
;;    wset,4
;; endif

end
