
function data_interp0,abscissa,data,abscissa1
  n_wl=n_elements(abscissa1)
  ref=make_array(n_wl,/double)
  for i=0,n_wl-1 do begin
    t=abscissa1[i]
    index1=where(abscissa le t,count1)
    index2=where(abscissa ge t,count2)
    if count1 gt 0 and count2 gt 0 then begin
      wvlen1_index=max(index1)
      wvlen2_index=min(index2)
      wvlen1=abscissa[wvlen1_index]
      wvlen2=abscissa[wvlen2_index]
      if wvlen2 eq wvlen1 then begin
        ref[i]=data[wvlen1_index]
      endif else begin
        ref[i]=(abs(t-wvlen2)*data[wvlen1_index]+abs(wvlen1-t)*data[wvlen2_index])/abs(wvlen2-wvlen1)
      endelse
    endif
    if count1 gt 0 and count2 le 0 then begin
      ref[i]=data[max(index1)]
    endif
    if count1 le 0 and count2 gt 0 then begin
      ref[i]=data[min(index2)]
    endif
  endfor
  return,ref
end
function mk_array,dim,data_type
  if isa(data_type) then begin
    case data_type of
      1:begin
        arr=make_array(dim,/byte)
      end
      2:begin
        arr=make_array(dim,/int)
      end
      12:begin
        arr=make_array(dim,/uint)
      end    
      3:begin
        arr=make_array(dim,/long)
      end
      13:begin
        arr=make_array(dim,/ulong)
      end
      14:begin
        arr=make_array(dim,/long64)
      end
      15:begin
        arr=make_array(dim,/unlong64)
      end
      4:begin
        arr=make_array(dim,/float)
      end
      5:begin
        arr=make_array(dim,/double)
      end
    endcase
  endif else begin
    arr=make_array(dim,/int)
  endelse
  return,arr
end
function read_hdr_file,file,interleave,data_type,ns,nl,nb, $
  srange=srange,lrange=lrange,brange=brange,BYTE_SWAP=BYTE_SWAP,$
  return_ptr=return_ptr
  sizeof=indgen(16)
  sizeof[1]=1
  sizeof[2]=2
  sizeof[3]=4
  sizeof[4]=4
  sizeof[5]=8
  sizeof[12]=2
  sizeof[13]=4
  sizeof[14]=8
  b0=0
  b1=nb-1
  if keyword_set(brange) then begin
    if n_elements(brange) eq 2 then begin
      if brange[0] ge 0 and brange[1] lt nb and brange[0] le brange[1] then begin
        b0=brange[0]
        b1=brange[1]
      endif
    endif
  endif
  s0=0
  s1=ns-1
  if keyword_set(srange) then begin
    if n_elements(srange) eq 2 then begin
      if srange[0] ge 0 and srange[1] lt ns and srange[0] le srange[1] then begin
        s0=srange[0]
        s1=srange[1]
      endif
    endif
  endif
  l0=0
  l1=nl-1
  if keyword_set(lrange) then begin
    if n_elements(lrange) eq 2 then begin
      if lrange[0] ge 0 and lrange[1] lt nl and lrange[0] le lrange[1] then begin
        l0=lrange[0]
        l1=lrange[1]
      endif
    endif
  endif
  data=mk_array([s1-s0+1,l1-l0+1,b1-b0+1],data_type)
  if interleave eq 1 then begin
    line=mk_array(s1-s0+1,data_type)
    if keyword_set(BYTE_SWAP) then begin
      if BYTE_SWAP eq 1 then begin
        openr,lun,file,/get_lun,/SWAP_IF_LITTLE_ENDIAN
      endif else begin
        openr,lun,file,/get_lun
      endelse
    endif else begin
      openr,lun,file,/get_lun
    endelse
    
    for j=l0,l1 do begin
      for k=b0,b1 do begin
        point_lun,lun,(long64(nb)*long64(ns)*j+k*long64(ns)+s0)*sizeof[data_type]
        readu,lun,line
        data[*,j-l0,k-b0]=line
      endfor
    endfor
    free_lun,lun
  endif else if interleave eq 2 then begin
    spec=mk_array(b1-b0+1,data_type)
    if keyword_set(BYTE_SWAP) then begin
      if BYTE_SWAP eq 1 then begin
        openr,lun,file,/get_lun,/SWAP_IF_LITTLE_ENDIAN
      endif else begin
        openr,lun,file,/get_lun
      endelse
    endif else begin
      openr,lun,file,/get_lun
    endelse
    for j=l0,l1 do begin
      for i=s0,s1 do begin
        point_lun,lun,((long64(ns)*long64(j)+i)*long64(nb)+long64(b0))*sizeof[data_type]
        readu,lun,spec
        data[i-s0,j-l0,*]=spec
      endfor
    endfor
    free_lun,lun
  endif else begin
    line=mk_array(s1-s0+1,data_type)
    if keyword_set(BYTE_SWAP) then begin
      if BYTE_SWAP eq 1 then begin
        openr,lun,file,/get_lun,/SWAP_IF_LITTLE_ENDIAN
      endif else begin
        openr,lun,file,/get_lun
      endelse
    endif else begin
      openr,lun,file,/get_lun
    endelse
    for k=b0,b1 do begin
      for j=l0,l1 do begin
        point_lun,lun,(long64(ns)*long64(nl)*k+j*long64(ns)+s0)*sizeof[data_type]
        readu,lun,line
        data[*,j-l0,k-b0]=line
      endfor
    endfor
    free_lun,lun
  endelse
  if keyword_set(return_ptr) then begin
    if return_ptr then begin

      return,ptr_new(data,/NO_COPY)
    endif
  endif else begin
    return,data
  endelse
end
function hdr_file_query,file_name,$ ;stirng
  ACQUISITION_TIME=ACQUISITION_TIME,$ ;variable
  BBL=BBL,$ ;array
  BNAMES=BNAMES,$ ;variable
  BYTE_SWAP=BYTE_SWAP,$ ;variable
  CLASS_NAMES=CLASS_NAMES,$ ;variable
  CLOUD_COVER=CLOUD_COVER,$ ;variable
  DATA_GAINS=DATA_GAINS,$ ;variable
  DATA_IGNORE_VALUE=DATA_IGNORE_VALUE,$ ;variable
  DATA_OFFSETS=DATA_OFFSETS,$ ;variable
  DATA_TYPE=DATA_TYPE,$ ;integer
  DEF_BANDS=DEF_BANDS,$ ;array
  DEF_ZRANGE=DEF_ZRANGE,$ ;variable
  DEF_STRETCH=DEF_STRETCH,$ ;variable
  DESCRIP=DESCRIP,$ ;variable
  DIMS=DIMS,$ ;integer
  FILE_TYPE=FILE_TYPE,$ ;variable
  FNAME=FNAME,$ ;variable
  FUNC_COMPLEX=FUNC_COMPLEX,$ ;{0|1|2|3|4}
  FWHM=FWHM,$ ;variable
  H_INFO=H_INFO,$ ;variable
  INTERLEAVE=INTERLEAVE,$ ;variable
  LOOKUP=LOOKUP,$ ;variable
  LUT_NAME=LUT_NAME,$ ;variable
  NB=NB,$ ;variable
  NL=NL,$ ;variable
  NS=NS,$ ;variable
  NUM_CLASSES=NUM_CLASSES,$ ;variable
  OFFSET=OFFSET,$ ;variable
  READ_PROCEDURE=READ_PROCEDURE,$ ;variable
  REFLECTANCE_SCALE_FACTOR=REFLECTANCE_SCALE_FACTOR,$ ;variable
  SENSOR_TYPE=SENSOR_TYPE,$ ;integer
  SNAME=SNAME,$ ;variable
  SOLAR_IRRADIANCE=SOLAR_IRRADIANCE,$ ;variable
  SPEC_NAMES=SPEC_NAMES,$ ;variable
  STA_NAME=STA_NAME,$ ;variable
  SUN_AZIMUTH=SUN_AZIMUTH,$ ;variable
  SUN_ELEVATION=SUN_ELEVATION,$ ;variable
  WAVELENGTH_UNITS=WAVELENGTH_UNITS,$ ;{0|1|2|3|4|5|6}
  WL=WL,$ ;variable
  XSTART=XSTART,$ ;variable
  YSTART=YSTART,$ ;variable
  map_info=map_info,$
  coordinate_system_string=coordinate_system_string
  H_INFO=dictionary('ACQUISITION_TIME',!NULL,$ ;variable
    'BBL',!NULL,$ ;array
    'BNAMES',!NULL,$ ;variable
    'BYTE_SWAP',!NULL,$ ;variable
    'CLASS_NAMES',!NULL,$ ;variable
    'CLOUD_COVER',!NULL,$ ;variable
    'DATA_GAINS',!NULL,$ ;variable
    'DATA_IGNORE_VALUE',!NULL,$ ;variable
    'DATA_OFFSETS',!NULL,$ ;variable
    'DATA_TYPE',!NULL,$ ;integer
    'DEF_BANDS',!NULL,$ ;array
    'DEF_ZRANGE',!NULL,$ ;variable
    'DEF_STRETCH',!NULL,$ ;variable
    'DESCRIP',!NULL,$ ;variable
    'DIMS',!NULL,$ ;integer
    'FILE_TYPE',!NULL,$ ;variable
    'FNAME',!NULL,$ ;variable
    'FUNC_COMPLEX',!NULL,$ ;{0|1|2|3|4}
    'FWHM',!NULL,$ ;variable
    'INTERLEAVE',!NULL,$ ;variable
    'LOOKUP',!NULL,$ ;variable
    'LUT_NAME',!NULL,$ ;variable
    'NB',!NULL,$ ;variable
    'NL',!NULL,$ ;variable
    'NS',!NULL,$ ;variable
    'NUM_CLASSES',!NULL,$ ;variable
    'OFFSET',!NULL,$  ;variable
    'READ_PROCEDURE',!NULL,$ ;variable
    'REFLECTANCE_SCALE_FACTOR',!NULL,$ ;variable
    'SENSOR_TYPE',!NULL,$ ;integer
    'SNAME',!NULL,$ ;variable
    'SOLAR_IRRADIANCE',!NULL,$ ;variable
    'SPEC_NAMES',!NULL,$ ;variable
    'STA_NAME',!NULL,$ ;variable
    'SUN_AZIMUTH',!NULL,$ ;variable
    'SUN_ELEVATION',!NULL,$ ;variable
    'WAVELENGTH_UNITS',!NULL,$ ;{0|1|2|3|4|5|6}
    'WL',!NULL,$ ;variable
    'XSTART',!NULL,$ ;variable
    'YSTART',!NULL,$ ;variable
    'map_info',!NULL,$
    'coordinate_system_string',!NULL,$
    'fname_hdr',!NULL,$
    'fname_data',!NULL)
  ;解析头文件名称
  n_part=n_elements(strsplit(file_name,'.',/extract))
  len_last=strlen((strsplit(file_name,'.',/extract))[n_part-1])
  file_pos=strmid(file_name,0,strlen(file_name)-len_last-1)
  fhdr=file_pos+'.hdr'
  is_exist=FILE_TEST(fhdr)
  ;读取头文件，并给头文件各个属性赋值
  if is_exist then begin
    openr,lun,fhdr,/get_lun
    ls=''
    readf,lun,ls
    if ls eq '' or strupcase(strtrim(ls,2)) ne 'ENVI' then begin
      dm=dialog_message('不是ENVI格式文件',title='错误',/ERROR)
      H_INFO=!NULL
      return,0
    endif else begin
      lines=list()
      while ~eof(lun) do begin
        readf,lun,ls
        if strpos(ls,'=') gt 0 then begin
          lines.add,ls
        endif else begin
          lines[lines.count()-1]=lines[lines.count()-1]+ls
        endelse
      endwhile
      free_lun,lun
      for i=0,lines.count()-1 do begin
        strs=strsplit(lines[i],'=',/extract)
        case strtrim(strs[0],2) of
          'acquisition time':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            ACQUISITION_TIME=tmp[0]
            delvar,tmp,strss
          end
          'bbl':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strsplit(tmp,',',/extract)
            tmp=strtrim(tmp,2)
            BBL=fix(tmp)
            delvar,tmp,strss
          end
          'band names':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strsplit(tmp,',',/extract)
            tmp=strtrim(tmp,2)
            BNAMES=tmp
            delvar,tmp,strss
          end
          'byte order':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            BYTE_SWAP=fix(tmp[0])
            delvar,tmp,strss
          end
          'class names':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            tmp=strsplit(tmp,',',/extract)
            CLASS_NAMES=tmp
            delvar,tmp,strss
          end
          'cloud cover':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            CLOUD_COVER=fix(tmp[0])
            delvar,tmp,strss
          end
          'data gain values':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            tmp=strsplit(tmp,',',/extract)
            DATA_GAINS=float(tmp)
            delvar,tmp,strss
          end
          'data ignore value':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            DATA_IGNORE_VALUE=fix(tmp[0])
            delvar,tmp,strss
          end
          'data offset values':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            tmp=strsplit(tmp,',',/extract)
            DATA_OFFSETS=float(tmp)
            delvar,tmp,strss
          end
          'data type':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            DATA_TYPE=fix(tmp[0])
            delvar,tmp,strss
          end
          'default bands':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            tmp=strsplit(tmp,',',/extract)
            DEF_BANDS=fix(tmp)
            delvar,tmp,strss
          end
          'z plot range':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            tmp=strsplit(tmp,',',/extract)
            DEF_ZRANGE=float(tmp)
            delvar,tmp,strss
          end
          'default stretch':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            DEF_STRETCH=tmp[0]
            delvar,tmp,strss
          end
          'description':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            DESCRIP=tmp[0]
            delvar,tmp,strss
          end
          'dimensions':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            DIMS=fix(tmp[0])
            delvar,tmp,strss
          end
          'file type':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            FILE_TYPE=tmp[0]
            delvar,tmp,strss
          end
          'fname':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            FNAME=tmp[0]
            delvar,tmp,strss
          end
          'function complex':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            FUNC_COMPLEX=fix(tmp[0])
            delvar,tmp,strss
          end
          'fwhm':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            tmp=strsplit(tmp,',',/extract)
            FWHM=float(tmp)
            delvar,tmp,strss
          end
          'interleave':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            INTERLEAVE=where(['bsq',$
              'bil',$
              'bip'] $
              eq tmp[0])
            delvar,tmp,strss
          end
          'class lookup':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            tmp=strsplit(tmp,',',/extract)
            LOOKUP=tmp
            delvar,tmp,strss
          end
          'lut_name':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            LUT_NAME=tmp[0]
            delvar,tmp,strss
          end
          'bands':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            NB=long(tmp[0])
            delvar,tmp,strss
          end
          'lines':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            NL=long(tmp[0])
            delvar,tmp,strss
          end
          'samples':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            NS=fix(tmp[0])
            delvar,tmp,strss
          end
          'classes':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            NUM_CLASSES=fix(tmp[0])
            delvar,tmp,strss
          end
          'header offset':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            OFFSET=fix(tmp[0])
            delvar,tmp,strss
          end
          'read procedure':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            READ_PROCEDURE=tmp[0]
            delvar,tmp,strss
          end
          'reflectance scale factor':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            REFLECTANCE_SCALE_FACTOR=float(tmp[0])
            delvar,tmp,strss
          end
          'sensor type':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            SENSOR_TYPE=tmp[0]
            delvar,tmp,strss
          end
          'sname':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            SNAME=tmp[0]
            delvar,tmp,strss
          end
          'solar irradiance':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            SOLAR_IRRADIANCE=float(tmp[0])
            delvar,tmp,strss
          end
          'spectrum names':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            tmp=strsplit(tmp,',',/extract)
            SPEC_NAMES=tmp
            delvar,tmp,strss
          end
          'sta_name':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            STA_NAME=tmp[0]
            delvar,tmp,strss
          end
          'sun azimuth':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            SUN_AZIMUTH=float(tmp[0])
            delvar,tmp,strss
          end
          'sun elevation':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            SUN_ELEVATION=float(tmp[0])
            delvar,tmp,strss
          end
          'wavelength units':begin ;{0|1|2|3|4|5|6}
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            WAVELENGTH_UNITS=where(['Micrometers',$
              'Nanometers',$
              'Wavenumber',$
              'GHz',$
              'MHz',$
              'Index',$
              'Unknown'] $
              eq tmp[0])
            delvar,tmp,strss
          end
          'wavelength':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            tmp=strsplit(tmp,',',/extract)
            WL=tmp
            delvar,tmp,strss
          end
          'x start':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            XSTART=fix(tmp[0])
            delvar,tmp,strss
          end
          'y start':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            YSTART=fix(tmp[0])
            delvar,tmp,strss
          end
          'map info':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            map_info=tmp[0]
            delvar,tmp,strss
          end
          'coordinate system string':begin
            strss=strs[1]
            if n_elements(strs) gt 2 then begin
              for k=2,n_elements(strs)-1 do begin
                strss=strss+'='+strs[k]
              endfor
            endif
            tmp=strss
            tmp=strsplit(strtrim(tmp,2),'{',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[1] : tmp
            tmp=strsplit(tmp,'}',/extract)
            tmp=n_elements(tmp) gt 1 ? tmp[0] : tmp
            tmp=strtrim(tmp,2)
            coordinate_system_string=tmp[0]
            delvar,tmp,strss
          end
          else:begin

          end
        endcase
      endfor
      H_INFO.ACQUISITION_TIME=ACQUISITION_TIME ;variable
      H_INFO.BBL=BBL ;array
      H_INFO.BNAMES=BNAMES ;variable
      H_INFO.BYTE_SWAP=BYTE_SWAP ;variable
      H_INFO.CLASS_NAMES=CLASS_NAMES ;variable
      H_INFO.CLOUD_COVER=CLOUD_COVER ;variable
      H_INFO.DATA_GAINS=DATA_GAINS ;variable
      H_INFO.DATA_IGNORE_VALUE=DATA_IGNORE_VALUE ;variable
      H_INFO.DATA_OFFSETS=DATA_OFFSETS ;variable
      H_INFO.DATA_TYPE=DATA_TYPE ;integer
      H_INFO.DEF_BANDS=DEF_BANDS ;array
      H_INFO.DEF_ZRANGE=DEF_ZRANGE ;variable
      H_INFO.DEF_STRETCH=DEF_STRETCH ;variable
      H_INFO.DESCRIP=DESCRIP ;variable
      H_INFO.DIMS=DIMS ;integer
      H_INFO.FILE_TYPE=FILE_TYPE ;variable
      H_INFO.FNAME=FNAME ;variable
      H_INFO.FUNC_COMPLEX=FUNC_COMPLEX ;{0|1|2|3|4}
      H_INFO.FWHM=FWHM ;variable
      H_INFO.INTERLEAVE=INTERLEAVE ;variable
      H_INFO.LOOKUP=LOOKUP ;variable
      H_INFO.LUT_NAME=LUT_NAME ;variable
      H_INFO.NB=NB ;variable
      H_INFO.NL=NL ;variable
      H_INFO.NS=NS ;variable
      H_INFO.NUM_CLASSES=NUM_CLASSES ;variable
      H_INFO.OFFSET=OFFSET ;variable
      H_INFO.READ_PROCEDURE=READ_PROCEDURE ;variable
      H_INFO.REFLECTANCE_SCALE_FACTOR=REFLECTANCE_SCALE_FACTOR ;variable
      H_INFO.SENSOR_TYPE=SENSOR_TYPE ;integer
      H_INFO.SNAME=SNAME ;variable
      H_INFO.SOLAR_IRRADIANCE=SOLAR_IRRADIANCE ;variable
      H_INFO.SPEC_NAMES=SPEC_NAMES ;variable
      H_INFO.STA_NAME=STA_NAME ;variable
      H_INFO.SUN_AZIMUTH=SUN_AZIMUTH ;variable
      H_INFO.SUN_ELEVATION=SUN_ELEVATION ;variable
      H_INFO.WAVELENGTH_UNITS=WAVELENGTH_UNITS ;{0|1|2|3|4|5|6}
      H_INFO.WL=WL ;variable
      H_INFO.XSTART=XSTART ;variable
      H_INFO.YSTART=YSTART ;variable
      H_INFO.map_info=map_info
      H_INFO.coordinate_system_string=coordinate_system_string
      H_INFO.fname_hdr=fhdr
      H_INFO.fname_data=file_name
      return,1
    endelse
  endif else begin
    dm=dialog_message('文件 '+file_name+' 的头文件不存在',title='错误',/ERROR)
    H_INFO=!NULL
    return,0
  endelse
end
pro hdr_file_setup,file_name,$ ;stirng
  ACQUISITION_TIME=ACQUISITION_TIME,$ ;variable
  BBL=BBL,$ ;array
  BNAMES=BNAMES,$ ;variable
  BYTE_SWAP=BYTE_SWAP,$ ;variable  BYTE_ORDER
  CLASS_NAMES=CLASS_NAMES,$ ;variable
  CLOUD_COVER=CLOUD_COVER,$ ;variable
  DATA_GAINS=DATA_GAINS,$ ;variable
  DATA_IGNORE_VALUE=DATA_IGNORE_VALUE,$ ;variable
  DATA_OFFSETS=DATA_OFFSETS,$ ;variable
  DATA_TYPE=DATA_TYPE,$ ;integer
  DEF_BANDS=DEF_BANDS,$ ;array
  DEF_ZRANGE=DEF_ZRANGE,$ ;variable   ZRANGE
  DEF_STRETCH=DEF_STRETCH,$ ;variable
  DESCRIP=DESCRIP,$ ;variable
  DIMS=DIMS,$ ;integer---
  FILE_TYPE=FILE_TYPE,$ ;variable
  FNAME=FNAME,$ ;variable
  FUNC_COMPLEX=FUNC_COMPLEX,$ ;{0|1|2|3|4}
  FWHM=FWHM,$ ;variable
  H_INFO=H_INFO,$ ;variable     INFO
  INTERLEAVE=INTERLEAVE,$ ;variable
  LOOKUP=LOOKUP,$ ;variable
  LUT_NAME=LUT_NAME,$ ;variable---
  NB=NB,$ ;variable
  NL=NL,$ ;variable
  NS=NS,$ ;variable
  NUM_CLASSES=NUM_CLASSES,$ ;variable
  OFFSET=OFFSET,$ ;variable
  READ_PROCEDURE=READ_PROCEDURE,$ ;variable
  REFLECTANCE_SCALE_FACTOR=REFLECTANCE_SCALE_FACTOR,$ ;variable
  SENSOR_TYPE=SENSOR_TYPE,$ ;integer
  SNAME=SNAME,$ ;variable---
  SOLAR_IRRADIANCE=SOLAR_IRRADIANCE,$ ;variable
  SPEC_NAMES=SPEC_NAMES,$ ;variable
  STA_NAME=STA_NAME,$ ;variable---
  SUN_AZIMUTH=SUN_AZIMUTH,$ ;variable
  SUN_ELEVATION=SUN_ELEVATION,$ ;variable
  WAVELENGTH_UNITS=WAVELENGTH_UNITS,$ ;{0|1|2|3|4|5|6}
  WL=WL,$ ;variable
  XSTART=XSTART,$ ;variable
  YSTART=YSTART,$ ;variable
  map_info=map_info,$
  coordinate_system_string=coordinate_system_string
  if keyword_set(H_INFO) then begin
    help,H_INFO,OUTPUT=hlp
    if (strsplit(hlp,/extract))[1] eq 'DICTIONARY' and H_INFO.count() eq 44 then begin
      ACQUISITION_TIME=H_INFO.ACQUISITION_TIME ;variable
      BBL=H_INFO.BBL ;array
      BNAMES=H_INFO.BNAMES ;variable
      BYTE_SWAP=H_INFO.BYTE_SWAP ;variable  BYTE_ORDER
      CLASS_NAMES=H_INFO.CLASS_NAMES ;variable
      CLOUD_COVER=H_INFO.CLOUD_COVER ;variable
      DATA_GAINS=H_INFO.DATA_GAINS ;variable
      DATA_IGNORE_VALUE=H_INFO.DATA_IGNORE_VALUE ;variable
      DATA_OFFSETS=H_INFO.DATA_OFFSETS ;variable
      DATA_TYPE=H_INFO.DATA_TYPE ;integer
      DEF_BANDS=H_INFO.DEF_BANDS ;array
      DEF_ZRANGE=H_INFO.DEF_ZRANGE ;variable   ZRANGE
      DEF_STRETCH=H_INFO.DEF_STRETCH ;variable
      DESCRIP=H_INFO.DESCRIP ;variable
      DIMS=H_INFO.DIMS ;integer---
      FILE_TYPE=H_INFO.FILE_TYPE ;variable
      FNAME=H_INFO.FNAME ;variable
      FUNC_COMPLEX=H_INFO.FUNC_COMPLEX ;{0|1|2|3|4}
      FWHM=H_INFO.FWHM ;variable
      INTERLEAVE=H_INFO.INTERLEAVE ;variable
      LOOKUP=H_INFO.LOOKUP ;variable
      LUT_NAME=H_INFO.LUT_NAME ;variable---
      NB=H_INFO.NB ;variable
      NL=H_INFO.NL ;variable
      NS=H_INFO.NS ;variable
      NUM_CLASSES=H_INFO.NUM_CLASSES ;variable
      OFFSET=H_INFO.OFFSET ;variable
      READ_PROCEDURE=H_INFO.READ_PROCEDURE ;variable
      REFLECTANCE_SCALE_FACTOR=H_INFO.REFLECTANCE_SCALE_FACTOR ;variable
      SENSOR_TYPE=H_INFO.SENSOR_TYPE ;integer
      SNAME=H_INFO.SNAME ;variable---
      SOLAR_IRRADIANCE=H_INFO.SOLAR_IRRADIANCE ;variable
      SPEC_NAMES=H_INFO.SPEC_NAMES ;variable
      STA_NAME=H_INFO.STA_NAME ;variable---
      SUN_AZIMUTH=H_INFO.SUN_AZIMUTH ;variable
      SUN_ELEVATION=H_INFO.SUN_ELEVATION ;variable
      WAVELENGTH_UNITS=H_INFO.WAVELENGTH_UNITS ;{0|1|2|3|4|5|6}
      WL=H_INFO.WL ;variable
      XSTART=H_INFO.XSTART ;variable
      YSTART=H_INFO.YSTART ;variable
      map_info=H_INFO.map_info
      coordinate_system_string=H_INFO.coordinate_system_string
    endif
  endif
  n_part=n_elements(strsplit(file_name,'.',/extract))
  len_last=strlen((strsplit(file_name,'.',/extract))[n_part-1])
  file_pos=strmid(file_name,0,strlen(file_name)-len_last-1)
  fhdr=file_pos+'.hdr'
  openw,lun_hdr,fhdr,/get_lun
  printf,lun_hdr,'ENVI'
  printf,lun_hdr,'description = {'
  if keyword_set(DESCRIP) and DESCRIP ne !NULL then begin
    printf,lun_hdr,'  '+DESCRIP+'}'
  endif else begin
    printf,lun_hdr,'  File Imported into ENVI ['+systime()+']}'
  endelse
  if keyword_set(NS) and NS ne !NULL then begin
    printf,lun_hdr,'samples = '+strtrim(string(NS),2)
  endif else begin
    printf,lun_hdr,'samples = 1'
  endelse
  if keyword_set(NL) and NL ne !NULL then begin
    printf,lun_hdr,'lines   = '+strtrim(string(NL),2)
  endif else begin
    printf,lun_hdr,'lines   = 1'
  endelse
  if keyword_set(NB) and NB ne !NULL then begin
    printf,lun_hdr,'bands   = '+strtrim(string(NB),2)
  endif else begin
    printf,lun_hdr,'bands   = 1'
  endelse
  if keyword_set(OFFSET) and OFFSET ne !NULL then begin
    printf,lun_hdr,'header offset = '+strtrim(string(OFFSET),2)
  endif else begin
    printf,lun_hdr,'header offset = 0'
  endelse
  if keyword_set(FILE_TYPE) and FILE_TYPE ne !NULL then begin
    printf,lun_hdr,'file type = '+strtrim(string(FILE_TYPE),2)
  endif else begin
    printf,lun_hdr,'file type = ENVI Standard'
  endelse
  if keyword_set(DATA_TYPE) and DATA_TYPE ne !NULL then begin
    printf,lun_hdr,'data type = '+strtrim(string(DATA_TYPE),2)
  endif else begin
    printf,lun_hdr,'data type = 1'
  endelse
  if keyword_set(INTERLEAVE) and INTERLEAVE ne !NULL then begin
    printf,lun_hdr,'interleave = '+(['bsq','bil','bip'])[INTERLEAVE]
  endif else begin
    printf,lun_hdr,'interleave = bsq'
  endelse
  if keyword_set(SENSOR_TYPE) and SENSOR_TYPE ne !NULL then begin
    printf,lun_hdr,'sensor type = '+SENSOR_TYPE
  endif else begin
    printf,lun_hdr,'sensor type = Unknown'
  endelse
  if keyword_set(BYTE_SWAP) and BYTE_SWAP ne !NULL then begin
    printf,lun_hdr,'byte order = '+strtrim(string(BYTE_SWAP),2)
  endif else begin
    printf,lun_hdr,'byte order = 0'
  endelse
  if keyword_set(XSTART) and XSTART ne !NULL then begin
    printf,lun_hdr,'x start = '+strtrim(string(XSTART),2)
  endif
  if keyword_set(YSTART) and YSTART ne !NULL then begin
    printf,lun_hdr,'y start = '+strtrim(string(YSTART),2)
  endif
  if keyword_set(map_info) and map_info ne !NULL then begin
    printf,lun_hdr,'map info = {'+map_info+'}'
  endif
  if keyword_set(coordinate_system_string) and coordinate_system_string ne !NULL then begin
    printf,lun_hdr,'coordinate system string = {'+coordinate_system_string+'}'
  endif
  if keyword_set(DEF_BANDS) and DEF_BANDS ne !NULL then begin
    printf,lun_hdr,'default bands = {'
    ns_line=3
    nb_DEF_BANDS=n_elements(DEF_BANDS)
    nl_DEF_BANDS=nb_DEF_BANDS mod ns_line eq 0 ? nb_DEF_BANDS/ns_line : nb_DEF_BANDS/ns_line+1
    if nl_DEF_BANDS le 1 then begin
      str_DEF_BANDS='  '
      for i=0,nb_DEF_BANDS-1 do begin
        str_DEF_BANDS=str_DEF_BANDS+strtrim(string(DEF_BANDS[i]),2)+','
      endfor
      str_DEF_BANDS=strmid(str_DEF_BANDS,0,strlen(str_DEF_BANDS)-1)
      str_DEF_BANDS=str_DEF_BANDS+'}'
      printf,lun_hdr,str_DEF_BANDS
    endif else begin
      for i=0,nl_DEF_BANDS-1 do begin
        str_DEF_BANDS='  '
        if i ne nl_DEF_BANDS-1 then begin
          for j=0,ns_line-1 do begin
            str_DEF_BANDS=str_DEF_BANDS+strtrim(string(DEF_BANDS[i*ns_line+j]),2)+','
          endfor
        endif else begin
          for k=i*ns_line,nb_DEF_BANDS-1 do begin
            str_DEF_BANDS=str_DEF_BANDS+strtrim(string(DEF_BANDS[k]),2)+','
          endfor
          str_DEF_BANDS=strmid(str_DEF_BANDS,0,strlen(str_DEF_BANDS)-1)
          str_DEF_BANDS=str_DEF_BANDS+'}'
        endelse
        printf,lun_hdr,str_DEF_BANDS
      endfor
    endelse
    delvar,nb_DEF_BANDS,nl_DEF_BANDS,str_DEF_BANDS,ns_line
  endif
  if keyword_set(ACQUISITION_TIME) and ACQUISITION_TIME ne !NULL then begin
    printf,lun_hdr,'acquisition time = '+ACQUISITION_TIME
  endif
  if keyword_set(WAVELENGTH_UNITS) and WAVELENGTH_UNITS ne !NULL then begin
    printf,lun_hdr,'wavelength units = '+(['Micrometers',$
      'Nanometers',$
      'Wavenumber',$
      'GHz',$
      'MHz',$
      'Index',$
      'Unknown'])[WAVELENGTH_UNITS]
  endif else begin
    printf,lun_hdr,'wavelength units = Unknown'
  endelse
  if keyword_set(REFLECTANCE_SCALE_FACTOR) and REFLECTANCE_SCALE_FACTOR ne !NULL then begin
    printf,lun_hdr,'reflectance scale factor = '+strtrim(string(REFLECTANCE_SCALE_FACTOR),2)
  endif
  if keyword_set(DEF_ZRANGE) and DEF_ZRANGE ne !NULL then begin
    printf,lun_hdr,'z plot range = {'+strtrim(string(DEF_ZRANGE[0]),2)+','+strtrim(string(DEF_ZRANGE[1]),2)+'}'
  endif
  if keyword_set(DATA_IGNORE_VALUE) and DATA_IGNORE_VALUE ne !NULL then begin
    printf,lun_hdr,'data ignore value = '+strtrim(string(DATA_IGNORE_VALUE),2)
  endif
  if keyword_set(DEF_STRETCH) and DEF_STRETCH ne !NULL then begin
    printf,lun_hdr,'default stretch = '+DEF_STRETCH
  endif
  if keyword_set(SOLAR_IRRADIANCE) and SOLAR_IRRADIANCE ne !NULL then begin
    printf,lun_hdr,'solar irradiance = '+strtrim(string(SOLAR_IRRADIANCE),2)
  endif
  if keyword_set(CLOUD_COVER) and CLOUD_COVER ne !NULL then begin
    printf,lun_hdr,'cloud cover = '+strtrim(string(CLOUD_COVER),2)
  endif
  if keyword_set(FUNC_COMPLEX) and FUNC_COMPLEX ne !NULL then begin
    printf,lun_hdr,'function complex = '+strtrim(string(FUNC_COMPLEX),2)
  endif
  if keyword_set(NUM_CLASSES) and NUM_CLASSES ne !NULL then begin
    printf,lun_hdr,'classes = '+strtrim(string(NUM_CLASSES),2)
  endif
  if keyword_set(SUN_AZIMUTH) and SUN_AZIMUTH ne !NULL then begin
    printf,lun_hdr,'sun azimuth = '+strtrim(string(SUN_AZIMUTH),2)
  endif
  if keyword_set(SUN_ELEVATION) and SUN_ELEVATION ne !NULL then begin
    printf,lun_hdr,'sun elevation = '+strtrim(string(SUN_ELEVATION),2)
  endif
  if keyword_set(SPEC_NAMES) and SPEC_NAMES ne !NULL then begin
    printf,lun_hdr,'spectrum names = {'
    ns_line=5
    nb_SPEC_NAMES=n_elements(LOOKUP)
    nl_SPEC_NAMES=nb_SPEC_NAMES mod ns_line eq 0 ? nb_SPEC_NAMES/ns_line : nb_SPEC_NAMES/ns_line+1
    if nl_SPEC_NAMES le 1 then begin
      str_SPEC_NAMES='  '
      for i=0,nb_SPEC_NAMES-1 do begin
        str_SPEC_NAMES=str_SPEC_NAMES+strtrim(string(SPEC_NAMES[i]),2)+','
      endfor
      str_SPEC_NAMES=strmid(str_SPEC_NAMES,0,strlen(str_SPEC_NAMES)-1)
      str_SPEC_NAMES=str_SPEC_NAMES+'}'
      printf,lun_hdr,str_SPEC_NAMES
    endif else begin
      for i=0,nl_SPEC_NAMES-1 do begin
        str_SPEC_NAMES='  '
        if i ne nl_SPEC_NAMES-1 then begin
          for j=0,ns_line-1 do begin
            str_SPEC_NAMES=str_SPEC_NAMES+strtrim(string(SPEC_NAMES[i*ns_line+j]),2)+','
          endfor
        endif else begin
          for k=i*ns_line,nb_SPEC_NAMES-1 do begin
            str_SPEC_NAMES=str_SPEC_NAMES+strtrim(string(SPEC_NAMES[k]),2)+','
          endfor
          str_SPEC_NAMES=strmid(str_SPEC_NAMES,0,strlen(str_SPEC_NAMES)-1)
          str_SPEC_NAMES=str_SPEC_NAMES+'}'
        endelse
        printf,lun_hdr,str_SPEC_NAMES
      endfor
    endelse
    delvar,nb_SPEC_NAMES,nl_SPEC_NAMES,str_SPEC_NAMES,ns_line
  endif
  if keyword_set(LOOKUP) and LOOKUP ne !NULL then begin
    printf,lun_hdr,'class lookup = {'
    ns_line=20
    nb_LOOKUP=n_elements(LOOKUP)
    nl_LOOKUP=nb_LOOKUP mod ns_line eq 0 ? nb_LOOKUP/ns_line : nb_LOOKUP/ns_line+1
    if nl_LOOKUP le 1 then begin
      str_LOOKUP='  '
      for i=0,nb_LOOKUP-1 do begin
        str_LOOKUP=str_LOOKUP+strtrim(string(LOOKUP[i]),2)+','
      endfor
      str_LOOKUP=strmid(str_LOOKUP,0,strlen(str_LOOKUP)-1)
      str_LOOKUP=str_LOOKUP+'}'
      printf,lun_hdr,str_LOOKUP
    endif else begin
      for i=0,nl_LOOKUP-1 do begin
        str_LOOKUP='  '
        if i ne nl_LOOKUP-1 then begin
          for j=0,ns_line-1 do begin
            str_LOOKUP=str_LOOKUP+strtrim(string(LOOKUP[i*ns_line+j]),2)+','
          endfor
        endif else begin
          for k=i*ns_line,nb_LOOKUP-1 do begin
            str_LOOKUP=str_LOOKUP+strtrim(string(LOOKUP[k]),2)+','
          endfor
          str_LOOKUP=strmid(str_LOOKUP,0,strlen(str_LOOKUP)-1)
          str_LOOKUP=str_LOOKUP+'}'
        endelse
        printf,lun_hdr,str_LOOKUP
      endfor
    endelse
    delvar,nb_LOOKUP,nl_LOOKUP,str_LOOKUP,ns_line
  endif
  if keyword_set(READ_PROCEDURE) and READ_PROCEDURE ne !NULL then begin
    printf,lun_hdr,'read procedure = {'
    ns_line=20
    nb_READ_PROCEDURE=n_elements(READ_PROCEDURE)
    nl_READ_PROCEDURE=nb_READ_PROCEDURE mod ns_line eq 0 ? nb_READ_PROCEDUREP/ns_line : nb_READ_PROCEDURE/ns_line+1
    if nl_READ_PROCEDURE le 1 then begin
      str_READ_PROCEDURE='  '
      for i=0,nb_READ_PROCEDURE-1 do begin
        str_READ_PROCEDURE=str_READ_PROCEDURE+strtrim(string(READ_PROCEDURE[i]),2)+','
      endfor
      str_READ_PROCEDURE=strmid(str_READ_PROCEDURE,0,strlen(str_READ_PROCEDURE)-1)
      str_READ_PROCEDURE=str_READ_PROCEDURE+'}'
      printf,lun_hdr,str_READ_PROCEDURE
    endif else begin
      for i=0,nl_READ_PROCEDURE-1 do begin
        str_READ_PROCEDURE='  '
        if i ne nl_READ_PROCEDURE-1 then begin
          for j=0,ns_line-1 do begin
            str_READ_PROCEDURE=str_READ_PROCEDURE+strtrim(string(READ_PROCEDURE[i*ns_line+j]),2)+','
          endfor
        endif else begin
          for k=i*ns_line,nb_READ_PROCEDURE-1 do begin
            str_READ_PROCEDURE=str_READ_PROCEDURE+strtrim(string(READ_PROCEDURE[k]),2)+','
          endfor
          str_READ_PROCEDURE=strmid(str_READ_PROCEDURE,0,strlen(str_READ_PROCEDURE)-1)
          str_READ_PROCEDURE=str_READ_PROCEDURE+'}'
        endelse
        printf,lun_hdr,str_READ_PROCEDURE
      endfor
    endelse
    delvar,nb_READ_PROCEDURE,nl_READ_PROCEDURE,str_READ_PROCEDURE,ns_line
  endif
  if keyword_set(BNAMES) and BNAMES ne !NULL then begin
    printf,lun_hdr,'band names = {'
    ns_line=5
    nb_BNAMES=n_elements(BNAMES)
    nl_BNAMES=nb_BNAMES mod ns_line eq 0 ? nb_BNAMES/ns_line : nb_BNAMES/ns_line+1
    if nl_BNAMES le 1 then begin
      str_BNAMES='  '
      for i=0,nb_BNAMES-1 do begin
        str_BNAMES=str_BNAMES+strtrim(string(BNAMES[i]),2)+','
      endfor
      str_BNAMES=strmid(str_BNAMES,0,strlen(str_BNAMES)-1)
      str_BNAMES=str_BNAMES+'}'
      printf,lun_hdr,str_BNAMES
    endif else begin
      for i=0,nl_BNAMES-1 do begin
        str_BNAMES='  '
        if i ne nl_BNAMES-1 then begin
          for j=0,ns_line-1 do begin
            str_BNAMES=str_BNAMES+strtrim(string(BNAMES[i*ns_line+j]),2)+','
          endfor
        endif else begin
          for k=i*ns_line,nb_BNAMES-1 do begin
            str_BNAMES=str_BNAMES+strtrim(string(BNAMES[k]),2)+','
          endfor
          str_BNAMES=strmid(str_BNAMES,0,strlen(str_BNAMES)-1)
          str_BNAMES=str_BNAMES+'}'
        endelse
        printf,lun_hdr,str_BNAMES
      endfor
    endelse
    delvar,nb_BNAMES,nl_BNAMES,str_BNAMES,ns_line
  endif
  if keyword_set(WL) and WL ne !NULL then begin
    printf,lun_hdr,'wavelength = {'
    ns_line=6
    nb_WL=n_elements(WL)
    nl_WL=nb_WL mod ns_line eq 0 ? nb_WL/ns_line : nb_WL/ns_line+1
    if nl_WL le 1 then begin
      str_WL='  '
      for i=0,nb_WL-1 do begin
        str_WL=str_WL+strtrim(string(WL[i]),2)+','
      endfor
      str_WL=strmid(str_WL,0,strlen(str_WL)-1)
      str_WL=str_WL+'}'
      printf,lun_hdr,str_WL
    endif else begin
      for i=0,nl_WL-1 do begin
        str_WL='  '
        if i ne nl_WL-1 then begin
          for j=0,ns_line-1 do begin
            str_WL=str_WL+strtrim(string(WL[i*ns_line+j]),2)+','
          endfor
        endif else begin
          for k=i*ns_line,nb_WL-1 do begin
            str_WL=str_WL+strtrim(string(WL[k]),2)+','
          endfor
          str_WL=strmid(str_WL,0,strlen(str_WL)-1)
          str_WL=str_WL+'}'
        endelse
        printf,lun_hdr,str_WL
      endfor
    endelse
    delvar,nb_WL,nl_WL,str_WL,ns_line
  endif
  if keyword_set(FWHM) and FWHM ne !NULL then begin
    printf,lun_hdr,'fwhm = {'
    ns_line=6
    nb_FWHM=n_elements(FWHM)
    nl_FWHM=nb_FWHM mod ns_line eq 0 ? nb_FWHM/ns_line : nb_FWHM/ns_line+1
    if nl_FWHM le 1 then begin
      str_FWHM='  '
      for i=0,nb_FWHM-1 do begin
        str_FWHM=str_FWHM+strtrim(string(FWHM[i]),2)+','
      endfor
      str_FWHM=strmid(str_FWHM,0,strlen(str_FWHM)-1)
      str_FWHM=str_FWHM+'}'
      printf,lun_hdr,str_FWHM
    endif else begin
      for i=0,nl_FWHM-1 do begin
        str_FWHM='  '
        if i ne nl_FWHM-1 then begin
          for j=0,ns_line-1 do begin
            str_FWHM=str_FWHM+strtrim(string(FWHM[i*ns_line+j]),2)+','
          endfor
        endif else begin
          for k=i*ns_line,nb_FWHM-1 do begin
            str_FWHM=str_FWHM+strtrim(string(FWHM[k]),2)+','
          endfor
          str_FWHM=strmid(str_FWHM,0,strlen(str_FWHM)-1)
          str_FWHM=str_FWHM+'}'
        endelse
        printf,lun_hdr,str_FWHM
      endfor
    endelse
    delvar,nb_FWHM,nl_FWHM,str_FWHM,ns_line
  endif
  if keyword_set(CLASS_NAMES) and CLASS_NAMES ne !NULL then begin
    printf,lun_hdr,'class names = {'
    ns_line=5
    nb_CLASS_NAMES=n_elements(CLASS_NAMES)
    nl_CLASS_NAMES=nb_CLASS_NAMES mod ns_line eq 0 ? nb_CLASS_NAMES/ns_line : nb_CLASS_NAMES/ns_line+1
    if nl_CLASS_NAMES le 1 then begin
      str_CLASS_NAMES='  '
      for i=0,nb_CLASS_NAMES-1 do begin
        str_CLASS_NAMES=str_CLASS_NAMES+strtrim(string(CLASS_NAMES[i]),2)+','
      endfor
      str_CLASS_NAMES=strmid(str_CLASS_NAMES,strlen(str_CLASS_NAMES)-1)
      str_CLASS_NAMES=str_CLASS_NAMES+'}'
      printf,lun_hdr,str_CLASS_NAMES
    endif else begin
      for i=0,nl_CLASS_NAMES-1 do begin
        str_CLASS_NAMES='  '
        if i ne nl_CLASS_NAMES-1 then begin
          for j=0,ns_line-1 do begin
            str_CLASS_NAMES=str_CLASS_NAMES+strtrim(string(CLASS_NAMES[i*ns_line+j]),2)+','
          endfor
        endif else begin
          for k=i*ns_line,nb_CLASS_NAMES-1 do begin
            str_CLASS_NAMES=str_CLASS_NAMES+strtrim(string(CLASS_NAMES[k]),2)+','
          endfor
          str_CLASS_NAMES=strmid(str_CLASS_NAMES,0,strlen(str_CLASS_NAMES)-1)
          str_CLASS_NAMES=str_CLASS_NAMES+'}'
        endelse
        printf,lun_hdr,str_CLASS_NAMES
      endfor
    endelse
    delvar,nb_CLASS_NAMES,nl_CLASS_NAMES,str_CLASS_NAMES,ns_line
  endif
  if keyword_set(BBL) and BBL ne !NULL then begin
    printf,lun_hdr,'bbl  = {'
    ns_line=26
    nb_BBL=n_elements(BBL)
    nl_BBL=nb_BBL mod ns_line eq 0 ? nb_BBL/ns_line : nb_BBL/ns_line+1
    if nl_BBL le 1 then begin
      str_BBL='  '
      for i=0,nb_BBL-1 do begin
        str_BBL=str_BBL+strtrim(string(BBL[i]),2)+','
      endfor
      str_BBL=strmid(str_BBL,strlen(str_BBL)-1)
      str_BBL=str_BBL+'}'
      printf,lun_hdr,str_BBL
    endif else begin
      for i=0,nl_BBL-1 do begin
        str_BBL='  '
        if i ne nl_BBL-1 then begin
          for j=0,ns_line-1 do begin
            str_BBL=str_BBL+strtrim(string(BBL[i*ns_line+j]),2)+','
          endfor
        endif else begin
          for k=i*ns_line,nb_BBL-1 do begin
            str_BBL=str_BBL+strtrim(string(BBL[k]),2)+','
          endfor
          str_BBL=strmid(str_BBL,0,strlen(str_BBL)-1)
          str_BBL=str_BBL+'}'
        endelse
        printf,lun_hdr,str_BBL
      endfor
    endelse
    delvar,nb_BBL,nl_BBL,str_BBL,ns_line
  endif
  if keyword_set(DATA_GAINS) and DATA_GAINS ne !NULL then begin
    printf,lun_hdr,'data gain values = {'
    ns_line=4
    nb_DATA_GAINS=n_elements(DATA_GAINS)
    nl_DATA_GAINS=nb_DATA_GAINS mod ns_line eq 0 ? nb_DATA_GAINS/ns_line : nb_DATA_GAINS/ns_line+1
    if nl_DATA_GAINS le 1 then begin
      str_DATA_GAINS='  '
      for i=0,nb_DATA_GAINS-1 do begin
        str_DATA_GAINS=str_DATA_GAINS+strtrim(string(DATA_GAINS[i]),2)+','
      endfor
      str_DATA_GAINS=strmid(str_DATA_GAINS,strlen(str_DATA_GAINS)-1)
      str_DATA_GAINS=str_DATA_GAINS+'}'
      printf,lun_hdr,str_DATA_GAINS
    endif else begin
      for i=0,nl_DATA_GAINS-1 do begin
        str_DATA_GAINS='  '
        if i ne nl_DATA_GAINS-1 then begin
          for j=0,ns_line-1 do begin
            str_DATA_GAINS=str_DATA_GAINS+strtrim(string(DATA_GAINS[i*ns_line+j]),2)+','
          endfor
        endif else begin
          for k=i*ns_line,nb_DATA_GAINS-1 do begin
            str_DATA_GAINS=str_DATA_GAINS+strtrim(string(DATA_GAINS[k]),2)+','
          endfor
          str_DATA_GAINS=strmid(str_DATA_GAINS,0,strlen(str_DATA_GAINS)-1)
          str_DATA_GAINS=str_DATA_GAINS+'}'
        endelse
        printf,lun_hdr,str_DATA_GAINS
      endfor
    endelse
    delvar,nb_DATA_GAINS,nl_DATA_GAINS,str_DATA_GAINS,ns_line
  endif
  if keyword_set(DATA_OFFSETS) and DATA_OFFSETS ne !NULL then begin
    printf,lun_hdr,'data offset values = {'
    ns_line=4
    nb_DATA_OFFSETS=n_elements(DATA_OFFSETS)
    nl_DATA_OFFSETS=nb_DATA_OFFSETS mod ns_line eq 0 ? nb_DATA_OFFSETS/ns_line : nb_DATA_OFFSETS/ns_line+1
    if nl_DATA_OFFSETS le 1 then begin
      str_DATA_OFFSETS='  '
      for i=0,nb_DATA_OFFSETS-1 do begin
        str_DATA_OFFSETS=str_DATA_OFFSETS+strtrim(string(DATA_OFFSETS[i]),2)+','
      endfor
      str_DATA_OFFSETS=strmid(str_DATA_OFFSETS,strlen(str_DATA_OFFSETS)-1)
      str_DATA_OFFSETS=str_DATA_OFFSETS+'}'
      printf,lun_hdr,str_DATA_OFFSETS
    endif else begin
      for i=0,nl_DATA_OFFSETS-1 do begin
        str_DATA_OFFSETS='  '
        if i ne nl_DATA_OFFSETS-1 then begin
          for j=0,ns_line-1 do begin
            str_DATA_OFFSETS=str_DATA_OFFSETS+strtrim(string(DATA_OFFSETS[i*ns_line+j]),2)+','
          endfor
        endif else begin
          for k=i*ns_line,nb_DATA_OFFSETS-1 do begin
            str_DATA_OFFSETS=str_DATA_OFFSETS+strtrim(string(DATA_OFFSETS[k]),2)+','
          endfor
          str_DATA_OFFSETS=strmid(str_DATA_OFFSETS,0,strlen(str_DATA_OFFSETS)-1)
          str_DATA_OFFSETS=str_DATA_OFFSETS+'}'
        endelse
        printf,lun_hdr,str_DATA_OFFSETS
      endfor
    endelse
    delvar,nb_DATA_OFFSETS,nl_DATA_OFFSETS,str_DATA_OFFSETS,ns_line
  endif
  free_lun,lun_hdr
end

pro write_hdr_file,lun,data,interleave,data_type,ns,nl,nb, $
  brange=brange,srange=srange,lrange=lrange,$
  is_ptr=is_ptr
  sizeof=indgen(16)
  sizeof[1]=1
  sizeof[2]=2
  sizeof[3]=4
  sizeof[4]=4
  sizeof[5]=8
  sizeof[12]=2
  sizeof[13]=4
  sizeof[14]=8
  b0=0
  b1=nb-1
  if keyword_set(brange) then begin
    if n_elements(brange) eq 2 then begin
      if brange[0] ge 0 and brange[1] lt nb and brange[0] le brange[1] then begin
        b0=brange[0]
        b1=brange[1]
      endif
    endif
  endif
  s0=0
  s1=ns-1
  if keyword_set(srange) then begin
    if n_elements(srange) eq 2 then begin
      if srange[0] ge 0 and srange[1] lt ns and srange[0] le srange[1] then begin
        s0=srange[0]
        s1=srange[1]
      endif
    endif
  endif
  l0=0
  l1=nl-1
  if keyword_set(lrange) then begin
    if n_elements(lrange) eq 2 then begin
      if lrange[0] ge 0 and lrange[1] lt nl and lrange[0] le lrange[1] then begin
        l0=lrange[0]
        l1=lrange[1]
      endif
    endif
  endif
  if interleave eq 1 then begin
    line=mk_array(s1-s0+1,data_type)
    for j=l0,l1 do begin
      for k=b0,b1 do begin
        line=data[*,j-l0,k-b0]
        point_lun,lun,(long64(nb)*long64(ns)*j+k*long64(ns)+s0)*sizeof[data_type]
        writeu,lun,line
      endfor
    endfor
  endif else if interleave eq 2 then begin
    spec=mk_array(b1-b0+1,data_type)
    for j=l0,l1 do begin
      for i=s0,s1 do begin
        spec=data[i-s0,j-l0,*]
        point_lun,lun,((long64(ns)*long64(j)+i)*long64(nb)+long64(b0))*sizeof[data_type]
        writeu,lun,spec
      endfor
    endfor
  endif else begin
    line=mk_array(s1-s0+1,data_type)
    for k=b0,b1 do begin
      for j=l0,l1 do begin
        line=data[*,j-l0,k-b0]
        point_lun,lun,(long64(ns)*long64(nl)*k+j*long64(ns)+s0)*sizeof[data_type]
        writeu,lun,line
      endfor
    endfor
  endelse
end

pro read_raw
  print,'This a set of codes for reading and writing ENVI file.'
end