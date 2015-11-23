pro solps_plots

    N = 38
    data = hash()

    ; temperature
    spawn, 'echo "te writ 95 f.y" | b2plot'
    data = data + hash('te',read_b2plot(N))
    
    ; particle flux
    spawn, 'echo "fnay sy m/ writ jxa f.y" | b2plot' 
    data = data + hash('flux',read_b2plot(N))
    
    ; electron heat flux
    spawn, 'echo "fhey sy m/ writ jxa f.y" | b2plot' 
    data = data + hash('eheatflux',read_b2plot(N))
    
    ; ion heat flux
    spawn, 'echo "fhiy sy m/ writ jxa f.y" | b2plot'
    data = data + hash('iheatflux',read_b2plot(N))
    
    ;ion density
    spawn, 'echo "na 1 zsel writ jxa f.y" | b2plot' 
    data = data + hash('idens',read_b2plot(N))
    
    ;atomic hydrogen density
    spawn, 'echo "nh 1 zsel writ jxa f.y" | b2plot' 
    data = data + hash('hdens',read_b2plot(N))

    ;poloidal diamagnetic current due to drift
    spawn, 'echo "cd 1 zsel writ jxa f.y" | b2plot'
    data = data + hash('poldia',read_b2plot(N))
    
    ;radial diamagnetic current due to drift
    spawn, 'echo "dd 1 zsel writ jxa f.y" | b2plot'
    data = data + hash('raddia',read_b2plot(N))
    
    ;molecular density
    spawn, 'echo "dmb2 1 zsel writ jxa f.y" | b2plot'
    data = data + hash('moldens',read_b2plot(N))
    
    ;atomic density
    spawn, 'echo "dab2 1 zsel writ jxa f.y" | b2plot'
    data = data + hash('atomdens',read_b2plot(N))
    
    ;radial atomic flux
    spawn, 'echo "rfla sy m/ writ jxa f.y" | b2plot'
    data = data + hash('atomflux',read_b2plot(N))
    
    ;radial molecular flux
    spawn, 'echo "rflm sy m/ writ jxa f.y" | b2plot'
    data = data + hash('molflux',read_b2plot(N))
  

    ; Write data to file

    FileName = 'ips-solps.nc'
    id = ncdf_create(FileName,/clobber) 
    ncdf_control, id, /fill
    nid = ncdf_dimdef( id, 'N', N)
    
    for i=0,data.count()-1 do begin
        ThisKey = (data.keys())[i]  
        print, 'Adding '+ThisKey+' to solps.nc'
        if i ne 0 then NCDF_CONTROL, id, /redef
        ThisId_r = ncdf_vardef(id,'r_'+ThisKey,[nid],/float)
        ThisId = ncdf_vardef(id,ThisKey,[nid],/float)
        NCDF_CONTROL, id, /ENDEF
        ncdf_varput, id, ThisId_r, ((data[ThisKey])[0,*])[*]
        ncdf_varput, id, ThisId, ((data[ThisKey])[1,*])[*]
    endfor 

    ncdf_close, id
    
end
