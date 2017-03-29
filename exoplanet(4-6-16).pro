;; Austin Bricker, 10 February 2016

;; This program is the first attempt to read a FITS file gathered of a
;; star, and translate that into a meaningful flux curve output.

;; Import a FITS file, and rename it for ease.

;; file0 = '~/idl/Senior_Research/KIC-9579641/kplr009579641-q0.fits'
;; file1 = '~/idl/Senior_Research/KIC-9579641/kplr009579641-q1.fits'
;; file2 = '~/idl/Senior_Research/KIC-9579641/kplr009579641-q2.fits'
;; file3 = '~/idl/Senior_Research/KIC-9579641/kplr009579641-q3.fits'
;; file4 = '~/idl/Senior_Research/KIC-9579641/kplr009579641-q4.fits'
;; file5 = '~/idl/Senior_Research/KIC-9579641/kplr009579641-q5.fits'

file0 = '~/idl/Senior_Research/KIC-11853905/kplr011853905-q0.fits'
file1 = '~/idl/Senior_Research/KIC-11853905/kplr011853905-q1.fits'
file2 = '~/idl/Senior_Research/KIC-11853905/kplr011853905-q2.fits'
file3 = '~/idl/Senior_Research/KIC-11853905/kplr011853905-q3.fits'
file4 = '~/idl/Senior_Research/KIC-11853905/kplr011853905-q4.fits'
file5 = '~/idl/Senior_Research/KIC-11853905/kplr011853905-q5.fits'

;; Extracting just one extention from the file
kepler0 = mrdfits(file0,1)
kepler1 = mrdfits(file1,1)
kepler2 = mrdfits(file2,1)
kepler3 = mrdfits(file3,1)
kepler4 = mrdfits(file4,1)
kepler5 = mrdfits(file5,1)

;; Gather time and flux data, and plot.
time0 = kepler0.TIME
flux0 = kepler0.PDCSAP_FLUX
time1 = kepler1.TIME
flux1 = kepler1.PDCSAP_FLUX 
time2 = kepler2.TIME
flux2 = kepler2.PDCSAP_FLUX
time3 = kepler3.TIME
flux3 = kepler3.PDCSAP_FLUX
time4 = kepler4.TIME
flux4 = kepler4.PDCSAP_FLUX
time5 = kepler5.TIME
flux5 = kepler5.PDCSAP_FLUX

time = [time0,time1,time2,time3,time4,time5]
flux = [flux0,flux1,flux2,flux3,flux4,flux5]
  
;; Ctrl-c-; comments large areas, BTW.

raw = plot(time,flux,COLOR='red',XTITLE='Time (Days)',YTITLE='Flux (Electrons/s)',$
           DIMENSIONS=[900,500])

;; To normalize, subtract each dataset by it's median
;; I really should go back and do this in a loop.


mean0 = MEAN(flux0,/NAN)
mean1 = MEAN(flux1,/NAN)
mean2 = MEAN(flux2,/NAN)
mean3 = MEAN(flux3,/NAN)
mean4 = MEAN(flux4,/NAN)
mean5 = MEAN(flux5,/NAN)

nor_flux0 = flux0 / mean0 - 1
nor_flux1 = flux1 / mean1 - 1
nor_flux2 = flux2 / mean2 - 1
nor_flux3 = flux3 / mean3 - 1
nor_flux4 = flux4 / mean4 - 1
nor_flux5 = flux5 / mean5 - 1

nor_flux = [nor_flux0,nor_flux1,nor_flux2,nor_flux3,nor_flux4,nor_flux5]

norplot = plot(time,nor_flux,COLOR='red',XTITLE='Time (Days)',YTITLE='Normalized Flux')

;; Remove noise from data via convolution, then fold over various
;; trial periods, and output the one with deepest dip.

delta = time0[1] - time0[0] ; Timestep of the data

bestFit = fltarr(3)             ; Output array
kernelMax = 5.
kernelWidth = 3
;; for kernelWidth = 1, kernelMax do begin
   for trialPeriod = 1., 10., 0.001 do begin
      if trialPeriod MOD 1 eq 0 then begin ; This section is entirely for fun.
         status = (kernelWidth - 1) / kernelMax
         print,status * 100,"% completed." ; It prints out how far along the analysis is :)
      endif
      kernel = fltarr(kernelWidth) + 1 ; Generates kernel of kernelWidth length
      conv = convol(nor_flux,kernel,/NAN) ; First convolve to reduce noise
      modTime = time MOD trialPeriod
      bin_scat, modTime, conv, XBINS=x, MEAN=convBinned
      currentMin = min(convBinned,minIndex)
      ;; test = plot(x, convBinned, TITLE=trialPeriod)
      if currentMin lt bestFit[0] then begin ; If found min is less than stored value, then execute.
         bestFit[0] = currentMin
         bestFit[1] = trialPeriod
         bestFit[2] = kernelWidth
      endif
   endfor
;; endfor

print,'Results: '
print,'Min value found as: ', bestFit[0]
print,'Best period found as: ', bestFit[1]
print,'Using a kernel width of: ', bestFit[2]

;; Now to fold the data to expose the decrease in flux.
t_0 = time[0] ; This is the initial time for the data plot
T = bestFit[1] ; This is the period at which the spike occurs
;; 5.4122 is correct period for KIC-9579641

;; phaseX = (time - t_0) / T MOD 2
phase = time MOD T

plot,phase,nor_flux,YRANGE=[-0.002,0.002],PSYM=3
     
end
