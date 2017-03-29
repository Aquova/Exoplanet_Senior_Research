;; Austin Bricker, 10 February 2016

;; This program is the first attempt to read a FITS file gathered of a
;; star, and translate that into a meaningful flux curve output.

;; Import a FITS file, and rename it for ease.

multi = 'off'                   ; If 'off', then only one dataset is analyzed.

file0 = '~/idl/Senior_Research/KIC-9579641/kplr009579641-q0.fits'
if multi ne 'off' then begin
  file1 = '~/idl/Senior_Research/KIC-9579641/kplr009579641-q1.fits'
  file2 = '~/idl/Senior_Research/KIC-9579641/kplr009579641-q2.fits'
  file3 = '~/idl/Senior_Research/KIC-9579641/kplr009579641-q3.fits'
  file4 = '~/idl/Senior_Research/KIC-9579641/kplr009579641-q4.fits'
  file5 = '~/idl/Senior_Research/KIC-9579641/kplr009579641-q5.fits'
endif

;; Extracting just one extention from the file
kepler0 = mrdfits(file0,1)
if multi ne 'off' then begin
  kepler1 = mrdfits(file1,1)
  kepler2 = mrdfits(file2,1)
  kepler3 = mrdfits(file3,1)
  kepler4 = mrdfits(file4,1)
  kepler5 = mrdfits(file5,1)
endif

;; Gather time and flux data, and plot.
time0 = kepler0.TIME
flux0 = kepler0.PDCSAP_FLUX

if multi ne 'off' then begin
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
endif
  
;; Ctrl-c-; comments large areas, BTW.

;; raw = plot(time0,flux0,COLOR='red',XTITLE='Time (Days)',YTITLE='Flux (Electrons/s)',$
;;            DIMENSIONS=[900,500])
if multi ne 'off' then begin
  raw = plot(time1,flux1,/OVERPLOT,COLOR='orange')
  raw = plot(time2,flux2,/OVERPLOT,COLOR='yellow')
  raw = plot(time3,flux3,/OVERPLOT,COLOR='green')
  raw = plot(time4,flux4,/OVERPLOT,COLOR='blue')
  raw = plot(time5,flux5,/OVERPLOT,COLOR='purple')
endif

;; To normalize, subtract each dataset by it's median
;; I really should go back and do this in a loop.

mean0 = MEAN(flux0,/NAN)
nor_flux0 = flux0 / mean0 - 1
if multi ne 'off' then begin
  mean1 = MEAN(flux1,/NAN)
  nor_flux1 = flux1 / mean1
  mean2 = MEAN(flux2,/NAN)
  nor_flux2 = flux2 / mean2
  mean3 = MEAN(flux3,/NAN)
  nor_flux3 = flux3 / mean3
  mean4 = MEAN(flux4,/NAN)
  nor_flux4 = flux4 / mean4
  mean5 = MEAN(flux5,/NAN)
  nor_flux5 = flux5 / mean5
endif

;; norplot = plot(time0,nor_flux0,COLOR='red',XTITLE='Time (Days)',YTITLE='Normalized Flux')
if multi ne 'off' then begin
  norplot = plot(time1,nor_flux1,/OVERPLOT,COLOR='orange')
  norplot = plot(time2,nor_flux2,/OVERPLOT,COLOR='yellow')
  norplot = plot(time3,nor_flux3,/OVERPLOT,COLOR='green')
  norplot = plot(time4,nor_flux4,/OVERPLOT,COLOR='blue')
  norplot = plot(time5,nor_flux5,/OVERPLOT,COLOR='purple')
endif

;; Now to find a way to search for the period automatically. 
;; Using a Fourier Transform (or similar procedure) gives an
;; approximate estimate to the period.

;; FFT function doesn't work with any NaN's present.
;; Following function replaces NaN's with value to the left.
;; no_NaN.pro is located in same folder as this program.

;; fixed_nor0 = NO_NAN(nor_flux0)
;; foo = LNP_TEST(time0, fixed_nor0, WK1=freq0 ,WK2=power0, JMAX=i0)
;; period0 = 1/freq0
;; plot,period0,power0

;; Remove noise from data via convolution, then fold over various
;; trial periods, and output the one with deepest dip.

delta = time0[1] - time0[0] ; Timestep of the data

bestFit = fltarr(3) ; Output array 
;; for kernelWidth = 1, 5 do begin
kernelWidth = 3
   for trialPeriod = 5., 6., 0.1 do begin
      print,'kernelWidth = ', kernelWidth
      print,'trialPeriod = ', trialPeriod
      kernel = fltarr(kernelWidth) + 1 ; Generates kernel of kernelWidth length
      conv = convol(nor_flux0,kernel,/NAN) ; First convolve to reduce noise
      modTime = time0 MOD trialPeriod
      bin_scat, modTime, conv, XBINS=x, MEAN=convBinned
      currentMin = min(convBinned,minIndex)
      test = plot(x, convBinned, TITLE=trialPeriod)
      print,'currentMin = ',currentMin
      if currentMin lt bestFit[0] then begin ; If found min is less than stored value, then execute.
         print,'New match!', trialPeriod
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
t_0 = time0[0] ; This is the initial time for the data plot
T = bestFit[1] ;This is the period at which the spike occurs. Currently a constant.
;; 5.4122 is correct period for KIC-9579641

;; phase0 = (time0 - t_0) / T MOD 2
phase0 = time0 MOD T
;; if multi ne 'off' then begin
;;   phase1 = (time1 - t_0) / T MOD 2
;;   phase2 = (time2 - t_0) / T MOD 2
;;   phase3 = (time3 - t_0) / T MOD 2
;;   phase4 = (time4 - t_0) / T MOD 2
;;   phase5 = (time5 - t_0) / T MOD 2
;; endif

plot,phase0,nor_flux0,YRANGE=[-0.002,0.002],PSYM=3
;; if multi ne 'off' then begin
;;   oplot,phase1,nor_flux1,PSYM=3
;;   oplot,phase2,nor_flux2,PSYM=3
;;   oplot,phase3,nor_flux3,PSYM=3
;;   oplot,phase4,nor_flux4,PSYM=3
;;   oplot,phase5,nor_flux5,PSYM=3
;; endif
      
      
end
