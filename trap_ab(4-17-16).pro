;; +
;;
;; NAME:
;;      trap_ab
;;
;; PURPOSE:
;;      To plot a trapezoidal waveform.
;;
;; CATEGORY:
;;      ???
;;
;; CALLING SEQUENCE:
;;      Result = trap_ab(x, b1, b2, h, c)
;;
;; INPUTS:
;;      x = Array of x-values to which trapezoidal waveform will be
;;      applied to.
;;      b1 = Length of shorter base of the trapezoid.
;;      b2 = Length of longer base of the trapezoid.
;;      h = Height of the trapezoid, positive if pulse is to be
;;      negative.
;;      c = x-value at which trapezoid is centered.
;;
;; OUTPUTS:
;;      Will output the y-values of a trapezoidal waveform, with the
;;      same number of elements as x. The trapezoid will be of height h,
;;      centered at c, and have bases of b1 and b2. 
;;
;; EXAMPLE:
;;      Given an array, the function will return y-values that
;;      correspond to a trapezoidal waveform. Intended to be used to fit a
;;      curve to an exoplanet flux curve, in conjuction with mpfitfun.pro.
;;
;; MODIFICATION HISTORY:
;;      Written by:     Austin Bricker, 15 April 2016
;;
;; -


function trap_ab, x, b1, b2, h, c
  
  partialStart = c - (b2 / 2)
  totalStart = c - (b1 / 2)
  totalEnd = c + (b1 / 2)
  partialEnd = c + (b2 / 2)

  lslope = h / (partialStart - totalStart)
  rslope = -lslope
  
  nMax = n_elements(x)
  y = fltarr(nMax)
  for i=0,nMax-1 do begin
     if (x[i] LT partialStart) OR (x[i] GE partialEnd) then begin
        y[i] = 1
     endif
     if (x[i] GE partialStart) AND (x[i] LT totalStart) then begin
        y[i] = 1 + lslope * (x[i] - partialStart)
     endif
     if (x[i] GE totalStart) AND (x[i] LT totalEnd) then begin
        y[i] = 1 - h
     endif
     if (x[i] GE totalEnd) AND (x[i] LT partialEnd) then begin
        y[i] = (1 - h) + rslope * (x[i] - totalEnd)
     endif
  endfor
  return, y
end

b1 = 20
b2 = 30
h = 0.1
c = 50

x = findgen(100)
y = trap_ab(x,b1,b2,h,c)
plot,x,y,YRANGE=[-0.1,1.1]

end
