;; +
;;
;; NAME:
;;      no_NaN
;;
;; PURPOSE:
;;      To replace any 'NaN' elements of an array with the element one
;;      index lower.
;;
;; CATEGORY:
;;      ???
;;
;; CALLING SEQUENCE:
;;      Result = no_NaN(arr)
;;
;; INPUTS:
;;      arr = An array.
;;
;; OUTPUTS:
;;      Will output a vector array, which has same number of elements as the input
;;      array, but with NaN instances replaced. 
;;
;; EXAMPLE:
;;      Given some array with sparatic NaN elements,
;;      Result = no_NaN(arr) yields the same array, but with the NaN's
;;      replaced with the element to the left. 
;;
;; MODIFICATION HISTORY:
;;      Written by:     Austin Bricker, 15 March 2016
;;
;; -

function NO_NAN, arr
  n = n_elements(arr)
  resultArr = fltarr(n)
  for i=0,n-1 do begin                ; Go along all elements of the array.
     if arr[i] eq arr[i] then begin   ; NaN's aren't equal to themselves.
        resultArr[i] = arr[i]
     endif else begin                 ; If it isn't equal to itself, it is a NaN.
        resultArr[i] = resultArr[i-1] ; Thus, replace it.
     endelse
  endfor
  return, resultArr                   ; Output the new array.
end


