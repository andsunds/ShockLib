local P={} -- List of functions to be returned when used as a package

local function erf_h(x)
   --[[ 
      Function from the picomath project (public domain):
      https://github.com/ghewgill/picomath
      Based on formula 7.1.26 in Abramowitz and Stegun.
      This function works best for high values of |x|.
   --]]
   -- constants
   a1 =  0.254829592
   a2 = -0.284496736
   a3 =  1.421413741
   a4 = -1.453152027
   a5 =  1.061405429
   p  =  0.3275911
   
   --x = math.abs(x)
   
   -- A&S formula 7.1.26
   t = 1.0/(1.0 + p*x)
   y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)
   return y
end
local function erf_l(x)
   --[[
      Function based on the Taylor series expansion of erf,
      formula 7.1.5 in Abramowitz and Stegun.
      This function works best for low values of |x|.
   --]]
   -- Constants
   b0 = 1.128379167095513 -- 2/math.sqrt(math.pi)
   b1 = 0.333333333333333 -- 1/3
   b2 = 0.10              -- 1/10
   b3 = 0.023809523809524 -- 1/42
   b4 = 0.004629629629630 -- 1/216
   xsq=x*x; -- pre-calc x^2
   y=b0*x*(1 -xsq*(b1 -xsq*(b2 -xsq*(b3 -b4*xsq))));
   return y
end

function P.erf(x)
   --[[
      An umbrella function that uses two different ways of 
      calcualting the error function, erf(x), depending on 
      the size of |x|. This function is continuous over the 
      cut-off point between the two differnt approximation 
      methods.
      The maximum absolute error is < 1.4e-7, and 
      the maximum relative error is < 2.7e-7.
      Both these error calculations have been made agains the
      built-in MATLAB erf(x).
   --]]
   -- Checks the sign of x
   neg = false
   if x < 0 then
      neg = true
      x=-x -- only a positive x is passed to erf_l and erf_h
   end
   
   -- Determines which function to use; the cut-off point between
   -- them is 0.358... and is choosen so that this implementation
   -- is continuous.
   if x>0.358645716186629 then
      y=erf_h(x)
   else
      y=erf_l(x)
   end
   -- Adjusts for the sign of x
   if neg then y=-y end
   
   return y
end



function P.erfc(x)
   -- Complementary error function
   return 1.000000000-P.erf(x)
end

-- Exporting the functions
return P
