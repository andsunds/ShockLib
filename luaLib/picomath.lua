local P={} -- List of functions to be returned when used as a package

function P.erf(x)
   --[[ 
      Function from the picomath project (public domain):
      https://github.com/ghewgill/picomath
      Based on formula 7.1.26 in Abramowitz and Stegun, with
      an absolute error <= 2.5e-7      
   --]]
   -- constants
   a1 =  0.254829592
   a2 = -0.284496736
   a3 =  1.421413741
   a4 = -1.453152027
   a5 =  1.061405429
   p  =  0.3275911
   
   -- Checks the sign of x
   neg = false
   if x < 0 then
      neg = true
      x=-x
   end
   --x = math.abs(x)
   
   -- A&S formula 7.1.26
   t = 1.0/(1.0 + p*x)
   y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)

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
