-- Lua script for reading in the approximation coefficients.
-- (c) Andréas Sundström, 2018

-- Loading the picomath package
local pmath = require "picomath"

local function getCoefs (filename)
   --[[
      This function reads the files containing the coefficients 
      exported from MATLAB.
      The format of the exported files must be:
      -------------
      |K      Cp_0|
      |Cm_0   Cp_1|
      |Cm_1   Cp_2|
      |...    ... |
      |0      Cp_n|
      ------------
      where K is the downstream oscillation wave number, Cm_k 
      the cosine transformation coefficient, and Cp_k the 
      polynomial coefficients for the numeratro polynomial, p(x),
      in the upstream approximation [see the function "rational"
      for more details]. 
   --]]

   -- opens the file with the coefficients
   local file = assert(io.open(filename, "r")) 
   local coefs = {} -- init list
   local N=1 -- init counter
   while true do
      -- reads the two columns into the coef list
      coefs[N],coefs[N+1]= file:read("*number", "*number")
      if not coefs[N] then break end -- break at end of file
      N=N+2 -- increase counter
   end
   file:close()
   N=N-1 -- N is now the length of coefs
   return coefs, N
end

local function lambda (filename)
   --[[
      This function returns the wavelength of the DS oscillations,
      by simply using the frequency that is given in the coefficient
      file.
   --]]
   local coefs = getCoefs(filename)
   return 2*math.pi/coefs[1]
end

local function rational (x, coefs, N)
   --[[
      Retruns the rational function used to approimate log(psi(x)),
      as well as its derivative, which is used to calcualte E(x).
      The rational function
      r(x) = p(x)/q(x),
      where
      q(x) = x^{n-1}+1, and
      p(x) = \sum_{m=0}^{n} Cp_{m}*x^{m}

      Note that coefs[i] and Cp_{m} does not line up, since coefs 
      begins with index 1 and also contains the CT coefficients.
   --]]
   local n       = 0 -- polynomial powers
   local p, q, r =0.0, 0.0, 0.0 -- init p,q,r
   local dp,dq,dr=0.0, 0.0, 0.0 -- init of their derivatives
   for i=2,N,2 do
      n=i/2-1
      p=p+coefs[i]*x^n
      if n>0 then
	 dp=dp+n*coefs[i]*x^(n-1)
      end
   end
   q=x^(n-1)+1
   dq=(n-1)*x^(n-2)
   r=p/q
   dr=dp/q-r*dq/q

   return r, dr
end

local function psimax (filename)
   --[[
      Returning the value of psimax, as determined by the CT
   --]]
   -- local coefs,N = getCoefs(filename)
   -- local fv=0.0
   -- for i=3,N-1,2 do
   --    n=(i-3)/2
   --    fv=fv+coefs[i]
   -- end
   -- return fv
   local coefs = getCoefs(filename)
   return math.exp(coefs[2])
end

local function psi (x, filename)
   --[[
      Function that calculates the approximation of psi(x), based on
      the coefficients given in the file: filename.
   --]]
   local coefs, N = getCoefs(filename)
   local fv, n=0.0
   if x<0 then
      --Use the CT in the downstream (x<0).
      local K=coefs[1]
      for i=3,N-1,2 do
	 n=(i-3)/2
	 fv=fv+coefs[i]*math.cos(n*K*x)
      end
   else
      --Otherwise use the exp(r(x)) approxiamtion.
      local r=rational(x, coefs,N)
      fv=math.exp(r)
   end
   return fv
end

local function Ex (x, filename)
   --[[
      Function that calculates the approximation of E(x), based on
      the coefficients given in the file: filename.
   --]]
   local coefs, N = getCoefs(filename) 
   local fv,n =0.0
   if x<0 then
      --Use the CT in the downstream (x<0).
      local K=coefs[1]
      for i=3,N-1,2 do
	 n=(i-3)/2
	 fv=fv+coefs[i]*n*K*math.sin(n*K*x)
      end
   else
      --Otherwise use the exp(r(x)) approxiamtion.
      local r,dr=rational(x, coefs,N)
      fv= -math.exp(r)*dr
   end
   return fv
end

local function fi_exp (v,ni0,taui,zetai,psi,M)
   --[[
      Helper function that calcualtes the exponential part of the 
      ion distribution function.
   --]]
   local ex=math.exp(-taui/(2*zetai)*(math.sqrt(v^2+2*zetai*psi)-M)^2)
   return ni0*math.sqrt(taui/(2*math.pi*zetai))*ex
end

local function fi (x,v, ni0,taui, zetai, M, filename, dv)
   --[[
      Calculates the ion distribution function, based on the approximated
      psi(x). This function uses the convention that x=0 is the location 
      of phimax, and hence the boundary between the US and DS. 

      As Gkyl does not hadle exact zero distribution functions, this 
      function implements a rapid exponential decay down to 1e-18 as the 
      "as good as zero" value.
   --]]
   local psi0, psimax0 = psi(x,filename),  psimax(filename)
   local v0
   -- a catch for slightly too big psi, otherwise v0 becomes nan and the
   -- whole resutls becomes nan
   if psi0<psimax0 then
      v0= math.sqrt(2*zetai*(psimax0-psi0))
   else
      v0=0.0
   end
   
   -- print("v0 = ",v0) --DEBUG
   local USDS
   if x<0 then USDS=-1 else USDS=1 end
   -- print("USDS = ",USDS) --DEBUG
   local fv
   if v<=USDS*v0 then
      fv= fi_exp(v,ni0,taui,zetai,psi0,M)
   else
      --[[
      local decay=10 -- Some large number to simulate a fast decay/discontinuity
      local ex=math.exp(-decay*taui*(v-USDS*v0))
      --fv= fi_exp(v0,ni0,taui,zetai,psi0,M)*ex
      fv= math.max(fi_exp(v0,ni0,taui,zetai,psi0,M)*ex, 1e-18)
      --]]
      if dv==nil or dv==nan then
	 dv=zetai/taui
      end
      local ex=math.exp(-10*(v-USDS*v0)^2/dv^2)
      fv= math.max(fi_exp(v0,ni0,taui,zetai,psi0,M)*ex, 1e-18)
   end
   return fv
end
local function ni1 (ni0,taui,zetai,M,filename)
   --[[
      Function that calculates the far upstream (psi=0) ion density,
      by using the erf functions provided in the picomath package.

      NOTE: this has a rather long evaluation time.
   --]]
   local sr_tz=math.sqrt(taui/zetai)
   local sr_zpmax=math.sqrt(2*zetai*psimax(filename))
   return 0.5*ni0*(1+2*pmath.erf(sr_tz*M)+pmath.erf(sr_tz*(sr_zpmax-M)))

end
--[[
A Decpricated version of ni1, using a poor numerical integration
method.
local function ni1 (ni0,taui,zetai,M,filename)
   --[[
      Function that calculates the far upstream (psi=0) ion density.
      This function numerically integrates the ion distribution function
      with psi=0.

      NOTE: this has a rather long evaluation time.
   --] ]
   local TZ = taui/zetai
   local dv = 1e-5/math.sqrt(TZ)
   local v0 = math.sqrt( 2*zetai*psimax(filename) )

   local I=0.0
   for v=-20/math.sqrt(TZ),v0,dv do
      I=I+fi_exp(v,ni0,taui,zetai,0,M)
   end
   return I*dv
end
--]]
local function fe_MB (x,v,ne1,zetae,M,filename)
   --[[
      Electron distribution function, here it is just a shifted
      Maxwellian.
   --]]
   local taue=-1 -- taue=Ze*Te/Te=Ze=-1, by normalization

   zetae=-math.abs(zetae) -- Must be negative
   local Te_Ze=math.abs(taue/zetae) -- Must be positive
   local psi0=psi(x,filename)
   local ex=math.exp( -0.5*Te_Ze*((v+M)^2 + 2*zetae*psi0) )
   return ne1*math.sqrt(Te_Ze/(2*math.pi))*ex
end


-- This part is for exporting the functions
local Approx={} -- This will be returned
-- These functions get exported
Approx.lambda = lambda
Approx.psi = psi
Approx.psimax = psimax
Approx.Ex = Ex
Approx.fi = fi
Approx.ni1= ni1
Approx.fe_MB = fe_MB
return Approx

 
