-- Gkyl ------------------------------------------------------------------------
--[[
   MAKE SURE TO RUN GKYL WITH AS MANY PROCESSES AS SPECIFIED BY
      decompCuts = {...}
   FURTHER DOWN IN THIS SCRIPT.
--]]

local Plasma = require "App.PlasmaOnCartGrid" -- Importing the Gkyl package

-- adding the path to my analytical approximation package
package.path = package.path..";"..os.getenv("SHOCKLIB").."/luaLib/?.lua"
-- NOTE: need to specify the right path depending on mashine.
-- With the definition below, you need to have an environment
-- variable "SHOCKLIB" leading to the Shock library under the Shock-project
-- in SVN [see the instruction notes for more info].

local Approx = require "shApprox" -- Importing the analytical approx pkg

-- Specifying the path to the file containing the apprimation coefficients
local DATA_PATH = os.getenv("SHOCKLIB").."/EXAMPLE/in-data/"
local filename = DATA_PATH.."COEFS-psi_M1-400_tau200_MB.tsv" -- filename with full path
--[[
   NOTE: You can also specify a relative path to your data file, e.g.
   local filename = "../in-data/COEFS-psi_M1-400_tau200_MB.tsv"
--]]

-- -----------------------------------------------------------------------------
--[[
This is where we specify the input paramters of the shock and 
calculate the normalized shock parameters.
The convention is to use _hat on dimensionful variables, to signify 
that they must be normalized (or are normalization parameters).

See my notes on normalization under:
<SVN-directory>/erc/andsunds/Shock-project/normalizations/
for more on this subject.
--]]


me_hat = 1 -- electron mass 
Z_e = -1 -- electron charge
mi_hat =  1836.153 -- ion mass
Z_i =  1 -- ion charge

Te_hat=200.0 -- electron temp, in some arbitrary unit
Ti_hat=1.0   -- ion temp, in the same arbitrary unit

-- The tau_x varaibles are electron-to-ion teperature ratios,
-- weighted with the charge of species x.
tau_i = Z_i*Te_hat/Ti_hat -- def of tau_i in our norm. scheme
tau_e = Z_e -- =Z_e*Te_hat/Te_hat

ni0_hat = 1.0      -- ion density in some arbitary unit
n0_hat=Z_i*ni0_hat -- =\sum_{j} Z_{j}*n_{j,0}, this is the norm. density
ni0=ni0_hat/n0_hat -- normalized ion density

cs_hat = math.sqrt( (Te_hat/n0_hat)*(Z_i^2*ni0_hat/mi_hat) ) -- Sound speed
M=1.400 -- Mach number

-- The zeta_x variables are normalized charge-to-mass-ratios of
-- species x, normalized against Te_hat/cs_hat^2, which is a
-- wheighted average over all the ion charge-to-mass ratios.
-- See my notes on normalization under:
-- SVN/erc/andsunds/Shock-project/normalizations/
-- for more on this subject.
zeta0_hat=cs_hat^2/Te_hat
zeta_i=(Z_i/mi_hat)/zeta0_hat
zeta_e=(Z_e/me_hat)/zeta0_hat

-- The far upstream densities
ni1 = Approx.ni1(ni0, tau_i, zeta_i,M,filename) -- This calcuation is slow,
                                                -- due to numerical integration
ne1 = Z_i*ni1 -- The far upstream electron density is \sum_{j} Z_{j}*n_{i,1}

-- The speed of the shock
uShock = 0 -- in this frame

-- Ish-thermal speeds
vTe=math.sqrt(math.abs(zeta_e/tau_e))
vTi=math.sqrt(math.abs(zeta_i/tau_i))

-- Pre-defined upper and lower limits for the ION velocities
-- (Needed to get dvIon)
vLimIon={-40.0*vTi-uShock, 40.0*vTi-uShock}
nCells=128
dvIon=(vLimIon[2]-vLimIon[1])/nCells
-- NOTE: dvIon is used to smooth out the discontinuity in the distribution
-- function. This is to prevent noise in the density due to the location
-- of the discontinutity with respect to the velocity grid. 


-- domain size and simulation time
lambdaD=1 -- in this normalization
LX = 40*lambdaD

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 1.0,             -- end time
   suggestedDt = 0.01,     -- this is just a first suggestion
   nFrame = 10,            -- number of output frames
   lower = {-LX/2},        -- configuration space lower left
   upper = {LX/2},         -- configuration space upper right
   cells = {256},          -- configuration space cells
   basis = "serendipity",  -- one of "serendipity" or "maximal-order"
   polyOrder = 2,          -- polynomial order
   timeStepper = "rk3",    -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {4},  -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- electrons
   elc = Plasma.VlasovSpecies {
      -- Since we only have a charge-to-mass ratio variable, and no
      -- normalized mass or charge individualy, we have to use zeta. 
      -- For electrons it is a little tricky, since zeta_e~-2000, Gkyl
      -- starts giving a lot of time-step warnings if we simply set
      --    charge = zeta_e, mass = 1,
      -- as we do for ions. Instead, we have to treat electrons a bit
      -- differntly by setting:
      charge = Z_e, mass = Z_e/zeta_e,
      -- velocity space grid
      lower = {-5.0*vTe-uShock},
      upper = {5.0*vTe-uShock},
      cells = {128},
      decompCuts = {1}, -- do not change, no parallelization in velocity space currently
      -- initial conditions
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 -- Maxwell-Boltzmann distributed electrons
	 return Approx.fe_MB (x,v, ne1, zeta_e, M, filename)
      end,
      -- boundary conditions
      bcx = { Plasma.VlasovSpecies.bcCopy, Plasma.VlasovSpecies.bcCopy},
      evolve = true, -- evolve species?
      -- write out density, flow, (total energy, and heat flux moments)
      -- diagnosticMoments = { "M0", "M1i" },
      -- coll = Plasma.VmLBOCollisions { collFreq = nuElc, },
   },
   -- protons
   ion = Plasma.VlasovSpecies {
      -- Since we only have the normalized charge-to-mass ratio, we
      -- can simply force the "mass" to be 1, and then set the
      -- "chrage" to be zeta_i:
      charge = Z_i, mass = Z_i/zeta_i, 
      -- velocity space grid
      lower = {vLimIon[1]},
      upper = {vLimIon[2]},
      cells = {nCells},
      decompCuts = {1}, -- do not change, no parallelization in velocity space currently
      -- initial conditions
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 -- The shock model ion distribution function.
	 return Approx.fi (x,v, ni0,tau_i, zeta_i, M, filename,dvIon)
      end,
      -- boundary conditions
      bcx = {Plasma.VlasovSpecies.bcCopy, Plasma.VlasovSpecies.bcCopy},
      evolve = true, -- evolve species?
      -- write out density, flow, (total energy, and heat flux moments)
      -- diagnosticMoments = { "M0", "M1i" },
      -- coll = Plasma.VmLBOCollisions { collFreq = nuIon, },
   },

   -- field solver
   field = Plasma.MaxwellField {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
	 local x=xn[1]
	 -- The electric fiels as specified by the approximation
	 -- coefficients in filename
	 local Ex = Approx.Ex (x, filename)
	 return Ex, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      -- boundary conditions
      bcx = {Plasma.MaxwellField.bcCopy, Plasma.MaxwellField.bcCopy},
      evolve = true, -- evolve field?
      hasMagneticField = {false},
   },
}
-- run application
plasmaApp:run()
