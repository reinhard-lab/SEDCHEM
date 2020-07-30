# **************************************************************************************************************** #
#
# SEDCHEM v0.1
# 
# **************************************************************************************************************** #
#
# SEDCHEM is a reaction-transport model for simulating the redox and acid-base chemistry
# of marine sediments, designed specifically to interrogate the mass flux and isotopic
# consequences of authigenic carbonate formation and benthic phosphorus recycling.
#
# Developed by:
#
#                 Chris Reinhard
#                 School of Earth and Atmospheric Sciences
#                 Georgia Institute of Technology
#           and
#                 Ming-Yu Zhao
#                 Department of Earth and Planetary Sciences
#                 Yale University
#
# **************************************************************************************************************** #
#
# **************************************************************************************************************** #
# LOAD REQUIRED PACKAGES
# **************************************************************************************************************** #
  library(ReacTran)
  library(seacarb)
  library(deSolve)
  library(marelac)
  library(NORMT3)
  library(dplyr)
  library(openxlsx)
# **************************************************************************************************************** #
#
#
# **************************************************************************************************************** #
# LOAD/SET PARAMETERS
# **************************************************************************************************************** #
#
# ---------------------------------------------------------------------------------------------------------------- #
# bottom water properties
# ---------------------------------------------------------------------------------------------------------------- #
#
  BW_Temp       <- 12.0                             # temperature [deg C] 
  BW_Sal        <- 28.4                             # salinity [o/oo]
  BW_Pr         <- 2.0                              # pressure (bar)
  D_sw          <- rho(BW_Sal,BW_Temp,BW_Pr)[[1]]   # density at {T,S,P} (kg/m^3)
  BW_pH         <- 7.64                             # pH 
  BW_H          <- 10^(-BW_pH)                      # mmol/cm^3
  BW_SumCO2     <- 2/1e3                            # [DIC] [mmol/cm^3] 
  BW_SumH2S     <- 0                                # [H2S] [mmol/cm^3]
  BW_Ca         <- 8.6/1e3                          # [Ca] [mmol/cm^3] 
  BW_Mg         <- 46/1e3                           # [Mg] [mmol/cm^3]  
  BW_Na         <- 380.56/1e3                       # [Na] [mmol/cm^3]
  BW_F          <- 0.07/1e3                         # [F] [mmol/cm^3] 
  BW_SumSO4     <- 22/1e3                           # [SO4] [mmol/cm^3]
  BW_O2         <- 0.15/1e3                         # [O2] [mmol/cm^3] 
  BW_NO3        <- 11.8/1e6                         # [NO3] [mmol/cm^3]                
  BW_MnII       <- 0/1e6                            # [Mn] [mmol/cm^3] 
  BW_FeII       <- 0/1e6                            # [Fe] [mmol/cm^3] 
  BW_SumNH4     <- 0/1e6                            # [NH4] [mmol/cm^3]
  BW_CH4        <- 0                                # [CH4] [mmol/cm^3]
  BW_SumH2PO4   <- 1/1e6                            # [PO4] [mmol/cm^3]
#
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# composition of organic C and Fe (oxy)hydroxides
# ---------------------------------------------------------------------------------------------------------------- #
#
  x             <- 10/106                           # molar ratio of N:C in highly reactive organic matter
  xx            <- 10/106                           # molar ratio of N:C in refractory organic matter
  y             <- 1.3/106                          # molar ratio of P:C in highly reactive organic matter 
  yy            <- 0.27/106                         # molar ratio of P:C in refractory organic matter
  r0            <- 0.24                             # molar ratio of P:Fe for Fe-bound P 
#
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# organic C remineralization from (12-G) continuum theory
# ---------------------------------------------------------------------------------------------------------------- #
#
  a             <- 0.15                             # avg lifetime of reactive C_org [yr] 
  v             <- 0.12                             # shape parameter for C_org distribution
#
  integrand <- function(x) {(x^(v-1))*exp(-x)}
#
  Total = integrate(integrand, lower = 0, upper = Inf, rel.tol =1e-8)$value
#
  k=1:12
  f=1:12
  k[1]= 1
  k[12]=1e-10
#
  for (i in 2:11) 
      {k[i]=10^(-i+1+0.5)}
#
  for (i in 2:11) 
      {f[i]=(integrate(integrand, lower = 0, upper = (a*10^(-i+2)), rel.tol =1e-8)$value-integrate(integrand, 
      	    lower = 0, upper = (a*10^(-i+1)), rel.tol =1e-8)$value)/Total}
#
  f[1] = (Total-integrate(integrand, lower = 0, upper = (a*1), rel.tol =1e-8)$value)/Total
  f[12]= integrate(integrand, lower = 0, upper = (a*1e-10), rel.tol =1e-8)$value/Total
#
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# solid fluxes at the sediment-water interface
# ---------------------------------------------------------------------------------------------------------------- #
  J.orgC        <- 0.3                                    #  mmol/(cm*cm*yr)  total C_org
  J.orgC1       <- J.orgC*f[1]                            #  mmol/(cm*cm*yr)  labile C_org
  J.orgC2       <- J.orgC*f[2]                            #  mmol/(cm*cm*yr)  labile C_org
  J.orgC3       <- J.orgC*f[3]                            #  mmol/(cm*cm*yr)  less reactive C_org
  J.orgC4       <- J.orgC*f[4]                            #  mmol/(cm*cm*yr)  less reactive C_org
  J.orgC5       <- J.orgC*f[5]                            #  mmol/(cm*cm*yr)  less reactive C_org
  J.orgC6       <- J.orgC*f[6]                            #  mmol/(cm*cm*yr)  less reactive C_org
  J.orgC7       <- J.orgC*f[7]                            #  mmol/(cm*cm*yr)  less reactive C_org
  J.orgC8       <- J.orgC*f[8]                            #  mmol/(cm*cm*yr)  less reactive C_org
  J.orgC9       <- J.orgC*f[9]                            #  mmol/(cm*cm*yr)  less reactive C_org
  J.orgC10      <- J.orgC*f[10]                           #  mmol/(cm*cm*yr)  less reactive C_org
  J.orgC11      <- J.orgC*f[11]                           #  mmol/(cm*cm*yr)  less reactive C_org
  J.orgC12      <- J.orgC*f[12]                           #  mmol/(cm*cm*yr)  less reactive C_org
#  
  J.orgCc       <- J.orgC7+J.orgC8+J.orgC9+J.orgC10+J.orgC11+J.orgC12
#
  J.FeOH3a      <- 0.213/100*2.5/56*1000*0.2*0.36         #  mmol/(cm*cm*yr)  amorphous Fe
  J.FeOH3b      <- 0                                      #  mmol/(cm*cm*yr)  crystalline Fe
  J.MnO2a       <- 1000/1e6*2.5/55*1000*0.2*0.36          #  mmol/(cm*cm*yr)  amorphous Mn
  J.MnO2b       <- 0                                      #  mmol/(cm*cm*yr)  crystalline Mn
  J.irPa        <- J.FeOH3a*r0                            #  mmol/(cm*cm*yr)  P bound by amorphous Fe
  J.irPb        <- J.FeOH3b*r0                            #  mmol/(cm*cm*yr)  P bound by crystalline Fe
  J.S0          <- 0                                      #  mmol/(cm*cm*yr)  S0
  J.FeS         <- 0                                      #  mmol/(cm*cm*yr)  FeS  
  J.FeS2        <- 0                                      #  mmol/(cm*cm*yr)  FeS2  
  J.CFA         <- 0.1*436/1e6*2.5/31*1000/4.8*0.2*0.36   #  mmol/(cm*cm*yr)  CFA
  J.calc        <- 2.88/100*2.5/12*1000*0.2*0.36          #  mmol/(cm*cm*yr)  calcite
  J.arag        <- 0                                      #  mmol/(cm*cm*yr)  aragonite
  J.biot        <- 0.004                                  #  mmol/(cm*cm*yr)  biotite
  J.magn        <- 0.062/100*2.5/56*1000*0.2*0.36         #  mmol/(cm*cm*yr)  magnetite
  J.SurfFe      <- 0                                      #  mmol/(cm*cm*yr)  adsorbed iron
#
# ---------------------------------------------------------------------------------------------------------------- #
# 
# ---------------------------------------------------------------------------------------------------------------- #
# equilibrium constants
# ---------------------------------------------------------------------------------------------------------------- #
#
  TK            <- 273.15 + BW_Temp
  K1CO2         <- K1(BW_Sal, BW_Temp)*(D_sw/1000)                      # mmol/(cm*cm*cm)
  K2CO2         <- K2(BW_Sal, BW_Temp)*(D_sw/1000)                      # mmol/(cm*cm*cm)
  K1P           <- K1p(BW_Sal, BW_Temp)*(D_sw/1000)                     # mmol/(cm*cm*cm)
  K2P           <- K2p(BW_Sal, BW_Temp)*(D_sw/1000)                     # mmol/(cm*cm*cm)
  K3P           <- K3p(BW_Sal, BW_Temp)*(D_sw/1000)                     # mmol/(cm*cm*cm)
  KS            <- Ks(BW_Sal, BW_Temp)*(D_sw/1000)                      # mmol/(cm*cm*cm)
  KN            <- Kn(BW_Sal, BW_Temp)*(D_sw/1000)                      # mmol/(cm*cm*cm)
  KW            <- Kw(BW_Sal, BW_Temp)*(D_sw/1000)*(D_sw/1000)          # mmol*mmol/cm^6
  KHS           <- Khs(BW_Sal,BW_Temp)[[1]]*(D_sw/1000)                 # mmol/(cm*cm*cm)
  Kspa          <- Kspa(BW_Sal, BW_Temp)[[1]]*(D_sw/1000)*(D_sw/1000)   # mmol*mmol/cm^6
  Kspc          <- Kspc(BW_Sal, BW_Temp)[[1]]*(D_sw/1000)*(D_sw/1000)   # mmol*mmol/cm^6
  Kspr          <- 10^(-9)                                              # mmol*mmol/cm^6  
  KspCFA        <- 10^(-83.231)                                         # mmol*mmol/cm^6 
#
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# bottom water solute concentrations derived from equilibrium
# ---------------------------------------------------------------------------------------------------------------- #
#
  BW_CO2        <- BW_H*BW_H/(BW_H*K1CO2 + BW_H*BW_H + K1CO2*K2CO2)*BW_SumCO2
  BW_HCO3       <- BW_H*K1CO2/(BW_H*K1CO2 + BW_H*BW_H + K1CO2*K2CO2)*BW_SumCO2
  BW_CO3        <- K1CO2*K2CO2/(BW_H*K1CO2 + BW_H*BW_H + K1CO2*K2CO2)*BW_SumCO2
  BW_HS         <- BW_SumH2S/(1+BW_H/KHS)
  BW_H2S        <- BW_SumH2S - BW_HS
  BW_PO4        <- 1/(1+BW_H/K3P+BW_H*BW_H/K2P/K3P+BW_H*BW_H*BW_H/K1P/K2P/K3P)*BW_SumH2PO4 
  BW_HPO4       <- (BW_H/K3P)/(1+BW_H/K3P+BW_H*BW_H/K2P/K3P+BW_H*BW_H*BW_H/K1P/K2P/K3P)*BW_SumH2PO4 
  BW_H2PO4      <- (BW_H*BW_H/K2P/K3P)/(1+BW_H/K3P+BW_H*BW_H/K2P/K3P+BW_H*BW_H*BW_H/K1P/K2P/K3P)*BW_SumH2PO4
  BW_H3PO4      <- (BW_H*BW_H*BW_H/K1P/K2P/K3P)/(1+BW_H/K3P+BW_H*BW_H/K2P/K3P+BW_H*BW_H*BW_H/K1P/K2P/K3P)*BW_SumH2PO4
  BW_SO4        <- BW_SumSO4/(1+BW_H/KS)
  BW_HSO4       <- BW_SumSO4-BW_SO4
  BW_NH3        <- BW_SumNH4/(1+BW_H/KN)
  BW_NH4        <- BW_SumNH4-BW_NH3
#
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# activity coefficients for bottom water and pore fluid solutes
# ---------------------------------------------------------------------------------------------------------------- #
#
  rCa           <- 0.21866266
  rNa           <- 0.65687384
  rMg           <- 0.2819323
  rPO4          <- 0.000037
  rCO3          <- 0.028861
  rF            <- 0.31540788
#
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# evaluation of solute diffusion coefficients
# ---------------------------------------------------------------------------------------------------------------- # 
#
  select        <- c("SO4","Ca","O2","HCO3","CO2","CO3","HS","H2S","H",
  	                 "NO3","NH4","H3PO4","H2PO4","HPO4","PO4","Mn","Fe",
  	                 "CH4","Na","Mg","F","HSO4","NH3")
  diffsal       <- diffcoeff(S=BW_Sal,t=BW_Temp,P=BW_Pr,species=select)
#
  diffsal       <- diffsal *100*100*31556926    # convert diffusion coefficients to cm*cm/yr   
#
  DSO4       	<- diffsal[1,1]                  
  DCa			<- diffsal[1,2]
  DO2           <- diffsal[1,3] 	
  DHCO3      	<- diffsal[1,4]
  DCO2     		<- diffsal[1,5]
  DCO3     		<- diffsal[1,6]
  DHS       	<- diffsal[1,7]
  DH2S       	<- diffsal[1,8]
  DH            <- diffsal[1,9]
  DNO3          <- diffsal[1,10]
  DNH4          <- diffsal[1,11]
  DH3PO4        <- diffsal[1,12]
  DH2PO4        <- diffsal[1,13]
  DHPO4         <- diffsal[1,14]
  DPO4          <- diffsal[1,15]
  DMn           <- diffsal[1,16]
  DFe           <- diffsal[1,17]
  DCH4          <- diffsal[1,18]
  DNa           <- diffsal[1,19]
  DMg           <- diffsal[1,20]
  DF            <- diffsal[1,21]
  DHSO4         <- diffsal[1,22]
  DNH3          <- diffsal[1,23]
#
# ---------------------------------------------------------------------------------------------------------------- # 
#
# ---------------------------------------------------------------------------------------------------------------- # 
# bioturbation parameters
# ---------------------------------------------------------------------------------------------------------------- # 
#
  Db            <- 1e-7*31556926   # biodiffusion coefficient at surface [cm*cm/yr]
  xbt           <- 3               # biodiffusion attenuation coefficient [cm]
  a0            <- 100             # bioirrigation coefficient at surface [yr-1]
  xbi           <- 0.8             # bioirrigation attenuation coefficient [cm]
#
# ---------------------------------------------------------------------------------------------------------------- # 
#
# ---------------------------------------------------------------------------------------------------------------- # 
# parameters for reaction kinetics
# ---------------------------------------------------------------------------------------------------------------- # 
#
# === primary redox reactions ==================================================================================== # 
#
  Ks_O2         <- 0.02/1e3      # limiting [O2] for aerobic respiration [mmol/cm^3]
  Ks_NO3        <- 0.004/1e3     # limiting [NO3] for nitrate reduction [mmol/cm^3]
  Ks_MnO2       <- 0.032*2.5     # limiting [MnO2] for Mn-oxide reduction [mmol/cm^3]
  Ks_FeOH3      <- 0.065*2.5     # limiting [FeOH3] for Fe-oxide reduction [mmol/cm^3]
  Ks_SO4        <- 1.6/1e3       # limiting [SO4] for sulfate reduction [mmol/cm^3]
  a_SO4         <- 0.2           # attenuation factor for sulfate reduction [dimensionless]   
#
# === secondary reactions ======================================================================================== #
#
  k7            <- 1e7           # rate constant for ammonium oxidation [cm^3/mmol*yr]
  k8            <- 6.94*1e5      # rate constant for Mn2+ oxidation [cm^3/mmol*yr]
  k9            <- 1.4*1e8       # rate constant for Fe2+ oxidation [cm^3/mmol*yr]
  k10           <- 3*1e5         # rate constant for FeS oxidation [cm^3/mmol*yr]
  k11           <- 1.89*1e4      # rate constant for FeS2 oxidation [cm^3/mmol*yr]
  k12           <- 1.6*1e5       # rate constant for H2S oxidation [cm^3/mmol*yr]
  k13           <- 1e10          # rate constant for aerobic CH4 oxidation [cm^3/mmol*yr]
  k14           <- 3*1e6         # rate constant for Mn-oxide reduction w/Fe2+ [cm^3/mmol*yr]
  k15           <- 2*1e4         # rate constant for Mn-oxide reduction w/H2S [cm^3/mmol*yr]
  k16           <- 8*1e3         # rate constant for Fe-oxide reduction w/H2S [cm^3/mmol*yr]
  k17           <- 1.48e6        # rate constant for FeS formation [cm^3/mmol*yr]
  k18           <- 1e4           # rate constant for AOM [cm^3/mmol*yr]
  k19           <- 3.16          # rate constant for S0 disproportionation [yr-1]
  k20           <- 7.26*1e3      # rate constant for FeS2 formation from S0 [cm^3/mmol*yr]
  k21           <- 0.57          # rate constant for Fe-oxide aging [yr-1]
  k22           <- 1.7           # rate constant for Mn-oxide aging [yr-1]
  k23_precip    <- 3e-6          # rate constant for aragonite precipitation [mmol/cm^3*yr]
  k23_diss      <- 0.5           # rate constant for aragonite dissolution [yr-1]
  k24_precip    <- 3e-7          # rate constant for calcite precipitation [mmol/cm^3*yr]
  k24_diss      <- 0.5           # rate constant for calcite dissolution [yr-1]
  k25_precip    <- 3e-6          # rate constant for rhodochrosite precipitation [mmol/cm^3*yr]
  k25_diss      <- 0.25          # rate constant for rhodochrosite dissolution [yr-1]
  k26           <- 2.7e-8        # rate constant for CFA formation [mmol/cm^3/yr]
  k27           <- 3.25*1e3      # rate constant for FeS2 formation from FeS [cm^3/mmol*yr]
  k28           <- 3e-4          # rate constant for biotite dissolution [yr-1] 
  k29           <- 0.95e-5*570   # rate constant for Fe3O4 dissolution by H2S [mM-0.5 yr-1]
  k30           <- 5e6           # rate constant for Fe2+ sorption 
  kFe           <- 500           # adsorption coefficient of Fe 
  KPO4          <- 10e-6         # limiting concentration for PO4 sorption on Fe-oxides [M]
#
# ---------------------------------------------------------------------------------------------------------------- #
#
#
# **************************************************************************************************************** #
# SET UP MODEL DOMAIN AND PROPERTIES
# **************************************************************************************************************** #
#
# ---------------------------------------------------------------------------------------------------------------- #
# boundary conditions for grid
# ---------------------------------------------------------------------------------------------------------------- #                     
#
  L             <- 300                       # size of sediment column [cm]
  N             <- 20                        # number of grid layers 
  dbl           <- 0.1                       # thickness of diffusive boundary layer [cm]
  por.0         <- 0.64                      # porosity at sediment seawater interface 
  por.inf       <- 0.64                      # porosity at infinite depth
  por_coeff     <- 5                         # attenuation coefficient for porosity profile [cm]
  rho_sed       <- 2.5                       # sediment density [g/cm^3]
  MAR           <- rho_sed*0.36*0.2          # sediment mass accumulation rate [g/cm3/yr]
  v.0           <- MAR/(1-por.0)/rho_sed     # solid-phase advection velocity at SWI [cm/y]
  v.inf         <- NULL                      # solid-phase advection velocity at infinite depth [cm/y]
#
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# set up overall grid
# ---------------------------------------------------------------------------------------------------------------- #  
#
  grid          <- setup.grid.1D(x.up = 0, L = L, N = N, dx.1 = 0.3)               # create finite element grid 
  por.grid      <- setup.prop.1D(func = p.exp, grid = grid, y.0 = por.0, 
                                 y.inf = por.inf, x.att = por_coeff)               # setup porosity profile
  svf.grid      <- setup.prop.1D(func = p.exp, grid = grid, y.0 = 1-por.0, 
                                 y.inf = 1-por.inf, x.att = por_coeff)             # setup solid volume fraction
  tort          <- 1 - 2*log(por.grid$int)                                         # calculate tortuosity from porosity
  adv           <- setup.compaction.1D(v.0 = v.0, v.inf = v.inf, por.0 = por.0, 
                                       por.inf = por.inf, por.grid = por.grid)     # setup advection grid
  u.grid        <- adv$u                                                           # porewater advection velocities [cm/y]
  v.grid        <- adv$v                                                           # solid-phase advection velocities [cm/y]
  f.grid        <- (1-por.grid$mid)/por.grid$mid                                   # solid-phase fraction
#
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# set up bioturbation grid
# ---------------------------------------------------------------------------------------------------------------- # 
#
  Dbfunc        <- function (x) 
                   return(Db*exp(-(x/xbt)^2))                        # biodiffusion coefficients [cm2 s-1]
  Db.grid       <- setup.prop.1D(func = Dbfunc, grid = grid)$int     # biodiffusion coefficient grid
  afunc         <- function (x) 
                   return(a0*exp(-x/xbi))                            # bioirrigation coefficients [yr-1]
  a.grid        <- setup.prop.1D(func = afunc, grid = grid)$mid      # bioirrigation coefficient grid
#
# ---------------------------------------------------------------------------------------------------------------- #
#
#
# **************************************************************************************************************** #
# MODEL FORMULATION
# **************************************************************************************************************** #
#
  Pmodel        <- function (t = 0, Conc, pars = NULL) 

  {
#
# ---------------------------------------------------------------------------------------------------------------- #
# pass in initial concentrations
# ---------------------------------------------------------------------------------------------------------------- #
#
  orgC1         <- Conc[1:N]                 # labile C_org
  orgC2    		<- Conc[(N+1):(2*N)]         # labile C_org
  orgC3    		<- Conc[(2*N+1):(3*N)]       # less reactive C_org
  orgC4    		<- Conc[(3*N+1):(4*N)]       # less reactive C_org
  orgC5    		<- Conc[(4*N+1):(5*N)]       # less reactive C_org
  orgC6    		<- Conc[(5*N+1):(6*N)]       # less reactive C_org
  CFA           <- Conc[(6*N+1):(7*N)]       # carbonate fluorapatite
  O2            <- Conc[(7*N+1):(8*N)]       # oxygen
  NO3		    <- Conc[(8*N+1):(9*N)]       # nitrate
  MnO2a   		<- Conc[(9*N+1):(10*N)]      # amorphous manganese oxides
  FeOH3a        <- Conc[(10*N+1):(11*N)]     # amorphous iron oxides
  SO4           <- Conc[(11*N+1):(12*N)]     # sulfate
  CH4           <- Conc[(12*N+1):(13*N)]     # methane
  SumNH4        <- Conc[(13*N+1):(14*N)]     # ammonium
  MnII          <- Conc[(14*N+1):(15*N)]     # dissolved manganese
  FeII          <- Conc[(15*N+1):(16*N)]     # dissolved iron
  FeS           <- Conc[(16*N+1):(17*N)]     # iron monosulfide
  FeS2          <- Conc[(17*N+1):(18*N)]     # pyrite
  SumH2S        <- Conc[(18*N+1):(19*N)]     # hydrogen sulfide
  MnO2b         <- Conc[(19*N+1):(20*N)]     # crystalline manganese oxides
  SumH2PO4      <- Conc[(20*N+1):(21*N)]     # phosphate
  S0            <- Conc[(21*N+1):(22*N)]     # elemental sulfur
  FeOH3b        <- Conc[(22*N+1):(23*N)]     # crystalline iron oxides
  F             <- Conc[(23*N+1):(24*N)]     # fluorine
  SumCO2        <- Conc[(24*N+1):(25*N)]     # dissolved inorganic carbon
  H             <- Conc[(25*N+1):(26*N)]     # proton
  Ca            <- Conc[(26*N+1):(27*N)]     # calcium
  arag          <- Conc[(27*N+1):(28*N)]     # aragonite
  calc          <- Conc[(28*N+1):(29*N)]     # calcite
  rhod          <- Conc[(29*N+1):(30*N)]     # rhodochrosite
  Mg            <- Conc[(30*N+1):(31*N)]     # magnesium
  magn          <- Conc[(31*N+1):(32*N)]     # magnetite
  irPa          <- Conc[(32*N+1):(33*N)]     # phosphate bound by amorphous iron oxides     
  irPb          <- Conc[(33*N+1):(34*N)]     # phosphate bound by crystalline iron oxides
#
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# dissociation terms
# ---------------------------------------------------------------------------------------------------------------- # 
#
  HS            <- SumH2S/(1+H/KHS)
  H2S           <- SumH2S - HS
  PO4           <- 1/(1+H/K3P+H*H/K2P/K3P+H*H*H/K1P/K2P/K3P)*SumH2PO4 
  HPO4          <- (H/K3P)/(1+H/K3P+H*H/K2P/K3P+H*H*H/K1P/K2P/K3P)*SumH2PO4 
  H2PO4         <- (H*H/K2P/K3P)/(1+H/K3P+H*H/K2P/K3P+H*H*H/K1P/K2P/K3P)*SumH2PO4 
  H3PO4         <- (H*H*H/K1P/K2P/K3P)/(1+H/K3P+H*H/K2P/K3P+H*H*H/K1P/K2P/K3P)*SumH2PO4 
  NH3           <- SumNH4/(1+H/KN)
  NH4           <- SumNH4-NH3
  CO2           <- H*H/(H*K1CO2 + H*H + K1CO2*K2CO2)*SumCO2
  HCO3          <- H*K1CO2/(H*K1CO2 + H*H + K1CO2*K2CO2)*SumCO2
  CO3           <- K1CO2*K2CO2/(H*K1CO2 + H*H + K1CO2*K2CO2)*SumCO2
  SurfFe        <- kFe*FeII/f.grid 
#
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# differential equations for tracking pH
# ---------------------------------------------------------------------------------------------------------------- # 
#
  deltaTP.deltaSumCO2     <-  2*H*H/(H*K1CO2 + H*H + K1CO2*K2CO2) + H*K1CO2/(H*K1CO2 + H*H + K1CO2*K2CO2)
  deltaTP.deltaSumH2S     <-  H/(H+KHS)
  deltaTP.deltaSumH2PO4   <-  2*K1P*H*H/(K1P*K2P*K3P+K1P*K2P*H+K1P*H*H+H*H*H)+K1P*K2P*H/(K1P*K2P*K3P+K1P*K2P*H+K1P*H*H+H*H*H)
  deltaTP.deltaSumNH4     <-  H/(H+KN)
  deltaCO2.deltaH         <- (2*H*(H*K1CO2 + H*H + K1CO2*K2CO2)-(2*H + K1CO2)*H*H)/(H*K1CO2 + H*H + K1CO2*K2CO2)^2 * SumCO2
  deltaHCO3.deltaH        <- (K1CO2*(H*K1CO2 + H*H + K1CO2*K2CO2)-(2*H + K1CO2)*H*K1CO2)/(H*K1CO2 + H*H + K1CO2*K2CO2)^2 * SumCO2 
  deltaH2PO4.deltaH       <- (2*K1P*H*(K1P*K2P*K3P+K1P*K2P*H+K1P*H*H+H*H*H)-(3*H*H+2*H*K1P+K1P*K2P)*K1P*H*H)/(K1P*K2P*K3P+K1P*K2P*H+K1P*H*H+H*H*H)^2*SumH2PO4
  deltaHPO4.deltaH        <- (K1P*K2P*(K1P*K2P*K3P+K1P*K2P*H+K1P*H*H+H*H*H)-(3*H*H+2*H*K1P+K1P*K2P)*H*K1P*K2P)/(K1P*K2P*K3P+K1P*K2P*H+K1P*H*H+H*H*H)^2*SumH2PO4
  deltaNH4.deltaH         <- KN/(H+KN)^2*SumNH4
  deltaH2S.deltaH         <- (KHS /(H*H + 2*H*KHS + KHS*KHS)) * SumH2S
  deltaTP.deltaH          <- 2*deltaCO2.deltaH + deltaHCO3.deltaH + 2*deltaH2PO4.deltaH + deltaHPO4.deltaH + deltaH2S.deltaH + deltaNH4.deltaH + 1 
#
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# functions for solute/solid transport [molecular diffusion, biodiffusion, bioirrigation, advection]
# ---------------------------------------------------------------------------------------------------------------- # 
#
  tranorgC1     <- tran.1D(C = orgC1, flux.up = J.orgC1, D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranorgC2     <- tran.1D(C = orgC2, flux.up = J.orgC2, D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranorgC3     <- tran.1D(C = orgC3, flux.up = J.orgC3, D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranorgC4     <- tran.1D(C = orgC4, flux.up = J.orgC4, D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranorgC5     <- tran.1D(C = orgC5, flux.up = J.orgC5, D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranorgC6     <- tran.1D(C = orgC6, flux.up = J.orgC6, D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranCFA       <- tran.1D(C = CFA, flux.up = J.CFA, D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranMnO2a     <- tran.1D(C = MnO2a, flux.up = J.MnO2a, D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranFeOH3a    <- tran.1D(C = FeOH3a, flux.up = J.FeOH3a, D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranFeS       <- tran.1D(C = FeS, flux.up = J.FeS, D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranFeS2      <- tran.1D(C = FeS2, flux.up = J.FeS2, D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranMnO2b     <- tran.1D(C = MnO2b, flux.up = J.MnO2b, D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranS0        <- tran.1D(C = S0, flux.up = J.S0, D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranFeOH3b    <- tran.1D(C = FeOH3b, flux.up = J.FeOH3b, D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranarag      <- tran.1D(C = arag, flux.up = J.arag, D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC          
  trancalc      <- tran.1D(C = calc, flux.up = J.calc, D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranrhod      <- tran.1D(C = rhod, flux.up = 0, D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC 
  tranmagn      <- tran.1D(C = magn, flux.up = J.magn,  D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranSurfFe    <- tran.1D(C = SurfFe, flux.up = J.SurfFe,  D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranirPa      <- tran.1D(C = irPa, flux.up = J.irPa,  D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranirPb      <- tran.1D(C = irPb, flux.up = J.irPb,  D = Db.grid, v = v.grid, VF = svf.grid, dx = grid)$dC
  tranO2        <- tran.1D(C = O2, C.up = BW_O2, D = DO2/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(O2-BW_O2)
  tranNO3       <- tran.1D(C = NO3, C.up = BW_NO3, D = DNO3/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(NO3-BW_NO3)
  tranSO4       <- tran.1D(C = SO4, C.up = BW_SO4, D = DSO4/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(SO4-BW_SO4)
  tranCH4       <- tran.1D(C = CH4, C.up = BW_CH4, D = DCH4/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(CH4-BW_CH4)
  tranNH3       <- tran.1D(C = NH3, C.up = BW_NH3, D = DNH3/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(NH3-BW_NH3)
  tranNH4       <- tran.1D(C = NH4, C.up = BW_NH4, D = DNH4/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(NH4-BW_NH4)
  tranMnII      <- tran.1D(C = MnII, C.up = BW_MnII, D = DMn/tort, v = u.grid, VF = por.grid, dx = grid)$dC - 0.2*a.grid*(MnII-BW_MnII)
  tranFeII      <- tran.1D(C = FeII, C.up = BW_FeII, D = DFe/tort, v = u.grid, VF = por.grid, dx = grid)$dC 
  tranHS        <- tran.1D(C = HS, C.up = BW_HS, D = DHS/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(HS-BW_HS)
  tranH2S       <- tran.1D(C = H2S, C.up = BW_H2S, D = DH2S/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(H2S-BW_H2S)
  tranPO4       <- tran.1D(C = PO4, C.up = BW_PO4, D = DPO4/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(PO4-BW_PO4)
  tranHPO4      <- tran.1D(C = HPO4, C.up = BW_HPO4, D = DHPO4/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(HPO4-BW_HPO4)
  tranH2PO4     <- tran.1D(C = H2PO4, C.up = BW_H2PO4, D = DH2PO4/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(H2PO4-BW_H2PO4)
  tranH3PO4     <- tran.1D(C = H3PO4, C.up = BW_H3PO4, D = DH3PO4/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(H3PO4-BW_H3PO4)
  tranF         <- tran.1D(C = F, C.up = BW_F, D = DF/tort, v = u.grid, VF = por.grid, dx = grid)$dC - 2*a.grid*(F-BW_F)
  tranHCO3      <- tran.1D(C = HCO3, C.up = BW_HCO3, D = DHCO3/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(HCO3-BW_HCO3)
  tranCO2     	<- tran.1D(C = CO2, C.up = BW_CO2, D = DCO2/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(CO2-BW_CO2)
  tranCO3     	<- tran.1D(C = CO3, C.up = BW_CO3, D = DCO3/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(CO3-BW_CO3)
  tranH         <- tran.1D(C = H, C.up = BW_H, D = DH/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(H-BW_H)
  tranCa        <- tran.1D(C = Ca, C.up = BW_Ca, D = DCa/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(Ca-BW_Ca)
  tranMg        <- tran.1D(C = Mg, C.up = BW_Mg, D = DMg/tort, v = u.grid, VF = por.grid, dx = grid)$dC - a.grid*(Mg-BW_Mg)
#
  tranSumNH4    <- tranNH3 + tranNH4
  tranSumH2S    <- tranHS + tranH2S
  tranSumH2PO4  <- tranPO4 + tranHPO4 + tranH2PO4 + tranH3PO4
  tranSumCO2    <- tranCO3 + tranHCO3 + tranCO2
  tranTP        <- 2*tranCO2 + tranHCO3 + 2*tranH2PO4 + tranHPO4 + tranH2S + tranNH4 + tranH 
#
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# functions for reaction kinetics
# ---------------------------------------------------------------------------------------------------------------- # 
#
# === organic carbon remineralization ============================================================================ # 
#
  reacOC1i      <- k[1]*orgC1
  reacOC2i      <- k[2]*orgC2
  reacOC3i      <- k[3]*orgC3
  reacOC4i      <- k[4]*orgC4
  reacOC5i      <- k[5]*orgC5
  reacOC6i      <- k[6]*orgC6
#
# === primary redox reactions ==================================================================================== #
#
  reacorgC11    <- reacOC1i*O2/(Ks_O2+O2)                           # aerobic respiration
  reacorgC21    <- reacOC2i*O2/(Ks_O2+O2)       
  reacorgC31    <- reacOC3i*O2/(Ks_O2+O2)        
  reacorgC41    <- reacOC4i*O2/(Ks_O2+O2)         
  reacorgC51    <- reacOC5i*O2/(Ks_O2+O2)        
  reacorgC61    <- reacOC6i*O2/(Ks_O2+O2)            
#
  reacorgC12    <- reacOC1i*(Ks_O2/(Ks_O2+O2))*                     # nitrate reduction
                            (NO3/(Ks_NO3+NO3))     
  reacorgC22    <- reacOC2i*(Ks_O2/(Ks_O2+O2))*
                            (NO3/(Ks_NO3+NO3))     
  reacorgC32    <- reacOC3i*(Ks_O2/(Ks_O2+O2))*
                            (NO3/(Ks_NO3+NO3))     
  reacorgC42    <- reacOC4i*(Ks_O2/(Ks_O2+O2))*
                            (NO3/(Ks_NO3+NO3))       
  reacorgC52    <- reacOC5i*(Ks_O2/(Ks_O2+O2))*
                            (NO3/(Ks_NO3+NO3))     
  reacorgC62    <- reacOC6i*(Ks_O2/(Ks_O2+O2))*
                            (NO3/(Ks_NO3+NO3))       
#
  reacorgC13    <- reacOC1i*(Ks_O2/(Ks_O2+O2))*                     # Mn-oxide reduction
                            (Ks_NO3/(Ks_NO3+NO3))*
                            (MnO2a/(Ks_MnO2+MnO2a))    
  reacorgC23    <- reacOC2i*(Ks_O2/(Ks_O2+O2))*
                            (Ks_NO3/(Ks_NO3+NO3))*
                            (MnO2a/(Ks_MnO2+MnO2a))   
  reacorgC33    <- reacOC3i*(Ks_O2/(Ks_O2+O2))*
                            (Ks_NO3/(Ks_NO3+NO3))*
                            (MnO2a/(Ks_MnO2+MnO2a))    
  reacorgC43    <- reacOC4i*(Ks_O2/(Ks_O2+O2))*
                            (Ks_NO3/(Ks_NO3+NO3))*
                            (MnO2a/(Ks_MnO2+MnO2a))   
  reacorgC53    <- reacOC5i*(Ks_O2/(Ks_O2+O2))*
                            (Ks_NO3/(Ks_NO3+NO3))*
                            (MnO2a/(Ks_MnO2+MnO2a))   
  reacorgC63    <- reacOC6i*(Ks_O2/(Ks_O2+O2))*
                            (Ks_NO3/(Ks_NO3+NO3))*
                            (MnO2a/(Ks_MnO2+MnO2a))   
#
  reacorgC14    <- reacOC1i*(Ks_O2/(Ks_O2+O2))*                     # Fe-oxide reduction
                            (Ks_NO3/(Ks_NO3+NO3))*
                            (Ks_MnO2/(Ks_MnO2+MnO2a))*
                            (FeOH3a/(Ks_FeOH3+FeOH3a))   
  reacorgC24    <- reacOC2i*(Ks_O2/(Ks_O2+O2))*
                            (Ks_NO3/(Ks_NO3+NO3))*
                            (Ks_MnO2/(Ks_MnO2+MnO2a))*
                            (FeOH3a/(Ks_FeOH3+FeOH3a))    
  reacorgC34    <- reacOC3i*(Ks_O2/(Ks_O2+O2))*
                            (Ks_NO3/(Ks_NO3+NO3))*
                            (Ks_MnO2/(Ks_MnO2+MnO2a))*
                            (FeOH3a/(Ks_FeOH3+FeOH3a))   
  reacorgC44    <- reacOC4i*(Ks_O2/(Ks_O2+O2))*
                            (Ks_NO3/(Ks_NO3+NO3))*
                            (Ks_MnO2/(Ks_MnO2+MnO2a))*
                            (FeOH3a/(Ks_FeOH3+FeOH3a))   
  reacorgC54    <- reacOC5i*(Ks_O2/(Ks_O2+O2))*
                            (Ks_NO3/(Ks_NO3+NO3))*
                            (Ks_MnO2/(Ks_MnO2+MnO2a))*
                            (FeOH3a/(Ks_FeOH3+FeOH3a))   
  reacorgC64    <- reacOC6i*(Ks_O2/(Ks_O2+O2))*
                            (Ks_NO3/(Ks_NO3+NO3))*
                            (Ks_MnO2/(Ks_MnO2+MnO2a))*
                            (FeOH3a/(Ks_FeOH3+FeOH3a))   
#
  reacorgC15    <- a_SO4*reacOC1i*(Ks_O2/(Ks_O2+O2))*               # sulfate reduction
                                  (Ks_NO3/(Ks_NO3+NO3))*
                                  (Ks_MnO2/(Ks_MnO2+MnO2a))*
                                  (Ks_FeOH3/(Ks_FeOH3+FeOH3a))*
                                  (SO4/(Ks_SO4+SO4))           
  reacorgC25    <- a_SO4*reacOC2i*(Ks_O2/(Ks_O2+O2))*
                                  (Ks_NO3/(Ks_NO3+NO3))*
                                  (Ks_MnO2/(Ks_MnO2+MnO2a))*
                                  (Ks_FeOH3/(Ks_FeOH3+FeOH3a))*
                                  (SO4/(Ks_SO4+SO4))          
  reacorgC35    <- a_SO4*reacOC3i*(Ks_O2/(Ks_O2+O2))*
                                  (Ks_NO3/(Ks_NO3+NO3))*
                                  (Ks_MnO2/(Ks_MnO2+MnO2a))*
                                  (Ks_FeOH3/(Ks_FeOH3+FeOH3a))*
                                  (SO4/(Ks_SO4+SO4))           
  reacorgC45    <- a_SO4*reacOC4i*(Ks_O2/(Ks_O2+O2))*
                                  (Ks_NO3/(Ks_NO3+NO3))*
                                  (Ks_MnO2/(Ks_MnO2+MnO2a))*
                                  (Ks_FeOH3/(Ks_FeOH3+FeOH3a))*
                                  (SO4/(Ks_SO4+SO4))           
  reacorgC55    <- a_SO4*reacOC5i*(Ks_O2/(Ks_O2+O2))*
                                  (Ks_NO3/(Ks_NO3+NO3))*
                                  (Ks_MnO2/(Ks_MnO2+MnO2a))*
                                  (Ks_FeOH3/(Ks_FeOH3+FeOH3a))*
                                  (SO4/(Ks_SO4+SO4))           
  reacorgC65    <- a_SO4*reacOC6i*(Ks_O2/(Ks_O2+O2))*
                                  (Ks_NO3/(Ks_NO3+NO3))*
                                  (Ks_MnO2/(Ks_MnO2+MnO2a))*
                                  (Ks_FeOH3/(Ks_FeOH3+FeOH3a))*
                                  (SO4/(Ks_SO4+SO4))           
#
  reacorgC16    <- a_SO4*reacOC1i*(Ks_O2/(Ks_O2+O2))*               # methanogenesis
                                  (Ks_NO3/(Ks_NO3+NO3))*
                                  (Ks_MnO2/(Ks_MnO2+MnO2a))*
                                  (Ks_FeOH3/(Ks_FeOH3+FeOH3a))*
                                  (Ks_SO4/(Ks_SO4+SO4))           
  reacorgC26    <- a_SO4*reacOC2i*(Ks_O2/(Ks_O2+O2))*
                                  (Ks_NO3/(Ks_NO3+NO3))*
                                  (Ks_MnO2/(Ks_MnO2+MnO2a))*
                                  (Ks_FeOH3/(Ks_FeOH3+FeOH3a))*
                                  (Ks_SO4/(Ks_SO4+SO4))          
  reacorgC36    <- a_SO4*reacOC3i*(Ks_O2/(Ks_O2+O2))*
                                  (Ks_NO3/(Ks_NO3+NO3))*
                                  (Ks_MnO2/(Ks_MnO2+MnO2a))*
                                  (Ks_FeOH3/(Ks_FeOH3+FeOH3a))*
                                  (Ks_SO4/(Ks_SO4+SO4))           
  reacorgC46    <- a_SO4*reacOC4i*(Ks_O2/(Ks_O2+O2))*
                                  (Ks_NO3/(Ks_NO3+NO3))*
                                  (Ks_MnO2/(Ks_MnO2+MnO2a))*
                                  (Ks_FeOH3/(Ks_FeOH3+FeOH3a))*
                                  (Ks_SO4/(Ks_SO4+SO4))           
  reacorgC56    <- a_SO4*reacOC5i*(Ks_O2/(Ks_O2+O2))*
                                  (Ks_NO3/(Ks_NO3+NO3))*
                                  (Ks_MnO2/(Ks_MnO2+MnO2a))*
                                  (Ks_FeOH3/(Ks_FeOH3+FeOH3a))*
                                  (Ks_SO4/(Ks_SO4+SO4))          
  reacorgC66    <- a_SO4*reacOC6i*(Ks_O2/(Ks_O2+O2))*
                                  (Ks_NO3/(Ks_NO3+NO3))*
                                  (Ks_MnO2/(Ks_MnO2+MnO2a))*
                                  (Ks_FeOH3/(Ks_FeOH3+FeOH3a))*
                                  (Ks_SO4/(Ks_SO4+SO4))           
#
  reacOC1       <- reacorgC11+reacorgC12+reacorgC13+reacorgC14+reacorgC15+reacorgC16
  reacOC2       <- reacorgC21+reacorgC22+reacorgC23+reacorgC24+reacorgC25+reacorgC26
  reacOC3       <- reacorgC31+reacorgC32+reacorgC33+reacorgC34+reacorgC35+reacorgC36
  reacOC4       <- reacorgC41+reacorgC42+reacorgC43+reacorgC44+reacorgC45+reacorgC46
  reacOC5       <- reacorgC51+reacorgC52+reacorgC53+reacorgC54+reacorgC55+reacorgC56
  reacOC6       <- reacorgC61+reacorgC62+reacorgC63+reacorgC64+reacorgC65+reacorgC66
#
  reacorgCa1    <- reacorgC11+reacorgC21
  reacorgCa2    <- reacorgC12+reacorgC22
  reacorgCa3    <- reacorgC13+reacorgC23
  reacorgCa4    <- reacorgC14+reacorgC24
  reacorgCa5    <- reacorgC15+reacorgC25
  reacorgCa6    <- reacorgC16+reacorgC26
  reacorgCb1    <- reacorgC31+reacorgC41+reacorgC51+reacorgC61
  reacorgCb2    <- reacorgC32+reacorgC42+reacorgC52+reacorgC62
  reacorgCb3    <- reacorgC33+reacorgC43+reacorgC53+reacorgC63
  reacorgCb4    <- reacorgC34+reacorgC44+reacorgC54+reacorgC64
  reacorgCb5    <- reacorgC35+reacorgC45+reacorgC55+reacorgC65
  reacorgCb6    <- reacorgC36+reacorgC46+reacorgC56+reacorgC66
#
  reacorgCa     <- reacOC1+reacOC2
  reacorgCb     <- reacOC3+reacOC4+reacOC5+reacOC6
#
  reacorgC4P    <- (reacOC1i+reacOC2i+reacOC3i+reacOC4i+reacOC5i+reacOC6i)*
                   (Ks_O2/(Ks_O2+O2))*
                   (Ks_NO3/(Ks_NO3+NO3))*
                   (Ks_MnO2/(Ks_MnO2+MnO2a))*
                   (irPa/(Ks_FeOH3+FeOH3a))   

# === secondary reactions ======================================================================================== #  
#
# --- saturation terms for mineral precipitation/dissolution kinetics -------------------------------------------- #
  omegaa        <- (CO3 * Ca)/Kspa
  omegac        <- (CO3 * Ca)/Kspc
  omegar        <- (CO3 * MnII)/Kspr
  omegaCFA      <- (Ca*rCa)^9.54*
                   (BW_Na*rNa)^0.33*
                   (Mg*rMg)^0.13*
                   (PO4*rPO4)^4.8*
                   (CO3*rCO3)^(1.2-2.3307)*
                   (F*rF)^2.48/KspCFA
# ---------------------------------------------------------------------------------------------------------------- #
#  reac7         <- k7*O2*SumNH4                            # R7
#  reac8         <- k8*O2*MnII                              # R8
#  reac9         <- k9*O2*FeII                              # R9
#  reac10        <- k10*O2*FeS                              # R10
#  reac11        <- k11*O2*FeS2                             # R11
#  reac12        <- k12*O2*SumH2S                           # R12
#  reac13        <- k13*O2*CH4                              # R13
#  reac14a       <- k14*MnO2a*FeII                          # R14a
#  reac14b       <- k14*MnO2b*FeII                          # R14b
#  reac14        <- reac14a+reac14b                         # R14
#  reac15a       <- k15*MnO2a*SumH2S                        # R15a
#  reac15b       <- k15*MnO2b*SumH2S                        # R15a
#  reac15        <- reac15a+reac15b                         # R15
#  reac16a       <- k16*FeOH3a*SumH2S                       # R16a
#  reac16b       <- k16*FeOH3b*SumH2S                       # R16b
#  reac16aP      <- k16*irPa*SumH2S                        # R16a
#  reac16bP      <- k16*irPb*SumH2S                        # R16b
#  reac16        <- reac16a+reac16b                         # R16
#  reac17        <- k17*FeII*SumH2S                         # R17
#  reac18        <- k18*SO4*CH4                             # R18
#  reac19        <- k19*S0                                  # R19
#  reac20        <- k20*FeS*S0                              # R20
#  reac21        <- k21*FeOH3a                              # R21
#  reac21P       <- k21*irPa                               # R21
#  reac22        <- k22*MnO2a                               # R22
#  reac23        <- ifelse (omegaa>1,k23_precip*(omegaa-1),
#  	                               -k23_diss*(1-omegaa)*arag) 
#  reac24        <- ifelse (omegac>1,k24_precip*(omegac-1),
#  	                               -k24_diss*(1-omegac)*calc)  
#  reac25        <- ifelse (omegar>1,k25_precip* (omegar-1),
#  	                               -k25_diss*(1-omegar)*rhod)
#  reac26        <- ifelse(omegaCFA>1,k26*(omegaCFA-1),0)       
#  reac27        <- k27*FeS*SumH2S                          # R28
#  reac28        <- k28*J.biot/(1-por.0)/v.0          # R29
#  reac29        <- magn*k29*((SumH2S*1000)^0.5)         # R30
#  reac30        <- k30*SurfFe*O2                           # R32
#
  reac7         <- k7*O2*SumNH4                                  # ammonium oxidation
  reac8         <- k8*O2*MnII                                    # Mn2+ oxidation
  reac9         <- k9*O2*FeII                                    # Fe2+ oxidation
  reac9P        <- k9*O2*FeII*r0*(SumH2PO4/(SumH2PO4+KPO4))      # P sorption during Fe2+ oxidation
  reac10        <- k10*O2*FeS                                    # FeS oxidation
  reac11        <- k11*O2*FeS2                                   # FeS2 oxidation
  reac12        <- k12*O2*SumH2S                                 # H2S oxidation
  reac13        <- k13*O2*CH4                                    # aerobic CH4 oxidtion
  reac14a       <- k14*MnO2a*FeII                                # amorphous Mn-oxide reduction w/Fe2+
  reac14aP      <- k14*MnO2a*FeII*r0*(SumH2PO4/(SumH2PO4+KPO4))  # P sorption during Fe2+ oxidation via amorphous Mn-oxides
  reac14b       <- k14*MnO2b*FeII                                # crystalline Mn-oxide reduction w/Fe2+
  reac14bP      <- k14*MnO2b*FeII*r0*(SumH2PO4/(SumH2PO4+KPO4))  # P sorption during Fe2+ oxidation via crystalline Mn-oxides
  reac14        <- reac14a+reac14b                               # total Mn-oxide reduction w/Fe2+
  reac14P       <- reac14aP+reac14bP                             # total P sorption during Fe2+ oxidation via Mn-oxides
  reac15a       <- k15*MnO2a*SumH2S                              # amorphous Mn-oxide reduction w/H2S
  reac15b       <- k15*MnO2b*SumH2S                              # crystalline Mn-oxide reduction w/H2S
  reac15        <- reac15a+reac15b                               # total Mn-oxide reduction w/H2S
  reac16a       <- k16*FeOH3a*SumH2S                             # amorphous Fe-oxide reduction w/H2S
  reac16b       <- k16*FeOH3b*SumH2S                             # crystalline Fe-oxide reduction w/H2S
  reac16aP      <- k16*irPa*SumH2S                               # P release from amorphous Fe-oxide reduction w/H2S
  reac16bP      <- k16*irPb*SumH2S                               # P release from crystalline Fe-oxide reduction w/H2S
  reac16        <- reac16a+reac16b                               # total Fe-oxide reduction w/H2S
  reac16P       <- reac16aP+reac16bP                             # total P release from Fe-oxide reduction
  reac17        <- k17*FeII*SumH2S                               # FeS formation
  reac18        <- k18*SO4*CH4                                   # AOM
  reac19        <- k19*S0                                        # S0 disproportionation
  reac20        <- k20*FeS*S0                                    # FeS2 formation from S0
  reac21        <- k21*FeOH3a                                    # Fe-oxide aging
  reac21P       <- k21*irPa                                      # P associated with aged Fe-oxides
  reac22        <- k22*MnO2a                                     # Mn-oxide aging
  reac23        <- ifelse(omegaa>1,k23_precip*(omegaa-1),        # aragonite precipitation/dissolution
  	                              -k23_diss*(1-omegaa)*arag)
  reac24        <- ifelse(omegac>1,k24_precip*(omegac-1),        # calcite precipitation/dissolution
  	                              -k24_diss*(1-omegac)*calc)
  reac25        <- ifelse(omegar>1,k25_precip*(omegar-1),        # rhodochrosite precipitation/dissolution
  	                              -k25_diss*(1-omegar)*rhod)
  reac26        <- ifelse(omegaCFA>1,k26*(omegaCFA-1),0)         # CFA precipitation/dissolution
  reac27        <- k27*FeS*SumH2S                                # FeS2 formation from FeS
  reac28        <- k28*J.biot/(1-por.0)/v.0                      # biotite dissolution
  reac29        <- magn*k29*((SumH2S*1000)^0.5)                  # Fe3O4 dissolution by H2S
  reac30        <- k30*SurfFe*O2                                 # Fe2+ sorption
  reac30P       <- k30*SurfFe*O2*r0*(SumH2PO4/(SumH2PO4+KPO4))   # P sorption
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# differential reaction-transport equations for solutes and solids
# ---------------------------------------------------------------------------------------------------------------- # 
#
  dorgC1        <- tranorgC1-reacOC1
  dorgC2        <- tranorgC2-reacOC2
  dorgC3        <- tranorgC3-reacOC3
  dorgC4        <- tranorgC4-reacOC4
  dorgC5        <- tranorgC5-reacOC5
  dorgC6        <- tranorgC6-reacOC6
#  
  dCFA          <- tranCFA+reac26
  dO2           <- tranO2-
                   (1+2*x)*f.grid*reacorgCa1-
                   (1+2*xx)*f.grid*reacorgCb1-
                   2*reac7-
                   0.5*reac8-
                   0.25*reac9-
                   2*f.grid*reac10-
                   3.5*f.grid*reac11-
                   2*reac12-
                   2*reac13-
                   0.25*f.grid*reac30
  dNO3          <- tranNO3+
                   x*f.grid*reacorgCa1+
                   xx*f.grid*reacorgCb1-
                   (0.8+0.6*x)*f.grid*reacorgCa2-
                   (0.8+0.6*xx)*f.grid*reacorgCb2+
                   reac7 
  dMnO2a        <- tranMnO2a-
                   2*reacorgCa3-
                   2*reacorgCb3+
                   1/f.grid*reac8-
                   reac14a-
                   reac15a-
                   reac22
  dFeOH3a       <- tranFeOH3a-
                   4*reacorgCa4-
                   4*reacorgCb4+
                   1/f.grid*reac9+
                   2*reac14-
                   2*reac16a-
                   reac21+
                   reac30
  dSO4          <- tranSO4-
                   0.5*f.grid*(reacorgCa5+reacorgCb5)+
                   f.grid*reac10+
                   2*f.grid*reac11+
                   reac12-
                   reac18+
                   f.grid*reac19
  dCH4          <- tranCH4+
                   0.5*f.grid*(reacorgCa6+reacorgCb6)-
                   reac13-
                   reac18
  dSumNH4       <- tranSumNH4+
                   f.grid*(reacorgCa3*x+reacorgCb3*xx+reacorgCa4*x+reacorgCb4*xx+reacorgCa5*x+reacorgCb5*xx+reacorgCa6*x+reacorgCb6*xx)-
                   reac7
  dMnII         <- tranMnII+
                   2*f.grid*(reacorgCa3+reacorgCb3)-
                   reac8+
                   f.grid*reac14+
                   f.grid*reac15-
                   f.grid*reac25
  RTFeII        <- tranFeII+
                   4*f.grid*(reacorgCa4+reacorgCb4)-
                   reac9+
                   f.grid*reac10+
                   f.grid*reac11-
                   2*f.grid*reac14+
                   2*f.grid*reac16-
                   reac17+
                   2*f.grid*reac28+
                   3*f.grid*reac29 
  RTSurfFe      <- tranSurfFe-reac30
  dFeII         <- 1/(1+kFe)*RTFeII + 1/(1+kFe)*f.grid*RTSurfFe 
  dFeS          <- tranFeS-
                   reac10+
                   1/f.grid*reac17-
                   reac20-
                   reac27
  dFeS2         <- tranFeS2-
                   reac11+
                   reac20+
                   reac27
  dSumH2S       <- tranSumH2S+
                   0.5*f.grid*(reacorgCa5+reacorgCb5)-
                   reac12-
                   f.grid*reac15-
                   f.grid*reac16-
                   reac17+
                   reac18+
                   3*f.grid*reac19-
                   f.grid*reac27-
                   f.grid*reac29 
  dMnO2b        <- tranMnO2b-
                   reac14b-
                   reac15b+
                   reac22 
  dSumH2PO4     <- tranSumH2PO4+
                   f.grid*(reacorgCa*y+reacorgCb*yy)-
                   4.8*f.grid*reac26+
                   4*f.grid*reacorgC4P-
                   reac9P-
                   2*f.grid*reac14P+
                   2*f.grid*reac16P-
                   f.grid*reac30P  
  dS0           <- tranS0+
                   reac15+
                   reac16-
                   4*reac19-
                   reac20+
                   reac29
  dFeOH3b       <- tranFeOH3b-
                   2*reac16b+
                   reac21
  dF            <- tranF-
                   2.48*f.grid*reac26
  dSumCO2       <- tranSumCO2+
                   f.grid*(reacorgCa + reacorgCb - 0.5*reacorgCa6 - 0.5*reacorgCb6)+
                   reac13+
                   reac18-
                   f.grid*reac23-
                   f.grid*reac24-
                   f.grid*reac25-
                   1.2*f.grid*reac26
  dTP           <- tranTP+
                   (2+x)*f.grid*reacorgCa1+
                   (2+xx)*f.grid*reacorgCb1+
                   1.2*f.grid*(reacorgCa2+reacorgCb2)-
                   0.6*f.grid*(reacorgCa2*x+reacorgCb2*xx)-
                   2*f.grid*(reacorgCa3+reacorgCb3)-
                   6*f.grid*(reacorgCa4+reacorgCb4)+
                   1.5*f.grid*(reacorgCa5+reacorgCb5)+
                   f.grid*(reacorgCa6+reacorgCb6)+
                   3*f.grid*(reacorgCa*y+reacorgCb*yy)+
                   12*f.grid*reacorgC4P+
                   reac7+
                   2*reac8+
                   (2-3*r0*SumH2PO4/(SumH2PO4+KPO4))*reac9+
                   2*f.grid*reac11+
                   reac12+
                   2*reac13+
                   (2-6*r0*SumH2PO4/(SumH2PO4+KPO4))*f.grid*reac14-
                   3*f.grid*reac15-
                   5*f.grid*(reac16a+reac16b)+
                   6*f.grid*(reac16aP+reac16bP)+
                   reac17+
                   reac18+
                   5*f.grid*reac19-
                   f.grid*reac26-
                   7*f.grid*reac28-
                   7*f.grid*reac29+
                   (1-3*r0*SumH2PO4/(SumH2PO4+KPO4))*f.grid*reac30+
                   kFe/(1+kFe)*RTFeII-1/(1+kFe)*f.grid*RTSurfFe
  dH            <- dTP/deltaTP.deltaH-
                   (dSumCO2*deltaTP.deltaSumCO2+
                   	dSumH2S*deltaTP.deltaSumH2S+
                   	dSumH2PO4*deltaTP.deltaSumH2PO4+
                   	dSumNH4*deltaTP.deltaSumNH4)/deltaTP.deltaH
  dCa           <- tranCa-
                   f.grid*reac23-
                   f.grid*reac24-
                   9.54*f.grid*reac26
  darag         <- tranarag + reac23
  dcalc         <- trancalc + reac24
  drhod         <- tranrhod + reac25
  dMg           <- tranMg - 0.13*f.grid*reac26  
  dmagn         <- tranmagn - reac29 
  dirPa         <- tranirPa-
                   4*reacorgC4P+
                   1/f.grid*reac9P+
                   2*reac14P-
                   2*reac16aP-
                   reac21P+
                   reac30P
  dirPb         <- tranirPb-
                   2*reac16bP+
                   reac21P
#
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# collect terms
# ---------------------------------------------------------------------------------------------------------------- # 
#
  return(list(c(
  	            dorgC1    = dorgC1,
  	            dorgC2    = dorgC2,
  	            dorgC3    = dorgC3,
  	            dorgC4    = dorgC4,
  	            dorgC5    = dorgC5,
  	            dorgC6    = dorgC6,
  	            dCFA      = dCFA,
  	            dO2       = dO2,
  	            dNO3      = dNO3,
  	            dMnO2a    = dMnO2a,
  	            dFeOH3a   = dFeOH3a,
  	            dSO4      = dSO4,
  	            dCH4      = dCH4,
  	            dSumNH4   = dSumNH4,
  	            dMnII     = dMnII,
  	            dFeII     = dFeII,
  	            dFeS      = dFeS,
  	            dFeS2     = dFeS2,
  	            dSumH2S   = dSumH2S,
  	            dMnO2b    = dMnO2b,
  	            dSumH2PO4 = dSumH2PO4,
  	            dS0       = dS0,
  	            dFeOH3b   = dFeOH3b,
  	            dF        = dF,
  	            dSumCO2   = dSumCO2,
  	            dH        = dH,
  	            dCa       = dCa,
  	            darag     = darag,
  	            dcalc     = dcalc,
  	            drhod     = drhod,
  	            dMg       = dMg,
  	            dmagn     = dmagn,
  	            dirPa     = dirPa,
  	            dirPb     = dirPb
  	           ),
              omegaa      = omegaa,
              omegac      = omegac,
              omegaCFA    = omegaCFA,
              PO4         = PO4,
              HPO4        = HPO4,
              H2PO4       = H2PO4,
              H3PO4       = H3PO4))

  }
#
# ---------------------------------------------------------------------------------------------------------------- #
#
#
# **************************************************************************************************************** #
# MODEL SOLUTION
# **************************************************************************************************************** #
#
# ---------------------------------------------------------------------------------------------------------------- #
# initialize state variables
# ---------------------------------------------------------------------------------------------------------------- #
#
  orgC1         <- rep(0,N)
  orgC2         <- rep(0,N)
  orgC3         <- rep(0,N)
  orgC4         <- rep(0,N)
  orgC5         <- rep(0,N)
  orgC6         <- rep(0,N)
  CFA           <- rep(0,N)
  O2            <- rep(0,N)
  NO3           <- rep(0,N)
  MnO2a         <- rep(0,N)
  FeOH3a        <- rep(0,N)
  SO4           <- rep(0,N)
  CH4           <- rep(0,N)
  SumNH4        <- rep(0,N)
  MnII          <- rep(0,N)
  FeII          <- rep(0,N)
  FeS           <- rep(0,N)
  FeS2          <- rep(0,N)
  SumH2S        <- rep(0,N)
  MnO2b         <- rep(0,N)
  SumH2PO4      <- rep(0,N) 
  S0            <- rep(0,N)
  FeOH3b        <- rep(0,N)
  F             <- rep(0,N)
  SumCO2        <- rep(BW_SumCO2,N)
  H             <- rep(BW_H,N)
  Ca            <- rep(0,N)
  arag          <- rep(0,N)
  calc          <- rep(0,N)
  rhod          <- rep(0,N)
  Mg            <- rep(0,N)
  magn          <- rep(0,N) 
  irPa          <- rep(0,N)
  irPb          <- rep(0,N)
#
  Conc          <- c(orgC1,
  	                 orgC2,
  	                 orgC3,
  	                 orgC4,
  	                 orgC5,
  	                 orgC6,
  	                 CFA,
  	                 O2,
  	                 NO3,
  	                 MnO2a,
  	                 FeOH3a,
  	                 SO4,
  	                 CH4,
  	                 SumNH4,
  	                 MnII,
  	                 FeII,
  	                 FeS,
  	                 FeS2,
  	                 SumH2S,
  	                 MnO2b,
  	                 SumH2PO4,
  	                 S0,
  	                 FeOH3b,
  	                 F,
  	                 SumCO2,
  	                 H,
  	                 Ca,
  	                 arag,
  	                 calc,
  	                 rhod,
  	                 Mg,
  	                 magn,
  	                 irPa,
  	                 irPb)
#
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# run model
# ---------------------------------------------------------------------------------------------------------------- #
#
  options(max.print=20000)
#
  times             <- seq(0, 20000, by = 1000) 
#
  system.time (sol <- ode.1D(y=Conc,
  	                         time=times,
  	                         func=Pmodel,
  	                         maxstep=10000,
  	                         parms=NULL,
  	                         nspec=34,
  	                         method="vode",
  	                         rtol=1e-14,
  	                         atol=1e-14,
  	                         verbose=TRUE,
  	                         mf=-22))
#
# ---------------------------------------------------------------------------------------------------------------- #
#
#
# **************************************************************************************************************** #
# PROCESS AND DISPLAY OUTPUT
# **************************************************************************************************************** #
#
# ---------------------------------------------------------------------------------------------------------------- #
# collect output
# ---------------------------------------------------------------------------------------------------------------- #
#
  orgC1         <- sol[21, (2):(N+1)]
  orgC2         <- sol[21, (N + 2):(2 * N + 1)]
  orgC3         <- sol[21, (2*N + 2):(3 * N + 1)]
  orgC4         <- sol[21, (3*N + 2):(4 * N + 1)]
  orgC5         <- sol[21, (4*N + 2):(5 * N + 1)]
  orgC6         <- sol[21, (5*N + 2):(6 * N + 1)]
  CFA           <- sol[21, (6*N + 2):(7 * N + 1)]
  O2            <- sol[21, (7*N + 2):(8 * N + 1)]
  NO3           <- sol[21, (8*N + 2):(9 * N + 1)]
  MnO2a         <- sol[21, (9*N + 2):(10 * N + 1)]
  FeOH3a        <- sol[21, (10*N + 2):(11 * N + 1)]
  SO4           <- sol[21, (11*N + 2):(12 * N + 1)]
  CH4           <- sol[21, (12*N + 2):(13 * N + 1)]
  SumNH4        <- sol[21, (13*N + 2):(14 * N + 1)]
  MnII          <- sol[21, (14*N + 2):(15 * N + 1)]
  FeII          <- sol[21, (15*N + 2):(16 * N + 1)]
  FeS           <- sol[21, (16*N + 2):(17 * N + 1)]
  FeS2          <- sol[21, (17*N + 2):(18 * N + 1)]
  SumH2S        <- sol[21, (18*N + 2):(19 * N + 1)]
  MnO2b         <- sol[21, (19*N + 2):(20 * N + 1)]
  SumH2PO4      <- sol[21, (20*N + 2):(21 * N + 1)]
  S0            <- sol[21, (21*N + 2):(22 * N + 1)]
  FeOH3b        <- sol[21, (22*N + 2):(23 * N + 1)]
  F             <- sol[21, (23*N + 2):(24 * N + 1)]
  SumCO2        <- sol[21, (24*N + 2):(25 * N + 1)]
  H             <- sol[21, (25*N + 2):(26 * N + 1)]
  Ca            <- sol[21, (26*N + 2):(27 * N + 1)]
  arag          <- sol[21, (27*N + 2):(28 * N + 1)]
  calc          <- sol[21, (28*N + 2):(29 * N + 1)]
  rhod          <- sol[21, (29*N + 2):(30 * N + 1)]
  Mg            <- sol[21, (30*N + 2):(31 * N + 1)]
  magn          <- sol[21, (31*N + 2):(32 * N + 1)]
  irPa          <- sol[21, (32*N + 2):(33 * N + 1)]
  irPb          <- sol[21, (33*N + 2):(34 * N + 1)]
  omegaa        <- sol[21, (34*N + 2):(35 * N + 1)]
  omegac        <- sol[21, (35*N + 2):(36 * N + 1)]
  omegaCFA      <- sol[21, (36*N + 2):(37 * N + 1)]
  PO4           <- sol[21, (37*N + 2):(38 * N + 1)]
  HPO4          <- sol[21, (38*N + 2):(39 * N + 1)]
  H2PO4         <- sol[21, (39*N + 2):(40 * N + 1)]
  H3PO4         <- sol[21, (40*N + 2):(41 * N + 1)]
#
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# convert output
# ---------------------------------------------------------------------------------------------------------------- #
#
  orgCat        <- (orgC1+orgC2)*12/1000/rho_sed*100                                  # labile C_org [wt%]
  orgCbt        <- (orgC3+orgC4+orgC5+orgC6)*12/1000/rho_sed*100                      # less reactive C_org [wt%]
  orgCct        <- (J.orgCc/0.2/0.36)*12/1000/rho_sed*100                             # refractory C_org [wt%]
  orgCt         <- orgCat + orgCbt + orgCct                                           # total C_org [wt%]
  orgNt         <- orgCat*x + orgCbt*xx + orgCct*xx                                   # total N_org [wt%]
  CFA_Pt        <- CFA*4.8*1000/rho_sed                                               # CFA phosphorus [umol/g]
  orgPt         <- ((orgC1+orgC2)*y+                                                  # total P_org [umol/g]
  	                (orgC3+orgC4+orgC5+orgC6+J.orgCc/0.2/0.36)*yy)*
                   1000/rho_sed
  O2t           <- O2*1e6                                                             # [O2] [umol/L]
  NO3t          <- NO3*1e6                                                            # [NO3] [umol/L]
  MnO2t         <- (MnO2a + MnO2b)*1000/rho_sed + 88/55                               # [MnO2] [umol/g]
  FeOH3t        <- (FeOH3a + FeOH3b)*1000/rho_sed                                     # [FeOH3] [umol/g]
  SO4t          <- SO4*1e3                                                            # [SO4] [mmol/L]
  CH4t          <- CH4*1e6                                                            # [CH4] [umol/L]
  SumNH4t       <- SumNH4*1e3                                                         # [NH4] [mmol/L]
  MnIIt         <- MnII*1e6                                                           # [Mn2+] [umol/L]
  FeIIt         <- FeII*1e6                                                           # [Fe2+] [umol/L]
  SumH2St       <- SumH2S*1e6                                                         # [H2S] [umol/L]
  SumH2PO4t     <- SumH2PO4*1e6                                                       # [PO4] [umol/L]
  S0t           <- S0*1000/rho_sed                                                    # [S0] [umol/g]
  SumCO2t       <- SumCO2*1e3                                                         # [DIC] [mmol/L]
  pH            <- -log10(H)                                                          # pH
  Cat           <- Ca*1e3                                                             # [Ca] [mmol/L]
  Mgt           <- Mg*1e3                                                             # [Mg] [mmol/L]
  Ft            <- F*1e6                                                              # [F] [umol/L]
  magnt         <- magn*56*100/1000/2.5                                               # [Fe3O4] [wt%]
  FeS.Fe        <- FeS*100/2.5*56/1000                                                # [FeS] [wt%]
  FeS2.Fe       <- FeS2*100/2.5*56/1000                                               # [FeS2] [wt%]
  irPat         <- irPa*1000/rho_sed                                                  # amorphous [Fe-P] [umol/g]
  irPbt         <- irPb*1000/rho_sed                                                  # crystalline [Fe-P] [umol/g]
  irPt          <- irPat + irPbt                                                      # total [Fe-P] [umol/g]
  aragt         <- arag*100/2.5*12/1000                                               # aragonite [mmol/cm3]
  calct         <- calc*100/2.5*12/1000                                               # calcite [mmol/cm3]
  rhodt         <- rhod*100/2.5*12/1000                                               # rhodochrosite [mmol/cm3]
  CP_org        <- (orgC1+orgC2+orgC3+orgC4+orgC5+orgC6+J.orgCc/0.2/0.36)/            # organic C/P
                   ((orgC1+orgC2)*y+(orgC3+orgC4+orgC5+orgC6+J.orgCc/0.2/0.36)*yy)
#
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# store and export output
# ---------------------------------------------------------------------------------------------------------------- #
#
  out_df <- as.data.frame(cbind(orgCt,         # organic carbon remineralization
  	                            orgNt,
  	                            orgPt,
  	                            SumNH4t,   
  	                            SumH2PO4t,
#
  	                            O2t,           # terminal electron acceptors
  	                            NO3t,
  	                            MnO2t,
  	                            FeOH3t,
  	                            SO4t,
#
  	                            CH4t,          # reduced solutes
  	                            MnIIt,
  	                            FeIIt,
  	                            SumH2St,
#
  	                            pH,            # charge/acid balance
  	                            SumCO2t,
  	                            Cat,
  	                            Mgt,
  	                            Ft,
#
                                omegaa,        # mineral saturation states
                                omegac,
                                omegaCFA,
#   
                                aragt,        # misc solid phases
                                calct,
                                rhodt,
                                CFA_Pt,
                                irPt,
                                magnt,
                                FeS.Fe,
                                FeS2.Fe,
                                CP_org)
                          ) 
#          
  date <- Sys.Date()
  filename <- paste("SEDCHEM_out",date,".csv",sep=" ")
  write.csv(out_df,file=filename)
#
# ---------------------------------------------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------------------------------------------- #
# plot selected output
# ---------------------------------------------------------------------------------------------------------------- #
#
  par(mfrow=c(1,5))
#
  matplot(orgCt,grid$x.mid,
  	      xlab="C_org [wt%]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(orgCt)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
  matplot(orgNt,grid$x.mid,
  	      xlab="N_org [wt%]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(orgNt)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
  matplot(orgPt,grid$x.mid,
  	      xlab="P_org [wt%]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(orgPt)),
  	      ylim=c(L,0),
          col=1,lty=1,lwd=3)
  matplot(SumNH4t,grid$x.mid,
  	      xlab="[NH4] [wt%]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(SumNH4t)),
  	      ylim=c(L,0),
          col=1,lty=1,lwd=3)
  matplot(SumH2PO4t,grid$x.mid,
  	      xlab="[PO4] [wt%]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(SumH2PO4t)),
  	      ylim=c(L,0),
          col=1,lty=1,lwd=3)
#
  par(mfrow=c(1,5))
#
  matplot(O2t,grid$x.mid,
  	      xlab="[O2] [umol/L]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(O2t)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
  matplot(NO3t,grid$x.mid,
  	      xlab="[NO3] [umol/L]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(NO3t)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
  matplot(MnO2t,grid$x.mid,
  	      xlab="[MnO2] [umol/g]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(MnO2t)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
  matplot(FeOH3t,grid$x.mid,
  	      xlab="[FeOH3] [umol/g]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(FeOH3t)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
  matplot(SO4t,grid$x.mid,
  	      xlab="[SO4] [mmol/L]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(SO4t)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
#
  par(mfrow=c(1,4))
#
  matplot(MnIIt,grid$x.mid,
  	      xlab="[Mn2+] [umol/L]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(MnIIt)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
  matplot(FeIIt,grid$x.mid,
  	      xlab="[Fe2+] [umol/L]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(FeIIt)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
  matplot(SumH2St,grid$x.mid,
  	      xlab="[Fe2+] [umol/L]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(FeIIt)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
  matplot(CH4t,grid$x.mid,
  	      xlab="[CH4] [umol/L]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(CH4t)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
#
  par(mfrow=c(1,5))
#
  matplot(pH,grid$x.mid,
  	      xlab="pH",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(pH)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
  matplot(SumCO2t,grid$x.mid,
  	      xlab="[DIC] [mmol/L]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(SumCO2t)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
  matplot(Cat,grid$x.mid,
  	      xlab="[Ca] [mmol/L]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(Cat)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
  matplot(Mgt,grid$x.mid,
  	      xlab="[Mg] [mmol/L]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(Mgt)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
  matplot(Ft,grid$x.mid,
  	      xlab="[F] [mmol/L]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(Ft)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
#
  par(mfrow=c(1,5))
#
  matplot(calct,grid$x.mid,
  	      xlab="[CaCO3] [mmol/cm3]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(calct)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
  matplot(rhodt,grid$x.mid,
  	      xlab="[MnCO3] [mmol/cm3]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(rhodt)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
  matplot(CFA_Pt,grid$x.mid,
  	      xlab="[CFA] [umol/g]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(CFA_Pt)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
  matplot(FeS.Fe,grid$x.mid,
  	      xlab="[FeS] [wt%]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(FeS.Fe)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
  matplot(FeS2.Fe,grid$x.mid,
  	      xlab="[FeS2] [wt%]",
  	      ylab="depth [cm]",
  	      type="l",
  	      xlim=c(0,max(FeS2.Fe)),
  	      ylim=c(L,0),
  	      col=1,lty=1,lwd=3)
#
# ---------------------------------------------------------------------------------------------------------------- # 
