def satvapP(Tair):
   # following constants taken from Vaisala humidity conversion forumulas
   # max error 0.083% for range -20 -> + 50C for constants A, m, Tn
   A = 6.116441     
   m = 7.591386
   Tn = 240.7263

   P_ws = A*10**(m*Tair/(Tair+Tn)) # in hPa
   
   return P_ws
   

def longwave_flux(Tair,SST,Fc,P_licor,H2O_ppm):
   emiss_lw = 0.985       # Emissivity of the Ocean (Dickey et al.)
   sigmaSB = 5.6697e-8    # Stefan-Boltzmann 
   
   TairK = Tair + 273.15
   SSTK = SST + 273.15
   
   ea = (P_licor*(H2O_ppm/1.0e+6))/(1.0 + (H2O_ppm/1.0e+6))
   
   # Berliand calculation of longwave flux
   lwest = -emiss_lw*sigmaSB*(TairK**4)*(0.39-0.05*sqrt(ea))*Fc - 
           4.0*emiss_lw*sigmaSB*(TairK**3)*(SSTK-TairK)

   return lwest
