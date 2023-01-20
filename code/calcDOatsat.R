# Calculate dissolved oxygen solubility at given temperature and salinity

calcDOatsat_GG <- function(t,s) {
  
  # Input:                         
  # From Garcia and Gordon 1992 (https://doi.org/10.4319/lo.1992.37.6.1307)
  # note: there is a typo in the equation in paper, correct equation used here
  #   t = temperature (ºC)
  #   s = salinity    (PSU)
  #
  # Output:
  #   DO_sat_GG = [oxygen] in µmol/kg
  #
  
  
  A0 =  5.80818;
  A1 =  3.20684;
  A2 =  4.11890;
  A3 =  4.93845;
  A4 =  1.01567;
  A5 =  1.41575;
  
  B0 = -7.01211e-03;
  B1 = -7.25958e-03;
  B2 = -7.93334e-03;
  B3 = -5.54491e-03;
  
  C0 = -1.32412e-07;
    
    #   Calculate Ts from T (ºC)  
    
    Ts = log((298.15 - t) / (273.15 + t));
    
    #   Calculate O2 solubility in µmol O2/kg in a few steps (expanded polynomial)
    
    A = ((((A5*Ts + A4)* Ts + A3)* Ts + A2)* Ts + A1)* Ts + A0;
    
    B = ((B3*Ts + B2)* Ts + B1)* Ts + B0;
    
    O2 = exp(A + s*(B + s*C0));
  
}
