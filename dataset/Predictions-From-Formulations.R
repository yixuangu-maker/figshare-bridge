#This script analyzes all the fourmulation information from the hydrogels in the 
#comprehensive Richbourg-Peppas dataset and fills in intermediate structural parameter
#and physical property predictions.

#call tidyverse library for table manipulations
library(tidyverse)

#call pracma library for error functions
library(pracma)

#hard-code set working directory
getwd()
homedir <- "C:/Users/nrr25b/OneDrive - Florida State University/NR-Research/Ideas-Other/RichbourgPeppasDatasetFix"
setwd(homedir)

#call in the formulation table
formdata <- read_csv("Form-Richbourg-Peppas-Hydrogel-Dataset.csv")

#function to calculate Stokes hydrodynamic radii from free diffusion coefficients
get_rs <- function (d0, tk, kb, visc) {
  a = tk*kb/(6*pi*d0*visc)*10^6*10^6*10^9
  a
}

#add each solute's rs to the main table
formdata %>% mutate(
  fl00_rs = get_rs(fl00_d0, tk, kb, visc),
  fd04_rs = get_rs(fd04_d0, tk, kb, visc),
  fd20_rs = get_rs(fd20_d0, tk, kb, visc),
  fd70_rs = get_rs(fd70_d0, tk, kb, visc),
  fp05_rs = get_rs(fp05_d0, tk, kb, visc),
  fp20_rs = get_rs(fp20_d0, tk, kb, visc),
  fp40_rs = get_rs(fp40_d0, tk, kb, visc)
) -> formdata

#function to predict swollen polymer volume fraction for each hydrogel formulation)
#uses nonlinear optimization via optim function
get_pspvf <- function(chi1, chi2, pd, v1, mr, nj, fced, fun, ipvf) {
  opter <- function(x, chi1, chi2, pd, v1, mr, nj, fced, fun, ipvf) {
    a = pd*v1/(mr*nj)*(1-fced)*(1-2/fun)*ipvf^(2/3)
    b = (log(1-x)+x+chi1*x^2+chi2*x^3)*x^(-1/3)
    return((a+b)^2)
  }
  res <- optim(ipvf, opter,
    chi1 = chi1,
    chi2 = chi2,
    pd = pd,
    v1 = v1,
    mr = mr,
    nj = nj,
    fced = fced,
    fun = fun,
    ipvf = ipvf,
    method = "Brent",
    lower = 0.001, #lower limit for spvf
    upper = 0.4 #upper limit for spvf
  )
  res[[1]]
}

#add p_spvf to table
formdata %>% rowwise() %>% mutate(
  p_spvf = get_pspvf(
    chi1 = chi1,
    chi2 = chi2,
    pd = pd,
    v1 = v1,
    mr = mr,
    nj = nj,
    fced = fced,
    fun = fun,
    ipvf = ipvf
  )
) ->formdata

#predict all the other relevant parameters!

#p_reference_ratio (p_rr) does not need its own function
#p_rr = ipvf/p_spvf

#function for shear modulus (kPa)
get_g <- function(rgas, tk, fun, fced, pd, nj, mr, ipvf, spvf) {
  a = rgas*tk*(1-2/fun)*(1-fced)*pd/(mr*nj)*ipvf^(2/3)*spvf^(1/3)
  a
}

#function for mesh size (nm)
get_ms <- function(spvf, fun, lbar, cinf, lambda, nj) {
  a = spvf^(-1/3)*sqrt((1-2/fun)*lbar^2*lambda*cinf*nj)
  a
}

#function for mesh radius, has to be applied after ms is calculated
#need to check if this works with mutate, might need to use rowwise()
get_rm <- function(ms, fun) {
  if(fun == 4) {
    a = sqrt(6)/3 * ms
  } else if (fun == 6) {
    a = 0.5 * ms
  } else if (fun == 8) {
    a = sqrt(2)/4 * ms
  } else {
    a = 0.01 * ms #this is an arbitrary low value 
                  #to avoid zero throwing an error
  }
  a
}

#function to predict a diffusion coefficient for a solute in a hydrogel formulation
get_d <- function(d0, rs, rfvw, spvf, rf, rm) {
  a = erf(rfvw/rs)*exp(-1*(rs/rfvw)^3*(spvf/(1-spvf)))
  b = erfc(rfvw/rs)*exp(-1*pi/4*((rs+rf)/rm)^2)
  a = d0*(a+b)
  a
}

#Now, implement adding all the predicted properties
formdata %>% mutate(
  #reference ratio
  p_rr = ipvf/p_spvf,
  p_g = get_g(
    rgas = rgas, 
    tk = tk, 
    fun = fun, 
    fced = fced, 
    pd = pd, 
    nj = nj, 
    mr = mr, 
    ipvf = ipvf, 
    spvf = p_spvf
  ),
  p_ms = get_ms(
    spvf = p_spvf, 
    fun = fun, 
    lbar = lbar, 
    cinf = cinf, 
    lambda = lambda, 
    nj = nj
  ),
  p_rm = get_rm(
    ms = p_ms,
    fun = fun
  ),
  #predict all solute diffusivities
  p_d_fl00 = get_d(
    d0 = fl00_d0, 
    rs = fl00_rs, 
    rfvw = rfvw, 
    spvf = p_spvf, 
    rf = rf,
    rm = p_rm
  ),
  p_d_fd04= get_d(
    d0 = fd04_d0, 
    rs = fd04_rs, 
    rfvw = rfvw, 
    spvf = p_spvf, 
    rf = rf,
    rm = p_rm
  ),
  p_d_fd20 = get_d(
    d0 = fd20_d0, 
    rs = fd20_rs, 
    rfvw = rfvw, 
    spvf = p_spvf, 
    rf = rf,
    rm = p_rm
  ),
  p_d_fd70 = get_d(
    d0 = fd70_d0, 
    rs = fd70_rs, 
    rfvw = rfvw, 
    spvf = p_spvf, 
    rf = rf,
    rm = p_rm
  ),
  p_d_fp05 = get_d(
    d0 = fp05_d0, 
    rs = fp05_rs, 
    rfvw = rfvw, 
    spvf = p_spvf, 
    rf = rf,
    rm = p_rm
  ),
  p_d_fp20 = get_d(
    d0 = fp20_d0, 
    rs = fp20_rs, 
    rfvw = rfvw, 
    spvf = p_spvf, 
    rf = rf,
    rm = p_rm
  ),
  p_d_fp40 = get_d(
    d0 = fp40_d0, 
    rs = fp40_rs, 
    rfvw = rfvw, 
    spvf = p_spvf, 
    rf = rf,
    rm = p_rm
  )
) ->formdata

#reformat column names for public sharing of the table
formdata %>% rename(
  "Formulation ID" = form,
  "Polymer" = poly,
  "Initial PVF" = ipvf,
  "Nj" = nj,
  "Jun. Fun." = fun,
  "FCED" = fced,
  "Precursor MW [g/mol]" = mw,
  "Chi-1" = chi1,
  "Chi-2" = chi2,
  "Dens [g/mL]" = pd,
  "V1 [mol/mL]" = v1,
  "Mr [g/mol]" = mr,
  "Gas Const." = rgas,
  "Temp [K]" = tk,
  "Avg. bond length [nm]" = lbar,
  "Bonds per RU" = lambda,
  "C-inf" = cinf,
  "Fiber radius [nm]" = rf,
  "Rad. free vol. water" = rfvw,
  "Boltz. const." = kb,
  "viscosity" = visc,
  "fl00_d0 [um^2/s]" = fl00_d0,
  "fl00_rs [nm]" = fl00_rs,
  "p_g [kPa]" = p_g,
  "p_ms [nm]" = p_ms,
  "p_rm [nm]" = p_rm,
  "p_d_fl00 [um^2/s]" = p_d_fl00
) -> formdata

#Write to a csv
write.csv(formdata, "Pred-RichPepp-Dataset.csv")