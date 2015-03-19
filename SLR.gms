* FUND Sea level rise
* AUTHOR: Alexey Shiklomanov

* Abbreviations:
*   eft = e folding time
*   SL = sea level

* "r" functions are draws from distributions:
*   rgamma(shape, scale)
*   rtriangle(lower limit, upper limit)

sets        t   "time"  / t0 * t3 /
            region / "Europe", "USA" /
            sector / "Agriculture", "Service" /
            vars   / "climate_sensitivity", 
                    "temp_eft_alpha", "temp_eft_beta_l", "temp_eft_beta_q", 
                    "eft_rho", "SL_temp_sensitivity"
                    /  
;

parameters          tol_param(region, sector, vars)

                    CS      "climate sensitivity"
                    alpha   "eft_alpha"
                    beta_l  "eft_beta_l"
                    beta_q  "eft_beta_q"
                    rho     "eft_rho"
                    gam     "SL_temp_sensitivity"
;

*** Temperature parameters 
tol_param(region, sector, "climate_sensitivity") = 3.0; 
tol_param(region, sector, "eft_alpha") = -42.7;
tol_param(region, sector, "eft_beta_l") = 29.1;
tol_param(region, sector, "eft_beta_q") = 0.001;
* tol_param(region, sector, "climate_sensitivity") = rgamma(6.48, 0.55);

*** Sea level parameters
tol_param(region, sector, "eft_rho") = 500;
tol_param(region, sector, "SL_temp_sensitivity") = 2;
* tol_param(region, sector, "eft_rho") = rtriangle(250, 1000);
* tol_param(region, sector, "SL_temp_sensitivity") = rgamma(6, 0.4);



*** Simplify parameter names for calls in functions
CS = tol_param(region, sector, "climate_sensitivity");
alpha = tol_param(region, sector, "eft_alpha");
beta_l = tol_param(region, sector, "eft_beta_l");
beta_q = tol_param(region, sector, "eft_beta_q");
rho = tol_param(region, sector, "eft_rho");
gam = tol_param(region, sector, "SL_temp_sensitivity");


*** Sea level rise equations
phi = max(alpha + beta_l*CS + beta_q*CS**2, 1);
temp(t+1) = (1 - 1/phi) * temp(t) + (1/phi) * (CS/(5.35*log(2))) * RF(t);
SL(t+1) = (1 - 1/rho) * SL(t) + gam*temp(t);

