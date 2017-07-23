; NAME: 
;      EXOMOP- EXOplanet MOdeling Package
;      See Pearson, Turner, Sagan. 2014. New Astronomy. Volume 27. 102-110, arXiv.1310.5397
;          Turner et al. 2016, MNRAS, arXiv.1603.02587 (Version 7)
;
;AUTHOR: 
;       Jake Turner, University of Virginia, jt6an@virginia.edu
;
;PURPOSE: 
; This program inputs a light curve of an exoplanet and fits a model fit to it.  
; This program performs a Levenberg-Marquardt least-squares fit to the data to find the best fit model
; Finds the best model fit using the Mandel & Agol (2002) prescription. 
; This program uses the Bootstrap Monte Carlo technique to find errors on the transit parameters. 
; This program also performs a Differential Evolution Markov Chain Monte Carlo (MCMC) fit to find parameters values and errors (ter Braak, 2006)
; Note: This program accounts for red noise using the Prayer Bead Method (Residual-Permutation Method) 
; 	There are two methods of the Prayer Bead Method used:      
;	1.) The method by Fernandez et al. 2009 but has updated to use the error bars from the residual (Like the Bootstrap Monte Carlo method). Assumes a Gaussian Distribution  
;	2.) The method in Todorov et al 2012 (APJ,746,111). Does not assume a Gaussian Distribution. 
; This program also uses the Time-Averaging Method (Beta Method) to further account for red noise (Pont et al 2006, 373,231)
;Also, uses the Wavelet method to account for red noise (Carter & Winn 2009, APJ, 704, 51) 
;
;WHAT IT DOES: 
; 1.)  Fits a linear or quadratic function to the OoT baseline simultaneously with the Mandel & Agol (2002) model for both the LM and MCMC methods
; 2.)  Fits your transit data using Levenberg-Marquardt (LM) non-linear least-squares fit to a Mandel & Agol (2002) model
; 3.)  Finds robust errors using the Bootstrap Monte Carlo technique 
; 4.)  Finds the asymmetry of the transit about the calculated mid-transit time (subtracts the transit by its mirror image)
; 5.)  Uses the Prayer Bead Method (Residual-Permutation Method) to access the influence of red noise on the transit and calculates new errors (Using LM residuals)
; 6.)  Uses the Time-Averaging Method (Beta Method) to further access the influence of red noise (Finds sigma_red_noise and sigma_white_noise) (Using LM residuals)
; 7.)  Uses the Wavelet Method to access the influence of red noise (finds red noise and white noise) (Using LM residuals)
; 8.)  Fits your de-trended transit data using a Differential Evolution Markov Chain Monte Carlo (MCMC) fit (ter Braak, 2006) and finds errors. 
;             Program stops when chains are well-defined using the Gelman-Rubin Statistic (Ford 2006)
; 9.)  Uses the Prayer Bead Method (Residual-Permutation Method) to access the influence of red noise on the transit and calculates new errors (Using MCMC residuals)
; 10.)  Uses the Time-Averaging Method (Beta Method) to further access the influence of red noise (Finds sigma_red_noise and sigma_white_noise) (Using MCMC residuals)
; 11.) Uses the Wavelet Method to access the influence of red noise (finds red noise and white noise) (Using MCMC residuals)
; 
;INPUTS: EXTERNAL FILES:
;       1.) Transit File (This file needs to contain: Time (HJD), Normalized Flux, Error)
;       2.) Input File (See example input file in folder and in README) 
;	NOTE (limitation): Time of transit needs to be inputed with only the ones place and the decimals (e.g. change 2454218.90003 to 8.90003) 
;
;OUTPUTS: 
;IN MAIN FOLDER: 
;Input_Value.dat       - A file containing the input values and what values were fixed. 
;plot_LM.ps 	         - Plot of the Best-fit model, data, and residuals from the Levenberg-Marquardt (LM) non-linear least-squares fit	
;plot_MCMC.ps 	       - Plot of the Best-fit model, data, and residuals from the Differential Evolution Markov Chain Monte Carlo (MCMC) fit 
;Best-Fit_LM.dat       - Time, Phase, original transit,detrended transit, error, best-fit model, and residuals from LM method 
;Best-Fit_MCMC.dat     - Time, Phase, original transit,detrended transit, error, best-fit model, and residuals from MCMC method
;FINAL_EVERYTHING.dat  - 
;          a. RMS of residuals from LM Method
;				   b. RMS of residuals from MCMC Method
;				   c.  Differential Evolution Markov Chain Monte Carlo Parameter Values and Errors 
;					  1. White,noise & red noise (from Wavelet and Time Av), Beta from Time-Averaging Method, Parameter Betas from Residual-Permutation Method, Beta from Wavelet
;				   d.  Monte Carlo Levenberg-Marquardt Parameter Values and Errors 
;				    1. White,noise & red noise (from Wavelet and Time Av), Beta from Time-Averaging Method, Parameter Betas from Residual-Permutation Method, Beta from Wavelet
;				   e. Other Useful Values
;					  1. Duration 
;					  2. Cadence
;					  3. Out of Transit STDEV
;					  4. Reduced Chi-Squared
;					  5. STDEV of Assmetry 
;
; IN MCMC FOLDER: 
;Final_Parameters_MCMC.dat        - Final parameters and errors from MCMC fit
;beta_plot.ps                     - Plot of relative rms vs bin size (line is with no red noise Original_RMS/sqrt(bin size) ) see Knutson et al 2012 (APJ, 754, 22) 
;gauss_prayer_mcmc.ps             - Plot of the distribution of fitted parameters by the Prayer Bead Method from the Monte Carlo simulations 
;Final_parameters_prayer_mcmc.dat - 
;                                    a. The fitted parameters of the planet with error bars derived from the Prayer Bead Red Noise method
;                                      (rp/R*,a/R*,Inclination, Mid-Transit Time, Linear and Quadratic LD) 
;                              	     b. Beta Values of fitted parameters (Red noise error/white noise error)
;Final_parameters_Beta_mcmc.dat  -   a. Sigma White Noise and error
;				                             b. Sigma Red Noise and errors
;				                             c. Beta Value calculated from the Time-Averaging Red Noise Technique 
;monte_value_prayer_mcmc.dat     - Each Monte Carlo run from the Prayer Bead Method (Best-fit model-monte carlo value)  
;probability_plot.ps              -Plot of Probablity Distribution Functions (PDFs)
;covar.ps                         -plot of covariances
;monte_value_prayer_method2.dat   -Each Monte Carlo run from the Prayer Bead Method 2 (Best-fit model-monte carlo value)  
;variablename_prayer_histogram.ps - histogram of monte_value_prayer_method2.dat for each parameter. Red-dashed lines are the error bars derived. 
;Gelman_Rubin.dat                 - The Gelman-Rubin Statitic and Number of independent draws from the MCMC analysis
;
; IN LM FOLDER:
;Final_parameters_prayer.dat         - 
;              a. The fitted parameters of the planet with error bars derived from the Prayer Bead Red Noise method
;                                    (rp/R*,a/R*,Inclination, Mid-Transit Time, Linear and Quadratic LD) 
;              b. Beta Values of fitted parameters (Red noise error/white noise error)
;Final_parameters_Beta.dat          - 
;              a. Sigma White Noise and error
;				       b. Sigma Red Noise and errors
;				       c. Beta Value calculated from the Time-Averaging Red Noise Technique 
;Final_parameters.dat -
;              a. The fitted parameters of the planet with error bars from the Levenberg-Marquardt Monte Carlo method
;                                    (rp/R*,a/R*,Inclination, Mid-Transit Time, Linear and Quadratic Limb Darkening) 
;			         b. Duration of transit in days and mins (with errors)
;			         c. Cadence of observations (days and secs)
;			         d. Out of Transit RMS 
;			         e. RMS of residuals 
;           	 f. Reduced Chi-Squared from best fit-model
;			         g. Lowest Monte Carlo Reduced Chi-Squared 
;			         h. STDEV of the transit subtracted by the mirror image of itself (from transitsymmetry.pro)
;				       i. Baseline Intercept calculated from Detrend_linfit.pro
;              j. Baseline slope calculated from Detrend_linfit.pro
;Best-Fit_LM.dat                     - Time, Phase, original transit,detrended transit, error, best-fit model, and residuals
;plot_LM.ps 	                       - Plot of the Best-fit model, data, and residuals	
;monte_value2.dat                    - Each Monte-Carlo run (Best-fit model-Monte Carlo value)  
;monte_value.dat                     - Each Monte-Carlo run (Monte Carlo values and Chi-Squared)
;gauss.ps 	                         - Plot of the distribution of fitted parameters from the Monte Carlo simulations
;Asymmetry_Plot.ps                   - 4 Plots showing the asymmetry process (the 4th plot shows the subtraction) 
;symmetry_test                       - time and flux used for the asymmetry test 
;monte_value_prayer.dat              - Each Monte Carlo run from the Prayer Bead Method (Best-fit model-monte carlo value)  
;gauss_prayer.ps                     - Plot of the distribution of fitted parameters by the Prayer Bead Method from the Monte Carlo simulations 
;beta_plot.ps                        - Plot of relative rms vs bin size (line is with no red noise Original_RMS/sqrt(bin size) ) see Knutson et al 2012 (APJ, 754, 22) 
;beta_value.dat                      - The Values plotted beta_plot.ps (Bin size,Relative RMS,Original_RMS/sqrt(bin size)) 
;monte_value_prayer_method2_mcmc.dat -  Each Monte Carlo run from the Prayer Bead Method 2 (Best-fit model-monte carlo value)  
;variablename_prayer_histogram.ps    - histogram of monte_value_prayer_method2_mcmc.dat for each parameter. Red-dashed lines are the error bars derived. 
;Wavelet_red_noise_LM.dat            - The Red and White Noise calcualted from the Wavelet Method 
;
; LIMITATIONS:
;      - Does not fit for period
;	     - Time of transit needs to be inputed with only the ones place and the decimals (e.g. change 245551.90003 to 1.90003) 
;
; ENSURING YOU HAVE RELIABLE VALUES AND ERROR BARS (See README): 
;	0. Check to make sure that the gauss.ps and gauss_prayer.ps plots are gaussian. If they are not, run your model with more/less fixed parameters 
;          A. Also, pay attention to the Bayesian information criterion value (lower values are usually better). 
;	1. The file FINAL_EVERYTHING.dat includes all parameter fits, respective errors, both red noise calculation results and residuals. 
;	2. It is recommend that you choose the method (LM or MCMC) that minimized the residuals
;	3. Multiple your error bars from FINAL_EVERYTHING.dat for each parameter by the largest beta (either the global value or individual parameter betas)
;	4. To get your final error bars multiple the error bars from step 3 by sqrt(Chi-Squared), where Chi-Squared is the Reduced Chi-Squared in       
;	     FINAL_EVERYTHING.dat file. Only do this if the Reduced Chi-Squared is greater than one. 
;
;
; SUBROUTINES CALLED (provided in /sub): beta.pro, bin.pro, beta_model.pro,model_fit.pro, 
;                                        Detrend_quadfit.pro, Detrend_linfit.pro,Useful_Values.pro,
;					                               transitsymetry.pro,makesymetric.pro,exomop_mcmc.pro,
;                                        exomop_chi2.pro,exomop_prayer.pro,exomop_wavelet.pro
;		                                     Subroutines from EXOFAST, Carter, Mandel & Agol (2002)
;
; MODIFICATION HISTORY 
; 2013/02/13   Jake Turner (LPL) V1.0 
; 	           Created
; 2013/02/17   Jake Turner (LPL) V2.0
;              Added linear detrending of Out of Transit baseline 
; 2013/02/20   Jake Turner (LPL) V2.1
;              Added quadratic detrending of Out of Transit baseline
; 2013/03/14   Jake Turner (LPL) V3.0 
;              Added the Time-Averaging Red Noise Method (Beta Method) to the code
;	              User now can choose whether to run Linear/Quadratic Detrending 
;2013/03/19    Jake Turner (LPL) V3.1
;              Added the option to include eccentricity and omega in the fit 
;2013/03/23     Jake Turner (LPL) V4.0 
;              Now includes an MCMC fit (adapted from Eastman, Gaudi & Agol, 2013,PASP, 125, 83)
;2013/04/29	Jake Turner (LPL) V5.0 
;		         Added another version of the residual permutation method (Todorov, 2012, APJ,746,111)
;		         Can now input whether to run Symmetry test
;		         Now has the option to use the full mcmc model with the mcmc residuals in the prayer bead method. (Very slow)
;2014/01/06     Jake Turner (UVA) V5.1
; 		          Fixed bug in BIC that affects the absolute value
;		             Period is now double precision
;2014/01/20  ake Turner (UVA) V6.0
;     		   The Wavelt Red noise method was added (Carter & Winn APJ, 704, 51)
;2014/03/29  Jake Turner (UVA) V6.1
;		         Time Averaging method was updated 
;2015/01/23  Jake Turner (UVA) V7.0
;              Changed input into text file
;              Implemented Simultaneous fitting of parameters and baseline functions	
;		           Updated time Averaging method 
;		           
pro EXOMOPv70,input_file


;------------------------------------------------
;------------------Version Number----------------
;------------------------------------------------
version = 7.0

;------------------------------------------------
;------------------Compile Subprograms-----------
;------------------------------------------------
RESOLVE_ALL,RESOLVE_EITHER=RESOLVE_EITHER

;------------------------------------------------
;-----Default Run--------------------------------
;------------------------------------------------
run_prayer      = 1		; Prayer-Method Red noise fit (1=Yes) & (0=NO)               (Default=1)
run_beta        = 1 	; Beta Red noise calculation  (1=Yes) & (0=NO) 			         (Default=1)
run_symmetry    = 1		; Symmetry calculation (Transit - Mirror image of itself)    (Default=1)  (0=Yes) & (1=NO)
run_MCMC        = 1 	; Run MCMC fit 		      (0=Yes) & (1=NO)                     (Default=1)
run_LM 	        = 0		; Run LM fit 		      (0=Yes) & (1=NO)                       (Default=0)
run_long_prayer = 1 	; Run prayer bead with full MCMC (0=Yes) & (1=NO) 		       (Default=1)  (Will take many hours!)
run_wavelet     = 1   ; Run Wavelet Red noise (1=YES) & (0=NO)                     (Default=1)

Common baseline_block,run_baseline,egress,ingress
;------------------------------------------------
;-----READ USER INPUTS --------------------------
;------------------------------------------------

;Read input file
read_input,input_file,transitfile=transitfile,foldername=foldername,run_LM=run_LM,run_MCMC=run_MCMC,maxsteps=maxsteps,Ninput=N,nchains=nchains,phase_indicator=phase_indicator, $
  rp=rp,a_rstar=a_rstar,incl=incl, mid_Transit=mid_Transit,period=period,ecc=ecc,omega=omega,linlimb=linlimb,quadlimb=quadlimb,                      $
  parinfo_rp=parinfo_rp,parinfo_ar=parinfo_ar,parinfo_inc=parinfo_inc,parinfo_mid=parinfo_mid,parinfo_lin=parinfo_lin,parinfo_quad=parinfo_quad,     $
  info_baseline=info_baseline,run_baseline=run_baseline,ingress=ingress,egress=egress,run_symmetry=run_symmetry,run_renormalize=run_renormalize
  
  file_mkdir, foldername
  file_mkdir, foldername+'/MCMC'
  file_mkdir, foldername+'/LM'

  FILE_COPY, [input_file], foldername    ;copy input file

;print user input file
openw,9,'Input_Values.dat'
printf,9,'Input Values into Model'
printf,9,'************************************'
printf,9, 'Transit File:',transitfile,FORMAT='(A13,2x,A30)'
printf,9,'Run MCMC (0=YES,1=NO)', run_MCMC
printf,9,'Run LM   (0=YES,1=NO)', run_LM
printf,9,'************************************'
printf,9, 'Number of Monte Carlo Simulations to Run: ', N
If run_mcmc eq 0 then begin 
	printf,9,'Maximum number of MCMC steps: ', maxsteps
	printf,9,'Number of Chains in MCMC fit: ', nchains
endif 
printf,9,'************************************'
printf,9,'Time in Phase? (0=YES,1=NO)', phase_indicator
printf,9, 'Rp/R*:',rp,FORMAT='(A6,9x,F9.7)
printf,9, 'a/R*:',a_rstar,FORMAT='(A5,10x,F9.6)
printf,9, 'Inclination:',incl,FORMAT='(A12,1x,A20)'
printf,9, 'Mid-Transit Time:   ', mid_Transit
printf,9, 'Period:', period,FORMAT='(A7,19x,F12.9)'
printf,9, 'Eccentricity:',ecc, FORMAT='(A13,2x,F10.7)'
printf,9, 'Omega (Degrees):',omega,FORMAT='(A16,2x,F9.4)'
printf,9, 'Linear Limb Darkening Coefficient:       ', linlimb
printf,9, 'Quadratic Limb Darkening Coefficient     ', quadlimb
printf,9,'************************************'
printf,9,'Fit for Rp/R*? (0=YES,1=NO)', parinfo_rp 
printf,9,'Fit for a/R*?  (0=YES,1=NO)', parinfo_ar 
printf,9,'Fit for I?     (0=YES,1=NO)', parinfo_inc 
printf,9,'Fit for Tmid?  (0=YES,1=NO)', parinfo_mid
printf,9,'Fit for Linear?(0=YES,1=NO)', parinfo_lin
printf,9,'Fit for Quad?  (0=YES,1=NO)', parinfo_quad
printf,9,'************************************'
printf,9,'Fit Baseline?  (0=YES,1=NO)', info_baseline
If info_baseline eq 0 then begin
	printf,9,'Function of baseline fit (1=linear, 2=quadratic):', run_baseline
	printf,9,'Number of Ingress Pts', round(ingress)
	printf,9,'Number of Egress Pts', round(egress)
endif 
printf,9,'Fit for Symmetry? (0=YES,1=NO)', run_symmetry  
printf,9,'Renormalize the light curve (0=YES, 1=NO):', run_renormalize
close,9 
file_move,['Input_Values.dat'],foldername
FILE_COPY, [transitfile], foldername		;copy transit file

;------------------------------------------------
;-----READ TRANSIT FILE--------------------------
;------------------------------------------------
COMMON chi2_block,transit		;set common

readcol,transitfile,time,flux_transit,error
size = n_elements(time) 

;save original transit
flux_transit_or = flux_transit

;Update Time if it is in phase 
If phase_indicator eq 0 then begin
  time = time*period		;only do if you input time is in phase
endif


;*************************   
;     Find baseline
;*************************
egress  = round(egress)
ingress = round(ingress)

time_base 		      = fltarr(size)
flux_transit_base 	= fltarr(size)
error_base		      = fltarr(size)
count               = 0 

for i=0,size-1 do begin
  If ((i lt ingress) OR (i gt egress)) then begin 
      time_base[i] 		      = time[i] 
      flux_transit_base[i]	= flux_transit[i]
      error_base[i] 		    = error[i]	
      count                 = count + 1
  endif
endfor
print,"Number of data points: ",count
time_base2 		      = fltarr(count)
flux_transit_base2 	= fltarr(count)
error_base2		      = fltarr(count)

time_base2[0:ingress-1]         	   = time_base[0:ingress-1]
time_base2[ingress:count-1] 	 	     = time_base[egress+1:size-1]
flux_transit_base2[0:ingress-1] 	   = flux_transit_base[0:ingress-1]
flux_transit_base2[ingress:count-1]	 = flux_transit_base[egress+1:size-1]
error_base2[0:ingress-1]	  	       = error_base[0:ingress-1]
error_base2[ingress:count-1] 		     = error_base[egress+1:size-1]

run_baseline = round(run_baseline)

;*************************
;       renormalize
;**************************
If run_renormalize eq 0 then begin
  base                  = fltarr(count)
  base[0:ingress-1]     = flux_transit[0:ingress-1]
  base[ingress:count-1] = flux_transit[egress+1:size-1]

  average      = MEAN(base)
  flux_transit = flux_transit/average    ;update
  error        = error/average
endif 

;If info_baseline eq 0 then begin
  ;----------------------------------------------------------------------------------------------------
  ;----------------------------------De-trend Baseline of Transit--------------------------------------
  ;----------------------------------------------------------------------------------------------------

 ;Linear Detrender
If run_baseline eq 1 then begin
	Detrend = detrend_linfit(time_base2,flux_transit_base2,error_base2)		;run detrender
endif 

;Quadratic Detrender
If run_baseline eq 2 then begin
  Detrend = detrend_quadfit(time_base2,flux_transit_base2,error_base2)		;run detrender
endif 


;------------------------------------------------
;-----Updated Structure with detrended transit---
;------------------------------------------------
;save structure of transit
trandata =  create_struct('time',time,'flux',flux_transit,'err',error)
transit=trandata

;------------------------------------------------
;-----Number of Fit Parameters----------
;------------------------------------------------

N_fit = 12

If run_LM eq 0 then begin
  ;------------------------------------------------
  ;-----Number of Monte Carlo Simulations----------
  ;------------------------------------------------
  N = N  ; Monte Carlo Simulations (make sure it is > than 1000)
  	    ; 10000 takes about 10 minutes to run 
  	    ; 100000 takes about 50 minutes to run
  
  ;------------------------------------------------
  ;-----initialize parameters----------------------
  ;------------------------------------------------
  rp_monte 	     = fltarr(N)
  rp_monte2	     = fltarr(N)
  a_rstar_monte	 = fltarr(N)
  a_rstar_monte2 = fltarr(N)
  incl_monte     = fltarr(N)
  incl_monte2    = fltarr(N)
  mid_monte      = fltarr(N)
  mid_monte2     = fltarr(N)
  lin_monte2		 = fltarr(N)
  quad_monte2	   = fltarr(N)
  lin_monte		   = fltarr(N)
  quad_monte		  = fltarr(N)
  fit_start      = dblarr(N_fit)
  
  flux_transit2  = fltarr(size)
  chi			       = fltarr(N)
  
  ;MCMC
  rp_monte_mcmc 	    	= fltarr(N)
  rp_monte2_mcmc	    	= fltarr(N)
  a_rstar_monte_mcmc	= fltarr(N)
  a_rstar_monte2_mcmc	= fltarr(N)
  incl_monte_mcmc    	= fltarr(N)
  incl_monte2_mcmc    = fltarr(N)
  mid_monte_mcmc     	= fltarr(N)
  mid_monte2_mcmc   	= fltarr(N)
  lin_monte2_mcmc		 = fltarr(N)
  quad_monte2_mcmc	 = fltarr(N)
  lin_monte_mcmc		 = fltarr(N)
  quad_monte_mcmc		 = fltarr(N)
  fit_start_mcmc     = dblarr(N_fit) 
  flux_transit2_mcmc = fltarr(size)
  chi_mcmc	     	   = fltarr(N)
  
  ;Prayer Bead Method 2
  rp_monte_prayer2	    	= fltarr(size)
  a_rstar_monte_prayer2		= fltarr(size)
  incl_monte_prayer2	  	= fltarr(size)
  mid_monte_prayer2	    	= fltarr(size)
  lin_monte_prayer2	    	= fltarr(size)
  quad_monte_prayer2	  	= fltarr(size)
  flux_transit3 		    	= fltarr(size)
  
  rp_monte_prayer2_mcmc      = fltarr(size)
  a_rstar_monte_prayer2_mcmc = fltarr(size)
  incl_monte_prayer2_mcmc    = fltarr(size)
  mid_monte_prayer2_mcmc     = fltarr(size)
  lin_monte_prayer2_mcmc     = fltarr(size)
  quad_monte_prayer2_mcmc    = fltarr(size)
  flux_transit4              = fltarr(size)
  
  ;------------------------------------------------
  ;-----Output Files------------------------------
  ;------------------------------------------------
  openW,1,'monte_value.dat',width=250
  openw,12,'monte_value2.dat'
  openw,5,'Best-Fit_LM.dat' , width=250
  
  ;------------------------------------------------
  ;-----Setup of Initial Guess----------------------
  ;------------------------------------------------
  incl = (!dPI/180.d)*incl	;Convert to radians !IMPORTANT
  omega = (!dPI/180.d)*omega ;Convert to radians !IMPORTANT
  
  fit_start[0] 	= rp
  fit_start[1] 	= a_rstar
  fit_start[2] 	= incl
  fit_start[3] 	= mid_transit 
  fit_start[4] 	= linlimb
  fit_start[5]  = quadlimb
  fit_start[6] 	= period
  fit_start[7] 	= ecc
  fit_start[8]  = omega
  
  ;No baseline fit
  If info_baseline eq 1.0 then begin
    fit_start[9]   = 0.0d      ;a  ;ax^2 + bx + c
    fit_start[10]  = 0.0d      ;b
    fit_start[11]  = 1.0d      ;c
  endif
  
  ;Linear
  If run_baseline eq 1 then begin
    fit_start[9]   = 0.0d            ;a 
    fit_start[10]  = Detrend[1]      ;b
    fit_start[11]  = Detrend[0]      ;c
  endif
  
  ;Quadratic Detrender
  If run_baseline eq 2 then begin
    fit_start[9]   = Detrend[2]      ;a 
    fit_start[10]  = Detrend[1]      ;b
    fit_start[11]  = Detrend[0]      ;c
  endif
  
  
  ;------------------------------------------------
  ;-----Create parameter info structure------------
  ;------------------------------------------------
  ;parinfo = replicate({fixed:0,MPSIDE:0,limited:[0,0], limits:[0.D,0.D]},9)	
  parinfo = replicate({fixed:0,MPSIDE:0,limited:[0,0], limits:[0.D,0.D]},N_fit)
  
  
  ;set limits on Rp/R* 
  parinfo[0].limited(0) = 1.0d
  parinfo[0].limits(0)  = 0.001d
  parinfo[0].limited(1) = 1.0d
  parinfo[0].limits(1)  = 0.4d
  
  ;set limits on a/R* 
  parinfo[1].limited(0) = 1.0d
  parinfo[1].limits(0)  = 1.5d
  
  ;flag inclinatoin
;  parinfo[2].limited(0) = 1.0d
;  parinfo[2].limited(1) = 1.0d
;  parinfo[2].limits(0) = 75.0d   ;lower
;  parinfo[2].limits(1) = 95.0d  


  
  ;Flag
  ;parinfo[1].limits(0)  = 7.88d
  ;parinfo[1].limited(1) = 1.0d
  ;parinfo[1].limits(1)  = 7.12d

  
  ;set limits for Linear Limb Darkening 
  parinfo[4].limited(0) = 1.0d
  parinfo[4].limits(0) = 0.0d
  parinfo[4].limited(1) = 1.0d
  parinfo[4].limits(1) = 2.0d
  
  ;set limits for Quadratic Limb Darkening 
  parinfo[5].limited(0) = 1.0d
  parinfo[5].limits(0) = -1.0d
  parinfo[5].limited(1) = 1.0d
  parinfo[5].limits(1) = 1.0d
  
  parinfo[0].fixed = byte(parinfo_rp)		; Update Fit for Rp/R*
  parinfo[1].fixed = byte(parinfo_ar)		; Update Fit for a/R*
  parinfo[2].fixed = byte(parinfo_inc)		; Update Fit for Inclination
  parinfo[3].fixed = byte(parinfo_mid)		; Update Fit for Mid-Transit
  parinfo[4].fixed = byte(parinfo_lin)		; Update Fit for linear limb darkening
  parinfo[5].fixed = byte(parinfo_quad)		; Update Fit for Quad limb darkening
  parinfo[6].fixed = 1			 	            ; Update Fit for Period (Never Fit!)
  parinfo[7].fixed = 1			 	            ; Update Fit for eccentricity (Never Fit!)
  parinfo[8].fixed = 1			 	            ; Update Fit for omega (Never Fit!)
  
  ;------------------------------------------------
  ;Baseline Fit
  ;------------------------------------------------
  ;No baseline fit
  If info_baseline eq 1.0 then begin
    parinfo[9].fixed  =  1                    ; Update Fit for a ; ax^2 + bx +c; don't run
    parinfo[10].fixed =  1                    ; Update Fit for b
    parinfo[11].fixed =  1                    ; Update Fit for c
  endif
  
  
  ;Linear
  If run_baseline eq 1 then begin
    parinfo[9].fixed  =  1                    ; Update Fit for a ; ax^2 + bx +c; don't run
    parinfo[10].fixed =  0                    ; Update Fit for b
    parinfo[11].fixed =  0                    ; Update Fit for c
  endif
  
  ;Quadratic Detrender
  If run_baseline eq 2 then begin
    parinfo[9].fixed  =  0                    ; Update Fit for a
    parinfo[10].fixed =  0                    ; Update Fit for b 
    parinfo[11].fixed =  0                    ; Update Fit for c 
  endif 
  
  
  for i=0,N_fit-1 do begin 		;set parameter structure 
    parinfo[i].mpside = 2		;  2 - two-sided derivative (f(x+h) - f(x-h))/(2*h) 
  endfor 				             ; the sidedness of the finite difference when computing numerical derivatives. 
  
  						
  ;---------------------------------------------------------------------------------------------------------
  ;--------------------------------------------Find Best Fit (Initial Guess)--------------------------------
  ;---------------------------------------------------------------------------------------------------------
  
  	flux_model = model_fitv2(time,fit_start)    ;produce model with initial parameters
  
  ;Perform Levenberg-Marquardt least-squares fit
  	FitValue   = mpfitfun('model_fitv2',time,flux_transit,error,fit_start,PERROR=perror,PARINFO=parinfo,/AUTODERIVATIVE,/quiet,BESTNORM= chi_squared2,DOF=DOF)
  	Best_fit   = FitValue
  	best_model = model_fitv2(time,Best_fit) 	;produce model with fitted parameters 
  
  ;find Reduced Chi_squared
  	print, chi_squared2/DOF
  	limit = chi_squared2/DOF
  
  ;repeat to get an initial guess that produces a reduced_chi_squared lt original. Attempts to see if first guess was stuck in local minimum. 
  repeat begin 
  	flux_model = model_fitv2(time,fit_start)	;produce model with initial parameters
  	FitValue   = mpfitfun('model_fitv2',time,flux_transit,error,fit_start,PERROR=perror,PARINFO=parinfo,/AUTODERIVATIVE,/quiet,BESTNORM= chi_squared2,DOF=DOF)
  	Best_fit   = FitValue
  	best_model = model_fitv2(time,Best_fit)  ;produce model with fitted parameters 
  
  ;find Chi_squared
  	print, chi_squared2/DOF
  	fit_start[0] = fit_start[0] + 0.001d		; small change to recalibrate LM fit 
  
  endrep until chi_squared2/DOF lt limit			; If Reduced Chi_squared is lt original  stop					
  
  	
  ;--------------------------------------
  ;----Output Best Fit Model--------------
  ;--------------------------------------
  
  best_model = model_fitv2(time,Best_fit)	
  
  best_mode_save = best_model 
  	
  ;********************************
  ;  Correct Transit/Model for baseline
  ;********************************
  
  ;Renormalize Linear
  If run_baseline eq 1 then begin
       flux_transit = flux_transit_or*(best_fit[11] + best_fit[10]*time) ;detrend
       error        = error[i]*(best_fit[11] + best_fit[10]*time )
       best_model  = best_model*(best_fit[11] + best_fit[10]*time )
  endif
  
  ;Renormalize Quadratic
  If run_baseline eq 2 then begin
     for i=0,size-1 do begin
       flux_transit[i] = flux_transit_or[i]*(best_fit[11] + best_fit[10]*time[i]+best_fit[9]*time[i]*time[i]) ;detrend
       error[i]        = error[i]*(best_fit[11] + best_fit[10]*time[i]+best_fit[9]*time[i]*time[i]) 
       best_model[i]   = best_model[i]*(best_fit[11] + best_fit[10]*time[i]+best_fit[9]*time[i]*time[i])
     endfor
  endif
  	
  	
  ;-------------------------------------------
  ;--Bayesian Information Criterion ---------
  ;------------------------------------------
  ;BIC = X^2 + k*ln(N), where k is the degrees of freedom and N is the number of data points (BIC; Schwarz 1978)
  ;This technique has been used to access over-fitting (Gibson et al, 2010, MNRAS, 404, L114-L118)
  
  fit_free = n_elements(where(parinfo.fixed eq 0))   ; find the number of free parameters fit
  BIC = chi_squared2  + fit_free*ALOG(size)	
  	
  print,'BIC: ',BIC	
  ;------------------------------
  ;--Reduced Chi-Squared---------
  ;------------------------------
  chi_reduced = chi_squared2/DOF	 ; chi/(DOF)
  
  	openw,8,'Initial_parameters'
  	printf,8,Best_fit[0], PERROR[0]
  	printf,8,Best_fit[1], PERROR[1]
  	printf,8,Best_fit[2]*(180/!PI), PERROR[2]*(180/!PI)	; convert to Degrees
  	printf,8,Best_fit[3], PERROR[3]
  	printf,8,Best_fit[4], PERROR[4]
  	printf,8,Best_fit[5], PERROR[5]
  	printf,8,Best_fit[9], PERROR[9]
  	printf,8,Best_fit[10], PERROR[10]
  	printf,8,Best_fit[11], PERROR[11]
  	close,8
  	file_move,['Initial_parameters'],foldername+'/LM
  
  	fit_start[0] 	= rp			;reset fit_start[0] 
  
  
  ;--------------------------------------
  ;----Save Best fit for MCMC------------
  ;--------------------------------------
  
  bestpars = fltarr(N_fit) 
  bestpars[0] = Best_fit[0]	;rp
  bestpars[1] = Best_fit[1]	;ap
  bestpars[2] = Best_fit[2]	;incl	
  bestpars[3] = Best_fit[3]	;mid
  bestpars[4] = Best_fit[4]	;lin limb
  bestpars[5] = Best_fit[5]	;quad limb
  bestpars[6] = fit_start[6] 	;period
  bestpars[7] = fit_start[7] 	;ecc
  bestpars[8] = fit_start[8]    	;omega
  bestpars[9] = Best_fit[9]      ;a
  bestpars[10] = Best_fit[10]      ;b
  bestpars[11] = Best_fit[11]       ;c
  
  ;------------------------------------------------
  ;-----Calculate Residuals and Phase--------------
  ;------------------------------------------------
  
  res = fltarr(1,size) 
  res = flux_transit - best_model
  
  phase = ( time - FitValue[3] ) / period
  
  ;-----------------------------------------------
  ;-----Print Transit and Best fit Model----------
  ;-----------------------------------------------
  
  print_best = fltarr(7,size)
  print_best[0,*] = time 
  print_best[1,*] = phase
  print_best[2,*] = flux_transit_or
  print_best[3,*] = flux_transit
  print_best[4,*] = error
  print_best[5,*] = best_model
  print_best[6,*] = res
  printf,5,print_best
  close,5
  file_move,['Best-Fit_LM.dat'],foldername
  
  ;--------------------------------------
  ;-----Plot Best Fit Model--------------
  ;-------------------------------------
  set_plot,'X'
  color_array=['navy blue','red','orange','green','blue','purple', 'white', 'yellow','black']
  ploterror, phase, flux_transit,error,yrange=[MIN(best_model)-14*STDEV(res),1.01],psym=SYM(1),symsize=0.5,color=Fsc_Color(color_array[8]),xtitle='Time since mid transit (Days)',ytitle='Relative Flux'
  oplot, phase, best_model,color=Fsc_Color(color_array[1]),thick = 3
  oplot,phase,res*0+MIN(best_model)-10*STDEV(res),color=Fsc_Color(color_array[1]),thick=3
  oplot,phase,res+(MIN(best_model)-10*STDEV(res)),color=Fsc_Color(color_array[8]),psym=SYM(1),symsize=.5
  
  
  set_plot,'ps'
  device,filename='plot_LM.ps',/color
  color_array=['navy blue','red','orange','green','blue','purple', 'white', 'yellow','black']
  ploterror, phase, flux_transit,error,yrange=[MIN(best_model)-14*STDEV(res),1.01],psym=SYM(1),symsize=0.5,color=Fsc_Color(color_array[8]),xtitle='Time since mid transit (Days)',ytitle='Relative Flux'
  oplot, phase, best_model,color=Fsc_Color(color_array[1]),thick = 8
  oplot,phase,res*0+MIN(best_model)-10*STDEV(res),color=Fsc_Color(color_array[1]),thick=8
  oplot,phase,res+(MIN(best_model)-10*STDEV(res)),color=Fsc_Color(color_array[8]),psym=SYM(1),symsize=.5
  device,/close
  file_move,['plot_LM.ps'],foldername
  
  
  ;------------------------------------------------------------------------------------------------------------------
  ;-------------Perform the Bootstrap Monte Carlo (find best value and error)----------------------------------------
  ;------------------------------------------------------------------------------------------------------------------
  print, '*******************************************************************************' 
  print, '*****************************MONTE CARLO STARTING******************************'
  print, '*******************************************************************************' 
  
  
  for i=0L, N-1 do begin			; begin Monte Carlo iterations
  	If i eq 5000 then begin 
  		print, 'Are We There Yet (You are at 5000)' 
  	endif
  	If i eq 9000 then begin 
  	   print, 'Just Keep Swimming (You are at 9000)'
  	endif
  	If i eq 50000 then begin 
  	   print, 'What were you thinking? (You are at 50000)'
  	endif
  
  	random = RANDOMN(seed, size)	;find random number with a mean of 0 and stdev of 1
  					; important (seed needs to change every iteration), set by time on computer. 
        for j=0, size-1 do begin
  		flux_transit2(j)=flux_transit(j)+random(j)*error(j)	
        endfor
  
  	FitValue = mpfitfun('model_fit',time,flux_transit2,error,fit_start,PERROR=perror,PARINFO=parinfo,yfit=yfit,/quiet,/AUTODERIVATIVE,BESTNORM= chi_squared,DOF=DOF)
  	
  ;------------------------------
  ;----Find Chi-Squared----------
  ;------------------------------
  	
  ;reduced Chi-Squared
  chi[i] = chi_squared/DOF
  
  	;find differences between assumed best fit value
  	rp_monte2[i] 		    = Best_fit[0]-FitValue[0]
  	a_rstar_monte2[i] 	= Best_fit[1]-FitValue[1]
  	incl_monte2[i]    	= Best_fit[2]-FitValue[2]
  	mid_monte2[i]   	= Best_fit[3]-FitValue[3]
  	lin_monte2[i] 		= Best_fit[4]-FitValue[4] 
  	quad_monte2[i] 		= Best_fit[5]-FitValue[5]
  
  	rp_monte[i] 		= FitValue[0]
  	a_rstar_monte[i] 	= FitValue[1]
  	incl_monte[i]    	= FitValue[2]
  	mid_monte[i]   		= FitValue[3]
  	lin_monte[i] 		= FitValue[4]
  	quad_monte[i] 		= FitValue[5]	
  	printf,12,rp_monte2[i],a_rstar_monte2[i],incl_monte2[i], mid_monte2[i],lin_monte[i],quad_monte[i]
  	printF,1,rp_monte[i],a_rstar_monte[i],incl_monte[i], mid_monte[i],lin_monte[i],quad_monte[i],chi[i]
  
  Endfor 											;End MC simulations 
  	close,1
  	close,12
  	file_move,['monte_value.dat','monte_value2.dat'],foldername+'/LM
  
  ;----------------------------------------------------------------------------------------
  ;------Fit Gaussian to Distribution And find SIGMA (Error) 
  ;----------------------------------------------------------------------------------------
  set_plot,'ps'
  device,filename='gauss.ps',/color
  openw,2,'Final_Parameters.dat'
  
  
  ;------------------------------
  ;------Fit Gaussians-----------
  ;------------------------------
  
  If parinfo_rp eq 0 then begin 
  	histogauss, rp_monte2, A 		; Fit a gaussian function to the distribution 
  	printf,2,Best_fit[0], A[2],'Rp/R*',FORMAT='(F9.7,2X,F9.7,2x,A5)' 	
  endif
  
  If parinfo_ar eq 0 then begin
  	histogauss, a_rstar_monte2, B 
  	printf,2,Best_fit[1], B[2],'a/R*',FORMAT='(F9.7,2X,F9.7,2x,A4)'	
  endif
  
  If parinfo_inc eq 0 then begin
  	histogauss, incl_monte2, C 
  	;convert inclination to degrees 
  	printf,2,Best_fit[2]*(180/!PI), C[2]*(180/!PI), 'Inclination (degrees)',FORMAT='(F9.6,2X,F9.7,2x,A21)'	
  endif
  
  If parinfo_mid eq 0 then begin
  	histogauss, mid_monte2, D
  	printf,2,Best_fit[3], D[2],'Mid-Transit (days)',FORMAT='(F10.7,2X,F9.7,2x,A18)'	
  endif 
  
  If parinfo_lin eq 0 then begin
  	histogauss, lin_monte, E
  	printf,2,Best_fit[4], E[2],'Linear Limb',FORMAT='(F9.7,2X,F9.7,2x,A15)'	
  endif 
  
  If parinfo_quad eq 0 then begin
  	histogauss, quad_monte, F
  	printf,2,Best_fit[5], F[2],'Quadratic Limb',FORMAT='(F9.6,2X,F9.7,2x,A15)'	
  endif 
  	
  
  ;------------------------------------------------------------------------------------
  ;------Output Useful Values (OoT Baseline, Cadence, Duration)------------------------
  ;------------------------------------------------------------------------------------
  
  best_model_input = best_model
  Value = Useful_Values(time,res,best_model_input,flux_transit,size)
  cadence_sec = 	Value[3]*24.d*60.d*60.d
  res_stdev   =	STDEV(res)*1000.d  
  
  If run_baseline eq 0 then begin 
      printf,2,Value[1],Value[2],' Duration (Days)',FORMAT='(F9.7,2X,F9.7,1x,A16)'  
      printf,2,Value[1]*24.d*60.d,Value[2]*24.d*60.d,' Duration (Min)',FORMAT='(F9.5,2X,F9.7,1x,A15)'
  endif 
  printf,2,Value[3],' Cadence (days)',FORMAT='(F9.7,12X,A15)'
  printf,2, cadence_sec, ' Cadence (secs)' ,FORMAT='(F9.5,3X,A24)'
  printf,2,Value[0],'OoT STDEV',FORMAT='(F9.7,13x,A9)'
  printf,2, res_stdev,' RMS of Residuals (mmag)',FORMAT='(F9.7,12X,A24)'
  printf,2,chi_reduced,'Reduced Chi-Squared',FORMAT='(F9.5,12X,A19)'
  printf,2,MIN(chi),' Lowest MC Reduced Chi^2',FORMAT='(F9.5,12X,A24)'
  If info_baseline eq 0 then begin
  	printf,2, Detrend[0],' Baseline Intercept', FORMAT='(F9.5,13X,A19)'
  	printf,2, Detrend[1],' Baseline Slope' , FORMAT='(F9.5,13x,A15)'
  		If run_baseline eq 2 then begin
  			printf,2, Detrend[2],' Baseline x^2 term' , FORMAT='(F9.5,13x,A15)'
  		endif
  endif 
  
  If run_symmetry eq 0 then begin
      ;------------------------------------------------------------------------------------------------------------------
      ;----------------------------Calculation of the Symmetry of the Transit--------------------------------------------
      ;----------This part of the program uses transitsymetry.pro and makesymetric.pro written by Tim Carleton (03/26/11)
      ;------------------------------------------------------------------------------------------------------------------
      
      openw,56,'symmetry_test'
      print_transit = fltarr(2,size) 
      print_transit[0,*] = time 
      print_transit[1,*] = flux_transit
      printf,56,print_transit
      close,56
      
      mid = Best_fit[3]
      subtransit = transitsymetry('symmetry_test',mid)
      
      printf,2,stdev(subtransit[1,*])*1000.d,'STDEV of Asymmetry (mmag)',FORMAT='(F9.7,13X,A25)'
      
      file_move,['symmetry_test','Asymmetry_Plot.ps'],foldername+'/LM
  
  endif
  
  close,/all
  ;device, /close
  file_move,['gauss.ps','Final_Parameters.dat'],foldername+'/LM

  If run_prayer eq 1 then begin
    ;-------------------------------------------Method 1---(There is Gaussian Distribution assumption in this method)---
    ;--------------------------------------(Red Nose Estimation)--------------------------------------------------------
    ;----------------------------Prayer Bead Method (Residual-Permutation Method)---------------------------------------
    ;-------Used to Find out how much correlated noise (Red noise) is a problem-----------------------------------------
    ;-----The general idea is to perform the Bootstrap MC and add the residuals (with errors to the transit) -----------
    ;-----See Southworth 2008 and Fernandez et al. 2009 for more information on PB method-------------------------------
    ;------------------------------------------------------------------------------------------------------------------
    ;----This routine does exactly what is described in these references with the exception that we --------------------
    ;----include the error from the residuals and add residual+evenly disturbed error to the transit, ------------------
    ;----thus allowing for more than d shifts (d is number of data points in the transit)-------------------------------
    ;-----We run this routine N*(size/size-1) times (N is # of Monte Carlo Sims, size # of data points) ----------------
    ;-------------------------------------------------------------------------------------------------------------------
    ;----- There is a Gaussian Distribution assumption in this method, so the true nature of the noise might not be----- 
    ;----- found. Method 2 of the red noise estimation allows for non-gaussian distributions ---------------------------
    ;-------------------------------------------------------------------------------------------------------------------
    
    print,'******************************************************************************' 
    print, '*****************************PRAYER BEAD STARTED*****************************'
    print,'******************************************************************************' 
    
    ;------------------------------------------------
    ;-----Calculate Residuals for Prayer Bead method
    ;------------------------------------------------
    
    prayer_number = N/size
    prayer_number = round(prayer_number) 
    
    res_new			= fltarr(size*(prayer_number+1))
    
    for j=0,prayer_number do begin
    for i = 0,size-1 do begin
    	res_new(i+size*j)     = res(i)
    endfor
    endfor
    
    runs = size*(prayer_number-1) 
    
    ;------------------------------------------------------------------------------------------------------------------
    ;-------------Find Error Through Monte Carlo-----------------------------------------------------------------------
    ;------------------------------------------------------------------------------------------------------------------
    openw,10,'monte_value_prayer.dat'
    openw,20,'Final_Parameters_prayer.dat'
    
    
    for i=0L, runs-1 do begin				; begin Monte Carlo iterations
    	If i eq 5000 then begin 
    		print, 'Where are we? (You are at 5000)' 
    	endif
    	If i eq 9000 then begin 
    		print, 'You still there? (You are at 9000)' 
    	endif
    
    	random = RANDOMN(seed, size)	;find random number with a mean of 0 and stdev of 1
    					; important (seed needs to change every iteration), set by time on computer. 
          for j=0, size-1 do begin
    		flux_transit2(j)=flux_transit(j)+ res_new(j+i)*random(j)
    
          endfor
    
    	FitValue = mpfitfun('model_fit',time,flux_transit2,error,fit_start,PERROR=perror,PARINFO=parinfo,yfit=yfit,/quiet,/AUTODERIVATIVE)
    	
    	;find differences between assumed best fit value
    	rp_monte[i] 		  = Best_fit[0]-FitValue[0]
    	a_rstar_monte[i] 	= Best_fit[1]-FitValue[1]
    	incl_monte[i]    	= Best_fit[2]-FitValue[2]
    	mid_monte[i]   		= Best_fit[3]-FitValue[3]
    	lin_monte[i] 		  = Best_fit[4]-FitValue[4] 
    	quad_monte[i] 		= Best_fit[5]-FitValue[5]
    	printf,10,rp_monte[i],a_rstar_monte[i],incl_monte[i], mid_monte[i],lin_monte[i],quad_monte[i]
    Endfor 
    
    
    ;---------------------------------------------------------
    ;------Fit Gaussian to Distribution And find SIGMA (Error) 
    ;---------------------------------------------------------
    set_plot,'ps'
    device,filename='gauss_prayer.ps',/color
    
    
    If parinfo_rp eq 0.0d then begin 
    	histogauss, rp_monte, AA 							; Fit a gaussian function to the distribution 
    	printf,20,Best_fit[0], AA[2],'Rp/R*',FORMAT='(F9.7,2X,F9.7,2x,A5)' 
    	printf,20,AA[2]/A[2],'Beta_RP',FORMAT='(F9.7,12x,A8)	
    endif
    
    If parinfo_ar eq 0.0d then begin
    	histogauss, a_rstar_monte, BB 
    	printf,20,Best_fit[1], BB[2],'a/R*',FORMAT='(F9.7,2X,F9.7,2x,A4)'
    	printf,20,BB[2]/B[2],'Beta_AR',FORMAT='(F9.7,2x,A8)	
    endif
    
    If parinfo_inc eq 0.0d then begin
    	histogauss, incl_monte, CC 
    
    	;convert inclination to degrees 
    	printf,20,Best_fit[2]*(180/!PI), CC[2]*(180/!PI), 'Inclination (degrees)',FORMAT='(F9.6,2X,F9.7,2x,A21)'
    	printf,20,CC[2]/C[2],'Beta_IN',FORMAT='(F9.6,2x,A8)	
    endif
    
    If parinfo_mid eq 0.0 then begin
    	histogauss, mid_monte, DD
    	printf,20,Best_fit[3], DD[2],'Mid-Transit (days)',FORMAT='(F10.7,2X,F9.7,2x,A18)'	
    	printf,20,DD[2]/D[2],'Beta_MID',FORMAT='(F9.7,12x,A8)'	
    endif 
    
    If parinfo_lin eq 0.0 then begin
    	histogauss, lin_monte, EE
    	printf,20,Best_fit[4], EE[2],'Linear Limb',FORMAT='(F9.7,2X,F9.7,2x,A15)'
    	printf,20,EE[2]/E[2],'Beta_lin',FORMAT='(F9.7,2x,A9)		
    endif 
    
    If parinfo_quad eq 0.0 then begin
    	histogauss, quad_monte, FF
    	printf,20,Best_fit[5], FF[2],'Quadratic Limb',FORMAT='(F9.7,2X,F9.7,2x,A15)'	
    	printf,20,FF[2]/F[2],'Beta_quad',FORMAT='(F9.6,2x,A10)		
    endif 
    
    printf,20,runs,'Number of Simulations',FORMAT='(I5,17X,A21)'
    
    close,10
    device,/close
    
    
    ;-------------------------------------------Method 2----------------------------------------------------
    ;--------------------------------------(Red Nose Estimation)--------------------------------------------
    ;----------------------------Prayer Bead Method (Residual-Permutation Method)---------------------------
    ;-------Used to Find out how much correlated noise (Red noise) is a problem-----------------------------
    ;-------Uses the techniques by Todorov et al 2012 (APJ,746,111) to find the error of a 
    ;-------non-Gaussian distribution.  --------------------------------------------------------------------
    print,'******************************************************************************' 
    print, '*****************************PRAYER BEAD METHOD 2 STARTED********************'
    print,'******************************************************************************' 
    
    openw,12,'monte_value_prayer_method2.dat'
    
    ;------------------------------------------------
    ;-----Calculate Residuals for Prayer Bead method
    ;------------------------------------------------
    
    res_new_prayer2			= res_new
    
    runs = size		
    
    ;just run through once with the residuals 
    for i=0L, runs-1 do begin				; begin Monte Carlo iterations
    	If i eq 1000 then begin 
    		print, 'Lots of Data Points :)' 
    	endif
    	
    	Fit_Value_3 = FitValue
    	
          for j=0, size-1 do begin
            flux_transit3(j)=flux_transit(j)+ res_new_prayer2(j+i)		; only add the residuals to the transit cyclically
          endfor
    
    	FitValue_3 = mpfitfun('model_fit',time,flux_transit3,error,fit_start,PERROR=perror,PARINFO=parinfo,yfit=yfit,/quiet,/AUTODERIVATIVE)
    	
    	;find differences between assumed best fit value
    	rp_monte_prayer2[i] 		  = Best_fit[0]-FitValue_3[0]
    	a_rstar_monte_prayer2[i] 	= Best_fit[1]-FitValue_3[1]
    	incl_monte_prayer2[i]    	= Best_fit[2]-FitValue_3[2]  ;radians
    	mid_monte_prayer2[i]   		= Best_fit[3]-FitValue_3[3]
    	lin_monte_prayer2[i] 		  = Best_fit[4]-FitValue_3[4] 
    	quad_monte_prayer2[i] 		= Best_fit[5]-FitValue_3[5]
    	printf,12,rp_monte_prayer2[i],a_rstar_monte_prayer2[i],incl_monte_prayer2[i], mid_monte_prayer2[i],lin_monte_prayer2[i],quad_monte_prayer2[i]
    Endfor 
    
    ;------------------------------------------------------------------
    ;-------Run EXOMOP_PRAYER (Method 2) for each fitted parameter-----
    ;------------------------------------------------------------------
    If parinfo_rp eq 0.0 then begin 
    	name = 'rp_prayer_histogram.ps'
    	exomop_prayer,rp_monte_prayer2,name,upper_error=upper_error_rp,lower_error=lower_error_rp
    	file_move,['rp_prayer_histogram.ps'],foldername+'/LM
    	printf,20, upper_error_rp/A[2],'Beta_RP_Upper (Method 2)',FORMAT='(F9.7,12x,A25)'	        ;LM
    	printf,20, abs(lower_error_rp)/A[2],'Beta_RP_Lower (Method 2)',FORMAT='(F9.7,12x,A25)'	   ;LM
    endif
    
    If parinfo_ar eq 0.0 then begin
    	name = 'ar_prayer_histogram.ps'
    	exomop_prayer, a_rstar_monte_prayer2,name,upper_error=upper_error_ar,lower_error=lower_error_ar
    	file_move,['ar_prayer_histogram.ps'],foldername+'/LM
    	printf,20, upper_error_ar/B[2],'Beta_AR_Upper (Method 2)',FORMAT='(F9.7,2x,A25)'	        ;LM
    	printf,20,abs(lower_error_ar)/B[2],'Beta_AR_Lower (Method 2)',FORMAT='(F9.7,2x,A25)'	    ;LM
    endif
    
    If parinfo_inc eq 0.0 then begin
    	name = 'inc_prayer_histogram.ps'
    	exomop_prayer, incl_monte_prayer2,name,upper_error=upper_error_inc,lower_error=lower_error_inc
    		file_move,['inc_prayer_histogram.ps'],foldername+'/LM
    	printf,20, upper_error_inc/C[2],'Beta_IN_Upper (Method 2)',FORMAT='(F9.6,2x,A25)'	      ;LM  
    	printf,20, abs(lower_error_inc)/C[2],'Beta_IN_Lower (Method 2)',FORMAT='(F9.6,2x,A25)'	 ;LM
    
    endif
    
    If parinfo_mid eq 0.0 then begin
    	name = 'mid_prayer_histogram.ps'
    	exomop_prayer,mid_monte_prayer2,name,upper_error=upper_error_mid,lower_error=lower_error_mid
    	file_move,['mid_prayer_histogram.ps'],foldername+'/LM
    	printf,20, upper_error_mid/D[2],'Beta_MID_Upper (Method 2)',FORMAT='(F9.7,12x,A25)'	    ;LM
    	printf,20,abs(lower_error_mid)/D[2],'Beta_MID_Lower (Method 2)',FORMAT='(F9.7,12x,A25)'	;LM
    
    endif 
    
    If parinfo_lin eq 0.0 then begin
    	name = 'lin_prayer_histogram.ps'
    	exomop_prayer,lin_monte_prayer2,name,upper_error=upper_error_lin,lower_error=lower_error_lin
    	file_move,['lin_prayer_histogram.ps'],foldername+'/LM
    	printf,20, upper_error_lin/E[2],'Beta_lin_Upper (Method 2)',FORMAT='(F9.7,2x,A24)'    ;LM
    	printf,20, abs(lower_error_lin)/E[2],'Beta_lin_Lower (Method 2)',FORMAT='(F9.7,2x,A24)'	;LM	
    
    endif 
    
    If parinfo_quad eq 0.0 then begin
    	name = 'quad_prayer_histogram.ps'
    	exomop_prayer,quad_monte_prayer2,name,upper_error=upper_error_quad,lower_error=lower_error_quad
    	file_move,['quad_prayer_histogram.ps'],foldername+'/LM
    	printf,20, upper_error_quad/F[2],'Beta_quad_Upper (Method 2)',FORMAT='(F9.6,2x,A26)'		      ;LM
    	printf,20, abs(lower_error_quad)/F[2],'Beta_quad_Lower (Method 2)',FORMAT='(F9.6,2x,A26)'   ;LM		
    
    endif 
    
    
    close,20
    close,12
    file_move,['Final_Parameters_prayer.dat','gauss_prayer.ps','monte_value_prayer.dat','monte_value_prayer_method2.dat'],foldername+'/LM
    


endif 	;close Prayer Bead Calculation

If run_beta eq 1 then begin
    ;-------------------------------------------------------------------------------------------------------
    ;--------------------------------------(Red Nose Estimation)--------------------------------------------
    ;----------------------------Time-Averaging Technique (Beta Method)--------------------------------------
    ;-------Used to Find out how much correlated noise (Red noise) is a problem-----------------------------
    ;-----See (Pont et al 2006, 373,231) for more information on this method--------------------------------
    ;------------------------------------------------------------------------------------------------------
    ;--------------------------------------------------------------------------------------------------------
    
    print,'******************************************************************************' 
    print, '*********************RED NOISE BETA METHOD STARTED***************************'
    print,'******************************************************************************' 
    
    
    beta_stuff = beta_exomop(time,res,error) 
    
    openw,99,'Final_Parameters_Beta.dat'
    printf,99,Beta_stuff[0], Beta_stuff[2],'White Noise',FORMAT='(F9.7,2X,F9.7,2x,A11)'
    printf,99,Beta_stuff[1], Beta_stuff[3],'Red Noise',FORMAT='(F9.7,2X,F9.7,2x,A9)'
    printf,99,Beta_stuff[4],'Beta',FORMAT='(F9.7,13x,A4)'
    close,99
    
     file_move,['beta_plot.ps','Beta_value.dat','Final_Parameters_Beta.dat'], foldername+'/LM
    
    endif 	;end run_beta
    
    save,/ALL,FILENAME ='everything_LM.sav'	;save everything
    file_move,['everything_LM.sav'],foldername+'/LM' ;move
    
    ;endif   ; END LM FIT
    ;close,/all
    
    If run_wavelet eq 1 then begin    ;Wavelet Red Noise LM
    ;-------------------------------------------------------------------------------------------------------
    ;--------------------------------------(Red Nose Estimation)--------------------------------------------
    ;--------------------------------------WAVELET-BASED METHOD --------------------------------------------
    ;---Used to find the amount of red and white noise, assumed 1/f^alpha red noise------------------------
    ;---See (Carter & Winn 2009, APJ, 704, 51) for more information on the method
    ;------------------------------------------------------------------------------------------------------
    ;--------------------------------------------------------------------------------------------------------
    
    print,'******************************************************************************' 
    print, '*********************WAVELET METHOD STARTED**********************************'
    print,'******************************************************************************' 
    
    ;Run Wavelet Red Noise and output Sigma Red Noise and Sigma White Noise 
      
     exomop_wavelet,res,rms_red=sigma_wavelet_red_LM, $
        rms_white=sigma_wavelet_white_LM
    
       beta_wavelet_LM = sqrt(1.d + (sigma_wavelet_red_LM/sigma_wavelet_white_LM)^2.d)
       
       openw,99,'Wavelet_red_noise_LM.dat'
       printf,99,sigma_wavelet_red_LM,'Red Noise'
       printf,99, sigma_wavelet_white_LM,'White Noise'
       printf,99,beta_wavelet_LM,'Beta (Wavelet)'
       close,99
      
       file_move,['Wavelet_red_noise_LM.dat'], foldername+'/LM
      
      
    endif ; END wavelet 
    
  endif   ; END LM FIT
    close,/all
    
    
    If run_MCMC eq 0 then begin
    ;--------------------------------------------------------------------------------------------------------
    ;-------------Differential Evolution Markov Chain Monte Carlo Fit (MCMC)-------------------------------
    ;--------------------------------------------------------------------------------------------------------
    ;---------Uses a Differential Evolution Markov Chain Monte Carlo fit (ter Braak, 2006)-------------------
    ;-------http://www.stat.columbia.edu/~gelman/stuff_for_blog/cajo.pdf.------------------------------------ 
    ;--------------------------------------------------------------------------------------------------------
    ;-----The program stops when the chains are well-mixed, as defined by------------------------------------
    ;-----Ford 2006 (http://adsabs.harvard.edu/abs/2006ApJ...642..505F)--------------------------------------
    ;-----using the Gelman-Rubin Statistic and number of independent draws,----------------------------------
    ;-----or when each chain has taken MAXSTEPS, whichever is first.-----------------------------------------
    ;--------------------------------------------------------------------------------------------------------
    ;--Uses exofast_demc and exofast_plotdist created by Eastman 2013 (http://arxiv.org/pdf/1206.5798v3.pdf)-
    ;--------------------------------------------------------------------------------------------------------
    
    
    print,'******************************************************************************' 
    print, '********************* MCMC Method Started ***********************************'
    print,'******************************************************************************' 
    
    If maxsteps lt 1000000 then begin
    	maxsteps = 1000000
    endif 
    If nchains lt 9 then begin 
    	nchains = 20
    endif 
    
    exomop_mcmcv2,bestpars,parinfo.fixed,maxsteps,nchains,tofit=tofit,res_mcmc=res_mcmc,$
    	    best_model_mcmc=best_model_mcmc,bestpars_mcmc=bestpars_mcmc,mcmcerrors=mcmcerrors
    
    ;-----------------------------------------------
    ;-----Print Transit and Best fit Model----------
    ;-----------------------------------------------
    openw,5,'Best-Fit_mcmc.dat' , width=250
    print_best_mcmc = fltarr(7,size)
    print_best_mcmc[0,*] = time 
    print_best_mcmc[1,*] = phase
    print_best_mcmc[2,*] = flux_transit_or
    print_best_mcmc[3,*] = flux_transit
    print_best_mcmc[4,*] = error
    print_best_mcmc[5,*] = best_model_mcmc
    print_best_mcmc[6,*] = res_mcmc
    printf,5,print_best
    close,5
    file_move,['Best-Fit_mcmc.dat'],foldername
    
    file_move,['probability_plot.ps','Final_Parameters_MCMC.dat','covar.ps'], foldername+'/MCMC'
    file_move,['plot_mcmc.ps'], foldername
    
    save,/ALL,FILENAME ='everything_MCMC.sav'		;save everything
    file_move,['everything_MCMC.sav'],foldername+'/MCMC' 	;move
    
    
    If run_prayer eq 1 then begin		;for MCMC fit
    ;-------------------------------------------Method 1---(There is a Gaussian Distribution assumption in this method)---
    ;--------------------------------------(Red Nose Estimation for MCMC FIT)-------------------------------------------
    ;----------------------------Prayer Bead Method (Residual-Permutation Method)---------------------------------------
    ;-------Used to Find out how much correlated noise (Red noise) is a problem-----------------------------------------
    ;-----The general idea is to perform the Bootstrap MC and add the residuals (with errors to the transit) ----------
    ;-----See Southworth 2008 and Fernandez et al. 2009 for more information on PB method-----------------------------
    ;-----------------------------------------------------------------------------------------------------------------
    ;----This routine does exactly what is described in these references with the exception that we ------------------
    ;----include the error from the residuals and add residual+evenly disturbed error to the transit, -----------------
    ;----thus allowing for more than d shifts (d is number of data points in the transit)------------------------------
    ;-----We run this routine N*(size/size-1) times (N is # of Monte Carlo Sims, size # of data points) ---------------
    ;-------------------------------------------------------------------------------------------------------------------
    ;----- There is a Gaussian Distribution assumption in this method, so the true nature of the noise might not be------- 
    ;----- found. Method 2 of the prayer bead red noise estimation allows for non-gaussian distributions ------------------
    ;-------------------------------------------------------------------------------------------------------------------
    
    print,'******************************************************************************' 
    print, '*****************************PRAYER BEAD STARTED (MCMC)**********************'
    print,'******************************************************************************' 
    
    ;-----------------------------------------------
    ;-----Update errors-------------------------
    ;-----------------------------------------------
    A_mcmc   = mcmcerrors[0,0] 	;rp/R*
    A_mcmc_2 = mcmcerrors[1,0]
    
    B_mcmc   = mcmcerrors[0,1] 	;a/R*
    B_mcmc_2 = mcmcerrors[1,1]
    
    C_mcmc   = mcmcerrors[0,2] 	;inclination in degrees
    C_mcmc_2 = mcmcerrors[1,2] 
    
    D_mcmc   = mcmcerrors[0,3] 	;mid-transit
    D_mcmc_2 = mcmcerrors[1,3]
    
    E_mcmc   = mcmcerrors[0,4] 	;lin limb 
    E_mcmc_2 = mcmcerrors[1,4]
    
    F_mcmc   = mcmcerrors[0,5] 	;quad limb 
    F_mcmc_2 = mcmcerrors[1,5]
    
    ;------------------------------------------------
    ;-----Calculate Residuals for Prayer Bead method
    ;------------------------------------------------
    
    res_new_mcmc			= fltarr(size*(prayer_number+1))
    
    for j=0,prayer_number do begin
    for i = 0,size-1 do begin
    	res_new(i+size*j)     = res_mcmc(i)
    endfor
    endfor
    
    runs = size*(prayer_number-1) 
    
    ;------------------------------------------------------------------------------------------------------------------
    ;-------------Find Error Through Monte Carlo-----------------------------------------------------------------------
    ;------------------------------------------------------------------------------------------------------------------
    openw,10,'monte_value_prayer_mcmc.dat'
    openw,20,'Final_Parameters_prayer_mcmc.dat'
    
    flux_transit2_mcmc = flux_transit
    
    for i=0L, runs-1 do begin				; begin Monte Carlo iterations
    	If i eq 5000 then begin 
    		print, 'On many long journeys have we gone? (You are at 5000)' 
    	endif
    	If i eq 9000 then begin 
    		print, 'Where you say? (You are at 9000)' 
    	endif
    
    	    random = RANDOMN(seed, size)	;find random number with a mean of 0 and stdev of 1
    					                          ; important (seed needs to change every iteration), set by time on computer. 
          for j=0, size-1 do begin
            flux_transit2_mcmc(j)=flux_transit(j)+ res_new(j+i)*random(j)
          endfor
    
    ;------------------------------------------------------------------------------
    ; Use the full MCMC model for the prayer bead execution------------------------
    ;------------------------------------------------------------------------------
    If run_long_prayer eq 0 then begin 		
    ;---------------------------------------------------
    ;-----Updated Structure with the permuted transit---
    ;---------------------------------------------------
     ;save structure of transit
    	trandata_mcmc_Res =  create_struct('time',time,'flux', flux_transit2_mcmc,'err',error)
    	transit=trandata_mcmc_Res
    
    	exomop_mcmcv2,bestpars,parinfo.fixed,maxsteps,nchains,tofit=tofit,res_mcmc=res_mcmc,$
    		    best_model_mcmc=best_model_mcmc,bestpars_mcmc=bestpars_mcmc_res,mcmcerrors=mcmcerrors
    	device,/close
    	
     ;find differences between assumed best fit value (for full mcmc fit)
    	 rp_monte_mcmc[i] 		      	= bestpars_mcmc[0]-bestpars_mcmc_res[0]
    	 a_rstar_monte_mcmc[i] 			= bestpars_mcmc[1]-bestpars_mcmc_res[1]
    	 incl_monte_mcmc[i]    			= bestpars_mcmc[2]-(bestpars_mcmc_res[2])
    	 mid_monte_mcmc[i]   		  	= bestpars_mcmc[3]-bestpars_mcmc_res[3]
    	 lin_monte_mcmc[i] 		   	 = bestpars_mcmc[4]-bestpars_mcmc_res[4] 
    	 quad_monte_mcmc[i] 	    		= bestpars_mcmc[5]-bestpars_mcmc_res[5]
    	 printf,10,rp_monte_mcmc[i],a_rstar_monte_mcmc[i],incl_monte_mcmc[i], mid_monte_mcmc[i],lin_monte_mcmc[i],quad_monte_mcmc[i]
    endif 
    
    
    ;------------------------------------------------------------------------------
    ;Only use MCMC residuals for prayer bead execution along with mpfitfun---------
    ;------------------------------------------------------------------------------
    If run_long_prayer eq 1 then begin 	
    	FitValue_mcmc = mpfitfun('model_fit',time,flux_transit2_mcmc,error,bestpars,PERROR=perror_mcmc,PARINFO=parinfo,yfit=yfit,/quiet,/AUTODERIVATIVE)
    
    
    	;find differences between assumed best fit value
    	rp_monte_mcmc[i] 		      	= bestpars_mcmc[0]-FitValue_mcmc[0]
    	a_rstar_monte_mcmc[i] 			= bestpars_mcmc[1]-FitValue_mcmc[1]
    ;	incl_monte_mcmc[i]    			= bestpars_mcmc[2]-(FitValue_mcmc[2])
    	incl_monte_mcmc[i]    			= bestpars_mcmc[2]-(FitValue_mcmc[2]*(180.d/(!PI)))
    	mid_monte_mcmc[i]   		  	= bestpars_mcmc[3]-FitValue_mcmc[3]
    	lin_monte_mcmc[i] 		   	= bestpars_mcmc[4]-FitValue_mcmc[4] 
    	quad_monte_mcmc[i] 	    		= bestpars_mcmc[5]-FitValue_mcmc[5]
    	printf,10,rp_monte_mcmc[i],a_rstar_monte_mcmc[i],incl_monte_mcmc[i], mid_monte_mcmc[i],lin_monte_mcmc[i],quad_monte_mcmc[i]
    endif 
    
    Endfor 	;end Monte Carlo Iterations
    
    ;---------------------------------------------------------
    ;------Fit Gaussian to Distribution And find SIGMA (Error) 
    ;---------------------------------------------------------
    set_plot,'ps'
    device,filename='gauss_prayer_mcmc.ps',/color
    
    
    If parinfo_rp eq 0.0 then begin 
    	histogauss, rp_monte_mcmc, AA_mcmc 							; Fit a gaussian function to the distribution 
    	printf,20, bestpars_mcmc[0], AA_mcmc[2],'Rp/R*',FORMAT='(F9.7,2X,F9.7,2x,A5)' 
    	Beta_rp_mcmc = AA_mcmc[2]/A_mcmc
    	printf,20,AA_mcmc[2]/A_mcmc,'Beta_RP',FORMAT='(F9.7,12x,A8)
    	printf,20,AA_mcmc[2]/A_mcmc_2,'Beta_RP',FORMAT='(F9.7,12x,A8)
    endif
    
    If parinfo_ar eq 0 then begin
    	histogauss, a_rstar_monte_mcmc, BB_mcmc 
    	Beta_ar_mcmc = BB_mcmc[2]/B_mcmc
    	printf,20, bestpars_mcmc[1], BB_mcmc[2],'a/R*',FORMAT='(F9.7,2X,F9.7,2x,A4)'
    	printf,20,BB_mcmc[2]/B_mcmc,'Beta_AR',FORMAT='(F9.7,2x,A8)	
    endif
    
    If parinfo_inc eq 0 then begin
    	histogauss, incl_monte_mcmc, CC_mcmc 
    	printf,20, bestpars_mcmc[2], CC_mcmc[2], 'Inclination (degrees)',FORMAT='(F9.6,2X,F9.7,2x,A21)'
    	printf,20,CC_mcmc[2]/C_mcmc,'Beta_IN',FORMAT='(F9.7,2x,A8)	
    endif
    
    If parinfo_mid eq 0 then begin
    	histogauss, mid_monte_mcmc, DD_mcmc
    	Beta_mid_mcmc = DD_mcmc[2]/DD_mcmc
    	printf,20, bestpars_mcmc[3], DD_mcmc[2],'Mid-Transit (days)',FORMAT='(F10.7,2X,F9.7,2x,A18)'	
    	printf,20,DD_mcmc[2]/D_mcmc,'Beta_MID',FORMAT='(F9.7,12x,A8)'
    	printf,20,DD_mcmc[2]/D_mcmc_2,'Beta_MID',FORMAT='(F9.7,12x,A8)'	
    endif 
    
    If parinfo_lin eq 0 then begin
    	histogauss, lin_monte_mcmc, EE_mcmc
    	printf,20, bestpars_mcmc[4], EE[2],'Linear Limb',FORMAT='(F9.7,2X,F9.7,2x,A15)'
    	printf,20,EE_mcmc[2]/E_mcmc,'Beta_lin',FORMAT='(F9.7,2x,A9)		
    endif 
    
    If parinfo_quad eq 0 then begin
    	histogauss, quad_monte_mcmc, FF_mcmc
    	printf,20, bestpars_mcmc[5], FF_mcmc[2],'Quadratic Limb',FORMAT='(F9.7,2X,F9.7,2x,A15)'	
    	printf,20,FF_mcmc[2]/F_mcmc,'Beta_quad',FORMAT='(F9.7,2x,A10)		
    endif 
    
    printf,20,runs,'Number of Simulations',FORMAT='(I5,17X,A21)'
    
    close,10
    device,/close
    
    ;-------------------------------------------Method 2------(Assumes a non-Gaussian error distribution)--
    ;--------------------------------------(Red Nose Estimation)--------------------------------------------
    ;----------------------------Prayer Bead Method (Residual-Permutation Method)---------------------------
    ;-------Used to Find out how much correlated noise (Red noise) is a problem-----------------------------
    ;-------Uses the techniques by Todorov et al 2012 (APJ,746,111) to find the error of a 
    ;-------non-Gaussian distribution.  --------------------------------------------------------------------
    print,'******************************************************************************' 
    print, '*****************************PRAYER BEAD METHOD 2 STARTED********************'
    print,'******************************************************************************' 
    
    openw,12,'monte_value_prayer_method2_mcmc.dat'
    
    ;------------------------------------------------
    ;-----Calculate Residuals for Prayer Bead method
    ;------------------------------------------------
    
    res_new_prayer2_mcmc     = res_new	; set new residual variable
    
    runs = size   				;runs is only the size of the array
    Fit_Value_4 = FitValue_mcmc		;set parameters
    
    ;just run through once with the residuals 
    for i=0L, runs-1 do begin       ; begin Monte Carlo iterations
      If i eq 1000 then begin 
        print, 'Lots of Data Points :)' 
      endif
      
          for j=0, size-1 do begin
            flux_transit4(j)=flux_transit(j)+ res_new_prayer2_mcmc(j+i)    ; only add the residuals to the transit cyclically
          endfor
    
    
    ;------------------------------------------------------------------------------
    ; Use the full MCMC model for the prayer bead execution------------------------
    ;------------------------------------------------------------------------------
    If run_long_prayer eq 0 then begin 		  ;will take a long time
    ;---------------------------------------------------
    ;-----Updated Structure with the permuted transit---
    ;---------------------------------------------------
     ;save structure of transit
    	trandata_mcmc_Res2 =  create_struct('time',time,'flux', flux_transit4,'err',error)
    	transit=trandata_mcmc_Res2
    
    	exomop_mcmcv2,bestpars,parinfo.fixed,maxsteps,nchains,tofit=tofit,res_mcmc=res_mcmc,$
    		    best_model_mcmc=best_model_mcmc,bestpars_mcmc=bestpars_mcmc_res2,mcmcerrors=mcmcerrors
    	device,/close
    	
     ;find differences between assumed best fit value (for full mcmc fit)
    	 rp_monte_prayer2_mcmc[i] 		      	= bestpars_mcmc[0]-bestpars_mcmc_res2[0]
    	 a_rstar_monte_prayer2_mcmc[i] 			= bestpars_mcmc[1]-bestpars_mcmc_res2[1]
    	 incl_monte_prayer2_mcmc[i]    			= bestpars_mcmc[2]-(bestpars_mcmc_res2[2])
    	 mid_monte_prayer2_mcmc[i]   		  	= bestpars_mcmc[3]-bestpars_mcmc_res2[3]
    	 lin_monte_prayer2_mcmc[i] 		   	 = bestpars_mcmc[4]-bestpars_mcmc_res2[4] 
    	 quad_monte_prayer2_mcmc[i] 	    		= bestpars_mcmc[5]-bestpars_mcmc_res2[5]
      printf,12,rp_monte_prayer2_mcmc[i],a_rstar_monte_prayer2_mcmc[i],incl_monte_prayer2_mcmc[i], mid_monte_prayer2_mcmc[i],lin_monte_prayer2_mcmc[i],quad_monte_prayer2_mcmc[i]
    endif 
    
    ;------------------------------------------------------------------------------
    ;Only use MCMC residuals for prayer bead execution along with mpfitfun---------
    ;------------------------------------------------------------------------------
    If run_long_prayer eq 1 then begin 	
      	FitValue_4 = mpfitfun('model_fit',time,flux_transit4,error,bestpars,PERROR=perror,PARINFO=parinfo,yfit=yfit,/quiet,/AUTODERIVATIVE)
    
    	;find differences between assumed best fit value
    	rp_monte_prayer2_mcmc[i] 		      	= bestpars_mcmc[0]-FitValue_4[0]
    	a_rstar_monte_prayer2_mcmc[i] 			= bestpars_mcmc[1]-FitValue_4[1]
    ;	incl_monte_prayer2_mcmc[i]    			= bestpars_mcmc[2]-(FitValue_4[2])
    	incl_monte_prayer2_mcmc[i]    			= bestpars_mcmc[2]-(FitValue_4[2]*(180.d/(!PI)))
    	mid_monte_prayer2_mcmc[i]   		  	= bestpars_mcmc[3]-FitValue_4[3]
    	lin_monte_prayer2_mcmc[i] 		   	= bestpars_mcmc[4]-FitValue_4[4] 
    	quad_monte_prayer2_mcmc[i] 	    		= bestpars_mcmc[5]-FitValue_4[5]
      printf,12,rp_monte_prayer2_mcmc[i],a_rstar_monte_prayer2_mcmc[i],incl_monte_prayer2_mcmc[i], mid_monte_prayer2_mcmc[i],lin_monte_prayer2_mcmc[i],quad_monte_prayer2_mcmc[i]
    endif 
    
    endfor ; end Monte Carlo 
    ;------------------------------------------------------------------
    ;-------Run EXOMOP_PRAYER (Method 2) for each fitted parameter-----
    ;------------------------------------------------------------------
    If parinfo_rp eq 0 then begin 
      name = 'rp_prayer_histogram.ps'
      exomop_prayer,rp_monte_prayer2_mcmc,name,upper_error=upper_error_rp_mcmc,lower_error=lower_error_rp_mcmc
      file_move,['rp_prayer_histogram.ps'],foldername+'/MCMC
      printf,20, upper_error_rp_mcmc/A_mcmc,'Beta_RP_Upper (Method 2)',FORMAT='(F9.7,12x,A25)'      ;MCMC
      printf,20, abs(lower_error_rp_mcmc)/A_mcmc_2,'Beta_RP_Lower (Method 2)',FORMAT='(F9.7,12x,A25)' ;MCMC 
    endif
    
    If parinfo_ar eq 0 then begin
      name = 'ar_prayer_histogram.ps'
      exomop_prayer, a_rstar_monte_prayer2_mcmc,name,upper_error=upper_error_ar_mcmc,lower_error=lower_error_ar_mcmc
      file_move,['ar_prayer_histogram.ps'],foldername+'/MCMC
      printf,20, upper_error_ar_mcmc/B_mcmc,'Beta_AR_Upper (Method 2)',FORMAT='(F9.7,2x,A25)'     ;MCMC
      printf,20,abs(lower_error_ar_mcmc)/B_mcmc_2,'Beta_AR_Lower (Method 2)',FORMAT='(F9.7,2x,A25)' ;MCMC
    endif
    
    If parinfo_inc eq 0 then begin
      name = 'inc_prayer_histogram.ps'
      exomop_prayer, incl_monte_prayer2_mcmc,name,upper_error=upper_error_inc_mcmc,lower_error=lower_error_inc_mcmc
      
      ;upper_error_inc_mcmc = upper_error_inc_mcmc*(180/!PI) ;convert to degrees
      ;lower_error_inc_mcmc = lower_error_inc_mcmc*(180/!PI)
    
      file_move,['inc_prayer_histogram.ps'],foldername+'/MCMC
      printf,20, upper_error_inc_mcmc/C_mcmc,'Beta_IN_Upper (Method 2)',FORMAT='(F9.6,2x,A25)'      ;MCMC
      printf,20, abs(lower_error_inc_mcmc)/C_mcmc_2,'Beta_IN_Lower (Method 2)',FORMAT='(F9.6,2x,A25)' ;MCMC
    endif
    
    If parinfo_mid eq 0 then begin
      name = 'mid_prayer_histogram.ps'
      exomop_prayer,mid_monte_prayer2_mcmc,name,upper_error=upper_error_mid_mcmc,lower_error=lower_error_mid_mcmc
      file_move,['mid_prayer_histogram.ps'],foldername+'/MCMC
      printf,20, upper_error_mid_mcmc/D_mcmc,'Beta_MID_Upper (Method 2)',FORMAT='(F9.7,12x,A25)'     ;MCMC
      printf,20,abs(lower_error_mid_mcmc)/D_mcmc_2,'Beta_MID_Lower (Method 2)',FORMAT='(F9.7,12x,A25)' ;MCMC 
    
    endif 
    
    If parinfo_lin eq 0 then begin
      name = 'lin_prayer_histogram.ps'
      exomop_prayer,lin_monte_prayer2,name,upper_error=upper_error_lin,lower_error=lower_error_lin
      file_move,['lin_prayer_histogram.ps'],foldername+'/MCMC
      printf,20, upper_error_lin/E_mcmc,'Beta_lin_Upper (Method 2)',FORMAT='(F9.7,2x,A24)'       ;MCMC
      printf,20, abs(lower_error_lin)/E_mcmc_2,'Beta_lin_Lower (Method 2)',FORMAT='(F9.7,2x,A24)'  ;MCMC  
    
    endif 
    
    If parinfo_quad eq 0 then begin
      name = 'quad_prayer_histogram.ps'
      exomop_prayer,quad_monte_prayer2_mcmc,name,upper_error=upper_error_quad_mcmc,lower_error=lower_error_quad_mcmc
      file_move,['quad_prayer_histogram.ps'],foldername+'/MCMC
      printf,20, upper_error_quad_mcmc/F_mcmc,'Beta_quad_Upper (Method 2)',FORMAT='(F9.6,2x,A26)'      ;MCMC   
      printf,20, abs(lower_error_quad_mcmc)/F_mcmc_2,'Beta_quad_Lower (Method 2)',FORMAT='(F9.6,2x,A26)' ;MCMC   
    
    endif 
    
    close,20
    close,12
    file_move,['Final_Parameters_prayer_mcmc.dat','gauss_prayer_mcmc.ps','monte_value_prayer_mcmc.dat','monte_value_prayer_method2_mcmc.dat','Gelman_Rubin.dat'],foldername+'/MCMC'


endif 	;end Red noise (Prayer Bead)


If run_beta eq 1 then begin	;MCMC
;-------------------------------------------------------------------------------------------------------
;--------------------------------------(Red Nose Estimation MCMC)---------------------------------------
;----------------------------Time-Averaging Technique (Beta Method)--------------------------------------
;-------Used to Find out how much correlated noise (Red noise) is a problem-----------------------------
;-----See (Pont et al 2006, 373,231) for more information on this method--------------------------------
;------------------------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------------------------

print,'******************************************************************************' 
print, '*********************RED NOISE BETA METHOD STARTED (MCMC) *******************'
print,'******************************************************************************' 


beta_stuff_mcmc = beta_exomop(time,res_mcmc,error)       ;Run Beta method

openw,lun,'Final_Parameters_Beta_mcmc.dat',/get_lun
printf,lun,Beta_stuff_mcmc[0], Beta_stuff_mcmc[2],'White Noise',FORMAT='(F9.7,2X,F9.7,2x,A11)'
printf,lun,Beta_stuff_mcmc[1], Beta_stuff_mcmc[3],'Red Noise',FORMAT='(F9.7,2X,F9.7,2x,A9)'
printf,lun,Beta_stuff_mcmc[4],'Beta',FORMAT='(F9.7,13x,A4)'
free_lun,lun

 file_move,['beta_plot.ps','Beta_value.dat','Final_Parameters_Beta_mcmc.dat'], foldername+'/MCMC

endif 	;end run_beta



If run_wavelet eq 1 then begin    ;Wavelet Red Noise LM
  ;-------------------------------------------------------------------------------------------------------
  ;--------------------------------------(Red Nose Estimation)--------------------------------------------
  ;--------------------------------------WAVELET-BASED METHOD --------------------------------------------
  ;---Used to find the amount of red and white noise, assumed 1/f^alpha red noise------------------------
  ;---See (Carter & Winn 2009, APJ, 704, 51) for more information on the method
  ;------------------------------------------------------------------------------------------------------
  ;--------------------------------------------------------------------------------------------------------
  
  print,'******************************************************************************' 
  print, '*********************WAVELET METHOD STARTED**********************************'
  print,'******************************************************************************' 
  
  ;Run Wavelet Red Noise and output Sigma Red Noise and Sigma White Noise 
    
  ;Run Wavelet Red Noise and output Sigma Red Noise and Sigma White Noise 
    
   exomop_wavelet,res_mcmc,rms_red=sigma_wavelet_red_MCMC, $
      rms_white=sigma_wavelet_white_MCMC
  
     beta_wavelet_MCMC = sqrt(1.d + (sigma_wavelet_red_MCMC/sigma_wavelet_white_MCMC)^2.d)
     
     openw,99,'Wavelet_red_noise_MCMC.dat'
     printf,99,sigma_wavelet_red_MCMC,'Red Noise'
     printf,99, sigma_wavelet_white_MCMC,'White Noise'
     printf,99,beta_wavelet_MCMC,'Beta (Wavelet)'
     close,99
    
     file_move,['Wavelet_red_noise_MCMC.dat'], foldername+'/MCMC
    
    
  endif ; END wavelet 
  
  close,/all


endif   ; END MCMC FIT

;------------------------------------------------------------------------------------------------------------------
;------------FINAL PARAMETERS and ERRORS---------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------------------

openw,20,'FINAL_EVERYTHING.dat'
printf,20, BIC,'Bayesian Information Criterion',FORMAT='(F9.3,3x,A30)'
printf,20,stdev(res)*1000.d,'Residuals from MCLM Method (mmag)',FORMAT='(F9.7,2X,A34)' 

If run_MCMC eq 0 then begin
	printf,20,stdev(res_mcmc)*1000.d,'Residuals from MCMC Method (mmag)',FORMAT='(F9.7,2X,A34)' 
	printf,20,''
	
	printf,20,'****** Differential Evolution Markov Chain Monte Carlo Values and Errors ********'
	If parinfo_rp eq 0 then begin 
		printf,20, bestpars_mcmc[0], mcmcerrors[0,0],mcmcerrors[1,0],'Rp/R*',FORMAT='(F9.7,2X,F9.7,2x,F9.7,2x,A5)' 	
	endif

	If parinfo_ar eq 0 then begin
		printf,20, bestpars_mcmc[1], mcmcerrors[0,1],mcmcerrors[1,1],'a/R*',FORMAT='(F9.7,2X,F9.7,2X,F9.7,2x,A4)'	
	endif

	If parinfo_inc eq 0 then begin
		printf,20, bestpars_mcmc[2]*(180/!PI), mcmcerrors[0,2]*(180/!PI),mcmcerrors[1,2]*(180/!PI), 'Inclination (degrees)',FORMAT='(F9.6,2X,F9.7,2x,F9.7,2x,A21)'	
	endif

	If parinfo_mid eq 0 then begin
		printf,20, bestpars_mcmc[3], mcmcerrors[0,3],mcmcerrors[1,3],'Mid-Transit (days)',FORMAT='(F10.7,2X,F9.7,2x,F9.7,1x,A18)'	
	endif 

	If parinfo_lin eq 0 then begin
		printf,20, bestpars_mcmc[4],mcmcerrors[0,4],mcmcerrors[1,4],'Linear Limb',FORMAT='(F9.7,2X,F9.7,2x,F9.7,A15)'	
	endif 

	If parinfo_quad eq 0 then begin
		printf,20, bestpars_mcmc[5], mcmcerrors[0,5],mcmcerrors[1,5],'Quadratic Limb',FORMAT='(F9.7,2X,F9.7,2X,F9.7,1x,A15)'	
	endif 
	printf,20,Beta_stuff_mcmc[0], Beta_stuff_mcmc[2],'White Noise (MCMC):Time',FORMAT='(F9.7,2X,F9.7,6x,A30)'
	printf,20,Beta_stuff_mcmc[1], Beta_stuff_mcmc[3],'Red Noise (MCMC):Time',FORMAT='(F9.7,2X,F9.7,4x,A30)'
	printf,20,Beta_stuff_mcmc[4],'Beta (MCMC):Time',FORMAT='(F9.7,24x,A15)'
	printf,20,sigma_wavelet_white_MCMC,'White Noise (MCMC):Wavelet',FORMAT='(F9.7,24x,A26)'
  printf,20,sigma_wavelet_red_MCMC,'Red Noise (MCMC):Wavelet',FORMAT='(F9.7,23x,A25)'
  printf,20,beta_wavelet_MCMC,'Beta (MCMC):Wavelet',FORMAT='(F9.7,24x,A19)'

;print Betas
	If parinfo_rp eq 0 then begin 
		printf,20,AA_mcmc[2]/MAX([mcmcerrors[0,0],mcmcerrors[1,0]]),'Beta_RP',FORMAT='(F9.7,22x,A9)'
		printf,20, upper_error_rp_mcmc/mcmcerrors[0,0],'Beta_RP_Upper (Method 2)',FORMAT='(F9.7,23x,A25)'      ;MCMC
    printf,20, abs(lower_error_rp_mcmc)/mcmcerrors[1,0],'Beta_RP_Lower (Method 2)',FORMAT='(F9.7,23x,A25)' ;MCMC 
	endif
	If parinfo_ar eq 0 then begin
		printf,20,BB_mcmc[2]/MAX([mcmcerrors[0,1],mcmcerrors[1,1]]),'Beta_AR',FORMAT='(F9.7,23x,A8)'	
		printf,20, upper_error_ar_mcmc/mcmcerrors[0,1],'Beta_AR_Upper (Method 2)',FORMAT='(F9.7,23x,A25)'     ;MCMC
    printf,20,abs(lower_error_ar_mcmc)/mcmcerrors[1,1],'Beta_AR_Lower (Method 2)',FORMAT='(F9.7,23x,A25)' ;MCMC
	endif
	If parinfo_inc eq 0 then begin
		printf,20,CC_mcmc[2]/MAX([mcmcerrors[0,2],mcmcerrors[1,2]]),'Beta_IN',FORMAT='(F9.7,23x,A8)'	
		printf,20, upper_error_inc_mcmc/(mcmcerrors[0,2]),'Beta_IN_Upper (Method 2)',FORMAT='(F9.6,23x,A25)'      ;MCMC
    printf,20, abs(lower_error_inc_mcmc)/(mcmcerrors[1,2]),'Beta_IN_Lower (Method 2)',FORMAT='(F9.6,23x,A25)' ;MCMC
	endif
	If parinfo_mid eq 0 then begin
		printf,20,DD_mcmc[2]/MAX([mcmcerrors[0,3],mcmcerrors[1,3]]),'Beta_MID',FORMAT='(F9.7,24x,A8)'	
		printf,20, upper_error_mid_mcmc/mcmcerrors[0,3],'Beta_MID_Upper (Method 2)',FORMAT='(F9.7,24x,A25)'     ;MCMC
    printf,20,abs(lower_error_mid_mcmc)/mcmcerrors[1,3],'Beta_MID_Lower (Method 2)',FORMAT='(F9.7,24x,A25)' ;MCMC 
	endif
	If parinfo_lin eq 0 then begin
		printf,20,EE_mcmc[2]/MAX([mcmcerrors[0,4],mcmcerrors[1,4]]),'Beta_lin',FORMAT='(F9.7,24x,A9)'	
		printf,20, upper_error_lin/mcmcerrors[0,4],'Beta_lin_Upper (Method 2)',FORMAT='(F9.7,24x,A25)'       ;MCMC
    printf,20, abs(lower_error_lin)/mcmcerrors[1,4],'Beta_lin_Lower (Method 2)',FORMAT='(F9.7,24x,A25)'  ;MCMC  	
	endif
	If parinfo_quad eq 0 then begin
		printf,20,FF_mcmc[2]/MAX([mcmcerrors[0,5],mcmcerrors[1,5]]),'Beta_quad',FORMAT='(F9.7,24x,A10)'
		printf,20, upper_error_quad_mcmc/mcmcerrors[0,5],'Beta_quad_Upper (Method 2)',FORMAT='(F9.6,24x,A26)'      ;MCMC   
    printf,20, abs(lower_error_quad_mcmc)/mcmcerrors[1,5],'Beta_quad_Lower (Method 2)',FORMAT='(F9.6,24x,A26)' ;MCMC   		
	endif
	
endif 

printf,20,''
printf,20,'**** Monte Carlo Levenberg-Marquardt Values and Errors******'
;values and errors
	If parinfo_rp eq 0.0 then begin 
		printf,20,Best_fit[0], AA[2],'Rp/R*',FORMAT='(F9.7,2X,F9.7,2x,A5)' 		
	endif

	If parinfo_ar eq 0.0 then begin
		printf,20,Best_fit[1], BB[2],'a/R*',FORMAT='(F9.7,2X,F9.7,2x,A4)'
	endif

	If parinfo_inc eq 0.0 then begin
		printf,20,Best_fit[2]*(180/!PI), CC[2]*(180/!PI), 'Inclination (degrees)',FORMAT='(F9.6,2X,F9.7,2x,A21)'	
	endif

	If parinfo_mid eq 0.0 then begin
		printf,20,Best_fit[3], DD[2],'Mid-Transit (days)',FORMAT='(F10.7,2X,F9.7,1x,A18)'	
	endif 

	If parinfo_lin eq 0 then begin
		printf,20,Best_fit[4], EE[2],'Linear Limb',FORMAT='(F9.7,2X,F9.7,1x,A15)'
	endif 

	If parinfo_quad eq 0 then begin
		printf,20,Best_fit[5], FF[2],'Quadratic Limb',FORMAT='(F9.7,2X,F9.7,2x,A15)'	
	endif 

;Beta's
printf,20,Beta_stuff[0], Beta_stuff[2],'White Noise (LM):Time',FORMAT='(F9.7,2X,F9.7,1x,A23)'
printf,20,Beta_stuff[1], Beta_stuff[3],'Red Noise (LM):Time',FORMAT='(F9.7,2X,F9.7,1x,A20)'
printf,20,Beta_stuff[4],'Beta (LM):Time',FORMAT='(F9.7,12x,A15)'
  printf,20,sigma_wavelet_white_LM,'White Noise (LM):Wavelet',FORMAT='(F9.7,11x,A26)'
  printf,20,sigma_wavelet_red_LM,'Red Noise (LM):Wavelet',FORMAT='(F9.7,10x,A25)'
  printf,20,beta_wavelet_LM,'Beta (LM):Wavelet',FORMAT='(F9.7,5x,A25)'

	If parinfo_rp eq 0 then begin 
		printf,20,AA[2]/A[2],'Beta_RP',FORMAT='(F9.7,12x,A8)	
		printf,20, upper_error_rp/A[2],'Beta_RP_Upper (Method 2)',FORMAT='(F9.7,12x,A25)'          ;LM
    printf,20, abs(lower_error_rp)/A[2],'Beta_RP_Lower (Method 2)',FORMAT='(F9.7,12x,A25)'     ;LM
	endif
	If parinfo_ar eq 0 then begin
		printf,20,BB[2]/B[2],'Beta_AR',FORMAT='(F9.7,12x,A8)	
		printf,20, upper_error_ar/B[2],'Beta_AR_Upper (Method 2)',FORMAT='(F9.7,12x,A25)'         ;LM
    printf,20,abs(lower_error_ar)/B[2],'Beta_AR_Lower (Method 2)',FORMAT='(F9.7,12x,A25)'     ;LM	
	endif
	If parinfo_inc eq 0 then begin
		printf,20,CC[2]/C[2],'Beta_IN',FORMAT='(F9.7,12x,A8)	
		printf,20, upper_error_inc/C[2],'Beta_IN_Upper (Method 2)',FORMAT='(F9.6,12x,A25)'        ;LM  
    printf,20, abs(lower_error_inc)/C[2],'Beta_IN_Lower (Method 2)',FORMAT='(F9.6,12x,A25)'  ;LM
	endif
	If parinfo_mid eq 0 then begin
		printf,20,DD[2]/D[2],'Beta_MID',FORMAT='(F9.7,13x,A8)'
		printf,20, upper_error_mid/D[2],'Beta_MID_Upper (Method 2)',FORMAT='(F9.7,13x,A25)'    ;LM
  printf,20,abs(lower_error_mid)/D[2],'Beta_MID_Lower (Method 2)',FORMAT='(F9.7,13x,A25)'  ;LM	
	endif
	If parinfo_lin eq 0 then begin
		printf,20,EE[2]/E[2],'Beta_lin',FORMAT='(F9.7,12x,A9)	
    printf,20, upper_error_lin/E[2],'Beta_lin_Upper (Method 2)',FORMAT='(F9.7,13x,A25)'      ;LM
  printf,20, abs(lower_error_lin)/E[2],'Beta_lin_Lower (Method 2)',FORMAT='(F9.7,13x,A25)'  ;LM
	endif
	If parinfo_quad eq 0 then begin
		printf,20,FF[2]/F[2],'Beta_quad',FORMAT='(F9.7,13x,A10)	
    printf,20, upper_error_quad/F[2],'Beta_quad_Upper (Method 2)',FORMAT='(F9.6,13x,A26)'         ;LM
    printf,20, abs(lower_error_quad)/F[2],'Beta_quad_Lower (Method 2)',FORMAT='(F9.6,13x,A26) '  ;LM  
	endif


;Other
printf,20,''
printf,20,'****Other Useful Values******'
If run_baseline eq 0 then begin
  printf,20,Value[1],Value[2],' Duration (Days)',FORMAT='(F9.7,2X,F9.7,1x,A16)'  
  printf,20,Value[1]*24.d*60.d,Value[2]*24.d*60.d,' Duration (Min)',FORMAT='(F9.5,2X,F9.7,1x,A15)' 
endif
printf,20,Value[3],' Cadence (days)',FORMAT='(F9.7,12X,A15)'
printf,20, cadence_sec, ' Cadence (secs)' ,FORMAT='(F9.5,3X,A24)'
printf,20,Value[0],'OoT STDEV (mmag)',FORMAT='(F9.7,13x,A16)'
printf,20,chi_reduced,'Reduced Chi-Squared',FORMAT='(F9.5,13X,A19)'
If info_baseline eq 0 then begin
	printf,20, Detrend[0],' Baseline Intercept', FORMAT='(F9.5,12X,A19)'
	printf,20, Detrend[1],' Baseline Slope' , FORMAT='(F9.5,12x,A15)'
		If run_baseline eq 2 then begin
			printf,20, Detrend[2],'Baseline x^2 term' , FORMAT='(F9.5,13x,A17)'
		endif
endif
If run_symmetry eq 0 then begin 
	printf,20,stdev(subtransit[1,*])*1000.d,'STDEV of Asymmetry (mmag)',FORMAT='(F9.7,12X,A26)'
endif 
printf,20,''
printf,20, version,'Version Number', FORMAT='(D5.2,20X,A15)'
printf,20, SYSTIME(),'  Date and Time of Model' 

printf,20,''
printf,20,'****** Errors are not final (i.e. The Betas have not been propagated )****' 
printf,20,'****See the README file for a description on how to properly do this ******'


save,/ALL,FILENAME ='everything.sav'		;save everything
file_move,['everything.sav','FINAL_EVERYTHING.dat'],foldername		;move

;To restore
;restore,FILENAME='everything.sav'
close,/all
print,'The End' 

END
