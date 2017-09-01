

//# make sure you add the ".0" to force float types

//# limits for the carrying capacity (limits of uninformative prior)
#define max_K (5700000)
#define min_K (10000)

//# limits for the linear factor alpha (limits of prior)
#define max_ecoFactorV	10.0
#define min_ecoFactorV	0.1
#define prior_ecoFactorV_meanGauss 3.0 //prior mean gaussian
#define prior_ecoFactorV_sdGauss 0.33 //prior SD gaussian

//# limits for the linear factor eta (limits of prior)
#define max_eipFactorV	10.0
#define min_eipFactorV	0.1
#define prior_eipFactorV_meanGauss 2.0 //prior mean gaussian
#define prior_eipFactorV_sdGauss 0.33  //prior SD gaussian

//# limits for non-linear factor rho (limits of uninformative prior)
#define max_eco2FactorV	3.0
#define min_eco2FactorV	0.001

//# limits for the human incubation period (limits of prior)
#define max_eipH	20.0
#define min_eipH	1.0
#define prior_eipH_meanGauss 6.5//prior mean gaussian
#define prior_eipH_sdGauss 0.33 //prior SD gaussian

//# limits for the human infectious period (limits of prior)
#define max_infpH	20.0
#define min_infpH	1.0
#define prior_infpH_meanGauss 6.0 //prior mean gaussian
#define prior_infpH_sdGauss 0.33  //prior SD gaussian

//# limits for first observation rate (limits of uninformative prior)
#define max_zeta	0.99
#define min_zeta	0.00001

//# limits for time of introduction (limits of uninformative prior)
#define max_T0	200.0
#define min_T0	-200.0

//# limits for second observation rate (limits of uninformative prior)
#define max_rho	0.99
#define min_rho	0.00001

//# fixed parameters
#define BITTING_RATE 0.25 //#adult mosquito bitting rate
#define INTRO_CASESH 2.0 //#cases introduced at time of introduction
#define muH (1.0/(365*70)) //#human life-span
#define R  1.987 //scoffield original article in cal deg-1 mol-1
#define F 0.5 //proportion of eggs hatched that are female (mosquito sex ratio)
#define effect_temp_epsHV (0.5) //probability of transmission per bite from human to vector
