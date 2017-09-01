

#include <vector>
#include <algorithm>
#include "new_random.h"
#include <limits>

using namespace std;


#define movAvgSDSeries true //save mean model behaviour?
#define ZERO_PROBABILITY (1e-200) //numerical limit for zero probability
#define NONSENSE -999999 //for general use on non-initialized variables
#define VIRTUALLY_ZERO (1.0/10000) //used to trim values to zero
#define VIRTUALLY_ONE (0.999) //used to trim values to one

//#definitions of climate driven ento-epidemiological parameters
#define temp_effect_muA (2.13 -0.3797*T +0.02457*pow(T,2) -0.0006778*pow(T,3) +0.000006794*pow(T,4))
#define temp_effect_epA (0.131 -0.05723*T +0.01164*pow(T,2) -0.001341*pow(T,3) +0.00008723*pow(T,4) -0.000003017*pow(T,5) +0.00000005153*pow(T,6) -0.000000000342*pow(T,7))
#define temp_effect_muV ( 0.8692 -0.159*T +(0.01116*pow(T,2)) -0.0003408*pow(T,3) +0.000003809*pow(T,4) )
#define temp_effect_gammaV (24.0*( 0.003359* (Tk/298.) * exp((15000./R)*(1/298.-1./Tk)) / (1.+ exp((6.203*pow(10.,21)/R)*(1./(-2.176*pow(10.,30))-1./Tk))) )) //its output is hours
#define effect_temp_epsVH (0.001044*T *(T-12.286)*pow((32.461-T),0.5))
#define temp_effect_c ((-1.848e+02 +2.794e+01*T -9.254e-01*pow(T,2) +9.226e-03*pow(T,3))/100.0 )
#define temp_effect_theta  (-5.4 +1.8*T -0.2124*pow(T,2) +0.01015*pow(T,3) -0.0001515*pow(T,4))

//#definitions of initial population sizes of the vector model classes
#define Q ( (epA_t/(epA_t+muA_t)) * (c_t*theta_t)/muV_t )
#define iVeq ((cK*(rain_t+1))*(1.-(1./Q))*(epA_t/muV_t))
#define iA ( iAeq<0 ? (NH0*1) : iAeq ) //it is assumed that the pop is not extinct, even if Temp0 would suggest it
#define iV ( iVeq<0 ? (NH0*1) : iVeq ) //it is assumed that the pop is not extinct, even if Temp0 would suggest it
#define iAeq ((cK*(rain_t+1))*(1.-(1./Q)))

//#weather driven formulations for entomological parameters
#define Za ((humi_t-meanHumi)/(sqrt(1+pow((humi_t-meanHumi),2))))
#define ZaWater (pow((1+Za),cEco2FactorV))
#define Zc ((rain_t-meanPrec)/(sqrt(1+pow((rain_t-meanPrec),2))))
#define ZcWater (pow((1+Zc),cEco2FactorV))
#define Ztheta ((rain_t-meanPrec)/(sqrt(1+pow((rain_t-meanPrec),2))))
#define ZthetaWater 1
#define ZmuA (meanPrec-(rain_t-meanPrec)/(sqrt(1+pow((rain_t-meanPrec),2))))
#define ZmuAWater 1
#define ZmuV (meanHumi-(humi_t-meanHumi)/(sqrt(1+pow((humi_t-meanHumi),2))))
#define ZmuVWater (pow((1+ZmuV),cEco2FactorV))

//# model's classes
#define TT  y[0]
#define SH  y[1]
#define EH  y[2]
#define IH  y[3]
#define RH  y[4] //#cummulative 'incidence' of recovereds
#define A   y[5]
#define SV  y[6]
#define EV  y[7]
#define IV  y[8]
#define RRH y[9] //#actual dynamic class for human recovereds

#define dSH dydx[1]
#define dEH dydx[2]
#define dIH dydx[3]
#define dRH dydx[4]
#define dA  dydx[5]
#define dSV dydx[6]
#define dEV dydx[7]
#define dIV dydx[8]
#define dRRH dydx[9]

#define NEQUATIONS 10 //count all the above

//#input data is stored in modelInputData, define easier names heres
#define sim_t (modelInputData[0])
#define temperature_t (modelInputData[1])
#define epi_cases_cm_t (modelInputData[2])
#define known_epi_points_t (modelInputData[3])
#define prec_t (modelInputData[4])
#define um_t (modelInputData[5])

//#define a set of variables; description is in bash script or function where used
bool static sampleGaussian;
double static SAMPLE_factor_T;
double static SAMPLE_factor_K;
double static SAMPLE_factor_fV;
double static SAMPLE_EipInfH;
double static SAMPLE_factor_zetaRho;
double static SAMPLE_iH0;
double time2run;
double timestart;
double stepsize;
int nsteps;
int nMCMCSteps;
double eta;
double NH;
double NV;
int introStep;
string tag;
string tagInPar;
int countIntros;
double* initialCondition;
bool noZeta;
bool fRho;
int tRho;
double meanPrec;
double meanHumi;
double NH0;
double psi;
double tpsi;
int maxSampleSims;


std::vector<std::vector<double> > modelInputData; //#keeps all input data
std::vector<double> data_odeVsdata; //#keeps all data points to be fitted
std::vector<double> timestep_odeVsdata; //#keeps all time steps to be fitted
std::vector<int> par_varying; //list of parameter indexes for the MCMC to fit
int par_nstep; //#number of parameter to co-jump per MCMC step

//#described in funtion 'main'
double* RH_t;
double* SH_t;
double* inc_t;
double* R0_t;
double* Re_t;
double* V_t;
double* IV_t;
double* A_t;
double* LEV_t;
double* IH_t;
double* LV_t;

//#described in funtion 'main'
double* prev_RH_t;
double* prev_SH_t;
double* prev_inc_t;
double* prev_R0_t;
double* prev_Re_t;
double* prev_V_t;
double* prev_IV_t;
double* prev_A_t;
double* prev_LEV_t;
double* prev_IH_t;
double* prev_LV_t;

//#described in funtion 'main'
double* mean_RH_t;
double* mean_SH_t;
double* mean_inc_t;
double* mean_R0_t;
double* mean_Re_t;
double* mean_V_t;
double* mean_IV_t;
double* mean_A_t;
double* mean_LEV_t;
double* mean_IH_t;
double* mean_LV_t;

//#described in funtion 'main'
double cT0;
double cK;
double cEcoFactorV;
double cEco2FactorV;
double cEipFactorV;
double cEipH;
double cInfpH;
double cZeta;
double cRho;
double cObsR;

//#described in funtion 'main'
double iT0;
double iK;
double iEcoFactorV;
double iEco2FactorV;
double iEipFactorV;
double iEipH;
double iInfpH;
double iZeta;
double iRho;

//#described in funtion 'main'
double* aProb_accept;
double* accept_accept;
double* T0_accept;
double* T0step_accept;
double* K_accept;
double* ecoFactorV_accept;
double* eco2FactorV_accept;
double* eipFactorV_accept;
double* eipH_accept;
double* infpH_accept;
double* zeta_accept;
double* rho_accept;

//#for each time step:
double T; //to be used as temperature in celsius
double Tk; //to be used as kelvin
double c_t; //to be used as hatching success
double muA_t;//to be used as death rate of aquatic vector
double epA_t; //to be used as development rate from aquatic to adult
double muV_t; //to be used as death rate of adult vector
double gammaV_t; //to be used as incubation period of vector
double epsilonVH_t; //to be used as prob of infection from vector to human
double epsilonHV_t; //to be used as prob of infection from human to vector
double theta_t; //to be used as oviposition rate
double rain_t; //to be used as precipitaiton
double humi_t; //to be used as humidity
double a_t; //to be used as bitting rate
double deltaH; //to be used as recovering rate of humans
double gammaH; //to be used as latency-to-infectious rate of humans


//#define file streams for output
ofstream fout_all;
ofstream fout_inc;
ofstream fout_time;
ofstream fout_raink;
ofstream fout_RH;
ofstream fout_SH;
ofstream fout_V;
ofstream fout_SV;
ofstream fout_a;
ofstream fout_r0;
ofstream fout_re;
ofstream fout_LEV;
ofstream fout_A2V;
ofstream fout_LV;
ostringstream fileNameBuilder;

//#function exports the mean model dynamics over all MCMC accepted states
void exportSeries(){

    //#get output files ready
    fileNameBuilder.str(""); //clean filename
    fileNameBuilder << "output/"<<tag<<".inc_t.data"; //define filename of output
    fout_inc.open(fileNameBuilder.str().c_str()); //open the file

    fileNameBuilder.str(""); //clean filename
    fileNameBuilder << "output/"<<tag<<".RH_t.data"; //define filename of output
    fout_RH.open(fileNameBuilder.str().c_str()); //open the file

    fileNameBuilder.str(""); //clean filename
    fileNameBuilder << "output/"<<tag<<".SH_t.data"; //define filename of output
    fout_SH.open(fileNameBuilder.str().c_str()); //open the file

    fileNameBuilder.str(""); //clean filename
    fileNameBuilder << "output/"<<tag<<".V_t.data"; //define filename of output
    fout_V.open(fileNameBuilder.str().c_str()); //open the file

    fileNameBuilder.str(""); //clean filename
    fileNameBuilder << "output/"<<tag<<".IV_t.data"; //define filename of output
    fout_SV.open(fileNameBuilder.str().c_str()); //open the file

    fileNameBuilder.str(""); //clean filename
    fileNameBuilder << "output/"<<tag<<".A_t.data"; //define filename of output
    fout_a.open(fileNameBuilder.str().c_str()); //open the file

    fileNameBuilder.str(""); //clean filename
    fileNameBuilder << "output/"<<tag<<".r0_t.data"; //define filename of output
    fout_r0.open(fileNameBuilder.str().c_str()); //open the file

    fileNameBuilder.str(""); //clean filename
    fileNameBuilder << "output/"<<tag<<".re_t.data"; //define filename of output
    fout_re.open(fileNameBuilder.str().c_str()); //open the file

    fileNameBuilder.str(""); //clean filename
    fileNameBuilder << "output/"<<tag<<".LEV_t.data"; //define filename of output
    fout_LEV.open(fileNameBuilder.str().c_str()); //open the file

    fileNameBuilder.str(""); //clean filename
    fileNameBuilder << "output/"<<tag<<".IH_t.data"; //define filename of output
    fout_A2V.open(fileNameBuilder.str().c_str()); //open the file

    fileNameBuilder.str(""); //clean filename
    fileNameBuilder << "output/"<<tag<<".LV_t.data"; //define filename of output
    fout_LV.open(fileNameBuilder.str().c_str()); //open the file

    //#cycle time and export each variable to a separate file
    for(int t=0; t< nsteps; t=t+1){
        fout_inc << mean_inc_t[t] << " ";
        fout_RH << mean_RH_t[t] << " ";
        fout_SH << mean_SH_t[t] << " ";
        fout_V << mean_V_t[t] << " ";
        fout_SV << mean_IV_t[t] << " ";
        fout_a << mean_A_t[t] << " ";
        fout_LEV << mean_LEV_t[t] << " ";
        fout_A2V << mean_IH_t[t] << " ";
        fout_LV << mean_LV_t[t] << " ";
        fout_r0 << mean_R0_t[t] << " ";
        fout_re << mean_Re_t[t] << " ";
    }

    //#close all files
    fout_RH endL; fout_RH.close();
    fout_inc endL; fout_inc.close();
    fout_raink endL; fout_raink.close();
    fout_V endL; fout_V.close();
    fout_SV endL; fout_V.close();
    fout_a endL; fout_a.close();
    fout_r0 endL; fout_r0.close();
    fout_re endL; fout_re.close();
    fout_LEV endL; fout_LEV.close();
    fout_A2V endL; fout_A2V.close();
    fout_LV endL; fout_LV.close();
}

//#for the ento-epidemiological formulations dependent on climate variables, some
//#of the polinomials are restricted in range; for example, for some temperature values
//#the polinomials may give negative parameter values which in reality should be zero, or give
//#values above one when the parameters should have a ceiling of 1; here, these ranges are trimmed
void fixRanges(){
        epsilonVH_t= (epsilonVH_t  <=0) ? VIRTUALLY_ZERO : epsilonVH_t;
        epsilonVH_t= (epsilonVH_t  >1) ? VIRTUALLY_ONE : epsilonVH_t;
        if(T>31) epsilonVH_t= VIRTUALLY_ZERO; //expression not derived for values >31 (lambrechts2011)

        epsilonHV_t= (epsilonHV_t  <=0) ? VIRTUALLY_ZERO : epsilonHV_t;
        epsilonHV_t= (epsilonHV_t  >1) ? VIRTUALLY_ONE : epsilonHV_t;

        c_t= (c_t<=0) ? VIRTUALLY_ZERO : c_t;
        c_t= (c_t>1) ? VIRTUALLY_ONE : c_t;

        muV_t= (muV_t<=0) ? VIRTUALLY_ZERO : muV_t;
        muV_t= (muV_t>1) ? VIRTUALLY_ONE : muV_t;

        gammaV_t= (gammaV_t<=0) ? VIRTUALLY_ZERO : gammaV_t;
        gammaV_t= (gammaV_t>1) ? VIRTUALLY_ONE : gammaV_t;

        muA_t= (muA_t<=0) ? VIRTUALLY_ZERO : muA_t;
        muA_t= (muA_t>1) ? VIRTUALLY_ONE : muA_t;

        epA_t= (epA_t<=0) ? VIRTUALLY_ZERO : epA_t;
        epA_t= (epA_t>1) ? VIRTUALLY_ONE : epA_t;
        if(T<10) epA_t= VIRTUALLY_ZERO; //expression not derived for values <10 (Yang2009)

        theta_t= (theta_t<=0) ? VIRTUALLY_ZERO : theta_t;
        theta_t= (theta_t>1) ? VIRTUALLY_ONE : theta_t;
}


//#here, the proposed parameter values are transformed / copied into the variables
//#that are actually used in the model's code / equations
void setStepParameters(int istep){

        //#transform / use temperatures
        T= temperature_t[istep];
        Tk= T+273.15;

        //#decide if the observation rate is the first or second
        cObsR= (istep< ((tRho+abs(timestart))/stepsize)) ? cZeta : cRho;

        rain_t= (prec_t[istep]);   //#get precipitation
        humi_t= (um_t[istep]); //get humidity
        a_t= (BITTING_RATE); //get bitting rate
        c_t= (temp_effect_c);  //get hatching
        theta_t= (temp_effect_theta);  //get overposition
        epA_t= (temp_effect_epA); //get rate of aquatic to adult
        epsilonVH_t= (effect_temp_epsVH);  //get trans prob vector to human
        epsilonHV_t= (effect_temp_epsHV);  //get trans prob human to vector
        gammaV_t= (temp_effect_gammaV * cEipFactorV);  //get EIP rate
        muV_t= (temp_effect_muV * cEcoFactorV);  //get mort rate adults
        muA_t= (temp_effect_muA * cEcoFactorV); //get mort rate aquatic

        a_t= (a_t) * ZaWater; //#correct bitting rate according to extra climate factors
        c_t= (c_t) * ZcWater; //#correct hatching according to extra climate factors
        theta_t= (theta_t) * ZthetaWater; //#correct oviposition according to extra climate factors
        muA_t= (muA_t) * ZmuAWater; //#correct aquatic death rate according to extra climate factors
        muV_t= (muV_t) * ZmuVWater; //#correct adult death rate according to extra climate factors

        //#correct polinomial numerical ranges
        fixRanges();

        c_t= F*c_t; //#effective female hatching
        deltaH= (1.0/cInfpH); //#transform period to rate of recovery
        gammaH= (1.0/cEipH); //#transform period to rate of latency-to-infectious

}

//#initialize time series (variables that will keep model's mean time series)
void initTimeSeries(){

    for(int t=0; t< nsteps; t=t+1){
        mean_V_t[t]= 0;
        mean_R0_t[t]= 0;
        mean_Re_t[t]= 0;
        mean_A_t[t]= 0;
        mean_inc_t[t]= 0;
        mean_RH_t[t]= 0;
        mean_SH_t[t]= 0;
        mean_LEV_t[t]= 0;
        mean_IH_t[t]= 0;
        mean_LV_t[t]= 0;

        prev_inc_t[t]= 0;
        prev_RH_t[t]= 0;
        prev_SH_t[t]= 0;
        prev_V_t[t]= 0;
        prev_LEV_t[t]= 0;
        prev_IH_t[t]= 0;
        prev_LV_t[t]=  0;
        prev_R0_t[t]= 0;
        prev_Re_t[t]= 0;
        prev_A_t[t]= 0;
    }

    //#also clean / reset model's classes
    for(int t=0; t< nsteps; t=t+1){
        inc_t[t]= 0;
        RH_t[t]= 0;
        SH_t[t]= 0;
        V_t[t]= 0;
        IV_t[t]= 0;
        LEV_t[t]= 0;
        IH_t[t]= 0;
        LV_t[t]=  0;
        R0_t[t]= 0;
        Re_t[t]= 0;
        A_t[t]= 0;
    }

}

//# save current time series as current mean dynamics and previously
//# accepted dynamics
void saveTimeSeriesOnly(){

    for(int t=0; t< nsteps; t=t+1){
        mean_V_t[t]= V_t[t];
        mean_IV_t[t]= IV_t[t];
        mean_R0_t[t]= R0_t[t];
        mean_Re_t[t]= Re_t[t];
        mean_A_t[t]= A_t[t];
        mean_inc_t[t]= inc_t[t];
        mean_RH_t[t]= RH_t[t];
        mean_SH_t[t]= SH_t[t];
        mean_LEV_t[t]= LEV_t[t];
        mean_IH_t[t]= IH_t[t];
        mean_LV_t[t]= LV_t[t];

        prev_inc_t[t]= inc_t[t];
        prev_RH_t[t]= RH_t[t];
        prev_SH_t[t]= SH_t[t];
        prev_V_t[t]= V_t[t];
        prev_IV_t[t]= IV_t[t];
        prev_LEV_t[t]= LEV_t[t];
        prev_IH_t[t]= IH_t[t];
        prev_LV_t[t]= LV_t[t];
        prev_R0_t[t]= R0_t[t];
        prev_Re_t[t]= Re_t[t];
        prev_A_t[t]= A_t[t];
    }
}

//#save the current model dynamics by calculating new mean model behaviour given the
//#previous mean and the current model dynamics
void saveTimeseries(double n, bool alsoCalcMovAvgs){

    double prevMean;

    //#do only save?
    if(n==0 | !alsoCalcMovAvgs) {
        saveTimeSeriesOnly();
        return;
    }

    //#also update mean behaviour
    n= n+1; //n starts at sim 0!
    for(int t=0; t< nsteps; t=t+1){

        prevMean= mean_inc_t[t];
        mean_inc_t[t]= (1./n)*(inc_t[t]+(n-1.)*prevMean);

        prevMean= mean_RH_t[t];
        mean_RH_t[t]= (1./n)*(RH_t[t]+(n-1.)*prevMean);

        prevMean= mean_SH_t[t];
        mean_SH_t[t]= (1./n)*(SH_t[t]+(n-1.)*prevMean);

        prevMean= mean_R0_t[t];
        mean_R0_t[t]= (1./n)*(R0_t[t]+(n-1.)*prevMean);

        prevMean= mean_Re_t[t];
        mean_Re_t[t]= (1./n)*(Re_t[t]+(n-1.)*prevMean);

        prevMean= mean_V_t[t];
        mean_V_t[t]= (1./n)*(V_t[t]+(n-1.)*prevMean);

        prevMean= mean_IV_t[t];
        mean_IV_t[t]= (1./n)*(IV_t[t]+(n-1.)*prevMean);

        prevMean= mean_LEV_t[t];
        mean_LEV_t[t]= (1./n)*(LEV_t[t]+(n-1.)*prevMean);

        prevMean= mean_IH_t[t];
        mean_IH_t[t]= (1./n)*(IH_t[t]+(n-1.)*prevMean);

        prevMean= mean_LV_t[t];
        mean_LV_t[t]= (1./n)*(LV_t[t]+(n-1.)*prevMean);

        prevMean= mean_A_t[t];
        mean_A_t[t]= (1./n)*(A_t[t]+(n-1.)*prevMean);

        //for next step, make current state the 'next' previous state
        prev_inc_t[t]= inc_t[t];
        prev_RH_t[t]= RH_t[t];
        prev_SH_t[t]= SH_t[t];
        prev_V_t[t]= V_t[t];
        prev_IV_t[t]= IV_t[t];
        prev_LEV_t[t]= LEV_t[t];
        prev_IH_t[t]= IH_t[t];
        prev_LV_t[t]= LV_t[t];
        prev_R0_t[t]= R0_t[t];
        prev_Re_t[t]= Re_t[t];
        prev_A_t[t]= A_t[t];
    }

}

//#save the previously accepted model dynamics by calculating new mean model behaviour given the
//#previous mean and the previously accepted model dynamics
void savePrevTimeseries(double n, bool alsoCalcMovAvgs){

    int t;
    double prevMean;

    //#do only save?
    if(n==0 | !alsoCalcMovAvgs) {
        saveTimeSeriesOnly();
        return;
    }

    //#also update mean behaviour
    n= n+1; //n starts at sim 0!
    for(int t=0; t< nsteps; t= t+1){

        prevMean= mean_inc_t[t];
        mean_inc_t[t]= (1./n)*(prev_inc_t[t]+(n-1.)*prevMean);

        prevMean= mean_RH_t[t];
        mean_RH_t[t]= (1./n)*(prev_RH_t[t]+(n-1.)*prevMean);

        prevMean= mean_SH_t[t];
        mean_SH_t[t]= (1./n)*(prev_SH_t[t]+(n-1.)*prevMean);

        prevMean= mean_R0_t[t];
        mean_R0_t[t]= (1./n)*(prev_R0_t[t]+(n-1.)*prevMean);

        prevMean= mean_Re_t[t];
        mean_Re_t[t]= (1./n)*(prev_Re_t[t]+(n-1.)*prevMean);

        prevMean= mean_IV_t[t];
        mean_IV_t[t]= (1./n)*(prev_IV_t[t]+(n-1.)*prevMean);

        prevMean= mean_V_t[t];
        mean_V_t[t]= (1./n)*(prev_V_t[t]+(n-1.)*prevMean);

        prevMean= mean_LEV_t[t];
        mean_LEV_t[t]= (1./n)*(prev_LEV_t[t]+(n-1.)*prevMean);

        prevMean= mean_IH_t[t];
        mean_IH_t[t]= (1./n)*(prev_IH_t[t]+(n-1.)*prevMean);

        prevMean= mean_LV_t[t];
        mean_LV_t[t]= (1./n)*(prev_LV_t[t]+(n-1.)*prevMean);

        prevMean= mean_A_t[t];
        mean_A_t[t]= (1./n)*(prev_A_t[t]+(n-1.)*prevMean);

    }

}


//#writes the current simulation into a file; this is used to save the first simulation as reference
//#useful to later compare to what mean behaviour the MCMC converged to, compared to the first
//#simulation with initial proposed parameter values
void exportSimulation( double**y ){

    printL "Exporting first sim." endL;
    fileNameBuilder.str(""); //clean filename
    fileNameBuilder << "output/"<<tag<<".firstsim.population.data"; //define filename of output
    fout_all.open(fileNameBuilder.str().c_str()); //open the file
    //#write the header
    fout_all << "time SH EH IH RH RRH A SV EV IV M" endL;
        for(int t=0; t< nsteps; t=t+1){ //cycle time
            for(int yi=0; yi< NEQUATIONS; yi++) fout_all << y[yi][t] << " ";
            fout_all << (y[6][t]+y[7][t]+y[8][t])/(y[1][t]+y[2][t]+y[3][t]+y[9][t]);
            fout_all endL;
        }
    fout_all.close();
}

//#gives used the feedback composed of all the currently proposed parameter values
void printProposedParams(){
    printL "Proposing: cT0 " << cT0 << " cK " << cK
        << " cEcoF " << cEcoFactorV
        << " cEco2F " << cEco2FactorV
        << " cEpiF " << cEipFactorV
        << " cEipH " << cEipH << " cInfpH " << cInfpH
        << " cZeta " << cZeta << " cRho " << cRho
        endL;
}

//#checks if the total population size in the model is correct; this is useful
//#if running the system with large time step for computational reasons, or when
//#compiling the code with math flags that speed computation but switch off
//#numerical rule checking
void testSimulation(double** y, int sim){
    double simNH= SH[nsteps-1] + EH[nsteps-1] + IH[nsteps-1] + RRH[nsteps-1];
    if(round(simNH)!=round(NH0)){
        printW "At sim: "<< sim << " maybe an error in the ODE? step size too small?" endL;
        printW "At the end Nh: " << simNH << " should be: " << NH0 endL;
        printProposedParams();
        exit(-3);
    }
}

//#given a proposed time of introduction (continuous), find the model's time step
//#closer to that proposal (model time step is more discrete)
void inline findAndSetIntroStep(bool show){
    introStep= distance( sim_t.begin(), upper_bound(sim_t.begin(),sim_t.end(), cT0) );
    if(show){
        printL  "Propose intro cT0: " << cT0 <<
                " will at: " << sim_t[introStep] <<
                " pos: " << introStep endL;
        getchar();
    }
    cT0= sim_t[introStep];
}

//#generally find the closest model step of timeV
int inline findStepOfTime(double timeV){
    return(distance( sim_t.begin(), upper_bound(sim_t.begin(),sim_t.end(), timeV) ));
}

//#here, MCMC jump of parameters is executed; N parameters are chosen randomly from
//#the total that co-jump and new values are sampled from respective distributions
void stepParameters(){
    bool stopFlag;
    int nTries;
    int ipar;

    //suffle the list of pars, effectively making sampling of parameters random
    random_shuffle(par_varying.begin(), par_varying.end());

    ipar= 0;
    stopFlag= false;
    for(int i=0; i< par_nstep; i++){

        nTries=100; //#safety, try 100 times to sample, if failed something is wrong
        ipar= par_varying[i];

        if(ipar==1) {
            cT0= NONSENSE;
            while(cT0<=min_T0 || cT0>=max_T0) {
                if(sampleGaussian) cT0= numGen->norGen(iT0, SAMPLE_factor_T);
                else cT0= numGen->cauGen(iT0, SAMPLE_factor_T);
                nTries--; if(nTries<0) { stopFlag=true; break; }
            }
            findAndSetIntroStep(false);
        }
        else if(ipar==2) {
            cK= NONSENSE;
            while(cK<=min_K || cK>=max_K) {
                if(sampleGaussian) cK= numGen->norGen(iK, SAMPLE_factor_K);
                else cK= numGen->cauGen(iK, SAMPLE_factor_K);
                nTries--; if(nTries<0) { stopFlag=true; break; }
            }
         }
        else if(ipar==3) {
            cEcoFactorV= NONSENSE;
            while(cEcoFactorV<=min_ecoFactorV || cEcoFactorV>=max_ecoFactorV){
                if(sampleGaussian) cEcoFactorV= numGen->norGen(iEcoFactorV, SAMPLE_factor_fV);
                else cEcoFactorV= numGen->cauGen(iEcoFactorV, SAMPLE_factor_fV);
                nTries--; if(nTries<0) { stopFlag=true; break; }
             }
         }
        else if(ipar==4) {
            cZeta= NONSENSE;
            while(cZeta<=min_zeta || cZeta>=max_zeta){
                if(sampleGaussian) cZeta= numGen->norGen(iZeta, SAMPLE_factor_zetaRho);
                else cZeta= numGen->cauGen(iZeta, SAMPLE_factor_zetaRho);
                nTries--; if(nTries<0) { stopFlag=true; break; }
             }
         }
        else if(ipar==5) {
            cRho= NONSENSE;
            while(cRho<=min_rho || cRho>=max_rho){
                if(sampleGaussian) cRho= numGen->norGen(iRho, SAMPLE_factor_zetaRho);
                else cRho= numGen->cauGen(iRho, SAMPLE_factor_zetaRho);
                nTries--; if(nTries<0) { stopFlag=true; break; }
             }
        }
        else if(ipar==6) {
            cEipH= NONSENSE;
            while(cEipH<=min_eipH || cEipH>=max_eipH){
                if(sampleGaussian) cEipH= numGen->norGen(iEipH, SAMPLE_EipInfH);
                else cEipH= numGen->cauGen(iEipH, SAMPLE_EipInfH);
                nTries--; if(nTries<0) { stopFlag=true; break; }
             }
         }
        else if(ipar==7) {
            cInfpH= NONSENSE;
            while(cInfpH<=min_infpH || cInfpH>=max_infpH){
                if(sampleGaussian) cInfpH= numGen->norGen(iInfpH, SAMPLE_EipInfH);
                else cInfpH= numGen->cauGen(iInfpH, SAMPLE_EipInfH);
                nTries--; if(nTries<0) { stopFlag=true; break; }
             }

         }
        else if(ipar==8) {
            cEipFactorV= NONSENSE;
            while(cEipFactorV<=min_eipFactorV || cEipFactorV>=max_eipFactorV){
                if(sampleGaussian) cEipFactorV= numGen->norGen(iEipFactorV, SAMPLE_factor_fV);
                else cEipFactorV= numGen->cauGen(iEipFactorV, SAMPLE_factor_fV);
                nTries--; if(nTries<0) { stopFlag=true; break; }
             }
         }
        else if(ipar==9) {
            cEco2FactorV= NONSENSE;
            while(cEco2FactorV<=min_eco2FactorV || cEco2FactorV>=max_eco2FactorV){
                if(sampleGaussian) cEco2FactorV= numGen->norGen(iEco2FactorV, SAMPLE_factor_fV);
                else cEco2FactorV= numGen->cauGen(iEco2FactorV, SAMPLE_factor_fV);
                nTries--; if(nTries<0) { stopFlag=true; break; }
             }
         }
        else{
            printW "got a random par that i dont recognize: " << ipar endL; exit(-1);
        }
        if(stopFlag){
            printW "Tried to find parameter step for X interations but failed." endL;
            printProposedParams();
            exit(-3);
        }
    }
}

//#if MCMC state was not accepted, roll back current proposal to previously
//#accepted parameter values / MCMC state
void inline rollBackCurrentParameters(){
    cT0= iT0;
    cK= iK;
    cEcoFactorV= iEcoFactorV;
    cEco2FactorV= iEco2FactorV;
    cEipFactorV= iEipFactorV;
    cEipH= iEipH;
    cInfpH= iInfpH;
    cZeta= iZeta;
    cRho= iRho;
}

//#save proposed parameter values as accepted
void inline saveParAndErrAP_acceptPar(int isave, double aProb, int isAccepted){

    //#save the current parameters as accepted
    T0_accept[isave]= cT0;
    T0step_accept[isave]= introStep;
    K_accept[isave]= cK ;
    ecoFactorV_accept[isave]= cEcoFactorV ;
    eco2FactorV_accept[isave]= cEco2FactorV ;
    eipFactorV_accept[isave]= cEipFactorV ;
    eipH_accept[isave]= cEipH ;
    infpH_accept[isave]= cInfpH ;
    zeta_accept[isave]= cZeta ;
    rho_accept[isave]= cRho ;

    //save error acceptance probs
    aProb_accept[isave]= aProb;
    accept_accept[isave]= isAccepted;

    //#accept the current parameters by setting current values to
    //#what will be seend as the previously accepted in the next MCMC step
    iT0= cT0;
    iK= cK;
    iEcoFactorV= cEcoFactorV;
    iEco2FactorV= cEco2FactorV;
    iEipFactorV= cEipFactorV;
    iEipH= cEipH;
    iInfpH= cInfpH;
    iZeta= cZeta;
    iRho= cRho;

}

//#reject proposed parameter values, accept previously accepted values and roll back propositions
void inline savePrevParAndErrAP_roolBackPar(int isave, double aProb, int isAccepted){

    //#save the current parameters
    T0_accept[isave]= T0_accept[isave-1];
    T0step_accept[isave]= T0step_accept[isave-1];
    K_accept[isave]= K_accept[isave-1] ;
    ecoFactorV_accept[isave]= ecoFactorV_accept[isave-1] ;
    eco2FactorV_accept[isave]= eco2FactorV_accept[isave-1] ;
    eipFactorV_accept[isave]= eipFactorV_accept[isave-1] ;
    eipH_accept[isave]= eipH_accept[isave-1] ;
    infpH_accept[isave]= infpH_accept[isave-1] ;
    zeta_accept[isave]= zeta_accept[isave-1] ;
    rho_accept[isave]= rho_accept[isave-1] ;

    //save error acceptance probs
    aProb_accept[isave]= aProb;
    accept_accept[isave]= isAccepted;

    //roll back parameters
    cT0= iT0;
    cK= iK;
    cEco2FactorV= iEco2FactorV;
    cEipFactorV= iEipFactorV;
    cEipH= iEipH;
    cInfpH= iInfpH;
    cZeta= iZeta;
    cRho= iRho;

}

//#calculate loglikelihood
long double inline likelihood_gaussPrior(){
    long double like= 0;
    long double prob;
    int tstep;

    //#look at probabilities of all wanted ode / data points
    for(int i=0; i< timestep_odeVsdata.size(); i=i+1){
        tstep= timestep_odeVsdata[i];
        prob= dpoi( RH_t[tstep], data_odeVsdata[i] );
        if(prob< 1e-300) prob= 1e-300; //avoids numerical under/overflows
        like= like + log(prob);
    }

    //take into account the defined priors
    like= like + log(dgaus(cInfpH, prior_infpH_meanGauss, prior_infpH_sdGauss))
               + log(dgaus(cEipH, prior_eipH_meanGauss, prior_eipH_sdGauss))
               + log(dgaus(cEcoFactorV, prior_ecoFactorV_meanGauss, prior_ecoFactorV_sdGauss))
               + log(dgaus(cEipFactorV, prior_eipFactorV_meanGauss, prior_eipFactorV_sdGauss));

    return like;
}

//#calculate loglikelihood ratio
long double inline MH_likelihood_ratio(long double logLp, long double logLa){
    return exp(logLp-logLa);
}

//initiliaze model population classes
void initPopulationStatus(bool checkflag){

    setStepParameters(0); //set current parameter values

    initialCondition[0]= timestart;
    initialCondition[1]= NH0; //iSH
    initialCondition[2]= 0; //iEH
    initialCondition[3]= 0; //iIH
    initialCondition[4]= 0; //iRH
    initialCondition[9]= 0; //iRRH

    initialCondition[5]= iA; //iA

    initialCondition[6]= iV; //iSV
    initialCondition[7]= 0; //iEV
    initialCondition[8]= 0; //iIV

    if(checkflag){
        printW "Population init: " endL;
        for(int x=0; x<NEQUATIONS; x++) printW initialCondition[x] endL;
        printW "-------" endL;
        getchar();
    }

}

//#function that executs the deterministic verion of model
void modelDeterministic(double x, double y[], double dydx[], int istep){

        double betaVH;
        double betaHV;
        double lamVH;
        double lamHV;
        double migrationI;

        //#population sizes
        NH= SH+EH+IH+RRH;
        NV= SV+EV+IV;

        setStepParameters(istep); //#set parameters to this time step values

        //#locally transform some parameters
        betaVH= (a_t*epsilonVH_t);
        betaHV= (a_t*epsilonHV_t);
        lamVH= (betaVH*IV/NH);
        lamHV= (betaHV*IH/NH);

        //#make decision on migration of infected individuals
        if(x<tpsi) migrationI= 0;
        else{
            migrationI= 0;
            if(numGen->uniGen()<=psi){
                migrationI= 1;
                countIntros++;
              }
        }

        //# note that introduction of the first cases at time of introduction
        //# is executed within the solver code, in file new_solvers_c.h

        //## system dynamics
        dSH= muH*(SH+EH+IH+RRH) -lamVH*SH -muH*SH -(migrationI);
        dEH= lamVH*SH - gammaH*EH -muH*EH;
        dIH= gammaH*EH - deltaH*IH -muH*IH +(migrationI);
        dRH= cObsR*(lamVH*SH);
        dRRH= deltaH*IH -muH*RRH;
        dA= (c_t*theta_t)*(1.0-A/(cK*(rain_t+1)))*NV-(epA_t+(muA_t))*A;
        dSV= (epA_t*A) -(lamHV*SV) -(muV_t*SV) ;
        dEV= (lamHV*SV) - (gammaV_t*EV) - (muV_t*EV);
        dIV= (gammaV_t*EV)- (muV_t*IV) ;

        //#save some realtime values of interest to be used in other code
        //#to export and/or calculate mean model behaviour
        inc_t[istep]= cObsR*(lamVH*SH);
        SH_t[istep]= SH;
        V_t[istep]= NV/NH;
        IV_t[istep]= IV;
        A_t[istep]= A;
        R0_t[istep]= (V_t[istep]*betaVH*betaHV*gammaV_t*gammaH)/(muV_t*(deltaH+muH)*(gammaH+muH)*(gammaV_t+muV_t));
        Re_t[istep]= ((SH/NH)*(SV/NH)*betaVH*betaHV*gammaV_t*gammaH)/(muV_t*(deltaH+muH)*(gammaH+muH)*(gammaV_t+muV_t));
        LEV_t[istep]= 1.0/muV_t;
        LV_t[istep]= 1.0/gammaV_t;
        IH_t[istep]= IH;
        RH_t[istep]= RH;

}

//#function that executs the deterministic verion of model
void modelStochastic(double x, double y[], double dydx[], int istep){

        double betaVH;
        double betaHV;
        double lamVH;
        double lamHV;
        double rate;
        double migrationI;

        //#population sizes
        NH= SH+EH+IH+RRH;
        NV= SV+EV+IV;

        setStepParameters(istep);//#set parameters to this time step values

        //#locally transform some parameters
        betaVH=(a_t*epsilonVH_t);
        betaHV=(a_t*epsilonHV_t);
        lamVH= (betaVH*IV/NH);
        lamHV= (betaHV*IH/NH);

        //#make decision on migration of infected individuals
        if(x<tpsi) migrationI= 0;
        else{
            migrationI= 0;
            if(numGen->uniGen()<=psi){
                migrationI= 1;
                countIntros++;
            }
        }

        //# note that introduction of the first cases at time of introduction
        //# is executed within the solver code, in file new_solvers_c.h

        //decisions on stochastic transitions
        int S_dec;
        int S_death;
        int S_2_exposed;
            rate=(lamVH+muH)*stepsize;
            if(rate>1) { printW "rate SH:" << rate << "you must decrement step size!" endL; exit(-1); }
            S_dec= numGen->binGen(SH,rate);
            if(S_dec>SH) S_dec= SH;
            rate= (lamVH/(lamVH+muH));
            if(rate>1) { printW "rate SH2:" << rate << "you must decrement step size!" endL; exit(-1); }
            S_2_exposed= numGen->binGen(S_dec,rate);
            S_death= (S_dec- S_2_exposed);

        int E_dec;
        int E_death;
        int E_2_infectious;
            rate=(gammaH+muH)*stepsize;
            if(rate>1) { printW "rate EH:" << rate << "you must decrement step size!" endL; exit(-1); }
            E_dec= numGen->binGen(EH,rate);
            if(E_dec>EH) E_dec= EH;
            rate= (gammaH/(gammaH+muH));
            if(rate>1) { printW "rate EH2:" << rate << "you must decrement step size!" endL; exit(-1); }
            E_2_infectious= numGen->binGen(E_dec,rate);
            E_death= (E_dec- E_2_infectious);

        int I_dec;
        int I_death;
        int I_2_recovered;
            rate=(deltaH+muH)*stepsize;
            if(rate>1) { printW "rate IH:" << rate << "you must decrement step size!" endL; exit(-1); }
            I_dec= numGen->binGen(IH,rate);
            if(I_dec>IH) I_dec= IH;
            rate= (deltaH/(deltaH+muH));
            if(rate>1) { printW "rate IH2:" << rate << "you must decrement step size!" endL; exit(-1); }
            I_2_recovered= numGen->binGen(I_dec,rate);
            I_death= (I_dec- I_2_recovered);

        int RR_death;
            rate=(muH)*stepsize;
            if(rate>1) { printW "rate RRH:" << rate << "you must decrement step size!" endL; exit(-1); }
            RR_death= numGen->binGen(RRH,rate);
            if(RR_death>RRH) RR_death= RRH;

        int A_dec;
        int A_2_adult;
        int A_death;
            rate=(epA_t+muA_t)*stepsize;
            if(rate>1) { printW "rate A1:" << rate << "you must decrement step size!" endL; exit(-1); }
            A_dec= numGen->binGen(A,rate);
            if(A_dec>A) A_dec= A;
            rate= ( epA_t / (epA_t+muA_t) );
            if(rate>1) { printW "rate A2:" << rate << "you must decrement step size!" endL; exit(-1); }
            A_2_adult= numGen->binGen(A_dec,rate);
            A_death= (A_dec- A_2_adult);

        int A_births;
            rate=(c_t*theta_t)*stepsize;
            if(rate>1) { printW "rate X:" << rate << "you must decrement step size!" endL; exit(-1); }
            A_births= numGen->binGen(NV,rate);

        int A_deathK;
            rate=(c_t*theta_t*NV/(cK*(rain_t+1)))*stepsize;
            if(rate>1) { printW "rate X:" << rate << "you must decrement step size!" endL; exit(-1); }
            A_deathK= numGen->binGen(A,rate);

        int SV_dec;
        int SV_death;
        int SV_to_exposed;
            rate=(lamHV+muV_t)*stepsize;
            if(rate>1) { printW "rate SV1:" << rate << "you must decrement step size!" endL; exit(-1); }
            SV_dec= numGen->binGen(SV,rate);
            if(SV_dec>SV) SV_dec= SV;
            rate= (lamHV/(lamHV+muV_t));
            if(rate>1) { printW "rate SV2:" << rate << "you must decrement step size!" endL; exit(-1); }
            SV_to_exposed= numGen->binGen(SV_dec,rate);
            SV_death= (SV_dec- SV_to_exposed);

        int EV_dec;
        int EV_death;
        int EV_to_infected;
            rate=(gammaV_t+muV_t)*stepsize;
            if(rate>1) { printW "rate EV1:" << rate << "you must decrement step size!" endL; exit(-1); }
            EV_dec= numGen->binGen(EV,rate);
            if(EV_dec>EV) EV_dec= EV;
            rate= (gammaV_t/(gammaV_t+muV_t));
            if(rate>1) { printW "rate EV2:" << rate << "you must decrement step size!" endL; exit(-1); }
            EV_to_infected= numGen->binGen(EV_dec,rate);
            EV_death= (EV_dec- EV_to_infected);


        int IV_death;
            rate=(muV_t)*stepsize;
            if(rate>1) { printW "rate IV:" << rate << "you must decrement step size!" endL; exit(-1); }
            IV_death= numGen->binGen(IV,rate);

        //model dynamics
        dSH= (S_death+E_death+I_death+RR_death) -S_2_exposed -S_death -migrationI;
        dEH= S_2_exposed - E_2_infectious -E_death;
        dIH= E_2_infectious - I_2_recovered -I_death +migrationI;
        dRH= cObsR*S_2_exposed;
        dRRH= I_2_recovered- RR_death;
        dA= A_births -A_deathK -A_2_adult -A_death;
        dSV= A_2_adult -SV_to_exposed -SV_death ;
        dEV= SV_to_exposed -EV_to_infected -EV_death;
        dIV= EV_to_infected -IV_death ;

        //force persistence of vectors as safety (has no impact on dynamics)
        if( (A+dA)<=0 ) dA= -A+1;
        if( (SV+dSV)<=0 ) dSV= -SV+1;

        //#save some realtime values of interest to be used in other code
        //#to export and/or calculate mean model behaviour
        inc_t[istep]= cObsR*(S_2_exposed)*1/stepsize;
        SH_t[istep]= SH;
        V_t[istep]= NV/NH;
        A_t[istep]= A;
        R0_t[istep]= (V_t[istep]*betaVH*betaHV*gammaV_t*gammaH)/(muV_t*(deltaH+muH)*(gammaH+muH)*(gammaV_t+muV_t));
        Re_t[istep]= ((SH/NH)*(SV/NH)*betaVH*betaHV*gammaV_t*gammaH)/(muV_t*(deltaH+muH)*(gammaH+muH)*(gammaV_t+muV_t));
        LEV_t[istep]= 1.0/muV_t;
        IV_t[istep]= IV;
        LV_t[istep]= 1.0/gammaV_t;
        IH_t[istep]= IH;
        RH_t[istep]= RH;
}
