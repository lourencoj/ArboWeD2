
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>
#include <vector>
#include <iterator>
#include <iomanip>

#include "new_2factors_kftemp_param_limits.h"
#include "new_printing.h"
#include "new_printing.h"
#include "new_2factors_kftemp_system.h"
#include "new_solvers_c.h"
#include "new_random.h"

using namespace std;

//## function that splits a verctor of strings given a delim character as separator
std::vector<double> &split(const string &s, char delim, std::vector<double> &elems) {
    stringstream ss(s);
    string item;
    while (std::getline(ss, item, delim)) elems.push_back(atof(item.c_str()));
    return elems;
}

//## function reads the input data and stores it in a vector of vectors
//## each added vector corresponds to one of the input series
void readDataInput(){

	string temp;
	double numbertemp;
	int itstep= 0;

	modelInputData.push_back(std::vector<double>()); //time series, A.K.A: sim_t
	modelInputData.push_back(std::vector<double>()); //temperature series, A.K.A: temperature_t
	modelInputData.push_back(std::vector<double>()); //cummulative incidence series, A.K.A: epi_cases_cm_t
	modelInputData.push_back(std::vector<double>()); //observed time points series, A.K.A: known_epi_points_t
	modelInputData.push_back(std::vector<double>()); //precipitation series, A.K.A: prec_t
	modelInputData.push_back(std::vector<double>()); //humidity series, A.K.A: um_t

	fileNameBuilder.str(""); //clean filename
	fileNameBuilder << "input/" << tag << ".4MCMC.data"; //define filename of output
	//could make a test if file was openned properly
	printL "Opening file " << fileNameBuilder.str() endL;
	fstream myfile1(fileNameBuilder.str().c_str(), std::ios_base::in);
  //#read file and populate each of the input series into separate vectors
	while (getline(myfile1, temp)) {
	    std::vector<double> elems;
	    split(temp, ' ', elems);
	    sim_t.push_back(elems[0]);
	    temperature_t.push_back(elems[1]);
	    epi_cases_cm_t.push_back(elems[2]);
	    known_epi_points_t.push_back(elems[3]);
	    prec_t.push_back(elems[4]);
	    um_t.push_back(elems[5]);
	    itstep++;
	}

  fileNameBuilder.str(""); //clean filename
	fileNameBuilder << "input/" << tag << ".4MCMC.precDist.data"; //define filename of output
	//could make a test if file was openned properly
	printL "Opening file " << fileNameBuilder.str() endL;
	fstream myfile2(fileNameBuilder.str().c_str(), std::ios_base::in);
  //#read file and populate each of the input series into separate vectors
	while (getline(myfile2, temp)) {
	    std::vector<double> elems;
	    split(temp, ' ', elems);
	    meanPrec= elems[0];
	    meanHumi= elems[1];
	    break;
	}
  printL "Read meanPrec " << meanPrec << " and meanHumi " << meanHumi endL;

  //#for the time points in which it is known that data was observed
  //#populate a new vector which memorizes the cummulative incidence
  //#at those time points; this data, not the interpolated is latter fitted
	for(int i=0; i<known_epi_points_t.size(); i++){
		if(known_epi_points_t[i]==1){
			//data point that exists in the original dataset (observed)
			data_odeVsdata.push_back(epi_cases_cm_t[i]);
			timestep_odeVsdata.push_back(i);
		}
	}

  //## time to run the model is determined by the input in the time series
	//## get time limits from the data and calculate number of time steps
	timestart= sim_t[0];
	time2run= sim_t[sim_t.size()-1];
	nsteps= sim_t.size();

  //## give message to the user on what was read about time
	printL "From the data and parameters, i think i have:" endL;
	printL "\t time: " << timestart << " to " << time2run <<" nsteps: " << nsteps << " step: " << sim_t[1]-sim_t[0] endL;

}

//# this function runs the alternative simulation mode (as opposed to MCMC fitting)
//# in here, a previous solution (posteriors) is read from a file, and with a sampling
//# probability, states of the MCMC chain (accepted into posteriors) are used to simulate
//# the model; simulations can be stochastic or deterministic; output files are named accordingly
//# to which type of simulation is performed; output files contain solutions (series) of all relevant
//# model classes, in which each row is a simulation using a sample of the posteiors
void doSimulationMode(int odeType, double sampleSims){

  string temp; //used for reading MCMC solutions

	printL "Entering sim mode with odeType " << odeType endL;

  //#create list of output files before simulation
  ofstream fout_incSim;
  fileNameBuilder.str(""); //clean filename
  if(odeType==1)fileNameBuilder << "output/det_"<<tag<<".inc_t_SIMS.csv"; //define filename of output
  if(odeType==2)fileNameBuilder << "output/sto_"<<tag<<".inc_t_SIMS.csv"; //define filename of output
  fout_incSim.open(fileNameBuilder.str().c_str()); //open the file

  ofstream fout_timeSim;
  fileNameBuilder.str(""); //clean filename
  if(odeType==1)fileNameBuilder << "output/det_"<<tag<<".time_SIMS.csv"; //define filename of output
  if(odeType==2)fileNameBuilder << "output/sto_"<<tag<<".time_SIMS.csv"; //define filename of output
  fout_timeSim.open(fileNameBuilder.str().c_str()); //open the file

  ofstream fout_RHSim;
  fileNameBuilder.str(""); //clean filename
  if(odeType==1)fileNameBuilder << "output/det_"<<tag<<".RH_t_SIMS.csv"; //define filename of output
  if(odeType==2)fileNameBuilder << "output/sto_"<<tag<<".RH_t_SIMS.csv"; //define filename of output
  fout_RHSim.open(fileNameBuilder.str().c_str()); //open the file

  ofstream fout_SHSim;
  fileNameBuilder.str(""); //clean filename
  if(odeType==1)fileNameBuilder << "output/det_"<<tag<<".SH_t_SIMS.csv"; //define filename of output
  if(odeType==2)fileNameBuilder << "output/sto_"<<tag<<".SH_t_SIMS.csv"; //define filename of output
  fout_SHSim.open(fileNameBuilder.str().c_str()); //open the file

  ofstream fout_r0Sim;
  fileNameBuilder.str(""); //clean filename
  if(odeType==1)fileNameBuilder << "output/det_"<<tag<<".r0_t_SIMS.csv"; //define filename of output
  if(odeType==2)fileNameBuilder << "output/sto_"<<tag<<".r0_t_SIMS.csv"; //define filename of output
  fout_r0Sim.open(fileNameBuilder.str().c_str()); //open the file

  ofstream fout_reSim;
  fileNameBuilder.str(""); //clean filename
  if(odeType==1)fileNameBuilder << "output/det_"<<tag<<".re_t_SIMS.csv"; //define filename of output
  if(odeType==2)fileNameBuilder << "output/sto_"<<tag<<".re_t_SIMS.csv"; //define filename of output
  fout_reSim.open(fileNameBuilder.str().c_str()); //open the file

  ofstream fout_VSim;
  fileNameBuilder.str(""); //clean filename
  if(odeType==1)fileNameBuilder << "output/det_"<<tag<<".V_t_SIMS.csv"; //define filename of output
  if(odeType==2)fileNameBuilder << "output/sto_"<<tag<<".V_t_SIMS.csv"; //define filename of output
  fout_VSim.open(fileNameBuilder.str().c_str()); //open the file

  ofstream fout_SVSim;
  fileNameBuilder.str(""); //clean filename
  if(odeType==1)fileNameBuilder << "output/det_"<<tag<<".IV_t_SIMS.csv"; //define filename of output
  if(odeType==2)fileNameBuilder << "output/sto_"<<tag<<".IV_t_SIMS.csv"; //define filename of output
  fout_SVSim.open(fileNameBuilder.str().c_str()); //open the file

  ofstream fout_LEVSim;
  fileNameBuilder.str(""); //clean filename
  if(odeType==1)fileNameBuilder << "output/det_"<<tag<<".LEV_t_SIMS.csv"; //define filename of output
  if(odeType==2)fileNameBuilder << "output/sto_"<<tag<<".LEV_t_SIMS.csv"; //define filename of output
  fout_LEVSim.open(fileNameBuilder.str().c_str()); //open the file

  ofstream fout_LVSim;
  fileNameBuilder.str(""); //clean filename
  if(odeType==1)fileNameBuilder << "output/det_"<<tag<<".LV_t_SIMS.csv"; //define filename of output
  if(odeType==2)fileNameBuilder << "output/sto_"<<tag<<".LV_t_SIMS.csv"; //define filename of output
  fout_LVSim.open(fileNameBuilder.str().c_str()); //open the file

  //# get ready to read from a particular MCMC solution file
  fileNameBuilder.str(""); //clean filename
	fileNameBuilder << "output/"<<tagInPar;
	fstream fin_accepted(fileNameBuilder.str().c_str(), std::ios_base::in);
	printL "opening file with MCMC solutions: " << fileNameBuilder.str().c_str() endL;

	//# retrict output amount by sampling series; necessary or else output is too big
	int sampleOut= 7*int(1.0/stepsize); //this gives output per X days (where X*1/stepsize)

  //# some auxiliar variables
  int read=0;
	int wroteTime=0;
  int simulationi=0;
	int runSims= 0;
	countIntros= 0;

  //# as defined by user, set a maximum number of sims to perform, independently of
  //# sampling probability given
  int NForcedBreak= maxSampleSims;

  //# read MCMC solution file, line by line
	while (getline(fin_accepted, temp)) {

		if(read==0) { read++; continue;} //first line is the header; reject
		read++;

    //# digest line into separate parameters in the MCMC solution
		if(numGen->uniGen() <= sampleSims){
        simulationi++;
        printW "read " << read << " simulation " << simulationi endL;
		    std::vector<double> elems;
		    split(temp, ' ', elems);
  			cT0= elems[0];
  			introStep= elems[1];
  			cK= elems[2];
		    cEcoFactorV= elems[3];
		    cEipFactorV= elems[4];
  			cZeta= elems[5];
  			//6= accept rate; ignored for simulation mode
  			//7= accept status; ignored for simulation mode
  			cRho= elems[8];
		    cEipH= elems[9];
		    cInfpH= elems[10];
		    cEco2FactorV= elems[11];

        //# initialize population before starting simulation
        initPopulationStatus(false);
        //# run simulation depending on sim type
        if(odeType==1) rkdumb(initialCondition, NEQUATIONS, timestart, time2run, nsteps);
        else stochBindumb(initialCondition, NEQUATIONS, timestart, time2run, nsteps);
        //write simulation solutions to output files; each class separately
		    for(int t=0; t< nsteps; t=t+sampleOut){
		        fout_incSim << inc_t[t] << " ";
		        fout_RHSim << RH_t[t] << " ";
		        fout_SHSim << SH_t[t] << " ";
		        fout_r0Sim << R0_t[t] << " ";
		        fout_reSim << Re_t[t] << " ";
		        fout_VSim << V_t[t] << " ";
		        fout_LEVSim << LEV_t[t] << " ";
		        fout_LVSim << LV_t[t] << " ";
		        fout_SVSim << IV_t[t] << " ";
		    }
        //# break row just written; each row is a simulation
        fout_incSim endL;
        fout_RHSim endL;
        fout_SHSim endL;
        fout_r0Sim endL;
        fout_reSim endL;
        fout_VSim endL;
        fout_SVSim endL;
        fout_LEVSim endL;
        fout_LVSim endL;
        //# if this is the first simulation, output a time series, shared by all simulations
        if(wroteTime==0){
            for(int t=0; t< nsteps; t=t+sampleOut) fout_timeSim << y[0][t] << " ";
            fout_timeSim endL;
            wroteTime=1;
        }
        NForcedBreak--;
        runSims++;
		}
    //# break simulation cycle if reached maximum of samples wanted
    if(NForcedBreak==0) break;
	}

  //# close output files
	fout_timeSim.close();
	fout_incSim.close();
	fout_RHSim.close();
	fout_SHSim.close();
	fout_r0Sim.close();
	fout_reSim.close();
	fout_VSim.close();
	fout_SVSim.close();
	fout_LEVSim.close();
	fout_LVSim.close();

  //# inform user of number of introductions per year (debuu)
	printW "intros per year (avg across all rec sims): " << (countIntros/runSims)/((time2run-tpsi)/365.0) endL;

  printL "All simulations are done. Exit." endL;
	exit(-1);
}

//# MAIX c function
int main(int argc, char *argv[]) {

	int nMCMCStepsAcceptedObs; //memorizes how many MCMC states have been accepted after burnin
	int nMCMCStepsAccepted; //memorizes how many MCMC states have been accepted
	int nMCMCStepsSaved; //memorizes how many MCMC states have been recorded
	int recordSeriesThres; //memorizes MCMC step at which burnin ends

	long double logLa; //loglikelihood of currently accepted chain state
	long double logLp; //loglikelihood of (new) proposed chain state
	long double randGuess; //auxiliar for flipping coin to assess probability events
	long double LH_ap; //likelihood of acceptance of proposed chain state
	long double aProb; //auxiliar for probability of acceptance of proposed chain state

  //# read parameters from the bash script
  //# parameters are described in the bash script

	stepsize= atof(argv[1]);
	double sampleSims= atof(argv[2]);
	nMCMCSteps= atoi(argv[3]);
	NH0= atof(argv[4]);
	eta= atof(argv[5]);
	numGen= new NumGenerator(atoi(argv[6]));
	tag= string(argv[7]);
  double sT0= atof(argv[8]);
  double sK= atof(argv[9]);
  double sEcoFactorV= atof(argv[10]);
  int odeType= atoi(argv[11]);
  double sZeta= atof(argv[12]);
  double sRho= atof(argv[13]);
  tRho= atoi(argv[14]);
	double sEipH= atof(argv[15]);
	noZeta= bool(atoi(argv[16]));
	double sInfpH= atof(argv[17]);
	double sEipFactorV= atof(argv[18]);
	double sEco2FactorV= atof(argv[19]);
	bool simulationMode= atoi(argv[20]);
	sampleGaussian= bool(atoi(argv[21]));
	SAMPLE_factor_T= atof(argv[22]);
	SAMPLE_factor_K= atof(argv[23]);
	SAMPLE_factor_fV= atof(argv[24]);
	SAMPLE_EipInfH= atof(argv[25]);
	SAMPLE_factor_zetaRho= atof(argv[26]);
	tagInPar= string(argv[27]);
  psi=  atof(argv[28])/(365.0/stepsize);
  tpsi= atof(argv[29]);
  maxSampleSims= atoi(argv[30]);

	//# give feedback to the user on critical choices made

  if(odeType==1) printW "will be using the deterministic system" endL;
  else if(odeType==2) printW "will be using the stochastic system" endL;
  else { printW "not sure what type of ODE i should use now? Abort." endL; exit(-3); }

  if(sampleGaussian) printW "Sampling is set to Gaussian." endL;
  else printW "Sampling is set to Cauchy." endL;

	readDataInput(); //see function for description
	printL "Have read the input data!" endL;

  //# make a decision on whether the model should use a second observation rate
	if(tRho/stepsize> nsteps){
		printW "Rho is not to be used: cautiously setting sRho=Zeta and fRho=false." endL;
		sRho= sZeta;
		fRho= false;
	} else fRho= true;

	//# vectors that will store the model solutions
	inc_t= new double[nsteps]; 		//incidence
	RH_t= new double[nsteps]; 		//recovered
	SH_t= new double[nsteps];    //susceptible
	R0_t= new double[nsteps]; 		//R0
	Re_t= new double[nsteps]; 		//Re
	A_t= new double[nsteps];     //aquatic vectors
	V_t= new double[nsteps]; 		//adult vectors
	IV_t= new double[nsteps]; 		//adult infected vectors
	LEV_t= new double[nsteps]; 	//life-span of adult vectors
	IH_t= new double[nsteps]; 	//time period from aquatic to adult
	LV_t= new double[nsteps]; 	//EIP of adult vectors

  //# these are used to calculate the mean behaviour of the model givel all
  //# accepted MCMC chain states; this avoids having output for all accepted
  //# states' simulations and only the mean is exported to output

    //series of the simulation last accepted by the MCMC
  	prev_inc_t= new double[nsteps]; 		//incidence
  	prev_RH_t= new double[nsteps]; 		//recovered
  	prev_SH_t= new double[nsteps];    //susceptible
  	prev_R0_t= new double[nsteps]; 		//R0
  	prev_Re_t= new double[nsteps]; 		//Re
  	prev_A_t= new double[nsteps];      //aquatic vectors
  	prev_LEV_t= new double[nsteps]; 	//adult vectors
  	prev_IH_t= new double[nsteps]; 		//adult infected vectors
  	prev_LV_t= new double[nsteps]; 		//life-span of adult vectors
  	prev_V_t= new double[nsteps]; 		//time period from aquatic to adult
  	prev_IV_t= new double[nsteps]; 		//EIP of adult vectors

    //actual mean of all accepted MCMC states' simulations
		mean_inc_t= new double[nsteps]; 			//incidence
		mean_RH_t= new double[nsteps]; 		//recovered
		mean_SH_t= new double[nsteps];     //susceptible
		mean_R0_t= new double[nsteps]; 		//R0
		mean_Re_t= new double[nsteps]; 		//Re
		mean_A_t= new double[nsteps];      //aquatic vectors
		mean_LEV_t= new double[nsteps]; 	//adult vectors
		mean_IH_t= new double[nsteps]; 		//adult infected vectors
		mean_LV_t= new double[nsteps]; 		//life-span of adult vectors
		mean_V_t= new double[nsteps]; 		//time period from aquatic to adult
		mean_IV_t= new double[nsteps]; 		//EIP of adult vectors


  //# create vector of initial conditions for the model
	initialCondition= new double[NEQUATIONS]; //keep ODE init pop conditions

  //# initialize the variables set above
	if(movAvgSDSeries) initTimeSeries();

  //# these vectors keep all the accepted states (to later get posteriors)
	accept_accept= new double[nMCMCSteps]; //acceptance rate along chain
	T0_accept= new double[nMCMCSteps]; //accepted times of introduction
	T0step_accept= new double[nMCMCSteps]; //accepted model step for times of introduction
	K_accept= new double[nMCMCSteps]; //accepted carrying capacity
	ecoFactorV_accept= new double[nMCMCSteps]; //accepted linear factors alpha
	eco2FactorV_accept= new double[nMCMCSteps]; //accepted non-linear factors rho
	eipFactorV_accept= new double[nMCMCSteps]; //accepted linear factors eta
	eipH_accept= new double[nMCMCSteps]; //accepted human incubation period
	infpH_accept= new double[nMCMCSteps]; //accepted human infectious period
	zeta_accept= new double[nMCMCSteps]; //accepted observation rate (first)
	aProb_accept= new double[nMCMCSteps]; //accepted acceptance probability along the chain
	rho_accept= new double[nMCMCSteps]; //accepted observation rate (second)

  //# initialize parameters for MCMC
  //# notation 'i' is used as currently accepted state
  //# notation 'c' is used as currently proposed state
  cT0= sT0;
 	cK= sK;
  cEcoFactorV= sEcoFactorV;
  cEco2FactorV= sEco2FactorV;
  cZeta= sZeta;
  cRho= sRho;
  cEipH= sEipH;
  cInfpH= sInfpH;
  cEipFactorV= sEipFactorV;

  //# the numbers in comment below are references for picking which parameters
  //# will later be used in the MCMC chains
  iT0= cT0; //#1
  iK= cK; //#2
  iEcoFactorV= cEcoFactorV; //#3
  iZeta= cZeta; //#4
  iRho= cRho; //#5
  iEipH= cEipH; //#6
  iInfpH= cInfpH; //#7
  iEipFactorV= cEipFactorV; //#8
  iEco2FactorV= cEco2FactorV; //#9

  //# make decisions on what parameters will be estimated by the MCMC
  //# if any parameter is not to be used by the MCMC, here is the place to skip it
  int ppns= 3; //sets the number of parameters to co-vary in every sampling MCMC jump
  if(fRho){
        //use second observation rate
        if(noZeta){
          //do not use first observation rate
            par_varying = {1,2,3, 5,6,7,8,9};
            par_nstep= ppns;
        }
        else{
          //use first observation rate
            par_varying = {1,2,3,4,5,6,7,8,9};
            par_nstep= ppns;
        }
	}else {
        //do not use second observation rate
    		if(noZeta){
          //do not use first observation rate
      			par_varying = {1,2,3,6,7,8,9};
            par_nstep= ppns;
  			}
  			else{
          //use first observation rate
    				par_varying = {1,2,3,4,6,7,8,9};
            par_nstep= ppns;
  			}
  }

  //# give feedback to user on decisions just made
  printProposedParams(); //getchar();
  if(par_nstep> par_varying.size()) { printW "Hey! par_nstep needs << par_varying.size() !" endL; exit(-1); }


  //# initialize general variables for the MCMC cycle; description of variables in their
  //# definition in the beginning of this function

	nMCMCStepsAcceptedObs= 0;
	nMCMCStepsAccepted= 0;
	nMCMCStepsSaved= 0;
	recordSeriesThres= round(nMCMCSteps*(1-eta));
	y= dmatrix(0, NEQUATIONS, 0, nsteps); //solution of the solver
  countIntros= 0;

  //# if doing simulation mode, call function, run and terminate there
	if(simulationMode) doSimulationMode(odeType, sampleSims);

  //# this will be MCMC mode (fitting)
	printL "starting MCMC cycle" endL;
	for(int sim=0; sim< nMCMCSteps; sim++){

		if(sim==0){

      //# make particular decisions on first run
			findAndSetIntroStep(false); //#find model step for proposed time of introduction
			printProposedParams(); //# give feedback to user
			initPopulationStatus(false); //#initialize model population
			rkdumb(initialCondition, NEQUATIONS, timestart, time2run, nsteps); //#run model
			exportSimulation(y); //#export this simulation / run
			testSimulation(y,sim); //#proxy-test if solver solved correctly
			logLa= likelihood_gaussPrior(); //#get loglikelihood
      //# proxy-text if something is wrong with first loglikelihood and abort if so
			if(logLa>1 | logLa==0) {printW "First step had invalid log-likelihoods" endL; exit(-1);}
			nMCMCStepsAcceptedObs++;
			if(movAvgSDSeries) saveTimeseries(nMCMCStepsSaved, false); //#record the this simulation
			nMCMCStepsSaved++;
			saveParAndErrAP_acceptPar(sim, aProb, 1); //#save MCMC chain state
			printL "Sim:" << sim << " started Markov Chain." endL;
		}
		else{

			stepParameters(); // jump parameters / MCMCM chain state
			initPopulationStatus(false); //initialize model population
			rkdumb(initialCondition, NEQUATIONS, timestart, time2run, nsteps); //# run model
			testSimulation(y,sim); //proxy-test if solver solved correctly
			logLp= likelihood_gaussPrior(); //get loglikelihood
			LH_ap= MH_likelihood_ratio(logLp, logLa); //get loglikelihood ratio
			if(isnan(LH_ap)){ printW "likelihood is NaN!" endL; exit(-11); } //if ratio went wrong, abort

			randGuess= numGen->uniGen(); //flipt 'uniform' coin
			aProb= (LH_ap>1) ? 1 : LH_ap; //set acceptance probability
			if(randGuess <= aProb){
				//proposed MCMC state is accepted
				nMCMCStepsAccepted++;
        //save model solutions while building / showing user feedbank for decision
				printL "(A)" << nMCMCSteps-sim << ((sim>= recordSeriesThres) ? "(YS)" : "(NS)");
				if(sim>= recordSeriesThres){
					nMCMCStepsAcceptedObs++;
					if(movAvgSDSeries) saveTimeseries(nMCMCStepsSaved, true); //if measureing mean model behaviour, record this solution
					nMCMCStepsSaved++;
					cout << "# (aR)" << round(100*nMCMCStepsAcceptedObs/(nMCMCStepsSaved+1)) << "%";
				} else{
					if(movAvgSDSeries) saveTimeseries(nMCMCStepsSaved, false); //if measureing mean model behaviour, record this solution
					cout << "# (aR)" << round(100*nMCMCStepsAccepted/(sim+1)) << "%";
				}

        //# continue to build / show user some feedback on decision made
				cout << " (T0)" << cT0 << "[" <<introStep<< "] (K)" << cK << " (eipH)" << cEipH << " (infpH)" << cInfpH endL;
				printL " (eco2f)" << cEco2FactorV << " (ecof)" << cEcoFactorV << " (epif)" << cEipFactorV ;
				if(fRho) cout << " (rho)" << cRho;
				if(!noZeta) cout << " (zeta)" << cZeta;
				cout << endl;
				printL " (LK)" << (LH_ap) endL endL;
        //save parameters and acceptance variables
				saveParAndErrAP_acceptPar(sim, aProb, 1);
        //set currently accepted likelihood
				logLa= logLp;
			}
			else{
        //#proposed MCMC state is rejected
        //# build / show no feedback to the user; save previous simulation as accepted
				if(sim>= recordSeriesThres){
					if(movAvgSDSeries) savePrevTimeseries(nMCMCStepsSaved, true);
					nMCMCStepsSaved++;
					//cout << " (aR)" << round(100*nMCMCStepsAcceptedObs/(nMCMCStepsSaved+1)) << "%";
				} else{
					if(movAvgSDSeries) savePrevTimeseries(nMCMCStepsSaved, false);
				}
				//# save previous MCMC state as accepted
				savePrevParAndErrAP_roolBackPar(sim, aProb, 0);
			}
		}

	} //#MCMC ends here

  //#free model memory
	free_dmatrix(y, 0, NEQUATIONS, 0, nsteps);

  //#give user some feedback on number of events
	printL "accepted sims: " << nMCMCStepsAcceptedObs << " observed sims: " << nMCMCStepsSaved endL;
	printL "acceptance rate (observed period) = " << 100*nMCMCStepsAcceptedObs/float(nMCMCStepsSaved) << " %" endL;

	//## export accepted parameters (posteriors)
	fileNameBuilder.str(""); //clean filename
	fileNameBuilder << "output/"<<tag<<".accepted_parameters.data"; //define filename of output
	fout_all.open(fileNameBuilder.str().c_str()); //open the file
  //#export header of output table
	fout_all << "T0 T0step K ecoFactorV eipFactorV zeta aProb accept rho eipH infpH eco2FactorV" endL;
  if(fRho){
    for(int i=0; i< nMCMCSteps; i=i+1)
  		fout_all 	<< T0_accept[i] << " "
  					<< T0step_accept[i] << " "
  					<< K_accept[i] << " "
  					<< ecoFactorV_accept[i] << " "
  					<< eipFactorV_accept[i] << " "
  					<< zeta_accept[i] << " "
  					<< aProb_accept[i] << " "
  					<< accept_accept[i] << " "
  					<< zeta_accept[i] << " "
  					<< eipH_accept[i] << " "
  					<< infpH_accept[i] << " "
  					<< eco2FactorV_accept[i]
            endL;
  }else{
    for(int i=0; i< nMCMCSteps; i=i+1)
  		fout_all 	<< T0_accept[i] << " "
  					<< T0step_accept[i] << " "
  					<< K_accept[i] << " "
  					<< ecoFactorV_accept[i] << " "
  					<< eipFactorV_accept[i] << " "
  					<< zeta_accept[i] << " "
  					<< aProb_accept[i] << " "
  					<< accept_accept[i] << " "
  					<< rho_accept[i] << " "
  					<< eipH_accept[i] << " "
  					<< infpH_accept[i] << " "
  					<< eco2FactorV_accept[i]
            endL;
  }

  //#close output file
	fout_all.close();

  //if measuring mean model behaviour, export that mean behaviour
	if(movAvgSDSeries) exportSeries();

	printL "done all " endL;
}
