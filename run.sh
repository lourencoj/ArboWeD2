#!/bin/bash

#####################################################################

# Evaluate a floating point number expression in this bash script
float_scale=2 #resolution for function below
function float_eval()
{
    local stat=0
    local result=0.0
    if [[ $# -gt 0 ]]; then
        result=$(echo "scale=$float_scale; $*" | bc -q 2>/dev/null)
        stat=$?
        if [[ $stat -eq 0  &&  -z "$result" ]]; then stat=1; fi
    fi
    echo $result
    return $stat
}

##make sure all R files used are executable
cd code
chmod +x *.R
cd ..

##to measure running time
START=$(date +%s)

##compiling instruction (current flags may be specific to development)
g++ -O3 code/new_2factors_kftemp_main.cpp code/nrutil.c -o model -std=c++11 -ffast-math -freciprocal-math -fno-math-errno -funroll-loops

stepsize=0.01 #ODE step size
NH0=612000 #total human population size
odeType=1 #1= deterministic, 2= stochastic; defaults to 1 if simMode=0
tag="ArboWeD2_dev" #all output files will start with this tag

simMode=0 #simMode (1) will try to read a list of accepted parameters and run the model with those parameters N parameter times
tagInPar="---" #file name in output containing MCMC posterior solutions to use in simMod
sampleSims=0.00625 #prob for sampling MCMC steps in sim mode
maxSampleSims=10000000000 #sample with freq sampleSims or/and stop at N=maxSampleSims

#input files
fileinEpi="dev_epi.csv" #epi data
fileinTemp="dev_temp.csv" #temperature data
fileinPrec="dev_prechum.csv" #humidity and precipitation data

psi=0 #probability per time step of migration of 1 human infection
tpsi=99999999999999 #time point at which psi gets activated

seed=-55521852 #seed for random number generator
nMCMCSteps=2000 #number of desired MCMC states
tRho=999999 #time point in which the observation rate switches from first to second
noZeta=0 #0= estimate observation rate; 1= do not estimate observation rate
eta=0.33 #last proportion of MCMC steps to consider for fitting / posteriors (burnin)

##initial guesses for possible estimated parameters
##note here that 'rho' is used for the second observation rate, not for the non-linear factor (*)
##this is for historical reasons and may differ from published material on this code
sT0=-124	#initial estimate of: introduction day
sK=2650000 #initial estimate of: aquatic carrying capacity
sEcoFactorV=1.5  #initial estimate of: linear factor for temperature on adult vector mortality (alpha)
sEco2FactorV=1.5 #initial estimate of: non-linear factor for humidity and precipitation (rho) (*)
sEpiFactorV=1.5 #initial estimate of: linear factor for temperature on EIP (eta)
sEipH=6.0 #initial estimate of: incubation period of humans
sInfpH=6.0 #initial estimate of: infectious period of humans
sZeta=0.01 #initial estimate of: first observation rate
sRho=0.01 #initial estimate of: second observation rate (activated after tRho)

precaccdays=7 #number of days to smooth climate data
precaccsides=2 #1=moving-forward average; 2=moving average
smoothtemp=0 #0=do not smooth temperature; 1=smooth temperature

#change these to guarantee acceptance RATES 0.15-0.25
sampleGaussian=1 #0=use Cauchy sampling; 1=use Gaussian sampling
SAMPLE_factor_T=0.5 #sampling jump for time of introduction
SAMPLE_factor_K=1000.0 #sampling jump for carrying capacity
SAMPLE_factor_fV=0.002 #sampling jump for linear and non-linear factors
SAMPLE_EipInfH=0.002 #sampling jump for infectious and incubation periods of humans
SAMPLE_factor_zetaRho=0.002 #sampling jump for observation rate

#interpolates climate data to model step size
./code/new_interpolate_data_IT.R $stepsize $tag $fileinEpi $fileinTemp $fileinPrec $precaccdays $precaccsides $smoothtemp

#runs c code that simulates or fits data
./model $stepsize $sampleSims $nMCMCSteps $NH0 $eta $seed $tag $sT0 $sK $sEcoFactorV $odeType $sZeta $sRho $tRho $sEipH $noZeta $sInfpH $sEpiFactorV $sEco2FactorV $simMode $sampleGaussian $SAMPLE_factor_T $SAMPLE_factor_K $SAMPLE_factor_fV $SAMPLE_EipInfH $SAMPLE_factor_zetaRho $tagInPar $psi $tpsi $maxSampleSims

#measures running time
echo "$(( $(date +%s) - $START )) seconds of execution."

##makes simple output png figures
if [ $simMode -eq 0 ]
	then
	./code/plot5.R $eta $tag $fileinEpi $tRho $NH0
fi
