#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>

//# the code in this file, as well as nrutil.h and nrutil.c is based on the code from:
//# Numerical recipes in C (2nd ed.): the art of scientific computing
//# Cambridge University Press New York, NY, USA ©1992
//# ISBN:0-521-43108-5

#include "nrutil.h"
#include "new_printing.h"

    void rk4(double [], double [], int , double , double , double [], double, int);
    void inline modelDeterministic(double x, double y[], double dydx[], int istep);
    void inline modelStochastic(double x, double y[], double dydx[], int istep);
    void inline auxiliarmodel(double y[], int k);

void rk4(double y[], double dydx[], int n, double x, double h, double yout[], int istep)
// Given values for the variables y[1..n] and their derivatives dydx[1..n] known at x, use the
// fourth-order Runge-Kutta method to advance the solution over an interval h and return the
// incremented variables as yout[1..n], which need not be a distinct array from y. The user
// supplies the routine model(x,y,dydx) , which returns derivatives dydx at x.
{
    int i;
    double xh,hh,h6;

    double dym[n+1];
    double dyt[n+1];
    double yt[n+1];
    hh=h*0.5;
    h6=h/6.0;
    xh=x+hh;

    for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];

    modelDeterministic(xh,yt,dyt,istep);

    for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
    modelDeterministic(xh,yt,dym,istep);

    for (i=1;i<=n;i++) {
      yt[i]=y[i]+h*dym[i];
      dym[i] += dyt[i];
    }

    modelDeterministic(x+h,yt,dyt,istep);

    for (i=1;i<=n;i++) yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);

}


double **y; //For communication back to main.

//#runs the deterministic version of the model
void rkdumb(double vstart[], int nvar, double x1, double x2, int nstep)
// Starting from initial values vstart[1..nvar] known at x1 use fourth-order Runge-Kutta
// to advance nstep equal increments to x2. The user-supplied routine model(x,v,dvdx)
// evaluates derivatives. Results are stored in the global variables y[1..nvar][1..nstep+1]
{
    int k, i;
    double x,h;

    double v[nvar];
    double vout[nvar];
    double dv[nvar];

    double importCases;

    // time variables
    x=x1;
    h=stepsize;


    //save initial point
    // Load starting values.
    for (k=0; k<nvar; k++) {
       v[k]= vstart[k];
       y[k][0]= v[k];
    }

    for (k=0; k<nstep; k++) { //Take nstep steps.

        //###################################################
        //address introduction
        importCases= (k == introStep) ? INTRO_CASESH : 0;
        v[1]= v[1];
        v[2]= v[2];
        v[3]= v[3];
        v[6]= v[6]- importCases;
        v[7]= v[7]+ round(importCases/2);
        v[8]= v[8]+ round(importCases/2);
        //###################################################

		modelDeterministic(x,v,dv,k);
		rk4(v,dv,nvar,x,h,vout,k);

        //###################################################

		if ((double)(x+h) == x)
            nrerror((char*)"Step size too small in routine rkdumb");

        x += h;

		for (i=0;i<nvar;i++) { //Store intermediate steps.
            if(i==0){
                v[i]= x;
                y[i][k+1]= x;
            }else{
                v[i]= vout[i];
                y[i][k+1]= v[i];
            }
		}



    }

}

void rk4Fake(double y[], double dydx[], int n, double x, double h, double yout[], int istep) {
    for (int i=1;i<=n;i++) yout[i]=y[i]+ dydx[i];
}

//#runs the stochastic version of the model; looks and behaves like a solver but does not use
//#a solver
void stochBindumb(double vstart[], int nvar, double x1, double x2, int nstep){

    double v[nvar];
    double vout[nvar];
    double dv[nvar];
    double importCases;
    double x;
    double h;
    int k, i;

    // time variables
    x= x1;
    h= stepsize;

    //save initial point
    // Load starting values.
    for (k=0; k<nvar; k++) {
       v[k]=vstart[k];
       y[k][0]=v[k];
    }

    for (k=0; k<nstep; k++) { //Take nstep steps.

        //###################################################
        //address introduction
        importCases= (k == introStep) ? INTRO_CASESH : 0;
        v[1]= v[1];
        v[3]= v[3];
        v[6]= v[6]- importCases;
        v[7]= v[7]+ round(importCases/2);
        v[8]= v[8]+ round(importCases/2);
        //###################################################

        modelStochastic(x,v,dv,k);
        rk4Fake(v,dv,nvar,x,h,vout,k);

        //###################################################

        x += h;

        for (i=0;i<nvar;i++) { //Store intermediate steps.
            if(i==0){
                v[i]= x;
                y[i][k+1]= x;
            }else{
                v[i]= vout[i];
                y[i][k+1]= v[i];
            }
        }

    }

}
