/***************************************************************************
@ F. Rousset 2005-2007

francois.rousset@umontpellier.fr

This file is part of Genepop'007
This software is a computer program whose purpose is to perform statistical analyses.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

 ***************************************************************************/

#include <iostream>
#include <sstream>
#include <limits>
#include "GenepopS.h"
#include "F_est.h"
#include "bootstrap.h"
#include "settings.h"
#include "proba.h"
#include "myutils.h"
#include "tools.h"
#ifdef COMPATIBILITYRCPP
#include <Rcpp.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#endif

using namespace std;

// This will be used across bootstrap methods
const double epsn_value=0.00002;

/** in no way specific to bootstrap**/
vector<double> bisection_search(double (*func)(double d),double x1,double x2,bool verbose) {
//R: function (x1, x2, f, tol = 1e-07, niter = 25, f.extra = NA, upcross.level = 0)
// rtbis, Numerical recipes
//Using bisection, find the root of a function func known to lie between x1 and x2. The root,
//returned as rtbis, will be re...ned until its accuracy is ...xacc.
    vector<double>v(0); // [0] : indicateur probleme eventuel; [1]: actual P value
    const double xacc=numeric_limits<double>::epsilon()*(abs(x1)+abs(x2))/2; //suggestion Num. Rec.
    const int maxit=-2*int(log(numeric_limits<double>::epsilon())/log(2.));
    int j;
    double dx,f,fmid,xmid,rtb;
    f=(*func)(x1);
    fmid=(*func)(x2);
    if (f*fmid >= 0.0) {
        if (verbose) {
                noR_cout<<"(!) From CKrig::bisection_search() : Root must be bracketed for bisection. "<<endl;
                noR_cout<<"   x1, f(x1), x2, f(x2) were "<<x1<<" "<<f<<" "<<x2<<" "<<fmid<<endl;
        }
        //if (cinGetOnError) cin.get();
        v.push_back(-1);
        return(v); // length 1
    }
    rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2); //Orient the search so that f>0
/*std::cout<<x1<<" ";
std::cout<<x2<<" ";
std::cout<<dx<<" ";
std::cout<<rtb;getchar();*/
    for (j=1;j<=maxit;j++) { //lies at x+dx.
        fmid=(*func)(xmid=rtb+(dx *= 0.5)); //Bisection loop.
        if (fmid <= 0.0) rtb=xmid;
        if (fabs(dx) < xacc || fmid == 0.0) {
            v.push_back(0);
            v.push_back(rtb);
            return(v);
        }
    }
    if (f*fmid >= 0.0) { //mainly means !=0...
        if (verbose) {
                noR_cout<<"(!) From CKrig::bisection_search() : Too many bisections. "<<endl;
        }
        //if (cinGetOnError) cin.get();
        v.push_back(-1);
        v.push_back(numeric_limits<double>::quiet_NaN());  // visible problem; otherwise perhps xmid ?
        return(v);
    }
// should never reach:
v.push_back(numeric_limits<double>::signaling_NaN());
return(v);
}

static CABCbootstrap* ABCptr;

CABCbootstrap::CABCbootstrap(size_t units) {
  nb_units=units;
  z=std::numeric_limits<double>::quiet_NaN();
}

int CABCbootstrap::bootstrapOverLociABC(// Function returning the statistic of interest
                                        double (*estimatingFnPtr)(vector<double> d), 
                                        // Simple textual info to add to progress bars
                                        string legend,
                                        // Where to write the bootstrap results 
                                        string bootOutfile,
                                        // Clear screen before running computations
                                        bool clear_screen) {
  ofstream bootOut;
  estimFnPtr=estimatingFnPtr;
  testLegend=legend; // copy to namespace variable used by other function
  // nb_units is a number of relevant loci (in particular with the relevant ploidy)
  // it is the function called through estimatingFnPtr that must map the 
  // weights to the chosen subset of loci using a set of indices for the 
  // mapping (see varForBootstrapGenepop::idxPloid in wrapper() )
  // Alex: here nb_units stems from the constructor above
  // 
  // Useful resources: 
  //   boot:::abc.ci()
  // 
  if (nb_units == 0) { 
#ifdef COMPATIBILITYRCPP
    Rcpp::Rcerr<<"(!) 0 replicates to bootstrap over";
#else
    cerr<<"(!) 0 replicates to bootstrap over";
#endif
    return(0);
  }
  
  delta.resize(nb_units);
  if (clear_screen) effacer_ecran();
  vector<double> locABCweight(nb_units);
  vector<double> dt(nb_units);
  vector<double> ddt(nb_units);
  double epsn,tm,tp,cq,bhat,curv,qnorm_arg_4_12;//ABC bootstrap variables
  double sigmahat=0.0;
  seuil_inf=(1.0-widthCI)/2.0;
  seuil_sup=1.0-seuil_inf;
  
  // point estimate for original data (ABCweight constant)
  for (size_t loc=0; loc<nb_units; loc++) {
    locABCweight[loc]=1.0/nb_units;
  }
  t0 = (*estimatingFnPtr)(locABCweight);  //first call to isoldeproc sets _first_of_repl=false (and if it was =true, also clears screen !)
  
  // does not delete the .MIG file here
  if (std::isnan(t0)) {
    tinf=numeric_limits<double>::quiet_NaN();
    tsup=numeric_limits<double>::quiet_NaN();
    noR_cout<<"\n No valid estimate for original data. Skipping further bootstrap computation.\n";
    if (pauseGP) {
      noR_cout<<"(Return) to continue";
      cin.get();
    }
    genepop_exit(-1, "No valid estimate for original data. Skipping further bootstrap computation.");
  }
  
  // Not enough loci for bootstrap
  if (nb_units<3) {
    bootOut.open(bootOutfile.c_str(),ios::app);
    if ( ! bootOut.is_open()) {
#ifdef COMPATIBILITYRCPP
      // Rcpp::Rcerr<<"(!) From bootstrapOverLociABC(): error while opening file "<<bootOutfile<<endl;
#else
      cerr<<"(!) From bootstrapOverLociABC(): error while opening file "<<bootOutfile<<endl;
      if (cinGetOnError) cin.get();
#endif
      
      genepop_exit(-1, "(!) From bootstrapOverLociABC(): error while opening file ");
    }
    noR_cout<<"\nNot enough loci to bootstrap over\n";
    bootOut<<"\nNot enough loci to bootstrap over\n";
    bootOut.close();
    if (pauseGP) {
      noR_cout<<"(Return) to continue"<<endl; getchar();
    }
    genepop_exit(-1, "Not enough loci to bootstrap over");
  }
  
  /* else: everything is in order, proceed with bootstrap */
  
  // Display information on terminal
  int consy=wherey();
  _gotoxy(0,consy+2);
  noR_cout<<" Computing [ "<<seuil_inf<<" - "<<seuil_sup<<" ] confidence interval"<<legend<<":";
  //pour borne inf
#ifdef COMPATIBILITYRCPP
  Rcpp::Rcerr<<"Computing confidence interval"<<legend<<":"<<endl;
#endif
  
  // Compute the upper and lower values for derivatives. Note that 
  // get_nearby_epsn handles the progress report
  vector<double> tminus = get_nearby_epsn(-1, nb_units, 
                                          // Starting at 1, 2*nb_units tasks 
                                          // on line +3
                                          1, 2 * nb_units, consy+3, 
                                          *estimatingFnPtr); 
  vector<double> tplus = get_nearby_epsn(1, nb_units, 
                                         // Starting at nb_units, 2*nb_units
                                         // tasks on line +3
                                         nb_units, 2 * nb_units, consy+3, 
                                         *estimatingFnPtr); 
  
  // Compute the derivatives of that function
  for (size_t loc=0; loc<nb_units; loc++) {
    // Approximation of derivative
    dt[loc]  = (tplus[loc] - tminus[loc]) / (2 * epsn_value);
    // Approximation of second derivative
    ddt[loc] = (tplus[loc] + tminus[loc] - 2*t0) / pow(epsn_value, 2);
    //cout<<tpinput[loc]<<" "<<dt[loc]<<" "<<ddt[loc]<<endl;
  }
//getchar();
  
  // sigmahat here defined to give root denominator of formula 4.9 for acceleration, up to 'nb_units' factor
  for (size_t loc=0; loc<nb_units; loc++) { 
    sigmahat += pow(dt[loc], 2);
  }
  sigmahat = sqrt(sigmahat) / nb_units;
  
  // Handle case where S.D. of bootstrap distribution is zero
  if ( sigmahat == 0 ) {
    bootOut.open(bootOutfile.c_str(),ios::app);
    if ( ! bootOut.is_open() ) {
#ifdef COMPATIBILITYRCPP
      // Rcpp::Rcerr<<"(!) From bootstrapOverLociABC(): error while opening file "<<bootOutfile<<endl;
#else
      cerr<<"(!) From bootstrapOverLociABC(): error while opening file "<<bootOutfile<<endl;
      if (cinGetOnError) { 
        cin.get();
      }
#endif
      genepop_exit(-1, "(!) From bootstrapOverLociABC(): error while opening file ");
    }
    noR_cout<<"\nNot enough data for constructing confidence interval\n";
    bootOut<<"\nNot enough data for constructing confidence interval\n";
    bootOut.close();
    if (pauseGP) {
      noR_cout<<"(Return) to continue"<<endl; getchar();
    }
    genepop_exit(-1, "Not enough data for constructing confidence interval");
  }
  
  // Estimate acceleration (equation 4.9)
  ahat = 0;
	for (size_t loc=0; loc<nb_units; loc++) { 
    ahat += pow(dt[loc], 3);
  }
	ahat = ahat / ( 6 * pow(double(nb_units), 3) * pow(sigmahat, 3) );
  
  // dhat in boot:::abc.ci()
	for (size_t loc=0; loc<nb_units; loc++) { 
    delta[loc] = dt[loc] / ( pow(double(nb_units),2) * sigmahat );
  }
  
  // 
  epsn = -epsn_value;
  for (size_t loc=0; loc<nb_units; loc++) { 
    locABCweight[loc] = delta[loc] * epsn + (1.0/nb_units);
  }
  tp = (*estimatingFnPtr)(locABCweight);
  
  for (size_t loc=0;loc<nb_units;loc++) { 
    locABCweight[loc] = - delta[loc] * epsn + (1.0/nb_units);
  }
  tm = (*estimatingFnPtr)(locABCweight);
  
  if (use_console_xy) {
    _gotoxy(0,consy+3);
    noR_cout<<" Computing confidence interval... about      % done              ";
    _gotoxy(41,consy+3);
    noR_cout << int(100*(2*double(nb_units)+2)/(2*double(nb_units)+5));
  }
  
  /* cq is chat in boot::abc.ci, where it is computed from  
   t(w.orig + eps * dhat)+t(w.orig + eps * dhat)-2 t0, where dhat corresponds to our vector delta 
   (and to k in DavisonH97 eq. 5.49). Thus 
   t(w.orig + eps * dhat) <-> tp, t(w.orig - eps * dhat) <-> tm, cq <-> chat */
  cq = (tp + tm - 2*t0 ) / ( 2 * sigmahat * pow(epsn, 2) );
  
  /* bhat is bhat in boot::abc.ci, where it is computed as below, as bhat= sum(L2)/(2 * n * n) where L2[i] <- (t1 - 2 * t0 + t2)/eps^2
   [with t1 <-> tplus, t2 <-> tminus, L <-> dt, L2 <-> ddt, bhat <->bhat] 
   Eq. 4.11 expresses it in terms of the eigensystem of a matrix of second derivatives, but DavisonH97 eq. 5.49 
   expresses it as the trace of this matrix, which corresponds to the present computation. */
  bhat = 0.0;
  for (size_t loc=0;loc<nb_units;loc++) { 
    bhat += ddt[loc];
  }
  bhat /= ( 2 * pow(double(nb_units), 2) );
  
  // Eq. 4.12
  curv = cq - ( bhat / sigmahat );
  // ndtr=pnorm
  qnorm_arg_4_12 = 2 * ndtr(ahat) * ndtr(curv);
  
  // ndtri=qnorm... ndtri(qnorm_arg_4_12) is a valid bound only if 'qnorm_arg_4_12' is a probability, which may fail to be the case.
  // In that case we need to bail. 
  if ( (qnorm_arg_4_12>=1) || (qnorm_arg_4_12<=0.0) ) { //z=INFINITY;
    bootOut.open(bootOutfile.c_str(),ios::app);
    if ( ! bootOut.is_open() ) {
#ifdef COMPATIBILITYRCPP
      // Rcpp::Rcerr<<"(!) From bootstrapOverLociABC(): error while opening file "<<bootOutfile<<endl;
#else
      cerr<<"(!) From bootstrapOverLociABC(): error while opening file "<<bootOutfile<<endl;
      if (cinGetOnError) { 
        cin.get();
      }
#endif
      genepop_exit(-1, "(!) From bootstrapOverLociABC(): error while opening file ");
    }
    
    noR_cout<<"Putative confidence interval"<<legend<<" has infinite range.\n Bootstrap computation aborted.";
    bootOut<<"Putative confidence interval"<<legend<<" has infinite range.\n Bootstrap computation aborted.";
    bootOut.close();
    if (pauseGP) {
      noR_cout<<"\n(Return) to continue"<<endl; getchar();
    }
    return 1; // exit() !
    
  // Correct bounds, actually retrieve estimate of bias
  } else { 
    z = ndtri(qnorm_arg_4_12);
  }
  
  // We now have estimates for bias and ahat, we can correct the bounds of 
  // the bootstrap CI. 
  
  // lambda in eq. 4.13, cf also lalpha in boot:::abc.ci()
  lambda_4_13 = ( z + ndtri(seuil_inf) ) /  
                 ( pow(1 - ahat * (z + ndtri(seuil_inf)), 2) ) ;
  for (size_t loc=0; loc<nb_units; loc++) { //pour borne inf
    locABCweight[loc] = delta[loc] * lambda_4_13 + 1.0/nb_units;
  }
  tinf = (*estimatingFnPtr)(locABCweight);
  
  if (use_console_xy) {
    _gotoxy(0,consy+3);
    noR_cout<<" Computing confidence interval... about      % done             ";
    _gotoxy(41,consy+3);
    noR_cout<< int(100*(2*double(nb_units)+4)/(2*double(nb_units)+5));
  }
  
  // Correct upper bound
  lambda_4_13 = ( z + ndtri(seuil_sup) ) / 
                 ( pow(1 - ahat * (z+ndtri(seuil_sup) ), 2) );
	for (size_t loc=0; loc<nb_units; loc++) { //pour borne sup
    locABCweight[loc] = delta[loc] * lambda_4_13 + 1.0/nb_units;
  }
	tsup = (*estimatingFnPtr)(locABCweight);
  
  // Invert tsup and tinf if bounds are inverted (why does this happen ?)
	if ( tsup < tinf ) { 
    double stat = tsup; 
    tsup=tinf; 
    tinf=stat;
  }
  
  _gotoxy(0,consy+3);
  noR_cout<<" Computing confidence interval... 100 % done                      ";
#ifdef COMPATIBILITYRCPP
  Rcpp::Rcerr << "\n ABC bootstrap results" << ":" << endl;
  Rcpp::Rcerr << "\n Point estimate and "<< 100 * widthCI << 
    "% confidence interval" << legend << ":" << endl;
  Rcpp::Rcerr << t0 << " [ " << tinf << " , " << tsup << " ]\n";
#else 
  noR_cout<<"\n ABC bootstrap results"<<":";
  noR_cout<<"\n Point estimate and "<< 100 * widthCI <<"% confidence interval" 
    << legend << ":" << endl;
  noR_cout << t0 << " [ " << tinf << " , " << tsup << " ]\n";
#endif
  
  bootOut.open(bootOutfile.c_str(), ios::app);
  if ( ! bootOut.is_open() ) {
#ifdef COMPATIBILITYRCPP
    // Rcpp::Rcerr<<"(!) From bootstrapOverLociABC(): error while opening file "<<bootOutfile<<endl;
#else
    cerr<<"(!) From bootstrapOverLociABC(): error while opening file "<<bootOutfile<<endl;
    if (cinGetOnError) { 
      cin.get();
    }
#endif
    
    genepop_exit(-1, "(!) From bootstrapOverLociABC(): error while opening file ");
  }
  
  if ( ! perf ) {	 //ecriture des resultats
    bootOut<<"\n ABC bootstrap results"<<legend<<":";
    bootOut<<"\n Point estimate and "<<100*widthCI<<"% confidence interval:\n"<<t0<<" [ "<<tinf<<" , "<<tsup<<" ]\n";
  }
  
  if ( std::strcmp(legend.c_str(), " for SLOPE")==0 && ! std::isnan(testPointslope)) {
    //FR->FR comme on n'est pas surde travailler sur des slopes ici, il faut un mecanisme genre migraine sur les testPoint
    testPointPvalue=Pvalue(testPointslope,true,false);
    if ( std::isnan(testPointPvalue) ) { 
      Pvalue(testPointslope,true,true); // verbose version
    }
#ifdef COMPATIBILITYRCPP
    Rcpp::Rcerr << "\nUnidirectional P-value for [ slope=" << testPointslope << 
      " ] : " << testPointPvalue << endl;
#else 
    noR_cout << "\nUnidirectional P-value for [ slope=" << testPointslope << 
      " ] : " << testPointPvalue << endl;
#endif
    bootOut<<"\nUnidirectional P-value for [ slope="<<testPointslope<<" ] : "<<testPointPvalue<<endl;
  }
    bootOut.close();
    ///// bug corrected 2015/12/8 : added lines:
  for (size_t loc=0;loc<nb_units;loc++) { 
    locABCweight[loc]=1.0/nb_units;
  }
  t0 = (*estimatingFnPtr)(locABCweight);  //need to reset any variable modified by estimatingFnPtr (slope() -> datamatrix::data used later by mantel test)
  
  return 0;
}

double norm_inter_sorted(vector<double>t_sorted, double alpha) {
  // based on the logic of boot:::norm.inter
  // assuming all bootstrap estimates are valid and that they are already sorted.
  size_t B = t_sorted.size();
  size_t B1 = B+1;
  double rk = alpha*B1;
  if (rk <= 1) {
    //   warning("extreme order statistics used as endpoints")
    return(t_sorted[0]);
  } else if (rk >= B) {
    //   warning("extreme order statistics used as endpoints")
    return(t_sorted[B-1]);
  } else {
    double k = trunc(rk);
    double temp1 = ndtri(alpha);
    double temp2 = ndtri(k/B1);
    double temp3 = ndtri((k+1)/B1);
    double tk = t_sorted[k-1];
    double tk1 = t_sorted[k];
    double resu = tk+ (tk1-tk)*(temp1-temp2)/(temp3-temp2);
    return(resu);
  }
  
  //     cbind(round(rk, 2), out)
      
}

int CABCbootstrap::bootstrapOverLociBCa(// Function returning the statistic of interest
                                        double (*estimatingFnPtr)(vector<double> d), 
                                        // Simple textual info to add to progress bars
                                        string legend,
                                        // Where to write the bootstrap results 
                                        string bootOutfile,
                                        // Whether to estimate acceleration 
                                        bool estim_accel, 
                                        // Clear screen before running computations
                                        bool clear_screen) {  
  
  ofstream bootOut;
  estimFnPtr=estimatingFnPtr;
  testLegend=legend; // copy to namespace variable used by other function
  double epsn_value=0.00002;
  string boot_type = estim_accel ? "BCa" : "BC"; 
  
  // nb_units is a number of relevant loci (in particular with the relevant
  // ploidy) it is the function called through estimatingFnPtr that must map the 
  // weights to the chosen subset of loci using a set of indices for the 
  // mapping (see varForBootstrapGenepop::idxPloid in wrapper() )
  if (nb_units == 0) { 
#ifdef COMPATIBILITYRCPP
    Rcpp::Rcerr<<"(!) 0 replicates to bootstrap over";
#else
    cerr<<"(!) 0 replicates to bootstrap over";
#endif
    return(0);
  }
  
  delta.resize(nb_units);
  if (clear_screen) effacer_ecran();
  vector<double> empinf(nb_units);
  seuil_inf=(1.0-widthCI)/2.0;
  seuil_sup=1.0-seuil_inf;
  
  // point estimate for original data (weights constant)
  vector<double> weights(nb_units, 1.0/nb_units); 
  t0 = (*estimatingFnPtr)(weights); 
  
  // does not delete the .MIG file here
  if (std::isnan(t0)) {
    tinf=numeric_limits<double>::quiet_NaN();
    tsup=numeric_limits<double>::quiet_NaN();
    noR_cout<<"\n No valid estimate for original data. Skipping further bootstrap computation.\n";
    if (pauseGP) {
      noR_cout<<"(Return) to continue";
      cin.get();
    }
    genepop_exit(-1, "No valid estimate for original data. Skipping further bootstrap computation.");
  }
  
  // Not enough loci for bootstrap
  if ( nb_units < 3 ) {
    bootOut.open(bootOutfile.c_str(),ios::app);
    if ( ! bootOut.is_open()) {
#ifdef COMPATIBILITYRCPP
      // Rcpp::Rcerr<<"(!) From bootstrapOverLociBCa(): error while opening file "<<bootOutfile<<endl;
#else
      cerr<<"(!) From bootstrapOverLociBCa(): error while opening file "<<bootOutfile<<endl;
      if (cinGetOnError) cin.get();
#endif
      
      genepop_exit(-1, "(!) From bootstrapOverLociBCa(): error while opening file ");
    }
    noR_cout<<"\nNot enough loci to bootstrap over\n";
    bootOut<<"\nNot enough loci to bootstrap over\n";
    bootOut.close();
    if (pauseGP) {
      noR_cout<<"(Return) to continue"<<endl; getchar();
    }
    genepop_exit(-1, "Not enough loci to bootstrap over");
  }
  
  /* else: everything is in order, proceed with bootstrap */
  
  // Display information on terminal
  int consy=wherey();
  _gotoxy(0,consy+2);
  noR_cout <<" Computing [ " << seuil_inf << " - " << seuil_sup <<
    " ] confidence interval"<<legend<<":";
  //pour borne inf
#ifdef COMPATIBILITYRCPP
  Rcpp::Rcerr << "Computing confidence interval" << legend << ":" << endl;
#endif
  
  // Start computations
  
  // We get the values of the simulation distribution, along with the number of
  // times the value is below the observed one.
#ifdef COMPATIBILITYRCPP
//   Rcpp::Rcerr << "Computing bootstrap distribution..." << endl;
  Progress progbar(nboot, true); 
#else 
  _gotoxy(0, consy+2);
  size_t progress = 1; 
#endif
  
  // Initialize random number generator 
  std::default_random_engine generator;
  std::uniform_int_distribution<size_t> distribution(0, nb_units-1);
  
  size_t nboot_below_obs = 0; 
  vector<double> simdistr(nboot); 
  double incr=1.0 / nb_units;
  for ( size_t n=0; n<nboot; n++ ) { 
    // Initialize weights to zero 
    std::fill(weights.begin(), weights.end(), 0); 
    for ( size_t bootn=0; bootn<nb_units; bootn++) { 
      // Get a random position in the vector. 
      size_t position = distribution(generator);
      weights[position] += incr; 
    }
    double current_value = (*estimatingFnPtr)(weights); 
    simdistr[n] = current_value; 
    nboot_below_obs += current_value < t0 ? 1 : 0; 
#ifdef COMPATIBILITYRCPP
    progbar.increment(); 
#else
    if (use_console_xy) {    
      _gotoxy(0, consy+3);
      noR_cout<<" Computing bootstrap distribution...       % done             ";
      _gotoxy(40, consy+3);
      noR_cout<< int(100 * (double)(progress) / ( nboot ) ) << flush;
    }    
progress++; 
#endif
  }
  
#ifdef COMPATIBILITYRCPP
  // We destroy the progress bar, so we are able to open a new one later 
  progbar.cleanup(); 
#endif
  
  // We compute the bias (eq. 2.8)
  z = ndtri( (double)(nboot_below_obs) / nboot ); 
  
  // Now we need to compute the acceleration. We first compute the empirical 
  // influence function using infinitesimal jacknife (sensu boot:::inf.jack())), 
  // and plug it in the estimate of a (eq. 6.6)
  double a = 0; 
  if ( estim_accel ) { // i.e., BCa, not BC
    
    // Compute t values with extra weight on observation.
    vector<double> tplus = get_nearby_epsn(1, nb_units, 
                                           // Starting at 1, nb_units tasks, 
                                           // on line consy+3
                                           1, nb_units, consy+3, 
                                           *estimatingFnPtr); 
    
    // Compute the actual influence function
    for (size_t loc=0; loc<nb_units; loc++) {
      empinf[loc]  = ( tplus[loc] - t0 ) / epsn_value; // correct scaling as explained in get_nearby_epsn()
    }
    
    // Estimate acceleration
    double uicube = 0; 
    double uisq   = 0; 
    for (size_t loc=0; loc<nb_units; loc++) { 
      uicube += pow(empinf[loc], 3); 
      uisq   += pow(empinf[loc], 2); 
    }
    a = uicube / (6. *pow(uisq, 1.5) ); 
    
    // a non-zero 'a' has no effect on CI as long as it does not move qhigh or qlow to another 1/B quantile step
    // => often no visible effect for small number B of bootstrap replicates.  
  } 
  
  // Okay, now we have z and a, we can use formula 2.3 to get quantiles 
  double zhigh = ndtri(seuil_sup); 
  double zlow  = ndtri(seuil_inf); 
  double qhigh = ndtr( z + ( z + zhigh ) / ( 1. - a * (z + zhigh)) ); 
  double qlow  = ndtr( z + ( z + zlow  ) / ( 1. - a * (z + zlow )) ); 
  
  std::sort(simdistr.begin(), simdistr.end());
  // BCa CI bounds:
  tinf = norm_inter_sorted(simdistr, qlow);
  tsup = norm_inter_sorted(simdistr, qhigh);

  // We only search for a P-value when dealing with slopes. Note that strcmp 
  // returns zero when true. 
  bool compute_pval = (std::strcmp(legend.c_str(), " for SLOPE") == 0);
  
  // next: not most elegant code, but there are other priorities...

  double pval = 0; 
  if ( compute_pval ) { 
    bool found_testpoint = false; 
    
    // If the test point is above or below the simulated distribution, then the 
    // P-value is 0 or 1. 
    if ( simdistr[0] > testPointslope ) { 
      pval = 0; 
      found_testpoint = true; 
    }
    if ( testPointslope > simdistr[nboot-1] ) { 
      pval = 1; 
      found_testpoint = true; 
    }
    
    // In other cases, we do a linear search to find the corresponding rank in 
    // the bootstrap distribution
    size_t i=1; 
    while ( (! found_testpoint) && i < nboot ) { 
      // Check if we are above the testPoint value. If that's true, 
      // then the observed value is between sim[i-1] and sim[i], so the rank 
      // of the observed value in the bootstrap distribution (merged with the 
      // observation) is i-1 + 1 = i
      if ( simdistr[i] > testPointslope ) { 
        pval = ( (double) i / (nboot + 1) ); 
        found_testpoint = true; 
      }
      i++; 
    }
  } // End of if ( compute_pval )
  
#ifdef COMPATIBILITYRCPP
  Rcpp::Rcerr << "Finished " << boot_type << " bootstrap" << endl; 
#else 
  _gotoxy(0,consy+3);
  noR_cout<<" Computing confidence interval... 100 % done                      ";
#endif
  
#ifdef COMPATIBILITYRCPP
  Rcpp::Rcerr << "\n " << boot_type << " bootstrap results" << ":" << endl;
  Rcpp::Rcerr << "\n Point estimate and "<< 100 * widthCI << 
    "% confidence interval" << legend << ":" << endl;
  Rcpp::Rcerr << t0 << " [ " << tinf << " , " << tsup << " ]\n";
#else 
  noR_cout << "\n " << boot_type << " bootstrap results" << ":";
  noR_cout << "\n Point estimate and " << 100 * widthCI << 
    "% confidence interval" << legend << ":" << endl;
  noR_cout << t0 << " [ " << tinf << " , " << tsup << " ]\n";
#endif
  
  bootOut.open(bootOutfile.c_str(), ios::app);
  if ( ! bootOut.is_open() ) {
#ifdef COMPATIBILITYRCPP
    // Rcpp::Rcerr<<"(!) From bootstrapOverLociBCa(): error while opening file "<<bootOutfile<<endl;
#else
    cerr<<"(!) From bootstrapOverLociBCa(): error while opening file "<<bootOutfile<<endl;
    if (cinGetOnError) { 
      cin.get();
    }
#endif
    
    genepop_exit(-1, "(!) From bootstrapOverLociBCa(): error while opening file ");
  }
  
  if ( ! perf ) {	 //ecriture des resultats
    bootOut << "\n " << boot_type << " bootstrap results" << legend << ":";
    bootOut << "\n Point estimate and " << 100*widthCI << 
      "% confidence interval:\n" << t0 << " [ " << tinf << " , " << tsup << 
      " ]\n";
  }
  
  if ( std::strcmp(legend.c_str(), " for SLOPE")==0 && ! std::isnan(testPointslope)) {
    
    // We need to format the P-value so when it is zero or one so we report 
    // something sensible. 
    std::ostringstream pvalstr;
    pvalstr.clear();
    // Note that this is floating-point comparison
    if ( pval == 0 ) { 
      pvalstr << "less than " << 1.0 / nboot; 
    } else if ( pval == 1 ) { 
      pvalstr << "greater than " << (double)(nboot - 1) / nboot; 
    } else { 
      pvalstr << pval; 
    }
    
#ifdef COMPATIBILITYRCPP
    Rcpp::Rcerr << endl << "Unidirectional P-value for [ slope=" << 
      testPointslope << " ] : " << pvalstr.str() << endl;
#else
    noR_cout << "\nUnidirectional P-value for [ slope=" << testPointslope << 
      " ] : " << pvalstr.str() << endl;
#endif
    bootOut << "\nUnidirectional P-value for [ slope=" << testPointslope << 
      " ] : " << pvalstr.str() << endl;
  }
  bootOut.close();
  
  ///// bug corrected 2015/12/8 : added lines:
  std::fill(weights.begin(), weights.end(), 1.0/nb_units); 
  t0 = (*estimatingFnPtr)(weights);  //need to reset any variabes modified by estimatingFnPtr (slope() -> datamatrix::data used later by mantel test)
  
  // Return
  return 0;
}


// Generic function to call the right bootstrap method above
int CABCbootstrap::bootstrapOverLoci(int method, 
                                     // Function returning the statistic of interest
                                     double (*estimatingFnPtr)(vector<double> d), 
                                     // Simple textual info to add to progress bars
                                     string legend,
                                     // Where to write the bootstrap results 
                                     string bootOutfile,
                                     // Clear screen before running computations
                                     bool clear_screen) {
  
// #ifdef COMPATIBILITYRCPP
//   Rcpp::Rcerr << "Using bootstrap method " << method << " (default: " 
//     << bootmethod << ")" << endl; 
// #endif
  
  if ( method == BOOT_METHOD_ABC ) { 
    CABCbootstrap::bootstrapOverLociABC(estimatingFnPtr, 
                                        legend, 
                                        bootOutfile, 
                                        clear_screen); 
    
  } else if ( method == BOOT_METHOD_BCA ) { 
    CABCbootstrap::bootstrapOverLociBCa(estimatingFnPtr, 
                                        legend, 
                                        bootOutfile, 
                                        true, // estimate alpha
                                        clear_screen); 
    
  } else if ( method == BOOT_METHOD_BC ) { 
    CABCbootstrap::bootstrapOverLociBCa(estimatingFnPtr, 
                                        legend, 
                                        bootOutfile, 
                                        false, // do not estimate alpha
                                        clear_screen); 
  
  } else { 
#ifdef COMPATIBILITYRCPP
    Rcpp::Rcerr << "Unknown bootstrap method (available methods are ABC (0) or BCa (1))";
#else
    cerr << "Unknown bootstrap method (available methods are ABC (0) or BCa (1))";
#endif
  }
  
  return 0; 
}

// Function that computes the finite differences
// sign: whether to compute -epsn or +epsn (set as -1 or 1)
// nloc: the number of loci 
// For progress report: 
//   starting_at_task: an integer between 1 and total_tasks
//   total_tasks: the number of total tasks to perform 
// estimatingFnPtr: a pointer to the function carrying out the estimation
vector<double> get_nearby_epsn(double sign, size_t nloc, 
                               size_t starting_at_task, 
                               size_t total_tasks, 
                               size_t consy, 
                               double (*estimatingFnPtr)(vector<double> d)) { 
  
  // This function will report progress. Note that progress here does not 
  // have to start at zero. When this function needs to be called twice 
  // (as in the ABC method), we jump straight to 50% completion. 
#ifdef COMPATIBILITYRCPP
  
  Progress progbar(total_tasks); 
  progbar.update(starting_at_task); 
#else 
  size_t progress = starting_at_task; 
#endif
  
  vector<double> weights(nloc, (1.0 - epsn_value)/nloc); 
  vector<double> results(nloc); 
  
  weights[0] += sign * epsn_value; // this defines the same weights as in boot:::inf.jack SO
  // following the same code, the empinf l_j are (<resulting t>-t0)/epsn_value
  results[0] = (*estimatingFnPtr)(weights);
	for (size_t current_loc=1; current_loc<nloc; current_loc++) {
    // Set previous weight to its original value
    weights[current_loc - 1] -= sign * epsn_value;
    // Adjust weight for current loc
    weights[current_loc] += sign * epsn_value;
    
    results[current_loc] = (*estimatingFnPtr)(weights);
#ifdef COMPATIBILITYRCPP
    progbar.increment();
#else
    if (use_console_xy) {
      _gotoxy(0, consy);
      noR_cout << 
        " Computing confidence interval... about      % done             "  <<
          flush;
      _gotoxy(41, consy);
      noR_cout << int(100 * (progress) / ( total_tasks ) ) << flush;
    } 
    progress++; 
#endif
  
  }
  
  return results; 
}

double CABCbootstrap::cancelland(double unidirPvalue) { //wrapper from (double) to (vector<double>)
    vector<double>weights(nb_units);
  lambda_4_13=(z+ndtri(unidirPvalue))/(pow(1-ahat*(z+ndtri(unidirPvalue)),2));
	for(size_t loc=0;loc<nb_units;loc++) //pour borne inf
		weights[loc]=delta[loc]*lambda_4_13+1.0/nb_units;
	double t=(*estimFnPtr)(weights);
	return(t-testPoint); // bisection_search must seek the zero of this. testPoint is global within the file...
}

double cancellandWrapper(double unidirPvalue) {
    //allows a member function to be called through an ordinary function pointer ie as first argument of bootstrapOverLoci
    return(ABCptr->cancelland(unidirPvalue)); // of course one needs a global pointer to the appropriate class
}

double CABCbootstrap::Pvalue(double testPt,bool unidir=false,bool verbose=true) { // gives (by default unidirectional) P-value associated with (currently) slope value
    /** returns a valid P value, or else
        NaN -> invalid call, no valid CI previouly computed
        -1 -> failure of bisection search
    **/

    vector<double> resu(1,-1);
    double guess=1;
    ABCptr=this;
    testPoint=testPt;  // for later uses outside Pvalue()
    double sig,level1=0.0,level2=0.0,pvalue;
    if(std::isnan(z)) {
        noR_cout<<"Attempt to compute P value by bootstrap\n    while confidence interval computation was not called, or failed.";
	     //(need to have bootOut open... )
	     //bootOut<<"Attempt to compute P value by ABC bootstrap\n    while confidence interval computation was not called, or failed.";
	     //bootOut.close();
	   	 if (pauseGP) {

            noR_cout<<"\n(Return) to continue"<<endl; getchar();

        }
	   	 return(numeric_limits<double>::quiet_NaN()); // not exit() !
    } /**ELSE **/

        noR_cout<<" Computing test"<<testLegend<<"= "<<testPoint<<"; beginning..";


    int it=1;
    while(resu.size()==1 && it <50) {  //resu.size()=1means bad starting values
        // note that resu.size()=2 either when solution is found or for bisection failure to converge (visible NaN)
        if (verbose && it==2) {
                noR_cout<<"(*) From Pvalue(): Problem finding starting values for bisection search"<<endl;
                noR_cout<<"tinf, t0, tsup were "<<tinf<<" "<<t0<<" "<<tsup<<endl;
                noR_cout<<"Initial levels were "<<level1<<" "<<level2<<endl;
        }
        guess*=10;
        /** just initial values... **/
        if (testPt<t0) {
            if (tinf<testPt) { // testtPt in ]tinf,t0[
                level1=seuil_inf;level2=0.5+0.01*it;
            } else { // testtPt in ]-inf,tinf]
                sig=(t0-tinf)/(-ndtri(seuil_inf)); // 1 if N(t0,1) CI => estimation sigma in N(t0,sigma2) CI
                //seuil_inf = 0.025 du CI would result in unidirPvalue=0.05
                level2=1-(1-0.01*it)*(1-2*seuil_inf);
                //level1 better not =0...
                level1=max(ndtr(guess*(testPt-t0)/sig),min(level2/2,pow(0.01,double(it)/5+1)));
            }
        } else {
            if (testPt<tsup) { // testtPt in [t0,tsup[
                level1=0.5-0.01*it;level2=seuil_sup;
            } else { // testtPt in [tsup,inf[
                //seuil_sup = 0.975 du CI would result in unidirPvalue=0.95
                level1=(1-0.01*it)*(1-2*(1-seuil_sup));
                sig=(tsup-t0)/(ndtri(seuil_sup)); // 1 if N(t0,1) CI => estimation sigma in N(t0,sigma2) CI
                //level2 better not =1...
                level2=min(ndtr(guess*(testPt-t0)/sig),max(1.-(1-level1)/2,1-pow(0.01,double(it)/5+1)));
            }
        }
        if (it>1 && verbose) {

            noR_cout<<"New initial levels "<<level1<<" "<<level2<<endl;

        }
        it++;
        /** main stuff **/
        resu=bisection_search(cancellandWrapper,level1,level2,verbose);
    }
    if (resu.size()==1) { // Exited while() without good starting values
        if (verbose) {
                noR_cout<<"(!) From Pvalue(): Failed to find starting values for bisection search";
                noR_cout<<"tinf, t0, tsup were "<<tinf<<" "<<t0<<" "<<tsup<<endl;
        }
        // now a reasonable guess as to why good starting values were not found...
        if (testPt<tinf) resu.push_back(0);
        else if (tsup<testPt) resu.push_back(1);
        else resu.push_back(numeric_limits<double>::quiet_NaN()); // no clear reason to reach this point => make it visible
    }
    if ( ! unidir) pvalue=2*(min(resu[1],1.-resu[1])); else pvalue=resu[1];
    return(pvalue);
}
