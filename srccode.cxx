#include <XSFunctions/Utilities/funcType.h>
#include <XSUtil/Numerics/Integrate.h>
#include <XSUtil/Numerics/Numerics.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSFunctions/functionMap.h>
#include "xsTypes.h"
#include <cmath>


extern "C" void NTsed(const RealArray& energyArray, 
                        const RealArray& params, 
                        int spectrumNumber, 
                        RealArray& fluxArray, 
                        RealArray& fluxErrArray, 
                        const string& initString)
{
    using namespace std;
    fluxArray.resize(energyArray.size()-1);
    fluxErrArray.resize(0);
    
    // Read in the model parameters
    Real mbh(pow(10., params[0]));
    Real mdot(params[1]);
    Real spin(params[2]);
    Real rin(params[3]);
    Real rout(params[4]);
    Real incl(params[5]);
    Real fcol(params[6]);
    Real dist(params[9]);

    // Set up relevant physical constants
    const Real gc(6.6743015e-11);
    const Real msol(1.98847e30);
    const Real hp(6.62607015e-34);
    const Real cl(2.99792458e8);
    const Real kb(1.380649e-23);
    const Real kbcgs(8.617333262e-5);
    const Real sb(5.670374419e-8);
    const Real ec(1.602176634e-19);
    const Real pc(3.0856775814913673e16);
    const Real mp(1.67262192369e-27);
    const Real st(6.6524587321e-29);

    // Compute the innermost stable circular orbit (ISCO)
    Real spinsign(1.0);
    if (signbit(spin) == true)
    {
        // if spin is < 0 then use a -ve in RMS calculation
        spinsign = -1.0;
    } else if (signbit(spin) == false)
    {
        // if spin is > 0 then use a +ve in RMS calculation
        spinsign = 1.0;
    }
    Real z1(1.+(pow(1.-pow(spin, 2.), 1./3.))*(pow(1+spin, 1./3.)+pow(1-spin, 1./3.)));
    Real z2(pow((3.*pow(spin, 2.)+pow(z1, 2.)), 1./2.));
    Real rms(3.+z2-spinsign*pow((3.-z1)*(3.+z1+2.*z2), 1./2.));

    // Compute the efficiency
    Real x0(pow(rms, 1./2.));
    Real fC(1.-3.*pow(x0, -2.)+2.*spin*pow(x0, -3.));
    Real fG(1.-2.*pow(x0, -2.)+spin*pow(x0, -3.));
    Real eta(1.-pow(fC, -1./2.)*fG);

    // Convert quantities into more useful units
    mbh = mbh*msol;
    Real LEdd(4.*M_PI*gc*mbh*mp*cl/st);
    mdot = (mdot*LEdd)/(pow(cl, 2.)*eta);
    dist = dist*1.e6*pc;
    Real rg((gc*mbh)/(pow(cl, 2.)));

    // If the rin parameter is -ve set it to the ISCO
    if (rin < 0.)
    {
        rin = rms;
    }

    // Set up the radial bins
    Real rstart(log10(rms));
    Real rstop(log10(rout));
    size_t rnum(1001);  // 1000 radial bins offers about dr/r of ~1% for Rin = 1 and Rout = 1e5
    Real rstep((rstop-rstart)/(rnum-1));
    RealArray redges(rnum);
    for (size_t n = 0; n < rnum; ++n){
        redges[n] = pow(10., rstart + n*rstep);
    }

    // Compute each radial bin mid-point and dr w.r.t. next bin
    RealArray rmidps(rnum-1);
    RealArray rdiff(rnum-1);
    for (size_t n = 0; n < rnum-1; ++n){
        rmidps[n] = pow(10, ((log10(redges[n+1])+log10(redges[n]))/2.));
        rdiff[n] = redges[n+1]-redges[n];
    }

    // Set up temperature profile and disc flux / luminosity
    RealArray tztr(rnum-1);
    Real tmax(0.);    

    // Compute the luminosity of the hot corona following eqn 2 of Kubota & Done (2018)
    // The hot corona is powered by the luminosity between the Rrms and Rin radii
    Real Ldiss(0.);
    int k = 0;
    while (redges[k+1] <= rin)
    {
        // Temperature profile for a Novikov-Thorne disc emissivity using the notation of Page & Thorne (1974)
        Real x(pow(rmidps[k], 1./2.));
        Real x0(pow(rms, 1./2.));
        Real x1(2.*cos((1./3.)*acos(spin)-M_PI/3.));
        Real x2(2.*cos((1./3.)*acos(spin)+M_PI/3.));
        Real x3(-2.*cos((1./3.)*acos(spin)));
        Real prefactor1((3.*mdot*pow(cl, 6.))/(8.*M_PI*pow(gc, 2.)*pow(mbh, 2.)));
        Real prefactor2(pow((pow(x, 7.)-3.*pow(x, 5.)+2.*spin*pow(x, 4.)), -1.));
        Real bigterm1(((3.*pow((x1-spin), 2.))/(x1*(x1-x2)*(x1-x3))) * log((x-x1)/(x0-x1)));
        Real bigterm2(((3.*pow((x2-spin), 2.))/(x2*(x2-x1)*(x2-x3))) * log((x-x2)/(x0-x2)));
        Real bigterm3(((3.*pow((x3-spin), 2.))/(x3*(x3-x1)*(x3-x2))) * log((x-x3)/(x0-x3)));
        Real fx(prefactor1*prefactor2*(x-x0-(3./2.)*spin*log(x/x0)-bigterm1-bigterm2-bigterm3));
        Real t(pow(fx/sb, 1./4.));
        tztr[k] = t;
        Real r(rmidps[k]*rg);
        Real dr(rdiff[k]*rg);
        Ldiss += 2*sb*pow(tztr[k], 4.)*2.*M_PI*r*dr;

        k++;
    }

    // Compute the accretion disc spectrum
    for (size_t j = k; j < rnum-1; ++j){
        // Temperature profile for a Novikov-Thorne disc emissivity using the notation of Page & Thorne (1974)
        Real x(pow(rmidps[j], 1./2.));
        Real x0(pow(rms, 1./2.));
        Real x1(2.*cos((1./3.)*acos(spin)-M_PI/3.));
        Real x2(2.*cos((1./3.)*acos(spin)+M_PI/3.));
        Real x3(-2.*cos((1./3.)*acos(spin)));
        Real prefactor1((3.*mdot*pow(cl, 6.))/(8.*M_PI*pow(gc, 2.)*pow(mbh, 2.)));
        Real prefactor2(pow((pow(x, 7.)-3.*pow(x, 5.)+2.*spin*pow(x, 4.)), -1.));
        Real bigterm1(((3.*pow((x1-spin), 2.))/(x1*(x1-x2)*(x1-x3))) * log((x-x1)/(x0-x1)));
        Real bigterm2(((3.*pow((x2-spin), 2.))/(x2*(x2-x1)*(x2-x3))) * log((x-x2)/(x0-x2)));
        Real bigterm3(((3.*pow((x3-spin), 2.))/(x3*(x3-x1)*(x3-x2))) * log((x-x3)/(x0-x3)));
        Real fx(prefactor1*prefactor2*(x-x0-(3./2.)*spin*log(x/x0)-bigterm1-bigterm2-bigterm3));
        Real t(pow(fx/sb, 1./4.));

        // Compute Done et al. (2012) colour-temperature correction if fcol is -ve
        Real ffcol(fcol);
        if (fcol < 0.)
        {
            if (t >= 1.e5)
            {
                ffcol = pow((72./(kbcgs*t/1.e3)), 1./9.);
            } else if (t > 3.e4 && t < 1.e5)
            {
                ffcol = pow((t/3.e4), 0.82);
            } else if (t <= 3.e4)
            {
                ffcol = 1.;
            }
        } 
        tztr[j] = t*ffcol;
        
        // Save the max disc temperature
        if (tztr[j] > tmax)
        {
            tmax = tztr[j];
        }

        // Compute the spectrum emitted by the current annulus of the disc, summing
        Real r(rmidps[j]*rg);
        Real dr(rdiff[j]*rg);
        for (size_t i = 0; i < fluxArray.size(); ++i){
            Real dE(energyArray[i+1]-energyArray[i]);
            Real midE(pow(10., (log10(energyArray[i])+log10(energyArray[i+1]))/2.));
            fluxArray[i] += ((2.*hp*pow((midE*1.e3*ec)/hp, 3.))/(pow(cl, 2.))) * (1./(exp((hp*(midE*1.e3*ec)/hp)/(kb*tztr[j]))-1.))*((2.*M_PI*r*dr*cos(incl*M_PI/180.))/(pow(dist, 2.)))*(1./(pow(ffcol, 4.))) * (1.e7/1.e4) * (1.51e26/midE)*dE;
        }
    }

    // Account for redshift on the disc emission
    const RealArray z = {params[10]};
    zashift(energyArray, z, spectrumNumber, fluxArray, fluxErrArray, initString);

    // Compute the hot corona emission using nthComp
    // ---> NOTE: z in nthcomp only controls energy shift, not flux
    // ---> NOTE: do not use zashift on nthcomp, use combination of 
    //      clumin and z(nthcomp) to achieve both energy and flux z
    //      correction
    RealArray nthfluxArray;
    RealArray nthfluxErrArray;
    const RealArray hotpars = {params[7], params[8], tmax*kbcgs/1.e3, 0, params[10]};
    nthcomp(energyArray, hotpars, spectrumNumber, nthfluxArray, nthfluxErrArray, initString);
    
    // Set the nthcomp flux to the dissipated luminosity between Rrms and Rin
    Real nthflx( log10((Ldiss*1.e7)/(4.*M_PI*(pow((1.+params[10])*(params[9]*1.e6*3.0856775814913673e16*1.e2), 2.)))) );
    const RealArray cflpars = {1.e-5, 1.e3, nthflx};
    cflux(energyArray, cflpars, spectrumNumber, nthfluxArray, nthfluxErrArray, initString);
    
    // Add the nthcomp component to the overall model
    for (size_t i = 0; i < fluxArray.size(); ++i){
        fluxArray[i] += nthfluxArray[i];
    }

    return;
}