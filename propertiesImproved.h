#ifndef PROPERTIESNEUTRAL_H
#define PROPERTIESNEUTRAL_H

#include "TMath.h"

namespace Eric
{
    constexpr Int_t returnNumberOfSteps(Int_t min, Int_t max, Double_t dStep)
    {
        return static_cast<Int_t>((max - min) / dStep + 1.0);
    }
}

namespace Constants
{
    constexpr double c{2.99792458e8};   //m/s
    constexpr double NA{6.022e23};        // mol⁻¹
    constexpr double me_c2 = 0.000511;     // GeV
    constexpr double prefactor = 256.0 / (9.0 * TMath::Pi());
    constexpr double angstromcmRatio{1.0e-8};
}

namespace Config {
    constexpr double thetaCrystal = 0.010;   // rad
    constexpr double lengthCrystal = 0.10;   // m

    constexpr double Ltwo{0.03}; //m
    constexpr double ltwo{0.0};

     constexpr double calculateXc(double Dp, double Atf)
    {
        return (Dp / 2 - Atf);
    }

    constexpr double calculateN (double rho, double A)
    {
        return Constants::NA * rho / A;    
    }

    constexpr double angstromToCm(double numberInAngstrom)
    {
        return numberInAngstrom * 1.0e8;
    }
}

namespace Silicium
    {
    // material properties Si(110)
    constexpr double d = 1.92;                 // Å 
    constexpr double aTF = 0.194;             // Å
    constexpr double Uxc = 16e-9;            // GeV
    constexpr double Uprime = 5.7;          // GeV/cm U'
    constexpr double I = 173e-9;             // GeV         
    constexpr double rho = 2.329;            // g/cm³
    constexpr double A = 28.0855;            // g/mol    
    constexpr double Z = 14.0;
    constexpr double aTF_cm{Config::angstromToCm(aTF)}; //cm
    constexpr double d_cm{Config::angstromToCm(d)}; //cm
    constexpr double N{Config::calculateN(Silicium::rho, Silicium::A)};
    constexpr double xc{Config::calculateXc(Silicium::d, Silicium::aTF)};
    }
    
namespace Germanium
    {
    //material properties Ge(110)
    constexpr double d = 2.00;   // Å 
    constexpr double aTF = 0.148;             // Å
    constexpr double Uxc = 27e-9;            // GeV
    constexpr double Uprime = 10.0;          // GeV/cm U'
    constexpr double I = 350e-9;             // GeV         
    constexpr double rho = 5.323;            // g/cm³
    constexpr double A = 72.630;            // g/mol    
    constexpr double Z = 32.0;
    constexpr double aTF_cm{Config::angstromToCm(aTF)}; //cm
    constexpr double d_cm{Config::angstromToCm(d)}; //cm
    constexpr double N{Config::calculateN(Germanium::rho, Germanium::A)};       
    constexpr double xc{Config::calculateXc(Germanium::d, Germanium::aTF)};
    }

    

#endif
