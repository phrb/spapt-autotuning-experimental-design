G4double G4hPairProductionModel::ComputeDMicroscopicCrossSection(
                                           G4double tkin,
                                           G4double Z,
                                           G4double pairEnergy)
//  differential cross section
{
  G4double bbbtf= 183. ;
  G4double bbbh = 202.4 ;
  G4double g1tf = 1.95e-5 ;
  G4double g2tf = 5.3e-5 ;
  G4double g1h  = 4.4e-5 ;
  G4double g2h  = 4.8e-5 ;

  G4double totalEnergy  = tkin + particleMass;
  G4double residEnergy  = totalEnergy - pairEnergy;
  G4double massratio    = particleMass/electron_mass_c2 ;
  G4double massratio2   = massratio*massratio ;
  G4double cross = 0.;

  SetElement(G4lrint(Z));

  G4double c3 = 0.75*sqrte*particleMass;
  if (residEnergy <= c3*z13) return cross;

  G4double c7 = 4.*electron_mass_c2;
  G4double c8 = 6.*particleMass*particleMass;
  G4double alf = c7/pairEnergy;
  G4double a3 = 1. - alf;
  if (a3 <= 0.) return cross;

  // zeta calculation
  G4double bbb,g1,g2;
  if( Z < 1.5 ) { bbb = bbbh ; g1 = g1h ; g2 = g2h ; }
  else          { bbb = bbbtf; g1 = g1tf; g2 = g2tf; }

  G4double zeta = 0;
  G4double zeta1 = 0.073*log(totalEnergy/(particleMass+g1*z23*totalEnergy))-0.26;
  if ( zeta1 > 0.)
  {
    G4double zeta2 = 0.058*log(totalEnergy/(particleMass+g2*z13*totalEnergy))-0.14;
    zeta  = zeta1/zeta2 ;
  }

  G4double z2 = Z*(Z+zeta);
  G4double screen0 = 2.*electron_mass_c2*sqrte*bbb/(z13*pairEnergy);
  G4double a0 = totalEnergy*residEnergy;
  G4double a1 = pairEnergy*pairEnergy/a0;
  G4double bet = 0.5*a1;
  G4double xi0 = 0.25*massratio2*a1;
  G4double del = c8/a0;

  G4double rta3 = sqrt(a3);
  G4double tmnexp = alf/(1. + rta3) + del*rta3;
  if(tmnexp >= 1.0) return cross;

  G4double tmn = log(tmnexp);
  G4double sum = 0.;

  // Gaussian integration in ln(1-ro) ( with 8 points)
  for (G4int i=0; i<8; i++)
  {
    G4double a4 = exp(tmn*xgi[i]);     // a4 = (1.-asymmetry)
    G4double a5 = a4*(2.-a4) ;
    G4double a6 = 1.-a5 ;
    G4double a7 = 1.+a6 ;
    G4double a9 = 3.+a6 ;
    G4double xi = xi0*a5 ;
    G4double xii = 1./xi ;
    G4double xi1 = 1.+xi ;
    G4double screen = screen0*xi1/a5 ;
    G4double yeu = 5.-a6+4.*bet*a7 ;
    G4double yed = 2.*(1.+3.*bet)*log(3.+xii)-a6-a1*(2.-a6) ;
    G4double ye1 = 1.+yeu/yed ;
    G4double ale=log(bbb/z13*sqrt(xi1*ye1)/(1.+screen*ye1)) ;
    G4double cre = 0.5*log(1.+2.25*z23*xi1*ye1/massratio2) ;
    G4double be;

    if (xi <= 1.e3) be = ((2.+a6)*(1.+bet)+xi*a9)*log(1.+xii)+(a5-bet)/xi1-a9;
    else            be = (3.-a6+a1*a7)/(2.*xi);

    G4double fe = (ale-cre)*be;
    if ( fe < 0.) fe = 0. ;

    G4double ymu = 4.+a6 +3.*bet*a7 ;
    G4double ymd = a7*(1.5+a1)*log(3.+xi)+1.-1.5*a6 ;
    G4double ym1 = 1.+ymu/ymd ;
    G4double alm_crm = log(bbb*massratio/(1.5*z23*(1.+screen*ym1)));
    G4double a10,bm;
    if ( xi >= 1.e-3)
    {
      a10 = (1.+a1)*a5 ;
      bm  = (a7*(1.+1.5*bet)-a10*xii)*log(xi1)+xi*(a5-bet)/xi1+a10;
    } else {
      bm = (5.-a6+bet*a9)*(xi/2.);
    }

    G4double fm = alm_crm*bm;
    if ( fm < 0.) fm = 0. ;

    sum += wgi[i]*a4*(fe+fm/massratio2);
  }

  cross = -tmn*sum*factorForCross*z2*residEnergy/(totalEnergy*pairEnergy);

  return cross;
}

