#ifndef TCONST_H
#define TCONST_H

//
//  Constants
//

const int MAXDIM =  200;              // Maximum no of elements in a card entry
const int MAXCENTRBIN = 10;           // Maximum no of centrality bins
const int MAXASSOCBIN = 10;           // Maximum no of associated pt bins
const int MAXTRIGGBIN = 10;           // Maximum no of trigger pt bins
const int MAXREACPLBIN = 7;           // Maximum no of reaction plane bins
const int MAXDPHIBIN  = 30;           // Maximum no of delta phi bins

const double EPS = 0.001;
const double pi           = 3.14159265358979; 
const double twopi        = 2*pi;

const double ElectronMass = .510999e-3;
const double MuonMass     = .105658;
const double PionMass     = .139570;
const double PizeroMass   = .134977;
const double KaonMass     = .493677;
const double KzeroMass    = .497672;
const double ProtonMass   = .938272;
const double NeutronMass  = .939565;

enum fillType {real, mixed};
enum corrType {trigg, assoc, centr, reacpl}; 
enum partType {unknown, hadron, proton, kaon, pion, photon, muon, electron, diphoton, dimuon, dielectron};
enum TEMC {PbSc, PbGl}; 

#endif
