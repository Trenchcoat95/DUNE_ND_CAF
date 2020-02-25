#include "CAF.C"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "EVGCore/EventRecord.h"
//#include "nusystematics/artless/response_helper.hh"
#include <stdio.h>

TRandom3 * rando;
const double mmu = 0.1056583745;
//TF1 * tsmear; // angular resolution function

// params will be extracted from command line, and passed to the reconstruction
struct params {
  double OA_xcoord;
  bool fhc, grid, IsGasTPC;
  int seed, run, subrun, first, n, nfiles;
  double trk_muRes, LAr_muRes, ECAL_muRes;
  double em_const, em_sqrtE;
  double michelEff;
  double CC_trk_length;
  double pileup_frac, pileup_max;
  double lartpc_len, gastpc_len, gastpc_B, gastpc_padPitch, gastpc_X0;
};

// Fill reco variables for muon reconstructed in magnetized tracker
void recoMuonTracker( CAF &caf, params &par )
{
  // smear momentum by resolution
  double p = sqrt(caf.LepE*caf.LepE - mmu*mmu);
  double reco_p = rando->Gaus( p, p*par.trk_muRes );
  caf.Elep_reco = sqrt(reco_p*reco_p + mmu*mmu);


  // assume perfect charge reconstruction
  caf.reco_q = (caf.LepPDG > 0 ? -1 : 1);

  // assume always muon for tracker-matched
  caf.reco_numu = 1; caf.reco_nue = 0; caf.reco_nc = 0;
  caf.muon_contained = 0; caf.muon_tracker = 1; caf.muon_ecal = 0; caf.muon_exit = 0;
  //caf.Ev_reco = caf.Elep_reco;
}

// Fill reco muon variables for muon contained in LAr
void recoMuonLAr( CAF &caf, params &par )
{
  // range-based, smear kinetic energy
  double ke = caf.LepE - mmu;
  double reco_ke = rando->Gaus( ke, ke*par.LAr_muRes );
  caf.Elep_reco = reco_ke + mmu;



  // assume negative for FHC, require Michel for RHC
  if( par.fhc ) caf.reco_q = -1;
  else {
    double michel = rando->Rndm();
    if( caf.LepPDG == -13 && michel < par.michelEff ) caf.reco_q = 1; // correct mu+
    else if( caf.LepPDG == 13 && michel < par.michelEff*0.25 ) caf.reco_q = 1; // incorrect mu-
    else caf.reco_q = -1; // no reco Michel
  }

  caf.reco_numu = 1; caf.reco_nue = 0; caf.reco_nc = 0;
  caf.muon_contained = 1; caf.muon_tracker = 0; caf.muon_ecal = 0; caf.muon_exit = 0;
  //caf.Ev_reco = caf.Elep_reco;
}

// Fill reco variables for muon reconstructed in magnetized tracker
void recoMuonECAL( CAF &caf, params &par )
{
  // range-based KE
  double ke = caf.LepE - mmu;
  double reco_ke = rando->Gaus( ke, ke*par.ECAL_muRes );
  caf.Elep_reco = reco_ke + mmu;


  // assume perfect charge reconstruction -- these are fairly soft and should curve a lot in short distance
  caf.reco_q = (caf.LepPDG > 0 ? -1 : 1);

  // assume always muon for ecal-matched
  caf.reco_numu = 1; caf.reco_nue = 0; caf.reco_nc = 0;
  caf.muon_contained = 0; caf.muon_tracker = 0; caf.muon_ecal = 1; caf.muon_exit = 0;
 

// Fill reco variables for true electron
void recoElectron( CAF &caf, params &par )
{
  caf.reco_q = 0; // never know charge
  caf.reco_numu = 0;
  caf.muon_contained = 0; caf.muon_tracker = 1; caf.muon_ecal = 0; caf.muon_exit = 0;

  // fake efficiency...threshold of 300 MeV, eff rising to 100% by 700 MeV
  if( rando->Rndm() > (caf.LepE-0.3)*2.5 ) { // reco as NC
    caf.Elep_reco = 0.;
    caf.reco_nue = 0; caf.reco_nc = 1;
  } else { // reco as CC
    caf.Elep_reco = rando->Gaus( caf.LepE, caf.LepE*(par.em_const + par.em_sqrtE/sqrt(caf.LepE)) );
    caf.reco_nue = 1; caf.reco_nc = 0;
  }


}

void decayPi0( TLorentzVector pi0, TVector3 &gamma1, TVector3 &gamma2 )
{
  double e = pi0.E();
  double mp = 134.9766; // pi0 mass

  double beta = sqrt( 1. - (mp*mp)/(e*e) ); // velocity of pi0
  double theta = 3.1416*rando->Rndm(); // theta of gamma1 w.r.t. pi0 direction
  double phi = 2.*3.1416*rando->Rndm(); // phi of gamma1 w.r.t. pi0 direction

  double p = mp/2.; // photon momentum in pi0 rest frame
  TLorentzVector g1( 0., 0., p, p ); // pre-rotation photon 1
  TLorentzVector g2( 0., 0., -p, p ); // pre-rotation photon 2 is opposite

  // rotate to the random decay axis in pi0 rest frame. choice of rotation about x instead of y is arbitrary
  g1.RotateX( theta );
  g2.RotateX( theta );
  g1.RotateZ( phi );
  g2.RotateZ( phi );

  // boost to lab frame with pi0 velocity. pi0 direction is z axis for this
  g1.Boost( 0., 0., beta );
  g2.Boost( 0., 0., beta );

  // make gamma1 the more energetic one
  if( g1.E() > g2.E() ) {
    gamma1 = g1.Vect();
    gamma2 = g2.Vect();
  } else {
    gamma1 = g2.Vect();
    gamma2 = g1.Vect();
  }

  // rotate from frame where pi0 is z' direction into neutrino frame
  TVector3 pi0dir = pi0.Vect().Unit(); // actually w.r.t. neutrino direction
  gamma1.RotateUz( pi0dir );
  gamma2.RotateUz( pi0dir );
}

// main loop function
void loop( CAF &caf, params &par, TTree * tree, std::string ghepdir, std::string fhicl_filename )
{
  // read in edep-sim output file
  int ifileNo, ievt, lepPdg, muonReco, nFS;
  float lepKE, muGArLen, hadTot, hadCollar;
  float hadP, hadN, hadPip, hadPim, hadPi0, hadOther;
  float p3lep[3], vtx[3], muonExitPt[3], muonExitMom[3];
  int fsPdg[100];
  float fsPx[100], fsPy[100], fsPz[100], fsE[100], fsTrkLen[100], fsTrkLenPerp[100];
  tree->SetBranchAddress( "ifileNo", &ifileNo );
  tree->SetBranchAddress( "ievt", &ievt );
  tree->SetBranchAddress( "lepPdg", &lepPdg );
  tree->SetBranchAddress( "muonReco", &muonReco );
  tree->SetBranchAddress( "lepKE", &lepKE );
  tree->SetBranchAddress( "muGArLen", &muGArLen );
  tree->SetBranchAddress( "hadTot", &hadTot );
  tree->SetBranchAddress( "hadCollar", &hadCollar );
  tree->SetBranchAddress( "hadP", &hadP );
  tree->SetBranchAddress( "hadN", &hadN );
  tree->SetBranchAddress( "hadPip", &hadPip );
  tree->SetBranchAddress( "hadPim", &hadPim );
  tree->SetBranchAddress( "hadPi0", &hadPi0 );
  tree->SetBranchAddress( "hadOther", &hadOther );
  tree->SetBranchAddress( "p3lep", p3lep );
  tree->SetBranchAddress( "vtx", vtx );
  tree->SetBranchAddress( "muonExitPt", muonExitPt );
  tree->SetBranchAddress( "muonExitMom", muonExitMom );
  tree->SetBranchAddress( "nFS", &nFS );
  tree->SetBranchAddress( "fsPdg", fsPdg );
  tree->SetBranchAddress( "fsPx", fsPx );
  tree->SetBranchAddress( "fsPy", fsPy );
  tree->SetBranchAddress( "fsPz", fsPz );
  tree->SetBranchAddress( "fsE", fsE );
  tree->SetBranchAddress( "fsTrkLen", fsTrkLen );
  tree->SetBranchAddress( "fsTrkLenPerp", fsTrkLenPerp );
  // Get GHEP file for genie::EventRecord from other file
  int current_file = -1;
  TFile * ghep_file = NULL;
  TTree * gtree = NULL;

  std::string mode = ( par.fhc ? "neutrino" : "antineutrino" );


  caf.pot = 0.;

  // Main event loop
  int N = tree->GetEntries();
  if( par.n > 0 && par.n < N ) N = par.n + par.first;
  for( int ii = par.first; ii < N; ++ii ) {

    tree->GetEntry(ii);
    if( ii % 100 == 0 ) printf( "Event %d of %d...\n", ii, N );

    caf.setToBS();

    // set defaults
    for( int j = 0; j < 100; ++j ) {
      caf.nwgt[j] = 7;
      caf.cvwgt[j] = ( caf.iswgt[j] ? 1. : 0. );
      for( unsigned int k = 0; k < 100; ++k ) {
        caf.wgt[j][k] = ( caf.iswgt[j] ? 1. : 0. );
      }
    }

    // make sure ghep file matches the current one, otherwise update to the current ghep file
    if( ifileNo != current_file ) {
      // close the previous file
      if( ghep_file ) ghep_file->Close();

      if( par.grid ) ghep_file = new TFile( Form("genie.%d.root", ifileNo) );
      else if( !par.IsGasTPC ) ghep_file = new TFile( Form("%s/LAr.%s.%d.ghep.root", ghepdir.c_str(),  mode.c_str(), ifileNo) );
      else ghep_file = new TFile( Form("%s/GAr.%s.%d.ghep.root", ghepdir.c_str(), mode.c_str(), ifileNo) );
      
      gtree = (TTree*) ghep_file->Get( "gtree" );
      caf.pot += gtree->GetWeight();
      printf( "New GHEP file with %g POT, total = %g\n", gtree->GetWeight(), caf.pot );

      // can't find GHepRecord
      if( gtree == NULL ) {
        printf( "Can't find ghep event record for file %d!!!\n", ifileNo );
        continue;
      }

      gtree->SetBranchAddress( "gmcrec", &caf.mcrec );
      current_file = ifileNo;
    }

    caf.vtx_x = vtx[0];
    caf.vtx_y = vtx[1];
    caf.vtx_z = vtx[2]; 
    caf.det_x = par.OA_xcoord;

    // configuration variables in CAF file; we don't use mvaresult so just set it to zero
    caf.run = par.run;
    caf.subrun = par.subrun;
    caf.event = ii;
    caf.isFD = 0;
    caf.isFHC = par.fhc;

    // get GENIE event record
    gtree->GetEntry( ievt );
    genie::EventRecord * event = caf.mcrec->event;
    genie::Interaction * in = event->Summary();
   
    
    // Get truth stuff out of GENIE ghep record
    caf.neutrinoPDG = in->InitState().ProbePdg();
    caf.neutrinoPDGunosc = in->InitState().ProbePdg(); // fill this for similarity with FD, but no oscillations
    caf.mode = in->ProcInfo().ScatteringTypeId();
    caf.Ev = in->InitState().ProbeE(genie::kRfLab);
    caf.interaction_type = in->ProcInfo().IsQuasiElastic();
    caf.LepPDG = in->FSPrimLeptonPdg();
    caf.isCC = (abs(caf.LepPDG) == 13 || abs(caf.LepPDG) == 11);
    TLorentzVector lepP4;
    TLorentzVector nuP4nuc = *(in->InitState().GetProbeP4(genie::kRfHitNucRest));
    TLorentzVector nuP4 = *(in->InitState().GetProbeP4(genie::kRfLab));

    caf.nP = 0;
    caf.nN = 0;
    caf.nipip = 0;
    caf.nipim = 0;
    caf.nipi0 = 0;
    caf.nikp = 0;
    caf.nikm = 0;
    caf.nik0 = 0;
    caf.niem = 0;
    caf.niother = 0;
    caf.nNucleus = 0;
    caf.nUNKNOWN = 0; // there is an "other" category so this never gets used
    caf.eP = 0.;
    caf.eN = 0.;
    caf.ePip = 0.;
    caf.ePim = 0.;
    caf.ePi0 = 0.;
    caf.eOther = 0.;
    caf.eRecoP = 0.;
    caf.eRecoN = 0.;
    caf.eRecoPip = 0.;
    caf.eRecoPim = 0.;
    caf.eRecoPi0 = 0.;
    caf.eOther = 0.;
    for( int i = 0; i < nFS; ++i ) {
      double ke = 0.001*(fsE[i] - sqrt(fsE[i]*fsE[i] - fsPx[i]*fsPx[i] - fsPy[i]*fsPy[i] - fsPz[i]*fsPz[i]));
      if( fsPdg[i] == caf.LepPDG ) {
        lepP4.SetPxPyPzE( fsPx[i]*0.001, fsPy[i]*0.001, fsPz[i]*0.001, fsE[i]*0.001 );
        caf.LepE = fsE[i]*0.001;
      }
      else if( fsPdg[i] == 2212 ) {caf.nP++; caf.eP += ke;}
      else if( fsPdg[i] == 2112 ) {caf.nN++; caf.eN += ke;}
      else if( fsPdg[i] ==  211 ) {caf.nipip++; caf.ePip += ke;}
      else if( fsPdg[i] == -211 ) {caf.nipim++; caf.ePim += ke;}
      else if( fsPdg[i] ==  111 ) {caf.nipi0++; caf.ePi0 += ke;}
      else if( fsPdg[i] ==  321 ) {caf.nikp++; caf.eOther += ke;}
      else if( fsPdg[i] == -321 ) {caf.nikm++; caf.eOther += ke;}
      else if( fsPdg[i] == 311 || fsPdg[i] == -311 || fsPdg[i] == 130 || fsPdg[i] == 310 ) {caf.nik0++; caf.eOther += ke;}
      else if( fsPdg[i] ==   22 ) {caf.niem++; caf.eOther += ke;}
      else if( fsPdg[i] > 1000000000 ) caf.nNucleus++;
      else {caf.niother++; caf.eOther += ke;}
    }

    // true 4-momentum transfer
    TLorentzVector q = nuP4-lepP4;

    // Q2, W, x, y frequently do not get filled in GENIE Kinematics object, so calculate manually
    caf.Q2 = -q.Mag2();
    caf.W = sqrt(0.939*0.939 + 2.*q.E()*0.939 + q.Mag2()); // "Wexp"
    caf.X = -q.Mag2()/(2*0.939*q.E());
    caf.Y = q.E()/caf.Ev;

    caf.theta_reco = -1.; // default value

    caf.NuMomX = nuP4.X();
    caf.NuMomY = nuP4.Y();
    caf.NuMomZ = nuP4.Z();
    caf.LepMomX = lepP4.X();
    caf.LepMomY = lepP4.Y();
    caf.LepMomZ = lepP4.Z();
    caf.LepE = lepP4.E();
    caf.LepNuAngle = nuP4.Angle( lepP4.Vect() );


    //--------------------------------------------------------------------------
    // Parameterized reconstruction
    //--------------------------------------------------------------------------
    if( !par.IsGasTPC ) {
	// Loop over final-state particles
	//caf.Ev_reco = 0;
	caf.nFSP = nFS;
	for( int i = 0; i < nFS; ++i ) {
	int pdg = fsPdg[i];
	double p = 0.001*sqrt(fsPx[i]*fsPx[i] + fsPy[i]*fsPy[i] + fsPz[i]*fsPz[i]);
	double KE = 0.001*fsE[i] - sqrt(0.001*0.001*fsE[i]*fsE[i] - p*p);
	double mass = 0.001*sqrt(fsE[i]*fsE[i] - fsPx[i]*fsPx[i] - fsPy[i]*fsPy[i] - fsPz[i]*fsPz[i]);
	caf.trkLen[i] = fsTrkLen[i];
    caf.trkLenPerp[i] = fsTrkLenPerp[i];
	if (fsTrkLen[i] > par.lartpc_len) { 
	double preco = rando->Gaus( p, 0.14 );
    double ereco = sqrt( preco*preco + mass*mass ) - mass;
	caf.preco[i] = preco;
	caf.partEvReco[i] = ereco;
    caf.truepdg[i] = pdg;
	printf( "pdg is %d", pdg );
	caf.truep[i] = p;
	caf.truepx[i] = fsPx[i];
	caf.truepy[i] = fsPy[i];
	caf.truepz[i] = fsPz[i];
	caf.mass[i] = mass;
	caf.ke[i] = KE;
	if( abs(fsPdg[i]) == 211 ) ereco += mass; // add pion mass
    caf.Ev_reco += ereco;
	if( abs(fsPdg[i]) == 211 ) caf.lartpc_pi_pl_mult++; }
    } // I think this closes the FSP loop

    } else {
    caf.nFSP = nFS;
	for( int i = 0; i < nFS; ++i ) {
	double ptrue = 0.001*sqrt(fsPx[i]*fsPx[i] + fsPy[i]*fsPy[i] + fsPz[i]*fsPz[i]);
	double mass = 0.001*sqrt(fsE[i]*fsE[i] - fsPx[i]*fsPx[i] - fsPy[i]*fsPy[i] - fsPz[i]*fsPz[i]);
	double KE = 0.001*fsE[i] - sqrt(0.001*0.001*fsE[i]*fsE[i] - ptrue*ptrue);
	caf.trkLen[i] = fsTrkLen[i];
    caf.trkLenPerp[i] = fsTrkLenPerp[i];
        
	// track length cut 2cm according to Thomas Campbell
	if( fsTrkLen[i] > par.gastpc_len ) {
	double pT = 0.001*sqrt(fsPy[i]*fsPy[i] + fsPz[i]*fsPz[i]); // transverse to B field, in GeV
	double nHits = fsTrkLen[i] / par.gastpc_padPitch; // doesn't matter if not integer as only used in eq
	// Gluckstern formula, sigmapT/pT, with sigmaX and L in meters
	double fracSig_meas = sqrt(720./(nHits+4)) * (0.01*par.gastpc_padPitch/sqrt(12.)) * pT / (0.3 * par.gastpc_B * 0.0001 * fsTrkLenPerp[i]*fsTrkLenPerp[i]);
	// multiple scattering term
	double fracSig_MCS = 0.052 / (par.gastpc_B * sqrt(par.gastpc_X0*fsTrkLenPerp[i]*0.0001));
	
	double sigmaP = ptrue * sqrt( fracSig_meas*fracSig_meas + fracSig_MCS*fracSig_MCS );
	double preco = rando->Gaus( ptrue, sigmaP );
	double ereco = sqrt( preco*preco + mass*mass ) - mass; // kinetic energy
	caf.erecon[i] = ereco; 
	if( abs(fsPdg[i]) == 211) ereco += mass; // add pion mass
	//else if( fsPdg[i] == 2212 && preco > 1.5 ) ereco += 0.1395; // mistake pion mass for high-energy proton
	caf.partEvReco[i] = ereco;
	int pdg = fsPdg[i];
	// threshold cut
	caf.Ev_reco += ereco;
	caf.truepdg[i] = pdg;
	printf("pdg is %d",pdg);
	caf.truep[i] = ptrue;
	caf.truepx[i] = fsPx[i];
	caf.truepy[i] = fsPy[i];
	caf.truepz[i] = fsPz[i];
	caf.mass[i] = mass;
	caf.ke[i] = KE;
	caf.preco[i] = preco;
	//if( fsPdg[i] == 211 || (fsPdg[i] == 2212 && preco > 1.5) ) caf.gastpc_pi_pl_mult++;
	if( fsPdg[i] == 211 ) caf.gastpc_pi_pl_mult++;
	else if( fsPdg[i] == -211 ) caf.gastpc_pi_min_mult++;
	
	if( abs(lepPdg) == 11 ) { // NC or numuCC reco as nueCC
	caf.reco_numu = 0; caf.reco_nue = 1; caf.reco_nc = 0;
	}
	
	if( (fsPdg[i] == 13 || fsPdg[i] == -13) && fsTrkLen[i] > 100. ) { // muon, don't really care about nu_e CC for now
	caf.partEvReco[i] += mass;
	caf.Elep_reco = sqrt(preco*preco + mass*mass);
	// assume perfect charge reconstruction
	caf.reco_q = (fsPdg[i] > 0 ? -1 : 1);
	caf.reco_numu = 1; caf.reco_nue = 0; caf.reco_nc = 0;
	caf.muon_tracker = 1;
	}
	} // this is to close the tracklength > 0 cut 
	else if( fsPdg[i] == 111 || fsPdg[i] == 22 ) {
          double ereco = 0.001 * rando->Gaus( fsE[i], 0.1*fsE[i] );
          caf.partEvReco[i] = ereco;
          //caf.Ev_reco += ereco;
        }
      } //this is to close the FSP loop
    }


    caf.fill();
  }

  // set POT
  caf.meta_run = par.run;
  caf.meta_subrun = par.subrun;

}

int main( int argc, char const *argv[] ) 
{

  if( (argc == 2) && ((std::string("--help") == argv[1]) || (std::string("-h") == argv[1])) ) {
    std::cout << "Help yourself by looking at the source code to see what the options are." << std::endl;
    return 0;
  }

  // get command line options
  std::string ghepdir;
  std::string outfile;
  std::string edepfile;
  std::string fhicl_filename;

  // Make parameter object and set defaults
  params par;
  par.IsGasTPC = false;
  par.OA_xcoord = 0.; // on-axis by default
  par.fhc = true;
  par.grid = false;
  par.seed = 7; // a very random number
  par.run = 1; // CAFAna doesn't like run number 0
  par.subrun = 0;
  par.n = -1;
  par.nfiles = 1;
  par.first = 0;
  par.trk_muRes = 0.02; // fractional muon energy resolution of HP GAr TPC
  par.LAr_muRes = 0.05; // fractional muon energy resolution of muons contained in LAr
  par.ECAL_muRes = 0.1; // fractional muon energy resolution of muons ending in ECAL
  par.em_const = 0.03; // EM energy resolution constant term: A + B/sqrt(E) (GeV)
  par.em_sqrtE = 0.1; // EM energy resolution 1/sqrt(E) term: A + B/sqrt(E) (GeV)
  par.michelEff = 0.75; // Michel finder efficiency
  par.CC_trk_length = 100.; // minimum track length for CC in cm
  par.pileup_frac = 0.1; // fraction of events with non-zero pile-up
  par.pileup_max = 0.5; // GeV
  par.gastpc_len = 2.; // track length cut in cm
  par.lartpc_len = 2.; // new track length cut in cm based on Thomas' study of low energy protons
  par.gastpc_B = 0.4; // B field strength in Tesla
  par.gastpc_padPitch = 0.1; // 1 mm. Actual pad pitch varies, which is going to be impossible to implement
  par.gastpc_X0 = 1300.; // cm = 13m radiation length

  int i = 0;
  while( i < argc ) {
    if( argv[i] == std::string("--edepfile") ) {
      edepfile = argv[i+1];
      i += 2;
    } else if( argv[i] == std::string("--ghepdir") ) {
      ghepdir = argv[i+1];
      i += 2;
    } else if( argv[i] == std::string("--outfile") ) {
      outfile = argv[i+1];
      i += 2;
    } else if( argv[i] == std::string("--fhicl") ) {
      fhicl_filename = argv[i+1];
      i += 2;
    } else if( argv[i] == std::string("--seed") ) {
      par.seed = atoi(argv[i+1]);
      par.run = par.seed;
      i += 2;
    } else if( argv[i] == std::string("--nevents") ) {
      par.n = atoi(argv[i+1]);
      i += 2;
    } else if( argv[i] == std::string("--nfiles") ) {
      par.nfiles = atoi(argv[i+1]);
      i += 2;
    } else if( argv[i] == std::string("--first") ) {
      par.first = atoi(argv[i+1]);
      i += 2;
    } else if( argv[i] == std::string("--oa") ) {
      par.OA_xcoord = atof(argv[i+1]);
      i += 2;
    } else if( argv[i] == std::string("--grid") ) {
      par.grid = true;
      i += 1;
    } else if( argv[i] == std::string("--rhc") ) {
      par.fhc = false;
      i += 1;
    } else if( argv[i] == std::string("--gastpc") ) {
      par.IsGasTPC = true;
      i += 1;
    } else i += 1; // look for next thing
  }

  printf( "Making CAF from edep-sim tree dump: %s\n", edepfile.c_str() );
  printf( "Searching for GENIE ghep files here: %s\n", ghepdir.c_str() );
  if( par.fhc ) printf( "Running neutrino mode (FHC)\n" );
  else printf( "Running antineutrino mode (RHC)\n" );
  if( par.grid ) printf( "Running grid mode\n" );
  else printf( "Running test mode\n" );
  printf( "Output CAF file: %s\n", outfile.c_str() );
  if( par.IsGasTPC ) printf( "Running gas TPC\n" );

  rando = new TRandom3( par.seed );
  CAF caf( outfile, par.IsGasTPC );


  TFile * tf = new TFile( edepfile.c_str() );
  TTree * tree = (TTree*) tf->Get( "tree" );

  loop( caf, par, tree, ghepdir, fhicl_filename );

  caf.version = 4;
  printf( "Run %d POT %g\n", caf.meta_run, caf.pot );
  caf.fillPOT();
  caf.write();

  caf.cafFile->Close();

  printf( "-30-\n" );


}

