import os
import os.path
import ROOT

from optparse import OptionParser

from array import array


gar_active_vols = ["GArTPC", "TPCChamber", "TPCGas", "cent_elec_shape", "cent_hc_shape", "TPC1_shape", "TPC1pad_shape", "TPC1fc_pvf_shape", "TPC1fc_kev_shape", "TPC2_shape", "TPC2pad_shape", "TPC2fc_pvf_shape", "TPC2fc_kev_shape", "TPC2fc_hc_shape"]

fname = "/pnfs/dune/persistent/users/mtanaz/edep/GAr/FHC/1/GAr.neutrino.0.edepsim.root"
tf = ROOT.TFile( fname, "OLD" )

tf.MakeProject( "EDepSimEvents", "*", "RECREATE++" )

events = ROOT.TChain( "EDepSimEvents", "main event tree" )
events.Add( fname )
event = ROOT.TG4Event()
events.SetBranchAddress("Event", ROOT.AddressOf(event))

ROOT.gROOT.SetBatch(1)

parser = OptionParser()
parser.add_option('--outfile', help='Output file name', default="gartpc_3_root_out.root")

(args, dummy) = parser.parse_args()
# a number of fiducial cuts in GAr
offset = [ 0., 300., 60. ]
fvLo = [ -150., -100., 50. ]
fvHi = [ 150., 100., 350. ]
collarLo = [ -170., -120., 30. ]
collarHi = [ 170., 120., 470. ]
# tree entries being created and populated
fout = ROOT.TFile( args.outfile, "RECREATE" )
tout = ROOT.TTree( "tree","tree" )
t_ifileNo = array('i',[0])
tout.Branch('ifileNo',t_ifileNo,'ifileNo/I')
t_ievt = array('i',[0])
tout.Branch('ievt',t_ievt,'ievt/I')
t_p3lep = array('f',3*[0.0])
tout.Branch('p3lep',t_p3lep,'p3lep[3]/F')
t_vtx = array('f',3*[0.0])
tout.Branch('vtx',t_vtx,'vtx[3]/F')
t_lepPdg = array('i',[0])
tout.Branch('lepPdg',t_lepPdg,'lepPdg/I')
t_lepKE = array('f',[0])
tout.Branch('lepKE',t_lepKE,'lepKE/F')
t_muonExitPt_active = array('f',3*[0.0])
tout.Branch('muonExitPt_active',t_muonExitPt_active,'muonExitPt_active[3]/F')
t_muonExitPt_passive = array('f',3*[0.0])
tout.Branch('muonExitPt_passive',t_muonExitPt_passive,'muonExitPt_passive[3]/F')
t_muonExitMom = array('f',3*[0.0])
tout.Branch('muonExitMom',t_muonExitMom,'muonExitMom[3]/F')
t_muonReco = array('i',[0])
tout.Branch('muonReco',t_muonReco,'muonReco/I')
t_muGArLen = array('f',[0])
tout.Branch('muGArLen',t_muGArLen,'muGArLen/F')
t_hadTot = array('f', [0.] )
tout.Branch('hadTot', t_hadTot, 'hadTot/F' )
t_hadCollar = array('f', [0.] )
tout.Branch('hadCollar', t_hadCollar, 'hadCollar/F' )
t_nFS = array('i',[0])
tout.Branch('nFS',t_nFS,'nFS/I')
t_fsPdg = array('i',100*[0])
tout.Branch('fsPdg',t_fsPdg,'fsPdg[nFS]/I')
t_fsPx = array('f',100*[0.])
tout.Branch('fsPx',t_fsPx,'fsPx[nFS]/F')
t_fsPy = array('f',100*[0.])
tout.Branch('fsPy',t_fsPy,'fsPy[nFS]/F')
t_fsPz = array('f',100*[0.])
tout.Branch('fsPz',t_fsPz,'fsPz[nFS]/F')
t_fsE = array('f',100*[0.])
tout.Branch('fsE',t_fsE,'fsE[nFS]/F')
t_fsTrkLen = array('f',100*[0.])
tout.Branch('fsTrkLen',t_fsTrkLen,'fsTrkLen[nFS]/F')
t_muGArEDep = array('f',[0])
tout.Branch('muGArEDep',t_muGArEDep,'muGArEDep/F')


N = events.GetEntries()
tgeo = None
if tgeo is None: # first OK file, get geometry
            tgeo = tf.Get("EDepSimGeometry")
tf.Close() # done with this one


for ient in range(N):

	events.GetEntry(ient)
	for ivtx,vertex in enumerate(event.Primaries): 
            t_vtx[0]=0.0; t_vtx[1]=0.0; t_vtx[2]=0.0;
	    t_p3lep[0]=0.0; t_p3lep[1]=0.0; t_p3lep[2]=0.0; 
            t_lepPdg[0] = 0
	    t_lepKE[0] = 0.
            t_muonExitPt_active[0] = 0.0; t_muonExitPt_active[1] = 0.0; t_muonExitPt_active[2] = 0.0;
	    t_muonExitPt_passive[0] = 0.0; t_muonExitPt_passive[1] = 0.0; t_muonExitPt_passive[2] = 0.0;
            t_muonExitMom[0] = 0.0; t_muonExitMom[1] = 0.0; t_muonExitMom[2] = 0.0;
            t_muonReco[0] = -1;
            t_muGArLen[0]=0.0;
	    t_muGArEDep[0]=0.0;
            t_hadTot[0] = 0.;
            t_hadCollar[0] = 0.;
            t_nFS[0] = 0;
	    reaction=vertex.Reaction
            for i in range(3):
                t_vtx[i] = vertex.Position[i] / 10. - offset[i] # cm       
            ileptraj = -1
            nfsp = 0
            fsParticleIdx = {}
 	    for ipart,particle in enumerate(vertex.Particles):
		e = particle.Momentum[3]
                p = (particle.Momentum[0]**2 + particle.Momentum[1]**2 + particle.Momentum[2]**2)**0.5
                m = (e**2 - p**2)**0.5
                t_fsPdg[nfsp] = particle.PDGCode
                t_fsPx[nfsp] = particle.Momentum[0]
                t_fsPy[nfsp] = particle.Momentum[1]
                t_fsPz[nfsp] = particle.Momentum[2]
                t_fsE[nfsp] = e
                fsParticleIdx[particle.TrackId] = nfsp
                nfsp += 1
                if abs(particle.PDGCode) in [11,12,13,14]:
                    ileptraj = particle.TrackId
                    t_lepPdg[0] = particle.PDGCode
                    for i in range(3): t_p3lep[i] = particle.Momentum[i]
                    t_lepKE[0] = e - m
            assert ileptraj != -1, "There isn't a lepton??"
            t_nFS[0] = nfsp
            muexit = 0
            exitKE = 0.
            exitP = None
            endVolIdx = -1 # where does the muon die
            if abs(t_lepPdg[0]) == 13:
		leptraj = event.Trajectories[ileptraj]

		for p in leptraj.Points:
                    pt = p.Position
                    node = tgeo.FindNode( pt.X(), pt.Y(), pt.Z() )
                    volName = node.GetName()
                    for v in gar_active_vols:
                        if v in volName:
				t_muonExitPt_active[0] = pt.X() / 10. - offset[0]
				t_muonExitPt_active[1] = pt.Y() / 10. - offset[1]
				t_muonExitPt_active[2] = pt.Z() / 10. - offset[2]		
			if v not in volName:
				t_muonExitPt_passive[0] = pt.X() / 10. - offset[0]
        	        	t_muonExitPt_passive[1] = pt.Y() / 10. - offset[1]
                	    	t_muonExitPt_passive[2] = pt.Z() / 10. - offset[2]
                    		t_muonExitMom[0] = p.Momentum.x()
                   	 	t_muonExitMom[1] = p.Momentum.y()
                    		t_muonExitMom[2] = p.Momentum.z()
                    		if abs(pt.X() / 10. - offset[0]) > 200.: muexit = 1 # side exit
                    		elif abs(pt.Y() / 10. - offset[1]) > 150.: muexit = 2 # top/bottom exit
                    		elif pt.Z() / 10. - offset[2] < 0.: muexit = 3 # upstream exit
                    		elif pt.Z() / 10. - offset[2] > 500.: muexit = 4 # downstream exit
                    		else:
                        		print "Hit in %s at position (%1.1f, %1.1f, %1.1f) unknown exit!" % (volName, pt.X()/10.-offset[0], pt.Y()/10.-offset[1], pt.Z()/10.-offset[2])
	                   	exitP = p.Momentum
      		            	exitKE = (exitP.x()**2 + exitP.y()**2 + exitP.z()**2 + 105.3**2)**0.5 - 105.3	
	        	       	endpt = leptraj.Points[-1].Position	
                	node = tgeo.FindNode( endpt.X(), endpt.Y(), endpt.Z() )
                	endVolName = node.GetName()	
			if "volWorld" in endVolName or "volDetEnclosure" in endVolName: endVolIdx = 0 # outside detector components
                	elif "volLArActive" in endVolName or "volPixelPlane" in endVolName: endVolIdx = 1 # active LAr
                	elif "volLAr" in endVolName or "DsPlane" in endVolName or "UsPlane" in endVolName: endVolIdx = 2 # Passive component of LAr
                	elif "volCylinder" in endVolName or "ArgonCube" in endVolName: endVolIdx = 2 # Passive component of LAr
                	elif "TPCChamber" in endVolName: endVolIdx = 6 # use "endcap yoke" idx for pressure vessel
                	elif "TPC" in endVolName: endVolIdx = 3 # very rare active tpc stopper
                	elif "ECALLeft" in endVolName or "ECALRight" in endVolName: endVolIdx = 4 # "endcap" ECALs
                	elif "ECAL" in endVolName or "SB" in endVolName: endVolIdx = 5 # "barrel ECALs
                	elif "Yoke" in endVolName: endVolIdx = 7
                	elif "Mag" in endVolName: endVolIdx = 8
                	hits = []
                	for key in event.SegmentDetectors:
                    		if key.first in ["TPC1", "TPC2"]:
                        		hits += key.second

 			tot_energy_deposit = 0.0
			tot_length = 0.0
                	for hit in hits:
                    		if hit.PrimaryId == ileptraj: # hit is due to the muon


 		                       	hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                	        	hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )
                        		tot_length += (hStop-hStart).Mag()
					tot_energy_deposit += hit.EnergyDeposit
                	t_muGArLen[0] = tot_length
			t_muGArEDep[0] = tot_energy_deposit

                # muon reconstruction method
                # 1 = contained in the avtive volume
                	if endVolIdx == 1:
                    		t_muonReco[0] = 1
                    		if muexit != 0: print "Muon exit %d but end vol is active" % muexit
                # 2 = gas TPC match
                	elif tot_length > 0.:
                    		t_muonReco[0] = 2
                # 3 = ECAL stopping
                	elif endVolIdx == 4 or endVolIdx == 5:
                    		t_muonReco[0] = 3
                # 4 = magnet/coil stopper
                	elif endVolIdx == 7 or endVolIdx == 8:
                    		t_muonReco[0] = 4
                # 5 = passive Ar stopper
                	elif endVolIdx == 2:
                    		t_muonReco[0] = 5
                # 6 = side-exiting
                	elif muexit == 1:
				t_muonReco[0] = 6
                # 7 = top/bottom-exiting
                	elif muexit == 2:
                    		t_muonReco[0] = 7
                # 8 = upstream-exiting
                	elif muexit == 3:
                    		t_muonReco[0] = 8

	    hits = []
            for key in event.SegmentDetectors:
                if key.first == "GArTPC":
                    hits += key.second

            collar_energy = 0.
            total_energy = 0.
            track_length = [0. for i in range(nfsp)]
            for hit in hits:
                hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )
                if event.Trajectories[hit.PrimaryId].ParentId == -1:
                    track_length[fsParticleIdx[hit.PrimaryId]] += (hStop-hStart).Mag()
	
                if hit.PrimaryId != ileptraj:
                    hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                    total_energy += hit.EnergyDeposit
                    # check if hit is in collar region
                    if hStart.x() < collarLo[0] or hStart.x() > collarHi[0] or hStart.y() < collarLo[1] or hStart.y() > collarHi[1] or hStart.z() < collarLo[2] or hStart.z() > collarHi[2]:
                        collar_energy += hit.EnergyDeposit

            t_hadTot[0] = total_energy
            t_hadCollar[0] = collar_energy
            for i in range(nfsp):
                t_fsTrkLen[i] = track_length[i]

            tout.Fill()
fout.cd()
tout.Write()
output_1.close()
print "reached the end of file"
