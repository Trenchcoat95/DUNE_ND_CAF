#!/usr/bin/env python

import sys
import os.path
import os
import ROOT
from optparse import OptionParser
from array import array


active = ["TPCGas_vol_PV_0","TPC2_PV_0", "TPC1_PV_0","TPC1fc_pvf_vol_PV_0","TPC2fc_pvf_vol_PV_0","TPC1pad_vol_PV_0","TPC2pad_vol_PV_0"] 


def loop( events, tgeo, tout, nfiles, okruns ):
 
    if len(okruns) == 0:
        #print "There are no runs in this TTree...skipping!"
        return

    #print "Inside event loop with %d files and first run %d" % (nfiles, okruns[0])

    # updated geometry with less steel
    offset = [ 0., 305., 5. ]
    fvLo = [ -300., -100., 50. ]
    fvHi = [ 300., 100., 450. ]
    collarLo = [ -320., -120., 30. ]
    collarHi = [ 320., 120., 470. ]

    event = ROOT.TG4Event()
    events.SetBranchAddress("Event",ROOT.AddressOf(event))
    #print "Set branch address"
 
    N = events.GetEntries()
    evt_per_file = N/nfiles
    #if N % nfiles:
        #print "Files don't all have the same number of events!!!"
        #print "AAAAAAAAAAAAAAAAAAAAAAAAHHHHHHHHHHHHHHHHhhhh that's bad stop"
        #print "\n\n\n\n\n\n\n\n\n\n\n"
        #print "AAAAAAAAAAAAAAAAAAAAAAAAHHHHHHHHHHHHHHHHhhhh seriously stop it"

    #print "Starting loop over %d entries" % N
    for ient in range(N):
        pi_pl_count = 0
        pi_min_count = 0
       
        nue_tot = 0.0

        #if ient % 100 == 0:
        #    print "Event %d" % (ient,N)
        events.GetEntry(ient)
        for ivtx,vertex in enumerate(event.Primaries):
  	    
            fileidx = ient/evt_per_file
            t_ifileNo[0] = okruns[fileidx]
            t_ievt[0] = ient%evt_per_file;
            t_vtx[0]=0.0; t_vtx[1]=0.0; t_vtx[2]=0.0;
	    t_true_numu[0] = 0
	    t_true_nue[0] = 0
	    t_ntrk[0] = 0
	    t_ntrkpoints[0] = 0
	    if 'nu:14' in vertex.Reaction and 'Weak[CC]' in vertex.Reaction: 
	    	#and '1000180400' in vertex.Reaction:
		t_true_numu[0]=1        
	    if 'nu:12' in vertex.Reaction and 'Weak[CC]' in vertex.Reaction: 
		#and '1000180400' in vertex.Reaction:
		t_true_nue[0]=1
	    if 'Weak[NC]' in vertex.Reaction: 
		#and '1000180400' in vertex.Reaction:
		t_true_nunc[0]=1
		#print vertex.Reaction
            for i in range(3): 
                t_vtx[i] = vertex.Position[i] / 10. # cm
          
            nfsp = 0
	    ntrk = 0
            ntrkpoints = 0
            for ipart,particle in enumerate(vertex.Particles):
                e = particle.Momentum[3]
                p = (particle.Momentum[0]**2 + particle.Momentum[1]**2 + particle.Momentum[2]**2)**0.5
                
		
		m = (e**2 - p**2)**0.5
                t_mass[nfsp] = m   
                t_fsPdg[nfsp] = particle.PDGCode
                t_fsPx[nfsp] = particle.Momentum[0]
                t_fsPy[nfsp] = particle.Momentum[1]
                t_fsPz[nfsp] = particle.Momentum[2]
                t_fsE[nfsp] = e
                nfsp += 1
	    for itrack,tracks in enumerate(event.Trajectories):
                t_trkpdg[ntrk] = tracks.PDGCode
		#print tracks.ParentId
		t_parentid[ntrk] = tracks.ParentId
		ntrk += 1
		for itrackpoints,trackpoints in enumerate(tracks.Points):
		    t_trkpositionx[ntrkpoints] = trackpoints.Position.X() / 10.
		    t_trkpositiony[ntrkpoints] = trackpoints.Position.Y() / 10.
		    t_trkpositionz[ntrkpoints] = trackpoints.Position.Z() / 10.
		  
		    t_trkpx[ntrkpoints] = trackpoints.Momentum.x()
                    t_trkpy[ntrkpoints] = trackpoints.Momentum.y() 
                    t_trkpz[ntrkpoints] = trackpoints.Momentum.z() 
		    t_proc[ntrkpoints] = trackpoints.Process
		    t_subproc[ntrkpoints] = trackpoints.Subprocess
		    ntrkpoints += 1
		    #print ntrkpoints
            	    pt = trackpoints.Position
                    startpt = tracks.Points[0].Position 
		    endpt = tracks.Points[-1].Position
                    
		    stnode = tgeo.FindNode( startpt.X(), startpt.Y(), startpt.Z() )
                    volstname = stnode.GetName()

		    endnode = tgeo.FindNode( endpt.X(), endpt.Y(), endpt.Z() )
		    if tracks.PDGCode == 14:
		    	continue
		    
		    trkptot = (trackpoints.Momentum.x()**2 + trackpoints.Momentum.y()**2 + trackpoints.Momentum.z()**2)**0.5 
		    fsptot = (particle.Momentum[0]**2 + particle.Momentum[1]**2 + particle.Momentum[2]**2)**0.5
		    pointtot = (trackpoints.Position.X()**2 + trackpoints.Position.Y()**2 + trackpoints.Position.Z()**2)**0.5 
		    vtxtot = (vertex.Position[0]**2 + vertex.Position[1]**2 + vertex.Position[2]**2)**0.5                  
		    
                    
		    for len_act in range(0,len(active)):
		    	if "Tile" in volendname and active[len_act] in volstname:
                        	t_p_loc[ntrk] = 1
                    for len_act in range(0,len(active)):
		    	if active[len_act] not in volendname and active[len_act] not in volstname:
                        	t_p_loc[ntrk] = 2
                    for len_act in range(0,len(active)):
		    	if active[len_act] not in volstname and "Tile" in volendname: 
		    		t_p_loc[ntrk] = 3
		    for len_act in range(0,len(active)):
		    	if active[len_act] not in volstname and active[len_act] in volendname:
				t_p_loc[ntrk] = 4
                    for len_act in range(0,len(active)):
                        if active[len_act] in volstname and active[len_act] in volendname:
                                t_p_loc[ntrk] = 5
		    #if "MPTYoke" in volstname and "MPTYoke" in volendname:
		      
		
	    t_ntrk[0] = ntrk
            t_nFS[0] = nfsp
            t_ntrkpoints[0] = ntrkpoints
            
	    tout.Fill()


if __name__ == "__main__":

    ROOT.gROOT.SetBatch(1)

    parser = OptionParser()
    parser.add_option('--outfile', help='Output file name', default="out.root")
    parser.add_option('--topdir', help='Input file top directory', default="")
    parser.add_option('--first_run', type=int, help='First run number', default=0)
    parser.add_option('--last_run', type=int, help='Last run number', default=0)
    parser.add_option('--rhc', action='store_true', help='Reverse horn current', default=False)
    parser.add_option('--grid', action='store_true', help='grid mode', default=False)

    (args, dummy) = parser.parse_args()

    # make an output ntuple
    fout = ROOT.TFile( args.outfile, "RECREATE" )
    tout = ROOT.TTree( "tree","tree" )
    
    runnumb = args.first_run
    out = open("true_nue_sample_"+str(runnumb)+".txt","w")
    out_2 = open("non_nue_sample_"+str(runnumb)+".txt","w") 
 


    t_ifileNo = array('i',[0])
    tout.Branch('ifileNo',t_ifileNo,'ifileNo/I')
    t_ievt = array('i',[0])
    tout.Branch('ievt',t_ievt,'ievt/I')
    t_vtx = array('f',3*[0.0])
    tout.Branch('vtx',t_vtx,'vtx[3]/F')
    
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
    t_fsTrkLenPerp = array('f',100*[0.])
    tout.Branch('fsTrkLenPerp',t_fsTrkLenPerp,'fsTrkLenPerp[nFS]/F')

    t_mass = array('f',100*[0.0])
    tout.Branch('mass',t_mass,'mass[nFS]/F')

    t_ntrk = array('i',[0])
    tout.Branch('ntrk',t_ntrk,'ntrk/I')

    t_ntrkpoints = array('i',[0])
    tout.Branch('ntrkpoints',t_ntrkpoints,'ntrkpoints/I')

    t_trkpdg = array('f',10000*[0.])
    tout.Branch('trkpdg',t_trkpdg,'trkpdg[ntrk]/F')
	
    t_parentid = array('f',10000*[0.])
    tout.Branch('parentid',t_parentid,'parentid[ntrk]/F')

    t_trkpositionx = array('f',100000*[0.])
    tout.Branch('trkpositionx',t_trkpositionx,'trkpositionx[ntrkpoints]/F')

    t_trkpositiony = array('f',100000*[0.])
    tout.Branch('trkpositiony',t_trkpositiony,'trkpositiony[ntrkpoints]/F')

    t_trkpositionz = array('f',100000*[0.])
    tout.Branch('trkpositionz',t_trkpositionz,'trkpositionz[ntrkpoints]/F')

    t_trkpx = array('f',100000*[0.])
    tout.Branch('trkpx',t_trkpx,'trkpx[ntrkpoints]/F')

    t_trkpy = array('f',100000*[0.])
    tout.Branch('trkpy',t_trkpy,'trkpy[ntrkpoints]/F')

    t_trkpz = array('f',100000*[0.])
    tout.Branch('trkpz',t_trkpz,'trkpz[ntrkpoints]/F')

    t_proc = array('f',100000*[0.])
    tout.Branch('proc',t_proc,'proc[ntrkpoints]/F')

    t_subproc = array('f',100000*[0.])
    tout.Branch('subproc',t_subproc,'subproc[ntrkpoints]/F')

    t_true_numu = array('i',[0])
    tout.Branch('true_numu',t_true_numu,'true_numu/I')
    
    t_true_nue = array('i',[0])
    tout.Branch('true_nue',t_true_nue,'true_nue/I')

    t_true_nunc = array('i',[0])
    tout.Branch('true_nunc',t_true_nunc,'true_nunc/I')

    t_p_loc = array('i',100000*[0])
    tout.Branch('p_loc',t_p_loc,'t_p_loc[ntrk]/I')

    tgeo = None

    events = ROOT.TChain( "EDepSimEvents", "main event tree" )

    neutrino = "neutrino"
    if args.rhc:
        neutrino = "antineutrino"

    #print "Building TChains for runs %d-%d..." % (args.first_run, args.last_run)
    nfiles = 0
    okruns = []
    for run in range( args.first_run, args.last_run+1 ):
        if args.grid:
            fname = "%s/edep.%d.root" % (args.topdir,run)
        else:
            fname = "%s/GAr.%s.%d.edepsim.root" % (args.topdir, neutrino, run)
        #print fname

        # see if it is an OK file
        if not os.access( fname, os.R_OK ):
            #print "Can't access file: %s" % fname
            continue
        tf = ROOT.TFile( fname )
        if tf.TestBit(ROOT.TFile.kRecovered): # problem with file
            #print "File is crap: %s" % fname
            continue
        nfiles += 1
        okruns.append( run )
	loaded = False
        if not loaded:
            loaded = True
            tf.MakeProject("EDepSimEvents","*","RECREATE++")

        # add it to the tchain
        events.Add( fname )

        if tgeo is None: # first OK file, get geometry
            tgeo = tf.Get("EDepSimGeometry")
        tf.Close() # done with this one

    #print "OK runs: ", sorted(okruns)
    #print "got %d events in %d files = %1.1f events per file" % (events.GetEntries(), nfiles, 1.0*events.GetEntries()/nfiles)
    loop( events, tgeo, tout, nfiles, sorted(okruns) )
    out.close()
    out_2.close()
    fout.cd()
    tout.Write()




