#!/opt/star/bin/perl

# Submit jobs on lxplus
use strict;
use integer;
use Switch; 

# chose something
my $dotype    = 1; # 0 for both, 1 signal only, 2 for mixed
my $dob       = 0; # 0 for both, 1 for minbias, 2 for central
my $dodimuon  = 0; # 0 for all, 1 for jpsi, 2 for ups, 3 for z0
my $tracktype = 0; # 0 both, 1 for global, 2 for sta
my $doalgo   = 0; # 0 all algos, 1 pp, 2 ppst , etc

# chose coverage
my $dobarrel = "True";
my $doendcap = "True";

# pick cuts:
my $ptcutsingle = 0.;
my $ptcutpair   = 0.;

my @massmin_sta  = ("", "2.0",  "7.0",  "60.0");
my @massmax_sta  = ("", "4.5", "12.0", "120.0");

my @massmin_glb  = ("", "2.9", "9.2", "82.0");
my @massmax_glb  = ("", "3.3", "9.8", "100.0");

my $JOBMAXEVENTS = -1;

my $CMSSW_BASE    = "/afs/cern.ch/user/m/mironov/scratch0/CMSSW_3_4_0/src/Correlation/DijetCorrelation";
#my $CMSSW_BASE = "/Users/eusmartass/Software/DijetCorrelation";
my $RUNDIR        = "$CMSSW_BASE/test/lxjobs/shop";
my $RUNCFG        = "$RUNDIR/trackplotter_cfg.py";

my $INFILEDIR  = "/castor/cern.ch/user/m/mironov/hltoutput/matching2/root/";
my $OUTFILEDIR = "/castor/cern.ch/user/m/mironov/hltoutput/plotting";

#-----------------------
my $logDir = "$OUTFILEDIR/log";
if(! -e $logDir){ system("rfmkdir $logDir"); }

my $pyDir = "$OUTFILEDIR/python";
if(! -e $pyDir){ system("rfmkdir $pyDir"); }

my $shDir = "$OUTFILEDIR/sh";
if(! -e $shDir){ system("rfmkdir $shDir"); }

my $rootDir = "$OUTFILEDIR/root";
if(! -e $rootDir){ system("rfmkdir $rootDir"); }

#-----------------------

my @datatype        = ("", "sgn","mix");

my @dimuon           = ("", "jpsi", "ups", "z0");
my @pdgid            = ("", "443",  "553", "23");

my @centrality      = ("", "minbias","central");


my @mytrktype       = ("", "Type1", "Type2");
my @realtrktype     = ("", "HltL3",   "HltL2Upd");

my @algos         = ("",
		     "std", "stdpp",
		     "oihit","oistate","oicombo", 
		     "cascade");

my $idimustart = 1;
my $idimuend   = 4;
switch ($dodimuon) {
  case (1) { $idimuend   = 2; }
  case (2) { $idimustart = 2; $idimuend=3; }
  case (3) { $idimustart = 3; }
  else { print "Running for all dimuon type! \n";}
}

# centrlaity or tracking type
my $ibstart   = 1;
my $ibend     = 3;
switch ($dob) {
  case (1) { $ibend  = 2;   }
  case (2) { $ibstart = 2;  }
  else { print "Running both minbias and centra/pptracking and hi tracking (if itype=sgn)l! \n";}
}

# ---- L3 or L2
my $itrkstart = 1;
my $itrkend   = 3;
switch ($tracktype) {
  case (1) { $itrkend = 2;   }
  case (2) { $itrkstart = 2; }
  else { print "Running both global and sta muons! \n";}
}

# ---- sgn or mixing
my $itypestart = 1;
my $itypeend   = 3;
switch ($dotype) {
  case (1) { $itypeend = 2; }
  case (2) { $itypestart = 2; }
  else { print "Running both signal and mixing! \n";}
}

# --- algo choice
my $ialgostart = 1;
my $ialgoend   = 7;
switch ($doalgo) {
  case (1) { $ialgoend=2; }
  case (2) { $ialgostart=2; $ialgoend=3; }
  case (3) { $ialgostart=3; $ialgoend=4; }
  case (4) { $ialgostart=4; $ialgoend=5; }
  case (5) { $ialgostart=5; $ialgoend=6; }
  case (6) { $ialgostart=6;  }
  else { print "Running all algos! \n";}
}

my $cov = "full";
if( ($dobarrel eq "True") && ($doendcap ne "True") )  {$cov = "barrel";}
if( ($dobarrel ne "True") && ($doendcap eq "True") )  {$cov = "endcap";}

my $centrking;
for( my $itype=$itypestart; $itype < $itypeend; $itype++)  
  {
    for ( my $idimu=$idimustart; $idimu<$idimuend; $idimu++)
      {
	for (my $ib=$ibstart; $ib<$ibend; $ib++)
	  {
	    for ( my $itrk=$itrkstart; $itrk<$itrkend; $itrk++ )
	      {
		for (my $ialgo=$ialgostart; $ialgo<$ialgoend; $ialgo++)
		  {
		    my $dimu       = $dimuon[$idimu];
		    my $algo       = $algos[$ialgo];
		    my $pdgpair    = $pdgid[$idimu];
		    my $type       = $datatype[1];
		    
		    my $massmax    = $massmax_glb[$idimu];
		    my $massmin    = $massmin_glb[$idimu];
		    if($itrk==2)
		      {
			$massmax    = $massmax_sta[$idimu];
			$massmin    = $massmin_sta[$idimu];
		      }
		    my $typenumber = $mytrktype[$itrk];
		    my $realtrk    = $realtrktype[$itrk];
		    
		    if($itype==1) {$centrking= "pp";}
		    else {$centrking= $centrality[$ib];}
		    
		    my $JOBINPUTFILESBASE  = $dimu."_".$type."_".$algo."_".$centrking;
		    # my $JOBINPUTFILES      = "file".$INFILEDIR.$JOBINPUTFILESBASE.".root";
		    my $JOBINPUTFILES      = $RUNDIR."/".$JOBINPUTFILESBASE.".root";
		    
		    my $JOBOUTPUTFILESBASE = $JOBINPUTFILESBASE."_".$realtrk."_".$cov;
		    my $JOBOUTPUTFILES     = $JOBOUTPUTFILESBASE."_hists.root" ;
		    
		    my $LOG                = $JOBOUTPUTFILESBASE.'.log';
		    my $CFG                = $RUNDIR.'/'.$JOBOUTPUTFILESBASE.'_cfg.py';
		    my $script             = $RUNDIR.'/'.$JOBOUTPUTFILESBASE.'.sh';
			    
		    print "input $JOBINPUTFILES and output $JOBOUTPUTFILES \n ";
		    # make unique configuration file for each job, starting	
		#     system("cat $RUNCFG | sed -e 's#JOB_INPUTFILENAME#$JOBINPUTFILES#' |\
#  sed -e s#JOB_OUTPUTFILENAME#$JOBOUTPUTFILES# |\
#  sed -e s#TYPE#$typenumber# |\
#  sed -e s#DOBARREL#$dobarrel# |\
#  sed -e s#DOENDCAP#$doendcap# |\
#  sed -e 's#PDGPAIR#$pdgpair#' |\
#  sed -e 's#PTCUTSINGLE#$ptcutsingle#' |\
#  sed -e 's#PTCUTPAIR#$ptcutpair#' |\
#  sed -e s#MASSMAX#$massmax# |\
#  sed -e s#MASSMIN#$massmin# |\
#  sed -e s#JOB_MAXEVENTS#$JOBMAXEVENTS# > $CFG ");	
		   
		   
# 		   system("cat > $script << END
#    !/bin/bash
#    #BSUB -q 8nm
#    # What to call this job
#    #BSUB -J plotting_mix

#    echo \"script executing cmssw\"
#    cd $RUNDIR
#    eval `scramv1 runtime -sh`
#    cd -

#    cmsRun $CFG >> $LOG 2>&1
#    echo \"script done executing\"
   
   
#    rfcp $JOBOUTPUTFILES $rootDir
#    rfcp $LOG $logDir
#    rfcp $CFG $pyDir

#    rm *

#    rm $CFG
#    rm $script
# END
#        ");
	    
	    
#  	    system("chmod u+x $script");
#  	    print `bsub < $script`;
		  }# algo type
		 }#track type loop
	      }#centrality sau tracking
	  }# dimu loop      
  }#singal or mixing
exit 0;



