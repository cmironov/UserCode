#!/opt/star/bin/perl

# Submit jobs on lxplus
use strict;
use integer;
use Switch; 

# chose something
my $dotype   = 2; # 0 for both, 1 signal only, 2 for mixed
my $dob      = 1; # 0 for both, 1 for mix_minbias, 2 for mix_central

my $dodimuon = 1; # 0 for all, 1 for jpsi, 2 for ups, 3 for z0
my $doalgo   = 1; # 0 all algos, 1 pp, 2 ppst , etc

my $JOBMAXEVENTS = -1;

my $CMSSW_BASE    = "/afs/cern.ch/user/m/mironov/scratch0/CMSSW_3_4_0/src/Correlation/DijetCorrelation";
#my $CMSSW_BASE = "/Users/eusmartass/Software/DijetCorrelation";
my $RUNDIR        = "$CMSSW_BASE/test/lxjobs/shop";
my $MIXRUNCFG     = "$RUNDIR/mixMatchHLT_cfg.py";
my $SGNRUNCFG     = "$RUNDIR/sgnMatchHLT_cfg.py";

my $OUTFILEDIR ="/castor/cern.ch/user/m/mironov/hltoutput/matching2";

# the 2 eta coverages
my @datatype     = ("", "sgn","mix");

my @dimuon       = ("", "jpsi", "ups", "z0");
my @pdgid         = ("", "443",  "553", "23");

my @centrality   = ("", "minbias","central");
my @algos         = ("",
		     "std", "stdpp",
		     "oihit","oistate","oicombo", 
		     "cascade");

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
# ---- dimu type
my $idimustart = 1;
my $idimuend   = 4;
switch ($dodimuon) {
  case (1) { $idimuend = 2; }
  case (2) { $idimustart = 2; $idimuend=3; }
  case (3) { $idimustart = 3; }
  else { print "Running for all dimuon type! \n";}
}

# ---- centrality
my $ibstart   = 1;
my $ibend     = 3;
switch ($dob) {
  case (1) { $ibend  = 2; }
  case (2) { $ibstart = 2;  }
  else { print "Running both minbias and central for mixing, or pptrk and hitrk for signal! \n";}
}

# ... sgn or mixing
my $itypestart = 1;
my $itypeend   = 3;
switch ($dotype) {
  case (1) { $itypeend = 2; }
  case (2) { $itypestart = 2; }
  else { print "Running both signal and mixing! \n";}
}
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
for( my $itype=$itypestart; $itype < $itypeend; $itype++)  
  {
    for ( my $idimu=$idimustart; $idimu<$idimuend; $idimu++)
      {
	for (my $ib=$ibstart; $ib<$ibend; $ib++)
	  {
	    for (my $ialgo=$ialgostart; $ialgo<$ialgoend; $ialgo++)
	      {
		my $algo       = $algos[$ialgo];
		my $dimu       = $dimuon[$idimu];
		
		my $pdgpair    = $pdgid[$idimu];
		my $type       = $datatype[$itype];
		my $centrking  = $centrality[$ib];
			
		my $JOBINPUTFILESBASE  = $dimu."_".$type."_".$algo;
		my $JOBINPUTFILES      = $JOBINPUTFILESBASE."_full_cfi";
		
		my $JOBOUTPUTFILES     = $JOBINPUTFILESBASE."_".$centrking.".root";
		if( $itype==1 ) {$JOBOUTPUTFILES     = $JOBINPUTFILESBASE."_pp.root";}
		print "Input: $JOBINPUTFILES; Output: $JOBOUTPUTFILES  \n ";

		my $CFG                = $RUNDIR.'/'.$JOBINPUTFILESBASE.'_cfg.py';
		my $LOG                = $JOBINPUTFILESBASE.'.log';
		my $script             = $RUNDIR.'/'.$JOBINPUTFILESBASE.'.sh';
			
		# make unique configuration file for each job, starting	
		my $RUNCFG = $SGNRUNCFG;
		if ($itype == 2) {$RUNCFG = $MIXRUNCFG;}
		
		print "input $JOBINPUTFILES and output $JOBOUTPUTFILES \n ";
		# make unique configuration file for each job, starting	
		system("cat $RUNCFG | sed -e 's#JOB_INPUTFILES#$JOBINPUTFILES#' | sed -e s#JOB_OUTPUTFILENAME#$JOBOUTPUTFILES# | sed -e s#PDGPAIR#$pdgpair# | sed -e s#JOB_MAXEVENTS#$JOBMAXEVENTS# > $CFG ");	
		
		
	 	system("cat > $script <<END
   !/bin/bash
   #BSUB -q 1nd
   # What to call this job
   #BSUB -J matching_mix

   echo \"script executing cmssw\"
   cd $RUNDIR
   eval `scramv1 runtime -sh`
   cd -

   cmsRun $CFG >> $LOG 2>&1
   echo \"script done executing\"
   
   rfcp $JOBOUTPUTFILES $rootDir
   rfcp $LOG $RUNDIR
   rfcp $CFG $pyDir
 
   rm *

   rm $CFG
   rm $script

  END
       ");
 	    system("chmod u+x $script");
 	    print `bsub < $script`;
	    
	      }#algorithm
	  }#centrality
      }# dimu loop      
  }#sgn or bck

exit 0;



