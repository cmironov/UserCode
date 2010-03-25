#!/opt/star/bin/perl -w
#
use strict;
use integer;
use Switch; 
#chose your settings:
my $move      = 1; # 0 for moving out, 1 for moving in


my $dotype   = 0; # 0 for both, 1 signal only, 2 for mixed
my $dodimuon = 3; # 0 for all, 1 for jpsi, 2 for ups, 3 for z0
my $doregion = 2; # 0 for both, 1 for |eta|<0.8, 2 for |eta|<2.4
my $doalgo   = 0; # 0 for all; 1 for ....


my @algos         = ("",
		     "std", "stdpp",
		     "oihit","oistate","oicombo", 
		     "cascade");
my @datatype      = ("","sgn","mix");
my @dimuon        = ("","jpsi","ups","z0");
my @eta           = ("","barrel","full");


my $CMSSW_BASE = "/afs/cern.ch/user/m/mironov/scratch0/CMSSW_3_4_0/src/UserCode/L3Switches";
#my $CMSSW_BASE = "/Users/eusmartass/Downloads/L3Switches";
my $RUNDIR     = "$CMSSW_BASE/test/jobs/macros";

my $indir          = "rfio:///castor/cern.ch/user/m/mironov/hltoutput/step3root/";
my $outdir         = "/tmp/mironov/step3root/";
#----------------

 
# --- set your dimuon
my $idimustart = 1;
my $idimuend   = 4;
switch ($dodimuon) {
  case (1) { $idimuend = 2; }
  case (2) { $idimustart = 2; $idimuend=3; }
  case (3) { $idimustart = 3; }
  else { print "Running for all dimuon type! \n";}
}
# --- mix or signal
my $itypestart = 1;
my $itypeend   = 3;
switch ($dotype) {
  case (1) { $itypeend = 2; }
  case (2) { $itypestart = 2; }
  else { print "Running both signal and mixing! \n";}
}
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


my $ietastart = 1;
my $ietaend   = 3;
switch ($doregion) {
  case (1) { $ietaend = 2; }
  case (2) { $ietastart = 2; }
  else { print "Running both barrel and full ! \n";}
}


if($move==0) 
  {
    print "Removing files from $outdir !!! \n";
    system("rm $outdir/*");
  }
else 
  {
    print "Copying files to $outdir \n";
    for ( my $idimu=$idimustart; $idimu<$idimuend; $idimu++)
      {
	for( my $itype=$itypestart; $itype < $itypeend; $itype++)
	  {
	    for (my $ialgo=$ialgostart; $ialgo<$ialgoend; $ialgo++)
	      {
		for (my $ieta=$ietastart;$ieta<$ietaend;$ieta++)
		  {
		    
		    my $dimu = $dimuon[$idimu];
		    my $type = $datatype[$itype];
		    my $algo = $algos[$ialgo];
		    my $cov  = $eta[$ieta];
	
		    my $rootfilein   = $indir."DQM_V0001_R000000001__".$dimu."_".$type."__".$algo."__".$cov.".root";
		    my $script       = $RUNDIR.'/'.$dimu."_".$type."_".$algo."_".$cov.".sh";
		    
		    print "Copying file $rootfilein to $outdir \n";
		    system("rfcp $rootfilein $outdir");
		  }
	      }
	  }
      }
  }# moving in


