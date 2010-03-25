#!/opt/star/bin/perl -w
#
use strict;
use integer;
use Switch; 
#chose your settings:
my $dotype       = 0; # 0 for both, 1 signal only, 2 for mixed
my $dodimuon     = 2; # 0 for all, 1 for jpsi, 2 for ups, 3 for z0

my $docentrality = 0; # 0 both, 1 minbias, 2 for central // or for different pp tracking
my $docoverage   = 0; # default full, 2 barrel
#=========================================================================

my @centrality   = ("", "minbias","central");
my @datatype     = ("",     "sgn",    "mix");
my @tracking     = ("",   "pptrk",   "hitrk");
my @dimuon       = ("",     "jpsi",   "ups","z0");

my @coverage     = ("","full","barrel");


my $workDir       = "/afs/cern.ch/user/m/mironov/scratch0/CMSSW_3_4_0/src/Correlation/DijetCorrelation/test/lxjobs/shop";
my $inDir          = "/castor/cern.ch/user/m/mironov/offoutput/matching2/";

#my $workDir        = "/Users/eusmartass/Downloads/L3Switches/test/jobs/macros";
my $outDir         = "$workDir/macros/";

#----------------
# --- centrality
my $icstart = 1;
my $icend   = 3;
switch ($docentrality) {
  case (1) { $icend = 2; }
  case (2) { $icstart = 2; }
  else { print "Running both central and minbias! \n";}
}

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


for ( my $idimu=$idimustart; $idimu<$idimuend; $idimu++)
  {
    for( my $itype=$itypestart; $itype < $itypeend; $itype++)
      {
	for(my $ic=$icstart; $ic<$icend;$ic++)
	  {
	   
		   
			my $dimu = $dimuon[$idimu];
			my $type = $datatype[$itype];

			my $centrking= $centrality[$ic];
			if($itype==1) {$centrking= $tracking[$ic];}
				
		
			my $JOBOUTPUTFILESBASE = $dimu."_".$type."_offline"."_".$centrking;
			my $JOBOUTPUTFILES     = $inDir.$JOBOUTPUTFILESBASE.".root" ;
			print "Moving file ... $JOBOUTPUTFILESBASE \n";
			
			system("rfcp $JOBOUTPUTFILES .");
		    
	  }#centrality
      }# itype
  }# dimu type



