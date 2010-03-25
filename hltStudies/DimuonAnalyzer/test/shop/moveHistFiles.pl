#!/opt/star/bin/perl -w
#
use strict;
use integer;
use Switch; 
#chose your settings:
my $dotype       = 2; # 0 for both, 1 signal only, 2 for mixed
my $dodimuon     = 0; # 0 for all, 1 for jpsi, 2 for ups, 3 for z0

my $docentrality = 0; # 0 both, 1 minbias, 2 for central // or for different pp tracking
my $dohisto      = 0; # 0 both, 1 for efficiency, 2 for fake
my $dobody       = 0; # 0 for both, 1 for single, 2 for pair

my $domuon       = 0; # 0 for both, 1 global, 2 for sta

my $docoverage   = 0; # default full, 2 barrel
#=========================================================================

my @centrality   = ("", "minbias","central");
my @datatype     = ("",     "sgn",    "mix");
my @tracking     = ("",   "pptrk",   "hitrk");
my @dimuon       = ("",     "jpsi",   "ups","z0");

my @coverage     = ("","full","barrel");
my @recomuons    = ("","GLB","STA");
my @histotype    = ("","Efficiency","Fake");

my $workDir       = "/afs/cern.ch/user/m/mironov/scratch0/CMSSW_3_4_0/src/Correlation/DijetCorrelation/test/lxjobs/shop";
my $inDir          = "/castor/cern.ch/user/m/mironov/offoutput/plotting/";

#my $workDir        = "/Users/eusmartass/Downloads/L3Switches/test/jobs/macros";
my $outDir         = "$workDir/macros/";

#----------------
my $cov = "full";
if($docoverage==2){$cov="barrel";}

# --- set histo type
my $ihstart = 1;
my $ihend   = 3;
switch ($dohisto) {
  case (1) { $ihend = 2; }
  case (2) { $ihstart = 2; }
  else { print "Running both efficiency and fake! \n";}
}

# --- set muon
my $imustart = 1;
my $imuend   = 3;
switch ($domuon) {
  case (1) { $imuend = 2; }
  case (2) { $imustart = 2; }
  else { print "Running both GLB and STA! \n";}
}

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

my $ibodystart = 1;
my $ibodyend   = 3;
switch ($dobody) {
  case (1) { $ibodyend = 2; }
  case (2) { $ibodystart = 2; }
  else { print "Running both single and dim- muon! \n";}
}


for ( my $idimu=$idimustart; $idimu<$idimuend; $idimu++)
  {
    for( my $itype=$itypestart; $itype < $itypeend; $itype++)
      {
	for(my $ic=$icstart; $ic<$icend;$ic++)
	  {
	    for(my $imu=$imustart; $imu<$imuend;$imu++)
	      {
		for(my $ih=$ihstart;$ih<$ihend;$ih++)
		  {
		    for(my $ibody=$ibodystart;$ibody<$ibodyend;$ibody++)
		      {
			my $dimu = $dimuon[$idimu];
			my $type = $datatype[$itype];

			my $centrking= $centrality[$ic];
			if($itype==1) {$centrking= $tracking[$ic];}
				
			my $recomu= $recomuons[$imu];
			my $histo = $histotype[$ih];
			my $body=$ibody;
		

			my $JOBOUTPUTFILESBASE = $dimu."_".$type."_offline"."_".$centrking."_".$recomu."_".$cov;
			my $JOBOUTPUTFILES     = $inDir.$JOBOUTPUTFILESBASE."_hists.root" ;
			print "Moving file ... $JOBOUTPUTFILESBASE \n";
			
			system("rfcp $JOBOUTPUTFILES .");
		      }# single or di-muon
		  }#histogram type
	      }#type of reco muon
	  }#centrality
      }# itype
  }# dimu type



