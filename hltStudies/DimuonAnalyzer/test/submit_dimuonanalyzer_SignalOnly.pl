#!/opt/star/bin/perl
#use strict;
#use warnings;

############
# Submit jobs on lxplus
#
# This scripts runs in $RUN_DIR       = "/afs/cern.ch/user/s/silvest/scratch0/CMSSW_3_5_0_pre3/src/my_code/DimuonAnalyzer/test"
# the folowwing cfg : $RUN_CFG       = "$RUN_DIR/batch_dimuonanalyzer_SignalOnly_cfg.py";
# which has the path of the input file writing in hard (cause I didnt find out who to write a whole path into a file replacing a string with 
# perl ( rfio:///castor/cern.ch/user/s/silvest/rootfiles/SignalOnly/PARTICLE/root/RECO/PARTICLE_SignalOnly_HIofflineReco_job1_RECO.root)
#
# Make sure to change the "TO BE ADAPTED" and create the related directories in castor and locally
#
# arguments  : perl submier_dimuonanalyzer_SignalOnly.pl particle_type nEvents acceptance_double muon_tag is_submit
#
############

## Checks to make sure your are using the script correctly
# particle
@possible=qw(upsilon jpsi z);
$isin=0;
foreach $possible (@possible) {
    if ( $possible eq $ARGV[0]) {$isin=1;}
}
($isin) || die "First argument must be one of (@possible)\n";
$what_particlue = "$ARGV[0]";# upsilon, jpsi, z

# number of arguments
if($#ARGV ==0 ||$#ARGV ==1|| $#ARGV >4){die "usage submier_dimuonanalyzer_SignalOnly.pl type-to-submit(@possible) nEvents acceptance[double] muon_tag submit\n$!";}

# is submit
$submit=$ARGV[4];
if($submit=~/submit/){$submit=true;}
else {$submit=false;} # by default false

# acceptance
$acceptance = $ARGV[2];# upsilon, jpsi, z
if($acceptance=~ /.8/){ $barrel="barrel";}
elsif($acceptance=~/2./){ $barrel="endcaps";}
else { die "pick a better acceptance or change the [endname]\n$!";}
@possibleTag=qw(STA GLB);
$isin=0;
foreach $possibleTag (@possibleTag) {
    if ( $possibleTag eq $ARGV[3]) {$isin=1;}
}
($isin) || die "Last argument must be one of (@possibleTag)\n";

# muon tag
$tag = "$ARGV[3]";
if($tag=~m/STA/){$muon_tag = "standAloneMuons:UpdatedAtVtx";}
elsif($tag=~m/GLB/){$muon_tag ="globalMuons";}
else{ die "pick a better muon tag or change the [endname]\n$!";}

print "submit $submit\n";

# number of events
$NUM_EVENTS_PER_JOB = $ARGV[1];


$nJobs = 1;
$starting_id = 1;

$pt_min_dimuon = 0;
$pt_min_muon = 0;
$pt_min_track = 0;

$eta_min_muon = -$acceptance;
$eta_max_muon = $acceptance;

$eta_min_track = -$acceptance;
$eta_max_track = $acceptance;

if($what_particlue=~ m/upsilon/)
  { 
    print "upsilon\n";
    $particle = "Upsilon";
    $pdgid = 553;
    $mass_max_dimuon = 14;
    $mass_min_dimuon=4;

    $outRootDir    = "/castor/cern.ch/user/s/silvest/rootfiles/SignalOnly/Upsilon";
    $BASE_NAME="DimuonAnalyzer_Upsilon_SignalOnly_2000events" .$end_name;

  }
if($what_particlue=~ m/jpsi/)
  {
    print "jpsi\n";
    $particle = "Jpsi";
    $pdgid = 443;
    $mass_max_dimuon = 5;
    $mass_min_dimuon=1;
    $outRootDir    = "/castor/cern.ch/user/s/silvest/rootfiles/SignalOnly/Jpsi";
    $BASE_NAME="DimuonAnalyzer_Jpsi_SignalOnly_1400events" .$end_name;
  }

if($what_particlue=~ m/z/)
  { 
    print "Z0\n";
    $particle = "Z";
    $pdgid = 23;
    $mass_max_dimuon = 120;
    $mass_min_dimuon=60;

    $BASE_NAME="DimuonAnalyzer_Z_SignalOnly_1000events" .$end_name;    
    $outRootDir    = "/castor/cern.ch/user/s/silvest/rootfiles/SignalOnly/Z";
  }


## TO BE ADAPTED
## naming files and output directory
# this name appends to the base name
$end_name = "_MassCut_".$barrel."_".$tag;
# my running directory
$RUN_DIR       = "/afs/cern.ch/user/s/silvest/scratch0/CMSSW_3_5_0_pre3/src/my_code/DimuonAnalyzer/test";
# the cfg
$RUN_CFG       = "$RUN_DIR/batch_dimuonanalyzer_SignalOnly_cfg.py";
# output directory that need to be created
$OUTPUT_LOG_DIR           = $outRootDir.'/log/';
$OUTPUT_ROOT_DIR          = $outRootDir.'/root/ANA/';
$OUTPUT_CFG_DIR          = $outRootDir.'/cfg/';
# final name of files    
$FINAL_PARAMETER_SET_NAME = $BASE_NAME;
$CFG                      = $RUN_DIR.'/cfg/'.$FINAL_PARAMETER_SET_NAME.'_cfg.py';
$OUTPUT_LOG               = $FINAL_PARAMETER_SET_NAME.'.log';
$OUTPUT_ROOTFILE          = $FINAL_PARAMETER_SET_NAME.'.root';
$script                   = 'sh/'.$FINAL_PARAMETER_SET_NAME.'.sh';
$skipEvts                 = 0;#($ijob-1)*$NUM_EVENTS_PER_JOB;#change this to 0 because not skipping any events

## Change in the cfg some TAGS to the values
system("cat $RUN_CFG | sed -e 's/CONDOR_MAXEVENTS/$NUM_EVENTS_PER_JOB/g' | sed -e 's/OUTPUTFILE/$OUTPUT_ROOTFILE/g'|  sed -e 's/CONDOR_SKIPEVENTS/$skipEvts/g' | sed -e 's/PARTICLE/$particle/g'| \
sed -e 's/MASS_MAX_DIMUON/$mass_max_dimuon/g'|\
 sed -e 's/MASS_MIN_DIMUON/$mass_min_dimuon/g'|\
 sed -e 's/PT_MIN_DIMUON/$pt_min_dimuon/g'|\
 sed -e 's/PDGID/$pdgid/g'|\
 sed -e 's/ETA_MIN_MU/$eta_min_muon/g'|\
 sed -e 's/ETA_MAX_MU/$eta_max_muon/g'|\
 sed -e 's/PT_MIN_MU/$pt_min_muon/g'|\
 sed -e 's/ETA_MAX_TRACK/$eta_max_track/g'|\
 sed -e 's/ETA_MIN_TRACK/$eta_min_track/g'|\
 sed -e 's/ETA_MIN_MU/$eta_min_muon/g'|\
 sed -e 's/MUON_TAG/$muon_tag/g'|\
 sed -e 's/PT_MIN_TRACK/$pt_min_track/g' > $CFG ");
print "writing $CFG\n";

## Write the .sh bash script (is will end up locally in your ./sh/ directory)
open (SUBMIT, ">$script") || die "Couldn't open $script";
print SUBMIT "#!/bin/bash
#BSUB -q 1nd

echo \"script executing cmssw\"
cd $RUN_DIR
eval `scramv1 runtime -sh`

cd -
echo \"cmsRun $CFG >> $OUTPUT_LOG 2>&1\"
cmsRun $CFG >> $OUTPUT_LOG 2>&1
echo \"script done executing\"
ls
rfcp $OUTPUT_LOG $OUTPUT_LOG_DIR
cp $OUTPUT_LOG $RUN_DIR/log/
rfcp $CFG $OUTPUT_CFG_DIR
rfcp $OUTPUT_ROOTFILE $OUTPUT_ROOT_DIR
cp $OUTPUT_ROOTFILE  $RUN_DIR/root/
rm $OUTPUT_LOG
rm $OUTPUT_ROOTFILE
";
close (SUBMIT);
chmod 0755, $script;

## submit jobs if true
if($submit=~ m/true/) {
  print "bsub -q 1nd -J $BASE_NAME < $script\n";
  print `bsub -q 1nd -J $BASE_NAME < $script`;
}

exit 0;





 
