#!/opt/star/bin/perl
#use strict;
#use warnings;

# Submit jobs on lxplus
#use strict;
@possible=qw(upsilon jpsi z);
$isin=0;
foreach $possible (@possible) {
    if ( $possible eq $ARGV[0]) {$isin=1;}
}
($isin) || die "First argument must be one of (@possible)\n";

if($#ARGV ==0 ||$#ARGV ==1|| $#ARGV >4){die "usage submit.pl type-to-submit(@possible) nEvents acceptance[double] muon_tag submit\n$!";}

$acceptance = $ARGV[2];# upsilon, jpsi, z

# cmssw output
$submit=$ARGV[4];
if($submit=~/submit/){$submit=true;}
else {$submit=false;}

if($acceptance=~ /.8/){ $barrel="barrel";}
elsif($acceptance=~/2./){ $barrel="endcaps";}
else { die "pick a better acceptance or change the [endname]\n$!";}
@possibleTag=qw(STA GLB);
$isin=0;
foreach $possibleTag (@possibleTag) {
    if ( $possibleTag eq $ARGV[3]) {$isin=1;}
}
($isin) || die "Last argument must be one of (@possibleTag)\n";

$what_particlue = "$ARGV[0]";# upsilon, jpsi, z
#
$tag = "$ARGV[3]";
if($tag=~m/STA/){$muon_tag = "standAloneMuons:UpdatedAtVtx";}
elsif($tag=~m/GLB/){$muon_tag ="globalMuons";}
else{ die "pick a better muon tag or change the [endname]\n$!";}

print "submit $submit\n";

$NUM_EVENTS_PER_JOB = $ARGV[1];
$end_name = "_MassCut_".$barrel."_".$tag;

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
    $pdgid = 553;
    $mass_max_dimuon = 20;
    $mass_min_dimuon=0;

    $import_source = "Upsilonembedded_RECO_cfi";
    $outRootDir    = "/castor/cern.ch/user/s/silvest/rootfiles/HI/Upsilon";
    $BASE_NAME="DimuonAnalyzer_Upsilon_embedding_500events" .$end_name;

  }
if($what_particlue=~ m/jpsi/)
  {
    print "jpsi\n";
    $pdgid = 443;
    $mass_max_dimuon = 6;
    $mass_min_dimuon=0;
    $outRootDir    = "/castor/cern.ch/user/s/silvest/rootfiles/HI/Jpsi";
    $BASE_NAME="DimuonAnalyzer_Jpsi_embedding_500events" .$end_name;

    $import_source = "Jpsinembedded_RECO_cfi";
  }
if($what_particlue=~ m/z/)
  { 
    print "Z0\n";
    $pdgid = 23;
    $mass_max_dimuon = 140;
    $mass_min_dimuon=0;

    $BASE_NAME="DimuonAnalyzer_Z_embedding_500events" .$end_name;    
    $outRootDir    = "/castor/cern.ch/user/s/silvest/rootfiles/HI/Z";
    $import_source = "Zembedded_RECO_cfi";
  }



# cmscode:
$RUN_DIR       = "/afs/cern.ch/user/s/silvest/scratch0/CMSSW_3_5_0_pre3/src/my_code/DimuonAnalyzer/test";
$RUN_CFG       = "$RUN_DIR/batch_dimuonanalyzer_embedding_cfg.py";

$rangeRandom = 10000000;


$OUTPUT_LOG_DIR           = $outRootDir.'/log/';
$OUTPUT_ROOT_DIR          = $outRootDir.'/root/ANA/';
$OUTPUT_CFG_DIR          = $outRootDir.'/cfg/';
    

$FINAL_PARAMETER_SET_NAME = $BASE_NAME;
$CFG                      = $RUN_DIR.'/cfg/'.$FINAL_PARAMETER_SET_NAME.'_cfg.py';
$OUTPUT_LOG               = $FINAL_PARAMETER_SET_NAME.'.log';
$OUTPUT_ROOTFILE          = $FINAL_PARAMETER_SET_NAME.'.root';
$script                   = 'sh/'.$FINAL_PARAMETER_SET_NAME.'.sh';
$skipEvts                 = 0;#($ijob-1)*$NUM_EVENTS_PER_JOB;#change this to 0

$seed = int(rand($rangeRandom));
print "final seed = $seed \n";

    system("cat $RUN_CFG | sed -e 's/CONDOR_MAXEVENTS/$NUM_EVENTS_PER_JOB/g' | sed -e 's/OUTPUTFILE/$OUTPUT_ROOTFILE/g'|  sed -e 's/CONDOR_SKIPEVENTS/$skipEvts/g' | sed -e 's/IMPORTSOURCE/$import_source/g'| \
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
#    system("chmod u+x $script");
#    print `bsub -q 1nh -J $BASE_NAME < $script`;
    if($submit=~ m/true/) {
        print "bsub -q 1nd -J $BASE_NAME < $script\n";
	print `bsub -q 1nd -J $BASE_NAME < $script`;
      }

  

exit 0;





 
