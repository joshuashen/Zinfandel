#!/bin/bash
#
# Wrapper for corona lite
# Borrowed from David Rio Deiros at BCM-HGSC
###  June 2008
# 

skip_matching=0; skip_pairing=1; skip_snips=1
to_lsf=0   # submit the jobs to the lsf cluster
to_local=1 # submit the jobs locally

# Things that have to be modified
main_dir="/hgsc/solid" # Where the cmaps_human are located
## PATH="$PATH:/hgsc/solid/hgsc_solid/bin"
project_dir="/hgsc/solid/analysis/Solid0100_20080514_1_Pilot2_YRI_B_library_MP" # where you'll run your project (most likely your current dir)
csfasta_files="/hgsc/solid/csfastas/Solid0100_20080514_1_Pilot2_YRI_B_library_MP" # where to search for the csfasta files

cmap_file="$main_dir/cmaps_human/cmap.txt"
full_ref="$main_dir/cmaps_human/full_hsap36.fasta"
cpus=4 # Number of cpus to use when calling run_local_jobs (should be #cpus/2)

logs="$project_dir/logs"
input_data="$project_dir/.input.txt"
corona_version="0.23"
bin="$main_dir/corona_lite_v$corona_version/bin"
local_tool="run_JOBS_locally.rb"
lsf_tool="submit_scripts_to_LSF.pl"
queue="allram"
wait_queue="login"
normal_queue="normal"

error()
{
  echo
  echo "ERROR: $1"
  exit
}

# Let's test everything is there
echo "making sure everything is there ..."
[ ! -f $cmap_file ] && error "$cmap_file not found"
[ ! -f $full_ref ] && error "$full_ref not found"
[ `which $local_tool` ] || error "$local_tool not found"
[ `which ruby` ] || error "ruby not found"
if [ $to_lsf == 1 ] 
then
  [ `which bsub` ] || error "You don't seem to have the lsf tools"
fi


# Let's get the color space reads from the run
echo "Finding csfasta files..."
if [ ! -f $input_data ]
then
  for c in `find $csfasta_files \( -name "*R3.csfasta" -o -name "*F3.csfasta" \) | xargs ls`
  do
    `echo $c | grep "R3" >/dev/null` && echo "r3_csfasta=\"$c\"" >> $input_data
    `echo $c | grep "F3" >/dev/null` && echo "f3_csfasta=\"$c\"" >> $input_data
    source $input_data
  done
else
  source $input_data
fi
if [[ ! -f $r3_csfasta || ! -f $f3_csfasta ]]
then
  error "ERROR: csfasta files not found."
fi
echo "r3: $r3_csfasta"
echo "f3: $f3_csfasta"

mkdir -p $logs

# Level 1: mapping
##################
matching_output="$project_dir/output/matching"
t_length="25"
n_errors="2"
max_hits="10"

if [ $skip_matching -ne 1 ]
then
  mkdir -p $matching_output/F3
  mkdir -p $matching_output/R3

  # a. generate lsf jobs (R3 reads)
  #
  $bin/matching_large_genomes_cmap_save_script.pl \
  -csfasta $r3_csfasta \
  -dir $matching_output/R3 \
  -cmap $cmap_file \
  -t $t_length \
  -e $n_errors \
  -z $max_hits

  # b. generate lsf jobs (F3 reads)
  #
  $bin/matching_large_genomes_cmap_save_script.pl \
  -csfasta $f3_csfasta \
  -dir $matching_output/F3 \
  -cmap $cmap_file \
  -t $t_length \
  -e $n_errors \
  -z $max_hits

  if [ $to_lsf -eq 1 ] 
  then
    $lsf_tool -q $queue -i $wait_queue -j $matching_output/R3/JOB_LIST.txt &> $logs/lsf_r3_jobs.log &
    $lsf_tool -q $queue -i $wait_queue -j $matching_output/F3/JOB_LIST.txt &> $logs/lsf_f3_jobs.log &
  elif [ $to_local -eq 1 ]
  then
    $local_tool "mapping_F3" $matching_output/F3/JOB_LIST.txt $project_dir $cpus &
    $local_tool "mapping_R3" $matching_output/R3/JOB_LIST.txt $project_dir $cpus &
  else
    error "Don't know where to run jobs."
  fi

  # Wait for the other jobs to finish
  wait
fi

# Level 2: pairing
##################
n_panels_per_group="10"
pairing_output="$project_dir/output/pairing"
isize_output="$project_dir/output/isize"

if [ $skip_pairing -ne 1 ]
then
  echo "Finding insert size"
  mkdir -p $isize_output
  max_isize=`find_insert_size.rb | ruby -ne 'print /\s([0-9]+)$/.match($_)[1] if /interval/.match($_)'`
  echo "isize = 100 $max_isize"
  mv *.png $isize_output/

  mkdir -p $pairing_output/split
  mkdir -p $pairing_output/tmp

  cd $pairing_output/tmp
  echo "Sorting R3 job..."
  echo "$bin/radsort $matching_output/R3/`basename $r3_csfasta`.ma.25.2 6 > $pairing_output/`basename $r3_csfasta`.ma.25.2.sorted" > sort1.sh
  echo -e "1\t./sort1.sh" > sort.lsf.txt
  
  echo "Sorting F3 job..."
  echo "$bin/radsort $matching_output/F3/`basename $f3_csfasta`.ma.25.2 6 > $pairing_output/`basename $f3_csfasta`.ma.25.2.sorted" > sort2.sh
  echo -e "2\t./sort2.sh" >> sort.lsf.txt
  echo -e "3\techo \"waiting for R3/F3 to finish\"\t1,2" >> sort.lsf.txt
  chmod 755 ./sort*.sh

  [ $to_lsf -eq 1 ] && $lsf_tool -q $normal_queue -i $wait_queue -j $pairing_output/tmp/sort.lsf.txt
  [ $to_local -eq 1 ] && $local_tool "pairing_sort" $pairing_output/tmp/sort.lsf.txt $project_dir $cpus

  cd ..
  # split
  echo "Splitting.."
  echo "$bin/split_sorted_match_files_into_groups_of_panels.pl -r3 $pairing_output/`basename $r3_csfasta`.ma.25.2.sorted -f3 $pairing_output/`basename $f3_csfasta`.ma.25.2.sorted -n $n_panels_per_group -dir $pairing_output/split" > split.sh
  echo -e "1\t./split.sh" > split.lsf.txt
  echo -e "2\techo \"waiting for split to finish\"\t1" >> split.lsf.txt
  chmod 755 ./split.sh

  [ $to_lsf -eq 1 ] && $lsf_tool -i $wait_queue -j $pairing_output/split.lsf.txt
  [ $to_local -eq 1 ] && $local_tool "pairing_split" $pairing_output/split.lsf.txt $project_dir $cpus

  # actual pairing
  echo "Pairing.."
  cd $pairing_output/split
  $bin/pairing_by_panel_save_script.pl \
  -dir_f3 F3_match_files -dir_r3 R3_match_files \
  -insert_start 100 -insert_end $max_isize -e 4 -dir runs \
  -ref $full_ref \
  -format multiple -mismatch_threshold 2

  [ $to_lsf -eq 1 ] && $lsf_tool -i $wait_queue -j $pairing_output/split/JOB_LIST.txt
  [ $to_local -eq 1 ] && $local_tool "pairing" $pairing_output/split/JOB_LIST.txt $project_dir $[$cpus * 2]
fi

# Level 3: SNPs
###############
snps_output="$project_dir/output/snps"

if [ $skip_snips -ne 1 ]
then
  mkdir -p $snps_output

  $bin/multi_chr_pairing_parser.pl \
  -mates $pairing_output/split/runs/F3_R3.mates.unique \
  -o_dir $snps_output/chr_unique

  cd $snps_output
  $bin/consensus_prep_and_wrapper_cmap_save_script.pl -mates_dir chr_unique \
  -ml_f3 25 -ml_r3 25 -e_f3 2 -e_r3 2 -insert_start 400 -insert_end 1450 \
  -cmap $cmap_file -o_dir runs #-12only

  [ $to_lsf -eq 1 ] && $lsf_tool -q $queue -i $wait_queue -j $snps_output/JOB_LIST.txt
  [ $to_local -eq 1 ] && $local_tool "snip_calling" $snps_output/JOB_LIST.txt $project_dir $[$cpus * 2]
fi
