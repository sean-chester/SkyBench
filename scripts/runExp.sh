#!/bin/bash
# -----------------------------------------------------------------------------
# command runner:
function docmd() {
  if [ $# -ne 2 ] ; then echo "do: $1" ; fi
  if( eval "(" $1 ")" )
  then 
    true
  else
	  echo "Command FAILED" ;
	  exit 1 ;
  fi
}

function printUsage() {
  cat << EOF
Highly configurable SkyBench experiment runner. Multiple values
for each option can be passed using double quotes, e.g.:
  -t "1 2 4 8"

Usage: ${0##*/} [-h -v -p] [-i dir] [-o dir] [-t num_threads] 
                [-s skylines] [-m T|D|P] [-d dims] [-c card] 
                [-x dist] [-a block_size] [-q pq_size]
       
  -h  display this help and exit
  -v  print all configuration variables and their values
  -i  input directory (by default: ./workloads). Note that input
      files follow the following naming convention: data-${dist}-${dim}-${card}.csv
  -o  output directory (by default: ./results)
  -t  number of threads (by default: 8)
  -s  skyline algorithms to test, e.g., -s "bskytree hybrid"
  -m  measure T=time (default), D=# ofdominance tests. Accepts only single option.
      Note that with D the running time increases significantly due to measuring
      overhead (esp. with many threads). The counts are accurate, though.
  -p  turn on profiler (i.e., compile with -DPROFILER=1 flag) that
      breaks down runtimes into 'pq-filter', 'select pivot', 
      'partition', 'phaseI', 'phaseII', and 'compress' phases.
  -x  distribution (by default: "C E A")

Only one of the following options can have multiple values. This
is the parameter that will be varied on the x-axis (the 1st column
in the output .csv file).
  -d  dimensionality (by default: 12)
  -c  cardinality (by default: 1000000)
  -a  alpha block size (by default: 1024)
  -q  priority queue size (used only by Hybrid, by default 16)

EXAMPLES:
To run our standard cardinality test:
./runExp.sh -p -i workloads/ -t 8 -c "500000 1000000 2000000 4000000 8000000" -d 12 -s "bskytree hybrid"
            
To run our standard dimensionality test:
./runExp.sh -p -i workloads/ -t 8 -c 1000000 -d "2 4 6 8 10 12 14 16 18 20 22 24" -s "bskytree hybrid"

EOF
}

function doTest() {
  OUTPUT=$(docmd "$program -f $datafile -t \"${threads}\" -s \"${algs}\" -a ${q_accum} -q ${pq}" 0) ;
  out_line="${out_line}${OUTPUT}" ;
  echo -e "$out_line";
  echo -e "$out_line" >> $res_file ;
}

function doPQ {
  for pq in $pq_sizes ; do
    if [[ $var_arg == "q" ]] ; then
      out_line="${pq}";
    fi
    doTest
  done
}

function doAlpha() {
  for q_accum in $alpha ; do
    if [[ $var_arg == "a" ]] ; then
      out_line="${q_accum}";
    fi
    doPQ ;
  done
}

function doDim() {
  for dim in $dims ; do
    if [[ $var_arg == "d" ]] ; then
      # Compile for each dimensionality
      make clean > /dev/null ;
      make -j4 all DT=${dt} DIMS=${dim} V=NVERBOSE \
                   PROFILER=${profiler} > /dev/null ;
      out_line="${dim}";
    fi
    # Input data filename can be identified already:
    datafile="${input_dir}/${file_prefix}-${dist}-${dim}-${card}.csv" ;
    # echo "DATAFILE: ${datafile}" ;
    doAlpha ;
  done
}

function doCard() {
  for card in $cards ; do
    if [[ $var_arg == "c" ]] ; then
      out_line="${card}";
    fi
    doDim ;
  done
}

# Print configuration
function printConfig() {
  echo "Configuration:"
  echo " input dir = \"$input_dir\"" ;
  echo " output dir = \"$output_dir\"" ;
  echo " threads = ${threads}" ;
  echo " dims = ${dims}" ;
  echo " cards = ${cards}" ;
  echo " dists = ${dists}" ;
  echo " skylines = ${algs}" ;
  echo " measure = ${measure}" ;
  echo " alpha sizes = ${alpha}" ;
  echo " pq sizes = ${pq_sizes}" ;
  echo " breakdown = ${profiler}" ;
  echo " multi-value option = -${var_arg}"
  echo "" ;
}
# -----------------------------------------------------------------------------

# Default values:
input_dir="./workloads" ;
output_dir="./results" ;
threads="8" ;
# hybrid is used by default:
algs="hybrid" ;
dims="12" ;
cards="1000000" ;
dists="C E A" ; # "core inde anti"
alpha="1024" ;
pq_sizes="16"
measure="T" ;
print_config=0 ;

# Constant values:
# breakdown time runs:
#sub_times="init phase1 phase2 sort/shift total" #hard-coded in PROFILER
sub_times="pq-filter pivot partition phase1 phase2 compress total"
st_algs="bskytree" # single-threaded algorithms
program="./bin/SkyBench" ;
file_prefix="data" ;
dt="0" ;
profiler="0" ;
header="" ;
var_arg="x" ; # variable argument

if [[ "$#" -eq 0 ]] ; then
  printUsage ;
  exit 0 ;
fi

while getopts ":hpvt:i:o:m:d:c:x:s:a:q:" opt; do
  case $opt in
    h)
      printUsage ;
      exit 0 ;
      ;;
    p)
      profiler=1 >&2 ;
      ;;
    v)
      print_config=1 >&2 ;
      ;;
    t)
      threads="$OPTARG" >&2 ;
      ;;
    m)
      measure="$OPTARG" >&2 ;
      num=$(echo "$measure" | wc -w)
      if [ $num -ne 1 ] ; then
        echo "ERROR: only one value is supported for -m" ;
        exit 1 ;
      fi
      ;;
    i)
      input_dir="$OPTARG" >&2 ;
      if [ ! -d "$input_dir" ]; then
        echo "ERROR: no such input directory \"$input_dir\" !" ;
        exit 1 ;
      fi
      ;;
    o)
      output_dir="$OPTARG" >&2 ;
      if [ ! -d "$output_dir" ]; then
        mkdir "$output_dir" ;
      fi
      ;;
    d)
      dims="$OPTARG" >&2 ;
      num=$(echo "$dims" | wc -w) ;
      if [ $num -ne 1 ] ; then
        if [ ! ${var_arg} == "x" ] ; then
          echo "ERROR: only one option can have multiple values!"
          exit 1 ;
        fi
        header="Dimensionality" ;
        var_arg="d" ;
      fi
      ;;
    c)
      cards="$OPTARG" >&2 ;
      num=$(echo "$cards" | wc -w) ;
      if [ $num -ne 1 ] ; then
        if [ ! ${var_arg} == "x" ] ; then
          echo "ERROR: only one option can have multiple values!"
          exit 1 ;
        fi
        header="Cardinality" ;
        var_arg="c" ;
      fi
      ;;
    x)
      dists="$OPTARG" >&2 ;
      ;;
    s)
      algs="$OPTARG" >&2 ;
      ;;
    a)
      alpha="$OPTARG" >&2 ;
      num=$(echo "$alpha" | wc -w) ;
      if [ $num -ne 1 ] ; then
        if [ ! ${var_arg} == "x" ] ; then
          echo "ERROR: only one option can have multiple values!"
          exit 1 ;
        fi
        header="Alpha" ;
        var_arg="a" ;
      fi
      ;;
    q)
      pq_sizes="$OPTARG" >&2 ;
      num=$(echo "$pq_sizes" | wc -w) ;
      if [ $num -ne 1 ] ; then
        if [ ! ${var_arg} == "x" ] ; then
          echo "ERROR: only one option can have multiple values!"
          exit 1 ;
        fi
        header="PQ" ;
        var_arg="q" ;
      fi
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2 ;
      printUsage ;
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2 ;
      exit 1
      ;;
  esac
done

STARTTIME=$(date +%s) ;

if [[ $print_config == 1 ]] ; then
  printConfig ;
fi

# Prepare header and print it to screen and to file
for a in $algs
do
  if [[ $st_algs == *${a}* ]] ; then
    if [[ $profiler == 1 ]] ; then
      for sub in $sub_times ; do
        header="${header} ${a}_${sub}" ;
      done
    else
      header="${header} ${a}" ;
    fi
  else
    for t in $threads ; do
      if [[ $profiler == 1 ]] ; then
	      for sub in $sub_times ; do
          header="${header} ${a}_t${t}_${sub}" ;
        done
      else
        header="${header} ${a}_t${t}" ;
      fi
		done
  fi
done #header is ready

if [[ ! ${var_arg} == "d" ]] ; then
  echo "Compiling with 'DIM=${dims}'.."
  make clean > /dev/null ;
  make -j8 all DT=${dt} DIMS=${dims} V=NVERBOSE \
               PROFILER=${profiler} > /dev/null ;
fi

# Special case when only dist (-x) is varied (prints all results 
# in the same file)
if [[ ${var_arg} == "x" ]] ; then
  header="Distribution${header}"
  res_file="data-${measure}-${cards}-${dims}-dist.csv";
  res_file="${output_dir}/${res_file}"
  echo "Output file: $res_file"
  echo -e "$header" > "$res_file" ;
  echo -e "$header" ;
  
  for dist in $dists ; do
    out_line="${dist}";
    doCard ;
  done
  
else
  # Print results for each distribution in separate files:
  for dist in $dists ; do
    # Generate result filename (following our naming conventions)
    case $var_arg in
      d)
        res_file="data-${measure}-${dist}-${cards}-dim.csv";
        ;;
      c)
        res_file="data-${measure}-${dist}-${dims}-card.csv";
        ;;
      a)
        res_file="data-${measure}-${dist}-${cards}-${dims}-alpha.csv";
        ;;
      q)
        res_file="data-${measure}-${dist}-${cards}-${dims}-pq.csv";
        ;;
    esac
    #res_file is ready:
    res_file="${output_dir}/${res_file}"
    echo "Output file: $res_file"
    echo -e "$header" > "$res_file" ;
    echo -e "$header" ;

    doCard ;
done
fi

ENDTIME=$(date +%s) ;
echo "" ;
echo "DONE [at $(date '+%H:%M:%S'), duration: $(($ENDTIME - $STARTTIME)) s.]" ;
