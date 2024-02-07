#!/bin/bash
CODE=${1/.py/}; shift

if [[ "$CODE" == "" || "$CODE" == "--help" || "$CODE" == "-h" || "$CODE" == "-?" ]]; then
    echo "scripts/prun.sh config.py [-j <N>] --<version>  <Sample> <OutputName> [ --noclean ] [ --inline-customize \"<options>\" ]"
    exit 0;
fi

N=8
if [[ "$1" == "-j" ]]; then N=$2; shift; shift; fi;

OPTS=""

if [[ "$1" == "--125X_v0" ]]; then
    shift;
    MAIN=/eos/cms/store/cmst3/group/l1tr/gpetrucc/12_5_X/NewInputs125X/150223/$1
    PREFIX="inputs125X_"
elif [[ "$1" == "--110X_v2" ]]; then
    shift;
    MAIN=/eos/cms/store/cmst3/group/l1tr/gpetrucc/11_1_0/NewInputs110X/110121.done/$1
    PREFIX="inputs110X_"
    OPTS="${OPTS} --replace Phase2C17I13M9=Phase2C9 "
    OPTS="${OPTS} --replace GeometryExtended2026D88=GeometryExtended2026D49 "
    OPTS="${OPTS} --replace 125X_mcRun4_realistic_v2=123X_mcRun4_realistic_v3 "
elif [[ "$1" == "--110X_v3" ]]; then
    shift;
    MAIN=/eos/cms/store/cmst3/group/l1tr/gpetrucc/12_3_X/NewInputs110X/220322/$1
    PREFIX="inputs110X_"
    OPTS="${OPTS} --replace Phase2C17I13M9=Phase2C9 "
    OPTS="${OPTS} --replace GeometryExtended2026D88=GeometryExtended2026D49 "
    OPTS="${OPTS} --replace 125X_mcRun4_realistic_v2=123X_mcRun4_realistic_v3 "
else 
    echo "You mush specify the version of the input samples to run on "
    echo "   --125X_v0 : 125X Phase2Fall22 MC in 12_5_X"
    echo "   --110X_v3 : 110X HLT MC inputs remade in 12_3_0_pre4; add --inline-customize 'oldInputs_12_3_X()'"
    echo "   --110X_v2 : 110X HLT MC inputs remade in 11_1_6;      add --inline-customize 'oldInputs_11_1_6()'"
    exit 1;
fi;
 
if [[ "$L1TPF_LOCAL_INPUT_DIR" != "" ]] && test -d $L1TPF_LOCAL_INPUT_DIR; then
    L1TPF_LOCAL_MAIN=$L1TPF_LOCAL_INPUT_DIR/$(basename $(dirname $MAIN))/$(basename $MAIN)
    if test -d $L1TPF_LOCAL_MAIN; then
        echo "Will use local directory $L1TPF_LOCAL_MAIN";
        MAIN=$L1TPF_LOCAL_MAIN;
    fi;
fi;

INPUT=$1; shift
if [[ "$1" != "" ]]; then
    OUTPUT="$1";
    shift;
else
    OUTPUT="${INPUT}";
fi;

clean=true
if [[ "$1" == "--noclean" ]]; then
    clean=false;
    shift
fi

PSCRIPT=$CMSSW_BASE/src/FastPUPPI/NtupleProducer/python/scripts/
$PSCRIPT/cmsSplit.pl --files "$MAIN/${PREFIX}*root" --label ${OUTPUT} ${CODE}.py --bash --n $N --rrb $OPTS  $* && bash ${CODE}_${OUTPUT}_local.sh 

if $clean; then
    REPORT=$(grep -o "tee \S\+_${OUTPUT}.report.txt" ${CODE}_${OUTPUT}_local.sh  | awk '{print $2}');
    if [[ "$REPORT" != "" ]] && test -f ${REPORT}; then
        if [[ "$N" == "1" ]] && grep -q "Job completed without exceptions" ${REPORT} || grep -q ', 0 failed' $REPORT; then
           bash ${CODE}_${OUTPUT}_cleanup.sh
        else
           echo "Failed jobs... will not delete anything."
           exit 1;
        fi;
    else
        bash ${CODE}_${OUTPUT}_cleanup.sh
    fi;
fi
