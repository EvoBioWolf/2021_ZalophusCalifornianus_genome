#These commands were repeated for each chromosome in VGPassembly aligned against SRassembly.v2 and SRassemby.v3 e.g., for i in *.fa; do sbatch nucmer_submission_script.sh $i GCF_900631625.1_zalCal2.2_genomic_masked.fna; done


name=`ls -1 ${1}|cut -f 1,2 -d '_'`


nucmer -c 100 -p $name"_nucmer_mumref_c100" ${1} ${2}

nucmer -p $name"nucmer_mumref" ${1} ${2}

nucmer -p $name"nucmer_mumref_g120" -g 120 ${1} ${2}

nucmer -c 100 -b 400 -p $name"nucmer_mumref_c100_b400" ${1} ${2}

nucmer --mum -p $name"nucmer_mum" ${1} ${2}

nucmer --mum -c 100 -p $name"nucmer_mum_c100" ${1} ${2}

nucmer --mum -c 100  -b 400 -p $name"nucmer_mum_c100" ${1} ${2}

nucmer --mum -p $name"nucmer_mum_g120" -g 120 ${1} ${2}


