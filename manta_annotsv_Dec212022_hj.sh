#! /bin/bash


##########
# hojinlee
# 20221221
# v5
# this script is written for variant calling and annotation.
# to do these things, Manta and AnnotSV are used respectively.
# follow usage, and then get DEL, DUP and INS ! 
#########

print_usage() {
	echo "$0 -i [bam or cram] -d [path to sample directory]"
	echo "ex ) $0 -i ../../P1.cram -d ./P1"
	echo "  -i : input file ex ) *.bam or *.cram"
	echo "  -d : directory with sampleID as name"
	echo "  -h : show this message"
}

while getopts i:d:h opts; do
	case ${opts} in
		i) input=$OPTARG
			;;
		d) outdir=$OPTARG
			;;
		h) print_usage
			exit;;
		?) print_usage # invalid options
			exit;;
	esac
done

module purge
module load Python/2.7.18-GCCcore-10.2.0


### make sample dir
#mkdir ${outdir}


### manta
python /home/jc2545/programs/manta-1.2.1.centos6_x86_64/bin/configManta.py \
	--bam ${input} \
	--referenceFasta /home/jc2545/ref_data/h_sapiens/1000genomes/2.5/b37/human_g1k_v37_decoy.fasta \
	--runDir ${outdir}

python2.7 ${outdir}/runWorkflow.py -m local


### extract sample name
sample=`basename ${outdir}`
cd ${outdir}
annotsv() {
	### pass
	bcftools view -i 'FILTER="PASS"' ./results/variants/diploidSV.vcf.gz > ./results/variants/diploidSV_edit.vcf
	
	### annotsv     
        $ANNOTSV/bin/AnnotSV/AnnotSV.tcl \
        -SVinputFile ./results/variants/diploidSV_edit.vcf \
        -genomeBuild GRCh37 \
        -rankFiltering 4,5 \
        -outputDir manta_annotsv \
        -typeOfAnnotation split
	
	rm -rf ./results/variants/diploidSV_edit.vcf
}

annotsv

manta_filtering() {	
	vartype=`echo $1`
	mode=`echo $2`

	head -n 1 ./manta_annotsv/diploidSV_edit.annotated.tsv > ./manta_annotsv/header.txt
	
        ##### PR SR
        # if there are not any SR evidence, filtering out the variant
        # 10th field means SR ref
        cat ./manta_annotsv/diploidSV_edit.annotated.tsv | grep -v '^Annot' |\
        awk 'BEGIN {FS="\t"}
        {
        n=split($14, a, "[:]")
        if (n == 6)
                {
                split($14, b, "[:,]")
                if ((b[8] > 0) && (b[10] > 0))
                        print $0
                } 
        }' > ./manta_annotsv/${sample}_manta.pass_temp

        ##### del ins del
        if [[ ${vartype} == "mantadel" ]]
        then
                grep_key="DEL"
        elif [[ ${vartype} == "mantains" ]]
        then
                grep_key="INS"
        elif [[ ${vartype} == "mantadup" ]]
        then
                grep_key="DUP"
        fi
	
	cat ./manta_annotsv/${sample}_manta.pass_temp |\
	awk 'BEGIN {FS="\t"}
	{
	if ($6 == "'${grep_key}'") {print $0}
	}' > ./manta_annotsv/${sample}_manta.pass.${vartype}_temp

	##### dom / rec
        if [[ ${mode} == "dom" ]]
        then
                cat ./manta_annotsv/${sample}_manta.pass.${vartype}_temp |\
		awk 'BEGIN {FS="\t"}
		{
		split($14, a, "[:,]")
		if ((a[1] == "0/1") || (a[1] == "0|1") || (a[1] == "1|0") || (a[1] == "1/0")) 
			{print $0}
		}' > ./manta_annotsv/${sample}_manta.pass.${vartype}.dom_temp 
	elif [[ ${mode} == "rec" ]]
	then
                cat ./manta_annotsv/${sample}_manta.pass.${vartype}_temp |\
                awk 'BEGIN {FS="\t"}
                {
                split($14, a, "[:,]")
                if ((a[1] == "1/1") || (a[1] == "1|1")) 
                        {print $0}
                }' > ./manta_annotsv/${sample}_manta.pass.${vartype}.rec_temp
	fi

	##### AF filtering
        if [[ ${mode} == "dom" ]]
        then
        # filtering 1000g_AF <= 0.001 (col 45), IMH_AF <= 0.001 (col 48), HI_DDDpercent <= 10 (col 70), pLI_ExAC >= 0.9 (col 73)
        cat ./manta_annotsv/${sample}_manta.pass.${vartype}.dom_temp |\
        awk -F '\t' '
        {
        if (($45 <= 0.001) && ($45 != "") && ($48 <= 0.001) && ($48 != "")) print $0
        } 
        ' > ./manta_annotsv/${sample}_manta.pass.${vartype}.dom.filtered_temp
	
        elif [[ ${mode} == "rec" ]]
        then
        # filtering 1000g_AF <= 0.01 (col 45), IMH_AF <= 0.01 (col 48), HI_DDDpercent <= 10 (col 70), pLI_ExAC >= 0.9 (col 73)
        cat ./manta_annotsv/${sample}_manta.pass.${vartype}.rec_temp |\
        awk -F '\t' '
        {
        if (($45 <= 0.01) && ($45 != "") && ($48 <= 0.01) && ($48 != "")) print $0
        }
        ' > ./manta_annotsv/${sample}_manta.pass.${vartype}.rec.filtered_temp
	fi

	cat ./manta_annotsv/header.txt ./manta_annotsv/${sample}_manta.pass.${vartype}.${mode}.filtered_temp > ./manta_annotsv/${sample}_manta.pass.${vartype}.${mode}.filtered.txt
	##### remove temp file
	rm -rf ./manta_annotsv/*temp
	rm -rf ./manta_annotsv/header.txt
}

# deletion
manta_filtering mantadel dom
manta_filtering mantadel rec
# duplication
manta_filtering mantadup dom
manta_filtering mantadup rec
# insertion
manta_filtering mantains dom
manta_filtering mantains rec

