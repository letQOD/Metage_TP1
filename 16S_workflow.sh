#! /bin/bash


###chmod +x 16S_workflow.sh to make it executable
###commandline : ./16S_workflow.sh dirfq outdir in tp_1 folder
###!!!change soft/databases paths if executed elsewhere!!!


dirfq=$1
outdir=$2


mkdir $outdir



#####################################
########## Quality check ############
#####################################



### fastqc ###
mkdir $outdir/fastqc
cd $dirfq
gunzip *.gz
ls *R1.fastq > tax_fqc.txt
N=($(cat tax_fqc.txt | wc -l))
for f in $(seq 1 $N);
do sample=($(sed -n "$f p" tax_fqc.txt));
ID=${sample%_*};
mkdir ../$2/fastqc/$ID;
fastqc -o $ID -f fastq $ID""_R1.fastq $ID""_R2.fastq; done

#####################################
############# Trimming ##############
#####################################




#### Alientrimmer ###
alien=databases/contaminants.fasta
cd ..
mkdir $2/alientrimmer
N=($(cat $1/tax_fqc.txt | wc -l))
for f in $(seq 1 $N);
do sample=($(sed -n "$f p" $1/tax_fqc.txt));
ID=${sample%_*};
java -jar soft/AlienTrimmer.jar -if fastq/$ID""_R1.fastq -ir fastq/$ID""_R2.fastq -c $alien -q 20 -of $2/alientrimmer/$ID""_R1.flt.fastq -or $2/alientrimmer/$ID""_R2.flt.fastq -os $2/alientrimmer/$ID""_RS.flt.fastq
done


#####################################
############# Vsearch ###############
#####################################
mkdir $2/vsearch
vsearch=soft/vsearch

N=($(cat $1/tax_fqc.txt | wc -l))
for f in $(seq 1 $N);
do sample=($(sed -n "$f p" $1/tax_fqc.txt));
ID=${sample%_*};
suffix=";sample=$ID"
$vsearch --fastq_mergepairs $2/alientrimmer/$ID""_R1.flt.fastq --reverse $2/alientrimmer/$ID""_R2.flt.fastq --label_suffix $suffix --fastaout $2/vsearch/$ID""_vsearch.fasta;
done

cat $2/vsearch/*.fasta > amplicon.fasta

sed -i -e "s/ //g" amplicon.fasta

### deduplication full length and prefix modes ###

$vsearch --derep_fulllength $2/vsearch/amplicon.fasta --log=vsearch_fl_log --sizeout --minuniquesize 10 --output $2/vsearch/amplicon_dedfl.fasta --uc $2/vsearch/amplicon_dedfl.uc


$vsearch --derep_prefix $2/vsearch/amplicon.fasta --log=vsearch_pref.log --sizeout --minuniquesize 10 --output $2/vsearch/amplicon_dedpref.fasta --uc $2/vsearch/amplicon_dedpref.uc


### remove chimeras ###

$vsearch --uchime_denovo $2/vsearch/amplicon_dedfl.fasta --nonchimeras $2/vsearch/amplicon_dedfluchim.fasta --sizein --chimeras $2/vsearch/chimeras_fl.fasta  ##on full length data

$vsearch --uchime_denovo $2/vsearch/amplicon_dedpref.fasta --nonchimeras $2/vsearch/amplicon_dedprefuchim.fasta --sizein --chimeras $2/vsearch/chimeras_pref.fasta  ##on prefix data


### clustering ###

$vsearch --cluster_size $2/vsearch/amplicon_dedfluchim.fasta --id 0.97 --sizein --sizeout --relabel OTU_ --centroids $2/vsearch/amplicon_dedfluchim_OTU.fasta   ##on full length data

$vsearch --cluster_size $2/vsearch/amplicon_dedprefuchim.fasta --id 0.97 --sizein --sizeout --relabel OTU_ --centroids $2/vsearch/amplicon_dedprefuchim_OTU.fasta   ##on prefix data


### table d'abondance ###

$vsearch-usearch_global $2/vsearch/amplicon.fasta --otutabout $2/vsearch/dedfluchim_OTU  -db $2/vsearch/amplicon_dedfluchim_OTU.fasta --id 0.97     ##on full length data


$vsearch -usearch_global $2/vsearch/amplicon.fasta --otutabout $2/vsearch/dedprefuchim_OTU -db $2/vsearch/amplicon_dedprefuchim_OTU.fasta --id 0.97     ##on prefix data



### annotation ###
$vsearch --usearch_global $2/vsearch/amplicon_dedfluchim_OTU.fasta --userout  $2/vsearch/amplicon_dedfluchim_annot.fasta --db databases/mock_16S_18S.fasta  --id 0.9 --top_hits_only --userfields    query+target    ##on full length data

$vsearch --usearch_global $2/vsearch/amplicon_dedprefuchim.fasta --userout  $2/vsearch/amplicon_dedprefuchim_annot.fasta --db databases/mock_16S_18S.fasta  --id 0.9 --top_hits_only --userfields query+target      ##on prefix data   

sed '1iOTU\tAnnotation' -i $2/vsearch/amplicon_dedfluchim_annot.fasta        ##on full length data
sed '1iOTU\tAnnotation' -i $2/vsearch/amplicon_dedprefuchim_annot.fasta      ##on prefix data




