# for easy referencing
$PWD = '/home/is5/assignment_1'
cd $PWD

TAGS=$(ls $PWD/*.fastq.gz | xargs -n 1 basename | sed 's/.fastq.gz//')

echo "Step 0: Fastqc"

for TAG in $TAGS; do
  echo 'Processing' $TAG
  OUTDIR="fastqc/$TAG"; mkdir -p "$OUTDIR" 
  fastqc -o "$OUTDIR" "$PWD/$TAG.fastq.gz" |& tee "$OUTDIR/$TAG.fastqc.log"
done

echo "Step 1: Trimmomatic"

######################################################
# Correction

TAGS=$(ls $PWD/*.fastq.gz | xargs -n 1 basename | sed 's/.fastq.gz//')
List=$TAGS
arr=($List)
for (( i=0; i<${#arr[@]} ; i+=2 )); do
  echo Doing trimmomatic on samples ${arr[i]}.fastq.gz ${arr[i+1]}.fastq.gz
  trimmomatic PE -phred33 ${arr[i]}.fastq.gz ${arr[i+1]}.fastq.gz  \
  ${arr[i]}.paired.fastq.gz ${arr[i]}.unpaired.fastq.gz ${arr[i+1]}.paired.fastq.gz \
  ${arr[i+1]}unpaired.fastq.gz SLIDINGWINDOW:4:17 AVGQUAL:28 |& tee trimmomatic.log; 
done


######################################################

echo "Step 0 once again, because we want to check the quality"

TAGS=$(ls $PWD/*.paired.fastq.gz | xargs -n 1 basename | sed 's/.fastq.gz//')
List=$TAGS
arr=($List)

for TAG in $TAGS; do
  echo 'Processing' $TAG
  OUTDIR="fastqc_new/$TAG"; mkdir -p "$OUTDIR" 
  fastqc -o "$OUTDIR" "$PWD/$TAG.fastq.gz" |& tee "$OUTDIR/$TAG.fastqc.log"
done

echo "Ya, I think it looks better now."
echo "Step 2: Merging reads"

for (( i=0; i<${#arr[@]} ; i+=2 )); do
fastq-join ${arr[i]}.fastq.gz ${arr[i+1]}.fastq.gz -o S1_ |& tee fastq-join.log;
done

echo "Step 3: Importing the dataset"

qiime tools import --type 'SampleData[JoinedSequencesWithQuality]' \
--input-path manifest.csv \
--output-path S1_demux.qza \
--input-format SingleEndFastqManifestPhred33

echo "Step 4: Dereplicate sequences"
  qiime vsearch dereplicate-sequences \
  --i-sequences S1_demux.qza \
  --o-dereplicated-table table.qza \
  --o-dereplicated-sequences rep-seqs.qza
  
echo "Step 5: Cluster sequences"

  qiime vsearch cluster-features-de-novo \
  --i-table table.qza \
  --i-sequences rep-seqs.qza \
  --p-perc-identity 0.99 \
  --o-clustered-table table-dn-99.qza \
  --o-clustered-sequences rep-seqs-dn-99.qza
  
echo "Step 6: Run de novo chimera detection"

  qiime vsearch uchime-denovo \
  --i-sequences rep-seqs-dn-99.qza \
  --i-table table-dn-99.qza \
  --o-chimeras chimeras \
  --o-nonchimeras nonchimeras \
  --o-stats stats

echo "Step 7: Exclude chimeras"

  qiime feature-table filter-features \
  --m-metadata-file nonchimeras.qza \
  --i-table table-dn-99.qza \
  --o-filtered-table table-nonchimeric.qza \

echo "Step 8: Export feature table"
  qiime tools export \
  --input-path  table-nonchimeric.qza \
  --output-path filtered_seqs_export
  
biom convert -i filtered_seqs_export/feature-table.biom -o filtered_seqs_export/OTU_table.txt --to-tsv
