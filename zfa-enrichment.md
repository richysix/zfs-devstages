## BioLayout 0.94 correlation network

Download required files
```bash
cd output/zfa/

wget http://www.berkeleybop.org/ontologies/zfa.obo
wget http://ontologizer.de/cmdline/Ontologizer.jar
wget http://zfin.org/downloads/phenoGeneCleanData_fish.txt
wget http://zfin.org/downloads/wildtype-expression_fish.txt

# get ZFIN ID to ZFA terms
cat <(cut -f 3,8 phenoGeneCleanData_fish.txt | grep ZFA) <(cut -f 1,4 wildtype-expression_fish.txt | grep ZFA) | sort -u \
| awk '{ print "ZFIN\t" $1 "\t" $1 "\t\t" $2 "\tRef\tND\t\tC\t" $1 "\t\tgene_product\ttaxon:7955\t20150515\tZFIN" }' > output/zfa/annotations.txt

cd ../../
```

The file with ZFIN IDs to Ensembl IDs is `danio_rerio_e85_zfin.txt`
It can be downloaded from the [topgo-wrapper](https://github.com/iansealy/topgo-wrapper/blob/master/data/danio_rerio_e85_go.txt) github repository.

The layout file is `dataFiles/zfs-grcz10.tpm_r-0.94_pearson.layout`

convert to cluster_info file
```bash
grep //NODECLASS dataFiles/zfs-grcz10.tpm_r-0.94_pearson.layout | sed -e 's|"||g' | \
perl -F"\t" -lane 'BEGIN{ %info_for = (); }
{next if $F[0] ne "//NODECLASS";
$info_for{ $F[1] }{ $F[3] } = $F[2];
}
END{ print join("\t", qw{ name gene_id biotype cluster chr start end } );
my @keys = ( "Gene ID", qw{ Biotype MCL_2.2_6 Chr Start End });
foreach $gene ( sort keys %info_for ){
# filter out genes not assigned to clusters
$info_for{$gene}{"MCL_2.2_6"} = $info_for{$gene}{"MCL_2.2_6"} eq "" ?
  "NA" : $info_for{$gene}{"MCL_2.2_6"};
print join("\t", $gene, map { $info_for{$gene}{$_} } @keys );
}
}' > dataFiles/zfs-grcz10.tpm_r-0.94_pearson.cluster.info.txt
```

```bash
# total genes = 12296
grep -v cluster dataFiles/zfs-grcz10.tpm_r-0.94_pearson.cluster.info.txt | wc -l
   13162

# genes NA to cluster = 1933
cut -f4 dataFiles/zfs-grcz10.tpm_r-0.94_pearson.cluster.info.txt | grep -v cluster | sort | uniq -c | grep NA
1723 NA

```

```bash
# get cluster sizes
mkdir -p output/zfa/0.94-e85
cut -f4 dataFiles/zfs-grcz10.tpm_r-0.94_pearson.cluster.info.txt | grep -v cluster | \
sort | uniq -c | sort -k2,2 > output/zfa/0.94-e85/cluster-sizes.txt
```

Make file of all genes.
```bash
cut -f2 dataFiles/zfs-grcz10.tpm_r-0.94_pearson.cluster.info.txt | sort | \
join -t$'\t' - output/zfa/danio_rerio_e85_zfin.txt | \
cut -f2 | sort -u > output/zfa/all-zfin-genes.txt

# check numbers
wc -l output/zfa/all-zfin-genes.txt
12360 output/zfa/all-zfin-genes.txt

# Genes that map to more than one ZFIN ID
cut -f2 dataFiles/zfs-grcz10.tpm_r-0.94_pearson.cluster.info.txt | sort | \
join -t$'\t' - output/zfa/danio_rerio_e85_zfin.txt | \
cut -f1 | sort | uniq -c | \
perl -lane 'if( $F[0] > 1 ){ print $_ }' | wc -l
78

# Genes that map to the SAME ZFIN ID
cut -f2 dataFiles/zfs-grcz10.tpm_r-0.94_pearson.cluster.info.txt | sort | \
join -t$'\t' - output/zfa/danio_rerio_e85_zfin.txt | cut -f2 | sort | \
uniq -c | perl -lane 'if( $F[0] > 1 ){ print $_ }' | wc -l
56

# Genes that don't map to a ZFIN ID
cut -f2 dataFiles/zfs-grcz10.tpm_r-0.94_pearson.cluster.info.txt | sort | \
join -t$'\t' -v1 - output/zfa/danio_rerio_e85_zfin.txt | \
wc -l
830
```

## Ontologizer
```bash
# run for all clusters sequentially
numClusters=254
ROOT=$(pwd)
ZFA_DIR=$ROOT/output/zfa
# ontologizerJar=$ZFA_DIR/Ontologizer.jar
ontologizerJar=~rw4/bin/Ontologizer.jar
for cluster in $( seq -w 1 $numClusters | sed -e 's|^|Cluster|' )
do
# make a separate directory for each cluster
mkdir -p $ZFA_DIR/0.94-e85/$cluster
cd $ZFA_DIR/0.94-e85/$cluster
# make a log file
touch $cluster.log
# get the Ensembl IDs of the genes in the cluster
grep $cluster $ROOT/data/zfs-grcz10.tpm_r-0.94_pearson.cluster.info.txt | cut -f2 | \
sort -u > ens-genes-in-cluster.txt
# convert to ZFIN IDs
join -t$'\t' ens-genes-in-cluster.txt $ZFA_DIR/danio_rerio_e85_zfin.txt | \
cut -f2 | sort -u > sig-zfin-genes.txt
# log numbers of IDs
echo $(wc -l ens-genes-in-cluster.txt) >> $cluster.log
echo $(wc -l sig-zfin-genes.txt) >> $cluster.log
# run Ontologizer
java -jar $ontologizerJar --mtc Benjamini-Hochberg --calculation Parent-Child-Union \
--association $ZFA_DIR/annotations.txt --go $ZFA_DIR/zfa.obo \
--population $ZFA_DIR/all-zfin-genes.txt \
--studyset $ZFA_DIR/0.94-e85/$cluster/sig-zfin-genes.txt 2>> $cluster.log
# make all and sig files
# Need to multiply the p-values by the number of clusters to account for testing multiple clusters
# All file
grep -v ^ZFS: table-sig-zfin-genes-Parent-Child-Union-Benjamini-Hochberg.txt | \
perl -F"\t" -lane '$num_clusters = "'$numClusters'";
if( $. == 1 ){ print join("\t", @F[0,9,10,12] ); }
else{
  print join("\t", $F[0], $F[9]*$num_clusters, $F[10]*$num_clusters, $F[12] );
}' > all.tsv
# sig file
grep -v ^ZFS: table-sig-zfin-genes-Parent-Child-Union-Benjamini-Hochberg.txt | \
perl -F"\t" -lane '$num_clusters = "'$numClusters'";
if( $. == 1 ){ print join("\t", @F[0,9,10,12] ); }
else{ # filter at 0.05 / number of clusters
  if ($F[10] * $num_clusters < 0.05 ){
    print join("\t", $F[0], $F[9]*$num_clusters, $F[10]*$num_clusters, $F[12] );
  }
}' > sig.tsv
# log numbers of genes
echo $(wc -l all.tsv) >> $cluster.log
echo $(wc -l sig.tsv) >> $cluster.log
cd $ROOT
echo "$cluster - DONE"
done
```

### Clusters with significant results
```bash
# List of clusters with significant results (excluding ZFA:0001093 - unspecified)
find $ZFA_DIR/0.94-e85/Cluster*/ | grep sig.tsv$ | \
xargs grep -cvE 'ZFA:0001093' | tr ':' '\t' | \
perl -lane 'if($F[1] > 1){ print $F[0]}' | sed -e 's|/sig\.tsv||; s|^.*/||'
Cluster001
Cluster004
Cluster005
Cluster007
Cluster012
Cluster013
Cluster016
Cluster018
Cluster022
Cluster036
Cluster041
Cluster044
Cluster050
Cluster053
Cluster088
Cluster103
Cluster214
```

Get genes for each significant ZFA term
```bash
# sort ENS gene to ZFIN ID by ZFIN ID
perl -F"\t" -lane 'print join("\t", @F[1,0]);' $ZFA_DIR/danio_rerio_e85_zfin.txt | \
sort -t$'\t' -k1,1 > $ZFA_DIR/zdb-gene-2-ens-gene.e85.txt

cd $ZFA_DIR/0.94-e85/
# variable for the total numbers of genes in the gene set
totalGenes=$( cat $ZFA_DIR/all-zfin-genes.txt | wc -l )
# go through clusters with significant results
for cluster in $( wc -l Cluster*/sig.tsv | \
perl -lane 'if($F[0] > 1){ print $F[1]}' | sed -e 's|/.*||' | grep -v total )
do
echo $cluster
cd $cluster
# make new sig file
echo -e "ID\tp\tp.adjusted\tname\tTotalInTerm\tSignificant\tExpected" > sig.nums.tsv
# variable for the size of this cluster
clusterSize=$(cat sig-zfin-genes.txt | wc -l)
for ZFAterm in $( grep ^ZFA sig.tsv | cut -f1 | grep -v 'ZFA:0001093')
do
# get all genes mapped to the ZFA term
grep $ZFAterm $ZFA_DIR/annotations.txt | cut -f2 | sort -u > $ZFAterm.allGenes.tmp
# count genes
zfaTotal=$(cat $ZFAterm.allGenes.tmp | wc -l)
# get genes mapped to the ZFA term in the cluster
join -t$'\t' $ZFAterm.allGenes.tmp sig-zfin-genes.txt > $ZFAterm.sigGenes.tmp
# count genes
zfaSig=$(cat $ZFAterm.sigGenes.tmp | wc -l)
# print info to file
zfaExpected=$(echo -e "$zfaTotal\t$totalGenes\t$clusterSize" | \
perl -lane 'printf("%.2f\n", $F[0]/$F[1]*$F[2])' )
if [ $( echo $zfaSig '>' $zfaExpected | bc -l ) == 1 ]
then
  grep $ZFAterm sig.tsv | sed -e "s|$|\t$zfaTotal\t$zfaSig\t$zfaExpected|" >> sig.nums.tsv
  grep $ZFAterm sig.tsv | cut -f1,4 | sed -e 's|\t|=|' > $ZFAterm.ens.e85.txt
  echo "Genes for ZFA term ($ZFAterm) = $zfaTotal" >> $ZFAterm.ens.e85.txt
  echo "Genes for ZFA term ($ZFAterm) in Cluster = $zfaSig" >> $ZFAterm.ens.e85.txt
  echo -e "$zfaTotal\t$totalGenes\t$clusterSize" | \
  perl -lane 'printf("%s%.2f\n", "Expected = ", $F[0]/$F[1]*$F[2]);' >> $ZFAterm.ens.e85.txt
  # get Ensembl IDs and info
  join -t$'\t' $ZFAterm.sigGenes.tmp $ZFA_DIR/zdb-gene-2-ens-gene.e85.txt | \
  cut -f2 | sort -u | \
  join -t$'\t' - $ROOT/dataFiles/annotation.txt | \
  cut -f1,7,8 >> $ZFAterm.ens.e85.txt
fi
done
# remove temp files
rm *tmp
# cat together
cat ZFA*.ens.e85.txt > zfa.ens.e85.txt
cd $ZFA_DIR/0.94-e85/
done
```

```bash
# count up lines for different criteria
for cluster in $( wc -l Cluster*/sig.tsv | perl -lane 'if($F[0] > 1){ print $F[1]}' | \
sed -e 's|/.*||' | grep -v total )
do
echo $cluster
perl -F"\t" -lane 'next if $. == 1; if( $F[5] > $F[6] ){ print $_ }' $cluster/sig.nums.tsv | wc -l
perl -F"\t" -lane 'next if $. == 1; next if $F[6] == 0;
if( $F[5]/$F[6] > 1.5 ){ print $_ }' $cluster/sig.nums.tsv | wc -l
done

# which ones are we throwing away with 1.5x cut off
for cluster in $( wc -l Cluster*/sig.tsv | perl -lane 'if($F[0] > 1){ print $F[1]}' | \
sed -e 's|/.*||' | grep -v total )
do
echo $cluster
perl -F"\t" -lane 'next if $. == 1; next if $F[6] == 0;
if( $F[5] > $F[6] && $F[5]/$F[6] <= 1.5 ){ print join("\t", @F, $F[5]/$F[6]); }' $cluster/sig.nums.tsv
done
```

Change enrichment cutoff to 1.3x

```bash
perl -e 'print join("\t", qw{Cluster ID p p.adjusted name TotalInTerm Expected Observed FoldEnrichment} ), "\n";' > zfa.e85.sig.tsv
for cluster in $(find ./ | grep 'sig.nums.tsv' | sed -e 's|/sig.nums.tsv||; s|\./||' | sort )
do
perl -F"\t" -lane 'if( $F[5] > $F[6] && $F[5]/$F[6] >= 1.3 ){ print join("\t", "'$cluster'", @F[0..4,6,5], sprintf("%.3f", $F[5]/$F[6]) ); }' $cluster/sig.nums.tsv >> zfa.e85.sig.tsv
done
```
