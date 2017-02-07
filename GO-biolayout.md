# Go Enrichment

topGO needs a file that maps between Ensembl gene IDs and GO terms. It is a tab-separated file with Ensembl gene IDs in the first column and a comma-separated  list of associated GO terms in the second column.

This can be made from the Ensembl annotation file (danio_rerio_e85_go.txt).
It can be downloaded from the [topgo-wrapper](https://github.com/iansealy/topgo-wrapper/blob/master/data/danio_rerio_e85_go.txt) github repository. This uses version 85. If newer annotation is required that isn't present in the repository, it can be produced using the get_ensembl_go_terms.pl script from the same repository (https://github.com/iansealy/topgo-wrapper).

Make the mapping files for each GO domain.
```bash
mkdir topgo
for domain in cellular_component molecular_function biological_process
do grep $domain dataFiles/danio_rerio_e85_go.txt | \
perl -F"\t" -lane 'BEGIN{%go_terms_for;}
{push @{ $go_terms_for{$F[0]} }, $F[1];}
END{foreach my $gene_id (sort keys %go_terms_for ){
print join("\t", $gene_id, join(",", @{ $go_terms_for{$gene_id} } ) ); } }' > topgo/$domain.e85.map
done
```

## Run topGO

topGO was run using a custom R script (topGO-biolayout-clusters.R). The clusterName option to the script runs the gene enrichment on the specified cluster. Otherwise it runs it on each cluster sequentially.
e.g.
```bash
mkdir topgo/0.94-e85/
# single cluster
Rscript topGO-biolayout-clusters.R \
--directory=topgo --OutputDirectory=0.94-e85 --clusterName=Cluster001 \
dataFiles/biolayout.GRCz10.publication_change_with_stage_r-0.94.cluster.info.txt

# or run all together
Rscript topGO-biolayout-clusters.R \
--directory=topgo --OutputDirectory=0.94-e85 \
dataFiles/biolayout.GRCz10.publication_change_with_stage_r-0.94.cluster.info.txt
```

## Pull out significantly enriched genes

A perl script was used to select which GO terms were enriched in which clusters.

```bash
dir='topgo/0.94-e85'
adjustedP=$(echo "0.05/252" | bc -l | sed -e 's|^|0|')
for cluster in $( seq -w 1 252 | sed -e 's|^|Cluster|' )
do
grep $cluster dataFiles/biolayout.GRCz10.publication_change_with_stage_r-0.94.cluster.info.txt \
> $dir/$cluster.gene-list.tsv
for ontology in BP MF CC
do
perl parse_GO_results.pl \
--enriched --output_file $dir/$cluster.GO.$ontology.sig.tsv \
$dir/$cluster.GO.$ontology.all.tsv \
$dir/$cluster.gene-list.tsv > $dir/$cluster.GO.$ontology.enriched.sig.tsv
done
done
```

join cluster results together
```bash
dir='topgo/0.94-e85'
for cluster in $( seq -w 1 252 | sed -e 's|^|Cluster|' )
do
for ontology in BP MF CC
do
domain=$( perl -e 'BEGIN{
%domains = ( BP => "biological_process",
              MF => "molecular_function",
              CC => "cellular_component" ); }
{ print $domains{'$ontology'}; }' )
sed -E "s|^|$cluster;$domain;|" $dir/$cluster.GO.$ontology.enriched.sig.tsv | \
tr ';' '\t' | sort -t$'\t' -k8,8nr
done
done | perl -F"\t" -lane 'BEGIN{ print join("\t", qw{Cluster GOTermID Description Domain AdjustedpValue Expected Observed log2FoldEnrichment})}
{print join("\t", @F[0,2,3,1,4..7]); }' > $dir/topgo-0.94-sigResults.tsv
```
