#!/usr/bin/env perl

# PODNAME: parse_GO_results.pl
# ABSTRACT: takes a file from the outputof out topGO analysis and an gene list file and creates a more readable output

use warnings;
use strict;
use Getopt::Long;
use autodie;
use Pod::Usage;
use Readonly;
use DateTime;

# get options
my %options;
my $SIG_LEVEL;
my $FOLD_ENRICHMENT;
my $MIN_GENES;
get_and_check_options();

# open gene list
my $gene_list_file = $ARGV[1];
my %gene_info;
if( !$gene_list_file ){
    pod2usage('A gene list file must be supplied!');
}
else{
    open my $gene_fh, '<', $gene_list_file;
    while(<$gene_fh>){
        chomp;
        my ( $gene_name, $gene_id, undef, ) = split /\t/, $_;
        if( exists $gene_info{$gene_id} ){
            die 'Found more than one entry in gene list for gene, ', $gene_id, "\n";
        }
        else{
            $gene_info{$gene_id} = $gene_name;
        }
    }
}

my $topGO_file = $ARGV[0];
if( !$topGO_file ){
    pod2usage('A topGO output file must be supplied!');
}
else{
    my %info_for;
    open my $topGO_fh, '<', $topGO_file;
    while(<$topGO_fh>){
        chomp;
        my ( $GOTerm, $description, $annotated,
            $significant, $expected, $pvalue, $gene_ids ) = split /\t/, $_;
        next if $GOTerm eq 'GO.ID'; # header line
        # check p value
        next if $pvalue eq "NA";
        $pvalue =~ s/< //xms;
        next if $pvalue >= $SIG_LEVEL;
        
        # check whether enriched or depleted
        my $consider = $options{enriched} && $options{depleted} ? 1
            : $options{enriched} && $significant > $expected ? 1
            : $options{depleted} && $significant < $expected;
        
        # get gene names
        if( $consider ){
            my $num_genes = 0;
            $info_for{$GOTerm} = {
                description => $description,
                expected => $expected,
                significant => $significant,
                pvalue => $pvalue,
            };
            foreach my $gene_id ( split /,/, $gene_ids ){
                if( !exists $gene_info{$gene_id} ){
                    if( $options{debug} ){
                        warn "MISSING: ", $gene_id, "\tCouldn't find an entry in the gene list file.\n";
                    }
                }
                else{
                    $num_genes++;
                    push @{ $info_for{$GOTerm}{genes} }, join(" = ", $gene_id, $gene_info{$gene_id} );
                }
            }
            if( $num_genes != $significant ){
                warn "GO TERM: ", $GOTerm, "\n",
                    "Numbers of gene ids matched does not equal the number of significant genes\n",
                    join(": ", 'Significant', $significant, ), "\n",
                    join(": ", 'Matched', $num_genes, ), "\n";
            }
        }
    }
    #output
    open my $out_fh, '>', $options{output_file};
    foreach my $term ( sort { $info_for{$a}{pvalue} <=> $info_for{$b}{pvalue} } keys %info_for ){
        print {$out_fh} "GOTerm: ", $term, "\n";
        print {$out_fh} "  description: ", $info_for{$term}{description}, "\n";
        print {$out_fh} "  p-value: ", $info_for{$term}{pvalue}, "\n";
        print {$out_fh} "  expected: ", $info_for{$term}{expected}, "\n";
        print {$out_fh} "  significant: ", $info_for{$term}{significant}, "\n";
        print {$out_fh} "  genes: \n";
        foreach my $gene_info ( sort @{ $info_for{$term}{genes} } ){
            print {$out_fh} "    - $gene_info", "\n";
        }
        
        if( $info_for{$term}{significant}/$info_for{$term}{expected} >= $FOLD_ENRICHMENT &&
            $info_for{$term}{significant} >= $MIN_GENES ){
            print join("\t", $term, $info_for{$term}{description}, $info_for{$term}{pvalue},
                       $info_for{$term}{expected}, $info_for{$term}{significant},
                       $info_for{$term}{significant}/$info_for{$term}{expected}, ), "\n";
        }
    }
    close( $out_fh );
}

sub get_and_check_options {
    
    GetOptions(
        \%options,
        'enriched',
        'depleted',
        'sig_level=f',
        'fold_enrichment=f',
        'min_genes=i',
        'output_file=s',
        'help',
        'man',
        'debug+',
        'verbose',
    ) or pod2usage(2);
    
    # Documentation
    if( $options{help} ) {
        pod2usage( -verbose => 0, -exitval => 1, );
    }
    elsif( $options{man} ) {
        pod2usage( -verbose => 2 );
    }
    
    if( !$options{enriched} && !$options{depleted} ){
        $options{enriched} = 1;
        $options{depleted} = 1;
    }
    
    Readonly $SIG_LEVEL => defined $options{sig_level}
        ? $options{sig_level}
        : 0.05;
        
    Readonly $FOLD_ENRICHMENT => defined $options{fold_enrichment}
        ? $options{fold_enrichment}
        : 1.5;
    
    Readonly $MIN_GENES => defined $options{min_genes}
        ? $options{min_genes}
        : 5;
    
    if( !$options{output_file} ){
        my $dt = DateTime->now;
        $options{output_file} = join('.', $dt->ymd, 'GO', 'sig', 'yml');
    }
    
    $options{debug} = !defined $options{debug} ? 0 : $options{debug};
    print "Settings:\n", map { join(' - ', $_, defined $options{$_} ? $options{$_} : 'off'),"\n" } sort keys %options if $options{verbose};
}

__END__

=pod

=head1 NAME

parse_GO_results.pl

=head1 DESCRIPTION

Takes the tsv file output from topGO and a gene list file with Ensembl ID to
gene name mappings and creates an output showing which genes are the ones
enriched/depleted.

=cut

=head1 SYNOPSIS

    parse_GO_results.pl [options] topGO_output gene_list_file
        --enriched              only consider terms that are enriched
        --depleted              only consider terms that are depleted
        --sig_level             significance level for reporting
        --fold_enrichment       fold enrichment cut-off
        --min_genes             minimum number of genes cut-off
        --output_file           name for the outputfile
        --help                  print this help message
        --man                   print the manual page
        --debug                 print debugging information
        --verbose               turn on verbose output


=head1 ARGUMENTS

=over

=item topGO output file

Output file from our topGO analsysis. Tab-separated

=item gene list

File of the genes that were tested for enrichment/depletion. Tab-separated.
Must contain these columns in this order:

=over

=item name - Gene Name

=item gene_id - Ensembl IDs

=back

=back

=head1 OPTIONS


=over

=item B<--enriched>

Only consider terms where there is an enrichment of genes for that term.
Default is both enrichment and depletion

=item B<--depleted>

Only consider terms where there is a depletion of genes for that term.
Default is both enrichment and depletion

=item B<--sig_level>

Only report term where the pvalue is below this level.
Default = 0.05

=item B<--fold_enrichment>

Only report term where the fold enrichment is above this level.
Default = 1.3

=item B<--min_genes>

Only report term where the number of observed genes annotated to the term is equal to or greater than this level.
Default = 5

=item B<--output_file>

Name for the output file. Defaults to the current date - DATE.GO.tsv

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=head1 DEPENDENCIES

None

=head1 AUTHOR

=over 4

=item *

Richard White <richard.white@sanger.ac.uk>

=back

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2016 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut