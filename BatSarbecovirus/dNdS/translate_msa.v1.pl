#!/usr/bin/perl
use strict;
use warnings;

# Hash for the genetic code
my %genetic_code = (
    'TTT' => 'F', 'TTC' => 'F', 'TTA' => 'L', 'TTG' => 'L',
    'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
    'ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I', 'ATG' => 'M',
    'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V',
    'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S',
    'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
    'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
    'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',
    'TAT' => 'Y', 'TAC' => 'Y', 'TAA' => '*', 'TAG' => '*',
    'CAT' => 'H', 'CAC' => 'H', 'CAA' => 'Q', 'CAG' => 'Q',
    'AAT' => 'N', 'AAC' => 'N', 'AAA' => 'K', 'AAG' => 'K',
    'GAT' => 'D', 'GAC' => 'D', 'GAA' => 'E', 'GAG' => 'E',
    'TGT' => 'C', 'TGC' => 'C', 'TGA' => '*', 'TGG' => 'W',
    'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',
    'AGT' => 'S', 'AGC' => 'S', 'AGA' => 'R', 'AGG' => 'R',
    'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G'
);

# Read FASTA file
sub read_fasta {
    my ($filename) = @_;
    my %sequences;
    my $header;
    
    open my $fh, '<', $filename or die "Cannot open file $filename: $!";
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^>(.*)/) {
            $header = $1;
            $sequences{$header} = '';
        } elsif ($line =~ /^[ATCGatgc-]+$/i && defined $header) {
            $sequences{$header} .= uc($line);
        }
    }
    close $fh;
    return %sequences;
}

# Translate nucleotide sequence to amino acid sequence, ignoring gaps
sub translate {
    my ($nucleotides) = @_;
    my $amino_acids = '';
    my $len = length($nucleotides);
    
    for (my $i = 0; $i < $len - 2; $i += 3) {
        my $codon = substr($nucleotides, $i, 3);
        if ($codon =~ /^[ATCG]{3}$/) { # Ensure complete and valid codons only
            $amino_acids .= $genetic_code{$codon} // '?'; # Use '?' for unknown codons
        }
    }
    return $amino_acids;
}

# Main subroutine to calculate mutation proportions
sub calculate_mutations {
    my %sequences = read_fasta('alignment.fasta');
    my ($ref_header) = keys %sequences;
    my $ref_seq = $sequences{$ref_header};
    my $ref_protein = translate($ref_seq);

    my $total_synonymous = 0;
    my $total_nonsynonymous = 0;

    foreach my $header (keys %sequences) {
        next if $header eq $ref_header;
        my $seq = $sequences{$header};
        my $min_length = length($ref_seq) < length($seq) ? length($ref_seq) : length($seq);

        for (my $i = 0; $i < $min_length - 2; $i += 3) {
            my $ref_codon = substr($ref_seq, $i, 3);
            my $codon = substr($seq, $i, 3);
            next if $ref_codon =~ /-/ || $codon =~ /-/ || length($ref_codon) < 3 || length($codon) < 3; # Skip gaps and incomplete codons

            if ($ref_codon ne $codon) {
                my $ref_aa = $genetic_code{$ref_codon} // '?';
                my $aa = $genetic_code{$codon} // '?';

                if ($ref_aa eq $aa) {
                    $total_synonymous++;
                } else {
                    $total_nonsynonymous++;
                }
            }
        }
    }

    if (($total_nonsynonymous + $total_synonymous) > 0) {
        my $proportion = $total_nonsynonymous / ($total_nonsynonymous + $total_synonymous);
        print "Proportion of Non-synonymous Mutations: $proportion\n";
    } else {
        print "No valid codon comparisons made.\n";
    }
}

# Run the analysis
calculate_mutations();
