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
        } elsif ($line =~ /^[ATCGatgcu-]+$/i && defined $header) {
            $sequences{$header} .= uc($line);
        }
    }
    close $fh;
    return %sequences;
}

# Count nucleotide differences
sub count_nucleotide_differences {
    my ($codon1, $codon2) = @_;
    my $differences = 0;
    for (my $i = 0; $i < 3; $i++) {
        if (substr($codon1, $i, 1) ne substr($codon2, $i, 1)) {
            $differences++;
        }
    }
    return $differences;
}


# Main subroutine to calculate mutation proportions
sub calculate_mutations {
    my ($filename, $ref_id) = @_;
    my %sequences = read_fasta($filename);

    unless (exists $sequences{$ref_id}) {
        die "Reference ID '$ref_id' not found in the provided FASTA file.\n";
    }

    my $ref_seq = $sequences{$ref_id};

    foreach my $header (keys %sequences) {
        next if $header eq $ref_id;  # Skip comparison with itself
        my $seq = $sequences{$header};
        my $total_mutations = 0;
        my $nonsynonymous_mutations = 0;
        my $miss_count = 0;

        for (my $i = 0; $i < length($ref_seq) - 2; $i += 3) {
            next if $i + 3 > length($ref_seq) || $i + 3 > length($seq);  # Ensure within bounds

            my $ref_codon = substr($ref_seq, $i, 3);
            my $codon = substr($seq, $i, 3);

            $miss_count ++ if $ref_codon =~ /-/ || $codon =~ /-/ || length($ref_codon) < 3 || length($codon) < 3;
            next if $ref_codon =~ /-/ || $codon =~ /-/ || length($ref_codon) < 3 || length($codon) < 3;

            my $nuc_diffs = count_nucleotide_differences($ref_codon, $codon);
            $total_mutations += $nuc_diffs;

            if ($nuc_diffs > 0) {
                my $ref_aa = $genetic_code{$ref_codon} // '?';
                my $aa = $genetic_code{$codon} // '?';

                if ($ref_aa ne $aa) {
                    $nonsynonymous_mutations += $nuc_diffs;
                }
            }
        }

        if ($total_mutations > 0) {
            my $proportion = $nonsynonymous_mutations / $total_mutations;
            print "Comparison $ref_id vs $header: Proportion of Nonsynonymous Mutations = $proportion ($nonsynonymous_mutations/$total_mutations) Miss count: $miss_count\n";
        } else {
            print "Comparison $ref_id vs $header: No nucleotide differences found.\n";
        }
    }
}

# File to be analyzed
my $fasta_file = shift;  # Change as necessary

my $reference_id = shift;
chomp $reference_id;

# Run the analysis
calculate_mutations($fasta_file, $reference_id);
