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
        } elsif ($line =~ /(\w)+/ && defined $header) {
            chomp $line;
            $sequences{$header} .= uc($line);
        }
    }
    close $fh;
    return %sequences;
}

# Evaluate individual nucleotide changes within codons
sub evaluate_nucleotide_changes {
    my ($codon1, $codon2) = @_;
    my $changes_causing_aa_change = 0;
    for (my $i = 0; $i < 3; $i++) {
        my $mutated_codon = $codon1;
        substr($mutated_codon, $i, 1, substr($codon2, $i, 1));  # Create mutated codon
        if ($genetic_code{$mutated_codon} && $genetic_code{$codon1} &&
            $genetic_code{$mutated_codon} ne $genetic_code{$codon1}) {
            $changes_causing_aa_change++;
        }
    }
    return $changes_causing_aa_change;
}


# Calculate mutation proportions
sub calculate_mutations {
    my ($filename, $ref_id) = @_;
    my %sequences = read_fasta($filename);
    
    unless (exists $sequences{$ref_id}) {
        die "Reference ID '$ref_id' not found in the provided FASTA file.\n";
    }

    foreach my $header (keys %sequences) {
      #  print $header, "\t", length ($sequences{$header}), "\n";
    }

    my $ref_seq = $sequences{$ref_id};

    foreach my $header (keys %sequences) {
        next if $header eq $ref_id;  # Skip comparison with itself
        my $seq = $sequences{$header};
        my $total_nucleotide_changes = 0;
        my $aa_changes_due_to_nucleotide_changes = 0;
        my $total_length = 0;

        for (my $i = 0; $i <= length($ref_seq) - 3; $i += 3) {
            next if $i + 3 > length($ref_seq) || $i + 3 > length($seq);

            my $ref_codon = substr($ref_seq, $i, 3);
            my $codon = substr($seq, $i, 3);
            next if $ref_codon =~ /-/ || $codon =~ /-/ || length($ref_codon) < 3 || length($codon) < 3;

            $total_length += 3;# unless $ref_codon =~ /-/ || length($ref_codon) < 3 || length($codon) < 3;
            
            for (my $j = 0; $j < 3; $j++) {
                if (substr($ref_codon, $j, 1) ne substr($codon, $j, 1)) {
                    $total_nucleotide_changes++;
                    if (evaluate_nucleotide_changes($ref_codon, $codon) > 0) {
                        $aa_changes_due_to_nucleotide_changes++;
                    }
                }
            }
        }

        if ($total_nucleotide_changes >= 0) {
            $total_nucleotide_changes ++;
            my $proportion = $aa_changes_due_to_nucleotide_changes / $total_nucleotide_changes;
            #print "Comparison $ref_id vs $header: Proportion of Nucleotide Changes Causing Amino Acid Changes = $proportion; Nucleotide changes: $total_nucleotide_changes among $total_length\n";
            print "$ref_id\t$header\t$total_length\t$total_nucleotide_changes\t$proportion\n";
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
