#!/usr/bin/perl
use strict;
use warnings;
use List::Util 'sum';

# Codon table for translating nucleotides to amino acids
my %codon_table = (
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
    'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G',
);

# Ensure two arguments are provided
die "Usage: $0 <alignment_file> <reference_id>\n" unless @ARGV == 2;

my ($filename, $reference_id) = @ARGV;
my %sequences;
my $ref_seq;

# Read the alignment file
open(my $fh, '<', $filename) or die "Cannot open file $filename: $!";
my $seq_id;
while (my $line = <$fh>) {
    chomp $line;
    if ($line =~ /^>(\S+)/) {  # Match sequence identifier after '>'
        $seq_id = $1;  # Store current sequence ID
        $sequences{$seq_id} = '';  # Initialize sequence string
    } elsif (defined $seq_id) {
        # Append to sequence string for the current ID
        $sequences{$seq_id} .= uc($line);
        $ref_seq .= uc($line) if $seq_id eq $reference_id;
    }
}
close $fh;

die "Reference sequence $reference_id not found in file.\n" unless defined $ref_seq;

# Debug: Output total number of sequences and reference sequence
#print "Loaded " . (keys %sequences) . " sequences.\n";
#print "Reference Sequence: $ref_seq\n";

# Analyze sequences
my %allele_counts;
my @positions;

for my $pos (0 .. length($ref_seq) - 3) {
    my $ref_codon = substr($ref_seq, $pos, 3);
    if ($ref_codon =~ /-/ || length($ref_codon) < 3) {  # Skip codons in the reference that contain gaps or are incomplete
        next;
    }
    if (!exists $codon_table{$ref_codon}) {  # Skip if codon is not in the table
        print "Skipping position $pos: Reference codon '$ref_codon' not recognized.\n";
        next;
    }

    my $ref_aa = $codon_table{$ref_codon};
    my %bases;

    foreach my $id (keys %sequences) {
        next if $id eq $reference_id;  # Skip the reference itself
        my $codon = substr($sequences{$id}, $pos, 3);
        if ($codon =~ /-/ || length($codon) < 3) {  # Skip codons that contain gaps or are incomplete
            next;
        }
        if (!exists $codon_table{$codon}) {  # Skip if codon is not in the table
            next;
        }

        my $aa = $codon_table{$codon};
        $bases{$codon}{'aa'} = $aa;
        $bases{$codon}{'count'}++;
    }

    if (keys %bases) {
        $allele_counts{$pos} = \%bases;
        push @positions, $pos;
    }
    else {
        print "No variable sites at position $pos.\n";
    }
}

# Output analysis results
if (@positions) {
    foreach my $pos (@positions) {
        my $ref_codon = substr($ref_seq, $pos, 3);
        my $ref_aa = $codon_table{$ref_codon};
        my $total = sum(map { $_->{'count'} } values %{$allele_counts{$pos}});

        foreach my $codon (keys %{$allele_counts{$pos}}) {
            my $count = $allele_counts{$pos}{$codon}{'count'};
            my $aa = $allele_counts{$pos}{$codon}{'aa'};
            my $change_type = ($aa eq $ref_aa) ? 'synonymous' : 'nonsynonymous';
            my $maf = $count / $total;

            #print "Position: $pos, Codon: $ref_codon->$codon, Amino Acid: $ref_aa->$aa, Change: $change_type, MAF: $maf, ($count / $total)\n";
            print "$pos\t$ref_codon->$codon\t$ref_aa->$aa\t$change_type\t$maf\t($count/$total)\n" if ($maf < 0.5);
        }
    }
} else {
    print "No variable sites found in any sequences compared to reference.\n";
}