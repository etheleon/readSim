#!/usr/bin/env perl

use Modern::Perl '2013';
use autodie;
use Bio::SeqIO;

die "$0 <chosen_completeGenomes> <chosen_woCompleteGenomes> <refseqNuclDIR>\n" unless $#ARGV == 2;

##################################################
#Part1: reads in the refseqIDs of the chosen genomes/scaffolds/contigs
##################################################
my $ref;
my %loc;
my ($completeGenome, $woCompleteGenome, $refseqDIR) = @ARGV;
mkdir "out/readSim.0300" unless -d "out/readSim.0300";

$refseqDIR =~ s/\/$//g;

open my $in, '<', $completeGenome;
while(<$in>)
{
#parentID        taxid   taxon   Chromosomes.RefSeq      total
#233191  926550  Caldilinea      NC_017079.1     25414.868178775
    unless($. == 1)
    {
        chomp;
        my ($genus, $taxid, $refseqid) = (split /\t/)[0,1,3];	#genus (not organism taxid)	$refseq sequence
        my @refseq                     = split ',', $refseqid;
        $ref->{$_}{'genus'}            = $genus for @refseq;
        $ref->{$_}{'taxid'}            = $taxid for @refseq;
    }
}
close $in;

say scalar keys %$ref, " refIDs of complete genomes stored";

#this part i can search from the sim.0102 output
open my $incomplete, '<', $woCompleteGenome;
while(<$incomplete>)
{
#taxid   parentID        Chromosome.RefSeq       total
#485913  363276  NR_042472.1     13662972
    unless ($. == 1)
    {
        chomp;
        my ($taxid,$genus,$refseqid) = (split /\t/)[0,1,2];	#KEY: actualTaxID, refseq
        $ref->{$refseqid}{'genus'}   = $genus;
        $ref->{$refseqid}{'taxid'}   = $taxid;
    }
}
close $incomplete;
say scalar keys %$ref, " refIDs of COMPLETE genomes and scaffold stored";

#mapping file
open my $mapping, ">", "out/readSim.0300.out.mapping";
say $mapping join "\t", qw/Genus Taxon RefSeq Start End/;    #HEADER

#concatenated sequences grouped by genus/leaftaxon
my $out = Bio::SeqIO->new(-file => ">out/readSim.0300.out.fna");

foreach (<"$refseqDIR/*genomic.fna.gz">)
{
    my $in  = Bio::SeqIO->new(-file => "zcat $_ |");
    while (my $seq = $in->next_seq ){
        my $seqid   = $seq->display_id;
        my ($refid) = (split /\|/, $seqid)[3];
        if(exists $ref->{$refid})
        {
            my %localref = %{$ref->{$refid}};
            $seqid       = join '|', 'genus', $localref{'genus'}, 'leaftaxa', $localref{'taxid'},$seqid;
            my $sequence = $seq->seq();
            my $length   = $seq->length;
            $seq->display_id($seqid);
            $out->write_seq($seq);
            joinNsplit($localref{'genus'}, $localref{'taxid'}, $sequence, $length, $refid);
            delete $ref->{$refid};
        }
    }
    $in->close;
}
$out->close;


#Which refseqs have not been found in the nr db removed
open my $report, '>', 'out/readSim.0300.out.notRemoved';
    say $report "#Following refseqIDs were not found";
    say $report join "\t", qw/refseq taxid/;
    say $report join("\t", $_, $ref->{$_}) for keys %$ref;
close $report;

sub joinNsplit  #concatenates the sequences together to single fastafile
{
    my ($genus, $taxid, $sequence, $length, $refid) = @_;
    my $location = "out/readSim.0300/$genus.fna";
    if (-e $location)
    {
        open my $taxaoutput, '>>', $location;
        print $taxaoutput $sequence;
        close $taxaoutput;

        my $start = $loc{$genus};
        my $end   = $start + $length - 1;
        say $mapping join "\t", $genus, $taxid, $refid, $start, $end;
    }else
    {
        $loc{$genus} = 0;
        open my $taxaoutput, '>', $location;
        say $taxaoutput join "|", '>genus', $genus, 'leaftaxa', $taxid;
        print $taxaoutput $sequence;
        close $taxaoutput;

        my $start = 0;
        my $end   = $start + $length - 1;
        say $mapping join "\t", $genus, $taxid, $refid, $start, $end;
    }
}
