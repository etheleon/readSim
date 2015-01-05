#!/usr/bin/env perl

use Modern::Perl '2013';
use Bio::SeqIO;
use autodie;
my $seqLength;

my $in  = Bio::SeqIO->new(-file => "out/readSim.0300.out.fna");
while(my $seq = $in->next_seq)
{
    my $id = $seq->display_id;
    my $length = $seq->length;
    my ($genus, $taxID) = (split /\|/, $id)[1,3];
    #say join "\t", $genus, $taxID;
    $seqLength->{$genus}{'length'} += $length;
    push @{$seqLength->{$genus}{'leaftaxa'}}, $taxID;
}

open my $out, ">", "out/readSim.0301.lengthtable.txt";
say $out join "\t", qw/genus length/;
say $out join("\t", $_, $seqLength->{$_}{'length'}) for keys %$seqLength;
