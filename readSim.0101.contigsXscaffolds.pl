#!/usr/bin/env perl

use Modern::Perl '2013';
use autodie;
use Bio::SeqIO;
use REST::Neo4p;
use REST::Neo4p::Query;
use List::MoreUtils 'each_arrayref';

die "usage: $0 <genera.names.file> <cypherurl> <refseqDB> <taxDB>\n" unless $#ARGV==3;

my ($genusTAXID,
    $cypher,
    $refseqDB,
    $taxDB
    ) = @ARGV;

$refseqDB =~ s/\/$//;
$taxDB    =~ s/\/$//;

my @genera;
my $missingChildren;
my %gihash;

open my $input, '<', $genusTAXID;
while(<$input>)
{
    chomp;
    push(@genera, $_) unless $. ==1; #skip header
}

#Find child nodes of genera
REST::Neo4p->connect($cypher);
my $querystmt = join "", <DATA>;
findLeaf($_, $querystmt) for @genera;

##################################################
say "Finding all leaf taxa belonging to list";
say "\t", scalar @genera, " genera with ", scalar keys %$missingChildren, " children";
##################################################

#$missingChildren->{9771}{'parent'} = 1000;
#open my $gi2tax, "$taxDB/testing.txt";
open my $gi2tax, "zcat $taxDB/gi_taxid_nucl.dmp.gz |";
while(<$gi2tax>)
{
    chomp;
    my($gi, $taxid) = split /\t/;
    if (exists $missingChildren->{$taxid})
    {
        $gihash{$gi} = $taxid;
    }
}
close $gi2tax;

say
"Stored GIs of missing Genera's children.
Total ", scalar keys %gihash, " GIs";    #should include a count

##################################################
#Output
#+------------------------------------------------
#open my $summedLength   , ">",  "out/readSim.01"; 	#Genus - length
open    my  $output,    ">",  "out/readSim.0101.output.txt";
        my  $seqout     = Bio::SeqIO->new(-file => '>out/readSim.0101.output.fna');
#+------------------------------------------------
##################################################

foreach my $fna (glob "$refseqDB/*")
#do
{
    #my $fna = "testing";
    #my $in  = Bio::SeqIO->new(-file => $fna, format =>'fasta');
    my $in  = Bio::SeqIO->new(-file => "zcat $fna |", format =>'fasta');
    while (my $seq = $in->next_seq )
    {
        my $sequenceID=$seq->display_id;
        my ($gi, $refseqID) = (split /\|/, $sequenceID)[1,3];
        if (exists $gihash{$gi})
        {
            my $taxid = $gihash{$gi};
            push @{$missingChildren->{$taxid}{'gi'}}, $gi;
            push @{$missingChildren->{$taxid}{'refseqID'}}, $refseqID;
            $missingChildren->{$taxid}{'length'}    += $seq->length;
            $seqout->write_seq($seq);
        }
    }
};

#Write Output
say $output join "\t", qw/taxid genus gi refseq combinedLength/;
for my $taxid (keys %$missingChildren)
{
    my %entry = %{$missingChildren->{$taxid}};
    if(exists $entry{'gi'})
    {
        my $genus = $entry{'parent'};
        my $combinedLength = $entry{'length'};
        my $ea = each_arrayref($entry{'gi'}, $entry{'refseqID'});
        while ( my ($giID, $refseqID) = $ea->() )
        {
            say $output join "\t", $taxid, $genus, $giID, $refseqID, $combinedLength;
        }
    }
}

sub findLeaf
{
    my ($taxid, $stmt) = @_;
    my $query = REST::Neo4p::Query->new($stmt,{ nameoftaxa => $taxid});
    $query->execute;
    while (my $row  = $query->fetch){
        my $genus   = $row ->[0];
        my $child   = $row ->[1];
        $missingChildren->{$child}{'parent'} =  $genus;
        $missingChildren->{$genus}{'parent'} =  $genus if !exists $missingChildren->{$genus};
    }
}

#$missingChildren = {
    #21917 => {      #the taxID
        #genus => 21916,      #the parentGenus
        #gi => [1234,5678, 91011],                        #GIs belonging to this taxID
        #refseq => ['NC_XXXX1', 'NC_XXXX2', 'NC_XXXX3']  #refseq belonging to this refseqID
    #}
#};


__DATA__
START
    genus = node:ncbitaxid(taxid={nameoftaxa})
MATCH
    p=(genus)<-[:childof*0..]-(lowerr:species)<-[:childof*0..]-(lowest)
RETURN
    genus.taxid,
    lowest.taxid
