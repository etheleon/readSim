#!/usr/bin/env perl

use Modern::Perl '2013';
use Cwd qw|abs_path cwd|;
use autodie;
use Bio::SeqIO;
use Parallel::ForkManager;

die "USAGE: $0 <path to nrDB> <path.to.taxonomyDB> <list.of.taxa> <number of threads>\n" unless $#ARGV == 3;

my %taxahash;
my %gihash;
my ($nrDIR, $taxonomyDIR, $list, $threads) = @ARGV;
$nrDIR      =~ s/\/$//g;
$taxonomyDIR=~ s/\/$//g;

my @home= split '/', abs_path($0);
pop @home;
my $homeDIR = join '/', @home;
my $pwd = cwd();

my $pm = new Parallel::ForkManager($threads);

unless(-e "$nrDIR/nr.gz")
{
system "gunzip -c $nrDIR/nr.gz > $nrDIR/nr";
}
system "$homeDIR/fasta-splitter.pl --n-parts=$threads $nrDIR/nr";
system "mv *part* $nrDIR/";

#Store condemned taxID
open my $tobedeleted, "<", $list;
while(<$tobedeleted>)
{
    unless($. == 1) #skip header
    {
        chomp;
        my ($parent, $taxa) =  (split(/\t/))[1,2];          #chosen genus underling
        $taxahash{$taxa}=$parent unless $parent =~ /33057/; #Keep Thauera
    }
}
close $tobedeleted;
say "#Finished storing condemned taxa. Total ", scalar keys %taxahash;

#Store condemned GI
-e "$taxonomyDIR/gi_taxid_prot.dmp.gz" ? gi2taxid("zcat $taxonomyDIR/gi_taxid_prot.dmp.gz|") : gi2taxid("$taxonomyDIR/gi_taxid_prot.dmp");

#Trimmed database
mkdir "out/db" unless -d "out/db";
my @splitFiles = <$nrDIR/nr.part*>;
foreach my $nrFile (@splitFiles)
{
#    ##################################################
     $pm->start and next;
#    #-------------------------------------------------

    my $in  = Bio::SeqIO->new(-file => "$nrFile");
    my $fileName = (split /\//, $nrFile)[-1];

    open my $trimmed, ">", "out/db/$fileName"."_trimmed";
    open my $removed, ">", "out/db/$fileName"."_removed";

    while (my $seq = $in->next_seq )
    {
        my $seqid        = $seq->display_id;
        my ($gi, $refid) = (split /\|/, $seqid)[1,3];
        if(!exists $gihash{$gi})
        {
            say $trimmed ">",$seqid;
            say $trimmed $seq->seq;
        }else
        {
            say "$nrFile"."::"."$gi\t$gihash{$gi}";
            say $removed ">",$seqid,"taxon::",$gihash{$gi};
            say $removed $seq->seq;
        }
    }
    close $trimmed;
    close $removed;

#    ------------------------------------------------#
    $pm->finish;                                     #
#    #################################################
}
$pm->wait_all_children;

system "cat out/db/*trimmed > out/db/nrtrimmed";

sub gi2taxid
{
    my ($input) = @_;
    open my $gitax, "<", $input;
    while(<$gitax>)
    {
        chomp;
        my($gi, $taxid) = split(/\t/);
        $gihash{$gi} = $taxahash{$taxid} if (exists $taxahash{$taxid});	#links GI with parentGenus
    }
    close $gitax;
    say "# Finished stored condemned gi. Total ",scalar keys %gihash;
}



__END__
species genus   underling       rank
7       6       7       species
7       6       438753  no rank
19      18      19      species
19      18      338963  no rank
41      40      41      species
41      40      675526  no rank
