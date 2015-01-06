#!/usr/bin/env perl

use Modern::Perl qw|2013|;
use autodie;
use Getopt::Long;
use Statistics::R;
use Parallel::ForkManager;

my ($abundanceProfile,
    $cypherurl,
    $taxDB,
    $refseqDB_nucl,
    $thread,
    $fastqDIR
    );

GetOptions (
    'abundance' => \$abundanceProfile,
    'cypherurl' => \$cypherurl,
    'refseqDB_nucl' => \$refseqDB_nucl,
    'taxDB'     =>   \$taxDB,
    'threads'   =>   \$thread,
    'fastQ dir' =>  \$fastqDIR
)   or die("Error in command line arguments\n");

##################################################
#+------------------------------------------------
#Init
#+------------------------------------------------
##################################################

my $R = Statistics::R->new() ;
my $pm = Parallel::ForkManager->new($threads);

my $fastqDIR =~ s/\/$//;
my $baseDIR = join '/', @path;

my @path = (split /\//, $0);
pop @path;

mkdir 'out' unless -d 'out';
mkdir 'figures' unless -d 'figures';

##################################################
#+------------------------------------------------
#Start
#+------------------------------------------------
##################################################

#0100 Choose genomes
system  "$baseDIR/readSim.0100.chooseGenomes.r $abundanceProfile $cypherurl";
system  "$baseDIR/readSim.0101.contigsXscaffolds.pl out/readSim.0100.woCompleteGenomes_taxidList $cypherurl $refseqDB_nucl $taxDB"
knit(   "$baseDIR/readSim.0102.chooseContigsXScaffolds.Rmd");

#0200
knit(   "$baseDIR/readSim.0201.whoTotrimmed.Rmd");
#$cypherurl

#0300
system  "$baseDIR/readSim.0300.genome_extraction.pl out/readSim.0100.chosen_completeGenomes out/readSim.0102.chosen_scaffolds.txt $refseqDB_nucl"

system  "$baseDIR/readSim.0301.getLength.pl"
knit(   "$baseDIR/readSim.0302.sequenceLengths.Rmd");

# 0400 Run simulation

#Concatenate all of the genus sequences into one
system "cat out/readSim.0300/* > out/readSim.0300.combined.fna"

#problem how am I going to feed all of this into readSim.0401
#Don't know don't care will fix this once xianghui asks if not save my time;
for my $lanes (<$fastqDIR/*>)
{

    $pm->start and next;
    system "$baseDIR/readSim.0401.readsimulation.pl out/readSim.0300.combined.fna $pair1 $pair2 $outputPrefix out/readSim.0100.abundance_NameTaxid.txt";
    $pm->finish;
}
$pm->wait_all_children;

sub knit
{
    my ($script) = @_;
    my $cmds = <<"EOF";
library(knitr);
knit2html($script, template="example/src/htmltemplate.html")
EOF
    $R->run($cmds);
    $R->stop();
}

__END__
CYPHERURL
http://192.168.100.1:7474
