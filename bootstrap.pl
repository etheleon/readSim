#!/usr/bin/env perl

use Modern::Perl qw|2013|;
use autodie;
use Getopt::Long;
use Statistics::R;

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

my $fastqDIR =~ s/\/$//;
my $R = Statistics::R->new() ;
my @path = (split /\//, $0);
pop @path;
my $baseDIR = join '/', @path;

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
open my $combinedFasta ">", 'out/readSim.0300.combined.fna';
foreach my $fastaFile (<out/readSim.0300/*>)
{
    open my $input, "<", $fasta;
    while(<$input>){$. == 2 ? print $combinedFasta "$_\n" : print $combinedFasta $_}
    close $input;
}

#prep
system "baseDIR/readSim.0400.presim.r"

#NOTE: This can be parallelised, have script to do this might want to consult me if you want to do this
for my $lanes (<$fastqDIR/*>)
{
    system "$baseDIR/readSim.0400.simu.shortgun.pl out/readSim.0300.out.fna $lanes $output data/abundanceProfile.txt";
}

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
