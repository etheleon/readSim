#!/usr/bin/env perl

use Modern::Perl qw|2013|;
use autodie;
use Statistics::R;
use JSON::Parse 'json_file_to_perl';
use Parallel::ForkManager;
use Cwd 'abs_path';

die "$0 <config.json> <threads> <restart>\n" unless $#ARGV == 1;

##################################################
#+------------------------------------------------
#Init
#+------------------------------------------------
##################################################
#
my ($configFile, $thread, $restart) = @ARGV;
my $config = json_file_to_perl($configFile);
my $pm = Parallel::ForkManager->new($thread);

my @path = (split /\//, abs_path($0));
pop @path;
my $baseDIR = join '/', @path;
say "## Base directory set at $baseDIR";

mkdir 'out' unless -d 'out';
mkdir 'figures' unless -d 'figures';
mkdir 'data' unless -d 'data';

if($restart)
{
    system "rm out/readSim.0300/* && rm out/readSim.0300*";
}

##################################################
#+------------------------------------------------
#Start
#+------------------------------------------------
##################################################

#0100 Choose genomes
knit    ("$baseDIR/readSim.0100.chooseGenomes.Rmd");
system  "$baseDIR/readSim.0101.contigsXscaffolds.pl out/readSim.0100.woCompleteGenomes_taxidList $config->{cypherurlperl} $config->{refseq} $config->{ncbiTaxonomy}";
knit    ("$baseDIR/readSim.0102.chooseContigsXScaffolds.Rmd");

##0200 Protein Database trimming
knit    ("$baseDIR/readSim.0201.whoTotrimmed.Rmd");
system  "$baseDIR/readSim.0202.trimDB.pl $config->{nr} $config->{ncbiTaxonomy} out/readSim.0201.output.txt $thread";

##0300 Genome filtering
system  "$baseDIR/readSim.0300.genome_extraction.pl out/readSim.0100.chosen_completeGenomes out/readSim.0102.chosen_scaffolds.txt $config->{refseq}";
system  "$baseDIR/readSim.0301.getLength.pl";
knit    ("$baseDIR/readSim.0302.sequenceLengths.Rmd");

# 0400 Run simulation
system "cat out/readSim.0300/* > out/readSim.0300.combined.fna";

foreach my $lane (@{$config->{fastq}{files}})
{
    $pm->start and next;
    ##################################################
    my $pair1 = $lane->[0];
    my $pair2 = $lane->[1];
    my $outputPrefix = $lane->[2];
    my $fastqDIR = $config->{fastq}{dir};
    $fastqDIR =~ s/\/$//;
    system "$baseDIR/readSim.0400.readsimulation.pl out/readSim.0300.combined.fna $fastqDIR/$pair1 $fastqDIR/$pair2 $outputPrefix out/readSim.0100.abundance_NameTaxid.txt 64";
    ##################################################
    $pm->finish;
};
$pm->wait_all_children;

sub knit
{
    my $R = Statistics::R->new() ;
    my ($script) = @_;
    my $cmds = <<"EOF";
library(knitr);
knit2html(input="$script", template="src/htmltemplate.html")
EOF
   $R->run($cmds);
   $R->stop();
}
