#!/usr/bin/env perl

use Modern::Perl qw|2014|;
use Getopt::Long;
use Statistics::R;
my $R = Statistics::R->new() ;
my $abundanceProfile= ''; # option variable with default value (false)
my $cypherurl;
my $taxDB;
my $refseqDB_nucl;
my $thread;

GetOptions (
    'abundance' => \$abundanceProfile,
    'cypherurl' => \$cypherurl,
    'refseqDB_nucl' => \$refseqDB_nucl,
    'taxDB'     =>   \$taxDB,
    'threads'   =>   \$thread
)   or die("Error in command line arguments\n");



#This script will control the execution of the whole pipeline
my @path = (split '/', $0);
pop @path;
my $baseDIR = join "/", @path;

mkdir 'out' unless -d 'out';
mkdir 'figures' unless -d 'figures';

#0100 Choose genomes
system $baseDIR."/readSim.0100.chooseGenomes.r $abundanceProfile $cypherurl";
system $baseDIR."/readSim.0101.contigsXscaffolds.pl out/readSim.0100.woCompleteGenomes_taxidList $cypherurl $refseqDB_nucl $taxDB"

my $cmds = <<"EOF";
library(knitr);
knit2html("$file.Rmd")
EOF
$R->run($cmds);
$R->stop();


# 0400 Run simulation

__END__
CYPHERURL
http://192.168.100.1:7474
