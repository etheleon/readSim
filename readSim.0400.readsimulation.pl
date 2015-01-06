#!/usr/bin/env perl

use Modern::Perl '2013';
use autodie;

die "usage: $0 <chosenGenomes.fna> <template.pair1.fq> <template.pair1.fq> <outputfileprefix> <abundanceInfo>
(indel-rate: proportion of indel over all errors; default=0.1)\n" unless $#ARGV==4;

my ($chosenGenomes, $templateFQ1,$templateFQ2, $outputFile, $abundanceInfo) = @ARGV;

##################################################
# Part0: Init
say "Initializing ...";
##################################################

my (
    $seq,   #hash ref
    %score,
    %ntrate,#nucleotide::the frequencies
    %globalnt,#nucleotide::frequency
);

my  $totalAbu;                                       # summed abundance
my  $id = 0;                                         # For adding to fastq header
my  $indelrate=0.01;                                 # the INDEl error rate:percentage of errors to be indel
my  $outputfile = $outputFile;
    $outputfile =~ s/.+\///; # the outputfile name

#Phred error probabilities
$score{$_} = phred($_) for 0..100;                  #Probability

##################################################
# Part1: read genomes
say "Reading Genomes ...";
##################################################

readFasta($chosenGenomes);
$seq->{$_}{length} = length($seq->{$_}{sequence}) for keys %$seq;
say "\tStored ",scalar keys %$seq, " sequences";

####################################################
## Part2: abundance
say "Reading abundances.. ";
###################################################

open my $abundance, '<', $abundanceInfo;
#UNFINISHED
#out/sim.0101.out2.txt	#this will read in full 413 taxa (bacteria and archea)
while(<$abundance>)
{
    unless ($. == 1)
    {
        chomp;
        my ($abu, $genus) = (split /\t/)[1,3];
        if (exists $seq->{$genus})
        {
            #this will only allow abundance to be stored if a sequence of that taxid exists
            $seq->{$genus}{abundance} = $abu;
        }
    }
}

##abundance summary
$totalAbu += $seq->{$_}{abundance}                              for keys %$seq;	#total reads per million abundance
$seq->{$_}{abundance} = ($seq->{$_}{abundance}/$totalAbu)       for keys %$seq;

#calculate Nucleotide Frequency
say "\tCalculating Nucleotide frequencies";
ntfreq($_) for keys %$seq;	#will push genus taxid into ntfreq

say "\tNucleotide frequencies";
say "\t\t$_: $globalnt{$_}" for keys %globalnt;

###taxid - range creation:: TAXID: (lower, upper)
#from the biggest to the smallest
foreach my $taxid (sort { $seq->{$b}{abundance} <=> $seq->{$a}{abundance} } keys %$seq)
{
    state $start = 0;
    my $abundance = $seq->{$taxid}{abundance};
    my $end = $start + $abundance;
    $seq->{$taxid}{start} = $start;
    $seq->{$taxid}{end}   = $end;
    $start = $end;
}

###################################################
## Simulate Shotgun sequence
say "Reading quality scores & simulating ...";
###################################################

open my $fastqOutputONE,    '>',    $outputFile."_1.fq";
open my $fastqOutputTWO,    '>',    $outputFile."_2.fq";

open my $fastqONE,          '<', $templateFQ1;
open my $fastqTWO,          '<', $templateFQ2;

while (!eof($fastqONE) and !eof($fastqTWO))
{
    #scrapping quality score only;
    <$fastqONE>;<$fastqONE>;<$fastqONE>;	#h1, sequence, h2
    <$fastqTWO>;<$fastqTWO>;<$fastqTWO>;	#h1, sequence, h2

#Process Quality

    my $qualONE =  <$fastqONE>;
    chomp $qualONE;

    my $qualTWO =  <$fastqTWO>;
    chomp $qualTWO;

    my $readlength = length($qualONE);	#template length
    my $qualitystringONE = processQual($qualONE);
    my $qualitystringTWO = processQual($qualTWO);

##Choosing taxa and loc to pluck sequence out from

    my $taxaofchoice = choosetaxa();
    processPairReads($taxaofchoice, $readlength, $qualitystringONE, $qualitystringTWO, $qualONE,$qualTWO);
}

say "All done. Please check $outputFile 1.fq/2.fq for results";

####################################################################################################
#Functions
####################################################################################################
#Calculates error-probability of mutation
sub phred
{
    my ($score) = @_;
    return 10**($score/(-10));
}

sub ntfreq
{
#count the number of ATCGs
    my ($taxid) = @_;
    my %ntrate = (
        'a' => ($seq->{$taxid}{'sequence'}=~tr/aA/AA/),
        't' => ($seq->{$taxid}{'sequence'}=~tr/tT/TT/),
        'g' => ($seq->{$taxid}{'sequence'}=~tr/gG/GG/),
        'c' => ($seq->{$taxid}{'sequence'}=~tr/cC/CC/)
    );
    $globalnt{$_} += ($ntrate{$_} / $seq->{$taxid}{length}) * $seq->{$taxid}{abundance} for (keys %ntrate); #correct one
}
##mutates the sequence based on frequency
sub mutate
{
    my ($readnt,$qual,$readLength) = @_;
    my @readnt = split(//, $readnt);
    my @qual = split(/\t/, $qual);		#stores the probability
    my $loc = 1;				#this is loc of buffered seqeunce in case of deletion event

    my @outputsequence;	#store sequence

    for(my $i=0; $i < $readLength; $i++) {
        my $qualityscore = $qual[$i];
        if(rand()<$qualityscore){
            #MUTATION
            #INDEL EVENT #####################################
            if(rand()<$indelrate)	{
                ##################################################
                if(int(rand(2))){ 	#INSERTION
                    my $extra = rand_nt();
                    $extra.= $readnt[$loc];
                    push @outputsequence, $extra;
                    $loc++
                }else{			#DELETION
                    $loc++;
                    push @outputsequence, $readnt[$loc];
                    $loc++
                }
                ##################################################

                #SUBSTITUTION EVENT ##############################
            }else	{
                ##################################################
                my $substitutedNT = rand_nt();
                push @outputsequence, $substitutedNT;
                $loc++
            }
            ##################################################
        }else{
            #NO MUTATION
            push @outputsequence, $readnt[$loc];
            $loc++;
        }
    }
    my $finalseq = join '', @outputsequence;
    return $finalseq;
}

##chooses the random ATGC to change
sub rand_nt
{
    my $rand = rand;
    return 'a' if($rand<$globalnt{'a'});
    return 't' if($rand<$globalnt{'a'}+$globalnt{'t'});
    return 'g' if($rand<$globalnt{'a'}+$globalnt{'t'}+$globalnt{'g'});
    return 'c' if($rand<$globalnt{'a'}+$globalnt{'t'}+$globalnt{'g'}+$globalnt{'c'});
    return 'n';
}

#chooses a random taxa
sub choosetaxa
{
    my $rand = rand();
    foreach (keys %$seq)
    {
        if($rand > $seq->{$_}{start} & $rand <= $seq->{$_}{end})
        {
            return $_;
        }
    }
}

sub processQual
{
    my ($qual) = @_;
    my @qual = map { ord($_) - 33 } split('',$qual);	#convert ASCII to indexNum
    my $qualitystring = join "\t", map {$score{$_}} @qual;
    return $qualitystring;
}

sub writeSequence
{
    my ($id, $sequence, $qual, $nameOfSequence, $start, $readLength,
        $outputfile, $filehandle, $pair) = @_;

    my $header = "READ_$id|taxID|$nameOfSequence|loc|$start-".$start+$readLength."|output|$outputfile/$pair";
        say $filehandle '@'.$header;
        say $filehandle $sequence;
        say $filehandle "+";
        say $filehandle $qual;
}

sub processPairReads
{
    my ($taxaofchoice,$readLength,$qualONE,$qualTWO,$asciiQualONE, $asciiQualTWO) = @_;

    my $insertSize = int(random_normal(150, 5)); #xc suggested to set min size in case probability kills you, not implemented
    my $genomelocation = int(rand($seq->{$taxaofchoice}{length}-$insertSize+1));
    ##################################################
    #read one
    ##################################################
    my $readnt = substr($seq->{$taxaofchoice}{sequence}, $genomelocation, $insertSize);
    my $newsequence = mutate($readnt, $qualONE,$readLength);

            writeSequence($id,$newsequence,$asciiQualONE,$taxaofchoice,$genomelocation,
            $readLength,$outputfile,$fastqOutputONE,1);
    ##################################################
    #read two
    ##################################################

    #a. reverse
    $readnt = scalar reverse $readnt;
    #b. complement
    $readnt =~ tr/ATGC/TACG/;
    $newsequence = mutate($readnt, $qualTWO, $readLength);

            writeSequence($id, $newsequence,$asciiQualTWO,$taxaofchoice,$genomelocation,
            $readLength,$outputfile,$fastqOutputTWO, 2);
    $id++;
}

sub readFasta
{
    my ($fastaFile) = @_;
    my $taxid;

    open my $in, '<', $fastaFile;
    while (<$in>)
    {
        my $genus;
        my $child;
        next if(m/^\s*$/); #remove blank lines
        if(m/>\s*?(\S+)/)
        {
            #>genus|1007|leaftaxa|984262
            chomp;
            ($genus, $child) = (split /\|/, $1)[1];  #the taxid of the sequence/genome
            $seq->{$genus}++;
        }else{
            chomp;
            $seq->{$genus}{'sequence'} .= $_;
        }
    }
    close $in;
}

sub gaussian_rand {
    my ($u1, $u2);  # uniformly distributed random numbers
    my $w;          # variance, then a weight
    my ($g1, $g2);  # gaussian-distributed numbers

    do {
        $u1 = 2 * rand() - 1;
        $u2 = 2 * rand() - 1;
        $w = $u1*$u1 + $u2*$u2;
    } while ( $w >= 1 );

    $w = sqrt( (-2 * log($w))  / $w );
    $g2 = $u1 * $w;
    $g1 = $u2 * $w;
    # return both if wanted, else just one
#    return wantarray ? ($g1, $g2) : $g1;
    return $g1;
}

sub random_normal
{
    my ($mean, $sd) = @_;
    my $num = gaussian_rand();
    my $newnum = ($num * $sd) + $mean;
}
