#!/usr/bin/env perl

use Modern::Perl '2013';
use autodie;

die "usage: $0 <chosenGenomes.fna> <template.pair1.fq> <template.pair1.fq> <outputfileprefix> <abundanceInfo>
    (indel-rate: proportion of indel over all errors; default=0.1)\n" unless $#ARGV==4;

my ($chosenGenomes, $templateFQ1,$templateFQ2, $outputFile, $abundanceInfo) = @ARGV;

##################################################
say "## Initializing ...";
##################################################

my $seq = {};
my %score;     #phred score
my %globalnt;  #nucleotide::frequency

my  $totalAbu;                                       # summed abundance
my  $id = 0;                                         # For adding to fastq header
my  $indelrate=0.01;                                 # the INDEl error rate:percentage of errors to be indel
    $outputFile =~ s/.+\///;                         # the outputfile name

#Phred error probabilities
$score{$_} = phred($_) for 0..100;                  #Probability

##################################################
# Part1: read genomes
say "## Reading Genomes/Scaffolds...";
##################################################

readFasta($chosenGenomes);
$seq->{$_}{length} = length($seq->{$_}{sequence}) for (keys %$seq);
say "##\tStored ",scalar keys %$seq, " sequences";

####################################################
say "## Reading abundances.. ";
###################################################

open my $abundance, '<', $abundanceInfo;
while(<$abundance>)
{
    #taxon   total   taxid
    #Acaryochloris   667.462796960323        155977
    unless ($. == 1)
    {
        chomp;
        my ($abu, $genus) = (split /\t/)[1,2];
        if (exists $seq->{$genus})
        {
            $seq->{$genus}{abundance} = $abu;
        }
    }
}

#Abundance as a percentage
$totalAbu += $seq->{$_}{abundance}                              for keys %$seq;	#total reads per million abundance
$seq->{$_}{abundance} = ($seq->{$_}{abundance}/$totalAbu)       for keys %$seq;

say "##\tCalculating Global nt frequencies";
ntfreq($_) for keys %$seq;

say "##\t\t$_: $globalnt{$_}" for keys %globalnt;
state $total += $globalnt{$_} for keys %globalnt;
say "##\tSums to: $total";

##taxid - range creation:: TAXID: (lower, upper)
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

####################################################
say "## Simulating reads NOW...";
####################################################

open my $fastqOutputONE,    '>',    "out/".$outputFile."_1.fq";
open my $fastqOutputTWO,    '>',    "out/".$outputFile."_2.fq";

open my $fastqONE,          '<', $templateFQ1;
open my $fastqTWO,          '<', $templateFQ2;

while (!eof($fastqONE) and !eof($fastqTWO))
{
    #ignoring rest of fastq, scrapping quality score only;
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

    my $chosenTaxon= choosetaxa();
    processPairReads($chosenTaxon, $readlength, $qualitystringONE, $qualitystringTWO, $qualONE,$qualTWO);
}

say '## ^Completed^';
say "## $id reads simulated";
say "## FastQ files at out/",$outputFile,"_1.fq and out/",$outputFile, "_2.fq";

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
    my @readnt = split '', $readnt;
    my @qual   = split /\t/, $qual;       #stores the probability
    my $loc    = 1;                        #this is loc of buffered seqeunce in case of deletion event

    my @outputsequence;

    for(my $i=0; $i < $readLength; $i++) {
        my $qualityscore = $qual[$i];
        if(rand()<$qualityscore){
            #MUTATION
            if(rand()<$indelrate)
            {
                #Indel event##################################
                if(int(rand(2)))
                {   #INSERTION
                    my $extra = rand_nt();
                    $extra.= $readnt[$loc];
                    push @outputsequence, $extra;
                    $loc++
                }else
                {   #DELETION
                    $loc++;
                    push @outputsequence, $readnt[$loc];
                    $loc++
                }
                ##############################################
            }else
                #Substitution event###########################
            {
                my $substitutedNT = rand_nt();
                push @outputsequence, $substitutedNT;
                $loc++
            }
                #############################################
        }else
        {
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
    foreach my $taxid (keys %$seq)
    {
        if($rand > $seq->{$taxid}{start} && $rand <= $seq->{$taxid}{end})
        {
            return $taxid;
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
        $outputFile, $filehandle, $pair) = @_;

    my $header = "READ_$id|$nameOfSequence|$start-".eval{$start+$readLength}."|$outputFile/$pair";
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
            $readLength,$outputFile,$fastqOutputONE,1);
    ##################################################
    #read two
    ##################################################

    #a. reverse
    $readnt = scalar reverse $readnt;
    #b. complement
    $readnt =~ tr/ATGC/TACG/;
    $newsequence = mutate($readnt, $qualTWO, $readLength);

            writeSequence($id, $newsequence,$asciiQualTWO,$taxaofchoice,$genomelocation,
            $readLength,$outputFile,$fastqOutputTWO, 2);
    $id++;
}

sub readFasta
{
    my ($fastaFile) = @_;
    my ($child, $genus);
    open my $in, '<', $fastaFile;
    while (<$in>)
    {
        next if(m/^\s*$/); #remove blank lines
        if(m/^>\s*?(\S+)/)
        {
            #>genus|1007|leaftaxa|984262
            chomp;
            ($genus, $child) = (split /\|/, $1)[1,3];  #the taxid of the sequence/genome
        }else{
            chomp;
            $seq->{$genus}{sequence} .= $_;
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
