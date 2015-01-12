#### Breakdown of the script


#### Data structure
```perl
#Example for one sequence
$seq = {
    genus => {
        sequence    => ‘ATCG’,      #sequence for that genus
        length      => 300,         #length of the sequence
        abundance   => 0.75         #relative abundance (as percentage)
        start       => 0            #the start position on the uniform distribution (badly phrased)
        end         => 0.75         #the end position on the uniform distribution
    }
}
```


##### Phred Score calculation

The phred and map functions stores 1.) Phred Quality Score and 2.) base-calling error probabilities 
as key-value pairs in hash `%score`.

```perl
sub phred
{
    my ($score) = @_;
    return 10**($score/(-10));
}

map {$score{$_}=phred($_)} 0..100
```


Function `gaussian_rand` generates a random number from a uniform distribution between 0 and 1
while `random_normal` uses the former to generate a number taken from a normal distribution.

The later is used to generate an insert size before truncation into read pairs.

```perl
sub gaussian_rand 
{
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
```

##### Substitution matrix

Substitution of template nucleotide upon failure of accurate base calling (based on phred score from template fastQ) is random (uniform distribution).

```perl
sub rand_nt
{
    my $rand = rand;
    return 'a' if($rand<$globalnt{'a'});
    return 't' if($rand<$globalnt{'a'}+$globalnt{'t'});
    return 'g' if($rand<$globalnt{'a'}+$globalnt{'t'}+$globalnt{'g'});
    return 'c' if($rand<$globalnt{'a'}+$globalnt{'t'}+$globalnt{'g'}+$globalnt{'c'});
    return 'n';
}
```

The mutation is based on ....
```perl
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
                    $loc++; #shouldnt this be recursive?
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
```
