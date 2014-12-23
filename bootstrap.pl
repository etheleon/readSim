#!/usr/bin/env perl

use Modern::Perl;

die "$0 <abundance profile>\n" unless $#ARGV == 0;

#This script will control the execution of the whole pipeline
my @path = (split '/', $0);
pop @path;
my $baseDIR = join "/", @path;

system $baseDIR."/"."readSim.0100.chooseGenomes.r";
