#!/usr/bin/perl

use warnings;
use strict;

my $filename = 'Mouse_hitdata.txt_filt2';
my $filename2 = 'CDD_channels_new_ids_list.txt';

open(FH, '<', $filename) or die $!;

while(<FH>){
   if ($_ =~ /^Q#/){  # Aquí se agrega el 'or' para 'NP'
        print $_;
   }else{
        system("grep -w '$_'  ../../CDD_channels_new_ids_list.txt");
   }
}

close(FH);
