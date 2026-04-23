#!/usr/bin/perl

use warnings;
use strict;
use POSIX 'locale_h';

# Local
setlocale(LC_ALL, 'en_US.UTF-8');

my $filename = 'Mouse_ids';

open(my $FH, '<', $filename) or die $!;

while (<$FH>) {
    chomp;  # Delete line breaks 
    my $accession_id = $_;

    # Use datasets to obtain IDs
    my $command = qq(../datasets summary gene accession "$accession_id" --report product | sed -n 's/.*"symbol":"\\([^"]*\\).*/\\1/p');

    # Get the output of command
    my $resultado = `$command`;

    # format "IDQUERY | IDresult"
    print("$accession_id | $resultado\n");
}

close($FH);
