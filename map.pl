#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper::Concise;

our $DIST         = 0.05;
our $SONIC        = 340.0;
our $FS           = 96000.0;
our $TS           = 0.0000104166667;
our $MAX_LAG      = 0.0001470588235;
our $LIMIT_SAMPLE = int(($MAX_LAG / $TS) + 1.5);
our $RANGE        = $LIMIT_SAMPLE * 2;

sub calc_answer;
sub pre_calc;
sub main_process;

&main_process();


sub main_process {

    my @lag_list;
    for (0 .. $RANGE) {
        push @lag_list, ($_ - $LIMIT_SAMPLE) * $TS;
    }
    &pre_calc(\@lag_list);

    print "DONE\n";
    return;
}

sub pre_calc {
    my $lag_list = shift;

    open(FH, ">>/Users/kosuke/Acoust/prog/map_p.csv");

    for my $i(0 .. $RANGE) {
        for my $j(0 .. $RANGE) {
            for my $k(0 .. $RANGE) {
                my $answers = calc_answer($lag_list->[$i],
                                          $lag_list->[$j],
                                          $lag_list->[$k]);
                # TODO:

            }
        }
    }

    close(FH);
}

sub calc_answer {

}
