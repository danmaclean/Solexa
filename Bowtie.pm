package Bowtie;
#################################
#
#Class for Bowtie Alignments .. 
#################################
use strict;
use Exporter;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()],
				ALL =>[qw()]);


sub new {
    my $class = shift;
    my $line = shift || die "No alignment\n\n";
    my $self = {};

    my @tmp = split(/\t/,$line);
    
    $$self{'_read_name'} = $tmp[0]; #|| die "missing data";
    $$self{'_strand'} = $tmp[1];# || die "missing data";
    $$self{'_target_name'} = $tmp[2];# || die "missing data";
    $$self{'_alignment_start'} = $tmp[3];# || die "missing data";
    $$self{'_read_seq'} = $tmp[4]; #|| die "missing data";
    $$self{'_read_qual'} = $tmp[5]; #|| die "missing data";#
#	warn $tmp[3], "\n";
    $$self{'_number_of_other_alignments'} = $tmp[6];  #|| die "missing data";

    if ($tmp[7] =~ m/\w+/){
	my %subs;

	my @tmp2 = split(/,/,$tmp[7]);
	foreach my $t (@tmp2){
	    my @tmp3 = split(/:/,$t);
	    my @tmp4 = split(/>/,$tmp3[1]);
	    $subs{$tmp3[0]}{$tmp4[0]}= $tmp4[1];
	    $$self{'_mismatches'} = \%subs;
        }
    }
    else{
	$$self{'_mismatches'} = 'None';
    }
    
    bless $self, $class;

}


##accessors
sub readname{ 
    my $self = shift;
    return $$self{'_readname'};

}

sub strand {
    my $self = shift;
    return $$self{'_strand'};

}
sub targetname {
    my $self = shift;
    return $$self{'_target_name'};
}
sub alignstart{
    my $self = shift;
    return $$self{'_alignment_start'};
}

sub readseq{
    my $self = shift;
    return $$self{'_read_seq'};
}

sub readqual{
    my $self = shift;
    return $$self{'_read_qual'};
}
sub other_alignments{
    my $self = shift;
    return $$self{'_number_of_other_alignments'};
}
sub mismatches{
    my $self = shift;
    if ($$self{'_mismatches'} eq 'None'){

	return '0';

    }
    else{
	return \%{$$self{'_mismatches'}};
    }
}


1;
