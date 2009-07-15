package Sequences;
#################################
#
#module of methods for dealing with raw solexa reads .. 
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

sub is_solexa_fastq{
## takes a quality score line and scans it to see if the line contains chars with ascii values of greater than 73.. which indicates solexa
	
	
	my @quals = split(/|/,$_[0]);
	foreach my $q (@quals){

		if (ord($q) > 73){
			
			return 1;

		}

	}

	
	return 	0;


}

sub sol2sanger{
	die "error in sol2sanger " unless $_[1];
	my $vers = shift; ## vers is the version of the solexa pipeline

	if ($vers =~ m/1.4/){
		my $quals = shift;
			# Added to eliminate carriage return conversion
		chomp $quals;
		my @quals = split( '', $quals );
		my $qual = '';
		foreach my $q (@quals){
			my $s = chr(ord($q) - 31);
			$qual = $qual . $s; 
			
		}
		return $qual;
	}
	else {

		die "No known version\n\n";

	}
}
sub get_read_lengths{
	my %hash;
	warn "@{$_[0]}";
	foreach my $file (@{$_[0]}){
	## returns reference to hash of lengths for every read
	open FILE, "<$file" || die "can't open $file\n\n";
		while (my $line = <FILE>){
			chomp $line;
			if ($line =~ m/^@/){
				my $seq = <FILE>;
				chomp $seq;
				my $l = length($seq);
				$hash{$l}{$file}{$line} =1;		

			}
		}
	close FILE;
	}
	return \%hash;
}
1;
