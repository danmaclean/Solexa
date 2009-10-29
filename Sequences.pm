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
sub split_line_in_two{
	my $line = shift;

	if (length $line % 2 eq 0 ){

		my $line1 = substr($line, 0, ($line / 2);
		my $line2 = substr($line, ($line / 2), $line / 2);

		return ($line1, $line2);
	}
	else {

		return 0;

	}

}


sub remove_adapter_3{ ##removes an adapter by exact match, adapter may have extra nt in read ... so returns chopped seq and  pos of adapter match  
    my $seq = shift;
    my $adapter = shift;
    my $match = 0;
   # warn "in sub, \n\n";
    while ($seq =~ m/$adapter/g){
	$match = length($`);     
    }
   # warn "match = $match\n\n";
    if($match){
	my $start = length($seq) - $match;
	$seq =  substr($seq, 0, (-1 * $start) );
#	warn $seq, "\n";
	return $seq, $start;
    }
    else{
	return 0;

    }

}

sub chop_n_3{ ## chops characters of the 3 
    my $seq = shift;
    my $length = shift; 
    return substr($seq, 0, -$length) 


}


1;
