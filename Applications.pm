package Applications;
#################################
#
#module of methods for sending solexa app jobs to the cluster .. 
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

my $maq_dir = '/home/studhold/maq-0.6.8_x86_64-linux';

sub make_bfq{ ## takes an input fastq and makes a bfq
	
	my $input = shift;
	my $output = $input;
	$output =~ s/\.fastq$/\.bfq/;
	system("bsub \"$maq_dir/maq fastq2bfq $input $output\" "); 
	return 0;

}

sub run_map {

	my $chr = shift;
	my $reads = shift;
	my $output = $reads;
	$output =~ s/\.bfq$/\.map/;
	
	system("bsub -q bigmem  \"$maq_dir/maq map $output $chr $reads\" "); 
	return 0;

}

1;
