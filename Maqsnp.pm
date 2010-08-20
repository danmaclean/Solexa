package Maqsnp;

use strict;
use Exporter;
use FileHandle;
use Data::Dumper;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()],
				ALL =>[qw()]);
#creates the Fastq sequence object,
sub new {
	my $self = {};
	my $class = shift;
	my ($chr, $pos, $ref_base, $consensus, $phred_consensus_qual, $read_depth, $average_hits, $highest_map_qual, $min_consensus_qual, $second_best, $log_likelihood, $third_best) = @_;
	
	$$self{'chromosome'} = $chr;
	$$self{'pos'} = $pos;
	$$self{'ref_base'} = $ref_base;
	$$self{'consensus'} = $consensus;
	$$self{'phred_consensus_qual'} = $phred_consensus_qual;
	$$self{'read_depth'} = $read_depth;
	$$self{'average_hits'} = $average_hits;
	$$self{'highest_map'} = $highest_map_qual;
	$$self{'min_consensus_qual'} = $min_consensus_qual;
	$$self{'second_best'} = $second_best;
	$$self{'log_likelihood'} = $log_likelihood;
	$$self{'third_best'} = $third_best;
	return bless($self, $class);
}


#getters, returns  information to you
sub chromosome {
	my $self = shift;
	$$self{'chromosome'} || 0;
}
sub position {
	my $self = shift;
	$$self{'pos'} || 0;
}
sub ref_base{
	my $self = shift;
	$$self{'ref_base'} || 0;
}
sub consensus{
	my $self = shift;
	$$self{'consensus'} || 0;
}
sub phred_consensus_qual {
	my $self = shift;
	$$self{'phred_consensus_qual'} || 0;
}
sub read_depth {
	my $self = shift;
	$$self{'read_depth'} || 0;
}
sub average_hits{
	my $self = shift;
	$$self{'average_hits'} || 0;
}
sub highest_map{
	my $self = shift;
	$$self{'highest_map'} || 0;
}
sub min_consensus_qual {
	my $self = shift;
	$$self{'min_consensus_qual'} || 0;
}
sub second_best {
	my $self = shift;
	$$self{'second_best'} || 0;
}
sub log_likelihood{
	my $self = shift;
	$$self{'log_likelihood'} || 0;
}
sub third_best{
	my $self = shift;
	$$self{'third_best'} || 0;
}
1;


=head1 NAME

Maqsnp - module that creates object representing maq snps from cns2snp

=head1 AUTHOR

Dan MacLean (dan.maclean@tsl.ac.uk)

=head1 SYNOPSIS

	use Solexa::Maqsnp
	my $snp = new Fastq($chr, $pos, $ref_base, $consensus, $phred_consensus_qual, $read_depth, $average_hits, $highest_map_qual, $min_consensus_qual, $second_best, $log_likelihood, $third_best);
	print $snp->consensus #prints out the consensus base
	print $snp->second_best; #prints the second best base call

=head1 DESCRIPTION

Fastq files are the sequence file variety created by many HTGS machines such as the Illumina GA2. This module creates an object
that allows you to access the individual parts of the sequence entry very easily. To create the object you must provide it with the 
sequence, qualities and ids. This module is really designed to be used with the Solexa::Parser module that will create the objects
automatically. Further methods will be created soon.

=head1 METHODS

=over

=item new(id,sequence,quality_id,quality_scores)

Creates a new Fastq object

	$id = '@LANE_ID_CO_ORD_NUMBER'
	$seq = 'ATGCATGCNNNATGCATGC';
	$qual_id = '+LANE_ID_CO_ORD_NUMBER';
	$qual = '!@£$%^&*(()_++_)*';
	
	$fastq = new Fastq($id, $seq, $qual_id, $qual);

=item id()

Returns the id of the Fastq object

	print $fastq->id;   #prints '@LANE_ID_CO_ORD_NUMBER'
	
=item seq()

Returns the sequence of the Fastq object

	print $fastq->seq;   #prints 'ATGCATGCNNNATGCATGC'

=item qual_id()

Returns the quality line id of the Fastq object

	print $fastq->qual_id;   #prints '+LANE_ID_CO_ORD_NUMBER'

=item quals()

Returns the quality score line of the Fastq object

	print $fastq->quals;   #prints '!@£$%^&*(()_++_)*'

=item quals_to_sanger

Alters the quality score of the line to the sanger-fastq standard ***only use if your read is already illumina-fastq standard***
That is from post pipeline version 1.3.

=item split()

Splits the read and quals in half and returns two new fastq objects with _1 and _2 appended to the ids

my ($split1,$split2) = $fastq->split;

=back

=head1 SEE ALSO
Solexa::Parser
