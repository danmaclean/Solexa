package Fastq;

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
	my ($id,$seq,$qual_id,$quals) = @_;
	chomp($id,$seq,$qual_id,$quals);
	$$self{'id'} = $id;
	$$self{'seq'} = $seq;
	$$self{'qual_id'} = $qual_id;
	$$self{'quals'} = $quals;
	return bless($self, $class);
}

#getters, returns sequence information to you
sub seq {
	my $self = shift;
	$$self{'seq'} || 0;
}
sub id {
	my $self = shift;
	$$self{'id'} || 0;
}
sub qual_id{
	my $self = shift;
	$$self{'qual_id'} || 0;
}
sub quals{
	my $self = shift;
	$$self{'quals'} || 0;
}


1;


=head1 NAME

Fastq - module that creates object representing fastq sequence files

=head1 AUTHOR

Dan MacLean (dan.maclean@tsl.ac.uk)

=head1 SYNOPSIS

	use Solexa::Fastq;
	my $fastq = new Fastq($id, $sequence, $qual_id, $quals);
	print $fastq->seq; #prints out the sequence line
	print $fastq->qual_id; #prints the quality line id

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

=back

=head1 SEE ALSO
Solexa::Parser