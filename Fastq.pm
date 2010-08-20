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
sub quals_to_sanger{
	my $self = shift;
	my $new_qual = _sol2sanger($self->quals);
	$$self{'quals'} = $new_qual;
}
sub split{
###splits the read in two and creates two new fastq objects	
	my $self = shift;
	return 0 unless length($self->seq) % 2 == 0;  
	my $string_1 = substr($self->seq, 0, length($self->seq) / 2) ;
	my $string_2 = substr($self->seq, length($self->seq) / 2, length($self->seq) / 2 );
	my $qual_1 = substr($self->quals, 0, length($self->quals) / 2);
	my $qual_2 = substr($self->quals, length($self->quals) / 2, length($self->quals) / 2 );
	#my $id_1 = $self->id.'_1';
	#my $id_2 = $self->id.'_2';
	#my $qid_1 = $self->qual_id.'_1';
	#my $qid_2 = $self->qual_id.'_2';
	my $id_1 = $self->id.'/1';
	my $id_2 = $self->id.'/2';
	my $qid_1 = $self->qual_id.'/1';
	my $qid_2 = $self->qual_id.'/2';
	my $fastq_1 = Fastq->new($id_1,$string_1,$qid_1,$qual_1);
	my $fastq_2 = Fastq->new($id_2,$string_2,$qid_2,$qual_2);
	return ($fastq_1, $fastq_2);		
}
sub write{

	my $self=shift;
	my $p = $self->id . "\n" .$self->seq . "\n" . $self->qual_id , "\n" , $self->quals, "\n";

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
sub _sol2sanger{
	my $quals = shift;
		my @quals = split( '', $quals );
		my $qual = '';
		foreach my $q (@quals){
			my $s = chr(ord($q) - 31);
			$qual = $qual . $s; 
			
		}
		return $qual;

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

=item quals_to_sanger

Alters the quality score of the line to the sanger-fastq standard ***only use if your read is already illumina-fastq standard***
That is from post pipeline version 1.3.

=item split()

Splits the read and quals in half and returns two new fastq objects with _1 and _2 appended to the ids

my ($split1,$split2) = $fastq->split;

=back

=head1 SEE ALSO
Solexa::Parser
