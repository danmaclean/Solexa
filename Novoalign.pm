package Novoalign;

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
#creates the Novoalign record ,

sub new {
	my $self = {};
	my $class = shift;
	my %info = %{$_[0]};#shift;

	$$self{'read_header'} = $info{'-read_header'} || undef;
	$$self{'endedness'} = $info{'-endedness'} || undef;
	$$self{'read_seq'} = $info{'-read_seq'} || undef;
	$$self{'base_quals'} = $info{'-base_quals'} || undef;
	$$self{'status'} = $info{'-status'} || undef;
	$$self{'alignment_score'} = $info{'-align_score'} || undef;
	$$self{'number_alignments'} = $info{'-number_alignments'} || undef;
	$$self{'alignment_quality'} = $info{'-align_qual'} || undef;
	$$self{'aligned_seq_header'} = $info{'-aligned_seq_header'} || undef;
	$$self{'aligned_offset'} = $info{'-aligned_offset'} || undef;
	$$self{'strand'} = $info{'-strand'} || undef;
	$$self{'pair'} = $info{'-pair'} || undef;
	$$self{'pair_offset'} = $info{'-pair_offset'} || undef;
	$$self{'pair_strand'} = $info{'-pair_strand'} || undef;
	$$self{'mismatches'} = $info{'-mismatches'} || undef;
	
	return bless($self, $class);
}

sub read_header{
	my $self = shift;
	return $$self{'read_header'};
}
sub endedness{
	my $self = shift;
	return $$self{'endedness'};
}
sub read_seq{
	my $self = shift;
	return $$self{'read_seq'};
}
sub base_quals{
	my $self = shift;
	return $$self{'base_quals'};	
}
sub status{
	my $self = shift;
	return $$self{'status'};	
}
sub alignment_score{
	my $self = shift;
	return $$self{'alignment_score'};
}
sub number_alignments{
	my $self = shift;
	return $$self{'number_alignments'};
}
sub alignment_quality{
	my $self = shift;
	return $$self{'alignment_quality'};
}
sub aligned_seq_header{
	my $self = shift;
	return $$self{'aligned_seq_header'};
}
sub aligned_offset{
	my $self = shift;
	return $$self{'aligned_offset'};	
}
sub strand{
	my $self = shift;
	return $$self{'strand'};	
}
sub pair{
	my $self = shift;
	return $$self{'pair'};
}
sub pair_offset{
	my $self = shift;
	return $$self{'pair_offset'};
}
sub pair_strand{
	my $self = shift;
	return $$self{'pair_strand'};
}
sub mismatches{
	my $self = shift;
	return %{$$self{'mismatches'}};
	
}