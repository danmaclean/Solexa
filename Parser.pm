package Parser;
#################################
#Dan MacLean Jun 19th 2009   dan.maclean@tsl.ac.uk
#module of methods for dealing with raw solexa reads .. 
#################################
use strict;
use Exporter;
use FileHandle;
use Data::Dumper;
use Solexa::Fastq;
use Solexa::Bowtie;
use Solexa::Soap;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()],
				ALL =>[qw()]);


sub new {
	my $class = shift;
	my %arg;
	print Dumper @_;
	###very messily go through the args array and set up an easier hash..... 
	for (my $i = 0; $i< scalar(@_); ++$i){
		if ($_[$i] eq '-format'){
			$arg{'-format'} = $_[$i+1];

		}
		if ($_[$i] eq '-file'){
			$arg{'-file'} = $_[$i+1];

		}
	}
	die "unknown format, allowed are bowtie fastq soap\n\n" unless $arg{'-format'} =~ m/bowtie|fastq|soap/i;
	die "no file provided\n\n" unless $arg{'-file'};
	my $self =  {};
	$$self{'format'} = $arg{'-format'};
	$$self{'file'} = $arg{'-file'};
	my $fh = new FileHandle;
	$$self{'handle'} = $fh;
	$fh->open($$self{'file'}, 'r') || die "couldn't open file $$self{'file'}\n\n";
	return bless($self, $class);
}


sub next{
	my $self = shift;
	if ($$self{'handle'}->eof){
		$$self{'handle'}->close;
		return 0;
	}
	if ($$self{'format'} =~ m/fastq/i){
		my $id = $$self{'handle'}->getline;
		if ($id =~ m/^@/){
			my $seq = $$self{'handle'}->getline;
			my $qual_id = $$self{'handle'}->getline;
			my $quals = $$self{'handle'}->getline;
			chomp($id,$seq,$qual_id,$quals);
			my $fastq = new Fastq($id, $seq, $qual_id, $quals);
			return $fastq;
		}
		else {
			die "Some error has occurred, was expecting an Illumina id line but got $id";
		}
	}
	elsif ($$self{'format'} =~ m/bowtie/){
		my $line = $$self{'handle'}->getline;
		my $bowt = new Bowtie($line);
		return $bowt;

	}
	elsif($$self{'format'} =~ m/soap/){
		my $line = $$self{'handle'}->getline;
		my $soap = new Soap($line);
		return $soap;


	}
	else{

		die "unknown format provided\n\n";

	}

}
1;	  



=head1 NAME

Solexa::Parser - module to allow easy reading of different file formats in use at TSL.

=head1 AUTHOR

Dan MacLean (dan.maclean@tsl.ac.uk)

=head1 SYNOPSIS

	use Solexa::Parser;
	my $parser = new Parser(-format=>'fastq', -file=>'my_solexa_data.fq');
	while(my $fastq = $parser->next){

		print $fastq->seq; #prints out the sequence line
		print $fastq->qual_id; #prints the quality line id
	}

=head1 DESCRIPTION

This module allows you to create a datastream of the individual records in different file formats, so that you dont need to know
how to parse them. You tell the parser in the new command what sort of file you want to use and the file name and it passes one record
at a time back to you. What you can do with that record depends on exactly what sort of file you had in the first place. You can
use Solexa::Fastq methods with a Fastq file, Solexa::Bowtie methods, with a bowtie file and so on.

=head1 REQUIRES

Solexa::Parser relies on the following modules being present
Solexa::Bowtie
Solexa::Fastq
Solexa::Soap

=head1 METHODS

=over

=item new(-file=> "somefile", -format=>'x')

Creates a new Parser object of format 'x', where x can be ONE of 'fastq' or 'bowtie' or 'soap'

	$file = '/home/myhome/fastq/my_solexa_data.fq';
	$parser = new Parser(-file=>$file, -format=>'fastq'); 

=item next()

Returns the next record from the file

	$fastq = $parser->next;

Use in a while loop to go through every record in the input file

	while ($fastq = $parser->next){

		#do something in here...

	}

=back

=head1 SEE ALSO

Solexa::Fastq
Solexa::Bowtie
Solexa::Soap
