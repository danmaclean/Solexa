package Parser;
#################################
#Dan MacLeanFeb 14th 2010 dan.maclean@tsl.ac.uk
#module of methods for dealing with raw solexa reads .. 
#################################

####Just adding a pointless comment####


use strict;
use Exporter;
use FileHandle;
use Data::Dumper;
use Solexa::Fastq;
use Solexa::Bowtie;
use Solexa::Novoalign;
use Solexa::Soap;
use Solexa::Maqsnp;
use Solexa::Novoalign;
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
	###very messily go through the args array and set up an easier hash..... 
	for (my $i = 0; $i< scalar(@_); ++$i){
		if ($_[$i] eq '-format'){
			$arg{'-format'} = $_[$i+1];

		}
		if ($_[$i] eq '-file'){
			$arg{'-file'} = $_[$i+1];

		}
	}

	die "unknown format, allowed are bowtie fastq soap maqsnp novoalign \n\n" unless $arg{'-format'} =~ m/bowtie|fastq|soap|maqsnp|novoalign/i;

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

	elsif($$self{'format'} =~ m/maqsnp/){
		my $line = $$self{'handle'}->getline;
		my @tmp = split(/\t/,$line);
		my $maqsnp = new Maqsnp($tmp[0],$tmp[1],$tmp[2],$tmp[3],$tmp[4],$tmp[5],$tmp[6],$tmp[7],$tmp[8],$tmp[9],$tmp[10],$tmp[11]);
		return $maqsnp;
	

	}

	elsif($$self{'format'} =~ m/novoalign/i){
		my $line = '#';
		while ($line =~ m/^#/){
			$line = $$self{'handle'}->getline;
		}
		my @els = split(/\t/, $line);
		my $read_header = $els[0];
		my $endedness = $els[1];
		my $read_seq = $els[2];
		my $base_quals = $els[3];
		my $status = $els[4];
		my ($align_score, $num_aligns, $align_qual, $aligned_seq_header,$aligned_offset, $strand, $pair, $pair_offset, $pair_strand, $mismatches);
		if ($status eq 'U'){
			$align_score = $els[5];
			$align_qual = $els[6];
			$aligned_seq_header = $els[7]; $aligned_seq_header =~ s/^>//;
			$aligned_offset = $els[8];
			$strand = $els[9];
			$pair = $els[10];
			$pair = $aligned_seq_header if $pair eq '.';
			$pair_offset = $els[11];
			$pair_strand = $els[12];
			$mismatches = $els[13];
			my %ms;
			if (defined $mismatches){
				my @t = split(/\s/, $mismatches);
				foreach my $m (@t){
					if ($m =~ m/>/){
						my @a = split(/>/,$m);
						$a[0] =~ m/(\d+)([ATCGN])/;
						my $offset = $1;
						my $orig = $2;
						$ms{'mismatches'}{$offset}{$orig} = $a[1];
					}
					elsif ($m =~ m/\+/){
						my @a = split(/\+/,$m);
						$ms{'inserts'}{$a[0]} = $a[1];						
					}
					elsif ($m=~ m/-/){
						my @a = split(/-/,$m);
						$ms{'deletes'}{$a[0]} = $a[1];						
					}
				}
			}
			my $novoalign = new Novoalign({'-read_header' => $read_header, '-endedness' => $endedness, '-read_seq' => $read_seq,
										'-base_quals' => $base_quals, '-status' => $status,  '-align_score' => $align_score, '-align_qual' => $align_qual,
										'-aligned_seq_header' => $aligned_seq_header, '-aligned_offset' => $aligned_offset, '-strand' => $strand, '-pair' => $pair,
										 '-pair_offset' => $pair_offset, '-pair_strand' => $pair_strand, '-mismatches' => \%ms}
			);
			return $novoalign;
		}
		elsif($status eq 'R'){
			$num_aligns = $els[5];
			$align_qual = $els[6];
			$aligned_seq_header = $els[7]; $aligned_seq_header =~ s/^>//;
			$aligned_offset = $els[8];
			$strand = $els[9];
			$pair = $els[10];
			$pair = $aligned_seq_header if $pair eq '.';
			$pair_offset = $els[11];
			$pair_strand = $els[12];
			$mismatches = $els[13];
			my %ms;
			if (defined $mismatches){
				my @t = split(/\s/, $mismatches);
				foreach my $m (@t){
					if ($m=~m/>/){
						my @a = split(/>/,$m);
						$a[0] =~ m/(\d+)([ATCGN])/;
						my $offset = $1;
						my $orig = $2;
						$ms{'mismatches'}{$offset}{$orig} = $a[1];
					}
					elsif ($m=~ m/\+/){
						my @a = split(/\+/,$m);
						$ms{'inserts'}{$a[0]} = $a[1];						
					}
					elsif ($m=~ m/-/){
						my @a = split(/-/,$m);
						$ms{'deletes'}{$a[0]} = $a[1];						
					}
				}
			}
			my $novoalign = new Novoalign({'-read_header' => $read_header, '-endedness' => $endedness, '-read_seq' => $read_seq,
										'-base_quals' => $base_quals, '-status' => $status,  '-number_alignments' => $num_aligns, '-align_qual' => $align_qual,
										'-aligned_seq_header' => $aligned_seq_header, '-aligned_offset' => $aligned_offset, '-strand' => $strand, '-pair' => $pair,
										 '-pair_offset' => $pair_offset, '-pair_strand' => $pair_strand, '-mismatches' => \%ms}
			);
			return $novoalign;	
		}
		else{
			my $novoalign = new Novoalign({'-read_header' => $read_header, '-endedness' => $endedness, '-read_seq' => $read_seq,
										'-base_quals' => $base_quals, '-status' => $status});
			return $novoalign;	
		}
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
Solexa::Novoalign

=head1 METHODS

=over

=item new(-file=> "somefile", -format=>'x')

Creates a new Parser object of format 'x', where x can be ONE of 'fastq' or 'bowtie' or 'soap' or 'novoalign

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
Solexa::Novoalign
Solexa::Soap
