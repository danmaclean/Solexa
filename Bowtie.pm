package Bowtie;
#################################
#
#Class for Bowtie Alignments .. 
#################################
use strict;
use Exporter;
use Digest::MD5 qw(md5 md5_hex md5_base64);
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
	chomp $line;
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

    if (defined $tmp[7] and $tmp[7] =~ m/\w+/){
		my %subs;

		my @tmp2 = split(/,/,$tmp[7]);
		foreach my $t (@tmp2){
	    	my @tmp3 = split(/:/,$t);
	    	my @tmp4 = split(/>/,$tmp3[1]);
			if ($tmp[1] eq '+'){
				my $offset = $tmp3[0] + $tmp[3];
	    			$subs{$offset}{$tmp4[0]}= $tmp4[1];
			}
        	else{
				my $offset = ($tmp[3] + length($tmp[4])) - ($tmp3[0] + 1);
				$subs{$offset}{$tmp4[0]}= $tmp4[1];
			}
		}	
	   	$$self{'_mismatches'} = \%subs;
	}
    else{
		$$self{'_mismatches'} = 'None';
    }
    
    bless $self, $class;

}


##accessors
sub readname{ 
    my $self = shift;
    return $$self{'_read_name'};

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
sub alignstop {
	my $self = shift;
	my $length = $self->alignstart + length( $self->readseq);
	return $length;
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
	return %{$$self{'_mismatches'}};
    }
}
sub gee_fu_gff{
	my $self = shift;
	my $id = md5_hex(sort %$self);
	my $seq = quote_string($self->readseq);
	my $qual = quote_string($self->readqual);
	my $group = "ID=$id;Note=<sequence>$seq</sequence><quality>$qual</quality>"; 
	my $gff = $self->targetname . "\t" . 'bowtie' . "\t" . 'read' . "\t" . $self->alignstart . "\t" . $self->alignstop . "\t" . '.' . "\t" . $self->strand . "\t" . '.' . "\t" . $group;
	return $gff;
}
sub quote_string{
	my $string = shift;
	$string =~ s/%/%25/g;
	$string =~ s/;/%3B/g;
	$string =~ s/=/%3D/g;
	$string =~ s/&/%26/g;
	$string =~ s/,/%2C/g;
	return $string;
}

1;


=head1 NAME

Bowtie - module that creates object representing Bowtie alignment files

=head1 AUTHOR

Dan MacLean (dan.maclean@tsl.ac.uk)

=head1 SYNOPSIS

	use Solexa::Bowtie;
	my $bwt = new bowtie($bowtie_line);
	print $bwt->readname; #prints out the reads name in this alignment
	print $bwt->readseq; #prints the reads sequence

=head1 DESCRIPTION

Bowtie is an NGS alignment program. This module creates an object that allows you to access the individual parts of the alignment
entry very easily. To create the object you must provide it with a single line from the Bowtie output file. This module is really 
designed to be used with the Solexa::Parser module that will create the objects automatically. Further methods will be created 
soon.

=head1 METHODS

=over

=item new(line)

Creates a new Bowtie object

	$line = 'line-from-a-bowtie-file';
	$bwt = new Bowtie($line);

=item readname()

Returns the readname of the alignment

	print $readname->readname;   #prints the id of the read in the alignment'

=item strand()

Returns the strand of the alignment

	print $bwt->strand;

=item targetname()

Returns the name of the target

	print $bwt->targetname;

=item alignstart()

Returns the start of the alignment on the reference

	print $bwt->alignstart();

=item alignstop()

Returns the stop position of the alignment on the reference

	print $bwt->alignstop();

=item readseq()

Returns the sequence of the read

	print $bwt->readseq()

=item readqual()

Returns the read qualities;

	print $bwt->readqual;

=item other_alignments()

Returns the number of other alignments this read is involved in

	print $bwt->other_alignments;

=item mismatches()

Returns a hash describing any mismatches that occured between the read and the reference in this alignment, if no mismatches occur
returns the string 'None'

	my %mismatches = $bwt->mismatches;

The hash returned has as keys the positions in which the read differs from the reference. The nucleotide number returned is that in the reference 
DNA, not the position in the read. There are some oddities in the Bowtie output, such as when a read matches the negative strand the mismatch read
must be read from the opposite end of the read. The hash looks like this...

	%mismatches{
		'137' => {'A' => 'G'},
		'124' => {'C => 'T'}
		}

=item gee_fu_gff()

Returns a gee_fu compatible gff string of the current object

	my $gff = $bwt->gee_fu_gff;

=back

=head1 SEE ALSO
Solexa::Parser
