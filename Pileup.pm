package Pileup;
#################################
#Dan MacLean Jun 19th 2009   dan.maclean@tsl.ac.uk
#module of methods for dealing with raw solexa reads .. 
#################################
use strict;
use Exporter;
use Storable qw(nstore retrieve);
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()],
				ALL =>[qw()]);


sub new {
    my $class = shift;
    my $file = shift || die;
    my $self = {};
    #warn $file, "\n";
    $$self{'_file'} = $file;

    open FILE, "<$file" || die "Can't get open $file\n\n";

    while (my $line = <FILE>){

	#chomp $line;
	my @tmp = split(/\s+/, $line);
	next if exists $$self{'_index'}{$tmp[0]};
    my $start = tell(FILE);
	$line = encode('UTF-8', $line);
	$$self{'_index'}{$tmp[0]} = $start - length($line); ##works through the pileup file and records the distance in bytes into the file that each contig starts
#	warn "$tmp[0]\n";		     
 
	    
    }
    close FILE;
   # warn "Index built for $file\n\n";
   # return $self;
    bless $self, $class;



}

sub _get_pile{
    my $self = shift;
    my $contig = shift;

   # warn "get pile called, getting pile for $contig\t";

    my %pileup;


    my $seekto = $$self{'_index'}{$contig};
    open FILE, "<$$self{'_file'}";
    seek(FILE, $seekto, 0);
    while (my $line = <FILE>){

	chomp $line;
	my @tmp = split(/\s+/, $line);
	last unless $contig eq  $tmp[0];
	$pileup{$tmp[0]}{$tmp[1]}{'_ref_base'} = $tmp[2];
	$pileup{$tmp[0]}{$tmp[1]}{'_coverage_depth'} = $tmp[3];

	

	my @bases = qw(A C G T N );
	foreach my $b (@bases){

	   $pileup{$tmp[0]}{$tmp[1]}{$b} = 0;

	}
	
	my %comp = ('A' => 'T', 'G' => 'C', 'T' => 'A', 'C' => 'G' );

	my @al = split(/|/, $tmp[4]);
	shift @al;	       

		       #### NB strand information is lost here...! 

	foreach my $b (@al){
	
	    if($b =~ m/\./){

		$pileup{$tmp[0]}{$tmp[1]}{$tmp[2]}++;

	    }

	    elsif ($b =~ m/,/){

		
		$pileup{$tmp[0]}{$tmp[1]}{$tmp[2]}++;

	    }
	    elsif ($b =~ m/n/i){

		$pileup{$tmp[0]}{$tmp[1]}{'N'}++;
	    }
	    elsif ($b eq 'A' or $b eq 'C' or $b eq 'G' or $b eq 'T' ){


		$pileup{$tmp[0]}{$tmp[1]}{$b}++;

	    }
	    else {
		my $u = uc($b);
		#my $n = $comp{$u};
		$pileup{$tmp[0]}{$tmp[1]}{$u}++;
	    }

	}


    }
    close FILE;
   # my @list = keys %{$$self{'_pileup'}};
   # warn "@list", "\n";
    $$self{'_pileup'}= \%pileup;
   # my @list2 = keys (%{$$self{'_pileup'}});
   # warn "@list2" , "\n";
    return $self;

}

sub name_file{
    my $self = shift;
    my $file =  $$self{'_file'};
    return $file;

}
sub pickle {
	
	my $self = shift;
	my $filename = shift;
	nstore($self, $filename); 
	return;

}

sub retrieve_pickle{
	my $file = shift;
	my $pileup = Storable::retrieve("$file");
	return $pileup;
}

sub get_coverage{  ## gets coverage for a single position on a particular contig
    my $self = shift;
    my $cont = shift;
    my $pos = shift;


    
    if (!defined $cont or !defined $pos){

	warn "Undefined contig or position";
	return 0;

    }
    if (!exists $$self{'_index'}{$cont}){

	warn "contig not in set ..!";
	return 0;
    }
    $self->_get_pile($cont) unless exists $$self{'_pileup'}{$cont};
    return $$self{'_pileup'}{$cont}{$pos}{'_coverage_depth'};

}

sub get_ref_base{

    my $self = shift;
    my $cont = shift;
    my $pos = shift;

    $self->_get_pile($cont) unless exists $$self{'_pileup'}{$cont};
    if (!defined $cont or !defined $pos){

	warn "Undefined contig or position";
	return 0;

    }
    if (!exists $$self{'_index'}{$cont}){

	warn "contig not in set ..!";
	return 0;
    }


    return $$self{'_pileup'}{$cont}{$pos}{'_ref_base'};

}

sub get_contigs{ ## returns list of contigs
    my $self = shift;
    #my $size = keys %{$$self;
    return keys %{$$self{'_index'}};
}

sub contig_count{ ## returns number of contigs
    my $self = shift;
    my $size = keys %{$$self{'_index'}};
    return $size;


}

sub contig_length{ ## returns a contig length
    my $self = shift;
    my $cont = shift;
    if (!defined $cont){

	warn "Undefined contig";

	return 0;

    }
    if (!exists $$self{'_index'}{$cont}){

	warn "contig not in set ..!";
	return 0;
    } 

    
    $self->_get_pile($cont) unless exists $$self{'_pileup'}{$cont};

    my @unsorted = keys %{$$self{'_pileup'}{$cont}};
    my  @sorted = sort {$b <=> $a } @unsorted;
    return shift @sorted;


}

sub nt_coverage{ ## returns hash of number of each nt called at a position
    my $self = shift;
    my $cont = shift;
    my $pos = shift;



    if (!defined $pos or !defined $cont){

	warn "Undefined position";
	return;

    }

    if (!exists $$self{'_index'}{$cont}){

	warn "contig not in set ..!";
	return 0;
    }

    $self->_get_pile($cont) unless exists $$self{'_pileup'}{$cont};

    my %t =  %{$$self{'_pileup'}{$cont}{$pos}};
    #delete $t{'_ref_base'};
    #delete $t{'_coverage_depth'};
    return %t;

}

sub average_depth{ ## returns average depth of coverage over all contigs
    my $self = shift;
    my $num = 0;

    if (exists $$self{'_average_depth'}){

	return $$self{'_average_depth'};

    }
    else{
	my $total = $self->total_nucleotides();
	my @contigs = $self->get_contigs();
	foreach my $c (@contigs){
	    $num += $self->contig_average_depth($c);
	
	}
	my $average_depth =  $num / scalar(@contigs);
	$$self{'_average_depth'} = $average_depth;
	return $average_depth;
    }

}

sub contig_average_depth{ ## returns single contig average depth
    my $self = shift;
    my $contig = shift;
  
    if (!defined $contig){

	warn "Undefined contig";
	return 0;

    }
    

    return $self->nucleotides_aligned_to_contig($contig) / $self->contig_length($contig);

}

sub total_nucleotides{ ## returns number of nucleotides in all contigs ie sum of lengths

    my $self = shift;
    my $length = 0;
    
    if (exists $$self{'_total_nucleotides'}){

	return $$self{'_total_nucleotides'};

    }
    else {
	my @contigs = $self->get_contigs;
	foreach my $c (@contigs) {

	    $length += $self->contig_length($c);

	}
	$$self{'_total_nucleotides'} = $length;
	return $length;
    }

}

sub nucleotides_aligned_to_contig{ ## returns the sum of aligned nucleotides at each position in a contig

    my $self = shift;
    my $contig = shift;
    my $num = 0;
    if (!defined $contig){

	warn "Undefined contig";
	return;

    }
    for (my $i = 1; $i <= $self->contig_length($contig); $i++){

	$num += $self->get_coverage($contig, $i);

    }
    return $num;

}

sub nucleotides_aligned{ ## returns the sum of nucleotides at each position over all contigs

    my $self = shift;
    my $num = 0;

    if (exists $$self{'_nucleotides_aligned'}){

	return $$self{'_nucleotides_aligned'};

    }
    else{

	foreach my $contig ($self->get_contigs){

	    $num += $self->nucleotides_aligned_to_contig($contig);

	}
	$$self{'_nucleotides_aligned'} = $num;
	return $num;
    }
}

sub _make_pileup_array {  ### makes a 2D array of the maq pileup with the reference in the first column

    my $self = shift;
    my $contig = shift;

    my $seekto = $$self{'_index'}{$contig};
    open FILE, "<$$self{'_file'}";
    seek(FILE, $seekto, 0);
    my @pile;
    my $max_cov = 0;
    while (my $line = <FILE>){

	chomp $line;
	my @tmp = split(/\s+/, $line);
	last unless $contig eq  $tmp[0];
	$max_cov = $tmp[3] if $tmp[3] > $max_cov;
	my $ref = $tmp[2];
 	my @al = split(/|/, $tmp[4]);
	#shift @al;	       
	
	$pile[$tmp[1] - 1 ][0] = $ref; 
	for (my $i = 0; $i <= $max_cov; ++$i){
	    my $read_base = $al[$i];
	    next unless defined $read_base;
	    if ($read_base eq '.'){
		$pile[$tmp[1]-1][$i] = $ref;
	    }
	    elsif ($read_base eq ','){

		$pile[$tmp[1]-1][$i] = lc($ref);

	    }
	    elsif($read_base =~ m/[ATCG]/i) {
		$pile[$tmp[1]-1][$i] = $read_base;
	    }
	}

    }
    close FILE;
    ###now fill out all the blank spaces in the array ...

    for (my $i = 0; $i < $self->contig_length($contig); $i++){

	for (my $j = 0; $j < $max_cov; $j++){

	    $pile[$i][$j] = '.' unless defined $pile[$i][$j];

	}
	

    }
    return \@pile;


}
sub _rotate_array { ##makes a rotation of a 2D array .. 

    my $self = shift;
    my $array = shift;##ref to array to rotate
    my @rotated_array;
    #foreach my 	

}

sub html_pile { ## converts  a 2d array in a pretty html format

    my $self = shift;
    my $contig = shift;
    my $orient = shift;
    $orient = 'up' unless defined $orient;
    my $pile = $self-> _make_pileup_array($contig);
    #my $rotated = $self->_rotate_array(\$pile);

    my $length = scalar(@{$$pile[0]});# get the length of the inner array


    my @html;

    if ($orient eq 'up'){
	for (my $i = $length; $i>= 0; $i--){  
	    for (my $j = 0; $j < scalar(@$pile); $j++){
	    
		my $h = '';
		my $el = $$pile[$j][$i];
		my $ref = $$pile[$j][0];
		my $p = $j + 1;
		if ($i == 0){

		    $h = "<span title=\"$p\" class='ref'>" . $el . '</span>';

		}
		else {
		    if (uc($el) =~ m/^$ref$/){
			if ($el =~ m/A/i){
			    $h =  "<span title=\"$p\" class='A'>" . $el . '</span>'; 

			}
			elsif ($el =~ m/T/i){
			    $h =  "<span title=\"$p\" class='T'>" . $el . '</span>';

			}
			elsif ($el =~ m/G/i){
			    $h = "<span title=\"$p\" class='G'>" . $el . '</span>';

			}
			elsif ($el =~ m/C/i){
			    $h = "<span title=\"$p\" class='C'>" . $el . '</span>';

			}
			else{
			    $h = '<span class=\'X\'>' . '&nbsp;'  . '</span>';

			}
		    }
		    else{
			if ($el =~ m/A/i){
			    $h =  "<span title=\"$p\" class='different_A'>" . $el . '</span>'; 

			}
			elsif ($el =~ m/T/i){
			    $h =  "<span title=\"$p\" class='different_T'>" . $el . '</span>';

			}
			elsif ($el =~ m/G/i){
			    $h = "<span title=\"$p\" class='different_G'>" . $el . '</span>';

			}
			elsif ($el =~ m/C/i){
			    $h = "<span title=\"$p\" class='different_C'>" . $el . '</span>';

			}
			else{
			    $h = "<span title=\"$p\" class=\'X\'>" . '&nbsp;'  . '</span>';

			}
			

		    }
	       
		}
		push @html, $h if defined $h;



		}
	    push @html, '<br>',"\n";

	}
    }
    elsif ($orient eq 'down'){
	for (my $i = 0; $i <= $length; $i++){
	    for (my $j = 0; $j < scalar(@$pile); $j++){
	    
		my $h = '';
		my $el = $$pile[$j][$i];
	    	my $ref = $$pile[$j][0];
		my $p = $j + 1;
		if ($i == 0){

		    $h = "<span title=\"$j\" class='ref'>" . $el . '</span>';

		}
		else{

		    if (uc($el) =~ m/$ref/){
			if ($el =~ m/A/i){
			    $h =  "<span title=\"$p\" class='A'>" . $el . '</span>'; 

			}
			elsif ($el =~ m/T/i){
			    $h =  "<span title=\"$p\" class='T'>" . $el . '</span>';

			}
			elsif ($el =~ m/G/i){
			    $h = "<span title=\"$p\" class='G'>" . $el . '</span>';

			}
			elsif ($el =~ m/C/i){
			    $h = "<span title=\"$p\" class='C'>" . $el . '</span>';

			}
			else{
			    $h = '<span title=\"$p\" class=\'X\'>' . '&nbsp;'  . '</span>';

			}
		    }
		    else{
			if ($el =~ m/A/i){
			    $h =  "<span title=\"$p\" class='different_A'>" . $el . '</span>'; 

			}
			elsif ($el =~ m/T/i){
			    $h =  "<span title=\"$p\" class='different_T'>" . $el . '</span>';

			}
			elsif ($el =~ m/G/i){
			    $h = "<span title=\"$p\" class='different_G'>" . $el . '</span>';

			}
			elsif ($el =~ m/C/i){
			    $h = "<span title=\"$p\" class='different_C'>" . $el . '</span>';

			}
			else{
			    $h = "<span title=\"$p\" class=\'X\'>" . '&nbsp;'  . '</span>';

			}
			


		    }
		}
		push @html, $h if defined $h;



	    }
	    push @html, '<br>',"\n";	    


	}


    }
    return join '', @html;
    #return $rotated;

}

sub html_style{
    my $sheet = $_[1];

    if ($sheet eq 'block'){

    return '
<!--
    body {

        background-color: black;
        font-family: Andale mono, monospace;
        font-size: 14px;
        color: black;  


    }
    .A {
        background-color: #00FFFF;
    }
    .C {
       background-color: #FFFF00;
    }
    .T {
       background-color: #00FF00;
    }
    .G {
       background-color: #FF0000;
    }
    .ref {

       background-color: white;
  
    }
   
-->

';
}
    elsif ($sheet eq 'simple'){

    return '
<!--
    body {

        background-color: black;
        font-family: Andale mono, monospace;
        font-size: 14px;
        color: white;  


    }
    .A {
        color: #00FFFF;
    }
    .C {
       color: #FFFF00;
    }
    .T {
       color: #00FF00;
    }
    .G {
       color: #FF0000;
    }
    .ref {

       color: white;
  
    }
   
-->

';
	

    }
    elsif ($sheet eq 'difference'){
    
	return '<!--
    body {

        background-color: black;
        font-family: Andale mono, monospace;
        font-size: 14px;
        color: black;  


    }
    .A {
        color: #00FFFF;
    }
    .C {
       color: #FFFF00;
    }
    .T {
       color: #00FF00;
    }
    .G {
       color: #FF0000;
    }


    .different_A{
      background-color: #00FFFF;

    }
    .different_C{

     background-color: #FFFF00;

    }
    .different_T {

     background-color: #00FF00;

    }

    .different_G  {

     background-color: #FF0000;


    }
    .ref{

     background-color: white;

    }
   
-->';
	


    }

}

1;
=head1 NAME

Pileup - module that creates object representing maq format pileup files. Will soon include methods for converting bowtie and SOAP
alignments into the pileup format for easy later analysis.

=head1 AUTHOR

Dan MacLean (dan.maclean@tsl.ac.uk)

=head1 SYNOPSIS

	use Solexa::Pileup;
	my $pileup = new Pileup(/home/macleand/Desktop/mypile);


=head1 DESCRIPTION

Pileup is the text based generic, easy to read output from the maq series of alignment programs. The format allows simple representation of
alignments, describing only the depth of coverage at each position on a reference sequence, the reference sequence itself and the 
and the identities of nucleotides that align over that position. It does not record anything to do with quality scores or anything to do
with the actual reads that make up the alignment, like identity or mate pairs.

=head1 METHODS

=over

=item new(pileup_file)

Creates a new Pileup object. The first time you run this method it creates an index of the file. The index allows it to access the different
bits describing each contig without running through the whole file. You can save this index by calling the pickle() method, which stores
This method runs for a long time. Try to use it once, the first time you make a file and then use the pickle() and retrieve_pickle().

	$file = '/home/macleand/Desktop/pileup_file.txt'
	$pileup = new Pileup($file);

=item pickle(outfile)

Creates a pickled (binarised version of the Pileup object) Pileup file and writes it to disk.

	my $outfile = '/home/macleand/Desktop/my_pileup.pickled';
	$pileup->pickle($outfile);

=item retrieve_pickle

Retrieves a pickled Pileup from disk, note that this is not an object method, call in the manner of subroutines, will return a
pileup object. Alternatively you can use the Perl standard module Storable::retrieve();

	my $pickled_file = '/home/macleand/Desktop/my_pileup.pickled';
	my $pileup = Pileup::retrieve_pickled($pickled_file);

	#or ..
	use Storable;
	my $pileup = Storable::retrieve($pickled_file);

=item name_file

Returns the name of the File you used to create this object

	print $pileup->name_file();   #prints '/home/macleand/Desktop/pileup_file.txt'

=item get_coverage(contig, position)

Returns coverage for a single position on a particular contig

	my $depth = $pileup->get_coverage('contig1', 124);

=item get_ref_base

Returns the reference_base at a particular position on a contig 

	my $ref = $pileup->get_ref_base('contig1', 345);

=item get_contigs

Returns list of contigs

	my @contigs = $pileup->get_contigs; 

=item contig_count

Returns the number of contigs in the pileup.

	print $pileup->contig_count();

=item contig_length(contig)

Returns the length of the contig

	print $pileup->contig_length('contig34');

=item nt_coverage(contig, position)

Returns hash with nucleotides as keys and number of times that nucleotide called at a position on a contig

	my %coverage = $pileup->nt_cpverage('contig1', 1234);

=item average_depth

Returns average depth of coverage over all contigs

	print $pileup->average_depth();

=item contig_average_depth(contig)

Returns single contig average coverage depth

	print $pileup->contig_average_depth('contig1234');

=item total_nucleotides

Returns number of nucleotides in all contigs ie sum of lengths

	print $pileup->total_nucleotides();

=item nucleotides_aligned_to_contig(contig) 

Returns the sum of aligned nucleotides at each position in a contig

	print $pileup->nucleotides_aligned_to_contig('contig12');

=item nucleotides_aligned

Returns the sum of nucleotides at each position over all contigs

	print $pileup->nucleotides_aligned();

=item html_pile(contig, direction)

Returns a stringified html version of a single contig for pretty printing. Provide with contig and direction. 'up' will return the
html so that the 'peaks' of the pileup will point up. 'down' will return the same but pointing down. Together these will allow you to
to print and compare two pileups over the same contig in one html document. Use this method in conjuction with the html_style() method
and Perl standard CGI module. 

	my $html = $pileup->html_pile('contig123', 'up');

=item html_style(style)

Returns a string that contains CSS to help render the html string from html_pile(). The style argument can be one of 'simple', 'block' or 'difference'
Simple colours the characters in the string only. Block renders the characters in the string in black and the background in block colour. Difference renders the nucleotides
that agree with the reference as simple, and those that disagree as block. Use in conjunction with the html_pile() method and CGI standard modules

	my $style = $pileup->html_style('difference');

=back

=head1 SEE ALSO
CGI;