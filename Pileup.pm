package Pileup;
#################################
#Dan MacLean Jun 19th 2009   dan.maclean@tsl.ac.uk
#module of methods for dealing with raw solexa reads .. 
#################################
use strict;
use Exporter;
use Storable qw(nstore);
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
	$$self{'_index'}{$tmp[0]} = $start - length($line);
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
		my $n = $comp{$u};
		$pileup{$tmp[0]}{$tmp[1]}{$n}++;
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

1;
