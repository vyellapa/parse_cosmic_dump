#! /gsc/bin/perl
################################################################################################################
## parse_cosmic_v70.pl - Parse the cosmic v70 file and create 5 column file
## based on /gscmnt/gc7201/info/medseq/bniu/pancan/HotSpot3D_01022014/cosmic_v67/bin/parse_cosmic.pl 
################################################################################################################

use strict;
use warnings;
use FileHandle;

#USAGE: perl script.pl infile
if(!(-e $ARGV[1]) || !(-e $ARGV[0]))
{
        die "USAGE: perl $0 <COMSIC_DUMP> <COSMIC_VCF>\n";
}

open(COSMIC, "$ARGV[0]") or die "Can't open pass file: $!\n";

while (<COSMIC>) {
	my $vcf_file=$ARGV[1];
	my $max_indel_length=100;
	my $fa="/gscuser/rmashl/ICGC_pilot_data/reference/hs37d5.fa";
        chomp;
        my $line = $_;
	my ($Gene_name,$Accession_Number,$Gene_CDS_length,$HGNC_ID,$Sample_name,$ID_sample,$ID_tumour,$Primary_site,$Site_subtype,$Primary_histology,$Histology_subtype,$Genome_wide_screen,$Mutation_ID,$Mutation_CDS,$Mutation_AA,$Mutation_Description,$Mutation_zygosity,$Mutation_GRCh37_genome_position,$Mutation_GRCh37_strand,$SNP,$FATHMM_score,$FATHMM_prediction,$Mutation_somatic_status,$Pubmed_PMID,$ID_STUDY,$Sample_source,$Tumour_origin,$Age,$Comments ) =split(/\t/, $line);

	# Skip variants of samples that were not Genome-wide screened
	next unless( $Genome_wide_screen eq 'y' );

	# Skip variants belonging to recurrent tumors, metastases, adjacent, etc.
	# We don't want to report the same variant from the same patient more than once
	next unless( $Tumour_origin =~ m/^(primary|secondary|surgery)/ );

	# Skip known germline variants or variants from samples annotated as normals
	next if( $Mutation_somatic_status =~ m/germline/  || $Sample_source eq 'cell-line');

	# Skip mutation types that are too hard to annotate
	next if( $Mutation_Description =~ m/^Whole gene deletion/ );

	# Skip variants that don't have both nucleotide change and loci available
	next unless( $Mutation_CDS and $Mutation_GRCh37_genome_position);
	
	# Skip variants that are not annotatable
	next if( $Mutation_CDS =~m/^c\.\?$/ || $Mutation_CDS =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*>[ACGTacgt]+$/);
	next if( $Mutation_CDS =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*>\d*$/ || $Mutation_CDS =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*[ACTGacgt]+>\d*$/);

	my ( $chr, $start, $stop ) = $Mutation_GRCh37_genome_position =~ m/^(\w+):(\d+)-(\d+)$/;
	die "Cannot identify locus for:\n$line\n" unless( $chr and $start and $stop );
	$chr =~ s/^(23|x)$/X/; $chr =~ s/^(24|y)$/Y/; $chr =~ s/^(25|M|m|mt)$/MT/;
	my $fetched_ref=`samtools faidx $fa $chr:$start-$stop| grep -v ^\\>|tr '\n' ' '|sed 's/ //g'`;
	chomp($fetched_ref);

	# Try to find out what kind of variant this is, and convert it to a minimal 5-column format
    	my ( $ref, $var, @tmp );
    	if( @tmp = $Mutation_CDS =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*([ACGTacgt]+)>([ACGTacgt]+)$/ ) { # SNP/DNP/ONP
        ( $ref, $var ) = ( uc( $tmp[0] ), uc( $tmp[1] ));
    	}
    	elsif( @tmp = $Mutation_CDS =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*ins([ACGTacgt]+)$/ ) { # Insertion
        ( $ref, $var ) = ( "-", uc( $tmp[0] ));
        
        if( length( $var ) > $max_indel_length ) {
            warn "Skipped: Long insertion >$max_indel_length bps in $Mutation_CDS at locus: $Mutation_CDS\n";
            next;
        }
    }
    	

	elsif( @tmp = $Mutation_CDS =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*ins(\d+)$/ ) { # Insertion without sequence specified
		my $var_len=$tmp[0];
		my $cosm_id="COSM".$Mutation_ID;
		my $vcf_line=`cat $vcf_file |grep -w "$chr\t$start\t$cosm_id"|awk -v a=$chr -v p=$start -v c=$cosm_id -F'\t' '{if(a==\$1 && p==\$2 && c=\$3 && length(\$5)>length(\$4) && length(\$4)==1) { print \$1,\$2,\$4,substr(\$5,2,length(\$5));}}'`;
		chomp($vcf_line);
		my $fetch=`samtools faidx $fa $chr:$start-$start| grep -v ^\\>|tr '\n' ' '|sed 's/ //g'`;
		my ($vcf_chr,$vcf_pos,$vcf_ref,$vcf_alt)=split(/ /, $vcf_line);
		if (!$vcf_chr or !$vcf_pos or !$vcf_ref or !$vcf_alt) {
            	warn "Skipped: Insertion unavialable in $vcf_file at locus: $Mutation_GRCh37_genome_position\n";
            	next;
			}
		# Skip if REF in VCF and Fasta dont match
		elsif ($vcf_chr eq $chr and $vcf_pos==$start and length($vcf_alt)==$var_len and $fetch ne $vcf_ref) {
		warn "Skipped: Insertion REF in VCF($vcf_ref) does not match Fatsa REF($fetch) at $chr:$start for $Mutation_CDS at locus: $Mutation_GRCh37_genome_position\n";
		next;
		}
		
	
		elsif ($vcf_chr eq $chr and $vcf_pos==$start and length($vcf_alt)==$var_len and $fetch eq $vcf_ref) {
			($stop,$ref,$var) = (($start+1),"-",$vcf_alt);
			$Mutation_GRCh37_strand="+"; #Since vcf has positive strand only
            warn "Success: Insertion avialable for $Mutation_CDS at locus: $Mutation_GRCh37_genome_position\n";
			}

		else {
		my $altlen=length($vcf_alt);
            	warn "Skipped: Insertion unavialable for(chr:$chr start:$start stop:$stop vcf_ref:$vcf_ref fetched:$fetch vcf_len:$altlen var_len:$var_len) $Mutation_CDS at locus: $Mutation_GRCh37_genome_position\n";
		next;
		}
	}


	elsif( @tmp = $Mutation_CDS =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*del([ACGTacgt]+)$/ ) { # Deletion
	( $ref, $var ) = ( uc( $tmp[0] ), "-" );
		if( length( $ref ) > $max_indel_length ) {
		warn "Skipped: Long deletion >$max_indel_length bps in $Mutation_CDS at locus: $Mutation_GRCh37_genome_position\n";
		next;
	}
	elsif( length( $ref ) != ( $stop - $start + 1 )) {
		my $seq_len=($stop-$start+1);
		warn "Skipped: Length of deleted sequence($seq_len) in $Mutation_CDS doesn't match locus: $Mutation_GRCh37_genome_position\n";
		next;
        }
    }

    elsif( @tmp = $Mutation_CDS =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*del(\d+)$/ ) { # Deletion without sequence specified
        my $del_length = ($tmp[0]);
        if( $del_length > $max_indel_length ) {
            warn "Skipped: Long deletion >$max_indel_length bps in $Mutation_CDS at locus: $Mutation_GRCh37_genome_position\n";
            next;
        }

        if( length( $fetched_ref ) != $del_length ) {
		my $seq_len=($stop-$start+1);
		my $seq_len_v2=length($fetched_ref);
            warn "Skipped: Length of deleted sequence in $Mutation_CDS(ref:$fetched_ref, length:$del_length, $seq_len) doesn't match locus($seq_len_v2): $Mutation_GRCh37_genome_position\n";
            next;
        }
        ( $ref, $var ) = ( uc( $fetched_ref ), "-" );
    }

	unless(( length($ref) == length($var) && length($ref) <= 3 ) ||
           ( $ref eq "-" && $var=~m/[ACGT]+/ ) || ( $ref=~m/[ACGT]+/ && $var eq "-" )) {
        warn "Skipped: Cannot annotate $ref/$var using $Mutation_CDS at locus: $Mutation_GRCh37_genome_position\n";
        next;
    }

    # Reverse complement the ref/var if the variant is from the negative strand
    my $rc_ref = revcomp($ref);
    # When possible, trust the fetched_ref to determine strand. COSMIC's strand info is sometimes wrong
    if( $rc_ref eq $fetched_ref or ( $ref eq "-" and $Mutation_GRCh37_strand eq "-" )) {
        if( $ref eq "-" && $var =~ m/[ACGT]+/ ) {
            $var = revcomp($var);
        }
        elsif( $ref =~ m/[ACGT]+/ && $var eq "-" ) {
            $ref = $rc_ref;
        }
        else {
            $ref = $rc_ref;
            $var = revcomp($var);
        }
    }
	print "$chr	$start	$stop	$ref	$var	$Sample_name	$ID_sample	$ID_tumour	$Mutation_ID\n";

}
close(COSMIC);


sub revcomp {

my $seq=shift;
my $rcomp=reverse($seq);
$rcomp=~ tr/ACGTacgt-/TGCAtgca-/;
return $rcomp;

}
 
