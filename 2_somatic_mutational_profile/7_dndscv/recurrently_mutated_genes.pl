#!/usr/bin/perl
#
use strict;
use warnings;

my $root = "/mnt/atgc-d1/drobles/ftalavera/vaping_lung_cancer/somatic_profile/variant_calling/tumour/annotation/";
my %files = ( 
	'MD6753a' => "$root/MD6753a_filtered_mutect2_PASS_selected_vep.vcf",
	'MD6754a' => "$root/MD6754a_filtered_mutect2_PASS_selected_vep.vcf",
	'MD6755a' => "$root/MD6755a_filtered_mutect2_PASS_selected_vep.vcf",
	'MD6756a' => "$root/MD6756a_filtered_mutect2_PASS_selected_vep.vcf",
	'MD6758a' => "$root/MD6758a_filtered_mutect2_PASS_selected_vep.vcf" );

my %mutations;
my @csqs_to_keep = ( 'missense', 'frameshift', 'stop', 'splice' );

#obtain recurrently mutated genes
foreach my $file ( keys %files ){
	open( FILE, $files{ $file } ) or die;
	while( <FILE> ){
		if( $_ =~ /^\#/ ){ next; }
		if( $_ =~ /CSQ\=([^\t]+)/ ){
			my @csqs = split( /\,/, $1 );
			foreach my $csq ( @csqs ){
				my @fields = split( /\|/, $csq );
				my ( $gene, $tr, $csq, $aa ) = ( 'undefined', 'undefined', 'undefined', 'undefined' );
				
				#check if we are interested in the consequence, if not, discard
				if( !grep{ $fields[1] =~ /$_/ } @csqs_to_keep ){ next; }
 
				if( defined( $fields[4] ) ){ #this is the gene name
					$gene = $fields[4];
					if( $fields[3] ne '' ){
						$gene .= " ($fields[3])";
					}
				}
				if( defined( $fields[6] ) ){ #this is the transcript
					$tr = $fields[6];
				}
				if( defined( $fields[1] ) ){ #this is the consequence
					$csq = $fields[1];
				}
				if( defined( $fields[14] ) and defined( $fields[15] ) ){ #this is aa info
					$aa = $fields[14] . ' ' . $fields[15];
				}
				#save in results hash: gene, transcript, sample, consequence and aa
				$mutations{ $gene }{ $tr }{ $file } = $csq . ", " . $aa;
			}
		}
	}
	close( FILE );
}

# Output file
my $output_file = 'tumour_recurrently_mutated_genes.txt';

# Open the file for writing
open my $fh, '>', $output_file or die "Cannot open $output_file: $!";

#print results
foreach my $gene ( sort keys %mutations ){
	foreach my $tr ( sort keys %{ $mutations{ $gene } } ){
		foreach my $sample ( sort keys %{ $mutations{ $gene }{ $tr } } ){
			print $fh ( $gene . "\t" . $tr . "\t" . $sample . "\t" . $mutations{ $gene }{ $tr }{ $sample } . "\n" );
		}
	}
}
	

# Close the file
close $fh;
