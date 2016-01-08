#!/usr/bin/perl -w
use DBI;
use FindBin;
use strict;

# Constants
package Constants;
use constant MM => 4;
use constant MM_TYPE => 'is';
use constant DB => 'patmatch_2013';
use constant PLANTDB => 'phytozome';

# MYSql
my $host   = "localhost";
my $userid;
my $passwd;

# Arguments
my $pattern	= $ARGV[0] || die "Must give pattern";
my $table_name = $ARGV[1] || die "Must give table name";

$pattern	= uc($pattern);
$pattern	=~ tr/uU/tT/;

my $mismatches = MM;
my $mismatch_types = MM_TYPE;

# Check pattern
my $SYNTAX_CHECKER_BIN = "perl " . $FindBin::Bin . '/patmatchPatternChecker.pl';
my $patStatus = `$SYNTAX_CHECKER_BIN 'dna' $pattern`;
chomp($patStatus);

# Variables
my ($gen_name,$hit_start,$hit_end,$target,$mirna) ='';
my $deltaG = 0;

# Files RNAhybrid target and microrna 5'3' y blast
my $target_file = $FindBin::Bin . "/extra_files/target_rnahybrid.txt";
my $mirna_file = $FindBin::Bin . "/extra_files/mirna_rnahybrid.txt"; 
my $blast_file_sequence = $FindBin::Bin . "/extra_files/" . "seq_blast.txt";
my $blast_database  = $FindBin::Bin . "/extra_files/blast/TAIR10_pep_20101214_updated";

my $sequence_file = '';
my $tab = $table_name . "_" .  $pattern;

hybrid_mirna_file($pattern,$mirna_file);
fill_table_mirnas($pattern,$table_name,$tab);

my $tab_db = create_table($tab);
my @specie_db = species(PLANTDB);
	
foreach my $file (@specie_db){
	
	my $fasta_file = $file->{'fasta'};
	my $specie = $file->{'specie'};
	
	$sequence_file = $FindBin::Bin . '/databases/' . PLANTDB . "/" . $fasta_file;
	print $fasta_file . " - ";
	print $specie . "\n";
	
	my ($gen_sv,$target_sv,$align_sv,$miR_sv,$deltaG_sv,$nro_mm_sv,$ins_sv) = '';
	my ($del_sv,$sust_mm_sv,$gu_mm_sv) = '';
	my ($filtro_mm,$family,$sub_family,$alias) = '';
	
	my %res_blast;
	my @res_ned,@res_mm;
	my @res_family;
	
	if ($patStatus eq "OK") { # syntax OK, run PatMatch
		my $SCAN_PIPELINE = "perl " . $FindBin::Bin . '/scan_pipeline.pl';
		open(OUTPUT, "$SCAN_PIPELINE -c '$pattern' '$sequence_file' '$mismatches' '$mismatch_types' |");
		open(INFO,$sequence_file);
		while (my $line = <OUTPUT>) {
			if ( $line =~ />(.*?):\[(\d*?),(\d*?)\]/ ) {
				#keep output parameters
				$gen_name  = $1;
				$hit_start = $2;
				$hit_end   = $3;
				$gen_sv    = $gen_name;
			}
			elsif ($line =~ /(.*)(\d*?)/ ) {
				$target = trim($1);
				$mirna = $pattern;
				$mirna = reverse complement($mirna);
				# Needleman
				@res_ned = Needleman($mirna,$target);
				$target_sv = $res_ned[0];
				$align_sv  = $res_ned[1];
				$miR_sv  = $res_ned[2];
				# Mismatches
				@res_mm = mismatches($target_sv,$align_sv,$miR_sv);
				$nro_mm_sv  = trim($res_mm[0]);
				$ins_sv  = trim($res_mm[1]);
				$del_sv  = trim($res_mm[2]);
				$sust_mm_sv = trim($res_mm[3]);
				$gu_mm_sv   = trim($res_mm[4]);

				hybrid_target_file($target,$target_file);

				#RNAhybrid				
				my $RNAhybrid_BIN = 'RNAhybrid -d 0 -m 12000';
				my $RNAhybrid_Status = `$RNAhybrid_BIN '-t' $target_file '-q' $mirna_file`;
				chomp($RNAhybrid_Status);
				if ( $RNAhybrid_Status =~ m/mfe: (.*)kcal\/mol/ ){
					$deltaG = $1;
					$deltaG =~ s/ //g;
					$deltaG_sv = $deltaG;
				}
					
				$filtro_mm = mm_position($align_sv);
				insert($tab_db,$specie,$gen_sv,$target_sv,$align_sv,$miR_sv,$nro_mm_sv,$ins_sv,$del_sv,$sust_mm_sv,$gu_mm_sv,$filtro_mm,$deltaG_sv);
			}	
		}
		close(OUTPUT);	
		close(INFO);
		blast_from_db_annotation($fasta_file,$specie,$tab_db);
		family_from_db($fasta_file,$tab_db);
	}
	
	else
	{
		print "Invalid pattern syntax";
	}
}

sub fill_table_mirnas{
	my ($sequence,$name,$table_reference) = @_;
	$target = reverse(complement($sequence));
	hybrid_target_file($target,$target_file);
	my $hyb_perf;
	
	#RNAhybrid				
	my $RNAhybrid_BIN = 'RNAhybrid -d 0 -m 12000';
	my $RNAhybrid_Status = `$RNAhybrid_BIN '-t' $target_file '-q' $mirna_file`;
	chomp($RNAhybrid_Status);
	if ( $RNAhybrid_Status =~ m/mfe: (.*)kcal\/mol/ ){
		$deltaG = $1;
		$deltaG =~ s/ //g;
		$hyb_perf = $deltaG;
	}

	my $sth;
	my $db = DB ;
	my $connectionInfo = "dbi:mysql:$db;$host";
	my $dbh = DBI->connect($connectionInfo,$userid,$passwd);
	
	my $query_insert = "INSERT INTO mirnas
		(name,sequence,table_reference,hyb_perf)
		values
		('$name',
		 '$sequence' ,
		 '$table_reference' ,
		 '$hyb_perf'
		);";

	$sth = $dbh->prepare($query_insert);
	$sth->execute();	
	
}

sub species{
	(my $db_plants) = @_;

	my $db = DB;
	my $connectionInfo="dbi:mysql:$db;$host";
	my $sth;
	my $dbh = DBI->connect($connectionInfo,$userid,$passwd);
	
	my $specie;
	my $fasta_file;
	my $specie_data = {};
	
	my $ref;
	my @res;

	my $query = "SELECT fasta,specie from plants where db = '$db_plants'";

	$sth = $dbh->prepare($query);
	$sth->execute();
	$sth->bind_columns(\$fasta_file,\$specie);

	if ($sth->rows > 0 ){
		while($ref = $sth->fetchrow_hashref() ) {
			push(@res, $ref);
		}
	}

	return @res;
}

sub insert{
	my ($table_db,$especie,$gen,$target,$align,$mirna,$mm,$ins,$del,$sust, $gu,$filtro_mm,$deltag) = @_;
	my $sth;
	my $db = DB ;
	my $connectionInfo = "dbi:mysql:$db;$host";
	my $dbh = DBI->connect($connectionInfo,$userid,$passwd);

	my $query_insert = "INSERT INTO $table_db
				(file,gen,target,align,mirna,mm,ins,del,sust,gu,filtro_mm,deltag)
				values			
				('$especie','$gen','$target','$align','$mirna','$mm','$ins','$del','$sust','$gu','$filtro_mm','$deltag');";

	$sth = $dbh->prepare($query_insert);
	$sth->execute();	

}
sub create_table{

	my ($tabmiRNA) = @_;
	my $sth;
	my $db = DB;
	my $connectionInfo = "dbi:mysql:$db;$host";
	my $dbh = DBI->connect($connectionInfo,$userid,$passwd);

	my $query_create = "CREATE TABLE $tabmiRNA (id INT AUTO_INCREMENT PRIMARY KEY, 
						file VARCHAR(60) NOT NULL,
						gen VARCHAR(30) NOT NULL,
						target  VARCHAR(30) NOT NULL,
						align  VARCHAR(30) NOT NULL,
						mirna  VARCHAR(30) NOT NULL,
						mm INT NOT NULL,
						ins INT NOT NULL,
						del INT NOT NULL,
						sust INT NOT NULL,
						gu INT NOT NULL,
						similar_ath VARCHAR(20) NOT NULL,
						similar_osa VARCHAR(20) NOT NULL,
						filtro_mm INT NOT NULL,
						family TEXT(1000) NOT NULL,
						sub_family VARCHAR(20) NOT NULL,
						alias VARCHAR(20) NOT NULL,
						deltag FLOAT NOT NULL
					);";

	$sth = $dbh->prepare($query_create);
	$sth->execute();
	$sth->finish();
	$dbh->disconnect;
	return $tabmiRNA;
}


sub blast_from_db_annotation{
	my ($file, $specie, $mirna) = @_;

	my $db= DB;
	my $connectionInfo="dbi:mysql:$db;$host";
	my $sth;
	my $dbh = DBI->connect($connectionInfo,$userid,$passwd);
	my $annotation_table = 'annotation_' . $file;
	my $annotation;
	my $query = "SELECT annotation from plants where fasta = '$file'";

	$sth = $dbh->prepare($query);
	$sth->execute();
	$sth->bind_columns(\$annotation);


	$sth->fetch();

	unless($annotation eq 'empty') {
		my $query_update = "UPDATE $mirna miR
							LEFT JOIN $annotation_table a
							ON 	miR.gen	= a.gen
							SET miR.similar_ath = SUBSTRING(a.similar_ath,1,9) ,
								miR.similar_osa = SUBSTRING(a.similar_osa,1,14) 
							WHERE file = '$specie'";
		
		$sth = $dbh->prepare($query_update);
		$sth->execute();
	}
}


sub family{
	(my $gen) = @_;

	my $db = DB;
	my $connectionInfo="dbi:mysql:$db;$host";
	my $sth;
	my $dbh = DBI->connect($connectionInfo,$userid,$passwd);
	my $family;
	my $sub_family;
	my $alias;
	my $query = "SELECT family,sub_family,gene_name from gene_families where locus_tag = '$gen';";

	$sth = $dbh->prepare($query);
	$sth->execute();
	$sth->bind_columns(\$family,\$sub_family,\$alias);

	if ($sth->rows > 0 ){
		while($sth->fetch()) {
			my @results = ($family,$sub_family,$alias);
			return @results;
		}
	}
	else{	
		my @results = ('','','');
		return @results;
	}
}

sub family_from_db{
	(my $file, $mirna) = @_;
	
	my $db = DB;
	my $connectionInfo="dbi:mysql:$db;$host";
	my $sth;
	my $dbh = DBI->connect($connectionInfo,$userid,$passwd);
	my $gen;

	my $query_update = "UPDATE $mirna miR
						LEFT JOIN gene_families f
						ON 	miR.similar_ath = f.locus_tag
						SET miR.family	  = f.family,
							miR.sub_family  = f.sub_family,
							miR.alias	   = f.gene_name" ;
	
	$sth = $dbh->prepare($query_update);
	$sth->execute();	
}

sub hybrid_mirna_file{
	my ($mirna,$file_mirna) = @_;
	open(MIRNA, ">$file_mirna") || die;
	print MIRNA ">mirna\n$mirna\n";
	close(MIRNA);
}

sub hybrid_target_file {
	my ($target,$file_target) = @_;
	open(OUT, ">$file_target") || die;
	print OUT ">target\n$target\n";
	close(OUT);
}

sub mm_position{
	(my $align) = @_;
	my @align = split (//,$align);
	my $i = 0;
	my $j = 0;
	@align = reverse(@align);
	foreach (@align) {
		if ($align[$i] eq "*"){
			if ($i+2 < 13 ){
				$j++;
			}
		}
		$i++;	
	}
	if ($j > 1 ){
		return 0;
	}
	else{
		return 1;
	}

}

sub mismatches{
	my ($target_mm,$align_mm,$mirna_mm) = @_;
	my $pos  = 0;
	my $nro_mm = 0;
	my $del_mm = 0;
	my $ins_mm = 0;
	my $sust_mm = 0;
	my $gu_mm = 0;
	my @align_mm = split(//, $align_mm);
	my @target_mm = split(//, $target_mm);
	my @mirna_mm = split(//, $mirna_mm);	
	# mismatches
	foreach (@align_mm) {
		if  ($_ =~ /\*/) {
			if (($target_mm[$pos] eq 'G' and $mirna_mm[$pos] eq 'T') or ($target_mm[$pos] eq 'T' and $mirna_mm[$pos] eq 'G')){
				$gu_mm++;
			}
			else{}
			$nro_mm ++;
		}
		$pos++;
	}
	# deletions
	foreach (@target_mm) {
		if  ($_ =~ /-/) {
				$del_mm ++;
		}
	}
	# insertions
	foreach (@mirna_mm) {
		if  ($_ =~ /-/) {
				$ins_mm ++;
		}
	
	}
	# sustitutions
	$sust_mm = $nro_mm - ($del_mm + $ins_mm + $gu_mm);
	
	my @results = ($nro_mm,$ins_mm,$del_mm,$sust_mm,$gu_mm);
	return @results;
	
}

sub Needleman {
	my ($seq1,$seq2) = @_;
	# scoring scheme
	my $MATCH	=  1; # +1 for letters that match
	my $MISMATCH = -1; # -1 for letters that mismatch
	my $GAP	  = -1; # -1 for any gap

	# initialization
	my @matrix;
	$matrix[0][0]{score}   = 0;
	$matrix[0][0]{pointer} = "none";

	for(my $j = 1; $j <= length($seq1); $j++) {
		$matrix[0][$j]{score}   = $GAP * $j ;
		$matrix[0][$j]{pointer} = "left";
	}
	for (my $i = 1; $i <= length($seq2); $i++) {
		$matrix[$i][0]{score}   = $GAP * $i ;
		$matrix[$i][0]{pointer} = "up";
	}

	# fill
	for(my $i = 1; $i <= length($seq2); $i++) {
		for(my $j = 1; $j <= length($seq1); $j++) {
			my ($diagonal_score, $left_score, $up_score);

			# calculate match score
			my $letter1 = substr($seq1, $j-1, 1);
			my $letter2 = substr($seq2, $i-1, 1);							
			if ($letter1 eq $letter2) {
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
			}
			else {
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
			}

			# calculate gap scores
			$up_score   = $matrix[$i-1][$j]{score} + $GAP;
			$left_score = $matrix[$i][$j-1]{score} + ($GAP*10);

			# choose best score
			if ($diagonal_score >= $up_score) {
				if ($diagonal_score >= $left_score) {
					$matrix[$i][$j]{score}   = $diagonal_score;
					$matrix[$i][$j]{pointer} = "diagonal";
				}
				else {
					$matrix[$i][$j]{score}   = $left_score;
					$matrix[$i][$j]{pointer} = "left";
                }
			} 
		else {
				if ($up_score >= $left_score) {
                    $matrix[$i][$j]{score}   = $up_score;
                    $matrix[$i][$j]{pointer} = "up";
                }
                else {
                    $matrix[$i][$j]{score}   = $left_score;
                    $matrix[$i][$j]{pointer} = "left";
                }
			}
		}
	}

	# trace-back

	my $align1 = "";
	my $align2 = "";
	my $align3;

	# start at last cell of matrix
	my $j = length($seq1);
	my $i = length($seq2);

	while (1) {
		last if $matrix[$i][$j]{pointer} eq "none"; # ends at first cell of matrix
		my $letter1 = substr($seq1, $j-1, 1);
			my $letter2 = substr($seq2, $i-1, 1);

		if ($matrix[$i][$j]{pointer} eq "diagonal") {
            if ($letter1 eq $letter2){
                $align3 .= "|";
            }
            else
            {
                $align3 .= "*";
            }
			$align1 .= substr($seq1, $j-1, 1);
			$align2 .= substr($seq2, $i-1, 1);
			$i--;
			$j--;


		}
		elsif ($matrix[$i][$j]{pointer} eq "left") {
			$align1 .= substr($seq1, $j-1, 1);
			$align2 .= "-";
            $align3 .= "*";
			$j--;

		}
		elsif ($matrix[$i][$j]{pointer} eq "up") {
			$align1 .= "-";
			$align2 .= substr($seq2, $i-1, 1);
            $align3 .= "*";
			$i--;

		}	
	}

	$align1 = reverse $align1;
	$align2 = reverse $align2;
	$align3 = reverse $align3;
	$align1 = complement($align1);

	my @results = ($align2,$align3,$align1);
	return @results;
}
sub complement
{
	my $sequence = shift;
	$sequence =~ tr/atcgATCG/tagcTAGC/;
	return $sequence;
}

sub trim
{
	(my $string) = @_;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	$string =~ s/\n//g;
	return $string;
}
