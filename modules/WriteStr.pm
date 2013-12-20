package WriteStr;
use strict;
use warnings;

require Exporter;
use vars qw(@ISA @EXPORT $VERSION);
our @ISA = qw(Exporter);
our @EXPORT = qw(pdbprint);
$VERSION=1.0;

sub pdbprint{
printf GOUT "reading pdb files...\n";
$katalog="$rootdir/$kat_pdb";
chdir "$katalog" or die "Nie moge wejsc do $kat_pdb: $!\n";
printf GOUT "changing directory to: ../$kat_pdb\n";

#check if all files are present, list of *.pdb files and *.dssp files should match
$matched=0;
 foreach $pdbn (@lipdb){#iteracja po plikach pdb
 $pat1=substr($pdbn,0,7);
 $pat1=~tr/A-Z/a-z/;
 $pat2=substr($dname,0,7);
 $pat2=~tr/A-Z/a-z/;
	if($pat1 eq $pat2){
	$matched=1;
	$pdbname=$pdbn;
	last;
	}
 }
  	if ($matched==0){
	printf GOUT "no PDB file found which matches dssp file: $dname: $!\n";
	next;
	}

printf GOUT "pdb file analyzed:\t$pdbname\n";
open(PDBINP, "<$pdbname") or die "Can’t open input file: $!";
my @pidibi=<PDBINP>;
chomp @pidibi;
close (PDBINP);

$pr=0;
$nr=0;
$na=0;
	foreach $dane(@pidibi){
		$info=substr($dane, 0, 4);
		if ($info eq "ATOM"){
		$reszta=substr($dane, 22, 8);
#		print "$pdb\t$pr\t\t$reszta\t$nr\t$na\t$at[$nr]\n";
#		$reszta=~s/(\d+)(\D+)/$1/}
			if($reszta ne $pr){
		 	$at[$nr]=$na;
#		print "$pdb\t$pr\t\t$reszta\t$nr\t$na\t$at[$nr]\n";
			$pr=$reszta;
			$nr++;
		 	$na=0;}
		$na++;			
		$pdbl[$nr][$na]=$dane;}
	}
#	$nr++;
	$at[$nr]=$na;#dla ostatniego
}
1;