package ReadStrDescr; #module reads structure descriptors: dssp secondary structures, vr
use strict;
use warnings;

require Exporter;
use vars qw(@ISA @EXPORT $VERSION);
our @ISA = qw(Exporter);
our @EXPORT = qw(dsspdig vrdig distdig);
$VERSION=1.0;

sub dsspdig{
#selection based on DSSP homepage
#    H = alpha helix
#    B = residue in isolated beta-bridge
#    E = extended strand, participates in beta ladder
#    G = 3-helix (3/10 helix)
#    I = 5 helix (pi helix)
#    T = hydrogen bonded turn
#    S = bend 
#	blanks stands for loops or irregular
#print "$#dssp\n";
my $nl=0;
my $nre=0;
my ($k,$line);
my($dssp,$loop,$ord,$struct)=@_;

	foreach  $line ($$dssp) {
	$nl++;
	if ($line=~/.*RESIDUE\s.*/){	
	for $k ($nl..$$dssp){
	 
 	 $line=$$dssp[$k];	
	 my @temp=split(/\s+/, $line);
	 chomp (@temp);
	 $temp[2]=~ s/(\d+)(\D+)/$1/;
#	 print "$pdb\t$temp[2]\n";
	 my $seqn=$temp[2];
	
	 $nre++;
#	 $chid=substr($line,11,1);
	 my $AA=substr($line,13,1);
	 my $sec=substr($line,16,1);
#wez pod uwage podzial dssp na lancuchy
#printf DSSPOUT "$nre\t$seqn\t$chid\t$AA\t$sec\n";
	 $$ord[$nre]=$seqn;
	 $$struct[$nre]=$AA;
	 if ($seqn =~/.*!.*/){$nre=$nre-1;}
	 if ($sec eq 'H' or $sec eq 'G' or $sec eq 'I'){
	 $$loop[$nre]="H";
	 }
	 elsif ($sec eq 'E'){
	 $$loop[$nre]="P";
	 }
	 elsif ($sec eq 'B' or $sec eq 'T' or $sec eq 'S'){
	 $$loop[$nre]="Z";
	 }
	 else{
	 $$loop[$nre]="L";
	 }
	} 	
	}  
	}
#close (DSSPOUT);
}

sub vrdig{
my (@vr,$j,$k,$line);
my ($fauer)=@_;
#zerowanie
 for $j (0..$$fauer+5){
 $vr[$j][1]=0;	
 $vr[$j][2]=0;	
 $vr[$j][3]=0;}	

 for $j (0..$$fauer){
	$line=$$fauer[$j];
 	my @temp2=split(/\s+/, $line);
 	chomp (@temp2);
 	$k=$temp2[1];	
 	$vr[$k][1]=$temp2[2];	
 	$vr[$k][2]=$temp2[3];	
 	$vr[$k][3]=$temp2[4];}
}

sub distdig{
my (@dist,$j,$k);
my ($rmdist)=@_;
 for $j (1..$$rmdist+1){#zerowanie
 $dist[$j][1]=0;	
 $dist[$j][2]=0;	
 $dist[$j][3]=0;
 $dist[$j][4]=0;
 }

 for $j (1..$$rmdist){
 my $line=$$rmdist[$j];
 chomp $line;
 my @temp2=split(/\s+/,$line);
 $k=$temp2[1];	
 $dist[$k][1]=$temp2[4];	
 $dist[$k][2]=$temp2[5];	
 $dist[$k][3]=$temp2[6];
 $dist[$k][4]=$temp2[7];		
 }
}
1;
