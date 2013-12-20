package LHBpred;
use strict;
use warnings;

require Exporter;
use vars qw(@ISA @EXPORT $VERSION);
our @ISA = qw(Exporter);
our @EXPORT = qw(lhbdesc);
$VERSION=1.0;

sub lhbdesc {#przypisanie wlasnosci
my ($j,$i,$k);
my ($ord,$aatab,$struct,$hfob_s1,$hfob_s2,$nhfob_s1,$lnorm,$hnorm,$bnorm,$lnum,$hnum,$bnum,$chftablr,$chftab,$nr,$brejk,$lbrejk,$hbrejk,$bbrejk,$motl,$lform,$bform,$hform)=@_;
my (@hfobia1,@hfobia2,@hfobia3,@w_hfobia1,@w_hfobia2,@w_hfobia3,@lscoef,@hscoef,@bscoef,@lsnorm,@bsnorm,@hsnorm);
my (@brejker,@lbrejker,@hbrejker,@bbrejker,@w_brejker,@w_lbrejker,@w_hbrejker,@w_bbrejker);
my (@former,@lformer,@hformer,@bformer,@w_lformer,@w_hformer,@w_bformer);
my (@w_wathbrr,@w_watbbrr,@llhb,@hlhb,@blhb,@p_lstate,@p_hstate,@p_bstate);
my (@lscoef2l,@hscoef2l,@bscoef2l,@lscoef2r,@hscoef2r,@bscoef2r,@lscoef2lr,@hscoef2lr,@bscoef2lr);
for $j(1..$$ord){# iteracja po bialku
$hfobia1[$j]=0;
$hfobia2[$j]=0;
$w_hfobia1[$j]=0;
$w_hfobia2[$j]=0;
$w_hfobia3[$j]=0;
$lscoef[$j]=0;
$bscoef[$j]=0;
$hscoef[$j]=0;
$lsnorm[$j]=0;
$bsnorm[$j]=0;
$hsnorm[$j]=0;
$brejker[$j]=0;
$lbrejker[$j]=0;
$hbrejker[$j]=0;
$bbrejker[$j]=0;
$lformer[$j]=0;
$hformer[$j]=0;
$bformer[$j]=0;
$w_brejker[$j]=0;
$w_lbrejker[$j]=0;
$w_hbrejker[$j]=0;
$w_bbrejker[$j]=0;
$w_lformer[$j]=0;
$w_hformer[$j]=0;
$w_bformer[$j]=0;
$w_wathbrr[$j]=0;
$w_watbbrr[$j]=0;
$llhb[$j]=0;
$hlhb[$j]=0;
$blhb[$j]=0;
$p_lstate[$j]=-2;
$p_hstate[$j]=-2;
$p_bstate[$j]=-2;
#hydrofobowo??
	for $i(0..19){
	if($$struct[$j] eq $$aatab[$i]){
	$hfobia1[$j]=$$hfob_s1[$i];
	$hfobia2[$j]=$$nhfob_s1[$i];
	$hfobia3[$j]=$$hfob_s2[$i];
	$lsnorm[$j]=$$lnorm[$i];
	$lscoef[$j]=$$lnum[$i];
	$hsnorm[$j]=$$hnorm[$i];
	$hscoef[$j]=$$hnum[$i];
	$bsnorm[$j]=$$bnorm[$i];
	$bscoef[$j]=$$bnum[$i];
	 if($j<$$ord){
	 for $k(0..19){
	  if($$struct[$j+1] eq $$aatab[$k]){
	  $nr=20*$i+($k+1);
	  $lscoef2lr[$j]=$$chftablr[$nr][1];
	  $hscoef2lr[$j]=$$chftablr[$nr][4];
	  $bscoef2lr[$j]=$$chftablr[$nr][7];
	  $lscoef2r[$j]=$$chftab[$nr][1];
	  $hscoef2r[$j]=$$chftab[$nr][4];
	  $bscoef2r[$j]=$$chftab[$nr][7];
#	  print "struct[$j],struct[$j+1],$lscoef2lr[$j],$hscoef2lr[$j],$bscoef2lr[$j],$lscoef2r[$j],$hscoef2r[$j],$bscoef2r[$j]\n";
	  }
	 }}
#print "$pdb\t$j\t$nr\t$struct[$j]\t$lscoef2lr[$j]\t$hscoef2lr[$j]\t$bscoef2lr[$j]\n";
	 if($j>1){
	 for $k(0..19){
	  if($$struct[$j-1] eq $$aatab[$k]){
	  $nr=20*$k+($i+1);
	  $lscoef2l[$j]=$$chftab[$nr][1];
	  $hscoef2l[$j]=$$chftab[$nr][4];
	  $bscoef2l[$j]=$$chftab[$nr][7];
#	  print "struct[$j],struct[$j-1],$lscoef2l[$j],$hscoef2l[$j],$bscoef2l[$j]\n";	  
	  }
	 }}
	last;}
	}
$lscoef2lr[$j]=$lscoef[$j];
$hscoef2lr[$j]=$hscoef[$j];
$bscoef2lr[$j]=$bscoef[$j];
$lscoef2l[$j]=$lscoef[$j];
$hscoef2l[$j]=$hscoef[$j];
$bscoef2l[$j]=$bscoef[$j];
$lscoef2r[$j]=$lscoef[$j];
$hscoef2r[$j]=$hscoef[$j];
$bscoef2r[$j]=$bscoef[$j];
#szukaj brejkerow po kolei
	for $i(0..1){
	if($$struct[$j] eq $$brejk[$i]){
	$brejker[$j]=1;
	last;}
	}
	for $i(0..7){
	if($$struct[$j] eq $$lbrejk[$i]){
	$lbrejker[$j]=1;
	last;}
	}
	for $i(0..4){
	if($$struct[$j] eq $$hbrejk[$i]){
	$hbrejker[$j]=1;
	last;}
	}
	for $i(0..7){
	if($$struct[$j] eq $$bbrejk[$i]){
	$bbrejker[$j]=1;
	last;}
	}
#szukaj formerów po kolei
	for $i(0..5){
	if($$struct[$j] eq $$lform[$i]){
	$lformer[$j]=1;
	last;}
	}
	for $i(0..7){
	if($$struct[$j] eq $$hform[$i]){
	$hformer[$j]=1;
	last;}
	}
	for $i(0..6){
	if($$struct[$j] eq $$bform[$i]){
	$bformer[$j]=1;
	last;}
	}
}#koniec iteracji po bialku

#analiza w oknie
my $st= int($motl/2);
my $stop=$$ord-$st;
my $lewy=$st-1;
my $prawy=$st;
if($st*2 < $motl){#dla nieparzystych dodaj jeden i bedzie srodek motywu
$st=$st+1;
$stop=$$ord-$st+1;
$lewy=$st-1;
$prawy=$st-1;}
for $j($st..$stop){# iteracja po bialku z uwzglednieniem szerokosci okna
	for $k($j-$lewy..$j+$prawy){#suma po oknie
#print "$j,$k,$hfobia1[$k]\n";
	$w_hfobia1[$j]=$w_hfobia1[$j]+$hfobia1[$k];#zageszczenie hydrofobowosci
	$w_hfobia2[$j]=$w_hfobia2[$j]+$hfobia2[$k];#zageszczenie hydrofobowosci
	$w_hfobia3[$j]=$w_hfobia3[$j]+$hfobia3[$k];#zageszczenie hydrofobowosci
	}
	#srednia
	$w_hfobia1[$j]=$w_hfobia1[$j]/$motl;
	$w_hfobia2[$j]=$w_hfobia2[$j]/$motl;
	$w_hfobia3[$j]=$w_hfobia3[$j]/$motl;
}

for $j(2..$$ord-1){
 #brejkery
 if($brejker[$j]==1){#zageszczenie brejkerow
 $w_brejker[$j]=$brejker[$j]+$brejker[$j+1]+$brejker[$j-1]}
 if($hbrejker[$j]==1){#zageszczenie helix brejkerow
 $w_hbrejker[$j]=$hbrejker[$j]+$hbrejker[$j+1]+$hbrejker[$j-1]}
 if($bbrejker[$j]==1){#zageszczenie beta brejkerow
 $w_bbrejker[$j]=$bbrejker[$j]+$bbrejker[$j+1]+$bbrejker[$j-1]}
 if($lbrejker[$j]==1){#zageszczenie loop brejkerow
 $w_lbrejker[$j]=$lbrejker[$j]+$lbrejker[$j+1]+$lbrejker[$j-1]}
 #formery
 if($hformer[$j]==1){#zageszczenie helix formerów
 $w_hbrejker[$j]=$hformer[$j]+$hformer[$j+1]+$hformer[$j-1]}
 if($bformer[$j]==1){#zageszczenie beta formerów
 $w_bformer[$j]=$bformer[$j]+$bformer[$j+1]+$bformer[$j-1]}
 if($lformer[$j]==1){#zageszczenie loop formerów
 $w_lformer[$j]=$lformer[$j]+$lformer[$j+1]+$lformer[$j-1]}
}
#funkcja przewiduj?ca
for $j(2..$$ord-1){
if($lformer[$j]==1){
$p_lstate[$j]=-0.7;}
if($bformer[$j]==1){
$p_bstate[$j]=-0.8;}
if($hformer[$j]==1){
$p_hstate[$j]=-0.9;}
#basic + normalized hydrophob
#$llhb[$j]=$w_hfobia2[$j]+($lscoef[$j]*$lscoef[$j]);
#$hlhb[$j]=$w_hfobia2[$j]+($hscoef[$j]*$hscoef[$j]);
#$blhb[$j]=$w_hfobia2[$j]+($bscoef[$j]*$bscoef[$j]);
#dependent probability + normalized hydrophob
#print "$pdb\t$j\t$struct[$j]\t$lscoef2lr[$j]\t$hscoef2lr[$j]\t$bscoef2lr[$j]\n";
#$llhb[$j]=$w_hfobia2[$j]+($lscoef2lr[$j]*$lscoef2lr[$j]);
#$hlhb[$j]=$w_hfobia2[$j]+($hscoef2lr[$j]*$hscoef2lr[$j]);
#$blhb[$j]=$w_hfobia2[$j]+($bscoef2lr[$j]*$bscoef2lr[$j]);
#dependent probability + normalized hydrophob
#v0
#$llhb[$j]=$lscoef2l[$j]*$lscoef2r[$j];
#$hlhb[$j]=$w_hfobia2[$j]+($hscoef2l[$j]*$hscoef2r[$j]);
#$blhb[$j]=$w_hfobia2[$j]+($bscoef2l[$j]*$bscoef2r[$j]);
#v1
#$llhb[$j]=(($lscoef2l[$j]*$lscoef2l[$j])+($lscoef2r[$j]*$lscoef2r[$j]))/2;
#$hlhb[$j]=$w_hfobia2[$j]+(($hscoef2l[$j]*$hscoef2l[$j])+($hscoef2r[$j]*$hscoef2r[$j]))/2;
#$blhb[$j]=$w_hfobia2[$j]+(($bscoef2l[$j]*$bscoef2l[$j])+($bscoef2r[$j]*$bscoef2r[$j]))/2;
#v2
#$llhb[$j]=($lscoef2l[$j]*$lscoef2l[$j])*($lscoef2r[$j]*$lscoef2r[$j]);
#$hlhb[$j]=$w_hfobia2[$j]+($hscoef2l[$j]*$hscoef2l[$j])*($hscoef2r[$j]*$hscoef2r[$j]);
#$blhb[$j]=$w_hfobia2[$j]+($bscoef2l[$j]*$bscoef2l[$j])*($bscoef2r[$j]*$bscoef2r[$j]);
#v4
#$llhb[$j]=$lscoef2l[$j]*$lscoef2r[$j]+$wathbrr[$j]+$watbbrr[$j];
#$hlhb[$j]=$w_hfobia2[$j]-$wathbrr[$j]+($hscoef2l[$j]*$hscoef2r[$j]);
#$blhb[$j]=$w_hfobia2[$j]-$watbbrr[$j]+($bscoef2l[$j]*$bscoef2r[$j]);
#v5
#$llhb[$j]=$w_hfobia1[$j]-($lscoef2l[$j]*$lscoef2r[$j]);
#$hlhb[$j]=$w_hfobia1[$j]-($hscoef2l[$j]*$hscoef2r[$j]);
#$blhb[$j]=$w_hfobia1[$j]-($bscoef2l[$j]*$bscoef2r[$j]);
#v6
#$llhb[$j]=sqrt (1/$w_hfobia3[$j])+($lscoef2l[$j]*$lscoef2r[$j]);
#$hlhb[$j]=sqrt (1/$w_hfobia3[$j])+($hscoef2l[$j]*$hscoef2r[$j]);
#$blhb[$j]=sqrt (1/$w_hfobia3[$j])+($bscoef2l[$j]*$bscoef2r[$j]);
#v7
$llhb[$j]=(1/$w_hfobia3[$j])+($lscoef2l[$j]*$lscoef2r[$j]);
$hlhb[$j]=(1/$w_hfobia3[$j])+($hscoef2l[$j]*$hscoef2r[$j]);
$blhb[$j]=(1/$w_hfobia3[$j])+($bscoef2l[$j]*$bscoef2r[$j]);
#$slhb[$j]=$llhb[$j]+
}

}
1;
