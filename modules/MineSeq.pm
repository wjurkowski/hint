package MineSeq;
use strict;
use warnings;

require Exporter;
use vars qw(@ISA @EXPORT $VERSION);
our @ISA = qw(Exporter);
our @EXPORT = qw(watpatt defchoufas sumchoufas seqprint);
$VERSION=1.0;

sub watpatt {#prints feature patterns in sequence
my (@proper,@sproper,@proper2,@wathbrr,@watbbrr,$j,$i);
my ($struct,$polar,$spolar,$plusk,$minusk,$hydrof,$shydrof,$ord,$nn1,$nn2,$nn3,$nn4,$nn21,$nn22,$nn23,$nn24)=@_;
for $j(1..$$ord){
 	$proper[$j]=0;
 	$proper2[$j]=0;
	$sproper[$j]=0;
 #znajdz AA naladowane i hydrofobowe
 	
	for $i(0..$$polar){
 	 if($$struct[$j] eq $polar[$i]){#jezeli naladowane ujemnie
 	 $proper[$j]='P';}
 	}
	for $i(0..$$spolar){
 	 if($$struct[$j] eq $spolar[$i]){#jezeli naladowane ujemnie
 	 $sproper[$j]='P';}
 	}
 	for $i(0..$$plusk){
	 if($$struct[$j] eq $plusk[$i]){#jezeli naladowane dodatnio 
	 $proper2[$j]='+';}
	}
	for $i(0..$$minusk){
	 if($$struct[$j] eq $minusk[$i]){#jezeli naladowane ujemnie
	 $proper2[$j]='-';}
 	}
 	for $i(0..$$hydrof){
 	 if($$struct[$j] eq $hydrof[$i]){#jezeli hydrofobowe
 	 $proper[$j]='H';}
 	}
	 for $i(0..$$shydrof){
 	 if($$struct[$j] eq $shydrof[$i]){#jezeli hydrofobowe
 	 $sproper[$j]='H';}
 	}
  }
	for $j($$ord+1..$$ord+4){#dodanie 4 pseudoreszt na koncu (do poszukiwan wzorcow wlasnosci)
	$proper2[$j]=' ';
	$proper[$j]=' ';
	$sproper[$j]=' ';}

print PATTERNS "PDB code\torder\tresid\tAA\tn--n+1\tn--n+2\tn--n+3\tn--n+4\tn--n+1\tn--n+2\tn--n+3\tn--n+4\n";
#print PATTERNS2 "PDB code\torder\tresid\tAA\tn--n+1\tn--n+2\tn--n+3\tn--n+4\n";
 for $j(1..$$ord){# iteracja po bialka
#wypelniaj pustkom
$$nn1[$j]=' ';
$$nn2[$j]=' ';
$$nn3[$j]=' ';
$$nn4[$j]=' ';
$$nn21[$j]=' ';
$$nn22[$j]=' ';
$$nn23[$j]=' ';
$$nn24[$j]=' ';
$wathbrr[$j]=0;
$watbbrr[$j]=0;
	if($proper[$j] eq 'P'){
	 if($proper[$j+1] eq 'H'){
	 $$nn1[$j]='PH';}
	 if($proper[$j+2] eq 'H'){
	 $$nn2[$j]='PH';}
	 if($proper[$j+3] eq 'H'){
	 $$nn3[$j]='PH';}
	 if($proper[$j+4] eq 'H'){
	 $$nn4[$j]='PH';}
	}
	if($proper[$j] eq 'H'){
	 if($proper[$j+1] eq 'P'){
	 $$nn1[$j]='HP';}
	 if($proper[$j+2] eq 'P'){
	 $$nn2[$j]='HP';}
	 if($proper[$j+3] eq 'P'){
	 $$nn3[$j]='HP';}
	 if($proper[$j+4] eq 'P'){
	 $$nn4[$j]='HP';}
	}
	if($proper2[$j] eq '+'){
#	print "check plus: $proper2[$j+1]\t$proper2[$j+2]\t$proper2[$j+3]\t$proper2[$j+4]\n";
	 if($proper2[$j+1] eq '+'){
	 $$nn21[$j]='++';}
	 if($proper2[$j+1] eq '-'){
	 $$nn21[$j]='+-';}
	 if($proper2[$j+2] eq '+'){
	 $$nn22[$j]='++';}
	 if($proper2[$j+2] eq '-'){
	 $$nn22[$j]='+-';}
	 if($proper2[$j+3] eq '+'){
	 $$nn23[$j]='++';}
	 if($proper2[$j+3] eq '-'){
	 $$nn23[$j]='+-';}
	 if($proper2[$j+4] eq '+'){
	 $$nn24[$j]='++';}
	 if($proper2[$j+4] eq '-'){
	 $$nn24[$j]='+-';}
	}
	if($proper2[$j] eq '-'){
#	print "check minus: $proper2[$j+1]\t$proper2[$j+2]\t$proper2[$j+3]\t$proper2[$j+4]\n";
	 if($proper2[$j+1] eq '-'){
	 $$nn21[$j]='--';}
	 if($proper2[$j+1] eq '+'){
	 $$nn21[$j]='-+';}
	 if($proper2[$j+2] eq '-'){
	 $$nn22[$j]='--';}
	 if($proper2[$j+2] eq '+'){
	 $$nn22[$j]='-+';}
	 if($proper2[$j+3] eq '-'){
	 $$nn23[$j]='--';}
	 if($proper2[$j+3] eq '+'){
	 $$nn23[$j]='-+';}
	 if($proper2[$j+4] eq '-'){
	 $$nn24[$j]='--';}
	 if($proper2[$j+4] eq '+'){
	 $$nn24[$j]='-+';}
	}
#water helix breakers
	if($sproper[$j] eq 'H'){
	 if($sproper[$j+3] eq 'P'){
	  if($sproper[$j+6] eq 'H'){
	 $wathbrr[$j]=$wathbrr[$j]+1.0;}
	 }
	 if($sproper[$j+4] eq 'P'){
	  if($sproper[$j+8] eq 'H'){
	 $wathbrr[$j]=$wathbrr[$j]+1.0;}
	 }
	}
	if($sproper[$j] eq 'P'){
	 if($sproper[$j+3] eq 'H'){
	  if($sproper[$j+6] eq 'P'){
	 $wathbrr[$j]=$wathbrr[$j]+1.0;}
	 }
	 if($sproper[$j+4] eq 'H'){
	  if($sproper[$j+8] eq 'P'){
	 $wathbrr[$j]=$wathbrr[$j]+1.0;}
	 }
	}
#water beta breakers
	if($sproper[$j] eq 'H'){
	 if($sproper[$j+2] eq 'P'){
	  if($sproper[$j+4] eq 'H'){
	 $watbbrr[$j]=$watbbrr[$j]+1.0;}
	 }
	}
	if($sproper[$j] eq 'P'){
	 if($sproper[$j+2] eq 'H'){
	  if($sproper[$j+4] eq 'P'){
	 $watbbrr[$j]=$watbbrr[$j]+1.0;}
	 }
	}
 }
}#end of watpatt

sub defchoufas{#calculate chou-fasman coefficients

my ($i,$j);
my (@l_c_peika,@h_c_peika,@b_c_peika,@l_b_peika,@h_b_peika,@b_b_peika,@lpeika,@hpeika,@bpeika);
my ($lsum,$hsum,$bsum,$lressum,$hressum,$bressum,$lr2sum,$hr2sum,$br2sum,$Plka,$Phka,$Pbka,$Plka_c,$Phka_c,$Pbka_c,$Plka_b,$Phka_b,$Pbka_b)=@_;
$nsum=$lsum+$hsum+$bsum;
for $i (0..19){#poczatkowe
$nressum[$i]=$lressum[$i]+$hressum[$i]+$bressum[$i];#suma zliczenia kolejnych reszt we wszystkich motywach
	for $j (0..19){#poczatkowe
	$nr2sum[$i][$j]=$lr2sum[$i][$j]+$hr2sum[$i][$j]+$br2sum[$i][$j];}	
$lpeika[$i]=0;
$l_c_peika[$i]=0;
$l_b_peika[$i]=0;
$hpeika[$i]=0;
$h_c_peika[$i]=0;
$h_b_peika[$i]=0;
$bpeika[$i]=0;
$b_c_peika[$i]=0;
$b_b_peika[$i]=0;
$Plka[$i]=0;
$Phka[$i]=0;
$Pbka[$i]=0;
$Plka_c[$i]=0;
$Phka_c[$i]=0;
$Pbka_c[$i]=0;
$Plka_b[$i]=0;
$Phka_b[$i]=0;
$Pbka_b[$i]=0;
}

for $i (0..19){#obliczanie prawdopodobienstw
$pei[$i]=$nressum[$i]/$nsum;#prawdopodobienstwo wystapienia danego aminokwasu w bialku
 for $j (0..19){
 $pei2[$i][$j]=$nr2sum[$i][$j]/$nsum;}#j.w. pod warunkiem znalezienia kolejnego aminokwasu
 if($lsum >0){#w petli
 $lpeika[$i]=$lressum[$i]/$lsum;#p znalezienia danej reszty we wszystkich fragmentach petlowych
  for $j (0..19){
  $lpeika2[$i][$j]=$lr2sum[$i][$j]/$lsum;}}#p znalezienia danej reszty pod warunkiem znalezienia kolejnej we wszystkich fragmentach petlowych
 if($l_c_sum>0){
 $l_c_peika[$i]=$l_c_ressum[$i]/$l_c_sum;# j.w. dla core 
  for $j (0..19){
  $l_c_peika2[$i][$j]=$l_c_r2sum[$i][$j]/$l_c_sum;}}
 if($l_b_sum>0){
 $l_b_peika[$i]=$l_b_ressum[$i]/$l_b_sum;#j.w. dla brzegow
  for $j (0..19){
  $l_b_peika2[$i][$j]=$l_b_r2sum[$i][$j]/$l_b_sum;}}
 if($hsum >0){#w helisie analogicznie jak w petli
 $hpeika[$i]=$hressum[$i]/$hsum;
  for $j (0..19){
  $hpeika2[$i][$j]=$hr2sum[$i][$j]/$hsum;}}
 if($h_c_sum>0){
 $h_c_peika[$i]=$h_c_ressum[$i]/$h_c_sum;
  for $j (0..19){	
  $h_c_peika2[$i][$j]=$h_c_r2sum[$i][$j]/$h_c_sum;}}
 if($h_b_sum>0){
 $h_b_peika[$i]=$h_b_ressum[$i]/$h_b_sum;
  for $j (0..19){
  $h_b_peika2[$i][$j]=$h_b_r2sum[$i][$j]/$h_b_sum;}}
 if($bsum>0){#w beta analogicznie jak w petli
 $bpeika[$i]=$bressum[$i]/$bsum;
  for $j (0..19){
  $bpeika2[$i][$j]=$br2sum[$i][$j]/$bsum;}}
 if($b_c_sum>0){
 $b_c_peika[$i]=$b_c_ressum[$i]/$b_c_sum;
  for $j (0..19){
  $b_c_peika2[$i][$j]=$b_c_r2sum[$i][$j]/$b_c_sum;}}
 if($b_b_sum>0){
 $b_b_peika[$i]=$b_b_ressum[$i]/$b_b_sum;
  for $j (0..19){
  $b_b_peika2[$i][$j]=$b_b_r2sum[$i][$j]/$b_b_sum;}}
}
for $i (0..19){#obliczanie wspolczynnikow
if($pei[$i]>0){
#print "ffff $hpeika[$i],$h_c_peika[$i],$h_b_peika[$i],$pei[$i]\n";
if($lsum >0){
$Plka[$i]=$lpeika[$i]/$pei[$i];
$Plka_b[$i]=$l_b_peika[$i]/$pei[$i];
$Plka_c[$i]=$l_c_peika[$i]/$pei[$i];
 for $j (0..19){
 if($pei2[$i][$j]>0){	
 $Plka2[$i][$j]=$lpeika2[$i][$j]/$pei2[$i][$j];
 $Plka_b2[$i][$j]=$l_b_peika2[$i][$j]/$pei2[$i][$j];
 $Plka_c2[$i][$j]=$l_c_peika2[$i][$j]/$pei2[$i][$j];}}
}
if($hsum >0){
$Phka[$i]=$hpeika[$i]/$pei[$i];
$Phka_b[$i]=$h_b_peika[$i]/$pei[$i];
$Phka_c[$i]=$h_c_peika[$i]/$pei[$i];
 for $j (0..19){
 if($pei2[$i][$j]>0){
 $Phka2[$i][$j]=$hpeika2[$i][$j]/$pei2[$i][$j];
 $Phka_b2[$i][$j]=$h_b_peika2[$i][$j]/$pei2[$i][$j];
 $Phka_c2[$i][$j]=$h_c_peika2[$i][$j]/$pei2[$i][$j];}}
}
if($bsum >0){
$Pbka[$i]=$bpeika[$i]/$pei[$i];
$Pbka_b[$i]=$b_b_peika[$i]/$pei[$i];
$Pbka_c[$i]=$b_c_peika[$i]/$pei[$i];
 for $j (0..19){
 if($pei2[$i][$j]>0){	
 $Pbka2[$i][$j]=$bpeika2[$i][$j]/$pei2[$i][$j];
 $Pbka_b2[$i][$j]=$b_b_peika2[$i][$j]/$pei2[$i][$j];
 $Pbka_c2[$i][$j]=$b_c_peika2[$i][$j]/$pei2[$i][$j];}}
}
}
}
# dla analizy zbiorczej wszystkich plikow
$nsumtot=$nsumtot+$nsum;
$lsumtot=$lsumtot+$lsum;
$l_c_sumtot=$l_c_sumtot+$l_c_sum;
$l_b_sumtot=$l_b_sumtot+$l_b_sum;
$hsumtot=$hsumtot+$hsum;
$h_c_sumtot=$h_c_sumtot+$h_c_sum;
$h_b_sumtot=$h_b_sumtot+$h_b_sum;
$bsumtot=$bsumtot+$bsum;
$b_c_sumtot=$b_c_sumtot+$b_c_sum;
$b_b_sumtot=$b_b_sumtot+$b_b_sum;

for $i (0..19){
$nressumtot[$i]=$nressumtot[$i]+$nressum[$i];
$lressumtot[$i]=$lressumtot[$i]+$lressum[$i];
$l_c_ressumtot[$i]=$l_c_ressumtot[$i]+$l_c_ressum[$i];
$l_b_ressumtot[$i]=$l_b_ressumtot[$i]+$l_b_ressum[$i];
$hressumtot[$i]=$hressumtot[$i]+$hressum[$i];
$h_c_ressumtot[$i]=$h_c_ressumtot[$i]+$h_c_ressum[$i];
$h_b_ressumtot[$i]=$h_b_ressumtot[$i]+$h_b_ressum[$i];
$bressumtot[$i]=$bressumtot[$i]+$bressum[$i];
$b_c_ressumtot[$i]=$b_c_ressumtot[$i]+$b_c_ressum[$i];
$b_b_ressumtot[$i]=$b_b_ressumtot[$i]+$b_b_ressum[$i];
 for $j (0..19){
 $nr2sumtot[$i][$j]=$nr2sumtot[$i][$j]+$nr2sum[$i][$j];
 $lr2sumtot[$i][$j]=$lr2sumtot[$i][$j]+$lr2sum[$i][$j];
 $l_c_r2sumtot[$i][$j]=$l_c_r2sumtot[$i][$j]+$l_c_r2sum[$i][$j];
 $l_b_r2sumtot[$i][$j]=$l_b_r2sumtot[$i][$j]+$l_b_r2sum[$i][$j];
 $hr2sumtot[$i][$j]=$hr2sumtot[$i][$j]+$hr2sum[$i][$j];
 $h_c_r2sumtot[$i][$j]=$h_c_r2sumtot[$i][$j]+$h_c_r2sum[$i][$j];
 $h_b_r2sumtot[$i][$j]=$h_b_r2sumtot[$i][$j]+$h_b_r2sum[$i][$j];
 $br2sumtot[$i][$j]=$br2sumtot[$i][$j]+$br2sum[$i][$j];
 $b_c_r2sumtot[$i][$j]=$b_c_r2sumtot[$i][$j]+$b_c_r2sum[$i][$j];
 $b_b_r2sumtot[$i][$j]=$b_b_r2sumtot[$i][$j]+$b_b_r2sum[$i][$j];
 }
}
}#end of defchoufas

sub seqprint{#printout sequences

my ($seqname,$llin)=@_;
my ($k1,$k2);
open(SEQ, ">$seqname") or die "Can’t open output file: $! $wyniki/$sequences/$seqname";
printf SEQ ">$pdb\n";

	for $k(1..$llin+1){# wydruk
	$k1=($k-1)*80+1;
	$k2=$k*80;
	for $j($k1..$k2){# wydruk
	if($j <= $#ord){
	printf SEQ ("%s",$struct[$j]);
#	print "$struct[$j]";
	}}
#	print "\n";
	print SEQ "\n";	
	for $m(0..7){
	$l=$k1+$m*10;
	if($l <= $#ord) {
#	print "$ord[$l]";
	printf SEQ ("|%-9d",$ord[$l]);}
	}
#	print "\n";
	print SEQ "\n";
#	for ($j=$k1;$j=$k2;$j=$j+5){# wydruk
#	for $j($k1..$k2;5){# wydruk
#	if($j <= $#ord){printf SEQ ("|%-9d",$ord[$j]);}}
		
#	printf SEQ ("|%-9d",$ord[]);
#	for $m(1..8){
#	$poz=$m*10*$k;
#print "$poz\n";
#	if($poz <= $#ord){printf SEQ ("|%-9d",$ord[$poz]);}}
	}

#	print SEQ "\n";
#	for $j(1..$#ord){# wydruk
#	 if($lstate[$j] !=0){
#	 print SEQ "$lstate[$j]";}
#	 if($hstate[$j] !=0){
#	 print SEQ "$hstate[$j]";}
#	 if($bstate[$j] !=0){
#	 print SEQ "$bstate[$j]";}
#	}
#	print SEQ "\n";
close (SEQ);	
}

sub sumchoufas{#final derivation of chou-fasman parameters

my $nsumtot=0;
for $i (0..19){#obliczanie prawdopodobienstw
$pei[$i]=$nressumtot[$i]/$nsumtot;
 for $j (0..19){
 $pei2[$i][$j]=$nr2sumtot[$i][$j]/$nsumtot;}

 if($lsumtot>0){#loops
 $lpeika[$i]=$lressumtot[$i]/$lsumtot;
  for $j (0..19){
  $lpeika2[$i][$j]=$lr2sumtot[$i][$j]/$lsumtot;}}
 if($l_c_sumtot>0){
 $l_c_peika[$i]=$l_c_ressumtot[$i]/$l_c_sumtot;
  for $j (0..19){
  $l_c_peika2[$i][$j]=$l_c_r2sumtot[$i][$j]/$l_c_sumtot;}}
 if($l_b_sumtot>0){
 $l_b_peika[$i]=$l_b_ressumtot[$i]/$l_b_sumtot;
  for $j (0..19){
  $l_b_peika2[$i][$j]=$l_b_r2sumtot[$i][$j]/$l_b_sumtot;}}
 if($hsumtot>0){#helices
 $hpeika[$i]=$hressumtot[$i]/$hsumtot;
  for $j (0..19){
  $hpeika2[$i][$j]=$hr2sumtot[$i][$j]/$hsumtot;}}
 if($h_c_sumtot>0){
 $h_c_peika[$i]=$h_c_ressumtot[$i]/$h_c_sumtot;
  for $j (0..19){
  $h_c_peika2[$i][$j]=$h_c_r2sumtot[$i][$j]/$h_c_sumtot;}}
 if($h_b_sumtot>0){
 $h_b_peika[$i]=$h_b_ressumtot[$i]/$h_b_sumtot;
  for $j (0..19){
  $h_b_peika2[$i][$j]=$h_b_r2sumtot[$i][$j]/$h_b_sumtot;}}
 if($bsumtot>0){#betas
 $bpeika[$i]=$bressumtot[$i]/$bsumtot;
  for $j (0..19){
  $bpeika2[$i][$j]=$br2sumtot[$i][$j]/$bsumtot;}}
 if($b_c_sumtot>0){
 $b_c_peika[$i]=$b_c_ressumtot[$i]/$b_c_sumtot;
  for $j (0..19){
  $b_c_peika2[$i][$j]=$b_c_r2sumtot[$i][$j]/$b_c_sumtot;}}
 if($b_b_sumtot>0){
 $b_b_peika[$i]=$b_b_ressumtot[$i]/$b_b_sumtot;
  for $j (0..19){
  $b_b_peika2[$i][$j]=$b_b_r2sumtot[$i][$j]/$b_b_sumtot;}}
}
print GOUT "chou-fasman like parms for given AA\n";
print GOUT "AA\tP hincz\tcore\tborder\tP helix\tcore\tborder\tP beta\tcore\tborder\n";
for $i (0..19){#obliczanie wspolczynnikow
if($pei[$i]>0){
 if($lsumtot>0){
 $Plka[$i]=$lpeika[$i]/$pei[$i];
 $Plka_b[$i]=$l_b_peika[$i]/$pei[$i];
 $Plka_c[$i]=$l_c_peika[$i]/$pei[$i];
  for $j (0..19){
  if($pei2[$i][$j]>0){
  $Plka2[$i][$j]=$lpeika2[$i][$j]/$pei2[$i][$j];
  $Plka_b2[$i][$j]=$l_b_peika2[$i][$j]/$pei2[$i][$j];
  $Plka_c2[$i][$j]=$l_c_peika2[$i][$j]/$pei2[$i][$j];}}
 }
 if($hsumtot>0){
 $Phka[$i]=$hpeika[$i]/$pei[$i];
 $Phka_b[$i]=$h_b_peika[$i]/$pei[$i];
 $Phka_c[$i]=$h_c_peika[$i]/$pei[$i];
  for $j (0..19){	
  if($pei2[$i][$j]>0){
  $Phka2[$i][$j]=$hpeika2[$i][$j]/$pei2[$i][$j];
  $Phka_b2[$i][$j]=$h_b_peika2[$i][$j]/$pei2[$i][$j];
  $Phka_c2[$i][$j]=$h_c_peika2[$i][$j]/$pei2[$i][$j];}}
 }
 if($bsumtot>0){
 $Pbka[$i]=$bpeika[$i]/$pei[$i];
 $Pbka_b[$i]=$b_b_peika[$i]/$pei[$i];
 $Pbka_c[$i]=$b_c_peika[$i]/$pei[$i];
  for $j (0..19){
  if($pei2[$i][$j]>0){
  $Pbka2[$i][$j]=$bpeika2[$i][$j]/$pei2[$i][$j];
  $Pbka_b2[$i][$j]=$b_b_peika2[$i][$j]/$pei2[$i][$j];
  $Pbka_c2[$i][$j]=$b_c_peika2[$i][$j]/$pei2[$i][$j];}} 
 }
printf GOUT ("%s\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",$aatab[$i],$Plka[$i],$Plka_c[$i],$Plka_b[$i],$Phka[$i],$Phka_c[$i],$Phka_b[$i],$Pbka[$i],$Pbka_c[$i],$Pbka_b[$i]);}
}

for $i (0..19){
print GOUT "chou-fasman like parms for given AA on n+1 position\n";
print GOUT "n+1AA\tnAA\tP hincz\tcore\tborder\tP helix\tcore\tborder\tP beta\tcore\tborder\n";
 for $j (0..19){#obliczanie ch-f dla p warunkowych 
# print GOUT "Who brejks all together assuming $aatab[$j] on n+1 position?\n";
# print GOUT "AA\tP hincz\tcore\tborder\tP helix\tcore\tborder\tP beta\tcore\tborder\n";
 if($pei2[$j]>0){
printf GOUT ("%s\t%s\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",$aatab[$j],$aatab[$i],$Plka2[$i][$j],$Plka_c2[$i][$j],$Plka_b2[$i][$j],$Phka2[$i][$j],$Phka_c2[$i][$j],$Phka_b2[$i][$j],$Pbka2[$i][$j],$Pbka_c2[$i][$j],$Pbka_b2[$i][$j]);}
 }
}

for $i (0..19){#sumowanie lewej z prawa
 for $j (0..19){
$plsum[$i][$j]=$pei2[$i][$j]+$pei2[$j][$i];
$plsum_l[$i][$j]=$lpeika2[$i][$j]+$lpeika2[$j][$i];
$plsum_l_c[$i][$j]=$l_c_peika2[$i][$j]+$l_c_peika2[$j][$i];
$plsum_l_b[$i][$j]=$l_b_peika2[$i][$j]+$l_b_peika2[$j][$i];
$plsum_h[$i][$j]=$hpeika2[$i][$j]+$hpeika2[$j][$i];
$plsum_h_c[$i][$j]=$h_c_peika2[$i][$j]+$h_c_peika2[$j][$i];
$plsum_h_b[$i][$j]=$h_b_peika2[$i][$j]+$h_b_peika2[$j][$i];
$plsum_b[$i][$j]=$bpeika2[$i][$j]+$bpeika2[$j][$i];
$plsum_b_c[$i][$j]=$b_c_peika2[$i][$j]+$b_c_peika2[$j][$i];
$plsum_b_b[$i][$j]=$b_b_peika2[$i][$j]+$b_b_peika2[$j][$i];
 }
}
for $i (0..19){#obliczanie wspolczynnikow dla lewej z prawa
 if($lsumtot>0){
  for $j (0..19){
  if($plsum[$i][$j]>0){
  $Plka2[$i][$j]=$plsum_l[$i][$j]/$plsum[$i][$j];
  $Plka_b2[$i][$j]=$plsum_l_b[$i][$j]/$plsum[$i][$j];
  $Plka_c2[$i][$j]=$plsum_l_c[$i][$j]/$plsum[$i][$j];}}
 }
 if($hsumtot>0){
  for $j (0..19){	
  if($plsum[$i][$j]>0){
  $Phka2[$i][$j]=$plsum_h[$i][$j]/$plsum[$i][$j];
  $Phka_b2[$i][$j]=$plsum_h_b[$i][$j]/$plsum[$i][$j];
  $Phka_c2[$i][$j]=$plsum_h_c[$i][$j]/$plsum[$i][$j];}}
 }
 if($bsumtot>0){
  for $j (0..19){
  if($plsum[$i][$j]>0){
  $Pbka2[$i][$j]=$plsum_b[$i][$j]/$plsum[$i][$j];
  $Pbka_b2[$i][$j]=$plsum_b_b[$i][$j]/$plsum[$i][$j];
  $Pbka_c2[$i][$j]=$plsum_b_c[$i][$j]/$plsum[$i][$j];}} 
 }
}
for $i (0..19){
print GOUT "chou-fasman like parms for given AA on n-1 and n+1 position\n";
print GOUT "adj AA\tnAA\tP hincz\tcore\tborder\tP helix\tcore\tborder\tP beta\tcore\tborder\n";
 for $j (0..19){#obliczanie ch-f dla p warunkowych 
 if($pei2[$j]>0){
printf GOUT ("%s\t%s\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",$aatab[$j],$aatab[$i],$Plka2[$i][$j],$Plka_c2[$i][$j],$Plka_b2[$i][$j],$Phka2[$i][$j],$Phka_c2[$i][$j],$Phka_b2[$i][$j],$Pbka2[$i][$j],$Pbka_c2[$i][$j],$Pbka_b2[$i][$j]);}
 }
}

}
1;
