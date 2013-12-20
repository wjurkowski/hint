package MotifWise;
use strict;
use warnings;

require Exporter;
use vars qw(@ISA @EXPORT $VERSION);
our @ISA = qw(Exporter);
our @EXPORT = qw(digmotifs takeloops takehelix takebeta);
$VERSION=1.0;

sub digmotifs {
my ($lsum,$bsum,$hsum)=_@;
my ($j);
my $nsum=0;
my $l_c_sum=0;
my $l_b_sum=0;
my $h_c_sum=0;
my $h_b_sum=0;
my $b_c_sum=0;
my $b_b_sum=0;
my $marg=2;
for $i(0..19){
$l_c_ressum[$i]=0;
$l_b_ressum[$i]=0;
$h_c_ressum[$i]=0;
$h_b_ressum[$i]=0;
$b_c_ressum[$i]=0;
$b_b_ressum[$i]=0;
 for $j(0..19){	
 $l_c_r2sum[$i][$j]=0;
 $l_b_r2sum[$i][$j]=0;
 $h_c_r2sum[$i][$j]=0;
 $h_b_r2sum[$i][$j]=0;
 $b_c_r2sum[$i][$j]=0;
 $b_b_r2sum[$i][$j]=0;}
}

	if($params{'print_motifs'}==1){ 
	print HINCHOUT "PDB code\tmotif\tresid\tAA\tlnR\tV\tD\tD aver\tphi_en\tpsi_en\n";
	print HELIXOUT "PDB code\tmotif\tresid\tAA\tlnR\tV\tD\tD aver\tphi_en\tpsi_en\n";
	print BETAOUT "PDB code\tmotif\tresid\tAA\tlnR\tV\tD\tD aver\tphi_en\tpsi_en\n";}
	if($params{'AA_stat'}==1){
	for $i (0..19){
	$laastat[$i]=0;
	$haastat[$i]=0;
	$baastat[$i]=0;
	$lressum[$i]=0;
	$hressum[$i]=0;
	$bressum[$i]=0;
	$lsum=0;
	$hsum=0;
	$bsum=0;
	 for $j (0..19){
	 $lr2sum[$i][$j]=0;
	 $hr2sum[$i][$j]=0;
	 $br2sum[$i][$j]=0;}}
	print HINCHSTAT "PDB code\tmotif\tAA\tcount\tlenght\tfreq\n";
	print HELIXSTAT "PDB code\tmotif\tAA\tcount\tlenght\tfreq\n";
	print BETASTAT "PDB code\tmotif\tAA\tcount\tlenght\tfreq\n";}
for $k (0..$mnr){
 if($mlen[$k]>=$params{'min_len'}){
  if($mottyp[$k] eq 'L'){#dla petli
  takeloop (\$lsum);	
  }#end of loop block	
  elsif($mottyp[$k] eq 'H'){
  takehelix (\$hsum);	
  }#end of helix block
  elsif($mottyp[$k] eq 'B'){
  takebeta (\$bsum);	
  }#end of beta block
 }#motif lenght limit
}#end of motif iteration

}

sub takeloop{
my ($lsum)=_@;
for $i (0..19){#zerowanie zliczen w obrebie danego motywu
	$laacnt[$i]=0;
	 for $j (0..19){
	 $laa2cnt[$i][$j]=0;}}
	 if($params{'printpdb'}==1){#zapis PDB poczatek i reszta oskrzydlajaca
	 $ind=rindex($pdbname,".");
	 $outpdb=substr($pdbname,0,$ind)."_loop_".$k.".pdb";
	 open(PDBOUT, ">$wyniki/$loops/$outpdb") or die "Can’t open output file: $!, $outpdb";
	 printf PDBOUT "REMARK File generated with perl script: hinge_alchemy.pl (W.Jurkowski) \n";
	 printf PDBOUT "REMARK User: %s Time: %s\n", $me, $czas;
	 	if($k>0){$dl=$mlen[$k-1];
  	 	for $i(1..$motif[$k-1][$dl][9]){
  	 	printf PDBOUT "$motpdb[$k-1][$dl][$i]\n";}}
 	 }
	 for $z (1..$mlen[$k]){
	  if($z >2 and $z <$mlen[$k]-1){#counting core and border residues for Chou-Fasman
		for $i (0..19){
	 	  if($aatab[$i] eq $motif[$k][$z][2]){#AA counted in loop cores
		  $l_c_ressum[$i]++;
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k][$z+1][2]){#AA counted in loop cores
		 	 $l_c_r2sum[$i][$j]++;	
			 last;}	
			}
		  }	
		}	
	  $l_c_sum++;}#overall core count
	  if($params{'print_motifs'}==1){ #prosty wydruk
	  print HINCHOUT "$motif[$k][$z][0]\t$k\t$motif[$k][$z][1]\t$motif[$k][$z][2]\t$motif[$k][$z][3]\t$motif[$k][$z][4]\t$motif[$k][$z][5]\t$motif[$k][$z][6]\t$motif[$k][$z][7]\t$motif[$k][$z][8]\n";}
	  if($params{'AA_stat'}==1){
	 	for $i (0..19){#zliczanie AA
	 	 if($aatab[$i] eq $motif[$k][$z][2]){
		 $laacnt[$i]++;}
			if($z<$mlen[$k]){
			for $j (0..19){#in pairs
	 	  	 if($aatab[$j] eq $motif[$k][$z+1][2]){
		  	 $laa2cnt[$i][$j]++;
			 last;}
			}}
		}
	  }
	  if($params{'printpdb'}==1){#zapis PDB
	   for $i(1..$motif[$k][$z][9]){
	   printf PDBOUT "$motpdb[$k][$z][$i]\n";}
	  }
	 }
#motifs borders
	 for $z (1..$marg){
		for $i (0..19){#begining of k
	 	  if($aatab[$i] eq $motif[$k][$z][2]){
		  $l_b_ressum[$i]++;
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k][$z+1][2]){
		 	 $l_b_r2sum[$i][$j]++;	
			 last;}	
			}	
		  }
		  if($params{'expband'}==1){#begining of k+1
		  if($k<$mnr){
		  if($mlen[$k+1]>=3){
		  if($aatab[$i] eq $motif[$k+1][$z][2]){
		  $l_b_ressum[$i]++;
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k+1][$z+1][2]){
		 	 $l_b_r2sum[$i][$j]++;	
			 last;}	
			}	
		  }}}}
		}
	 $l_b_sum++;}		
	 for $z ($mlen[$k]-$marg+1..$mlen[$k]-1){#ending of k
		for $i (0..19){
	 	  if($aatab[$i] eq $motif[$k][$z][2]){
		  $l_b_ressum[$i]++;
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k][$z+1][2]){
		 	 $l_b_r2sum[$i][$j]++;	
			 last;}	
			}	
		  }
		}
	 $l_b_sum++;}
	 if($params{'expband'}==1){#ending of k-1
	 if($k>0){
	 if($mlen[$k-1]>=2){
	 for $z ($mlen[$k-1]-$marg+1..$mlen[$k-1]-1){
		for $i (0..19){
		  if($aatab[$i] eq $motif[$k-1][$z][2]){
		  $l_b_ressum[$i]++;
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k-1][$z+1][2]){
		 	 $l_b_r2sum[$i][$j]++;}	
			}	
		  }
		}
	 $l_b_sum++;}}}}
	 	for $i (0..19){#for terminal residue of k
	 	  if($aatab[$i] eq $motif[$k][$mlen[$k]][2]){
		  $l_b_ressum[$i]++;
			if($k<$mnr){
			if($mlen[$k+1]>=2){
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k+1][1][2]){
		 	 $l_b_r2sum[$i][$j]++;	
			 last;}	
			}}}	
		  }
		  if($params{'expband'}==1){#for terminal residue of k-1
		  if($k>0){
		  if($mlen[$k-1]>=2){
		  if($aatab[$i] eq $motif[$k-1][$mlen[$k-1]][2]){
		  $l_b_ressum[$i]++;
			if($k<$mnr){
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k][1][2]){
		 	 $l_b_r2sum[$i][$j]++;	
			 last;}	
			}}	
		  }}}}
		}	

	 if($params{'printpdb'}==1){#zapis PDB reszta oskrzydlajaca i zakonczenie
	  if($k lt $mnr){for $i(1..$motif[$k+1][1][9]){
	  printf PDBOUT "$motpdb[$k+1][1][$i]\n";}}
	 printf PDBOUT "TER\nEND\n";
	 close (PDBOUT);}
	 if($params{'AA_stat'}==1){#warunki formerow i brejkerow
	 $lhform=$laacnt[0]+$laacnt[5]+$laacnt[6]+$laacnt[8]+$laacnt[9]+$laacnt[10]+$laacnt[11]+$laacnt[12]+$laacnt[13]+$laacnt[17]+$laacnt[19];
	 $lhbrej=$laacnt[2]+$laacnt[7]+$laacnt[14]+$laacnt[18];
	 $lbform=$laacnt[4]+$laacnt[5]+$laacnt[9]+$laacnt[10]+$laacnt[12]+$laacnt[13]+$laacnt[16]+$laacnt[17]+$laacnt[18]+$laacnt[19];
	 $lbbrej=$laacnt[2]+$laacnt[6]+$laacnt[8]+$laacnt[11]+$laacnt[14]+$laacnt[15];
	 $lreszt=$laacnt[1]+$laacnt[3];
	  for $i (0..19){
	  $lressum[$i]=$lressum[$i]+$laacnt[$i];
	  $laastat[$i]=$laacnt[$i]/$mlen[$k];
	  printf HINCHSTAT ("%s\t%4d\t%s\t%5.3f\t%5.3f\t%5.3f\n",$motif[$k][1][0],$k,$aatab[$i],$laacnt[$i],$mlen[$k],$laastat[$i]);
		for $j (0..19){
		$lr2sum[$i][$j]=$lr2sum[$i][$j]+$laa2cnt[$i][$j];
		}
	  }
	  $lhfstat=$lhform/$mlen[$k];
	  $lhbstat=$lhbrej/$mlen[$k];
	  $lbfstat=$lbform/$mlen[$k];
	  $lbbstat=$lbbrej/$mlen[$k];
	  $lrestat=$lreszt/$mlen[$k];
	  printf HINCHSTAT ("%s\t%4d\t%s\t%4d\t%4d\t%5.3f\n",$motif[$k][1][0],$k,HF,$lhform,$mlen[$k],$lhfstat);
	  printf HINCHSTAT ("%s\t%4d\t%s\t%4d\t%4d\t%5.3f\n",$motif[$k][1][0],$k,HB,$lhbrej,$mlen[$k],$lhbstat);
	  printf HINCHSTAT ("%s\t%4d\t%s\t%4d\t%4d\t%5.3f\n",$motif[$k][1][0],$k,BF,$lbform,$mlen[$k],$lbfstat);
	  printf HINCHSTAT ("%s\t%4d\t%s\t%4d\t%4d\t%5.3f\n",$motif[$k][1][0],$k,BB,$lbbrej,$mlen[$k],$lbbstat);
	  printf HINCHSTAT ("%s\t%4d\t%s\t%4d\t%4d\t%5.3f\n",$motif[$k][1][0],$k,RE,$lreszt,$mlen[$k],$lrestat);}
	 $lsum=$lsum+$mlen[$k]; #sumaryczna dlugosc motywow

}

sub takehelix {
my ($hsum)=_@;
for $i (0..19){
	$haacnt[$i]=0;
	 for $j (0..19){
	 $haa2cnt[$i][$j]=0;}}
	 if($params{'printpdb'}==1){#zapis PDB
	 $ind=rindex($pdbname,".");
	 $outpdb=substr($pdbname,0,$ind)."_helix_".$k.".pdb";
	 open(PDBOUT, ">$wyniki/$helices/$outpdb") or die "Can’t open output file: $!, $outpdb";
	 printf PDBOUT "REMARK File generated with perl script: hinge_alchemy.pl (W.Jurkowski) \n";
	 printf PDBOUT "REMARK User: %s Time: %s\n", $me, $czas;
	 	if($k>0){ $dl=$mlen[$k-1];
  	 	for $i(1..$motif[$k-1][$dl][9]){
  	 	printf PDBOUT "$motpdb[$k-1][$dl][$i]\n";}}
 	 }
	 for $z (1..$mlen[$k]){
	  if($z >2 and $z <$mlen[$k]-1){#counting core and border residues for Chou-Fasman
		for $i (0..19){
	 	  if($aatab[$i] eq $motif[$k][$z][2]){#AA counted in helix cores
		  $h_c_ressum[$i]++;
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k][$z+1][2]){#AA counted in helix cores
		 	 $h_c_r2sum[$i][$j]++;	
			 last;}	
			}
		  }	
		}	
	  $h_c_sum++;}#overall core count

	  if($params{'print_motifs'}==1){ #prosty wydruk
	  print HELIXOUT "$motif[$k][$z][0]\t$k\t$motif[$k][$z][1]\t$motif[$k][$z][2]\t$motif[$k][$z][3]\t$motif[$k][$z][4]\t$motif[$k][$z][5]\t$motif[$k][$z][6]\t$motif[$k][$z][7]\t$motif[$k][$z][8]\n";}
	  if($params{'AA_stat'}==1){
	 	 for $i (0..19){#zliczanie AA
	 	  if($aatab[$i] eq $motif[$k][$z][2]){
		  $haacnt[$i]++;}
			if($z<$mlen[$k]){
			for $j (0..19){#in pairs
	 	  	 if($aatab[$j] eq $motif[$k][$z+1][2]){
		  	 $haa2cnt[$i][$j]++;
			 last;}
			}}	
		 }
	  }
	  if($params{'printpdb'}==1){#zapis PDB
	  	for $i(1..$motif[$k][$z][9]){
	  	printf PDBOUT "$motpdb[$k][$z][$i]\n";}
	  }
	 }

	 #motifs borders
	 for $z (1..$marg){
		for $i (0..19){#begining of k
	 	  if($aatab[$i] eq $motif[$k][$z][2]){
		  $h_b_ressum[$i]++;
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k][$z+1][2]){
		 	 $h_b_r2sum[$i][$j]++;	
			 last;}	
			}	
		  }
		  if($params{'expband'}==1){#begining of k+1
		  if($k<$mnr){
		  if($mlen[$k+1]>=3){
		  if($aatab[$i] eq $motif[$k+1][$z][2]){
		  $h_b_ressum[$i]++;
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k+1][$z+1][2]){
		 	 $h_b_r2sum[$i][$j]++;	
			 last;}	
			}	
		  }}}}
		}
	 $h_b_sum++;}		
	 for $z ($mlen[$k]-$marg+1..$mlen[$k]-1){#ending of k
		for $i (0..19){
	 	  if($aatab[$i] eq $motif[$k][$z][2]){
		  $h_b_ressum[$i]++;
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k][$z+1][2]){
		 	 $h_b_r2sum[$i][$j]++;	
			 last;}	
			}	
		  }
		}
	 $h_b_sum++;}
	 if($params{'expband'}==1){#ending of k-1
	 if($k>0){
	 if($mlen[$k-1]>=2){
	 for $z ($mlen[$k-1]-$marg+1..$mlen[$k-1]-1){
		for $i (0..19){
		  if($aatab[$i] eq $motif[$k-1][$z][2]){#AA counted in loop borders
		  $h_b_ressum[$i]++;
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k-1][$z+1][2]){
		 	 $h_b_r2sum[$i][$j]++;	
			 last;}	
			}	
		  }
		}
	 $h_b_sum++;}}}}
	 	for $i (0..19){#for terminal residue of k
	 	  if($aatab[$i] eq $motif[$k][$mlen[$k]][2]){#AA counted in loop borders
		  $h_b_ressum[$i]++;
			if($mlen[$k+1]>=2){
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k+1][1][2]){
		 	 $h_b_r2sum[$i][$j]++;	
			 last;}	
			}}	
		  }
		  if($params{'expband'}==1){#for terminal residue of k-1
		  if($k>0){
		  if($mlen[$k-1]>=2){
		  if($aatab[$i] eq $motif[$k-1][$mlen[$k-1]][2]){#AA counted in loop borders
		  $h_b_ressum[$i]++;
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k][1][2]){
		 	 $h_b_r2sum[$i][$j]++;	
			 last;}	
			}	
		  }}}}
		}	

	 if($params{'printpdb'}==1){#zapis PDB reszta oskrzydlajaca i zakonczenie
	  if($k lt $mnr){for $i(1..$motif[$k+1][1][9]){
		printf PDBOUT "$motpdb[$k+1][1][$i]\n";}}
	 printf PDBOUT "TER\nEND\n";
	 close (PDBOUT);}
	 if($params{'AA_stat'}==1){#warunki formerow i brejkerow
	 $hhform=$haacnt[0]+$haacnt[5]+$haacnt[6]+$haacnt[8]+$haacnt[9]+$haacnt[10]+$haacnt[11]+$haacnt[12]+$haacnt[13]+$haacnt[17]+$haacnt[19];
	 $hhbrej=$haacnt[2]+$haacnt[7]+$haacnt[14]+$haacnt[18];
	 $hbform=$haacnt[4]+$haacnt[5]+$haacnt[9]+$haacnt[10]+$haacnt[12]+$haacnt[13]+$haacnt[16]+$haacnt[17]+$haacnt[18]+$haacnt[19];
	 $hbbrej=$haacnt[2]+$haacnt[6]+$haacnt[8]+$haacnt[11]+$haacnt[14]+$haacnt[15];
	 $hreszt=$haacnt[1]+$haacnt[3];
	  for $i (0..19){
	  $hressum[$i]=$hressum[$i]+$haacnt[$i];
	  $haastat[$i]=$haacnt[$i]/$mlen[$k];
	  printf HELIXSTAT ("%s\t%4d\t%s\t%5.3f\t%5.3f\t%5.3f\n",$motif[$k][1][0],$k,$aatab[$i],$haacnt[$i],$mlen[$k],$haastat[$i]);
		for $j (0..19){
		$hr2sum[$i][$j]=$hr2sum[$i][$j]+$haa2cnt[$i][$j];
		}
	  }
	 $hhfstat=$hhform/$mlen[$k];
	 $hhbstat=$hhbrej/$mlen[$k];
	 $hbfstat=$hbform/$mlen[$k];
	 $hbbstat=$hbbrej/$mlen[$k];
	 $hrestat=$hreszt/$mlen[$k];
	 printf HELIXSTAT ("%s\t%4d\t%s\t%4d\t%4d\t%5.3f\n",$motif[$k][1][0],$k,HF,$hhform,$mlen[$k],$hhfstat);
	 printf HELIXSTAT ("%s\t%4d\t%s\t%4d\t%4d\t%5.3f\n",$motif[$k][1][0],$k,HB,$hhbrej,$mlen[$k],$hhbstat);
	 printf HELIXSTAT ("%s\t%4d\t%s\t%4d\t%4d\t%5.3f\n",$motif[$k][1][0],$k,BF,$hbform,$mlen[$k],$hbfstat);
	 printf HELIXSTAT ("%s\t%4d\t%s\t%4d\t%4d\t%5.3f\n",$motif[$k][1][0],$k,BB,$hbbrej,$mlen[$k],$hbbstat);
	 printf HELIXSTAT ("%s\t%4d\t%s\t%4d\t%4d\t%5.3f\n",$motif[$k][1][0],$k,RE,$hreszt,$mlen[$k],$hrestat);}	
	 $hsum=$hsum+$mlen[$k];#sumaryczna dlugosc motywow
}

sub takebeta{
my ($bsum)=_@;
for $i (0..19){
	$baacnt[$i]=0;
	 for $j (0..19){
	 $baa2cnt[$i][$j]=0;}}
	 if($params{'printpdb'}==1){#zapis PDB
	 $ind=rindex($pdbname,".");
	 $outpdb=substr($pdbname,0,$ind)."_beta_".$k.".pdb";
	 open(PDBOUT, ">$wyniki/$betas/$outpdb") or die "Can’t open output file: $!, $outpdb";
	 printf PDBOUT "REMARK File generated with perl script: hinge_alchemy.pl (W.Jurkowski) \n";
	 printf PDBOUT "REMARK User: %s Time: %s\n", $me, $czas;
	 	if($k>0){$dl=$mlen[$k-1];
  	 	for $i(1..$motif[$k-1][$dl][9]){
  	 	printf PDBOUT "$motpdb[$k-1][$dl][$i]\n";}}
	 }
	 for $z (1..$mlen[$k]){
	  if($z >2 and $z <$mlen[$k]-1){#counting core and border residues for Chou-Fasman
		for $i (0..19){
	 	  if($aatab[$i] eq $motif[$k][$z][2]){#AA counted in loop cores
		  $b_c_ressum[$i]++;
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k][$z+1][2]){#AA counted in loop cores
		 	 $b_c_r2sum[$i][$j]++;	
			 last;}	
			}
		  }	
		}	
	  $b_c_sum++;}#overall core count
		
	  if($params{'print_motifs'}==1){ #prosty wydruk
	  print BETAOUT "$motif[$k][$z][0]\t$k\t$motif[$k][$z][1]\t$motif[$k][$z][2]\t$motif[$k][$z][3]\t$motif[$k][$z][4]\t$motif[$k][$z][5]\t$motif[$k][$z][6]\t$motif[$k][$z][7]\t$motif[$k][$z][8]\n";}
	  if($params{'AA_stat'}==1){
	 	for $i (0..19){#zliczanie AA
	 	 if($aatab[$i] eq $motif[$k][$z][2]){
		 $baacnt[$i]++;}
			if($z<$mlen[$k]){
			for $j (0..19){#in pairs
	 	  	 if($aatab[$j] eq $motif[$k][$z+1][2]){
		  	 $baa2cnt[$i][$j]++;
			 last;}
			}}	
		}
	  }
	  if($params{'printpdb'}==1){#zapis PDB
		for $i(1..$motif[$k][$z][9]){
		printf PDBOUT "$motpdb[$k][$z][$i]\n";}
	  }
	 }

	 #motifs borders count
	 for $z (1..$marg){
		for $i (0..19){#begining of k
	 	  if($aatab[$i] eq $motif[$k][$z][2]){
		  $b_b_ressum[$i]++;
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k][$z+1][2]){
		 	 $b_b_r2sum[$i][$j]++;	
			 last;}	
			}	
		  }
		  if($params{'expband'}==1){#begining of k+1
		  if($k<$mnr){
		  if($mlen[$k+1]>=3){
		  if($aatab[$i] eq $motif[$k+1][$z][2]){
		  $b_b_ressum[$i]++;
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k+1][$z+1][2]){
		 	 $b_b_r2sum[$i][$j]++;	
			 last;}	
			}	
		  }}}}
		}
	 $b_b_sum++;}		
	 for $z ($mlen[$k]-$marg+1..$mlen[$k]-1){#ending of k
		for $i (0..19){
	 	  if($aatab[$i] eq $motif[$k][$z][2]){
		  $b_b_ressum[$i]++;
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k][$z+1][2]){
		 	 $b_b_r2sum[$i][$j]++;	
			 last;}	
			}	
		  }
		}
	 $b_b_sum++;}
	 if($params{'expband'}==1){#ending of k-1
	 if($k>0){
	 if($mlen[$k-1]>=2){ 
	 for $z ($mlen[$k-1]-$marg+1..$mlen[$k-1]-1){
		for $i (0..19){
		  if($aatab[$i] eq $motif[$k-1][$z][2]){#AA counted in loop borders
		  $b_b_ressum[$i]++;
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k-1][$z+1][2]){
		 	 $b_b_r2sum[$i][$j]++;	
			 last;}	
			}	
		  }
		}
	 $b_b_sum++;}}}}
	 	for $i (0..19){#for terminal residue of k
	 	  if($aatab[$i] eq $motif[$k][$mlen[$k]][2]){#AA counted in loop borders
		  $b_b_ressum[$i]++;
			if($mlen[$k+1]>=2){
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k+1][1][2]){
		 	 $b_b_r2sum[$i][$j]++;	
			 last;}	
			}}	
		  }
		  if($params{'expband'}==1){#for terminal residue of k-1
		  if($k>0){
		  if($mlen[$k-1]>=2){
		  if($aatab[$i] eq $motif[$k-1][$mlen[$k-1]][2]){#AA counted in loop borders
		  $b_b_ressum[$i]++;
			for $j (0..19){
	 	  	 if($aatab[$j] eq $motif[$k][1][2]){
		 	 $b_b_r2sum[$i][$j]++;	
			 last;}	
			}	
		  }}}}
		}	

	 if($params{'printpdb'}==1){#zapis PDB reszta oskrzydlajaca i zakonczenie
	  if($k lt $mnr){for $i(1..$motif[$k+1][1][9]){
	  printf PDBOUT "$motpdb[$k+1][1][$i]\n";}}
	 printf PDBOUT "TER\nEND\n";
	 close (PDBOUT);}
	 if($params{'AA_stat'}==1){
	 $bhform=$baacnt[0]+$baacnt[5]+$baacnt[6]+$baacnt[8]+$baacnt[9]+$baacnt[10]+$baacnt[11]+$baacnt[12]+$baacnt[13]+$baacnt[17]+$baacnt[19];
	 $bhbrej=$baacnt[2]+$baacnt[7]+$baacnt[14]+$baacnt[18];
	 $bbform=$baacnt[4]+$baacnt[5]+$baacnt[9]+$baacnt[10]+$baacnt[12]+$baacnt[13]+$baacnt[16]+$baacnt[17]+$baacnt[18]+$baacnt[19];
	 $bbbrej=$baacnt[2]+$baacnt[6]+$baacnt[8]+$baacnt[11]+$baacnt[14]+$baacnt[15];
	 $breszt=$baacnt[1]+$baacnt[3];
	  for $i (0..19){
	  $bressum[$i]=$bressum[$i]+$baacnt[$i];
	  $baastat[$i]=$baacnt[$i]/$mlen[$k];
	  printf BETASTAT ("%s\t%4d\t%s\t%4d\t%4d\t%5.3f\n",$motif[$k][1][0],$k,$aatab[$i],$baacnt[$i],$mlen[$k],$baastat[$i]);
		for $j (0..19){
		$br2sum[$i][$j]=$br2sum[$i][$j]+$baa2cnt[$i][$j];
		}		
	  }
	$bhfstat=$bhform/$mlen[$k];
	$bhbstat=$bhbrej/$mlen[$k];
	$bbfstat=$bbform/$mlen[$k];
	$bbbstat=$bbbrej/$mlen[$k];
	$brestat=$breszt/$mlen[$k];
	printf BETASTAT ("%s\t%4d\t%s\t%4d\t%4d\t%5.3f\n",$motif[$k][1][0],$k,HF,$bhform,$mlen[$k],$bhfstat);
	printf BETASTAT ("%s\t%4d\t%s\t%4d\t%4d\t%5.3f\n",$motif[$k][1][0],$k,HB,$bhbrej,$mlen[$k],$bhbstat);
	printf BETASTAT ("%s\t%4d\t%s\t%4d\t%4d\t%5.3f\n",$motif[$k][1][0],$k,BF,$bbform,$mlen[$k],$bbfstat);
	printf BETASTAT ("%s\t%4d\t%s\t%4d\t%4d\t%5.3f\n",$motif[$k][1][0],$k,BB,$bbbrej,$mlen[$k],$bbbstat);
	printf BETASTAT ("%s\t%4d\t%s\t%4d\t%4d\t%5.3f\n",$motif[$k][1][0],$k,RE,$breszt,$mlen[$k],$brestat);}
	$bsum=$bsum+$mlen[$k];#sumaryczna dlugosc motywow

}
1;

