#!/usr/bin/perl -w
#parametry programu:[lista plikow vr] [lista plikow dssp] [nazwa pliku wynikowego]

if ($#ARGV != 0) {die "Program uzywany z parametrami! [plik z parametrami] \n";}
my $rootdir=`pwd`;
chomp $rootdir;

use warnings;
use strict;
use lib '.';
use HingeAlchemy::ReadStrDescr;
use HingeAlchemy::ConstHAlch;
use HingeAlchemy::LHBclass;
use HingeAlchemy::LHBpred;
use HingeAlchemy::MineSeq;
use HingeAlchemy::MineStr;
use HingeAlchemy::MotifWise;
use HingeAlchemy::WriteStr;

my ($defwyn, $lin);
#wczytanie opcji 
open(OPCJE, "< $ARGV[0]") or die "Can not open an input file: $!";
my @param=<OPCJE>;
close (OPCJE);
chomp @param;
my %params = ();
my @para=();
foreach $lin(@param){
my @para=split(/\s+/,$lin);
$params{$para[0]}=$para[1];
}
my $czas=localtime(time());
my $me=getlogin();	

my $roll=int(rand 10000000) +1;
my $wyniki="run".$roll;
#directories with input and output data
my $kat_dssp=$params{'dssp_dir'};
my $kat_pdb=$params{'inppdb_dir'};
if($params{'run_mode'} ==1){
my $kat_vr=$params{'vr_dir'};
$defwyn="$kat_vr-"."$kat_dssp";
mkdir("$wyniki/$defwyn", 0755) if (! -d "$wyniki/$defwyn");
my $wyn_dssp_vr='vr_dssp';
mkdir("$wyniki/$defwyn/$wyn_dssp_vr", 0755) if (! -d "$wyniki/$defwyn/$wyn_dssp_vr");
}
elsif($params{'run_mode'} == 2){
my $kat_dist=$params{'dist_dir'};
$defwyn="$kat_dist-"."$kat_dssp";
mkdir("$wyniki/$defwyn", 0755) if (! -d "$wyniki/$defwyn");
}
elsif($params{'run_mode'} == 3){
my $kat_dist=$params{'dist_dir'};
my $kat_vr=$params{'vr_dir'};
$defwyn="$kat_vr-"."$kat_dist-"."$kat_dssp";
mkdir("$wyniki/$defwyn", 0755) if (! -d "$wyniki/$defwyn");
my $wyn_dssp_vr='vr_dssp';
mkdir("$wyniki/$defwyn/$wyn_dssp_vr", 0755) if (! -d "$wyniki/$defwyn/$wyn_dssp_vr");
}


if($params{'lhb_coeff'} ==1){#poszukiwanie wzorcow
my $lhbwyn=$params{'lhbdir'};
mkdir("$wyniki/$lhbwyn", 0755) if (! -d "$wyniki/$lhbwyn");
}
if($params{'printseq'}==1){#drukowanie sekwencji
my $sequences=$params{'seqdir'};
mkdir("$wyniki/$sequences", 0755) if (! -d "$wyniki/$sequences");
}
if($params{'printpdb'}==1){#drukowanie pdb
my $loops=$params{'lpdbdir'};
my $helices=$params{'hpdbdir'};
my $betas=$params{'bpdbdir'};
mkdir("$wyniki/$loops", 0755) if (! -d "$wyniki/$loops");
mkdir("$wyniki/$helices", 0755) if (! -d "$wyniki/$helices");
mkdir("$wyniki/$betas", 0755) if (! -d "$wyniki/$betas");
}

#pliki wynikowe
my $gwyn=$params{'output'};
my $awyn=$params{'averout'};
my $wwyn=$params{'backout'};
my $hwyn=$params{'alfaout'};
my $bwyn=$params{'betaout'};
my $twyn=$params{'turnout'};
my $lwyn=$params{'loopout'};

#ogolny
open(GOUT,"> $wyniki/$defwyn/$gwyn") or die "Can’t write output file: $!";

if($params{'average'} ==1){#wyniki usrednione po strukturach drugorzedowych
open(AOUT,"> $wyniki/$defwyn/$awyn") or die "Can’t write output file: $!";}
if($params{'sec_sort'} ==1){#wyniki VR i dist dla kolejnych posortowanych typow struktury drugorzedowej bez usrednienia
open(WOUT,"> $wyniki/$defwyn/$wwyn") or die "Can’t write output file: $!";
open(HOUT,"> $wyniki/$defwyn/$hwyn") or die "Can’t write output file: $!";
open(BOUT,"> $wyniki/$defwyn/$bwyn") or die "Can’t write output file: $!";
open(TOUT,"> $wyniki/$defwyn/$twyn") or die "Can’t write output file: $!";
open(LOUT,"> $wyniki/$defwyn/$lwyn") or die "Can’t write output file: $!";}
if($params{'motifs'}==1){#wyniki dla analizy motywow
if($params{'print_motifs'}==1){
my $motwyn1=$params{'motifout'};
open(HINCHOUT, "> $wyniki/$defwyn/hinge_$motwyn1") or die "Can’t write output file: $!";
open(HELIXOUT, "> $wyniki/$defwyn/helix_$motwyn1") or die "Can’t write output file: $!";
open(BETAOUT, "> $wyniki/$defwyn/beta_$motwyn1") or die "Can’t write output file: $!";}
if($params{'AA_stat'}==1){
my $motwyn2=$params{'motifstatout'};
open(HINCHSTAT, "> $wyniki/$defwyn/hinge_$motwyn2") or die "Can’t write output file: $!";
open(HELIXSTAT, "> $wyniki/$defwyn/helix_$motwyn2") or die "Can’t write output file: $!";
open(BETASTAT, "> $wyniki/$defwyn/beta_$motwyn2") or die "Can’t write output file: $!";}
}
if($params{'patterns'} ==1){#poszukiwanie wzorcow
my $pattwyn=$params{'patternout'};
open(PATTERNS, ">$wyniki/$defwyn/1-$pattwyn") or die "Can’t open output file: $!";
#open(PATTERNS2, ">$wyniki/$defwyn/2-$pattwyn") or die "Can’t open output file: $!";
}

printf GOUT "Main ouput file generated with perl script: VR_DSSP_filter8.pl (W.Jurkowski) \n";
printf GOUT "User: %s Time: %s\n", $me, $czas;
printf GOUT "Selected parameters were used: \n";
while ( my ($key, $value) = each(%params) ) {
printf GOUT "$key = $value\n";
}
#opis plikow outputowych
printf GOUT "description of output files:\n";
if($params{'run_mode'}==1){
printf GOUT "VROUT *.dssp.vr - dssp vs VR for each chain separately\n";
printf GOUT "[Order #]\t[Residue #]\t[Structure Code]\t[lnR]\t[V]\n";
printf GOUT "backout - dssp vs VR for all chains\n";
printf GOUT "[Order #]\t[Residue #]\t[Structure Code]\t[lnR]\t[V]\n";
printf GOUT "loopout - dssp vs VR for loop residues, all chains\n";
printf GOUT "[Order #]\t[Residue #]\t[Structure Code]\t[lnR]\t[V]\n";
printf GOUT "alfaout - dssp vs VR for helical residues, all chains\n";
printf GOUT "[Order #]\t[Residue #]\t[Structure Code]\t[lnR]\t[V]\n";
printf GOUT "betaout - dssp vs VR for beta-strand residues, all chains\n";
printf GOUT "[Order #]\t[Residue #]\t[Structure Code]\t[lnR]\t[V]\n";
printf GOUT "turnout - dssp vs VR for turn residues, all chains\n";
printf GOUT "[Order #]\t[Residue #]\t[Structure Code]\t[lnR]\t[V]\n";
}
elsif($params{'run_mode'}==2){
printf GOUT "backout - dssp vs VR for all chains\n";
printf GOUT "[Order #]\t[Residue #]\t[Structure Code]\t[D:Ell-Nat distance]\t[D averaged]\t[ellphi-phi]\t[ellpsi-psi]\n";
printf GOUT "loopout - dssp vs VR for loop residues, all chains\n";
printf GOUT "[Order #]\t[Residue #]\t[Structure Code]\t[D:Ell-Nat distance]\t[D averaged]\t[ellphi-phi]\t[ellpsi-psi]\n";
printf GOUT "alfaout - dssp vs VR for helical residues, all chains\n";
printf GOUT "[Order #]\t[Residue #]\t[Structure Code]\t[D:Ell-Nat distance]\t[D averaged]\t[ellphi-phi]\t[ellpsi-psi]\n";
printf GOUT "betaout - dssp vs VR for beta-strand residues, all chains\n";
printf GOUT "[Order #]\t[Residue #]\t[Structure Code]\t[D:Ell-Nat distance]\t[D averaged]\t[ellphi-phi]\t[ellpsi-psi]\n";
printf GOUT "turnout - dssp vs VR for turn residues, all chains\n";
printf GOUT "[Order #]\t[Residue #]\t[Structure Code]\t[D:Ell-Nat distance]\t[D averaged]\t[ellphi-phi]\t[ellpsi-psi]\n";
}
elsif($params{'run_mode'}==3){
printf GOUT "backout - dssp vs VR for all chains\n";
printf GOUT "[Order #]\t[Residue #]\t[Structure Code]\t[lnR]\t[V]\t[D:Ell-Nat distance]\t[D averaged]\t[ellphi-phi]\t[ellpsi-psi]\n";
printf GOUT "loopout - dssp vs VR for loop residues, all chains\n";
printf GOUT "[Order #]\t[Residue #]\t[Structure Code]\t[lnR]\t[V]\t[D:Ell-Nat distance]\t[D averaged]\t[ellphi-phi]\t[ellpsi-psi]\n";
printf GOUT "alfaout - dssp vs VR for helical residues, all chains\n";
printf GOUT "[Order #]\t[Residue #]\t[Structure Code]\t[lnR]\t[V]\t[D:Ell-Nat distance]\t[D averaged]\t[ellphi-phi]\t[ellpsi-psi]\n";
printf GOUT "betaout - dssp vs VR for beta-strand residues, all chains\n";
printf GOUT "[Order #]\t[Residue #]\t[Structure Code]\t[lnR]\t[V]\t[D:Ell-Nat distance]\t[D averaged]\t[ellphi-phi]\t[ellpsi-psi]\n";
printf GOUT "turnout - dssp vs VR for turn residues, all chains\n";
printf GOUT "[Order #]\t[Residue #]\t[Structure Code]\t[lnR]\t[V]\t[D:Ell-Nat distance]\t[D averaged]\t[ellphi-phi]\t[ellpsi-psi]\n";
}

#wczytaj pliki z nazwami plikow
printf GOUT "reading in data...\n";
#wspolny
if($params{'printpdb'}==1){#dla zapisania pdb
my $pdb_inp=$params{'pdb_input'};
open(INPUT4, "< $pdb_inp") or die "Can’t open input file: $!";
my @lipdb=<INPUT4>;
close (INPUT4);
chomp @lipdb;
}
my $dssp_inp=$params{'dssp_input'};
open(INPUT2, "< $dssp_inp") or die "Can’t open input file: $!";
my @lidssp=<INPUT2>;
close (INPUT2);
chomp @lidssp;
#wybrane
if($params{'run_mode'}==1){
my $vr_inp=$params{'vr_input'};
open(INPUT1, "< $vr_inp") or die "Can’t open input file: $!";
my @livr=<INPUT1>;
close (INPUT1);
chomp @livr;
}
elsif($params{'run_mode'}==2){
my $dist_inp=$params{'dist_input'};
open(INPUT3, "< $dist_inp") or die "Can’t open input file: $!";
my @lidist=<INPUT3>;
close (INPUT3);
chomp @lidist;
}
elsif($params{'run_mode'}==3){
my $vr_inp=$params{'vr_input'};
open(INPUT1, "< $vr_inp") or die "Can’t open input file: $!";
my @livr=<INPUT1>;
close (INPUT1);
chomp @livr;
my $dist_inp=$params{'dist_input'};
open(INPUT3, "< $dist_inp") or die "Can’t open input file: $!";
my @lidist=<INPUT3>;
close (INPUT3);
chomp @lidist;
}

#zmienne wielobialkowe

my $lsumtot=0;
my $l_c_sumtot=0;
my $l_b_sumtot=0;
my $hsumtot=0;
my $h_c_sumtot=0;
my $h_b_sumtot=0;
my $bsumtot=0;
my $b_c_sumtot=0;
my $b_b_sumtot=0;
my (@nressumtot,@lressumtot,@l_c_ressumtot,@l_b_ressumtot,@hressumtot,@h_c_ressumtot,@h_b_ressumtot);
my (@bressumtot,@b_c_ressumtot,@b_b_ressumtot,@nr2sumtot,@lr2sumtot,@l_c_r2sumtot,@l_b_r2sumtot);
my (@hr2sumtot,@h_c_r2sumtot,@h_b_r2sumtot,@br2sumtot,@b_c_r2sumtot,@b_b_r2sumtot);
my ($i,$j,$l,$dname);

for $i (0..19){
 $nressumtot[$i]=0;
 $lressumtot[$i]=0;
 $l_c_ressumtot[$i]=0;
 $l_b_ressumtot[$i]=0;
 $hressumtot[$i]=0;
 $h_c_ressumtot[$i]=0;
 $h_b_ressumtot[$i]=0;
 $bressumtot[$i]=0;
 $b_c_ressumtot[$i]=0;
 $b_b_ressumtot[$i]=0;
 for $j (0..19){
 $nr2sumtot[$i][$j]=0;
 $lr2sumtot[$i][$j]=0;
 $l_c_r2sumtot[$i][$j]=0;
 $l_b_r2sumtot[$i][$j]=0;
 $hr2sumtot[$i][$j]=0;
 $h_c_r2sumtot[$i][$j]=0;
 $h_b_r2sumtot[$i][$j]=0;
 $br2sumtot[$i][$j]=0;
 $b_c_r2sumtot[$i][$j]=0;
 $b_b_r2sumtot[$i][$j]=0;
}}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if($params{'runtype'}==1){
foreach $dname (@lidssp){#iterates DSSP files
#reads in structural descriptors 
my ($kat_dssp,$dname);
my $katalog="$rootdir/$kat_dssp";
my $ind=index($dname,".");
my $pdb=substr($dname,0,$ind);
chdir "$katalog" or die "Nie moge wejsc do $kat_dssp: $!\n";
printf GOUT "changing directory to: $katalog\n";
printf GOUT "dssp file analyzed:\t$dname\n";
open(DSSPINP, "< $dname") or die "Can’t open input file: $!";
#open(DSSPOUT, "> ../$wyniki/$wyn_dssp/$dname".".o") or die "Can’t write output file: $!";	
my @dssp=<DSSPINP>;
close (DSSPINP);
#chomp @dssp;
printf GOUT "digesting dssp files...\n";#digest dssp files
dsspdig (\@dssp,\@loop,\@ord,\@struct);

if($params{'run_mode'} == 1 or $params{'run_mode'} == 3){
chdir	"$rootdir/$kat_vr" or die "Nie moge wejsc do $kat_vr: $!\n";	
printf GOUT "changing directory to: ../$kat_vr\n";
printf GOUT "digesting vr files...\n";
#check if all files are present, list of *.vr files and *.dssp files should match
foreach $vrdata (@livr){
 my $pat1=substr($vrdata,0,7);
 $pat1=~tr/A-Z/a-z/;
 my $pat2=substr($dname,0,7);
 $pat2=~tr/A-Z/a-z/;
	if($pat1 eq $pat2){
	$matched=1;
	$vrname=$vrdata;
	last;
	}
 }
  	if ($matched==0){
	printf GOUT "no VR file found matching dssp file: $dname: $!\n";
	last;
	}
#if matched proceedes vr vs. dssp analysis
 open(VRINP, "< $vrname") or (printf GOUT "Can’t open input file: $vrname\n" and next);
 my $ind=index($vrname,".");
 my $fname=substr($vrname,0,$ind).".dssp.vr";
 my @fauer=<VRINP>;
 close (VRINP);
vrdig (\@fauer);#analiza VR
}

if($params{'run_mode'} == 2 or $params{'run_mode'} == 3){
chdir	"$rootdir/$kat_dist" or die "Nie moge wejsc do $kat_dist: $!\n";	
printf GOUT "changing directory to: ../$kat_dist\n";
printf GOUT "digesting fipsi distances...\n";

#check if all files are present, list of *.vr files and *.dssp files should match
my $matched=0;
 foreach $distdata (@lidist){
 my $pat1=substr($distdata,0,7);
 $pat1=~tr/A-Z/a-z/;
 my $pat2=substr($dname,0,7);
 $pat2=~tr/A-Z/a-z/;
	if($pat1 eq $pat2){
	$matched=1;
	$distname=$distdata;
	last;
	}
 }
  	if ($matched==0){
	printf GOUT "no RM Distances file found which matches dssp file: $dname: $!\n";
	last;}

#if matched proceedes dist vs. dssp analysis
 open(DISTINP, "< $distname") or (printf GOUT "Can’t open input file: $distname\n" and next);
 my @rmdist=<DISTINP>;
 close (DISTINP);
 
distdig (\@rmdist);# analiza odleglosci
}
 
if($params{'printpdb'}==1){
pdbread ();#reads input pdb files for later modified pdb print
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#poczatek analiz zbiorczych
printf GOUT "back to parent directory $rootdir/ \n";
chdir	"$rootdir/" or die "Nie moge wejsc do $!\n";	

if($params{'sec_sort'}==1){#testowanie i zliczanie reszt petlowych, nie petlowych, helikalnych,beta, zwrotnych 
 if($params{'run_mode'} == 1){
 open (VROUT, "> ../$defwyn/$wyn_dssp_vr/$fname") or die "Can’t open output file: $!";}
secsort ();	
}#koniec sortowania

if($params{'average'} == 1){#watrosci srednie	
avsecsort();	
} # koniec bloku obliczania srednich i odchylen standardowych

if($params{'patterns'}==1){#analiza wzorcow oddzialywan	
watpatt ();  
}#koniec analizy wzorcow

if($params{'motifs'}==1){#definicja motywow
motdef (); 
}#koniec definicji motywow

if($params{'motifs'}==1){#analiza motywow
digmotifs ();
}# koniec bloku analizy motywow

if($params{'chou_fasman'}==1){#definicje Chou-Fasmana
print GOUT "Who brejks?\n";
print GOUT "AA\tP hincz\tcore\tborder\tP helix\tcore\tborder\tP beta\tcore\tborder\n";
defchoufas ();
printf GOUT ("%s\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",$aatab[$i],$Plka[$i],$Plka_c[$i],$Plka_b[$i],$Phka[$i],$Phka_c[$i],$Phka_b[$i],$Pbka[$i],$Pbka_c[$i],$Pbka_b[$i]);
}#koniec

if($params{'lhb_coeff'}==1){# wyznaczanie lhb
$lhbname="$pdb"."_lhb.dat";
$w_lhbname="$pdb"."_w_lhb.dat";
open(LHB, ">$wyniki/$lhbwyn/$lhbname") or die "Can’t open output file: $! $wyniki/$lhbwyn/$lhbname";
open(WLHB, ">$wyniki/$lhbwyn/$w_lhbname") or die "Can’t open output file: $! $wyniki/$lhbwyn/$w_lhbname";
printf LHB ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n","resid","AA","loop","helix","beta","p loop","p helix","p beta","H","break","lbreak","hbreak","bbreak","lform","hform","bform","wato h","wato b","llhb","hlhb","blhb");
printf WLHB ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n","resid","AA","loop","helix","beta","p loop","p helix","p beta","H","break","lbreak","hbreak","bbreak","lform","hform","bform","wato h","wato b","llhb","hlhb","blhb");
my $motl=$params{'lhb_window'};
lhbdescr (\@ord,\@aatab,\@struct,\@hfob_s1,\@hfob_s2,\@nhfob_s1,\@lnorm,\@hnorm,\@bnorm);
for $j(1..$#ord){# wydruk
printf LHB ("%s\t%s\t%3.1f\t%3.1f\t%3.1f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",$j,$struct[$j],$lstate[$j],$hstate[$j],$bstate[$j],$p_lstate[$j],$p_hstate[$j],$p_bstate[$j],$hfobia1[$j],$brejker[$j],$lbrejker[$j],$hbrejker[$j],$bbrejker[$j],$lformer[$j],$hformer[$j],$bformer[$j],$wathbrr[$j],$watbbrr[$j],$llhb[$j],$hlhb[$j],$blhb[$j]);
printf WLHB ("%s\t%s\t%3.1f\t%3.1f\t%3.1f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",$j,$struct[$j],$lstate[$j],$hstate[$j],$bstate[$j],$p_lstate[$j],$p_hstate[$j],$p_bstate[$j],$w_hfobia1[$j],$w_brejker[$j],$w_lbrejker[$j],$w_hbrejker[$j],$w_bbrejker[$j],$w_lformer[$j],$w_hformer[$j],$w_bformer[$j],$wathbrr[$j],$watbbrr[$j],$llhb[$j],$hlhb[$j],$blhb[$j]);
}
close (LHB);
}#koniec

if($params{'VR_class'}==1){#motifs classifications basing on VR values clusters
my $ns=$params{'VR_nsigm'};
my $minl=$params{'min_len'};
vrclass(\$ns,\$mnr,\@mlen,\@mottyp,\$minl,\@motif,\$lsum,\$bsum,\$hsum,\@l_pcl,\@h_pcl,\@b_pcl,\@l_totccl,\@h_totccl,\@b_totccl);
print GOUT "p of finding given VR class (1-9) for loops, helices and beta in $pdb\n";
print GOUT "1\t2\t3\t4\t5\t6\t7\t8\t9\n";
printf GOUT ("%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",$l_pcl[1],$l_pcl[2],$l_pcl[3],$l_pcl[4],$l_pcl[5],$l_pcl[6],$l_pcl[7],$l_pcl[8],$l_pcl[9]);
printf GOUT ("%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",$h_pcl[1],$h_pcl[2],$h_pcl[3],$h_pcl[4],$h_pcl[5],$h_pcl[6],$h_pcl[7],$h_pcl[8],$h_pcl[9]);
printf GOUT ("%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",$b_pcl[1],$b_pcl[2],$b_pcl[3],$b_pcl[4],$b_pcl[5],$b_pcl[6],$b_pcl[7],$b_pcl[8],$b_pcl[9]);
}
if($params{'printseq'}==1){
$seqname="$wyniki/$sequences/$pdb".".fasta";
$llin=int $#ord/80.0;
seqprint (\$seqname,\$llin);
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

}#koniec iteracji po plikach dssp
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if($params{'chou_fasman'}==1){#definicje Chou-Fasmana
sumchoufas ();
}# koniec chou-fasmana

if($params{'VR_class'}==1){#motifs classifications basing on VR values clusters
sumvrclass (\$lsumtot,\$hsumtot,\$bsumtot,\@l_totccl,\@h_totccl,\@b_totccl,\@l_totpcl,\@h_totpcl,\@b_totpcl);
print GOUT "p of finding given VR class (1-9) for loops, helices and beta\n";
print GOUT "1\t2\t3\t4\t5\t6\t7\t8\t9\n";
printf GOUT ("%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",$l_totpcl[1],$l_totpcl[2],$l_totpcl[3],$l_totpcl[4],$l_totpcl[5],$l_totpcl[6],$l_totpcl[7],$l_totpcl[8],$l_totpcl[9]);
printf GOUT ("%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",$h_totpcl[1],$h_totpcl[2],$h_totpcl[3],$h_totpcl[4],$h_totpcl[5],$h_totpcl[6],$h_totpcl[7],$h_totpcl[8],$h_totpcl[9]);
printf GOUT ("%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",$b_totpcl[1],$b_totpcl[2],$b_totpcl[3],$b_totpcl[4],$b_totpcl[5],$b_totpcl[6],$b_totpcl[7],$b_totpcl[8],$b_totpcl[9]);
}
}#end of runtype=1

if($params{'runtype'}==2){
my ($line,$fname,@struct,@sasekw);

#reads sequence
my $inplist=$params{'inpseq_list'};
open(INPUT1, "< $inplist") or die "Can not open an input file: m2 $inplist $!";
my @lista=<INPUT1>;#wczytanie pliku wejsciowego, nazwy plikow z sekwencjami w tablicy @lista
close (INPUT1);
chomp @lista;

foreach $fname (@lista){#po plikach
my $sekw=$params{'inpseq_dir'};
#$katalog="$sekw";
#chdir "$katalog" or die "Nie moge wejsc do $sekw: $!\n";
open(INPUT2, "< $sekw/$fname") or die "Can’t open input file: m3 $sekw/$fname $!";
printf GOUT "file analyzed: $fname\n";

my @data=<INPUT2>;
close (INPUT2);
chomp @data;
my $nl=-1;
my $cc=0;
LINIA:	foreach $line (@data){
	$nl++;
	 if ($line=~/^>/){
	 $cc++;
	 goto READ}
	 else{
	 next LINIA}
READ:		if($cc==1){
		@struct=split //,$data[$nl+1];}
		if($cc==2){
		@sasekw=split //,$data[$nl+1];} 
	}
my $naa=$#struct;#dlugosc sekwencji AA
my $nsa=$#sasekw;#dlugosc sekwencji SA
#print "@struct\n, $naa";

if($params{'patterns'}==1){#analiza wzorcow oddzialywan	
watpatt (\@struct,\@polar,\@spolar,\@plusk,\@minusk,\@hydrof,\@shydrof,\@ord,\@nn1,\@nn2,\@nn3,\@nn4,\@nn21,\@nn22,\@nn23,\@nn24);
  for $j(1..$$ord){
	print PATTERNS "$pdb\t$j\t$ord[$j]\t$struct[$j]\t$nn1[$j]\t$nn2[$j]\t$nn3[$j]\t$nn4[$j]\t$nn21[$j]\t$nn22[$j]\t$nn23[$j]\t$nn24[$j]\n";
  }
}#koniec analizy wzorcow

if($params{'lhb_coeff'}==1){# wyznaczanie lhb

$lhbname="$pdb"."_lhb.dat";
$w_lhbname="$pdb"."_w_lhb.dat";
open(LHB, ">$wyniki/$lhbwyn/$lhbname") or die "Can’t open output file: $! $wyniki/$lhbwyn/$lhbname";
open(WLHB, ">$wyniki/$lhbwyn/$w_lhbname") or die "Can’t open output file: $! $wyniki/$lhbwyn/$w_lhbname";
printf LHB ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n","resid","AA","loop","helix","beta","p loop","p helix","p beta","H","break","lbreak","hbreak","bbreak","lform","hform","bform","wato h","wato b","llhb","hlhb","blhb");
printf WLHB ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n","resid","AA","loop","helix","beta","p loop","p helix","p beta","H","break","lbreak","hbreak","bbreak","lform","hform","bform","wato h","wato b","llhb","hlhb","blhb");
lhbdescr ();
for $j(1..$#ord){# wydruk
printf LHB ("%s\t%s\t%3.1f\t%3.1f\t%3.1f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",$j,$struct[$j],$lstate[$j],$hstate[$j],$bstate[$j],$p_lstate[$j],$p_hstate[$j],$p_bstate[$j],$hfobia1[$j],$brejker[$j],$lbrejker[$j],$hbrejker[$j],$bbrejker[$j],$lformer[$j],$hformer[$j],$bformer[$j],$wathbrr[$j],$watbbrr[$j],$llhb[$j],$hlhb[$j],$blhb[$j]);
printf WLHB ("%s\t%s\t%3.1f\t%3.1f\t%3.1f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",$j,$struct[$j],$lstate[$j],$hstate[$j],$bstate[$j],$p_lstate[$j],$p_hstate[$j],$p_bstate[$j],$w_hfobia1[$j],$w_brejker[$j],$w_lbrejker[$j],$w_hbrejker[$j],$w_bbrejker[$j],$w_lformer[$j],$w_hformer[$j],$w_bformer[$j],$wathbrr[$j],$watbbrr[$j],$llhb[$j],$hlhb[$j],$blhb[$j]);
}
close (LHB);
}

for $j(0..$naa){# wydruk
#print "$wathbrr[$j]\n";

printf ("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%6.3f\t%7.3f\t%7.3f\t%7.3f\n",$fname,$j,$struct[$j],$brejker[$j],$lbrejker[$j],$hbrejker[$j],$bbrejker[$j],$lformer[$j],$hformer[$j],$bformer[$j],$wathbrr[$j],$watbbrr[$j],$hfobia1[$j],$llhb[$j],$hlhb[$j],$blhb[$j]);
printf LHB ("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%6.3f\t%7.3f\t%7.3f\t%7.3f\n",$fname,$j,$struct[$j],$brejker[$j],$lbrejker[$j],$hbrejker[$j],$bbrejker[$j],$lformer[$j],$hformer[$j],$bformer[$j],$wathbrr[$j],$watbbrr[$j],$hfobia1[$j],$llhb[$j],$hlhb[$j],$blhb[$j]);
printf WLHB ("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%6.3f\t%7.3f\t%7.3f\t%7.3f\n",$fname,$j,$struct[$j],$w_brejker[$j],$w_lbrejker[$j],$w_hbrejker[$j],$w_bbrejker[$j],$w_lformer[$j],$w_hformer[$j],$w_bformer[$j],$wathbrr[$j],$watbbrr[$j],$w_hfobia1[$j],$llhb[$j],$hlhb[$j],$blhb[$j]);
}
close (LHB);
}#koniec
}
}

$koniec=localtime(time());
printf GOUT "Run completed. Time: $koniec\n";
close (AOUT);
close (GOUT);

