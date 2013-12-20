package LHBclass;
use strict;
use warnings;

require Exporter;
use vars qw(@ISA @EXPORT $VERSION);
use HingeAlchemy::ConstHAlch;
our @ISA = qw(Exporter);
our @EXPORT = qw(vrclass sumvrclass);
$VERSION=1.0;

sub vrclass{
my (@Ra1_r,@Ra1_l,@Va1_r,@Va1_l,@l_ccl,@b_ccl,@h_ccl,@l_pcl,@b_pcl,@h_pcl,$R,$V);
my (@Ra1_1,@Ra1_2,@Ra1_3,@Ra1_4,@Ra1_5,@Ra1_6,@Ra1_7,@Ra1_8,@Ra1_9);
my (@Va1_1,@Va1_2,@Va1_3,@Va1_4,@Va1_5,@Va1_6,@Va1_7,@Va1_8,@Va1_9);
my ($k,$l,$z);
my ($ns,$mnr,$mlen,$mottyp,$minl,$motif,$lsum,$bsum,$hsum,$l_pcl,$h_pcl,$b_pcl,$l_totccl,$h_totccl,$b_totccl)=@_;
$Ra1_l[1]=$Ra1_1[0]-$ns*$Ra1_1[1];
$Ra1_r[1]=$Ra1_1[0]+$ns*$Ra1_1[1];
$Ra1_l[2]=$Ra1_2[0]-$ns*$Ra1_2[1];
$Ra1_r[2]=$Ra1_2[0]+$ns*$Ra1_2[1];
$Ra1_l[3]=$Ra1_3[0]-$ns*$Ra1_3[1];
$Ra1_r[3]=$Ra1_3[0]+$ns*$Ra1_3[1];
$Ra1_l[4]=$Ra1_4[0]-$ns*$Ra1_4[1];
$Ra1_r[4]=$Ra1_4[0]+$ns*$Ra1_4[1];
$Ra1_l[5]=$Ra1_5[0]-$ns*$Ra1_5[1];
$Ra1_r[5]=$Ra1_5[0]+$ns*$Ra1_5[1];
$Ra1_l[6]=$Ra1_6[0]-$ns*$Ra1_6[1];
$Ra1_r[6]=$Ra1_6[0]+$ns*$Ra1_6[1];
$Ra1_l[7]=$Ra1_7[0]-$ns*$Ra1_7[1];
$Ra1_r[7]=$Ra1_7[0]+$ns*$Ra1_7[1];
$Ra1_l[8]=$Ra1_8[0]-$ns*$Ra1_8[1];
$Ra1_r[8]=$Ra1_8[0]+$ns*$Ra1_8[1];
$Ra1_l[9]=$Ra1_9[0]-$ns*$Ra1_9[1];
$Ra1_r[9]=$Ra1_9[0]+$ns*$Ra1_9[1];
$Va1_l[1]=$Va1_1[0]-$ns*$Va1_1[1];
$Va1_r[1]=$Va1_1[0]+$ns*$Va1_1[1];
$Va1_l[2]=$Va1_2[0]-$ns*$Va1_2[1];
$Va1_r[2]=$Va1_2[0]+$ns*$Va1_2[1];
$Va1_l[3]=$Va1_3[0]-$ns*$Va1_3[1];
$Va1_r[3]=$Va1_3[0]+$ns*$Va1_3[1];
$Va1_l[4]=$Va1_4[0]-$ns*$Va1_4[1];
$Va1_r[4]=$Va1_4[0]+$ns*$Va1_4[1];
$Va1_l[5]=$Va1_5[0]-$ns*$Va1_5[1];
$Va1_r[5]=$Va1_5[0]+$ns*$Va1_5[1];
$Va1_l[6]=$Va1_6[0]-$ns*$Va1_6[1];
$Va1_r[6]=$Va1_6[0]+$ns*$Va1_6[1];
$Va1_l[7]=$Va1_7[0]-$ns*$Va1_7[1];
$Va1_r[7]=$Va1_7[0]+$ns*$Va1_7[1];
$Va1_l[8]=$Va1_8[0]-$ns*$Va1_8[1];
$Va1_r[8]=$Va1_8[0]+$ns*$Va1_8[1];
$Va1_l[9]=$Va1_9[0]-$ns*$Va1_9[1];
$Va1_r[9]=$Va1_9[0]+$ns*$Va1_9[1];
for $l (1..9){
$l_ccl[$l]=0;
$h_ccl[$l]=0;
$b_ccl[$l]=0;
$l_pcl[$l]=0;
$h_pcl[$l]=0;
$b_pcl[$l]=0;
}
for $k (0..$mnr){
 if($$mlen[$k]>=$minl){
  if($$mottyp[$k] eq 'L'){#loops
	for $z (1..$$mlen[$k]){
	$R=$$motif[$k][$z][3];
	$V=$$motif[$k][$z][4];
	 for $l (1..9){
#	 print "$R $V $Ra1_l[$l] $Ra1_r[$l] $Va1_l[$l] $Va1_r[$l]\n";
  	 if($R >= $Ra1_l[$l] and  $R <= $Ra1_r[$l] and $V >= $Va1_l[$l] and $V <= $Va1_r[$l]){$l_ccl[$l]++;}
	 }
	}
  }
  if($$mottyp[$k] eq 'H'){#helices
	for $z (1..$$mlen[$k]){
	$R=$$motif[$k][$z][3];
	$V=$$motif[$k][$z][4];
	 for $l (1..9){
	 if($R >= $Ra1_l[$l] and  $R <= $Ra1_r[$l] and $V >= $Va1_l[$l] and $V <= $Va1_r[$l]){$h_ccl[$l]++;}
	 }
	}
  }
  if($$mottyp[$k] eq 'B'){#betas
	for $z (1..$$mlen[$k]){
	$R=$$motif[$k][$z][3];
	$V=$$motif[$k][$z][4];
	 for $l (1..9){
  	 if($R >= $Ra1_l[$l] and  $R <= $Ra1_r[$l] and $V >= $Va1_l[$l] and $V <= $Va1_r[$l]){$b_ccl[$l]++;}
	 }
	}
  }
 }
}
#counts probabilities
for $l (1..9){
$$l_totccl[$l]=0;
$$h_totccl[$l]=0;
$$b_totccl[$l]=0;
}

for $l (1..9){
 if($$lsum>0){
 $$l_pcl[$l]=$l_ccl[$l]/$$lsum;}
 if($$hsum>0){
 $$h_pcl[$l]=$h_ccl[$l]/$$hsum;}
 if($$bsum>0){
 $$b_pcl[$l]=$b_ccl[$l]/$$bsum;}
$$l_totccl[$l]=$$l_totccl[$l]+$l_ccl[$l];
$$h_totccl[$l]=$$h_totccl[$l]+$h_ccl[$l];
$$b_totccl[$l]=$$b_totccl[$l]+$b_ccl[$l];
}
}

sub sumvrclass{
my ($l);
my ($lsumtot,$hsumtot,$bsumtot,$l_totccl,$h_totccl,$b_totccl,$l_totpcl,$h_totpcl,$b_totpcl)=@_;

  for $l (1..9){
	$$l_totpcl[$l]=0;
	$$h_totpcl[$l]=0;
	$$b_totpcl[$l]=0;
  }
  for $l (1..9){
	if($$lsumtot>0){#loops
	  $$l_totpcl[$l]=$$l_totccl[$l]/$$lsumtot;}
	if($$hsumtot>0){#helices
	  $$h_totpcl[$l]=$$h_totccl[$l]/$$hsumtot;}
	if($$bsumtot>0){#betas
	  $$b_totpcl[$l]=$$b_totccl[$l]/$$bsumtot;}
  }
}
1;