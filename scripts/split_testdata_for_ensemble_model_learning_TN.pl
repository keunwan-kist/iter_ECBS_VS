use strict;
use warnings;
use List::Util qw/shuffle/;

 
if(@ARGV<2){
	die "[].pl [directory for target.data] [new_mol1_2.txt] [CV]\n";
}

my $targ_dir = $ARGV[0];
my $new_mol = $ARGV[1];
my $numberofarrays = $ARGV[2] || 7;

my $targ_data_fn = "$targ_dir/target.data";
my $seed_data_fn = "$targ_dir/seed_target.data";

my %direct_chems;
open(IN,$seed_data_fn) or die "can't open file:$_\n";
while(<IN>){
	chomp;
	my ($chem,@vals)=split;
	$direct_chems{$chem}=1;
}
close IN;

my %exp_Pos;
my %exp_Neg;
open(IN,$new_mol) or die "can't open file:$_\n";
while(<IN>){
  chomp;
  #BTB04770 93.47858529 0.9275 1 P14679
  my ($chem,$activity,$dump,$class,$targ)=split;
  if($class == 1){
    $exp_Pos{$chem}=$_;
  }else{
    $exp_Neg{$chem}=$_;
  }
}
close IN;

for my $c (keys %exp_Pos){
	if(not defined $direct_chems{$c}){
		die "$c is not defined in direct_chems!!\n";;
	}
}

my %T;
open(IN,$targ_data_fn) or die "can't open file:$_\n";
while(<IN>){
	chomp;
	my @tmp=split;
	$T{"$tmp[0] $tmp[1]"}=$_;
}
close IN;

my @list_ary = keys %direct_chems;
my @shuffled_ary = shuffle @list_ary;
my @accts_split;
my $length =  sprintf("%5d",@shuffled_ary / $numberofarrays + 0.5);
my $start  =  0;

my @list_ary_neg = keys %exp_Neg;
my @shuffled_ary_neg = shuffle @list_ary_neg;
my $length_neg =  sprintf("%5d",@shuffled_ary_neg / $numberofarrays + 0.5);
my $start_neg  =  0;


for my $i (0 .. ($numberofarrays-1))
{
	my $end = ($i == $numberofarrays-1) ? $#shuffled_ary : $start + $length - 1;
	my $end_neg = ($i == $numberofarrays-1) ? $#shuffled_ary_neg : $start_neg + $length_neg - 1;
	@{$accts_split[$i]} = (@shuffled_ary[$start .. $end], @shuffled_ary_neg[$start_neg .. $end_neg]);
	#print "s-e $start $end\n";
	$start += $length;
	$start_neg += $length_neg;
}

my $cnt=@shuffled_ary;
my $cnt_neg=@shuffled_ary_neg;
my $tot_n = @accts_split;
my $set_num = sprintf("%5d",($cnt+$cnt_neg)/$tot_n + 0.5);
print STDERR "cnt/numberofarrays ($cnt+$cnt_neg)/$tot_n = $set_num\n";

my $splitn=1;
for my $p (@accts_split){
	my @arys = @{$p};			# --> num of shuffled data array for each split
	my $ary_n = @arys;
	#print "$p $ary_n\n";
	
	my %test_chems;
	for my $ch (@arys){
		$test_chems{$ch}=1;
	}

	my %test_set;
	my %train_set;

	for my $dp (keys %T){
		my ($d1,$d2)=split(/\s+/,$dp);
		if(defined $test_chems{$d1} or defined $test_chems{$d2}){
			$test_set{$T{$dp}} =1;
		}else{
			$train_set{$T{$dp}} =1;
		}
	}

	open(Otest,">$targ_dir/testset$splitn.target.data") or die "can't open file:\n";
	open(Otest_desc,">$targ_dir/testset$splitn.blind_chem.list") or die "can't open file:\n";
	open(Otrain,">$targ_dir/trainset$splitn.target.data") or die "can't open file:\n";

	for my $ch (keys %test_chems){
		my $cla = 1;
		$cla=0 if(defined $exp_Neg{$ch});
		print Otest_desc "$ch $cla\n";
	}

	for my $line (keys %test_set){		
		print Otest "$line\n";
	}
	for my $line (keys %train_set){		
		print Otrain "$line\n";
	}

	close IN;
	close Otest;
	close Otest_desc;
	close Otrain;
		
	$splitn++;
}

		



