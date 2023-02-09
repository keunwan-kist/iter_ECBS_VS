use strict;
use Data::Dumper;
use Storable;
use List::Util qw(shuffle);

my $DAT_FILES="../data_files";		# data file directory 

## files used in data_files dir ...
# $DAT_FILES/mat.storable.data
# $DAT_FILES/uniprot_fam_map_fwd_v2.dat
# $DAT_FILES/uniprot_fam_map_rev_v2.dat

# parameters to define chemical pairs  
my $N_cut_target = 25;			# upper limit (# of chems) per class
my $N_cut_family = 100;			# upper limit (# of chems) per fam_id
my $neg_multiplier = 2;			#  multiplier for true negatives
my $per_neg_cnt = 10;				# parameter to control (# of chem) for targets 
my $new_pos_overweight = 5;	# weight factor for new actives to newP - randN 

my $cur=`pwd`;chomp $cur;

if(@ARGV<6){
	print "This script produces the filtered-negatives and positives chem pairs by evo groups\n";
	die "[].pl [int.drug_target_map.txt] [new mol txt] [BindingDB_All_2D_v2.1uM_Affinity] [target ID] [family ID] [out dir name] [update scheme: e.g. PP_NP or PP_NP_NN]\n";
}

my $dt_map = $ARGV[0];	# simplified chem-tar int. data, BindingDB + DrugBank 
my $new_mol = $ARGV[1];	# new exp. data 
my $aff = $ARGV[2];			# parsed binding affinity data from BindingDB
my $target_name = $ARGV[3];	
my $family_info = $ARGV[4];	# family names if needed 
my @family_names=split(/,/,$family_info);
my $out_dir = $ARGV[5];
my $update_scheme = $ARGV[6];		# PP, NP, NN, or combination of them (e.g. PP_NP_NN) 

my %AF;	
open(IN,$aff) or die "can't open file:$aff\n";
while(<IN>){
  chomp;
  #my ($targ,$ch,$assay,$type,$val_parsed,$val,$unit)=split(/\t/,$_);
	#next if($val_parsed eq "na");
  #$AF{$targ}{$ch}{"$type;$unit"}{$val_parsed}=1;		# can handle multiple values for same aff. 
	my ($r_id,$d_id,$class,$val)=split;
  next if($val =~ />/);
  $val=~s/[><]//g;
  $AF{$d_id}=$val;  # only high affinity ones
}
close IN;

my %FeatVec;
print STDERR "\nPre-made feature file loading..\n";
## loading final filtered mat file 
my $hashref = retrieve("$DAT_FILES/mat.storable.data");
%FeatVec = %{$hashref};
print STDERR "\nMat file loading finished...\n";


my $promiscous_cut_targ = 60;		# excluding promiscous targets, cut defined by pval 0.0001
my $promiscous_cut_fam = 17;		# excluding promiscous family info
my %promiscous_chems;
my %target_chems;

my (%DT_f,%DT_r);
open(IN,$dt_map) or die "can't open file:$_\n";
while(<IN>){
	chomp;
	my ($uniprot,$chemid)=split;
	next if(not defined $FeatVec{$chemid}); #print STDERR "WARN $chemid is not defined in Mat file!\n";
	
	$DT_f{$uniprot}{$chemid}=1;
	$DT_r{$chemid}{$uniprot}=1;
	
	$target_chems{$chemid} =1 if($uniprot eq $target_name);
}
close IN;

my %exp_Pos;
my %exp_Neg;
open(IN,$new_mol) or die "can't open file:$_\n";
while(<IN>){
	chomp;
	#BTB04770 93.47858529 0.9275 1 P14679
	my ($chem,$activity,$dump,$class,$targ)=split;
	if(not defined $FeatVec{$chem}){
		die "New mol $chem is not defined in FeatVec, please add featvec to mat.sto file\n";
	}
	if($class == 1){			# new active
		$exp_Pos{$chem}=$_;
		$DT_f{$targ}{$chem}=1;	# update of chem-targ interactions
		$DT_r{$chem}{$targ}=1;  # update of chem-targ interactions
		$target_chems{$chem}=1;
	}else{								# new inactive 
		$exp_Neg{$chem}=$_;	
	}
}
close IN;

my @target_chems_ary = keys %target_chems;
my @exp_pos_chems_ary = keys %exp_Pos;
my @exp_neg_chems_ary = keys %exp_Neg;

print STDERR "------------------------------------------------------\n";
print STDERR "Direct $target_name target-binding chems: @target_chems_ary\n";
print STDERR "New exp. $target_name target-binding chems: @exp_pos_chems_ary\n";
print STDERR "New exp. $target_name NOT-binding chems: @exp_neg_chems_ary\n";
print STDERR "------------------------------------------------------\n";


my %gDrugMap;
my %gGroupMap;

for my $chemid (keys %DT_r){		# DT_r is updated with exp data 
	my @uniprots = keys %{$DT_r{$chemid}};
	my $uniprots_n = @uniprots;
	if($uniprots_n > $promiscous_cut_fam){		 # tight
		$promiscous_chems{family}{$chemid}=$uniprots_n;
	}	
	if($uniprots_n > $promiscous_cut_targ){		# loose 
		print STDERR "[Promiscous Outlier] $chemid - $uniprots_n targets .. excluded $promiscous_cut_targ/$promiscous_cut_fam\n";
		$promiscous_chems{target}{$chemid}=$uniprots_n;
		next;
	}
	for my $uniprot (@uniprots){
		$gDrugMap{target}{$uniprot}{$chemid}=1;	
		$gGroupMap{target}{$chemid}{$uniprot}=1;	
	}
}

print STDERR "Family Map loading....\n";
my $fm_f = retrieve("$DAT_FILES/uniprot_fam_map_fwd_v2.dat");
my $fm_r = retrieve("$DAT_FILES/uniprot_fam_map_rev_v2.dat");

#### def. for evolutionary information from diff DB .. 
my @categories = ("target","G3DSA","TIGR","smart","supfam","pfam","fam","print");		
#############################################################

my %cat_hash;
for my $c (@categories){ $cat_hash{$c}=1;	}
my @categories_nt = @categories; shift @categories_nt;	# w.o. TARGET


for my $uniprot (keys %DT_f){
	for my $class (keys %{$fm_r->{$uniprot}}){		# family mapping 
		next if(not defined $cat_hash{$class});
		for my $fam_id (keys %{$fm_r->{$uniprot}{$class}}){
			for my $chem (keys %{$DT_f{$uniprot}}){
				$gDrugMap{$class}{$fam_id}{$chem}=1;
				$gGroupMap{$class}{$chem}{$fam_id}=1;
			}	# chem..
		}		# fam id
	}	# categories
}

##################################################

my @aff_cut = (50000,10000,5000,1000,500,100,50,10,5,1,0.5,0.1,0.05,0.01);		# 50uM ~ 1nM

print STDERR "------------------------------------------------------\n";
print STDERR "	Writing data to $cur/$out_dir\n";
print STDERR "------------------------------------------------------\n";
mkdir("$cur/$out_dir");
chdir("$cur/$out_dir");

# making.. uniprot-related FAMILY variables
my %related_chems;
my %related_chems_reduced;
my %drugbank_n;
my %evo_related_targets;

for my $class (@categories_nt){
	
	my @fam_ids = keys %{$fm_r->{$target_name}{$class}};	# only for the target-related evo..
	print STDERR "! WARNING $class has no info for $target_name, will use @family_names\n" if(@fam_ids<1);
	
	push @fam_ids,@family_names if($family_info ne "NA");
	print STDERR "Evolutionary Infomation used - @fam_ids\n";

	for my $fam_id (@fam_ids){
			
		for my $t (keys %{$fm_f->{$class}{$fam_id}}){
			$evo_related_targets{$t}=1;		# evo-related targets!
		}

		my @chs = keys %{$gDrugMap{$class}{$fam_id}};
		for my $ch (@chs){
			$drugbank_n{$class}{$fam_id}{$ch}=1 if($ch=~/DB/);		# drugbank, no aff info
			$related_chems{$class}{$fam_id}{$ch} = 1;	
		}
		my @drugbank_fam_ch = keys %{$drugbank_n{$class}{$fam_id}};
		my $drugbank_fam_n = @drugbank_fam_ch;

		my @cur_fam_ch = keys %{$related_chems{$class}{$fam_id}};
		my $cur_fam_ch_n = @cur_fam_ch;
		my $original_cur_fam_ch_n = $cur_fam_ch_n;

		if($drugbank_fam_n > $N_cut_family){		# when drugbank alrady has enough chems .. 

			@drugbank_fam_ch = shuffle @drugbank_fam_ch;
			my @rand_chems = @drugbank_fam_ch[0..$N_cut_family];
			
			for my $ch (@rand_chems){
				$related_chems_reduced{$class}{$fam_id}{$ch}=1;	# contain only reduced info.. 
			}
			for my $ch (keys %exp_Pos){		# add exp new pos chems, okay for all target-related evo 
				$related_chems_reduced{$class}{$fam_id}{$ch}=1;
			}

			@cur_fam_ch = keys %{$related_chems_reduced{$class}{$fam_id}};
			$cur_fam_ch_n =  @cur_fam_ch;
			print STDERR "=> Evolutionary $class $fam_id $cur_fam_ch_n enough randomly sampled from $drugbank_fam_n drugbank\n";	

		}else{	# when drugbank is not enough .. 
			
			if($cur_fam_ch_n > $N_cut_family){		# but when the binding db exceeds the cut
				
				for my $ch (@drugbank_fam_ch){	# include drugbank chems first 
					$related_chems_reduced{$class}{$fam_id}{$ch}=1;
				}

				my $limit_flag=0;	

				for my $af_cut (reverse @aff_cut){
					for my $ch (@cur_fam_ch){
						next if($ch=~/DB/);
						next if(not defined $AF{$ch});
						if($AF{$ch} < $af_cut){
							$related_chems_reduced{$class}{$fam_id}{$ch}=1;
							my $cnt = keys %{$related_chems_reduced{$class}{$fam_id}};
							if($cnt >= $N_cut_family){
								$limit_flag=1;last;
							}
						}	
					}
					last if($limit_flag==1);	# only include the fixed number
				}
				
				for my $ch (keys %exp_Pos){	# add exp new pos chems, anyway
					$related_chems_reduced{$class}{$fam_id}{$ch}=1;
				}

				@cur_fam_ch = keys %{$related_chems_reduced{$class}{$fam_id}}; 
				$cur_fam_ch_n = @cur_fam_ch;
				print STDERR "=> Evolutionary $class $fam_id $cur_fam_ch_n randomly sampled with Affcut, including $drugbank_fam_n drugbank\n";	

			}else{		# when binding db number is not enough 

				for my $ch (keys %{$related_chems{$class}{$fam_id}}){		# include all
					$related_chems_reduced{$class}{$fam_id}{$ch}=1;
				}

				for my $ch (keys %exp_Pos){		# add exp new pos chems, okay for all target-related evo 
					$related_chems_reduced{$class}{$fam_id}{$ch}=1;
				}

				@cur_fam_ch = keys %{$related_chems_reduced{$class}{$fam_id}};
				$cur_fam_ch_n =  @cur_fam_ch;
				print STDERR "=> Evolutinary $class $fam_id $cur_fam_ch_n < $N_cut_family\n";
			}
		}	
		print STDERR "$class $fam_id Map - now $cur_fam_ch_n (start $original_cur_fam_ch_n) added, including DrugBank $drugbank_fam_n\n";
		print STDERR "[Final Evolutionary $class Chems]: @cur_fam_ch\n";
	}		# fam id

}	# categories


print STDERR "\n-------------Evolutionarily-related Target Update starting from here -----------\n"; 
$evo_related_targets{$target_name}=1;		# to make sure 

for my $uniprot (keys %evo_related_targets){		# only for all evo-related targets... 

	for my $ch (keys %{$DT_f{$uniprot}}){
		$related_chems{'target'}{$uniprot}{$ch} = 1;
		$drugbank_n{'target'}{$uniprot}{$ch}=1 if($ch=~/DB/);
	} 
	
	my @drugbank_targ_ch = keys %{$drugbank_n{'target'}{$uniprot}};
	my $drugbank_targ_ch_n = @drugbank_targ_ch;

	my @cur_targ_ch = keys %{$related_chems{'target'}{$uniprot}}; 
	my $cur_targ_ch_n = @cur_targ_ch; 
	my $original_cur_targ_ch_n = $cur_targ_ch_n; 

	if($drugbank_targ_ch_n > $N_cut_target){
		@drugbank_targ_ch = shuffle @drugbank_targ_ch;
		my @rand_chems = @drugbank_targ_ch[0..$N_cut_target];
	
		for my $ch (@rand_chems){
			$related_chems_reduced{'target'}{$uniprot}{$ch}=1;
		}
		if($uniprot eq $target_name){
			for my $ch (keys %exp_Pos){	# add new chems only for the specified targete 
				$related_chems_reduced{'target'}{$uniprot}{$ch}=1;
			}
		}

		@cur_targ_ch = keys %{$related_chems_reduced{'target'}{$uniprot}}; 
		$cur_targ_ch_n =  @cur_targ_ch;
		print STDERR "=> Target $uniprot $cur_targ_ch_n enough randomly sampled from $drugbank_targ_ch_n drugbank\n";
	
	}else{		# when drugbank has not enough number 
	
		if($cur_targ_ch_n > $N_cut_target){			# but binding db exceeds the number
			
			for my $ch (@drugbank_targ_ch){
				$related_chems_reduced{'target'}{$uniprot}{$ch}=1;
			}
					
			my $limit_flag=0;	
		
			for my $af_cut (reverse @aff_cut){	# from high aff
				for my $ch (@cur_targ_ch){
					next if($ch=~/DB/);
					next if(not defined $AF{$ch});
					if($AF{$ch} < $af_cut){
						$related_chems_reduced{'target'}{$uniprot}{$ch}=1;
						my $cnt = keys %{$related_chems_reduced{'target'}{$uniprot}};
						if($cnt >= $N_cut_target){
							$limit_flag=1;
							last;
						}
					}	
				}
				last if($limit_flag==1);
			}
			
			if($uniprot eq $target_name){
				for my $ch (keys %exp_Pos){			# add new mol info anyway 
					$related_chems_reduced{'target'}{$uniprot}{$ch}=1;
				}
			}

			@cur_targ_ch = keys %{$related_chems_reduced{'target'}{$uniprot}}; 
			$cur_targ_ch_n = @cur_targ_ch;
			print STDERR "=> Target $uniprot $cur_targ_ch_n random-sampled with Affinity Cutoff, DrugBank $drugbank_targ_ch_n\n";

		}else{	# when number is not enough ..

			for my $ch (keys %{$related_chems{target}{$uniprot}}){
				$related_chems_reduced{'target'}{$uniprot}{$ch}=1;
			}
			if($uniprot eq $target_name){
				for my $ch (keys %exp_Pos){			# add new mol info anyway 
					$related_chems_reduced{'target'}{$uniprot}{$ch}=1;
				}
			}

			@cur_targ_ch = keys %{$related_chems_reduced{'target'}{$uniprot}}; 
			$cur_targ_ch_n = @cur_targ_ch;
			print STDERR "=> Target $uniprot $cur_targ_ch_n < $N_cut_target\n";
		}
	}# else 

	print STDERR "[Final Target $uniprot Stat] $cur_targ_ch_n ($original_cur_targ_ch_n) added, including DrugBank $drugbank_targ_ch_n : @cur_targ_ch\n";

} # Now, "related_chems_reduced" is ready for use 


###### Delete promiscous chems ########### 
for my $class (keys %related_chems_reduced){
	for my $fam_id (keys %{$related_chems_reduced{$class}}){
		my %deleted_ch;
		for my $ch  (keys %{$related_chems_reduced{$class}{$fam_id}}){
			if($class eq 'target' and defined $promiscous_chems{target}{$ch}){
				delete $related_chems_reduced{$class}{$fam_id}{$ch};
				$deleted_ch{$ch}=1;
			}
			if($class ne 'target' and defined $promiscous_chems{family}{$ch}){
				delete $related_chems_reduced{$class}{$fam_id}{$ch};
				$deleted_ch{$ch}=1;
			}
		}
		if(keys %{$related_chems_reduced{$class}{$fam_id}} == 0){
			print STDERR "Nothing left... so Considering promiscous chems more again...\n";
			for my $ch (keys %deleted_ch){
				$related_chems_reduced{$class}{$fam_id}{$ch}=1;
			}
		}
	}
}	

###### NOW it is not considering target unrelevant info in DrugMap
###### That is, it only considers the chemical info in related_chems_reduced var. 
print STDERR "\n---------------------------------------------------\n";
print STDERR "               Staring to write [X].data\n";
print STDERR "--------------------------------------------------------\n";
my @all_chems = keys %DT_r;

for my $c (@categories){

	if(not defined $related_chems_reduced{$c}){
		print STDERR "! WARNING $c has no info for $target_name\n";
		next;
	}	
	print STDERR "Processing $c feat to ./$c.data : pos and neg:[$per_neg_cnt x pos] \n";
	open(O,">$c.data") or die "can't open file:$c\n";
	
	my %negpair;
	my %pospair;

	my %direct_QSAR_data;

	if($c eq "target"){
		open(O4,">seed_$c.data") or die "can't open file:seed_$c.data\n";
	}


	my $max_per_neg_cnt=0;
	
	for my $group (keys %{$related_chems_reduced{$c}}){
		my @ds = keys %{$related_chems_reduced{$c}{$group}};	

		my $direct_target_flag = 0;
		$direct_target_flag = 1 if($c eq "target" and $group eq $target_name);
		
		my @tmp_ds;
		for my $d (@ds){		# no exp. chems for now 
			next if(defined $exp_Pos{$d} or defined $exp_Neg{$d});
			next if(not defined $FeatVec{$d});
			push @tmp_ds,$d;		
		}
		@ds = @tmp_ds;

		print STDERR "[$c default data w.o. new mol] $group : @ds\n";	
		
		my %pos_chems;	# for knowns in DB
		my $pos_pair_cnt=0;

		for(my $i=0;$i<@ds;$i++){	# for all drug-pairs 
			for(my $j=$i;$j<@ds;$j++){		# contain self! 
				next if(defined $pospair{$ds[$i]}{$ds[$j]});
				next if(defined $negpair{$ds[$i]}{$ds[$j]});
				
				my @bits = &AB_on_bit(\%FeatVec,$ds[$i],$ds[$j]);	
				print O "$ds[$i] $ds[$j] @bits 1\n";
			
				@{$direct_QSAR_data{Pos}{$ds[$i]}} = @{$FeatVec{$ds[$i]}} if($direct_target_flag==1);
        @{$direct_QSAR_data{Pos}{$ds[$j]}} = @{$FeatVec{$ds[$j]}} if($direct_target_flag==1);

				$pospair{$ds[$i]}{$ds[$j]}=1;
				$pospair{$ds[$j]}{$ds[$i]}=1;
				$pos_chems{$ds[$i]}=1;
				$pos_chems{$ds[$j]}=1;
				$pos_pair_cnt++;
			} # for j
		}	# for i	
		
		my @pos_chems_ary = keys %pos_chems;
		my $pos_cnt = @pos_chems_ary;
		
		print STDERR "=> Positive chems for $c $group ($pos_cnt):  @pos_chems_ary\n";

		$per_neg_cnt = int($pos_pair_cnt*$neg_multiplier/$pos_cnt) if($pos_cnt>0);	
		$max_per_neg_cnt = $per_neg_cnt if($c eq "target" and $per_neg_cnt > $max_per_neg_cnt);
		print STDERR "=> $c $group PER_NEG_CNT $per_neg_cnt (max $max_per_neg_cnt)\n";	

		my %neg_all_chems;
		for my $pch (keys %pos_chems){	# default Pos-Rand negative samples 
			my $ncnt=0;	
			while($ncnt <= $per_neg_cnt){
				my $rand_ch = $all_chems[ rand(@all_chems) ];
				next if(not defined $FeatVec{$rand_ch});
				next if(defined $pospair{$rand_ch}{$pch});
				next if(defined $negpair{$rand_ch}{$pch});
				next if( &all_group_matched(\%gGroupMap,$pch,$rand_ch) );
				
				$neg_all_chems{$rand_ch}=1;
				my @bits = &AB_on_bit(\%FeatVec,$pch,$rand_ch);		
				print O "$pch $rand_ch @bits 0\n";

				@{$direct_QSAR_data{Neg}{$rand_ch}} = @{$FeatVec{$rand_ch}} if($direct_target_flag==1);
				@{$direct_QSAR_data{Pos}{$pch}} = @{$FeatVec{$pch}} if($direct_target_flag==1);

				$negpair{$pch}{$rand_ch}=1;
				$negpair{$rand_ch}{$pch}=1;
				$ncnt++;
			}
		}		# pos chem - random neg pairs  	
		
		# for $c, $group .. for fam info, all are pos-samples 
		if($update_scheme =~ /PP/){			# Update P-TP, TP-TP
			print STDERR "-----> UPDATE [newP:knownP] by $update_scheme \n"; 
			my @exp_pos_chems = keys %exp_Pos;

			for(my $k1=0;$k1<@exp_pos_chems;$k1++){		# TODO- exp can be from diff targets..
				my $ch1 = $exp_pos_chems[$k1];
			
				for my $pch (keys %pos_chems){	 
					next if(defined $pospair{$pch}{$ch1});
					next if(defined $negpair{$pch}{$ch1});
					my @bits = &AB_on_bit(\%FeatVec,$pch,$ch1);		
					
					if($c ne "target"){		# for fam, all are positives 
						
						print O "$pch $ch1 @bits 1\n";
						$pospair{$pch}{$ch1}=1;
						$pospair{$ch1}{$pch}=1;

					}elsif($c eq "target"){
						
						if($direct_target_flag==1){	
							@{$direct_QSAR_data{Pos}{$pch}} = @{$FeatVec{$pch}};
							@{$direct_QSAR_data{Pos}{$ch1}} = @{$FeatVec{$ch1}};

							print O "$pch $ch1 @bits 1\n";  # for target, only consider a given uniprot for PP 

								$pospair{$pch}{$ch1}=1;
							$pospair{$ch1}{$pch}=1;
						}
					}
				} # pos chem 
			} # k1
		}	# PP end
		
		if($direct_target_flag==1 and ($update_scheme =~ /NP/ or $update_scheme =~ /PN/)){
			print STDERR "-----> UPDATE [trueNeg:knownP] TRUE NEG by $update_scheme \n"; 
			for my $nch (keys %exp_Neg){
				for my $pch (keys %pos_chems){	# P-TN, known actives for a target!!! 

					next if(defined $pospair{$pch}{$nch});
					next if(defined $negpair{$pch}{$nch});
					next if( &all_group_matched(\%gGroupMap,$pch,$nch) );

					my @bits = &AB_on_bit(\%FeatVec,$pch,$nch);		
					print O "$pch $nch @bits 0\n";
					@{$direct_QSAR_data{Pos}{$pch}} = @{$FeatVec{$pch}};
					@{$direct_QSAR_data{Neg}{$nch}} = @{$FeatVec{$nch}};

					$negpair{$pch}{$nch}=1;
					$negpair{$nch}{$pch}=1;
				}
			} 
		}	# TrueNeg <-> known pos negative samples 

	} # group end 

	if($update_scheme =~ /PP/){			# Update P-TP, TP-TP
		my @exp_pos_chems = keys %exp_Pos;
		print STDERR "-----> UPDATE [newP:newP] Pos by $update_scheme \n"; 
		print STDERR "-----> UPDATE [newP:random Neg] Neg by $update_scheme \n"; 

		for(my $k1=0;$k1<@exp_pos_chems;$k1++){
			my $ch1 = $exp_pos_chems[$k1];
			my $ncnt=0;
			while($ncnt <= $max_per_neg_cnt * $new_pos_overweight ){		# !!!!weight factor for exp pos!!!!
				my $rand_ch = $all_chems[ rand(@all_chems) ];
				next if(not defined $FeatVec{$rand_ch});

				next if(defined $pospair{$rand_ch}{$ch1});
				next if(defined $negpair{$rand_ch}{$ch1});
				next if( &all_group_matched(\%gGroupMap,$ch1,$rand_ch) );

				my @bits = &AB_on_bit(\%FeatVec,$ch1,$rand_ch);
				print O "$ch1 $rand_ch @bits 0\n";
				@{$direct_QSAR_data{Pos}{$ch1}} = @{$FeatVec{$ch1}};
				@{$direct_QSAR_data{Neg}{$rand_ch}} = @{$FeatVec{$rand_ch}};

				$negpair{$ch1}{$rand_ch}=1;
				$negpair{$rand_ch}{$ch1}=1;
				$ncnt++;
			}
			
			for(my $k2=$k1;$k2<@exp_pos_chems;$k2++){
				my $ch2 = $exp_pos_chems[$k2];
				next if(defined $pospair{$ch1}{$ch2});
				next if(defined $negpair{$ch1}{$ch2});

				my @bits = &AB_on_bit(\%FeatVec,$ch1,$ch2);		
				print O "$ch1 $ch2 @bits 1\n";
				
				@{$direct_QSAR_data{Pos}{$ch1}} = @{$FeatVec{$ch1}};
				@{$direct_QSAR_data{Pos}{$ch2}} = @{$FeatVec{$ch2}};

				$pospair{$ch1}{$ch2}=1;
				$pospair{$ch2}{$ch1}=1;
			} # k2
		}	# k1

	}	# PP end

#####################################################
# from here, can't update evo negatives for fam, since the concept is valid only for targets  			
#####################################################
	if($c eq "target"){	

		if($update_scheme =~ /NP/ or $update_scheme =~ /PN/){
			print STDERR "-----> UPDATE [TrueNeg:newP] by $update_scheme \n"; 
			for my $nch (keys %exp_Neg){
				for my $pch (keys %exp_Pos){	# TP-TN negative samples 

					next if(defined $pospair{$pch}{$nch});
					next if(defined $negpair{$pch}{$nch});
					next if( &all_group_matched(\%gGroupMap,$pch,$nch) );

					my @bits = &AB_on_bit(\%FeatVec,$pch,$nch);		
					@{$direct_QSAR_data{Pos}{$pch}} = @{$FeatVec{$pch}};
					@{$direct_QSAR_data{Neg}{$nch}} = @{$FeatVec{$nch}};

					print O "$pch $nch @bits 0\n";

					$negpair{$pch}{$nch}=1;
					$negpair{$nch}{$pch}=1;
				}
			}
		}		# NP end

		if($update_scheme =~ /NN/){
			print STDERR "-----> UPDATE [TrueNeg:random N] by $update_scheme \n"; 
			for my $nch (keys %exp_Neg){  # Rand-TN negative samples
				my $ncnt=0;
				while($ncnt <= $max_per_neg_cnt * $new_pos_overweight){		# weight factor for exp Neg
					my $rand_ch = $all_chems[ rand(@all_chems) ];
					next if(not defined $FeatVec{$rand_ch});
					next if(defined $pospair{$rand_ch}{$nch});
					next if(defined $negpair{$rand_ch}{$nch});
					next if( &all_group_matched(\%gGroupMap,$nch,$rand_ch) );

					my @bits = &AB_on_bit(\%FeatVec,$nch,$rand_ch);
					@{$direct_QSAR_data{Neg}{$nch}} = @{$FeatVec{$nch}};
					@{$direct_QSAR_data{Neg}{$rand_ch}} = @{$FeatVec{$rand_ch}};
					
					print O "$nch $rand_ch @bits 0\n";

					$negpair{$nch}{$rand_ch}=1;
					$negpair{$rand_ch}{$nch}=1;
					$ncnt++;
				}
			}   # chem end
		}	# NN end

		for my $ch (keys %{$direct_QSAR_data{Pos}}){
			my @feat = @{$direct_QSAR_data{Pos}{$ch}};
			print O4 "$ch @feat\n";
		}
		close O4;
	}

	close O;

}	# categories


#################################

sub all_group_matched {
	my ($GroupMap,$d1,$d2) = @_;
	
	for my $c (keys %{$GroupMap}){
		my $p = &is_group_matched($GroupMap,$d1,$d2,$c); 
		
		if( $p eq "0" ){
			#print "=== d1 $d1 d2 $d2 [$c:$p]\n" if($d1 eq "DB00528" and $d2 eq "DB01054") ;
		
		}else{
			#print "=== d1 $d1 d2 $d2 [$c:$p]\n" if($d1 eq "DB00528" and $d2 eq "DB01054") ;
			#print STDERR "$d1 $d2 matched by $p / $c : EXCLUDED\n";
			return 1;			
		}
	}
	return 0;
}

sub is_group_matched {
	my ($P,$d1,$d2,$cat) = @_;
	

	my @g1 = keys %{$P->{$cat}{$d1}};
	my @g2 = keys %{$P->{$cat}{$d2}};
	
	for my $p1 (@g1){
		for my $p2 (@g2){
			return $p1 if($p1 eq $p2);
		}
	}
	return 0;
}	

sub AB_on_bit {
  my ($D,$d1,$d2) = @_;
  my @d1_val= @{$D->{$d1}};
  my @d2_val= @{$D->{$d2}};
  my @d;
  for(my $k=0;$k<@d1_val;$k++){
		push @d, $d1_val[$k]+$d2_val[$k]; 	# A + B   
	}
  return @d;
}


