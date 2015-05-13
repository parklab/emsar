#!/usr/bin/perl

if(@ARGV<2) { print "usage: $0 fpkm_file_dir g2t_file\n"; exit; }
$dir = shift @ARGV;
$g2t = shift @ARGV;

@fileprefix = split/\n/,`ls -1 $dir/*.fpkm |sed s/\.fpkm\$//g`;

my ($EMSARdir) = $0 =~/(.+)\/[^\/]+$/;


for my $pref (@fileprefix){
  `$EMSARdir/FPKM2gFPKM.pl $g2t $pref.fpkm > $pref.gfpkm`;
}

`$EMSARdir/merge_gReadcount.pl $dir/*.gfpkm > $dir/gReadcount.all`;

`$EMSARdir/merge_gTPM.pl $dir/*.gfpkm > $dir/TPM.all`;

