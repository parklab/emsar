#!/usr/bin/perl
if(@ARGV<2) { print "usage: $0 g2tfile FPKMfile\n"; exit; }
my $g2tfile = shift @ARGV;
my $FPKMfile = shift @ARGV;

open $G2T, $g2tfile or die "Can't open g2tfile $g2tfile\n";
while(<$G2T>){
  chomp;
  my ($g,$t) = split/\t/;
  $t2g{$t}=$g;
}
close $G2T;


open $FPKM, $FPKMfile or die "Can't open FPKMfile $FPKMfile\n";
<$FPKM>;
while(<$FPKM>){
  chomp;
  my ($t,$fpkm,$readcount,$tpm) = (split/\t/)[0,1,4,6];
  $gFPKM{$t2g{$t}}+=$fpkm;
  $gReadcount{$t2g{$t}}+=$readcount;
  $gTPM{$t2g{$t}}+=$tpm;
}
close $FPKM;

print "geneID\tFPKM\tiReadcount\tiReadcount.int\tTPM\n";
for my $g (keys %gFPKM){
  $gReadcount_int = roundoff($gReadcount{$g});
  print "$g\t$gFPKM{$g}\t$gReadcount{$g}\t$gReadcount_int\t$gTPM{$g}\n";
}

sub roundoff {
  ($_[0]-int($_[0])>=0.5)?int($_[0])+1:int($_[0]);
}

