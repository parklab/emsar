#!/usr/bin/perl

if(@ARGV<3) { print "usage: $0 fastafile segmentfile g2tfile\n"; exit; }
my $fastafile = shift @ARGV;
my $segmentfile = shift @ARGV;
my $g2tfile = shift @ARGV;

open IN, $g2tfile or die "Can't open g2t file $g2tfile\n";
while(<IN>){
  chomp;
  my ($g,$t) = split/\t/;
  $t2g{$t}=$g;  ## transcript-to-gene mapping
  $nIsoforms{$g}++;  ## number of isoforms for each gene
}
close IN;



$header='';
$seq='';
open IN, $fastafile or die "Can't open fasta file $fastafile\n";
while(<IN>){
  chomp;
  if(/^>/) { 
    if($header ne '' && $seq ne '') { 
      $tlen{$header}=length($seq);  # transcript length
      $GC{$header} = ($seq=~tr/CG/CG/) / length($seq); # GC content
    }
    ($header)=split/\s+/,$'; $seq=''; 
  }
  else {
     $seq.=$_;
  }
}
close IN;

# last entry
$tlen{$header}=length($seq);  # transcript length
$GC{$header} = ($seq=~tr/CG/CG/) / length($seq); # GC content





open SEG, $segmentfile or die "Can't open segments file $segmentfile\n";
<SEG>;
while(<SEG>){
  my ($tnames,$EUMA) = (split/\t/)[3,4];
  my @tarray = split/\+/,$tnames;
  my $multi=0;
  if(@tarray>1){
    my $g1 = $t2g{$tarray[0]};
    for($i=1;$i<=$#tarray;$i++){
       if($t2g{$tarray[$i]} ne $g1) { $multi=1; last; }
    }
  }
  if($multi==1) { for($i=0;$i<=$#tarray;$i++){ $MULTI{$tarray[$i]}+=$EUMA; } }
  else { for($i=0;$i<=$#tarray;$i++){ $SINGLE{$tarray[$i]}+=$EUMA; } }
  if($multi==0 && @tarray==$nIsoforms{$t2g{$tarray[0]}}) { $gEUMA{$t2g{$tarray[0]}}+=$EUMA; }

  if(@tarray==1){
     $unique_len{$tarray[0]} = $EUMA;
  }

}
close SEG;




print "transcript_id\tgene\ttranscript_length\tGC_content\tnIsoforms\ttotal_effective_length\tisoform_unique_length\tgene_unique_length\tmulti_gene_length\tgene_unique_isoform_common_length\tisoform_unique_proportion\tgene_unique_proportion\tgene_unique_isoform_common_proportion\n";
for my $t (keys %t2g){
  if(!exists $unique_len{$t}) { $unique_len{$t}=0; }
  if(!exists $MULTI{$t}) { $MULTI{$t}=0; }
  if(!exists $SINGLE{$t}) { $SINGLE{$t}=0; }
  if(!exists $gEUMA{$t2g{$t}}) { $gEUMA{$t2g{$t}}=0; }
  my $totalEUMA = $MULTI{$t}+$SINGLE{$t};
  my $isoform_uniq_p = $totalEUMA>0?$unique_len{$t}/$totalEUMA:'NA';
  my $gene_uniq_p = $totalEUMA>0?$SINGLE{$t}/$totalEUMA:'NA';
  my $gEUMA_p = $totalEUMA>0?$gEUMA{$t2g{$t}}/$totalEUMA:'NA';
  print "$t\t$t2g{$t}\t$tlen{$t}\t$GC{$t}\t$nIsoforms{$t2g{$t}}\t$totalEUMA\t$unique_len{$t}\t$SINGLE{$t}\t$MULTI{$t}\t$gEUMA{$t2g{$t}}\t$isoform_uniq_p\t$gene_uniq_p\t$gEUMA_p\n";
}


