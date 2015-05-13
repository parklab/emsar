#!/usr/bin/perl

my @files = @ARGV; #gReadcount files or other files in the same format.

for my $i (0..$#files){
 my $file = $files[$i];
 open $IN, $file or die "Can't open file $file\n";
 <$IN>;
 while(<$IN>){
   chomp;
   my ($id,$val) = (split/\t/)[0,3];
   $hash{$id}[$i]=$val;
 }
 close $IN;
}


# find longest common suffix for files
$minlen=1000000;
for my $file (@files){
  $minlen = length($file) if length($file)<$minlen;
}

my %suffs;
my $k=0;
do {
 %suffs=();
 $k++;
 for my $file (@files){
   $suff = substr($file,length($file)-$k);
   $suffs{$suff}=1;
 }
} while ( scalar(keys(%suffs))==1 &&  $k<$minlen);

my ($commonsuffix) = substr($files[0],length($files[0])-($k-1));
if(scalar(@files) == 1) { $commonsuffix = ''; }

# title is between "/" and common suffix.
my @title=();
for my $file (@files){
  if($file=~ /([^\/]+)$commonsuffix$/) { push @title, $1; }
  else { push @title, $file; }
}


## printing
print "ID\t@title\n";
$"="\t";
for my $id (keys %hash){
  for my $j (0..$#files){
    if(!defined $hash{$id}[$j]) {$hash{$id}[$j]=0; }
  }
  print "$id\t@{$hash{$id}}\n";
}

