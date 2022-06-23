### Find chromosome specific sgRNAs of hg19 (fasta);
### Usage: perl FindChrSpecificRepSeq.pl hg19.fa hg19.chrN.sgRNA chrN

use warnings;
use strict;
use List::MoreUtils qw(any);
use Parallel::ForkManager;

my ($fasta,$outfile,$target) = @ARGV;   ## $target=chrNAME
open (FASTA, $fasta) or die "error(input):$!";  ### hg19(mm10).fasta
open (OUTFILE, ">", $outfile) or die "error(output):$!";

my $MAX_PROCESSES=8;
my %seq;
my $header;
my $seq="";
foreach(<FASTA>){
	chomp;
	if(/^>/){
	   $header=(split/>/,$_)[1];
  }else{
       $seq=$_;
       $seq{$header} = uc($seq);
  }
}

close(FASTA);

foreach my $elm(sort keys %seq){
	if ($elm =~ /$target/){  ### when just want to search for one chromosome;
     my @sgRNAs;
	   for (my $n = 0; $n < length($seq{$elm})-23; $n++){
	     my $tmpstr=substr($seq{$elm},$n,23);
       if ( any { $_ eq $tmpstr } @sgRNAs ) {
           next;
       }
	     if($tmpstr !~ /.*N.*/){
            if($tmpstr =~ /.*GG$/ || $tmpstr =~ /^CC.*/){
               my ($count,@pos) = NumofStr($seq{$elm},$tmpstr);
               if($count > 10){
               	  my $hit=0;
                  my $pm = new Parallel::ForkManager($MAX_PROCESSES);
               	  for my $chr(sort keys %seq){
               	  	  if($chr !~ $target){
                        my $pid = $pm->start and next;
               	  	  	if( NumofHit($seq{$chr},$tmpstr) !=0){
                            $hit++;
                         }
                         $pm->finish;
               	  	  }
               	  }
                  $pm->wait_all_children;
               	  if ($hit==0){
               	  	  my $chrpos="";
               	  	  foreach my $pos(sort @pos){
                         $chrpos.=$pos.";";
               	  	  }
               	  	  print OUTFILE "$elm\t$tmpstr\t$count\t$chrpos\n";
                      push @sgRNAs, $tmpstr;
               	  }
               }      	
	        }
	     }
	   }
 }
}

sub NumofStr{
	my $t=0;
	my ($sequence, $query)=@_;
	my $loci=index($sequence,$query);
	my @loci;
    while($loci != -1){
          $t=$t+1;
          push @loci,$loci;
          $loci=index($sequence,$query,$loci+1);
	}     
    return ($t,@loci);
}
  
sub NumofHit{
	my $t=0;
	my ($sequence, $query)=@_;
	my $loci=index($sequence,$query);
  if($loci != -1){
          $t=$t+1;
#          $loci=index($sequence,$query,$loci+1);
	}     
  return $t;
}      

