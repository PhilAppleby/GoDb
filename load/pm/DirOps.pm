package DirOps;

require Exporter;

use strict;
use warnings;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(print_file);

our %qsub_tmplt = (
  qsub_impute_cmd => "qsub -q c6100-long.q -o %s -e %s -t 1-%d:1 %s/single_impute_chr_process.sh %s %d",
);
sub print_file {
  my $dir = shift;
  my $file = shift;
  my @params = @_;
  print "$dir, $file\n";
  for my $parm (@params) {
    print "$parm\n";
  }
}

#-------------------------------------------------------------------
# Make a standard chunk file record
# chunk number, start and end fields forced to a standard precision
#-------------------------------------------------------------------
sub make_chunk_rec {
  my $ckno = shift;
  my $start = shift;
  my $end = shift;
  my $numsnps = shift;
  my $interval = shift;
  my $chunks_ref = shift; # assoc array
  my $chunknum_ref = shift; # list

  $ckno = sprintf("%.2d", $ckno);
  $start = sprintf("%.10d", $start);
  $end = sprintf("%.10d", $end);
  my $chunkrec =  "$ckno,$start,$end,$numsnps,$interval";
  $$chunks_ref{$ckno} = $chunkrec;
  push(@$chunknum_ref, $ckno);

  return $chunkrec;
}
#-------------------------------------------------------------------
# Remove small chunks - from a chunks dict
#-------------------------------------------------------------------
sub remove_small_chunks {
  my $chunknum_ref = shift;
  my $chunks_ref = shift;

  for my $ckno (@$chunknum_ref) {
    my ($cknum, $start, $end, $numsnps) = split(/,/, $$chunks_ref{$ckno});
    if (($numsnps < 200) and ($cknum > 1)) {
      my $prvckno = sprintf("%.2d", ($ckno - 1));
      my ($prvcknum, $prvstart, $prvend, $prvnumsnps) = split(/,/, $$chunks_ref{$prvckno});
      my $newnumsnps = $prvnumsnps + $numsnps;
      my $newinterval = $end - $prvstart;
      $$chunks_ref{$prvckno} = "$prvcknum,$prvstart,$end,$newnumsnps,$newinterval";
      $$chunks_ref{$ckno} = "$cknum,0,0,0";
    }
  }
  return $chunks_ref;
}
#-------------------------------------------------------------------
# Write chunk index - from a chunks dict
#-------------------------------------------------------------------
sub write_chunk_index {
  my $dir = shift;
  my $chr = shift;
  my $chunks_ref = shift;

  my $outfile = sprintf("%s/chr%.2d_chunk_index.txt", $dir, $chr);
  open(my $ofh, ">", $outfile) || die "Can't open $outfile for output";
  my $fileckno = 0;
  for my $ckno (sort keys %$chunks_ref) {
    my ($cknum, $start, $end, $numsnps, $intrvl) = split(/,/, $$chunks_ref{$ckno});
    if ($numsnps > 0) {
      $fileckno++;
      $fileckno = sprintf("%.2d", $fileckno);
      print $ofh "$fileckno,$start,$end,$numsnps,$intrvl\n";
    }   
  }  
# filehandle goes out of scope ...
  return $fileckno;
}
#-------------------------------------------------------------------
# Exported fn to create a chunk index files
#-------------------------------------------------------------------
sub create_file_chunks {
  my $dir = shift;
  my $file = shift;
  my ($logdir, $conffile, $shdir, $centromere_end_ref, $chnksize, $printonly) = @_;

  my $ckcount = 0;
  my $fullpath = $dir . "/" . $file;
  if ($file =~ /chr(\d+)/) {
    my $chr = $1;
    my $chrkey = sprintf("chr%d", $1); # cytobank file has chr1, chr2, etc NOT chr01, chr02

    my %chunks = ();
    my @chunknum = ();
    print "\n$fullpath, $chr\n";
    open(my $fh, "<", $fullpath) || die "Can't open $fullpath for input";
    my $count = 0;
    my $min = 9999999999;
# $lastpos, $lastcount hold last chunk intervals numbers
    my $lastpos = $min;
    my $lastcount = $count;
# $prevpos, $prevcount hold previous record numbers
    my $prevpos = $min;
    my $prevcount = $count;
#
    my $max = 0;
    my $chnknum = 1;
    my $past_centromere = 0;

    while (<$fh>) {
      $count++;
      chomp;
      my  @record = split();
      $max = $record[3];
      if ($record[3] < $min) {
        print "LT test: $record[3]\n";
        $lastpos = $prevpos = $min = $record[3];
      }
      if ($past_centromere == 0) {
        if ($record[3] >= $$centromere_end_ref{$chrkey}) {
          $past_centromere = 1;
          my $posdiff = $prevpos - $lastpos;
          my $diff = $prevcount - $lastcount;
          my $chunkrec = make_chunk_rec($chnknum, $lastpos, $prevpos, $diff, $posdiff, \%chunks, \@chunknum);
          $lastpos = $record[3];
          $lastcount = $count - 1;
          $chnknum++;
          next;
        }
      }
      if ($record[3] >= ($lastpos + $chnksize)) {
        my $posdiff = $record[3] - $lastpos;
        my $diff = $count - $lastcount;
        if ($diff >= 200) {
          my $chunkrec = make_chunk_rec($chnknum, $lastpos, $record[3], $diff, $posdiff, \%chunks, \@chunknum);
          $lastpos = $record[3] + 1;
          $lastcount = $count;
          $chnknum++;
        }
      }
      $prevpos = $record[3];
      $prevcount = $count;
    } # while loop over file records

    my $interval = $max - ($min - 1);
    my $numchunks = $interval / $chnksize;
    my $posdiff = $max - $lastpos;
    my $diff = $count - $lastcount;
    my $chunkrec = make_chunk_rec($chnknum, $lastpos, $max, $diff, $posdiff, \%chunks, \@chunknum);
    print "\nFILE SNPS: $count\n";

# reprocess to eliminate small chunks
    my $chunks_ref = remove_small_chunks(\@chunknum, \%chunks);
    $ckcount = write_chunk_index($dir, $chr, $chunks_ref);
    my $cmd = sprintf($qsub_tmplt{qsub_impute_cmd}, $logdir, $logdir, $ckcount, $shdir, $conffile, $chr);
    if ($printonly eq "false") {
      system($cmd);
    }
    print "CMD: $cmd\n";
  }
  else {
    print "Filename: $fullpath doess not contain chromosome #\n";
  }
  print "\nFILE CHUNKS: $ckcount\n";
}



1; # end of module
