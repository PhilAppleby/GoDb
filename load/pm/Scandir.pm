package Scandir;

use strict;
use warnings;

sub new() {
  my $class = shift;

  return bless {}, $class;
}

sub scan_dir {
  my $self = shift;
  my $dir = shift;
  my $suffix = shift;
  my $oper = shift;
  my $params = shift; # can be empty

  opendir(my $dh, $dir) or die "Can't open dir [$dir]";

  my @files= grep{/${suffix}$/} readdir $dh;

  for my $file (@files) {
    my $fullpath = $dir . "/" . $file;
    if (($file eq "..") or ($file eq ".")) {
      next;
    }
    if (-l $fullpath) {
      next;
    }
    &$oper($dir, $file, @$params);
  }
}
1; # end of module
