use strict;
use warnings;

use Scandir;
use Getopt::Long;

my %options = (
  bindir => q{},
  dirname => q{},
  cfgfile => q{},
  printonly => q{},
);

$options{printonly} = "false";

GetOptions(
  'tmplt=s' => \$options{tmplt},
  'fsfx=s' => \$options{fsfx},
  'bindir=s' => \$options{bindir},
  'dirname=s' => \$options{dirname},
  'cfgfile=s' => \$options{cfgfile},
  'printonly=s' => \$options{printonly},
);

sub load_markers {
  my $dir = shift;
  my $file = shift;
  my ($script_tmplt, $bindir, $cfgfile, $printonly) = @_;
  print $cfgfile . "\n";
  if ($file =~ /\.*chr(\d+)/) {
    print $script_tmplt . "\n";
    print $1 . "\n";
    my $cmd = sprintf $script_tmplt, $bindir, $cfgfile, $1;
    if ($printonly ne "false") {
      print $cmd . "\n";
    }
    else {
      print "EXECUTE: " . $cmd . "\n";
      system($cmd);
    }
  }
}

sub process_dir {
  my $dir = shift;
  my $script_tmplt = shift;
  my $file_sfx = shift;
  my $bindir = shift;
  my $cfgfile = shift;
  my $printonly = shift;

  my $scanner = Scandir->new();

  my @params = ($script_tmplt, $bindir, $cfgfile, $printonly);

  $scanner->scan_dir($dir, $file_sfx, \&load_markers, \@params);
}

my $script_tmplt = "%s/load_markers_from_single_info_file.sh %s %s";
process_dir($options{dirname}, $options{tmplt}, $options{fsfx}, $options{bindir}, $options{cfgfile}, $options{printonly})


