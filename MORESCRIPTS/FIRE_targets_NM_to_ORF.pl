#
# a basic script that converts a .target file to ORF names from NM
#


#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

my $ta = Table->new;
$ta->loadFile("/Users/olivier/PROGRAMS/CANCERGENES/SANGER/refLink.simplified.human");
my $h_ref = $ta->getIndexKV(0,1);

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
foreach my $r (@$a_ref) {
  my $n = shift @$r;
  my @a = ();
  foreach my $s (@$r) {
    if (defined($h_ref->{$s})) {
      push @a, $h_ref->{$s} if (!Sets::in_array($h_ref->{$s}, @a));
    } else {
      push @a, $s;
    }
  }
  print "$n\t" . join("\t", @a) . "\n";
}

