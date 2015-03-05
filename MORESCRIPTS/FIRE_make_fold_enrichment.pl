#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Getopt::Long;
use Sets;
use Table;
use strict;

my $motiffile = undef;
my $expfile   = undef;

GetOptions("expfile=s"   => \$expfile,
           "motiffile=s" => \$motiffile);

my $ff     = Sets::filename($expfile);

if (!defined($motiffile)) {
  $motiffile = "$expfile\_FIRE/DNA/$ff.summary";
}

my $ta = Table->new;
$ta->loadFile($motiffile);
my $a_ref_mot = $ta->getColumn(0);
my $h_ref_name = $ta->getIndexKV(0,1);

# load signif
my $signif = "$expfile\_FIRE/DNA/$ff.signif.motifs.rep";
$ta->loadFile($signif);
my $a_ref_sig = $ta->getArray();
my %H = ();
foreach my $r (@$a_ref_sig) {
  $H{$r->[0]}{$r->[1]} = $r->[4];
}

# load pvalues
my $matrix = "$expfile\_FIRE/DNA/$ff.matrix";
$ta->loadFile($matrix);
my $a_ref_mat = $ta->getArray();
my %M = ();
foreach my $r (@$a_ref_mat) {
  $M{$r->[0]} = $r->[2];
}

print "Motif\tFold-enrichment\tp-value\tMotif name\n";
foreach my $m (@$a_ref_mot) {
  if (!defined($H{$m})) {
    die "$m not defined in $signif ?.\n";
  }
  my $fo = sprintf("%3.2f", $H{$m}{1} / $H{$m}{0});
  my $lp = $M{$m};
  if ($lp > 0) {
    $lp = -1 * $lp;
  }
  my $pv = Sets::round_pvalue_up( exp( $lp * log(10) ) );
  my $na = $h_ref_name->{$m};
  print "$m\t$fo\t$pv\t$na\n";
  
}
