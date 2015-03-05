use lib "$ENV{FIREDIR}/SCRIPTS";
use Table;
use Sets;
use Fire;

use lib qw(PostScript-Simple-0.07/lib);
use PostScript::Simple;

use strict;


if (@ARGV == 0) {
  die "Usage: perl motif_table__to_eps.pl table\n";

}

#
# load clusters to get order
#
my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $head  = shift @$a_ref;

#
#  START DRAWING
#
#

my $xbase = 50;
my $ybase = 50;
my $h     = 30;
my $xsize = 595;

my $nm    = @$a_ref;
my $ysize = $nm * $h + $ybase + 10;

my $p = new PostScript::Simple(xsize     => $xsize,
			       ysize     => $ysize,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");

$p->setlinewidth(0.5);
$p->setcolour("black");
$p->setfont("Arial", 12);

my @offset = (25, 105, 190, 260);

$p->text({align => "center"}, $xbase + $offset[0],  $ysize - ( $ybase - 15),  $head->[0] );
$p->text({align => "center"}, $xbase + $offset[1],  $ysize - ( $ybase - 15 ), $head->[1] );
$p->text({align => "center"}, $xbase + $offset[2],  $ysize - ( $ybase - 15 ), $head->[2] );
$p->text({align => "center"}, $xbase + $offset[3],  $ysize - ( $ybase - 15 ), $head->[3] );

#$p->line( $xbase, $ysize - ( $ybase - 10 ), $xbase + 400, $ysize - ( $ybase - 10 ));

my $i = 0;
foreach my $r (@$a_ref) {
  print join("\t", @$r) . "\n";
  my $e1 = &get_eps_logo_object($r->[0], 0, $i);
  my $eh = $e1->height;
  my $ew = $e1->width;
  my $ew_new = int(0.5 + $ew * $h / $eh); #print "e=$ew_new\n";
  $e1->scale($h / $eh);
  
  $p->_add_eps($e1, $xbase,  $ysize - ( $ybase + $i*$h+$h ) ); 

  $p->text({align => "center"}, $xbase + $offset[1], $ysize - ( $ybase + $i*$h+$h/2 ), $r->[1] );
  $p->text({align => "center"}, $xbase + $offset[2], $ysize - ( $ybase + $i*$h+$h/2 ), $r->[2] );
  $p->text({align => "center"}, $xbase + $offset[3], $ysize - ( $ybase + $i*$h+$h/2 ), $r->[3] );

  $i ++;
}


$p->output("$ARGV[0].eps");
system("ps2pdf  -dEPSCrop -dAutoRotatePages=/None $ARGV[0].eps $ARGV[0].pdf");



sub get_eps_logo_object {
  
  my ($re, $rna, $cnt) = @_;
  
  #
  # make motif logo
  #
  my $mo = Sets::myre2wm($re);
  if ($rna == 1) {
    $mo =~ s/T/U/g;
  }
  
  my $d           = "TMP";
  
  if (!defined($cnt)) {
    $cnt         = 0;
  }

  mkdir $d if (! -e $d);   # to be cleaned

  open OUT, ">$d/$cnt.txt" or die "cannot open $d/$cnt.txt\n";
  print OUT $mo;
  close OUT;
  system("$ENV{FIREDIR}/SCRIPTS/weblogo/seqlogo -f $d/$cnt.txt -F EPS  -a -c -M -n -Y -w 5 -h 3 > $d/$cnt.eps");

  #
  # START: integrate Motif LOGO
  #  
  my $e  = new PostScript::Simple::EPS(file => "$d/$cnt.eps");
 
  return $e;  
}

