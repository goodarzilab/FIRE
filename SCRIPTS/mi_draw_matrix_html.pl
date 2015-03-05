use lib "$ENV{FIREDIR}/SCRIPTS";
use lib "$ENV{FIREDIR}/SCRIPTS/PostScript-Simple-0.07/lib";

use Table;
use Sets;
use Getopt::Long;
use PostScript::Simple;
use strict;
use Data::Dumper;

my $scriptdir      = "$ENV{FIREDIR}/SCRIPTS";
my $progdir        = "$ENV{FIREDIR}/PROGRAMS";

my $matfile        = undef;
my $densityfile    = undef;
my $motiffile      = undef;
my $quantized      = 0;
my $clustfile      = undef;
my $rna            = 0;
my $sortcolsbysize = 0;
my $sortrowsbymax  = 0;
my $columnsfile    = undef;
my $distfile       = undef;
my $limit          = undef;
my $ps2pdf         = 0;
my $w              = undef;
my $h              = undef;
my $every          = 1;
my $consfile       = undef;
my $gofile         = undef;
my $colmap         = undef;
my $dopng          = 0;
my $rootdir        = '.';;
my $showstatlabels = 1;
my $showscale      = 1;
my $ybase          = undef;
my $redoweblogo    = 1;
my $motifnames     = undef;
my $outeps         = undef;
my $outfile        = undef;
my $lp_t_draw      = 10;
my $expfile        = undef;

if (@ARGV == 0) {
  die "Usage: perl mi_draw_matrix.pl  --matfile=FILE --summaryfile=FILE --clustfile=FILE --columnsfile=FILE --gofile=FILE --ps2pdf=INT --every=1 --quantized=INT\n";
}

GetOptions (
	    'expfile=s'        => \$expfile,
	    'outfile=s'        => \$outfile,
	    'matfile=s'        => \$matfile,
	    'densityfile=s'    => \$densityfile,
	    'outeps=s'         => \$outeps,
	    'summaryfile=s'    => \$motiffile,
	    'clustfile=s'      => \$clustfile,
	    'motifnames=s'     => \$motifnames,
	    'columnsfile=s'    => \$columnsfile,
	    'distfile=s'       => \$distfile,
	    'consfile=s'       => \$consfile,
	    'sortcolsbysize=s' => \$sortcolsbysize,
	    'sortrowsbymax=s'  => \$sortrowsbymax,
	    'showstatlabels=s' => \$showstatlabels,
	    'showscale=s'      => \$showscale,
	    'w=s'              => \$w,
	    'h=s'              => \$h,
	    'gofile=s'         => \$gofile,
	    'every=s'          => \$every,
	    'limit=s'          => \$limit,
	    'ps2pdf=s'         => \$ps2pdf,
	    'colmap=s'         => \$colmap,
	    'dopng=s'          => \$dopng,
	    'ybase=s'          => \$ybase,
	    'rootdir=s'        => \$rootdir,
	    'redoweblogo=s'    => \$redoweblogo,
	    'quantized=s'      => \$quantized,
	    'lp_t_draw=s'      => \$lp_t_draw);


if (! -e $clustfile) {
  print "$clustfile does not exist, ignoring.\n";
  $clustfile = undef;
}

my $ta = Table->new;


#
#  load GO (for clusters)
#
my %GO = ();
if (defined($gofile) && (-e $gofile)) {
  $ta->loadFile($gofile);
  %GO = %{ $ta->getIndexKV(0,1) };
  print "GO file loaded.\n";
}


my @MOTIFS = ();
my $a_ref_clust = undef;


#
#  read in the matrix file
#
$ta->loadFile($matfile);
my $a_ref_M      = $ta->getArray();

if (!defined($a_ref_clust)) {
  @MOTIFS = @{ $ta->getColumn(0) };
  shift @MOTIFS;
}
my $a_ref_H = shift @$a_ref_M; shift @$a_ref_H;

my %MATRIX       = ();
foreach my $r (@$a_ref_M) {
  my $m = shift @$r;
  $MATRIX{ $m } = $r;
}



#
#  read in the density file
#
my %DENSITIES       = ();
if (defined($densityfile)) {
  my $a_ref_D = undef;
  $ta->loadFile($densityfile);
  $a_ref_D      = $ta->getArray();
  shift @$a_ref_D;
  foreach my $r (@$a_ref_D) {
    my $m = shift @$r;
    $DENSITIES{ $m } = $r;
  }
}


#
# read in motifs names
#
$ta->loadFile($motifnames);
my $h_ref_names   = $ta->getIndexKV(0,1);


#
# read in the clusters
#
if (defined($clustfile)) {
  $ta->loadFile($clustfile);
  $a_ref_clust = $ta->getColumn(1);
  my $col      = $ta->getColumn(0);
  @MOTIFS      = @$col;
}


#
#  read in the summary file
#
$ta->loadFile($motiffile);
my $a_ref_mo = $ta->getArray();
my %STAT         = ();
my %MOTIF_NUMBER = ();

my %motif_name = ();

my $cnt = 0;
foreach my $r (@$a_ref_mo) {

  my %a_tmp = ( "RNA"    => $r->[1],
		"COPIES" => $r->[2],
		"MI"     => $r->[3],
		"RANK"   => $r->[4],
		"Z"      => $r->[5],
		"R"      => $r->[6],
		"S"      => $r->[7],
		"SEED"   => $r->[8],
		"DIST"   => $r->[9],
		"ORIE"   => $r->[10],
		"CONS"   => $r->[11],
		"NAME"   => $h_ref_names->{ $r->[0] });

  $STAT{ $r->[0] }         = \%a_tmp;
  $MOTIF_NUMBER{ $r->[0] } = $cnt ++;

  $motif_name{$r->[0]} = $r->[8];

  if (($r->[8] eq '') || ($r->[8] eq '0')) {
   $motif_name{$r->[0]} = $r->[0];
  }
}



$ta->loadFile($columnsfile);
my $a_ref_cols  = $ta->getColumn(0);

#
# load color map
#
my $A_REF_COLMAP = undef;
if (defined($colmap)) {

  $ta->setDelim(" ");
  $ta->loadFile($colmap);
  $A_REF_COLMAP = $ta->getArray();
  $ta->setDelim('\t');
}


#$a_ref_M = sort2DArrayByMax($a_ref_M);



print "Now doing the graphical display.\n";

#
#  START DRAWING
#
#

die "No motifs in \@MOTIFS ...\n" if (@MOTIFS == 0);

my $xsize = 600 ;
my $xbase = 120 ;
my $ybase = 180 ;
my $w = int( 0.5 + $xsize / @$a_ref_H ) ;
if ($w<6){
  $w = 6 ;
  $xsize = $w*@$a_ref_H + 1200 ;
}else{
  $xsize = $w*@$a_ref_H + 1200 ;
}

my $h = 50 ;

my $d = "$motiffile" . "_OUT";
my $fn      = Sets::filename($motiffile);
my $dname = $fn . "_OUT";

if (! -e $d) {
  mkdir $d;
}

my $fnOne = Sets::filename($matfile);

my $min  = undef;
my $max  = undef;
if (defined($densityfile)) {
  $min = 0.0;
  $max = 1.0;
} else {
  $min = -$lp_t_draw;
  $max =  $lp_t_draw;
}

my $fsize = $h/4;

open O, "> $outfile";
print O (
	"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\"
   \"http://www.w3.org/TR/html4/strict.dtd\">\n"
);
print O ("<head>\n");
print O ("<title>FIRE results</title>\n");
print O ("<style type=\"text/css\">");
print O ("
* {
padding: 0;
margin: 0;
}
body {
font-family: Helvetica, Arial, sans-serif;
}
h1, p, p.head {
width:960px;
padding: 10px;
}
td.grid {
width: $w\px;
height: $h\px;
font-size:$fsize\px;
}
td.gridA {
font-size:$fsize\px;
font-family: Helvetica, Arial, sans-serif;
text-align: center;
padding: 0px 2px;
}
td.gridheader {
height: auto;
width: $w\px;
border: 1px solid #FFF;
text-align: center;
background-color: #000;
}
td.gridheaderA {
height: auto;
width: $w\px;
border: 1px solid #FFF;
text-align: center;
font-size: 10px;
}
td.gridheaderB {
height: auto;
width: $w\px;
border: 1px solid #FFF;
text-align: center;
background-color: #FFF;
}
td.gridlabel {
white-space: nowrap;
height: $w\px;
padding-left: 5px;
}
table.legend {
margin-left: 10px;
top: 0px;
left: 0px;
position: relative;
}
table.legend td {
text-align: center;
}
table.legend tr {
font-size:0px;
height: 2px;
}
table.main {
position: absolute;
top: $ybase\px;
left: $xbase\px;
padding: 10px;
}
");
print O ("</style>\n");
print O ("</head>\n");
print O ("<body>\n\n");

print O ("<h1><abbr title=\"Finding Informative Regulatory Elements\">FIRE</abbr> results</h1>\n");

print O ("<p class=\"head\">These results are also available as <a href=\"./$fn.pdf\"><abbr title=\"Portable Document Format\">PDF</abbr></a> and <a href=\"./$fn.eps\"><abbr title=\"Encapsulated PostScript\">EPS</abbr></a> documents.</p>\n");

$fnOne =~ s/matrix/fullmimatrix/;

print O ("<p>The <a href=\"$fnOne.html\">motif interaction matrix</a> may also be informative.</p>");

print O ("<p>Depending on your display resolution, scrolling or zooming may be necessary.</p>\n");


############ Scale ############
my $scalex = $xbase/3 ;
my $hscale = $h*2/3 ;
my $W = $hscale ;

print O ("<table class=\"legend\" width=\"$hscale\px\" cellpadding=\"0\" cellspacing=\"0\">\n");
my $p = new PostScript::Simple(xsize => $W/2,
			       ysize => $h+10,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("Helvetica", 10);
$p->text( { rotate => 90 }, 10, 5, "over-rep") ;
$p->output("$d/over.eps") ;
system("convert -density 100 $d/over.eps $d/over.png");
print O "<tr><td style=\"border:1px #FFFFFF solid;\">" ;
print O "<img src=\"./$dname/over.png\" alt=\"over-representation\" /></td></tr>\n" ;

print O "<tr><td style=\"font-size: 12px; height: 14px; \">$max</td></tr>\n" ;

my $t = $max;
my $res = 70 ;
for (my $i=0; $i<=$res; $i++) {
  my @col = () ;
  if (!defined($colmap)) {
    if ($i>$res/2){
      @col = Sets::interp_general( $t, [0, 0, 0], [255, 0, 0], $min, 0);
    }else{
      @col = Sets::interp_general( $t, [255, 0, 0], [255, 255, 0], 0, $max);
    }
  }else{
    @col = Sets::interp_from_matlab_colormap( $t, $A_REF_COLMAP, $min, $max);
  }
  if ($i==$res/2){
    print O "<tr>\n" ;
    print O "<td style=\"font-size: 12px; height: 14px; \">0</td>\n" ;
    print O "</tr>\n" ;
  }else{
    my $color = &RGBDecToHex($col[0], $col[1], $col[2]);
    print O "<tr>\n" ;
    print O "<td style=\"background-color:$color\">&nbsp;</td>\n" ;
    print O "</tr>\n" ;
    $t -= ($max - $min) / $res;
  }
}

print O "<tr><td style=\"font-size: 12px; height: 14px; \">$min</td></tr>\n" ;
my $p = new PostScript::Simple(xsize => $W/2,
			       ysize => $h+10,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("Courier", 10);
$p->text( { rotate => 90 }, 10, 5, "under-rep") ;
$p->output("$d/under.eps") ;
system("convert -density 100 $d/under.eps $d/under.png");
print O "<tr><td>" ;
print O "<img src=\"./$dname/under.png\" alt=\"under-representation\"/></td></tr>\n" ;

print O ("\n</table>\n");

my $min_i1 =  1000000;
my $max_i2 = -1000000;
my @bins   = ();
if ($quantized == 0){
  foreach my $c (@$a_ref_H) {
    my ($i1, $i2) = $c =~ /\[(.+?);(.+?)\]/;
    $min_i1 = $i1 if ($i1 < $min_i1);
    $max_i2 = $i2 if ($i2 > $max_i2);
    my @a_tmp = ($i1, $i2); push @bins, \@a_tmp;
  }
  my $l = $xbase-30 ;
  my $t = $ybase+10;
  $l = sprintf("%2.2f", $l) ;
  print O "<div style=\"font-size:9px ; top:$t\px ; left:$l\px ; position:absolute ;\">$max_i2</div>" ;
  my $t = $ybase+$h*1.5-8 ;
  $l = sprintf("%2.2f", $l) ;
  print O "<div style=\"font-size:9px ; top:$t\px ; left:$l\px ; position:absolute ;\">$min_i1</div>" ;
}

print O ("<table width=\"$xsize\px\" cellpadding=\"0\" cellspacing=\"0\" class=\"main\">\n");

if ($quantized == 1){
  print O "<tr>\n" ;
  for (my $j=0; $j<@$a_ref_H; $j++) {
    if ($GO{ $a_ref_cols->[$j] } ne "") {
      my $s = "Enriched in: ".$GO{ $a_ref_cols->[$j]} ;
      print O "<td class=\"gridheaderA\" onclick=\"alert('$s')\" title=\"$s\">" ;
    }else{
      print O "<td class=\"gridheaderA\" style=\"background-color: #FFFFFF;\">" ;
    }

    my $c = $a_ref_H->[$j] ;
    my $W = 20 ;

    my $p = new PostScript::Simple(xsize => $w,
				   ysize => $W+10,
				   colour    => 1,
				   eps       => 1,
				   units     => "pt");
    $p->setcolour("black");
    $p->setfont("Courier", 6);
    $p->text( { rotate => 90 }, $w*2/3, 3, $c) ;
    $p->output("$d/$c.eps") ;
    system("convert -density 100 $d/$c.eps $d/$c.png");
    print O "<img width=\"$w\" src=\"./$dname/$c.png\" alt=\"$c\" /></td>\n" ;
  }
}else{
  my $th = $max_i2 - $min_i1;
	$th = 1e-3 if ($th < 1e-3) ;

  my $H = $h*1.5 ;
  print O "<tr style=\"height:$H\px\">\n" ;
  $H-=3 ;
  my $j = 0;
  foreach my $c (@$a_ref_H) {
    my $h1 = $H * ($bins[$j]->[0] - $min_i1 ) / $th ;
    my $h2 = $H * ($bins[$j]->[1] - $min_i1 ) / $th ;

    my $s = "lower bound = ".$bins[$j]->[0]." and upper bound = ".$bins[$j]->[1] ;
    if ($GO{ $a_ref_cols->[$j] } ne "") {
      $s.=" Enriched in: ".$GO{ $a_ref_cols->[$j]} ;
    }

    print O "<td class=\"gridheader\" onclick=\"alert('$s')\" title=\"$s\">" ;
    my $dw = $w ;
    my $dh = $h2-$h1 ;
    $dh = 1 if ($dh<1) ;
    my $y = ($h1+$h2)/2-($H-2)/2 ;
    print "$y\t$H\t$h1\n" ;
#    print O "<div style=\"width:$dw\px; height:$dh\px; background-color:#F00 ; left:0\px ; bottom:$y\px; position:relative\">" ;
    print O "<div style=\"width:100\%; height:$dh\px; background-color:#F00 ; left:0\px ; bottom:$y\px; position:relative\">" ;

    print O "</div>\n" ;
    print O "</td>\n" ;
    $j++ ;
  }
}

my $W = 80 ;

my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("courier", 10);
$p->text( { rotate => 45 }, 10, 5, "motif") ;
$p->output("$d/motif.eps") ;
system("convert -density 100 $d/motif.eps $d/motif.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/motif.png\" style=\"top:5px ; left:25px; position:relative;\" alt=\"motif\"></td>\n" ;

my $W = 50 ;
my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("courier", 10);
$p->text( { rotate => 45 }, 10, 5, "location") ;
$p->output("$d/location.eps") ;
system("convert -density 100 $d/location.eps $d/location.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/location.png\" style=\"top:5px ; left:5px; position:relative;\" alt=\"location\"></td>\n" ;

my $W = 80 ;
my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("courier", 10);
$p->text( { rotate => 45 }, 10, 5, "MI") ;
$p->output("$d/mi.eps") ;
system("convert -density 100 $d/mi.eps $d/mi.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/mi.png\" style=\"top:5px ; left:30px; position:relative;\" alt=\"mutual information\"></td>\n" ;

my $W = 80 ;
my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("courier", 10);
$p->text( { rotate => 45 }, 10, 5, "z-score") ;
$p->output("$d/z.eps") ;
system("convert -density 100 $d/z.eps $d/z.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/z.png\" style=\"top:5px ; left:15px; position:relative;\" alt=\"z-score\"></td>\n" ;

my $W = 80 ;
my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("courier", 10);
$p->text( { rotate => 45 }, 10, 0, "robustness") ;
$p->output("$d/robustness.eps") ;
system("convert -density 100 $d/robustness.eps $d/robustness.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/robustness.png\" style=\"top:5px ; left:15px; position:relative;\" alt=\"robustness\"></td>\n" ;

my $W = 60 ;
my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("courier", 10);
$p->text( { rotate => 45 }, 10, 0, "position") ;
$p->text( { rotate => 45 }, 38, 6, "bias") ;
$p->output("$d/pos.eps") ;
system("convert -density 100 $d/pos.eps $d/pos.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/pos.png\" style=\"top:5px ; left:3px; position:relative;\" alt=\"position bias\"></td>\n" ;

my $W = 60 ;
my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("courier", 10);
$p->text( { rotate => 45 }, 10, 0, "orientation") ;
$p->text( { rotate => 45 }, 38, 6, "bias") ;
$p->output("$d/ori.eps") ;
system("convert -density 100 $d/ori.eps $d/ori.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/ori.png\" style=\"top:5px ; left:3px; position:relative;\" alt=\"orientation bias\"></td>\n" ;

my $W = 60 ;
my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("courier", 9);
$p->text( { rotate => 45 }, 10, 0, "conservation") ;
$p->text( { rotate => 45 }, 38, 6, "index") ;
$p->output("$d/cons.eps") ;
system("convert -density 100 $d/cons.eps $d/cons.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/cons.png\" style=\"top:5px ; left:10px; position:relative;\" alt=\"conservation index\"></td>\n" ;

my $W = 80 ;
my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("courier", 10);
$p->text( { rotate => 45 }, 10, 4, "seed") ;
$p->output("$d/seed.eps") ;
system("convert -density 100 $d/seed.eps $d/seed.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/seed.png\" style=\"top:5px ; left:30px; position:relative;\" alt=\"seed\"></td>\n" ;

my $W = 500 ;
$W = 80;
my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("courier", 10);
$p->text( { rotate => 45 }, 5, 4, "motif") ;
$p->text( { rotate => 45 }, 23, 5, "name") ;
$p->output("$d/name.eps") ;
system("convert -density 100 $d/name.eps $d/name.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/name.png\" style=\"top:5px ; left:20px; position:relative;\" alt=\"motif name\"></td>\n" ;

$W = 80;
my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("courier", 10);
$p->text( { rotate => 45 }, 5, 4, "protein") ;
$p->text( { rotate => 45 }, 23, 5, "array") ;
$p->output("$d/parray.eps") ;
system("convert -density 100 $d/parray.eps $d/parray.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/parray.png\" style=\"top:5px ; left:20px; position:relative;\" alt=\"protein array\"></td>\n" ;

print O "</tr>\n" ;

my $cntm = 0;
my $NN   = @$a_ref_H;
my $i    = 0;


foreach my $re (@MOTIFS) {
  print O "<tr class=\"motif\">\n" ;

  my $myre = $re;
  if ($STAT{$re}->{RNA} == 1) {
    $myre =~ s/T/U/g;
  }

  my $r = undef;
  if (!defined($densityfile)) {
    $r = $MATRIX{ $re };  # get the data (row)
  } else {
    $r = $DENSITIES{ $re };
  }

  print "Processing $re ... ";

  #
  #  MOTIF NAME
  #
  my $myii = $MOTIF_NUMBER{ $re };

  print "Outputing motif $myii.png ... \n";

  my $mo = Sets::myre2wm($myre);
  open OUT, ">$d/$myii.txt" or die "cannot open $d/$myii.txt\n";
  print OUT $mo;
  close OUT;

  system("$scriptdir/weblogo/seqlogo -f $d/$myii.txt -F PNG  -a -c -M -n -Y -w 5 -h 3 > $d/$myii.png");

  my $j = 0;
  foreach my $c (@$r) {

    my @col = ();
    if (!defined($colmap)) {
      @col = Sets::interp_general( $c, [0, 0, 204], [255, 255, 0], $min, $max);
    } else {
      @col = Sets::interp_from_matlab_colormap( $c, $A_REF_COLMAP, $min, $max);
    }

    my $color = &RGBDecToHex($col[0], $col[1], $col[2]);
    my $lp = - Sets::log10( 0.05 / $NN );
    my $pv ;
    if ($c>0){
      $pv = 10**(-1*$c) ;
    }else{
      $pv = 10**$c ;
    }
    $pv = sprintf("p = %1.4e", $pv) ;

    if (abs($c) < $lp) {
      print O "<td class=\"grid\" style=\"background-color:$color\" onclick=\"alert('$pv')\" title=\"$pv\">&nbsp;</td>\n" ;
    }else{
      if ($c>0){
	print O "<td class=\"grid\" style=\"background-color: $color; border:1px #FF0000 solid\" onclick=\"alert('over-representation: $pv')\" title=\"over-representation: $pv\">&nbsp;</td>\n" ;
      }else{
	print O "<td class=\"grid\" style=\"background-color: $color; border:1px #0000FF solid\" onclick=\"alert('under-representation: $pv')\" title=\"under-representation: $pv\">&nbsp;</td>\n" ;
      }
    }
    $j ++;
  }

  my $myseed = $STAT{$re}->{SEED};
  if ($STAT{$re}->{RNA} == 1) {
    $myseed =~ s/T/U/g;
  }

  print O "<td class=\"gridA\">" ;
  print O "<a href=\"./$dname/$motif_name{$re}\.pdf\"><img src=\"./$dname/$myii.png\" height=\"$h\" style=\"top:5px ; position:relative; border: 0px;\" alt=\"$myre\" title=\"$myre\"></a></td>\n" ;
  print O "<td class=\"gridA\">" ;
  print O ($STAT{$re}->{RNA}==1?"3'UTR":"5'")."</td>\n" ;
  print O "<td class=\"gridA\">" ;
  print O sprintf("%4.3f", $STAT{$re}->{"MI"})."</td>\n" ;
  print O "<td class=\"gridA\">" ;
  print O sprintf("%3.1f", $STAT{$re}->{"Z"})."</td>\n" ;
  print O "<td class=\"gridA\">" ;
  print O $STAT{$re}->{R}."</td>\n" ;
  print O "<td class=\"gridA\">" ;
  print O  ($STAT{$re}->{DIST}==1?"Y":"-")."</td>\n" ;
  print O "<td class=\"gridA\">" ;
  if ($STAT{$re}->{ORIE} == 0) {
    print O "-</td>\n" ;
  }else{
    if ($STAT{$re}->{ORIE} == 1){
      print O "&#8594;</td>\n" ;
    }else{
      print O "&#8592;</td>\n" ;
    }
  }
  print O "<td class=\"gridA\">" ;
  print O (defined($STAT{$re}->{CONS})?$STAT{$re}->{CONS}:'-')."</td>\n" ;

  my $gaps ;
  my $format=undef;
  if ($myseed =~ /\./)
  {
      $myseed =~ /(\.+)/;
      $gaps=$1;
      $gaps=length($gaps);
      $format = ".($gaps)";
  }
  $myseed =~ s/\.+/$format/ if (defined $format and $gaps>2);
  print O "<td class=\"gridA\"><span class=\"seed\">" ;
  print O $myseed."</span></td>\n" ;

  my $na = $STAT{$re}->{NAME};  # dirty fix

  my @curmotifs = split(/\;/, $na);
  my $ht_curmotifs = $curmotifs[0];
  my @ht_curmotiflist = split(/\,/, $ht_curmotifs);
  my @ht_names = ();
  my @ht_pa = ();

  foreach my $curmotif (@ht_curmotiflist) {
  	$curmotif =~ s/^M\d{5}\_//;
  	$curmotif =~ s/^J\_...?\d+\_//;
	$curmotif =~ s/^P\_...?\d+\_//;
  	$curmotif =~ s/\.txt$//;

  	if ($curmotif =~/^PA\_(\d+)$/) {
  		push(@ht_pa, $curmotif);
  	} else {
  		push(@ht_names, $curmotif);
  	}
  }

  if (scalar(@ht_pa) == 0) {
  	push(@ht_pa, '-');
  }

  if (scalar(@ht_names) == 0) {
  	push(@ht_names, '-');
  }

# for normal hits at 0.8 threshold
  print O "<td class=\"gridA\" style=\"white-space:nowrap; text-align: left;\">" ;
  my $firstitemname = 0;
  #foreach my $item (@ht_names) {
  for(my $i = 0; $i < 3; $i++) {
  	  		if (!defined $ht_names[$i]) {
  			last;
  		}
  	my $item = $ht_names[$i];
  	if ($firstitemname == 0) {
  		print O $item;
  	} else {
  		print O ', ', $item;
  	}
  	$firstitemname = 1;
  }
  if (scalar(@ht_names) > 3) {
  	print O ' ...';
  }

  print O "</td>\n" ;
#  print O "</tr>\n" ;

# for protein array hits

  print O "<td class=\"gridA\" style=\"white-space:nowrap; text-align: left;\">" ;
    my $firstitemname = 0;
  #foreach my $item (@ht_pa) {
  	for(my $i = 0; $i < 3; $i++) {
  		if (!defined $ht_pa[$i]) {
  			last;
  		}
  	my $item = $ht_pa[$i];
if ($item =~ m/^PA\_(\d+)$/) {
  	$item = '<a href="http://bioinfo.wilmer.jhu.edu/PDI/motif.cgi?id='.$1.'"><!--<abbr title="Protein Array">PA</abbr> match (-->ID: '.$1.'<!--)--></a>';
  }
  	if ($firstitemname == 0) {
  		print O $item;
  	} else {
  		print O ', ', $item;
  	}
  	$firstitemname = 1;
  }

  if (scalar(@ht_pa) > 3) {
  	print O ' ...';
  }

  print O "</td>\n" ;
  print O "</tr>\n" ;

  if ($a_ref_clust->[$i] ne $a_ref_clust->[$i+1]) {
    print O "<tr style=\"height:5\px\">\n" ;
    print O "<td style=\"font-size:0px\"></td>\n" ;
    print O "</tr>\n" ;
  }

  $i ++;
}

print O ("\n</table>\n");
print O ("\n</body>\n");
print O ("</html>\n");

sub RGBDecToHex {
  my ($red, $green, $blue) = @_;
  return sprintf("#%02lx%02lx%02lx", $red, $green, $blue);
}
