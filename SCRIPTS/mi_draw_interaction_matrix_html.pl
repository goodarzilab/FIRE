use lib "$ENV{FIREDIR}/SCRIPTS";
use lib "$ENV{FIREDIR}/SCRIPTS/PostScript-Simple-0.07/lib";

use PostScript::Simple;

# use Data::Dumper;

use Table;
use Sets;

use strict;


my $scriptdir = "$ENV{FIREDIR}/SCRIPTS";
my $progdir   = "$ENV{FIREDIR}/PROGRAMS";

my $colmap    = "$scriptdir/HEATMAPS/cmap_hot.txt";
my $showscale = 1;
my $ps2pdf    = 1;
my $rootdir   = ".";

if (Sets::exist_parameter(\@ARGV, "-colmap") == 1) {
  $colmap   = Sets::get_parameter(\@ARGV, "-colmap");
}
     
my $summaryfile   = Sets::get_parameter(\@ARGV, "-summaryfile");
Sets::defined_exists_or_die($summaryfile, "Please define -summaryfile.\n");

my $matrixfile   = Sets::get_parameter(\@ARGV, "-matrixfile");
Sets::defined_exists_or_die($matrixfile, "Please define -matrixfile.\n");

my $resmatrixfile   = Sets::get_parameter(\@ARGV, "-resmatrixfile");
Sets::defined_exists_or_die($resmatrixfile, "Please define -resmatrixfile.\n");

my $clustfile = undef;
if (Sets::exist_parameter(\@ARGV, "-clustfile") == 1) {
  $clustfile   = Sets::get_parameter(\@ARGV, "-clustfile");
}

my $orderfile = undef;
if (Sets::exist_parameter(\@ARGV, "-orderfile") == 1) {
  $orderfile   = Sets::get_parameter(\@ARGV, "-orderfile");
}

my $outfile = undef;
if (Sets::exist_parameter(\@ARGV, "-outfile") == 1) {
  $outfile   = Sets::get_parameter(\@ARGV, "-outfile");
}

my $motifdir    = Sets::get_parameter(\@ARGV, "-motifdir");
Sets::defined_exists_or_die($motifdir, "Please define -motifdir.\n");

my $fn      = Sets::filename($matrixfile);
my $dname      = Sets::filename($summaryfile);
$dname .= '_OUT';

my @MOTIFS = ();

my $d = $motifdir ;
if (! -e $d) {
  mkdir $d;
}

my $ta = Table->new;

#
#  read in the summary file
#
print "Read in summary file ... ";
$ta->loadFile($summaryfile);
my $a_ref_mo = $ta->getArray();
my %STAT      = ();
my %MOTIF_NUMBER = ();
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
		"CONS"   => $r->[11]);

  $STAT{ $r->[0] } = \%a_tmp;   
  $MOTIF_NUMBER{ $r->[0] } = $cnt ++;
  push @MOTIFS, $r->[0];
}
print "Done.\n";


#
# read in the clusters
#
my $a_ref_clust = undef;
my $h_ref_clust = undef;
if (defined($clustfile)) {
  $ta->loadFile($clustfile);
  $a_ref_clust = $ta->getColumn(1);
  $h_ref_clust = $ta->getIndexKV(0,1);
  my $col         = $ta->getColumn(0);
  @MOTIFS      = @$col;  # over-write
}


if (defined($orderfile)) {
  print "Using orderfile\n";
  $ta->loadFile($orderfile);
  my $col         = $ta->getColumn(0);
  shift @$col if ($col->[0] eq "");
  @MOTIFS      = @$col;  # over-write
}



#
#  read in resmatrix file
#
$ta->loadFile($resmatrixfile);
my $a_ref_resmat = $ta->getArray();

my %RESMATRIX = ();
foreach my $r (@$a_ref_resmat) {
  my $a = shift @$r;
  my $b = shift @$r;
  #print "$a -- $b\n";
  $RESMATRIX{ $a }{ $b } = $r;
  $RESMATRIX{ $b }{ $a } = $r;
  
}


#
#  read in matrix file
#
$ta->loadFile($matrixfile);
my $a_ref_mat = $ta->getArray();

my %MATRIX = ();
foreach my $r (@$a_ref_mat) {
  my $a = shift @$r;
  my $b = shift @$r;
  #print "$a -- $b\n";
  $MATRIX{ $a }{ $b } = $r;
  $MATRIX{ $b }{ $a } = $r;
  
}



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



print "Now doing the graphical display.\n";

#
#  START DRAWING
#
#

my $xbase = 75;
my $ybase = 130;
my $w     = 20;
my $h     = 20;

my $xsize = $xbase + $w * scalar(@MOTIFS) + 135;
my $ysize = $ybase + $h * scalar(@MOTIFS) + 135; 


print "xsize = $xsize, ysize = $ysize, xbase = $xbase, ybase = $ybase\n";

my $p = new PostScript::Simple(xsize => $xsize, #papersize => "A4",
			       ysize => $ysize, #papersize => "A4",
			       #$papersize => "A4",
			       colour    => 1,
			       #landscape    => 1,
			       #coordorigin => "LeftTop",
			       #direction   => "RightDown",
			       eps       => 1,
			       units     => "pt");

$p->setlinewidth(0.5);
print $outfile;
open O, "> $outfile" or die;
my $mainwidth = $xbase+1000;
my $fontsize = $w*2/3;
print O (
	"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" 
   \"http://www.w3.org/TR/html4/strict.dtd\">\n"
);
print O ("<head>\n");
print O ("<title>FIRE results (motif interaction matrix)</title>\n");
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
text-align: center;
}
td.gridA {
font-family: courier, monospace;
text-align: center;
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
position: relative;
left: 50px;
top: -240px;
padding: 10px;
}
table.main tr {
font-size: $fontsize\px;
}
");
print O ("</style>\n");
print O ("</head>\n");
print O ("<body>\n\n");

print O ("<h1><abbr title=\"Finding Informative Regulatory Elements\">FIRE</abbr> results (motif interaction matrix)</h1>\n");

print O ("<p class=\"head\">These results are also available as <a href=\"./$fn.pdf\"><abbr title=\"Portable Document Format\">PDF</abbr></a> and <a href=\"./$fn.eps\"><abbr title=\"Encapsulated PostScript\">EPS</abbr></a> documents.</p>\n");

print O ("<p>Depending on your display resolution, scrolling or zooming may be necessary.</p>\n");

print O ('<p> In terms of cell color, a yellow cell represents positive correlation while a red cell represents negative correlation. In terms of cell borders, a black border represents functional module separation, a blue border represents significant interaction (p&lt;1e-4) between two DNA motifs, a pink border represents significant interaction between two RNA motifs, and a green border represents significant interaction between a DNA and RNA motif. The + symbol represents significant co-localization.</p>'."\n");

my $min  = -0.01;
my $max  =  0.01;
my $cntm = 0;
my $n    = scalar(@MOTIFS);

############ Scale ############
my $scalex = $xbase/3 ;
my $hscale = $h/2 ;
my $W = $hscale ;

print O ("<table class=\"legend\" cellpadding=\"0\" cellspacing=\"0\" style=\"width:$hscale\px\">\n");

print O "<tr><td style=\"font-size: 12px; height: 14px;\">(Pos)</td></tr>\n" ;
print O "<tr><td style=\"font-size: 12px; height: 14px;\">$max</td></tr>\n" ;

my $t = $max;
my $res = 70 ;
for (my $i=0; $i<=$res/2; $i++) {
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
  my $color = &RGBDecToHex($col[0], $col[1], $col[2]);
  print O "<tr>\n" ;
  print O "<td style=\"background-color:$color\">&nbsp;</td>\n" ;
  print O "</tr>\n" ;
  $t -= ($max - $min) / $res;
}
print O "<tr><td style=\"font-size: 12px; height: 14px;\">0.0</td></tr>\n" ;
print O "<tr><td style=\"background-color:#FFFFFF\">&nbsp;</td></tr>\n" ;
print O "<tr><td style=\"font-size: 12px; height: 14px;\">(Neg)</td></tr>\n" ;
print O "<tr><td style=\"font-size: 12px; height: 14px;\">$min</td></tr>\n" ;

my $t = $min;
for (my $i=$res ; $i>$res/2; $i--) {
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
  my $color = &RGBDecToHex($col[0], $col[1], $col[2]);
  print O "<tr>\n" ;
  print O "<td style=\"background-color:$color\">&nbsp;</td>\n" ;
  print O "</tr>\n" ;
  $t += ($max - $min) / $res;
}
print O "<tr><td style=\"font-size: 12px; height: 14px;\">0.0</td></tr>\n" ;

print O ("\n</table>\n");

print O ("<table class=\"main\" cellpadding=\"0\" cellspacing=\"0\">\n") ;

my $i = 0;
print O ("<tr>") ;
foreach my $re1 (@MOTIFS) {
  my $ty = ($STAT{$re1}->{RNA}==0?"5'":"3'UTR");

  my $p = new PostScript::Simple(xsize => $w/2,
				 ysize => 50,
				 colour    => 1,
				 eps       => 1,
				 units     => "pt");
  $p->setcolour("black");
  $p->setfont("Courier", 8);
  $p->text( { rotate => 90 }, 5, 5, "$ty") ;

  $ty =~ s/\ /_/gi ;
  $ty =~ s/'//gi ;
  $p->output("$d/$ty\_pos.eps") ;

  system("convert -density 100 $d/$ty\_pos.eps $d/$ty\_pos.png"); 
  print O "<td class=\"grid\" style=\"border:1px #FFF solid; background-color: #FFF;\" >" ;
  print O "<img src=\"./$dname/$ty\_pos.png\" alt=\"$ty\" /></td>\n" ;

  $i++;
}
print O ("</tr>") ;
print O ("<tr>") ;
my $i = 0;
foreach my $re1 (@MOTIFS) {
  my $myii = $MOTIF_NUMBER{ $re1 };
  print "$re1\t$myii\n";

  my $d  = $motifdir;  
  my $e  = new PostScript::Simple::EPS(file => "$d/$myii.eps");

  # get height
  my $eh = $e->height;
  my $ew = $e->width;
  $e->rotate(90) ;

  my $p = new PostScript::Simple(xsize => $eh,
				 ysize => $ew+5,
				 colour    => 1,
				 eps       => 1,
				 units     => "pt");
  $p->_add_eps($e,$ew*2/3,15);
  $p->output("$d/$myii\_r.eps") ;
  system("convert -density 100 $d/$myii\_r.eps $d/$myii\_r.png"); 

  print O "<td class=\"grid\"style=\"border:1px #FFF solid; background-color: #FFF;\">" ;
  print O "<img src=\"./$dname/$myii\_r.png\" width=\"$w\" alt=\"$re1\" title=\"$re1\" /></td>\n" ;
  
  $i++;
}
print O ("</tr>") ;
$i    = 0;
foreach my $re1 (@MOTIFS) {
  print O ("<tr>") ;
  my $rna1  = $STAT{$re1}->{RNA};

  my $d       = $motifdir;
  my $j       = 0;
  foreach my $re2 (@MOTIFS) {
    my $rna2  = $STAT{$re2}->{RNA};
    my $c = undef;
    
    my $r   = $RESMATRIX{ $re1 }{ $re2 }; 
    my $t   = undef;    

    if (defined($r)) {
      $t = $r->[1];
    } else {
      $t = 10000;
    }

    if ($re1 ne $re2) {
      my $r = $MATRIX{ $re1 }{ $re2 }; 
      if (abs($r->[3]) < 1e-10) {
	$c = 0.0;
      } else {
	$c = $r->[0] * $r->[3] / abs($r->[3]);
      }
    } else {
      $c = 1.0;
    }

    $c = 0 if ($c eq "nan") ;
      
    my @col = ();
    if (!defined($colmap)) {
      @col = Sets::interp_general( $c, [0, 0, 204], [255, 255, 0], $min, $max);
    } else {
      @col = Sets::interp_from_matlab_colormap( $c, $A_REF_COLMAP, $min, $max);
    }

    my $color = &RGBDecToHex($col[0], $col[1], $col[2]);
    my $style = "" ;
    if ($j==0){
      $style.="border-left:1px #000 solid; " ;
    }
    if ($h_ref_clust->{$re2} != $h_ref_clust->{$MOTIFS[$j+1]}) {
      $style.="border-right:1px #000 solid; " ;
    }
    if ($i==0){
      $style.="border-top:1px #000 solid; " ;
    }
    if ($h_ref_clust->{$re1} != $h_ref_clust->{$MOTIFS[$i+1]}) {
      $style.="border-bottom:1px #000 solid; " ;
    }

    my $content = "&nbsp;" ;
    my $bcol = "" ;
    if ($t == 0) {
      $style="" ;
      if ($rna1 != $rna2) {
	$bcol = "#5DFC0A" ;
      } else {
	if ($rna1 == 0) {
	  $bcol = "#0000FF" ;
	} else {
	  $bcol = "#FF66FF" ;
	}
      }
      $style="border:1px $bcol solid; " ;

      my $c1  = 0;
      if (($r->[4] ne "") && ($r->[4] ne "nan") && ($r->[4] <= 100)) { 
	$c1 = 1;
      }
      
      my $c2 = 0;
      if (($r->[6] ne "") && ($r->[6] ne "nan") && ($r->[6] <= 100)) { 
	$c2 = 1;
      }
  
      my $c3 = 0;
      if (($r->[8] ne "") && ($r->[8] ne "nan") && ($r->[8] <= 100)) { 
	$c3 = 1;
      }
      
      my $c4 = 0;
      if (($r->[10] ne "") && ($r->[10] ne "nan") && ($r->[10] <= 100)) { 
	$c4 = 1;
      }

      
      my $c5 = 0;
      if (($r->[12] ne "") && ($r->[12] ne "nan") && ($r->[12] <= 100)) { 
	$c5 = 1;
      }
      
      my $c6 = 0;
      if (($r->[14] ne "") && ($r->[14] ne "nan") && ($r->[14] <= 100)) { 
	$c6 = 1;
      }
      
      if ($c1 == 1) { # || ($c2 == 1)) {
	$content = "+" ;
      }
    }

    if ($c!=0){
      print O "<td class=\"grid\" style=\"$style background-color:$color;\" onclick=\"alert('$c')\" title=\"$c\">$content</td>\n" ;
    }else{
      print O "<td class=\"grid\" style=\"$style background-color:$color;\" >$content</td>\n" ;
    }

    $j ++;
  }
  my $myii = $MOTIF_NUMBER{ $re1 };
  #system("convert -density 100 $d/$myii.eps $d/$myii.png"); 

  print O "<td class=\"grid\" style=\"background-color: #FFF;\">" ;
  print O "<img src=\"./$dname/$myii.png\" height=\"$h\" alt=\"$re1\" title=\"$re1\" /></td>\n" ;

  print O "<td class=\"grid\" style=\"background-color:#FFFFFF\">" ;
  print O ($STAT{$re1}->{RNA}==1?"3'UTR":"5'")."</td>\n" ;  
  
  $i ++;
  print O ("</tr>") ;
}



$i    = 0;
foreach my $re1 (@MOTIFS) {
  my $j     = 0;
  my $rna1  = $STAT{$re1}->{RNA};

  foreach my $re2 (@MOTIFS) {
    my $rna2  = $STAT{$re2}->{RNA};

    #print "$re1 -- $re2\n";

    my $r   = $RESMATRIX{ $re1 }{ $re2 }; 
    my $t   = undef;    

    if (defined($r)) {
      $t = $r->[1];
    } else {
      $t = 10000;
    }

    
    
    
    #print "t = $t\n";
    
    if ($t == 0) {

      if ($rna1 != $rna2) {
	$p->setcolour("green");
      } else {
	if ($rna1 == 0) {
	  $p->setcolour("blue");
	} else {
	  $p->setcolour((255, 102, 255));	  
	}
      }

     
      $p->box({filled => 0}, 
	      $xbase + $j * $w ,      $ysize - ($ybase + $i*$h ) , 
	      $xbase + $j * $w + $w , $ysize - ($ybase + ($i*$h+$h) ));


      # BEWARE: used the shift function twice above !

      my $c1  = 0;
      if (($r->[4] ne "") && ($r->[4] ne "nan") && ($r->[4] <= 100)) { 
	$c1 = 1;
      }
      
      my $c2 = 0;
      if (($r->[6] ne "") && ($r->[6] ne "nan") && ($r->[6] <= 100)) { 
	$c2 = 1;
      }
  
      my $c3 = 0;
      if (($r->[8] ne "") && ($r->[8] ne "nan") && ($r->[8] <= 100)) { 
	$c3 = 1;
      }
      
      my $c4 = 0;
      if (($r->[10] ne "") && ($r->[10] ne "nan") && ($r->[10] <= 100)) { 
	$c4 = 1;
      }

      
      my $c5 = 0;
      if (($r->[12] ne "") && ($r->[12] ne "nan") && ($r->[12] <= 100)) { 
	$c5 = 1;
      }
      
      my $c6 = 0;
      if (($r->[14] ne "") && ($r->[14] ne "nan") && ($r->[14] <= 100)) { 
	$c6 = 1;
      }
            
      if ($c1 == 1) { # || ($c2 == 1)) {
	
	$p->setcolour("black");
	$p->line(  $xbase + $j * $w + $w/4,
		   $ysize - ($ybase +  $i*$h + $h/2    ) ,
		   $xbase + $j * $w + 3*$w/4,
		   $ysize - ($ybase +  $i*$h + $h/2    ));
	
	$p->line(  $xbase + $j * $w + $w/2,
		   $ysize - ($ybase +  $i*$h + $h/4    ) ,
		   $xbase + $j * $w + $w/2,
		   $ysize - ($ybase +  $i*$h + 3*$h/4    ));
      }
      
      
    }
    $j ++;
    
  }
  $i++;
}

print O ("\n</table>\n");
print O ("\n</body>\n");
print O ("</html>\n");


sub RGBDecToHex {
  my ($red, $green, $blue) = @_;
  return sprintf("#%02lx%02lx%02lx", $red, $green, $blue);
}
