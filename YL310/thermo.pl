#!/usr/bin/perl
# thermo.pl
 
$pi = 3.1415926535;
$avogad = 6.0221417930E23;
@specialt = ( 298.15 );
if ( $#ARGV < 0 ) {
	printf( "Usage:  perl thermo.pl [g03freq.out|datafile.dat] <freq scale factor> <temperature>\n" );
	exit;
}
open( FILE, $ARGV[0] );
printf( "opening %s for input\n\n", $ARGV[0] );
$scale = 1.0;
$symno = 1;
# arguments can be scale factor (in range 0.5-1.5) or temperature
for ( $i = 1; $i <= $#ARGV; $i++ ) {
	$tt = $ARGV[ $i ];
	if ( abs( $tt - 1.0 ) <= 0.5 ) {
		# vibrational scaling factor -- don't allow too far from 1.0
		$scale = $tt;
	} elsif ( $tt >= 10. ) {
		# a temperature -- don't allow them too ridiculously low
		push( @specialt, $tt );
	}
}
$atom = $linear = 0;
$g98 = $datfile = 0;	# type of input data file
printf( "Scaling vibrational frequencies by %.4f\n", $scale );
$ifreq = $ielec = $ixvib = 0;
$zpe = $xzpe = 0;
while ( <FILE> ) {
  #
  # this block is for G98 frequency output file
  #
  if ( /Molecular mass:/ ) {
    @line = split();
	$mass = $line[2];
	$g98 = 1;	# flag for g98 input file
  }
  if ( /ROTATIONAL CONSTANT/i ) {
    @line = split();
	$rot[0] = $line[3];
	$rot[1] = $line[4];
	$rot[2] = $line[5];
  }
  if ( /Full point group/ ) {
    @line = split();
	$group = $line[3];
	if ( $group =~ /\*/ ) {
		# molecule is linear 
		$linear = 1;
	}
	if ( $group eq "KH" ) {
		# molecule is an atom
		$atom = 1;
	}
  }
  if ( /Frequencies/ ) {
    @line = split();
	for ( $i = 2; $i < 5; $i++ ) {
		$tt = $line[$i];
		if ( $tt > 0.0 ) {
			$freq[$ifreq++] = $tt;
		} else {
			printf "\tDiscarding vibrational frequency of %.3f\n", $tt;
		}
	}
  }
  if ( /Charge/ && /Multiplicity/ ) {
	# electronic degeneracy
	$energy[0] = 0.0;
	$degen[0] = ( split() )[5];
	$ielec = 1;
  }
  if ( /ROTATIONAL SYMMETRY NUMBER/i ) {
    @line = split();
	$symno = $line[3];
  }
  #
  # this block is for keyword file
  #
  if ( /^MASS/i ) {
	@line = split();
	$mass = $line[1];
	$datfile = 1;	# flag for keyword input file
  }
  if ( /^GHZ/i ) {
	# rotational constant(s) in GHz
	@line = split();
	if ( $#line < 3 ) {
		$linear = 1;
		$rot[0] = $line[1];
		if ( $rot[0] == 0 ) {
			# molecule is an atom
			$linear = 0;
			$atom = 1;
		}
	}
	else {
		for ( $i = 0; $i < 3; $i++ ) {
			$rot[$i] = $line[$i + 1];
		}
	}
  }
  if ( /^VIB/i ) {
	@line = split();
	for ( $i = 1; $i <= $#line; $i++ ) {
		$freq[$ifreq++] = $line[$i];
	}
  }
  if ( /^XVIB/i ) {
  	# explicit vibrational energy levels (e.g., from anharmonic calc.)
	# note: must include the zero!
	( $i, @{$xvib[$ixvib]} ) = split();	# ignore i
	$ixvib++;
  }
  if ( /^AVIB/i ) {
  	# anharmonic oscillator: E(v) = v[w - (v+1)x]
	# two or three fields:	w, x, and wavenumber limit
	@line = split();
	push( @w, $line[1] );
	push( @x, $line[2] );
	push( @lim, $line[3] );
  }
  if ( /^ELEC/i ) {
	# electronic level info:  degeneracy, then energy in cm^-1
	$degen[$ielec] = 1;
	$energy[$ielec] = 0.0;
	@line = split();
	$degen[$ielec] = $line[1];
	$energy[$ielec++] = $line[2];
  }
  if ( /^SYMNO/i ) {
	# symmetry number
	@line = split();
	$symno = $line[1];
  }
}
printf( "Molecular mass = %.3f u\n", $mass );
print "External symmetry number = $symno\n";
if ( $linear ) {
	$rot[0] = $rot[2] unless ( $rot[0] > 0 );	# use last if first is zero
	printf( "Rotational constant (GHz) = %.5f\n", $rot[0] );
	$#rot = 0;	# truncate any useless elements
} elsif ( $atom ) {
	print "Molecule is an atom.\n";
} else {
	printf( "Rotational constants (GHz) = %.5f %.5f %.5f\n", $rot[0], $rot[1], $rot[2] );
}
if ( $ifreq > 0 ) {
	print "Harmonic vibrational frequencies (cm^-1):\n";
	print "\tunscaled\tscaled\n";
} else {
	print "Harmonic vibrational frequencies (cm^-1):  none\n";
}
for ( $i = 0; $i < $ifreq; $i++ ) {
	$tt = $freq[$i];
	$freq[$i] *= $scale;
	printf( "%3d\t%6.3f\t\t%6.3f\n", $i + 1, $tt, $freq[$i] );
	$zpe += $freq[$i] / 2;	# scaled harmonic contribution to vibrational zero-point energy
}

unless ( $datfile ) {
	$reply = "yes";
	if ( -e "thermo.dat" ) {
		print "\n***Overwrite file 'thermo.dat'? ";
		chomp( $reply = yes );
	}
	if ( $reply =~ /^y/i ) {
		open( DAT, ">thermo.dat" ) or warn "Unable to open 'thermo.dat' for output!\n";
		printf DAT "MASS\t%.6f\n", $mass;
		printf DAT "GHZ" . ( "\t%.6f" x @rot ) . "\n", @rot;
		if ( $ifreq > 0 ) {
			print DAT "VIB\t";
			for ( $i = 0; $i < $ifreq; $i++ ) {
				printf DAT "%.4f ", $freq[$i];	# these are scaled frequencies
			}
			print DAT "\n";
		}
		printf DAT "SYMNO\t%d\n", $symno;
		if ( $ielec > 0 ) {
			for ( $i = 0; $i < $ielec; $i++ ) {
				printf DAT "ELEC\t%d\t%.4f\n", $degen[$i], $energy[$i];
			}
		}
		close( DAT );
	}
}
if ( $ixvib > 0 ) {
	print "Explicit vibrational ladders (cm^-1):\n";
	for ( $i = 0; $i < $ixvib; $i++ ) {
		$fmt = ( " %.4f" ) x ( $#{$xvib[$i]}+1 );
		printf( "%3d  $fmt\n", $i+1, @{$xvib[$i]} );
		$zpe += $xvib[$i][0];	# usually will be zero
		$xzpe += $xvib[$i][0];	# keep track of the explicit ZPE
	}
}
if ( $#w > -1 ) {
	print "Anharmonic oscillators (cm^-1):\n";
	for ( $i = 0; $i <= $#w; $i++ ) {
		$lim[$i] = 10000 if ( $lim[$i] <= 0 );	# no limit requested
		$vmax = $w[$i] / ( 2 * $x[$i] ) - 1;		# v at turnover point
		$t = $vmax * ( $w[$i] - ( $vmax + 1 ) * $x[$i] );
		# don't exceed turnover point and don't go negative
		$lim[$i] = $t if ( $t < $lim[$i] and $t > $w[$i] );
		printf( "%3d\tw = %.4f\tx = %.4f\tlimit = %.4f\n", $i+1, $w[$i], $x[$i], $lim[$i] );
		$zpe += $w[$i] / 2 - $x[$i] / 4;	# anharmonic ZPE
	}
}
if ( $ielec > 0 ) {
	print "Electronic states:\n";
	printf( "\tdegen\tcm^-1\n" );
	for ( $i = 0; $i < $ielec; $i++ ) {
		printf( "%3d\t%3d\t%.3f\n", $i+1, $degen[$i], $energy[$i] );
	}
}
printf( "\n" );

########## begin calculations ###########
$R = 8.314510;		# J/(mol K)
$h = 6.626076E-34;	# J s
$k = 1.38066E-23;	# J/K
$p = 1.0e5;			# Pa
$c = 2.99792458E08;	# m/s

printf "Vibrational zero-point energy = %.4f cm^-1 ", $zpe;
$zpe *= $h * $avogad * $c / 10;		# convert from cm^-1 to kJ/mol
printf "= %.4f kJ/mol.\n", $zpe;
if ( $xzpe > 0 ) {
	printf "(includes %.4f cm^-1 = ", $xzpe;
	$xzpe *= $h * $avogad * $c / 10;	# convert from cm^-1 to kJ/mol
	printf "%.4f kJ/mol from explicit vibrational levels)\n";
	print "***WARNING: explicit vibrational ladders should generally start at zero***\n";
} 

# array indices: 0 = translation, 1 = rotation, 2 = vibration, 3 = electronic
for ( $i = 0; $i < 4; $i++ ) {
	$s[$i] = $cp[$i] = $ddh[$i] = 0.0;
}
# translational contributions
# convert mass from u to kg
$mass /= $avogad * 1000.0;
$s[0] = $R * ( 1.5 * log( 2 * $pi * $mass / $h / $h ) - log( $p ) + 2.5 );
$cp[0] = 2.5 * $R;
$ddh[0] = 2.5 * $R;
# rotational contributions
for ( $i = 0; $i < 3; $i++ ) {
	# convert from GHz to s^-1
	$rot[$i] *= 1.0e09;
}
if ( $rot[0] > 0.0 ) {
	# don't try to include rotations for atoms
	if ( $linear == 1 ) {
		# linear molecule
		$s[1] = $R * ( log( $k / ( $symno * $h * $rot[0] ) ) + 1.0 );
		$cp[1] = $R;
		$ddh[1] = $R;
	}
	else {
		# non-linear molecule
		$s[1] = log( $rot[0] ) + log( $rot[1] ) + log( $rot[2] ) - log( $pi );
		$s[1] *= -0.5;
		$s[1] += 1.5 - log( $symno );
		$s[1] *= $R;
		$cp[1] = 1.5 * $R;
		$ddh[1] = 1.5 * $R;
	}
}
$convert = 100.0 * $c * $h / $k;	# conversion factor from cm^-1 to K
foreach ( @freq ) { $_ *= $convert; }	# convert harmonic frequencies
foreach ( @energy ) { $_ *= $convert; }	# convert electronic energies
# convert explicit vibrational ladders
for ( $i = 0; $i < $ixvib; $i++ ) {
	foreach ( @{ $xvib[$i] } ) { $_ *= $convert; }
}
# convert anharmonic frequencies
for ( $i = 0; $i <= $#w; $i++ ) {
	$w[$i] *= $convert;
	$x[$i] *= $convert;
	$lim[$i] *= $convert;
}
# loop over temperatures and print totals
printf( "\nT (K)\tS (J/mol.K)\tCp (J/mol.K)\tddH (kJ/mol)\n" );
$tstep = 20;
@tlist = ();
for ( $t = 20; $t <= 4000; $t += $tstep ) {
	push( @tlist, $t );
}
foreach $t ( @specialt ) {
	for ( $i = 0; $i < $#tlist; $i++ ) {
		if ( $tlist[$i] < $t and $t < $tlist[$i + 1] ) {
			# insert the special temperature here
			splice( @tlist, $i+1, 0, $t );
			last;
		}
	}
	push( @tlist, $t ) if ( $t > $tlist[$#tlist] );
	unshift( @tlist, $t ) if ( $t < $tlist[0] );
}
foreach $t ( @tlist ) {
	# translational contribution
	$stot = $s[0] + $R * 2.5 * log( $k * $t );
	#printf( "Strans = %.2f\t", $stot );
	if ( $linear == 1 ) {
		$tt = $s[1] + $R * log ( $t );
		$stot += $tt;
	}
	elsif ( $rot[0] > 0.0 ) {
		$tt = log( $k * $t ) - log( $h );
		$tt *= 1.5 * $R;
		$tt += $s[1];
		$stot += $tt;
	}
	#printf( "Srot = %.2f\t", $tt );
	$s[2] = $cp[2] = $ddh[2] = 0.0;
	# harmonic vibrational contributions
	for ( $i = 0; $i < $ifreq; $i++ ) {
		next if ( $freq[$i] <= 0 );	# skip "negative" and zero frequencies
		$tt = $freq[$i] / $t;
		$ttt = exp( - $tt );
		$s[2] -= log( 1.0 - $ttt );
		$s[2] += $tt * $ttt / ( 1.0 - $ttt );
		$cp[2] += $tt * $tt * $ttt / ( 1.0 - $ttt ) / ( 1.0 - $ttt );
		$ddh[2] += $tt * $ttt / ( 1.0 - $ttt );
	}
	#printf( "Svib = %.2f\t", $s[2] * $R );
	$s[3] = $cp[3] = $ddh[3] = 0.0;
	# electronic level contributions
	if ( $ielec > 0 ) {
		$s1 = $s2 = $s3 = 0.0;
		for ( $i = 0; $i < $ielec; $i++ ) {
			$tt = $energy[$i] / $t;
			$ttt = exp( -$tt );
			$s0 = $degen[$i] * $ttt;
			$s1 += $s0;
			$s0 *= $tt;
			$s2 += $s0;
			$s0 *= $tt;
			$s3 += $s0;
		}
		$s[3] = log( $s1 ) + $s2 / $s1;
		#printf( "Selec = %.2f\t", $s[3] * $R );
		#$cp[3] = $s3 + ( $s2 * $s2 / $s1 / $s1 );
		$cp[3] = $s3 / $s1 - $s2 * $s2 / ( $s1 * $s1 );
		$ddh[3] = $s2 / $s1;
	}
	# explicit vibrational ladders
	for ( $j = 0; $j < $ixvib; $j++ ) {
		$s1 = $s2 = $s3 = 0.0;
		for ( $i = 0; $i <= $#{ $xvib[$j] }; $i++ ) {
			$tt = $xvib[$j][$i] / $t;
			$s0 = exp( -$tt );
			$s1 += $s0;
			$s0 *= $tt;
			$s2 += $s0;
			$s0 *= $tt;
			$s3 += $s0;
		}
		$s[3] += log( $s1 ) + $s2 / $s1;
		$cp[3] += $s3 / $s1 - $s2 * $s2 / ( $s1 * $s1 );
		$ddh[3] += $s2 / $s1;
	}
	# anharmonic vibrations
	for ( $i = 0; $i <= $#w; $i++ ) {
		$e = $s1 = $s2 = $s3 = 0.0;
		for ( $v = 1; $e <= $lim[$i]; $v++ ) {
			$tt = $e / $t;
			$s0 = exp( -$tt );
			$s1 += $s0;
			$s0 *= $tt;
			$s2 += $s0;
			$s0 *= $tt;
			$s3 += $s0;
			$e = $v * ( $w[$i] - ( $v + 1 ) * $x[$i] );
		}
		$s[3] += log( $s1 ) + $s2 / $s1;
		$cp[3] += $s3 / $s1 - $s2 * $s2 / ( $s1 * $s1 );
		$ddh[3] += $s2 / $s1;
	}
	$stot += $R * ( $s[2] + $s[3] );
	$cptot = $cp[0] + $cp[1] + $R * ( $cp[2] + $cp[3] );
	$ddhtot = $ddh[0] * $t;
	$ddhtot += $ddh[1] * $t;
	$ddhtot += $R * $t * ( $ddh[2] + $ddh[3] );
	$ddhtot /= 1000.0;
	printf( "%.2f\t\t%.3f\t\t%.3f\t\t%.3f\n", $t, $stot, $cptot, $ddhtot );
}
if ( $reply eq "yes" ) {
	print "\nInput data (including scaled vibrational frequencies) written to\n";
	print "file 'thermo.dat'.\n";
}
# If appropriate, print warning about undetected spatial degeneracy
if ( $degen[0] > 1 and ( $linear or $rot[0] == 0.0 or $symno > 2 ) ) {
	# electronic degeneracy AND
	#  linear molecule, atom, or possible non-Abelian point group
	print "\n**WARNING**\nThis program does not detect electronic spatial degeneracies.\n";
	print "If they apply to this molecule, save and edit 'thermo.dat' and then run:\n";
	print "\tperl thermo.pl thermo.dat\n";
}
