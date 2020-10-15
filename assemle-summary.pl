#! /usr/bin/perl

#summary the contig status of an assembly

$scaffold = shift;
$large = shift;
if($large eq "") {
    $large = 1000;
}

open (SCAF, "<$scaffold") or die;


while (<SCAF>) {
    if($_ =~ /^>(\S+)/) {
	$name = $1;
	$names[@names] = $name;
	$id = @names -1;
    } else {
	chomp $_;
	$seq{$name} .= $_;
    }
}

@names = sort {length($seq{$b}) <=> length($seq{$a})} @names;

$largest_scaf = $names[0];
$largest_length = length($seq{$names[0]});

for(local $i=0; $i < @names; ++$i) {

    @contseq = split /N+/, $seq{$names[$i]};
    for(local $n=0; $n<@contseq; ++$n) {
	$seq2{"$names[$i]-$n"} = $contseq[$n];
	$names2[@names2] = "$names[$i]-$n";
    }

    local $at = 0; $gc = 0; $basen = 0;
    for(local $j=0; $j< length($seq{$names[$i]}); ++$j) {
	$base = substr($seq{$names[$i]}, $j, 1);
	if($base =~ /[ATat]/) {
	    $at ++;
	} elsif ($base =~ /[GCgc]/) {
	    $gc ++;
	} else {
	    $basen ++;
	}
    }

    if(length($seq{$names[$i]}) >= $large) {
	$no_large_scaf ++;
	$len_large_scaf += length($seq{$names[$i]});
	$at_large_scaf += $at;
	$gc_large_scaf += $gc;
	$basen_large_scaf += $basen;
    }
    $no_scaf ++;
    $len_scaf += length($seq{$names[$i]});
    $at_scaf += $at;
    $gc_scaf += $gc;
    $basen_scaf += $basen;    
}

for(local $i=0; $i<@names; ++$i) {
    $length += length($seq{$names[$i]});
    if($length >= 0.5 * $len_large_scaf && $N50_large eq "") {
	$N50_large = length($seq{$names[$i]});
	$N50_large_num = $i+1;
    }
    if($length >= 0.9 * $len_large_scaf && $N90_large eq "") {
	$N90_large = length($seq{$names[$i]});
	$N90_large_num = $i+1;
    }
    if($length >= 0.5 * $len_scaf && $N50 eq "") {
	$N50 = length($seq{$names[$i]});
	$N50_num = $i+1;
    }
    if($length >= 0.9 * $len_scaf && $N90 eq "") {
	$N90 = length($seq{$names[$i]});
	$N90_num = $i+1;
    }
}


print "$scaffold\nLarge scaffolds (>$large bps):\n  Largest scaffold: $largest_scaf\n  Largest length: $largest_length\n  No. of large scaffolds: $no_large_scaf\n  Bases in large scaffolds: $len_large_scaf\n  N50 scaffold: $N50_large_num\n  N50 length: $N50_large\n  N90 scaf: $N90_large_num\n  N90 length: $N90_large\n  GC content: ".(int(($gc_large_scaf*100000)/($gc_large_scaf+$at_large_scaf))/1000)."%\n  N rate: ".(int(($basen_large_scaf*100000)/$len_large_scaf)/1000)."%\n\n";

print "All scaffolds:\n  Number: $no_scaf\n  Total base: $len_scaf\n  GC content: ".(int(($gc_scaf*100000)/($gc_scaf+$at_scaf))/1000)."%\n  N rate: ".(int($basen_scaf*100000/$len_scaf)/1000)."%\n";

=pod
@names2 = sort {length($seq2{$b}) <=> length($seq2{$a})} @names2;

$largest_contig = $names2[0];
$largest_length = length($seq2{$names2[0]});

$N50_large_num="";
$N50_large = "";
$N50_length = "";
$N50_num="";
$N50 = "";
$N90_large_num="";
$N90_large = "";
$N90_length="";
$N90_num="";
$N90 = "";
$length=0;
$at = 0;
$gc = 0;
$basen = 0;

for(local $i=0; $i < @names2; ++$i) {

    local $at = 0; $gc = 0; $basen = 0;
    for(local $j=0; $j< length($seq2{$names2[$i]}); ++$j) {
	$base = substr($seq2{$names2[$i]}, $j, 1);
	if($base =~ /[ATat]/) {
	    $at ++;
	} elsif ($base =~ /[GCgc]/) {
	    $gc ++;
	} else {
	    $basen ++;
	}
    }

    if(length($seq2{$names2[$i]}) >= $large) {
	$no_large_contig ++;
	$len_large_contig += length($seq2{$names2[$i]});
	$at_large_contig += $at;
	$gc_large_contig += $gc;
	$basen_large_contig += $basen;
    }
    $no_contig ++;
    $len_contig += length($seq2{$names2[$i]});
    $at_contig += $at;
    $gc_contig += $gc;
    $basen_contig += $basen;
}

for(local $i=0; $i<@names2; ++$i) {
    $length += length($seq2{$names2[$i]});
    if($length >= 0.5 * $len_large_contig && $N50_large eq "") {
	$N50_large = length($seq2{$names2[$i]});
	$N50_large_num = $i+1;
    }
    if($length >= 0.9 * $len_large_contig && $N90_large eq "") {
	$N90_large = length($seq2{$names2[$i]});
	$N90_large_num = $i+1;
    }
    if($length >= 0.5 * $len_contig && $N50 eq "") {
	$N50 = length($seq2{$names2[$i]});
	$N50_num = $i+1;
    }
    if($length >= 0.9 * $len_contig && $N90 eq "") {
	$N90 = length($seq2{$names2[$i]});
	$N90_num = $i+1;
    }
}

print "\nLarge contigs (>$large bps):\n  Largest contig: $largest_contig\n  Largest length: $largest_length\n  No. of large contigs: $no_large_contig\n  Bases in large contigs: $len_large_contig\n  N50 contig: $N50_large_num\n  N50 length: $N50_large\n  N90 contig: $N90_large_num\n  N90 length: $N90_large\n  GC content: ".(int(($gc_large_contig*100000)/($gc_large_contig+$at_large_contig))/1000)."%\n\n";

print "All contigs:\n  Number: $no_contig\n  Total base: $len_contig\n  GC content: ".(int(($gc_contig*100000)/($gc_contig+$at_contig))/1000)."%\n\n";
=cut
