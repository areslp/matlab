# convert meshlab OBJ that is actually
# pointset (only v/vn records) with estimated normals 
# to .pts format (which is just .obj v/vn format). 
# Basically MeshLab leaves out some normals, so we
# fill them in assuming that the normal is like a state
# (ie vn 'sets' the current normal)

open FILE,$ARGV[0];

$lastwasvn = "";
$lastvn = "";
while (<FILE>) {
    @data = split;
    next if $data[0] eq '#';
    if ($data[0] eq "vn") {
	$lastvn = "$data[1] $data[2] $data[3]";
	print "vn $lastvn \n";
	$lastwasvn = 1;
    } elsif ($data[0] eq "v") {
	if (! $lastwasvn) {
	    print "vn $lastvn \n";
	}
	print "v $data[1] $data[2] $data[3]\n";
	$lastwasvn = 0;
    }
}    
