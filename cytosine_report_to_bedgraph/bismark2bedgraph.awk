#!/usr/bin/awk

BEGIN {
	FS = "\t";
	if ( coverage == "true" ) {
		printf("#<Chr>\t<Start>\t<End>\t<Coverage>\n")
	} else { 
		printf("#<Chr>\t<Start>\t<End>\t<Strand;Meth_ratio>\n")
	}
}

{
	if ( $6 ~ context && ($4 > 0 || $5 > 0) ) { 

		chr_name = $1;
		chr_pos = $2;
		strand = $3;
		c_meth_count = $4;
		c_unmeth_count = $5;

		if ( coverage == "true" ) {
			nb_reads = c_meth_count + c_unmeth_count
			printf("%s\t%s\t%s\t%s\n", chr_name, chr_pos, chr_pos, nb_reads)
		} else {
			if ( strand == "-") { 
				s = "-"; 
			} else { 
				s = "";
			}
			meth_ratio = c_meth_count / (c_meth_count + c_unmeth_count); 
			printf("%s\t%s\t%s\t%s%s\n", chr_name, chr_pos, chr_pos, s, meth_ratio)	
		}
	}
}