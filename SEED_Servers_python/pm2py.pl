#__perl__
#


use strict;

my $p=0;
my $method = undef;
while (<>) {
	if (/^\=/) {
		$p=0;
		if ($method) {
			print "        '''\n\n";
			print<<EOF;
        self.function('$method')
        return self.retrieve(data)

EOF
		undef $method;
		}
	}
	if (s/^\=head3\s+//) {
		$p=1; 
		$method = $_;
		chomp($method);
		print "    def $method(self,data):\n";
		print "        '''\n";
		next;
	}
	print "        $_" if ($p && $_ !~ /^\s+/);

}


