#__perl__
#


use strict;

my $p=0;
my $method = undef;
while (<>) {
	if (/^\=/) {
		$p=0;
		if ($method eq "classification_of") {
			print "        '''\n\n";
			print <<EOF;
        self.function('classification_of')
        ss = self.retrieve(data)
        for k in ss:
            if (len(ss[k]) == 1):
                print(k)
                ss[k].append(None)
        return ss
EOF
		undef $method;
		}
		elsif ($method) {
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


