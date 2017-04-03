use 5.008001;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'Schedule::SGE',
    VERSION_FROM      => 'SGE.pm', # finds $VERSION
    PREREQ_PM         => {}, # e.g., Module::Name => 1.1
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT       => 'Schedule::SGE is a suite of modules for interacting with the Sun Grid Engine. The base module Schedule::SGE handles locating the executables and making sure everything works fine. The three modules Schedule::SGE::Run, Schedule::SGE::Control, and Schedule::SGE::Status are for different interactions with the queues',
       AUTHOR         => 'rob <Rob@nonet>') : ()),
);
