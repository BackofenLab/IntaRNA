
# steps when releasing a new version:

- increase counter in configure.ac

- create distribution source code tar ball
  - extract tar ball and check if compiling and tests working

- prepare release on github:
  - fill changelog
  - upload source code tar ball
  - (upload win binary)
  - (upload linux binary)
	- (via cygwin) to get cygwin dlls used for compilation on Windows run
	  - `for f in `ldd ./IntaRNA | grep "/usr/bin/" | awk '{print $3}'`; do cp $f .; done`
	  - zip these together with `IntaRNA.exe` (and additional script files etc.)
  - (upload API docu pdf and html.zip)
  
  
- publish release on github
  
- update [bioconda recipe](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/intarna)
- update [galaxy tool shed information](https://github.com/bgruening/galaxytools/tree/master/tools/rna_tools/intarna/intarna.xml)

- update webserver

- tweet release 
