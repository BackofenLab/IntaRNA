
# steps when releasing a new version:

- increase counter in configure.ac

- create distribution source code tar ball
- prepare release on github:
  - upload source code tar ball
  - fill changelog
  - (upload win binary)
  - (upload linux binary)
  
- tag release on github
  
- update bioconda recipe
- update galaxy tool shed information

- update binary on webserver

- tweet release 
