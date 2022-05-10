# Turbospectrum2019
synthetic spectrum calculation companion to MARCS

Please acknowledge use of this software through citation of 
Plez, B., 2012, Astrophysics Source Code Library, record ascl:1205.004

see: http://adsabs.harvard.edu/abs/2012ascl.soft05004P

In December 2020, Stark broadening was added, and a new release added: v19.1.3

New release on May 28, 2021 with a few small bug corrections: v19.1.4

2022-Mar-24: Modification of a couple of read statements in bsyn.f avoiding a crash with older ifort versions. 
No new release created. You may update by downloading the latest version of the routine.

2022-Apr-20: Update of Utilities/vald3line-BPz-freeformat.f, a code that translates line list from the VALD (http://vald.astro.uu.se) database to the TS format

2022-May-10: The modifications of March 24 were removed, due to problems arising with modern compilers. A couple of small changes added elsewhere.

Further development of the code with NLTE capability (v20) is currently made in a private repository, soon to become public (2022). 
The v19 public repository has not been archived, and is open to improvements by anyone.
