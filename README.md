Policy-controlled Signature From NTRU Lattice
===========

This software is a proof-of-concept implementation of policy-controlled signature scheme over NTRU lattices, described in the paper "Policy-controlled Signature from NTRU Lattice", of Zi-Yuan Liu, Jen-Chieh Hsu, Raylin Tso

Warning
=======
This code is not to be considered secure, efficient or fully portable. Its purpose is not to be used for actual encryption, but to provide the research community a tool to verify, analyze and reproduce the statements made in our paper.

How to use?
===========

To modify the parameters, edit the values N0 and q0 in params.h.

To run on an Unix machine with g++:
```
$ make
$ ./PCS_NTRU
```

If GMP and NTL are not in a standard directory, you have to modify the CCFLAGS and LDFLAGS in the Makefile to indicate where they are.


Reference
=========

This project is based on Thomas Léo Ducas, Vadim Lyubashevsky and Thomas Prest's work "Efficient Identity-Based Encryption over NTRU Lattices". For more information, see https://github.com/tprest/Lattice-IBE

